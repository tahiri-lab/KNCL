# Evaluation of the effects of the parameter k in the k-NCL tree completion algorithm

# Computes BSD(k–NCL) distances on kNCL-completed trees and produces

import os
import math
import pickle
import numpy as np
from ete3 import Tree
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
import multiprocessing as mp

# Import kNCL from a Python module
# Ensure kncl.py containing the kNCL implementation is in the same folder
from kncl import kNCL

# CONFIG
# If True, compute exact BSD values (no early exit).
# If False, BSD loop uses an early-exit bound for speed.
STRICT_IDENTICAL_OUTPUT = True

# Parallelism:
#   None or 0  -> use os.cpu_count()
#   >0         -> exact number of workers
#   <0         -> (cpu_count + PARALLEL_JOBS), e.g., -1 uses all but one core
PARALLEL_JOBS = -1


# Plotting setup

def _setup_matplotlib():
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True  # LaTeX rendering
    import matplotlib.pyplot as plt
    return plt


# Helpers

def resolve_n_jobs(n_jobs):
    """Resolve n_jobs allowing negatives (relative to CPU count)."""
    cpu = os.cpu_count() or 1
    if n_jobs is None or n_jobs == 0:
        return cpu
    if n_jobs < 0:
        return max(1, cpu + n_jobs)
    return n_jobs


def _running_in_notebook():
    # Heuristic: if running under IPython kernel in Jupyter
    return ("IPYKERNEL" in os.environ) or ("JPY_PARENT_PID" in os.environ)


def _effective_workers(n_jobs):
    """
    Return the actual number of worker processes to use.
    In Windows+Jupyter we explicitly compute a positive 'all but one' number
    instead of falling back to single-process.
    """
    if os.name == "nt" and _running_in_notebook():
        # force parallelism in notebook: all but one logical core
        cpu = os.cpu_count() or 1
        return max(1, cpu - 1)
    # otherwise, respect the normal resolution (supports negatives like -1)
    return resolve_n_jobs(n_jobs)


# Data loading

def load_datasets(filenames):
    datasets = []
    for file in filenames:
        with open(file, 'r') as f:
            trees = [Tree(line.strip(), format=1) for line in f if line.strip()]  # Newick per line
            datasets.append(trees)
    return datasets


# Function to compute overlap between two trees (Jaccard)
def calculate_jaccard_overlap(tree1, tree2):
    taxa1 = set(tree1.get_leaf_names())
    taxa2 = set(tree2.get_leaf_names())
    intersection = len(taxa1 & taxa2)
    union = len(taxa1 | taxa2)
    return intersection / union if union > 0 else 0.0



# For fast BSD(k–NCL) calculations

def squared_distance_sum_capped(t1, t2, leaves, cap=None):
    """
    Sum_{(i<j)} (d1(i,j) - d2(i,j))^2 with optional early exit.
    If `cap` is provided (BSD upper bound), abort as soon as the running
    squared sum exceeds cap**2.
    """
    cap2 = None if cap is None or not np.isfinite(cap) else cap * cap
    total = 0.0
    leaf_list = list(leaves)
    n = len(leaf_list)
    for idx_a in range(n - 1):
        a = leaf_list[idx_a]
        for idx_b in range(idx_a + 1, n):
            b = leaf_list[idx_b]
            d1 = t1.get_distance(a, b)
            d2 = t2.get_distance(a, b)
            diff = d1 - d2
            total += diff * diff
            if cap2 is not None and total > cap2:
                return total  # early exit: can't beat current best
    return total


def compute_bsd_kncl(T1, T2, k, cap=None):
    """
    Run kNCL(T1, T2, k), then compute BSD(k–NCL) distance between the completed trees.
    Returns: bsd_distance (float) or None on failure.
    """
    try:
        # Copy to avoid in-place modifications upstream
        T1_completed, T2_completed = kNCL(T1.copy(), T2.copy(), k)
    except Exception as e:
        print(f"Error in kNCL function: {e}")
        try:
            print(f"Tree 1 with error: {T1.write(format=1)}")
            print(f"Tree 2 with error: {T2.write(format=1)}")
        except Exception:
            pass
        return None

    try:
        leaves_completed = set(leaf.name for leaf in T1_completed.get_leaves())
        cap_to_use = None if STRICT_IDENTICAL_OUTPUT else cap
        bsd_sq = squared_distance_sum_capped(
            T1_completed, T2_completed, leaves_completed, cap=cap_to_use
        )
        if bsd_sq < 0:
            print("Warning: Squared distance sum is negative. Returning None for BSD.")
            return None
        return math.sqrt(bsd_sq)
    except Exception as e:
        print(f"Error computing BSD: {e}")
        return None


# Worker for parallel execution

def _evaluate_pair_worker(args):
    """
    Worker function run in a separate process.
    args: (newick1, newick2, N_cl, overlap)
    Returns:
      (overlap_bin, per_k_bsd, opt_k_bsd, overlap) or None
    """
    newick1, newick2, N_cl, overlap = args
    try:
        T1 = Tree(newick1, format=1)
        T2 = Tree(newick2, format=1)

        # --- Candidate k's ---
        k_values = {
            '$k=2$': 2,
            '$k=3$': 3,
            '$k=\\lfloor\\sqrt{N_{cl}}\\rfloor$': int(round(math.sqrt(N_cl))),
            '$k=\\lfloor\\frac{N_{cl} + 2}{2}\\rfloor$': int(round((2 + N_cl) / 2)),
            '$k=N_{cl}-1$': N_cl - 1,
            '$k=N_{cl}$': N_cl,
        }
        # Keep only valid k in [2, N_cl]
        k_values = {lab: k for lab, k in k_values.items() if 2 <= k <= N_cl}

        # Evaluate in an order that tends to tighten a BSD bound quickly.
        sqrtk = int(round(math.sqrt(N_cl)))

        def sort_key(item):
            lab, kval = item
            if kval == sqrtk:
                return (0, 0, 0)
            if kval in (2, 3):
                return (1, kval, 0)
            return (2, abs(kval - sqrtk), kval)

        ordered_k_items = sorted(k_values.items(), key=sort_key)

        min_bsd = float('inf')
        optimal_k_bsd = None
        per_k_bsd = {}

        for k_label, k in ordered_k_items:
            cap = None if STRICT_IDENTICAL_OUTPUT else min_bsd
            bsd_distance = compute_bsd_kncl(T1, T2, k, cap=cap)

            # record BSD
            if bsd_distance is not None:
                per_k_bsd[k_label] = bsd_distance
                if bsd_distance < min_bsd:
                    min_bsd = bsd_distance
                    optimal_k_bsd = k

        if not per_k_bsd:
            return None

        overlap_bin = round(overlap, 1)
        return (overlap_bin, per_k_bsd, optimal_k_bsd, overlap)

    except Exception:
        return None


# Core computation (parallel)

def compute_and_save_distances(dataset, dataset_name, save_file, n_jobs=None):
    """
    Compute BSD(k–NCL) distances for a single dataset and save to file.
    Stores BSD(k–NCL) results only.

    n_jobs:
      - None or 0: auto (use os.cpu_count())
      - >0: run with that many worker processes
      - <0: relative to CPU count (e.g., -1 => use all but one core)
    """
    print(f"\nProcessing dataset: {dataset_name}")
    # {overlap_bin: {k_label: [values]}}
    bsd_results = {}
    # [(overlap, optimal_k)]
    optimal_k_data_bsd = []

    num_trees = len(dataset)

    # Precompute leaf-name sets once per tree to avoid repeated get_leaf_names() calls
    leaf_sets = [set(t.get_leaf_names()) for t in dataset]
    
    # Precompute Newick strings
    newicks = [t.write(format=1) for t in dataset]

    # Prepare tasks after applying pair-level filters
    tasks = []
    for i in range(num_trees):
        leaves1 = leaf_sets[i]
        for j in range(i + 1, num_trees):
            leaves2 = leaf_sets[j]
            inter = len(leaves1 & leaves2)
            if inter < 3:
                continue
            union = len(leaves1 | leaves2)
            overlap = inter / union
            if overlap < 0.1 or overlap > 0.9:
                continue
            N_cl = inter
            tasks.append((newicks[i], newicks[j], N_cl, overlap))

    if not tasks:
        print(f"No qualifying pairs in {dataset_name}.")
        with open(save_file, 'wb') as f:
            pickle.dump((bsd_results, optimal_k_data_bsd), f)
        return bsd_results, optimal_k_data_bsd

    n_workers = _effective_workers(n_jobs)
    if n_workers == 1:
        print("Running in single-process mode.")
        for idx, args in enumerate(tasks, 1):
            out = _evaluate_pair_worker(args)
            if out is None:
                continue
            overlap_bin, per_k_bsd, opt_k_bsd, overlap = out

            for k_label, dist in per_k_bsd.items():
                bsd_results.setdefault(overlap_bin, {}).setdefault(k_label, []).append(dist)

            if opt_k_bsd is not None:
                optimal_k_data_bsd.append((overlap, opt_k_bsd))

            if idx % 100 == 0:
                print(f"Processed {idx}/{len(tasks)} pairs in {dataset_name}")
    else:
        print(f"Running with {n_workers} worker processes.")
        try:
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futures = [ex.submit(_evaluate_pair_worker, args) for args in tasks]
                for idx, fut in enumerate(as_completed(futures), 1):
                    try:
                        out = fut.result()
                    except Exception as e:
                        # Log and skip failed tasks; don't crash the whole run
                        print(f"Worker failed on a pair: {repr(e)}")
                        continue
                    if out is None:
                        continue

                    overlap_bin, per_k_bsd, opt_k_bsd, overlap = out

                    for k_label, dist in per_k_bsd.items():
                        bsd_results.setdefault(overlap_bin, {}).setdefault(k_label, []).append(dist)

                    if opt_k_bsd is not None:
                        optimal_k_data_bsd.append((overlap, opt_k_bsd))

                    if idx % 100 == 0:
                        print(f"Processed {idx}/{len(tasks)} pairs in {dataset_name}")
        except BrokenProcessPool:
            print("Process pool crashed — falling back to single-process mode.")
            for idx, args in enumerate(tasks, 1):
                out = _evaluate_pair_worker(args)
                if out is None:
                    continue
                overlap_bin, per_k_bsd, opt_k_bsd, overlap = out

                for k_label, dist in per_k_bsd.items():
                    bsd_results.setdefault(overlap_bin, {}).setdefault(k_label, []).append(dist)

                if opt_k_bsd is not None:
                    optimal_k_data_bsd.append((overlap, opt_k_bsd))

                if idx % 100 == 0:
                    print(f"Processed {idx}/{len(tasks)} pairs in {dataset_name}")

    # Save BSD(k–NCL) results for this dataset
    with open(save_file, 'wb') as f:
        pickle.dump((bsd_results, optimal_k_data_bsd), f)
    print(f"BSD(k–NCL) results for {dataset_name} saved to {save_file}")

    return bsd_results, optimal_k_data_bsd


# Plotting (BSD(k–NCL))

def plot_bsd_vs_k_per_dataset(all_bsd_results, dataset_names):
    plt = _setup_matplotlib()
    # Use the exact same keys as the compute step
    k_labels_order = [
        '$k=2$', '$k=3$', '$k=\\lfloor\\sqrt{N_{cl}}\\rfloor$',
        '$k=\\lfloor\\frac{N_{cl} + 2}{2}\\rfloor$', '$k=N_{cl}-1$', '$k=N_{cl}$'
    ]
    # Display labels
    label_map = {
        '$k=2$': r'$2$',
        '$k=3$': r'$3$',
        '$k=\\lfloor\\sqrt{N_{cl}}\\rfloor$': r'$\lfloor\sqrt{N_{cl}}\rfloor$',
        '$k=\\lfloor\\frac{N_{cl} + 2}{2}\\rfloor$': r'$\lfloor\frac{N_{cl} + 2}{2}\rfloor$',
        '$k=N_{cl}-1$': r'$N_{cl}-1$',
        '$k=N_{cl}$': r'$N_{cl}$'
    }

    fig, axs = plt.subplots(2, 2, figsize=(15, 15))
    axs = axs.flatten()
    for idx, results in enumerate(all_bsd_results):
        ax = axs[idx]
        for overlap_bin in sorted(results.keys()):
            bsd_means = []
            label_used = False
            for k_label in k_labels_order:
                bsd_values = results.get(overlap_bin, {}).get(k_label, [])
                if bsd_values:
                    label_used = True
                bsd_mean = np.mean(bsd_values) if bsd_values else np.nan
                bsd_means.append(bsd_mean)

            if label_used:
                ax.plot(k_labels_order, bsd_means, marker='o', label=f'{overlap_bin}')

        ax.set_xlabel(r'$k$', fontsize=14)
        ax.set_ylabel('Average BSD($k$–NCL)', fontsize=14)
        ax.set_title(dataset_names[idx], loc='left', fontsize=16)
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xticks(range(len(k_labels_order)))
        ax.set_xticklabels([label_map[k] for k in k_labels_order], fontsize=14)
        ax.grid(True)

    # Shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title=r'Overlap level $p$',
                   loc='lower center', bbox_to_anchor=(0.5, -0.035),
                   ncol=len(labels), fontsize=12, title_fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('bsd_vs_k_per_dataset.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('bsd_vs_k_per_dataset.svg', format='svg', bbox_inches='tight')
    #plt.show()

# Optional
def plot_optimal_k_vs_overlap_per_dataset(all_optimal_k_results, dataset_names):
    plt = _setup_matplotlib()
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))
    axs = axs.flatten()
    for idx, optimal_k_data in enumerate(all_optimal_k_results):
        overlaps = [data[0] for data in optimal_k_data]
        ks = [data[1] for data in optimal_k_data]
        ax = axs[idx]
        ax.scatter(overlaps, ks, alpha=0.6)
        ax.set_xlabel('Overlap level $p$')
        ax.set_ylabel('Optimal $k$ (BSD)')
        ax.set_title(f'{dataset_names[idx]}')
        ax.set_xlim(0.1, 0.9)
        ax.set_ylim(1, None)
        ax.grid(True)
    plt.tight_layout()
    plt.savefig('optimal_k_vs_overlap_per_dataset_bsd.pdf', format='pdf')
    plt.savefig('optimal_k_vs_overlap_per_dataset_bsd.svg', format='svg')
    plt.show()


# Main

if __name__ == "__main__":
    
    mp.freeze_support()
    if os.name == "nt":
        try:
            mp.set_start_method("spawn")
        except RuntimeError:
            # already set by the parent process
            pass

    datasets = load_datasets(['amphibians100.txt', 'birds100.txt', 'mammals100.txt', 'sharks100.txt'])
    dataset_names = ["(a) Amphibians", "(b) Birds", "(c) Mammals", "(d) Sharks"]

    # Resolve requested parallelism with Windows/Jupyter safeguards
    n_workers = _effective_workers(PARALLEL_JOBS)

    all_bsd_results = []
    all_optimal_k_bsd_results = []

    for i, dataset_name in enumerate(dataset_names):
        # Make a filesystem-safe cache filename
        save_file = f"{dataset_name.replace(' ', '_').replace('(', '').replace(')', '').replace('.', '').lower()}_results.pkl"
        if os.path.exists(save_file):
            # Load previously saved results
            with open(save_file, 'rb') as f:
                loaded = pickle.load(f)

            if isinstance(loaded, (list, tuple)) and len(loaded) == 2:
                bsd_results, optimal_k_bsd = loaded
                print(f"Loaded cached BSD(k–NCL) results for {dataset_name}")
            elif isinstance(loaded, (list, tuple)) and len(loaded) == 4:
                bsd_results, optimal_k_bsd, *_ = loaded
                print(f"Loaded cached BSD+ results for {dataset_name}")
            else:
                print(f"Unrecognized cache format for {dataset_name}; recomputing.")
                bsd_results, optimal_k_bsd = compute_and_save_distances(
                    datasets[i], dataset_name, save_file, n_jobs=n_workers
                )
        else:
            # Process and save results (parallel where possible)
            bsd_results, optimal_k_bsd = compute_and_save_distances(
                datasets[i], dataset_name, save_file, n_jobs=n_workers
            )

        all_bsd_results.append(bsd_results)
        all_optimal_k_bsd_results.append(optimal_k_bsd)

    # Combine all datasets (optional aggregate analysis if needed later)
    combined_bsd_results = {}
    combined_optimal_k_bsd_results = []
    for bsd_result in all_bsd_results:
        for overlap_bin, k_data in bsd_result.items():
            if overlap_bin not in combined_bsd_results:
                combined_bsd_results[overlap_bin] = {}
            for k_label, values in k_data.items():
                combined_bsd_results[overlap_bin].setdefault(k_label, []).extend(values)
    for optimal_k_data in all_optimal_k_bsd_results:
        combined_optimal_k_bsd_results.extend(optimal_k_data)

    print("All datasets combined. Proceeding to plots...")

    # Visualizations: BSD(k–NCL)
    plot_bsd_vs_k_per_dataset(all_bsd_results, dataset_names)
    #plot_optimal_k_vs_overlap_per_dataset(all_optimal_k_bsd_results, dataset_names) #optional