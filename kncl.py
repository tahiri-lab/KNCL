#!/usr/bin/env python3
"""
This script implements the k-NCL algorithm for completing phylogenetic trees
that are defined on different but overlapping taxon sets.

Author: Aleksandr Koshkarov

Usage:
    python kncl.py -i input.newick k -o output.txt

Requires:
    - ete3 (pip install ete3)

Implementation note (performance):    
    This implementation prioritizes clarity and determinism, so it may run
    somewhat slower in practice than the theoretical analysis suggests.
"""

import argparse
from functools import lru_cache
from math import floor, log2
from typing import Dict, Iterable, List, Sequence, Tuple
from ete3 import Tree

# ──────────────────────────────────────────────────────────────────────
#  Helper functions
# ──────────────────────────────────────────────────────────────────────
def _refresh_name_cache(tree: Tree):
    tree._nmap = {n.name: n for n in tree.traverse() if n.name}


def ensure_unique_internal_node_names(tree: Tree, prefix: str):
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and (not node.name or node.name.startswith(prefix)):
            node.name = f"{prefix}_{counter}"
            counter += 1
    _refresh_name_cache(tree)


def _is_descendant(node: Tree, ancestor: Tree) -> bool:
    """True iff 'ancestor' lies on the path from 'node' up to the root."""
    return node == ancestor or ancestor in node.iter_ancestors()


# ──────────────────────────────────────────────────────────────────────
#  File IO helpers
# ──────────────────────────────────────────────────────────────────────
def parse_input_file(filename: str) -> List[Tree]:
    with open(filename, "r") as f:
        return [Tree(line.strip(), format=1) for line in f if line.strip()]


def write_output_file(filename: str, results: Sequence[str]):
    with open(filename, "w") as f:
        for entry in results:
            f.write(entry + "\n\n")


# ──────────────────────────────────────────────────────────────────────
#  Find common / distinct leaves and maximal distinct-leaf subtrees
# ──────────────────────────────────────────────────────────────────────
def find_common_leaves(T1: Tree, T2: Tree) -> List[str]:
    return list(set(T1.get_leaf_names()).intersection(T2.get_leaf_names()))


def find_distinct_leaves(T: Tree, common: Iterable[str]) -> List[str]:
    commons = set(common)
    return [n for n in T.get_leaf_names() if n not in commons]


def findSD(tree: Tree, dl_set: Iterable[str]) -> List[Tree]:
    dl = set(dl_set)
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            n.allDistinct = n.name in dl
        else:
            n.allDistinct = all(ch.allDistinct for ch in n.children)
    return [n for n in tree.traverse("postorder")
            if n.allDistinct and (n.up is None or not n.up.allDistinct)]


# ──────────────────────────────────────────────────────────────────────
#  LCA / distance oracle (Euler tour + RMQ)
# ──────────────────────────────────────────────────────────────────────
class DistOracle:
    __slots__ = ("dist_to_root", "_euler", "_depth", "_first", "_st", "_log2", "_dist_cached")

    def __init__(self, tree: Tree, cache_size: int = 20_000):
        self.dist_to_root: Dict[Tree, float] = {}
        self._annotate_depth(tree)
        self._build_lca_struct(tree)

        @lru_cache(maxsize=cache_size)
        def _dist_cached(a: Tree, b: Tree) -> float:
            lca = self._lca(a, b)
            return self.dist_to_root[a] + self.dist_to_root[b] - 2 * self.dist_to_root[lca]

        self._dist_cached = _dist_cached

    def dist(self, a: Tree, b: Tree) -> float:
        return self._dist_cached(a, b)

    def _annotate_depth(self, tree: Tree):
        for n in tree.traverse("preorder"):
            self.dist_to_root[n] = 0.0 if n.up is None else self.dist_to_root[n.up] + n.dist

    def _build_lca_struct(self, tree: Tree):
        euler, depth, first = [], [], {}

        def dfs(v: Tree, d: int):
            first.setdefault(v, len(euler))
            euler.append(v)
            depth.append(d)
            for ch in v.children:
                dfs(ch, d + 1)
                euler.append(v)
                depth.append(d)

        dfs(tree, 0)
        self._euler, self._depth, self._first = euler, depth, first

        m = len(euler)
        kmax = floor(log2(m)) + 1
        st = [[0] * m for _ in range(kmax)]
        st[0] = list(range(m))
        for k in range(1, kmax):
            half = 1 << (k - 1)
            for i in range(m - (1 << k) + 1):
                a, b = st[k - 1][i], st[k - 1][i + half]
                st[k][i] = a if depth[a] < depth[b] else b
        self._st = st
        self._log2 = [0] * (m + 1)
        for i in range(2, m + 1):
            self._log2[i] = self._log2[i >> 1] + 1

    def _lca(self, a: Tree, b: Tree) -> Tree:
        ia, ib = self._first[a], self._first[b]
        if ia > ib:
            ia, ib = ib, ia
        span = ib - ia + 1
        k = self._log2[span]
        left = self._st[k][ia]
        right = self._st[k][ib - (1 << k) + 1]
        return self._euler[left] if self._depth[left] < self._depth[right] else self._euler[right]

    def dist_leaf_to_node(self, leaf: Tree, node: Tree) -> float:
        return self.dist(leaf, node)


# ──────────────────────────────────────────────────────────────────────
#  Adjustment-rate functions
# ──────────────────────────────────────────────────────────────────────
def compute_global_adjustment_rate(T1: Tree, T2: Tree, common_names: Sequence[str],
                                   D1: DistOracle, D2: DistOracle) -> float:
    L1 = [T1 & n for n in common_names]
    L2 = [T2 & n for n in common_names]
    num = den = 0.0
    for i in range(len(L1)):
        for j in range(i + 1, len(L1)):
            num += D1.dist(L1[i], L1[j])
            den += D2.dist(L2[i], L2[j])
    return 1.0 if abs(den) < 1e-15 else num / den


def compute_leaf_based_adjustment_rate(Dtgt: DistOracle, Dsrc: DistOracle,
                                       leaf_tgt: Tree, leaf_src: Tree,
                                       common_names: Sequence[str], tree_tgt: Tree, tree_src: Tree) -> float:
    num = den = 0.0
    for name in common_names:
        # Including the identical leaf adds zero; we skip it
        if name == leaf_tgt.name:
            continue
        num += Dtgt.dist(leaf_tgt, tree_tgt & name)
        den += Dsrc.dist(leaf_src, tree_src & name)
    return 1.0 if abs(den) < 1e-15 else num / den


# ──────────────────────────────────────────────────────────────────────
#  Ordering helpers: idx(l) for common leaves, idx(e) for original edges
# ──────────────────────────────────────────────────────────────────────
def common_order_and_index(T1: Tree, T2: Tree) -> Tuple[List[str], Dict[str, int]]:
    common = find_common_leaves(T1, T2)
    ordered = sorted(common)  # ascending lexicographic order of labels
    idx = {name: i for i, name in enumerate(ordered)}
    return ordered, idx


def assign_edge_ranks(tree: Tree, feature_name: str = "edge_rank") -> None:
    """Assign a fixed depth-first order rank to each original branch.
    The rank is stored on the child node as 'feature_name'.
    """
    rank = 1
    for node in tree.traverse("preorder"):
        for ch in node.children:
            ch.add_feature(feature_name, rank)
            rank += 1


# ──────────────────────────────────────────────────────────────────────
#  Misc helpers
# ──────────────────────────────────────────────────────────────────────
def scale_subtree(root: Tree, factor: float) -> None:
    for n in root.traverse():
        n.dist *= factor


def get_k_nearest_common_leaves(Dsrc: DistOracle, root: Tree,
                                common_names: Sequence[str], tree_src: Tree,
                                k: int, idx_map: Dict[str, int]) -> List[str]:
    pairs = [(name, Dsrc.dist_leaf_to_node(tree_src & name, root)) for name in common_names]
    # Lexicographic ordering: first by distance, then by fixed idx(l)
    pairs.sort(key=lambda t: (t[1], idx_map[t[0]]))
    return [name for name, _ in pairs[:min(k, len(pairs))]]


# ──────────────────────────────────────────────────────────────────────
#  Objective-function optimization
# ──────────────────────────────────────────────────────────────────────
def find_optimal_insertion_point(target_tree: Tree, Dtgt: DistOracle, subtree_root: Tree,
                                 ncl_names: Sequence[str], d_p: Dict[str, float],
                                 dl_names: Iterable[str], tol: float = 1e-12) -> Tuple[Tuple[Tree, Tree], float, float]:
    """Return ((parent, child), x_opt, f_min).
    Tie-breaking rules:
       1) smallest depth of candidate point v_e(x*)
       2) smallest original branch rank idx(e)
    """
    best_edge: Tuple[Tree, Tree] = None
    best_x: float = 0.0
    best_val: float = float("inf")
    best_depth: float = float("inf")
    best_rank: int = 10**18

    dl_set = set(dl_names)
    m = len(ncl_names)
    if m == 0:
        return None, 0.0, float("inf")

    root = target_tree

    for child in target_tree.traverse("postorder"):
        parent = child.up
        if parent is None:
            continue

        # Skip edges fully inside an already inserted subtree
        if getattr(child, "inInserted", False):
            continue

        # Avoid placing inside a branch that still contains target-tree distinct leaves
        if set(child.get_leaf_names()) & dl_set:
            continue

        e_len = child.dist
        if e_len < 1e-15:
            continue

        # The quadratic objective f_e(x) = sum_c (A_c + eps_c * x * e_len - d_p[c])^2
        # The minimizer (unconstrained) is:
        #   x* = (1 / (e_len * m)) * sum_c eps_c * (d_p[c] - A_c)
        numerator = 0.0
        for name in ncl_names:
            leaf = target_tree & name
            A_c = Dtgt.dist_leaf_to_node(leaf, parent)
            eps = -1 if _is_descendant(leaf, child) else 1
            numerator += eps * (d_p[name] - A_c)
        x_opt = numerator / (e_len * m)

        # Constrain x to [0, 1) (or [0, 1 - eps] for leaves to avoid zero-length child)
        if child.is_leaf():
            # keep a small terminal branch
            x_opt = max(0.0, min(x_opt, 1.0 - 1e-3 / e_len))
        else:
            x_opt = max(0.0, min(x_opt, 1.0))

        # Evaluate the objective at x*= and at boundary x=0 (node attach)
        for x in (x_opt, 0.0):
            of_val = 0.0
            for name in ncl_names:
                leaf = target_tree & name
                A_c = Dtgt.dist_leaf_to_node(leaf, parent)
                eps = -1 if _is_descendant(leaf, child) else 1
                obs = A_c + eps * x * e_len
                diff = obs - d_p[name]
                of_val += diff * diff

            # Candidate point depth from the root (for tie-breaking rule 1)
            cand_depth = Dtgt.dist_to_root[parent] + x * e_len

            # Edge rank (on the original branch) — stored on the child node and inherited on splits
            edge_rank = getattr(child, "edge_rank", 10**18)

            better = False
            if of_val < best_val - tol:
                better = True
            elif abs(of_val - best_val) <= tol:
                if cand_depth < best_depth - tol:
                    better = True
                elif abs(cand_depth - best_depth) <= tol and edge_rank < best_rank:
                    better = True

            if better:
                best_val, best_edge, best_x = of_val, (parent, child), x
                best_depth, best_rank = cand_depth, edge_rank

    return best_edge, best_x, best_val


def insert_subtree_at_point(target_tree: Tree, edge: Tuple[Tree, Tree], x_opt: float, subtree_copy: Tree):
    parent, child = edge
    orig_len = child.dist
    d_up, d_down = x_opt * orig_len, (1 - x_opt) * orig_len

    if x_opt < 1e-15:
        # Attach directly at the existing node (no split)
        parent.add_child(subtree_copy)
    else:
        # Split the original branch into parent->mid and mid->child
        mid = Tree()
        mid.dist = d_up

        # Inherit the original branch rank idx(e) for both new sub-branches
        if hasattr(child, "edge_rank"):
            mid.add_feature("edge_rank", child.edge_rank)

        parent.remove_child(child)
        parent.add_child(mid)
        child.dist = d_down
        mid.add_child(child)
        mid.add_child(subtree_copy)

    _refresh_name_cache(target_tree)


# ──────────────────────────────────────────────────────────────────────
#  Main k-NCL
# ──────────────────────────────────────────────────────────────────────
def kNCL(T1_raw: Tree, T2_raw: Tree, k: int = None) -> Tuple[Tree, Tree]:
    T1, T2 = T1_raw.copy(), T2_raw.copy()
    ensure_unique_internal_node_names(T1, "T1IN")
    ensure_unique_internal_node_names(T2, "T2IN")

    # Fixed common-leaf ordering and idx(l) map
    common_ordered, idx_map = common_order_and_index(T1, T2)
    if len(common_ordered) < 2:
        raise ValueError("Need at least two common leaves")

    # Enforce 2 ≤ k ≤ |CL|
    if k is None:
        k = (len(common_ordered) + 2) // 2   # default value
    if k < 2:
        raise ValueError("k must be ≥ 2")
    k = min(k, len(common_ordered))

    # Distinct leaves and maximal distinct-leaf subtrees
    DL1, DL2 = find_distinct_leaves(T1, common_ordered), find_distinct_leaves(T2, common_ordered)
    SD1, SD2 = findSD(T1, DL1), findSD(T2, DL2)

    # Assign fixed original-branch ranks idx(e) for each tree (stored on child nodes)
    assign_edge_ranks(T1)
    assign_edge_ranks(T2)

    # Distance oracles
    D1, D2 = DistOracle(T1), DistOracle(T2)

    # Global adjustment rates
    r12 = compute_global_adjustment_rate(T1, T2, common_ordered, D1, D2)
    r21 = compute_global_adjustment_rate(T2, T1, common_ordered, D2, D1)

    def _insert_all(target_tree: Tree, source_tree: Tree, subtrees: Sequence[Tree],
                    global_r: float, Dtgt: DistOracle, Dsrc: DistOracle,
                    idx_map_local: Dict[str, int]):
        for S in subtrees:
            # Copy and globally scale subtree (including its root edge length)
            S_copy = S.copy()
            scale_subtree(S_copy, global_r)
            for n in S_copy.traverse():
                n.add_feature("inInserted", True)

            # k nearest common leaves in the source tree
            ncl = get_k_nearest_common_leaves(Dsrc, S, common_ordered, source_tree, k, idx_map_local)
            if not ncl:
                continue

            # Attachment node in the source tree
            att_node_src = S.up if S.up else S

            # Position distances in the target tree
            d_p = {}
            for name in ncl:
                leaf_src = source_tree & name
                leaf_tgt = target_tree & name
                r_lc = compute_leaf_based_adjustment_rate(Dtgt, Dsrc, leaf_tgt, leaf_src,
                                                          common_ordered, target_tree, source_tree)
                d_p[name] = Dsrc.dist_leaf_to_node(leaf_src, att_node_src) * r_lc

            # Find optimal insertion point in the target tree (with deterministic tie-breaking)
            edge, x_opt, _ = find_optimal_insertion_point(
                target_tree, Dtgt, S, ncl, d_p, find_distinct_leaves(target_tree, common_ordered)
            )

            if edge:
                insert_subtree_at_point(target_tree, edge, x_opt, S_copy)
                # Rebuild the oracle (structure and depths changed)
                Dtgt.__init__(target_tree)
                _refresh_name_cache(target_tree)

    _insert_all(T1, T2, SD2, r12, D1, D2, idx_map)
    _insert_all(T2, T1, SD1, r21, D2, D1, idx_map)

    # Clear internal node names for output cleanliness
    for t in (T1, T2):
        for n in t.traverse():
            if not n.is_leaf():
                n.name = ""

    return T1, T2


# ──────────────────────────────────────────────────────────────────────
#  Main CLI logic
# ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Run k-NCL algorithm on Newick trees." )
    parser.add_argument("-i", "--input", required=True, help="Input file with Newick trees (one per line)")
    parser.add_argument("k", type=int, help="k value (≥ 2 and ≤ |CL|) for the k-NCL algorithm. If not specified, the default value is ⌊(|CL|+2)/2⌋ per tree pair.")
    parser.add_argument("-o", "--output", required=True, help="Output file for results")
    args = parser.parse_args()

    trees = parse_input_file(args.input)
    if len(trees) < 2:
        raise ValueError("Input must contain at least two trees.")

    results = []
    for i in range(len(trees)):
        for j in range(i + 1, len(trees)):
            try:
                T1_plus, T2_plus = kNCL(trees[i], trees[j], args.k)
                res = (
                    f"Tree pair {i+1} & {j+1}:\n"
                    f"Completed Tree 1:\n{T1_plus.write(format=1).strip()}\n"
                    f"Completed Tree 2:\n{T2_plus.write(format=1).strip()}"
                )
            except Exception as e:
                res = f"Tree pair {i+1} & {j+1}: cannot complete (error: {e})"
            results.append(res)

    write_output_file(args.output, results)


if __name__ == "__main__":
    main()
