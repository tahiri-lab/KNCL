#!/usr/bin/env python
# coding: utf-8

# ## Cluster separation analysis by method
# Tested: RF(+), RF(k-NCL), BSD(k-NCL)


import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class RenameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith("numpy._core"):
            module = module.replace("numpy._core", "numpy.core")
        return super().find_class(module, name)

def renamed_load(file_obj):
    return RenameUnpickler(file_obj).load()

PKL_PATH = Path("distance_matrices_clusters.pkl")
with open(PKL_PATH, "rb") as fh:
    matrices = renamed_load(fh)

metrics = list(matrices.keys())
datasets = list(next(iter(matrices.values())).keys())  # ['Amphibians','Birds','Mammals','Sharks']
cluster_size = 5

# Helpers

def within_between(mat: np.ndarray, cluster_size: int = 5):
    n = mat.shape[0]
    idx = np.arange(n)
    clusters = [idx[i:i+cluster_size] for i in range(0, n, cluster_size)]
    within = []
    for cl in clusters:
        for a in range(len(cl)):
            for b in range(a+1, len(cl)):
                within.append(mat[cl[a], cl[b]])
    between = []
    for i in range(len(clusters)):
        for j in range(i+1, len(clusters)):
            for a in clusters[i]:
                for b in clusters[j]:
                    between.append(mat[a, b])
    return np.asarray(within, float), np.asarray(between, float)

def auc_between_gt_within(within: np.ndarray, between: np.ndarray):
    # AUROC-like probability P(between > within), ties count 0.5
    w = within.reshape(-1, 1)
    b = between.reshape(1, -1)
    gt = (b > w).sum()
    eq = (b == w).sum()
    total = w.size * b.size
    return (gt + 0.5 * eq) / total

def cliffs_delta_from_auc(auc: float):
    # delta = P(b>w) - P(b<w) = 2*AUC - 1
    return 2*auc - 1

def silhouette_score_from_dist(mat: np.ndarray, cluster_size: int = 5):
    # mean silhouette using the distance matrix directly
    n = mat.shape[0]
    labels = np.repeat(np.arange(n // cluster_size), cluster_size)
    K = n // cluster_size
    clusters = [np.where(labels == k)[0] for k in range(K)]
    a = np.zeros(n, float)
    b = np.zeros(n, float)
    for i in range(n):
        ci = labels[i]
        same = clusters[ci]
        a[i] = mat[i, same[same != i]].mean() if len(same) > 1 else 0.0
        means = [mat[i, clusters[cj]].mean() for cj in range(K) if cj != ci]
        b[i] = min(means) if means else 0.0
    s = (b - a) / np.maximum(a, b)
    s = np.where(np.isfinite(s), s, 0.0)
    return float(np.mean(s))

def dunn_index(mat: np.ndarray, cluster_size: int = 5):
    n = mat.shape[0]
    labels = np.repeat(np.arange(n // cluster_size), cluster_size)
    K = n // cluster_size
    clusters = [np.where(labels == k)[0] for k in range(K)]
    # max intra-cluster distance (diameters)
    max_diam = 0.0
    for cl in clusters:
        sub = mat[np.ix_(cl, cl)]
        iu = np.triu_indices_from(sub, k=1)
        if iu[0].size:
            max_diam = max(max_diam, float(sub[iu].max()))
    # min inter-cluster distance
    min_inter = np.inf
    for i in range(K):
        for j in range(i+1, K):
            sub = mat[np.ix_(clusters[i], clusters[j])]
            min_inter = min(min_inter, float(sub.min()))
    return (min_inter / max_diam) if max_diam > 0 else (np.inf if np.isfinite(min_inter) else 0.0)

# Compute all stats

rows = []
for metric in metrics:
    for dataset in datasets:
        mat = matrices[metric][dataset]
        within, between = within_between(mat, cluster_size=cluster_size)
        auc = auc_between_gt_within(within, between)
        delta = cliffs_delta_from_auc(auc)
        med_w = float(np.median(within))
        med_b = float(np.median(between))
        median_ratio = (med_b / med_w) if med_w != 0 else np.inf
        sil = silhouette_score_from_dist(mat, cluster_size=cluster_size)
        dunn = dunn_index(mat, cluster_size=cluster_size)
        rows.append({
            "Metric": metric,
            "Dataset": dataset,
            "Within_N": within.size,
            "Between_N": between.size,
            "Median_within": med_w,
            "Median_between": med_b,
            "Median_ratio_btw_over_within": median_ratio,
            "AUROC_P(between>within)": auc,
            "Cliffs_delta": delta,
            "Mean_silhouette": sil,
            "Dunn_index": dunn,
        })

df = pd.DataFrame(rows)

# Per-dataset ranks (higher is better)
rank_cols = ["AUROC_P(between>within)", "Cliffs_delta",
             "Median_ratio_btw_over_within", "Mean_silhouette", "Dunn_index"]
for col in rank_cols:
    df[f"Rank_{col}"] = df.groupby("Dataset")[col].rank(ascending=False, method="average")
df["Mean_rank_across_metrics"] = df[[f"Rank_{c}" for c in rank_cols]].mean(axis=1)

# Overall (avg across datasets)
overall = (
    df.groupby("Metric")[["AUROC_P(between>within)", "Cliffs_delta",
                          "Median_ratio_btw_over_within", "Mean_silhouette",
                          "Dunn_index", "Mean_rank_across_metrics"]]
      .mean()
      .sort_values("Mean_rank_across_metrics", ascending=True)
      .reset_index()
)

# Save outputs

df_path = Path("cluster_separation_metrics_by_dataset.csv")
overall_path = Path("cluster_separation_overall.csv")
df.round(6).to_csv(df_path, index=False)
overall.round(6).to_csv(overall_path, index=False)

# Overall silhouette bar
fig = plt.figure()
ax = fig.add_subplot(111)
ax.bar(overall["Metric"], overall["Mean_silhouette"])
ax.set_ylabel("Mean silhouette")
ax.set_title("Overall cluster separation by method (mean silhouette across species)")
fig.savefig("overall_silhouette_bar.png", bbox_inches="tight")