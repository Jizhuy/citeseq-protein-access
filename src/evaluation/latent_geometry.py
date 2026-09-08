"""Exploratory pairwise latent-geometry comparisons. Not a superiority test."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
from sklearn.neighbors import NearestNeighbors


def compare_latent_geometries(
    latents: dict[str, np.ndarray],
    *,
    seed: int = 0,
    n_distance_subset: int = 2000,
    knn_ks: tuple[int, ...] = (15, 30, 50),
) -> pd.DataFrame:
    """Compare every pair with one shared cell subset and one seed."""
    names = list(latents)
    if len(names) < 2:
        raise ValueError("Need at least two latent matrices.")
    n_cells = next(iter(latents.values())).shape[0]
    for name, matrix in latents.items():
        if matrix.ndim != 2 or matrix.shape[0] != n_cells:
            raise ValueError(f"{name} has incompatible shape {matrix.shape}.")
        if not np.isfinite(matrix).all():
            raise ValueError(f"{name} contains NaN or Inf.")

    rng = np.random.default_rng(seed)
    subset_n = min(n_distance_subset, n_cells)
    subset = rng.choice(n_cells, size=subset_n, replace=False)
    subset.sort()

    neighbor_index: dict[str, dict[int, np.ndarray]] = {}
    for name, matrix in latents.items():
        neighbor_index[name] = {}
        for k in knn_ks:
            nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(matrix)
            neighbor_index[name][k] = nn.kneighbors(return_distance=False)[:, 1:]

    rows: list[dict[str, Any]] = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            _, _, disparity = procrustes(latents[a], latents[b])
            d_a = pdist(latents[a][subset])
            d_b = pdist(latents[b][subset])
            corr, p_value = pearsonr(d_a, d_b)
            for k in knn_ks:
                overlaps = [
                    len(set(row_a) & set(row_b)) / float(k)
                    for row_a, row_b in zip(neighbor_index[a][k], neighbor_index[b][k])
                ]
                rows.append(
                    {
                        "representation_a": a,
                        "representation_b": b,
                        "n_cells": n_cells,
                        "n_latent_a": int(latents[a].shape[1]),
                        "n_latent_b": int(latents[b].shape[1]),
                        "seed": seed,
                        "pairwise_distance_subset_n": subset_n,
                        "pairwise_distance_subset_shared": True,
                        "pairwise_distance_pearson_r": float(corr),
                        "pairwise_distance_pearson_p": float(p_value),
                        "knn_k": int(k),
                        "mean_knn_overlap": float(np.mean(overlaps)),
                        "median_knn_overlap": float(np.median(overlaps)),
                        "procrustes_disparity": float(disparity),
                        "label": "exploratory representation comparison",
                        "not_a_superiority_test": True,
                    }
                )
    return pd.DataFrame(rows)
