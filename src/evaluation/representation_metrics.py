"""Metrics that can be computed without biological cell-type labels."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

CELLTYPE_UNAVAILABLE_REASON = "Ground-truth cell-type labels unavailable"


def celltype_metrics_unavailable() -> dict[str, Any]:
    """Return NA cell-type metrics. Do not fill zeros."""
    return {
        "ARI": np.nan,
        "NMI": np.nan,
        "celltype_silhouette": np.nan,
        "celltype_asw": np.nan,
        "celltype_metrics_reason": CELLTYPE_UNAVAILABLE_REASON,
    }


def batch_silhouette(latent: np.ndarray, batch_labels: pd.Series | np.ndarray) -> float:
    """Silhouette of the latent space with respect to batch labels.

    This is a batch-structure diagnostic, not a claim that mixing is better.
    """
    labels = np.asarray(batch_labels).astype(str)
    if latent.ndim != 2:
        raise ValueError(f"Expected 2D latent matrix, got shape {latent.shape}")
    if latent.shape[0] != len(labels):
        raise ValueError("Latent rows and batch labels have different lengths.")
    if len(np.unique(labels)) < 2:
        raise ValueError("Batch silhouette requires at least two batch categories.")
    return float(silhouette_score(latent, labels, metric="euclidean"))


def latent_diagnostics(latent: np.ndarray) -> dict[str, Any]:
    if not np.isfinite(latent).all():
        n_nan = int(np.isnan(latent).sum())
        n_inf = int(np.isinf(latent).sum())
        raise ValueError(f"Latent representation contains NaN ({n_nan}) or Inf ({n_inf}).")
    means = latent.mean(axis=0)
    sds = latent.std(axis=0, ddof=1)
    return {
        "n_cells": int(latent.shape[0]),
        "n_latent": int(latent.shape[1]),
        "n_nan": 0,
        "n_inf": 0,
        "dim_mean": means.tolist(),
        "dim_sd": sds.tolist(),
        "global_mean": float(latent.mean()),
        "global_sd": float(latent.std(ddof=1)),
    }


def exploratory_latent_comparison(
    latent_a: np.ndarray,
    latent_b: np.ndarray,
    seed: int = 0,
    n_distance_subset: int = 2000,
    knn_k: int = 15,
) -> dict[str, Any]:
    """Exploratory scVI vs totalVI geometry. Not a superiority test."""
    if latent_a.shape != latent_b.shape:
        raise ValueError(
            f"Latent shapes differ: {latent_a.shape} vs {latent_b.shape}."
        )
    n_cells = latent_a.shape[0]
    from scipy.spatial import procrustes
    from scipy.spatial.distance import pdist
    from scipy.stats import pearsonr
    from sklearn.neighbors import NearestNeighbors

    _, _, disparity = procrustes(latent_a, latent_b)

    nn_a = NearestNeighbors(n_neighbors=knn_k + 1, metric="euclidean").fit(latent_a)
    nn_b = NearestNeighbors(n_neighbors=knn_k + 1, metric="euclidean").fit(latent_b)
    idx_a = nn_a.kneighbors(return_distance=False)[:, 1:]
    idx_b = nn_b.kneighbors(return_distance=False)[:, 1:]
    overlaps = [
        len(set(row_a) & set(row_b)) / float(knn_k) for row_a, row_b in zip(idx_a, idx_b)
    ]

    rng = np.random.default_rng(seed)
    subset_n = min(n_distance_subset, n_cells)
    subset = rng.choice(n_cells, size=subset_n, replace=False)
    corr, p_value = pearsonr(pdist(latent_a[subset]), pdist(latent_b[subset]))
    return {
        "label": "exploratory representation comparison",
        "not_a_superiority_test": True,
        "n_cells": n_cells,
        "n_latent": int(latent_a.shape[1]),
        "procrustes_disparity": float(disparity),
        "knn_k": knn_k,
        "mean_knn_overlap": float(np.mean(overlaps)),
        "median_knn_overlap": float(np.median(overlaps)),
        "pairwise_distance_subset_n": subset_n,
        "pairwise_distance_pearson_r": float(corr),
        "pairwise_distance_pearson_p": float(p_value),
    }


def compute_representation_metrics(
    latent: np.ndarray,
    batch_labels: pd.Series | np.ndarray,
    cell_type_labels: pd.Series | np.ndarray | None = None,
) -> dict[str, Any]:
    """Batch silhouette only. Cell-type metrics are NA until labels exist."""
    metrics = {
        "batch_silhouette": batch_silhouette(latent, batch_labels),
        **celltype_metrics_unavailable(),
        **latent_diagnostics(latent),
    }
    if cell_type_labels is not None:
        unique = pd.Series(cell_type_labels).astype(str).unique().tolist()
        if unique and set(unique) - {"unknown"}:
            metrics["celltype_metrics_reason"] = (
                CELLTYPE_UNAVAILABLE_REASON
                + "; PHASE 3 still withholds cell-type scores until an independent annotation source is added."
            )
    return metrics
