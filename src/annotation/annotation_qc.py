"""Sanity checks for transferred labels. Does not rewrite labels."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

MARKER_ALIASES = {
    "CD3": ("CD3", "CD3_TotalSeqB"),
    "CD4": ("CD4", "CD4_TotalSeqB"),
    "CD8a": ("CD8a", "CD8", "CD8a_TotalSeqB"),
    "CD14": ("CD14", "CD14_TotalSeqB"),
    "CD16": ("CD16", "CD16_TotalSeqB"),
    "CD19": ("CD19", "CD19_TotalSeqB"),
    "CD56": ("CD56", "CD56_TotalSeqB"),
}


def resolve_protein_column(columns: list[str], marker: str) -> str | None:
    cleaned = {str(c): str(c).replace("_TotalSeqB", "") for c in columns}
    for original, short in cleaned.items():
        if short == marker or original == marker:
            return original
    for alias in MARKER_ALIASES.get(marker, ()):
        if alias in columns:
            return alias
    return None


def protein_marker_validation(
    adata: ad.AnnData,
    labels: pd.Series,
    label_level: str,
    protein_key: str = "protein_expression",
) -> pd.DataFrame:
    protein = adata.obsm[protein_key]
    if not isinstance(protein, pd.DataFrame):
        raise TypeError("protein_expression must be a DataFrame.")
    labels = pd.Series(pd.Series(labels).to_numpy(), index=adata.obs_names).astype(str)
    rows = []
    for marker in MARKER_ALIASES:
        col = resolve_protein_column(list(protein.columns), marker)
        if col is None:
            continue
        values = protein[col].to_numpy()
        for cell_type in sorted(labels.unique()):
            arr = values[labels.to_numpy() == cell_type]
            rows.append(
                {
                    "label_level": label_level,
                    "predicted_cell_type": cell_type,
                    "protein_marker": marker,
                    "protein_column": col,
                    "n_cells": int(len(arr)),
                    "median": float(np.median(arr)),
                    "q25": float(np.quantile(arr, 0.25)),
                    "q75": float(np.quantile(arr, 0.75)),
                    "mean": float(np.mean(arr)),
                    "detection_fraction": float(np.mean(arr > 0)),
                    "source": "original observed protein counts",
                    "note": "sanity check only; labels were not edited",
                }
            )
    return pd.DataFrame(rows)


def celltype_by_batch(labels: pd.Series, batch: pd.Series, label_level: str) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "cell_type": pd.Series(labels).astype(str).to_numpy(),
            "batch": pd.Series(batch).astype(str).to_numpy(),
        }
    )
    counts = frame.value_counts(["cell_type", "batch"]).rename("n_cells").reset_index()
    totals = counts.groupby("batch")["n_cells"].transform("sum")
    counts["proportion_within_batch"] = counts["n_cells"] / totals
    type_tot = counts.groupby("cell_type")["n_cells"].transform("sum")
    counts["proportion_within_type"] = counts["n_cells"] / type_tot
    counts["label_level"] = label_level
    return counts.sort_values(["label_level", "cell_type", "batch"])


def independent_rna_umap(adata: ad.AnnData, seed: int = 0) -> np.ndarray:
    """Neutral RNA visualization. Not scVI/totalVI UMAP."""
    tmp = ad.AnnData(X=adata.layers["counts"].copy() if "counts" in adata.layers else adata.X.copy())
    tmp.obs_names = adata.obs_names
    tmp.var_names = adata.var_names
    import scanpy as sc

    sc.pp.highly_variable_genes(tmp, n_top_genes=2000, flavor="seurat_v3")
    tmp = tmp[:, tmp.var["highly_variable"]].copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.scale(tmp, max_value=10)
    sc.tl.pca(tmp, n_comps=30, random_state=seed)
    sc.pp.neighbors(tmp, n_neighbors=15, random_state=seed)
    sc.tl.umap(tmp, random_state=seed)
    return np.asarray(tmp.obsm["X_umap"])
