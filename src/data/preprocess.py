"""Validate CITE-seq AnnData objects without overwriting raw counts."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from src.data.load_citeseq import summarize_matrix
from src.utils.logging_utils import get_logger

logger = get_logger(__name__)


def add_library_metrics(adata: ad.AnnData) -> ad.AnnData:
    rna = adata.layers.get("counts", adata.X)
    if hasattr(rna, "sum"):
        rna_counts = np.asarray(rna.sum(axis=1)).ravel()
    else:
        rna_counts = np.asarray(rna).sum(axis=1)
    adata.obs["rna_counts"] = rna_counts
    adata.obs["n_genes_by_counts"] = np.asarray((rna > 0).sum(axis=1)).ravel()

    protein = adata.obsm["protein_counts"]
    observed = adata.obsm.get("protein_observed_mask")
    if observed is None:
        observed = pd.DataFrame(True, index=protein.index, columns=protein.columns)
    numeric = protein.where(observed, other=np.nan)
    adata.obs["protein_counts_sum"] = numeric.sum(axis=1, skipna=True).to_numpy()
    adata.obs["n_proteins_detected"] = (numeric > 0).sum(axis=1, skipna=True).to_numpy()
    adata.obs["n_proteins_measured"] = observed.sum(axis=1).to_numpy()
    return adata


def store_optional_log1p(adata: ad.AnnData) -> ad.AnnData:
    """Store a float32 visualization layer without modifying X or counts.

    The full gene matrix becomes dense after library-size normalization, so
    this is opt-in. Phase 2 defaults to leaving it off.
    """
    counts = adata.layers["counts"]
    if hasattr(counts, "toarray"):
        dense = counts.toarray()
    else:
        dense = np.asarray(counts)
    lib = dense.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    adata.layers["log1p_norm"] = np.log1p(dense / lib * 1e4).astype(np.float32)
    return adata


def prepare_processed_adata(
    adata: ad.AnnData,
    store_log1p_layer: bool = True,
) -> ad.AnnData:
    out = adata.copy()
    if "counts" not in out.layers:
        out.layers["counts"] = out.X.copy()
    add_library_metrics(out)
    if store_log1p_layer:
        store_optional_log1p(out)
    out.uns["x_is_raw_counts"] = True
    out.uns["raw_counts_preserved"] = True
    return out


def inventory_metadata(adata: ad.AnnData, dataset: str) -> pd.DataFrame:
    rows = []
    for col in adata.obs.columns:
        series = adata.obs[col]
        n_unique = int(series.nunique(dropna=False))
        sample_vals = series.dropna().astype(str).unique()[:8].tolist()
        rows.append(
            {
                "dataset": dataset,
                "slot": "obs",
                "column": col,
                "dtype": str(series.dtype),
                "n_unique": n_unique,
                "n_missing": int(series.isna().sum()),
                "example_values": "; ".join(sample_vals),
            }
        )
    for col in adata.var.columns:
        series = adata.var[col]
        rows.append(
            {
                "dataset": dataset,
                "slot": "var",
                "column": col,
                "dtype": str(series.dtype),
                "n_unique": int(series.nunique(dropna=False)),
                "n_missing": int(series.isna().sum()),
                "example_values": "; ".join(series.dropna().astype(str).unique()[:8].tolist()),
            }
        )
    return pd.DataFrame(rows)


def dataset_summary_row(adata: ad.AnnData, dataset: str) -> dict[str, Any]:
    rna_stats = summarize_matrix(adata.layers.get("counts", adata.X), "rna")
    protein = adata.obsm["protein_counts"]
    observed = adata.obsm.get("protein_observed_mask")
    if observed is None:
        observed = pd.DataFrame(True, index=protein.index, columns=protein.columns)
    protein_numeric = protein.where(observed, other=0.0)
    protein_stats = summarize_matrix(protein_numeric.to_numpy(), "protein")
    n_unobserved = int((~observed.to_numpy().astype(bool)).sum())
    celltype_col = adata.uns.get("celltype_column_detected")
    cell_types = adata.obs["cell_type"].astype(str)
    unknown_only = set(cell_types.unique()) <= {"unknown"}
    batches = sorted(adata.obs["batch"].astype(str).unique().tolist())
    return {
        "dataset": dataset,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_proteins": int(protein.shape[1]),
        "n_proteins_fully_observed": int(observed.all(axis=0).sum()),
        "n_protein_entries_unobserved": n_unobserved,
        "batches": "; ".join(batches),
        "n_batches": len(batches),
        "celltype_column_detected": celltype_col,
        "cell_type_labels_available": bool(celltype_col) and not unknown_only,
        "n_cell_types": int(cell_types.nunique()),
        "cell_types": "; ".join(sorted(cell_types.unique().tolist())[:20]),
        "rna_sparsity": rna_stats["sparsity"],
        "protein_sparsity": protein_stats["sparsity"],
        "rna_mean_counts_per_cell": rna_stats["mean_counts_per_cell"],
        "protein_mean_counts_per_cell": protein_stats["mean_counts_per_cell"],
        "rna_approx_integer": rna_stats["approx_integer"],
        "protein_approx_integer": protein_stats["approx_integer"],
        "rna_has_negative": rna_stats["has_negative"],
        "protein_has_negative": protein_stats["has_negative"],
        "rna_layout": rna_stats["layout"],
        "raw_counts_preserved": bool(adata.uns.get("raw_counts_preserved", False)),
        "protein_join": adata.uns.get("protein_join", "single_dataset"),
        "n_genes_dropped_pbmc10k": adata.uns.get("n_genes_dropped_pbmc10k", 0),
        "n_genes_dropped_pbmc5k": adata.uns.get("n_genes_dropped_pbmc5k", 0),
        "n_proteins_dropped_pbmc10k": adata.uns.get("n_proteins_dropped_pbmc10k", 0),
        "n_proteins_dropped_pbmc5k": adata.uns.get("n_proteins_dropped_pbmc5k", 0),
        "obs_columns": "; ".join(map(str, adata.obs.columns)),
        "obsm_keys": "; ".join(map(str, adata.obsm.keys())),
        "layer_keys": "; ".join(map(str, adata.layers.keys())),
    }


def make_train_test_split(
    adata: ad.AnnData,
    test_fraction: float,
    seed: int,
    stratify_by: str | None = "batch",
) -> pd.DataFrame:
    index = np.arange(adata.n_obs)
    stratify = None
    if stratify_by and stratify_by in adata.obs:
        stratify = adata.obs[stratify_by].astype(str).to_numpy()
        counts = pd.Series(stratify).value_counts()
        if (counts < 2).any():
            logger.warning(
                "Stratification column %s has classes with <2 cells; using an unstratified split.",
                stratify_by,
            )
            stratify = None
    train_idx, test_idx = train_test_split(
        index,
        test_size=test_fraction,
        random_state=seed,
        stratify=stratify,
    )
    split = pd.Series("train", index=adata.obs_names, name="split")
    split.iloc[test_idx] = "test"
    frame = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "split": split.to_numpy(),
            "seed": seed,
            "test_fraction": test_fraction,
        }
    )
    return frame
