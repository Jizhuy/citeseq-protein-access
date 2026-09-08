"""Exploratory protein-marker summaries. Do not convert thresholds into labels."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

DEFAULT_MARKERS = ("CD3", "CD4", "CD8a", "CD14", "CD16", "CD19", "CD56")


def protein_marker_summary(
    adata: ad.AnnData,
    *,
    protein_key: str = "protein_expression",
    markers: tuple[str, ...] = DEFAULT_MARKERS,
) -> pd.DataFrame:
    if protein_key not in adata.obsm:
        raise KeyError(f"obsm['{protein_key}'] is missing.")
    protein = adata.obsm[protein_key]
    if not isinstance(protein, pd.DataFrame):
        raise TypeError("protein_expression must be a DataFrame with protein names.")
    cleaned = {str(c): str(c).replace("_TotalSeqB", "") for c in protein.columns}
    rows = []
    for original, short in cleaned.items():
        if short not in markers:
            continue
        values = protein[original].to_numpy()
        batch_labels = adata.obs["batch"].astype(str).to_numpy()
        if len(batch_labels) != len(values):
            raise ValueError("Protein matrix and obs['batch'] have different lengths.")
        for batch in sorted(set(batch_labels)):
            arr = values[batch_labels == batch]
            rows.append(
                {
                    "protein_name_original": original,
                    "protein_name_cleaned": short,
                    "batch": batch,
                    "n_cells": int(len(arr)),
                    "min": float(arr.min()),
                    "q25": float(pd.Series(arr).quantile(0.25)),
                    "median": float(pd.Series(arr).median()),
                    "q75": float(pd.Series(arr).quantile(0.75)),
                    "max": float(arr.max()),
                    "mean": float(arr.mean()),
                    "fraction_zero": float((arr == 0).mean()),
                    "note": "exploratory marker summary only; not a cell-type label",
                }
            )
        rows.append(
            {
                "protein_name_original": original,
                "protein_name_cleaned": short,
                "batch": "all",
                "n_cells": int(len(values)),
                "min": float(values.min()),
                "q25": float(pd.Series(values).quantile(0.25)),
                "median": float(pd.Series(values).median()),
                "q75": float(pd.Series(values).quantile(0.75)),
                "max": float(values.max()),
                "mean": float(values.mean()),
                "fraction_zero": float((values == 0).mean()),
                "note": "exploratory marker summary only; not a cell-type label",
            }
        )
    return pd.DataFrame(rows)


def write_protein_marker_sanity(adata: ad.AnnData, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    table = protein_marker_summary(adata)
    path = output_dir / "protein_marker_quantiles_by_batch.csv"
    table.to_csv(path, index=False)
    readme = output_dir / "README.md"
    readme.write_text(
        "# Exploratory protein-marker summaries\n\n"
        "These tables summarize raw protein counts for commonly used PBMC markers.\n"
        "They are **not** validated cell-type labels.\n\n"
        "Do not threshold these values into ground-truth annotations unless a later\n"
        "research decision explicitly documents and approves that procedure.\n",
        encoding="utf-8",
    )
    return path
