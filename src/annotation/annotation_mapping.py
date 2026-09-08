"""Attach validated external labels onto AnnData without inventing values."""

from __future__ import annotations

from typing import Any

import anndata as ad
import pandas as pd

from src.annotation.annotation_validator import (
    AnnotationValidationError,
    assert_valid_annotations,
)


def attach_annotations(
    adata: ad.AnnData,
    table: pd.DataFrame,
    *,
    validated_key: str = "cell_type_validated",
    coarse_key: str = "cell_type_coarse",
    fine_key: str = "cell_type_fine",
    overwrite_placeholder_cell_type: bool = False,
    copy: bool = True,
) -> ad.AnnData:
    """Add validated labels as new obs columns.

    By default the PHASE 2 placeholder `obs['cell_type']` is left unchanged.
    Unmatched cells are not dropped; validation fails instead.
    """
    report = assert_valid_annotations(table, adata.obs_names, require_full_coverage=True)
    if not report.ok:
        raise AnnotationValidationError("; ".join(report.errors))

    out = adata.copy() if copy else adata
    indexed = table.drop_duplicates(subset=["cell_id"], keep="first").set_index("cell_id")
    aligned = indexed.reindex(out.obs_names.astype(str))
    out.obs[validated_key] = aligned["cell_type"].to_numpy()
    if "cell_type_coarse" in aligned.columns:
        out.obs[coarse_key] = aligned["cell_type_coarse"].to_numpy()
    if "cell_type_fine" in aligned.columns:
        out.obs[fine_key] = aligned["cell_type_fine"].to_numpy()
    if "annotation_source" in aligned.columns:
        out.obs["annotation_source"] = aligned["annotation_source"].to_numpy()
    if "annotation_confidence" in aligned.columns:
        out.obs["annotation_confidence"] = aligned["annotation_confidence"].to_numpy()
    if overwrite_placeholder_cell_type:
        out.obs["cell_type"] = out.obs[validated_key]
    out.uns["annotation_attached"] = {
        "validated_key": validated_key,
        "overwrite_placeholder_cell_type": overwrite_placeholder_cell_type,
        "n_cells": int(out.n_obs),
        "n_classes": int(pd.Series(out.obs[validated_key]).nunique()),
    }
    return out


def select_label_level(obs: pd.DataFrame, level: str) -> pd.Series:
    """Return coarse, fine, or generic validated labels. No hard-coded class list."""
    mapping = {
        "coarse": ["cell_type_coarse", "cell_type_validated", "cell_type"],
        "fine": ["cell_type_fine", "cell_type_validated", "cell_type"],
        "validated": ["cell_type_validated", "cell_type"],
    }
    if level not in mapping:
        raise KeyError(f"Unknown label level {level}. Use coarse, fine, or validated.")
    for column in mapping[level]:
        if column in obs.columns:
            return obs[column].astype(str)
    raise KeyError(f"No annotation column available for level={level}.")


def annotation_summary(labels: pd.Series) -> dict[str, Any]:
    counts = labels.astype(str).value_counts(dropna=False)
    return {
        "n_cells": int(len(labels)),
        "n_classes": int(counts.shape[0]),
        "label_frequencies": {str(k): int(v) for k, v in counts.items()},
    }
