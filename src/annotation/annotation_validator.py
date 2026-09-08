"""Validate external annotations against AnnData cell IDs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

PLACEHOLDER_LABELS = {"unknown", "", "nan", "none", "na", "n/a"}


class AnnotationValidationError(ValueError):
    """Raised when an annotation table is unsafe to attach or evaluate."""


@dataclass
class AnnotationValidationReport:
    n_table_rows: int
    n_adata_cells: int
    n_unique_table_ids: int
    n_matched: int
    n_unmatched_in_table: int
    n_unmatched_in_adata: int
    n_duplicate_table_ids: int
    n_missing_labels: int
    unmatched_in_table: list[str] = field(default_factory=list)
    unmatched_in_adata: list[str] = field(default_factory=list)
    duplicate_cell_ids: list[str] = field(default_factory=list)
    label_frequencies: dict[str, int] = field(default_factory=dict)
    placeholder_labels_present: bool = False
    errors: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_table_rows": self.n_table_rows,
            "n_adata_cells": self.n_adata_cells,
            "n_unique_table_ids": self.n_unique_table_ids,
            "n_matched": self.n_matched,
            "n_unmatched_in_table": self.n_unmatched_in_table,
            "n_unmatched_in_adata": self.n_unmatched_in_adata,
            "n_duplicate_table_ids": self.n_duplicate_table_ids,
            "n_missing_labels": self.n_missing_labels,
            "unmatched_in_table": self.unmatched_in_table,
            "unmatched_in_adata": self.unmatched_in_adata,
            "duplicate_cell_ids": self.duplicate_cell_ids,
            "label_frequencies": self.label_frequencies,
            "placeholder_labels_present": self.placeholder_labels_present,
            "errors": self.errors,
            "ok": self.ok,
        }


def labels_are_validated(labels: pd.Series | np.ndarray) -> bool:
    """True only if there are at least two non-placeholder classes."""
    values = pd.Series(labels).astype(str)
    values = values[values.notna()]
    unique = {v.strip().lower() for v in values}
    unique -= PLACEHOLDER_LABELS
    return len(unique) >= 2


def validate_annotation_table(
    table: pd.DataFrame,
    adata_obs_names: pd.Index,
    require_full_coverage: bool = True,
) -> AnnotationValidationReport:
    """Check uniqueness, overlap, missing labels, and frequencies.

    Does not silently drop unmatched cells. Failures are collected in `errors`.
    """
    errors: list[str] = []
    cell_ids = table["cell_id"].astype(str)
    labels = table["cell_type"].astype(str)
    adata_ids = pd.Index(adata_obs_names.astype(str))

    missing_mask = labels.isna() | labels.str.strip().eq("") | labels.str.lower().isin(
        PLACEHOLDER_LABELS
    )
    n_missing = int(missing_mask.sum())
    if n_missing:
        errors.append(f"{n_missing} rows have missing or placeholder cell_type values.")

    duplicated = cell_ids[cell_ids.duplicated(keep=False)]
    duplicate_ids = sorted(duplicated.unique().tolist())
    if duplicate_ids:
        errors.append(
            f"{len(duplicate_ids)} duplicated cell_id values. Unmatched/duplicate rows were not dropped."
        )

    table_set = set(cell_ids)
    adata_set = set(adata_ids)
    unmatched_table = sorted(table_set - adata_set)
    unmatched_adata = sorted(adata_set - table_set)
    if unmatched_table:
        errors.append(
            f"{len(unmatched_table)} annotation cell_id values are absent from AnnData.obs_names. "
            "They were not dropped."
        )
    if require_full_coverage and unmatched_adata:
        errors.append(
            f"{len(unmatched_adata)} AnnData cells have no annotation row. "
            "They were not silently excluded."
        )

    if not labels_are_validated(labels):
        errors.append(
            "Validated ground-truth cell-type labels unavailable: "
            "fewer than two non-placeholder classes."
        )

    frequencies = labels.value_counts(dropna=False).to_dict()
    frequencies = {str(k): int(v) for k, v in frequencies.items()}
    return AnnotationValidationReport(
        n_table_rows=int(len(table)),
        n_adata_cells=int(len(adata_ids)),
        n_unique_table_ids=int(cell_ids.nunique()),
        n_matched=int(len(table_set & adata_set)),
        n_unmatched_in_table=len(unmatched_table),
        n_unmatched_in_adata=len(unmatched_adata),
        n_duplicate_table_ids=len(duplicate_ids),
        n_missing_labels=n_missing,
        unmatched_in_table=unmatched_table[:50],
        unmatched_in_adata=unmatched_adata[:50],
        duplicate_cell_ids=duplicate_ids[:50],
        label_frequencies=frequencies,
        placeholder_labels_present=bool(
            any(str(v).strip().lower() in PLACEHOLDER_LABELS for v in frequencies)
        ),
        errors=errors,
    )


def assert_valid_annotations(
    table: pd.DataFrame,
    adata_obs_names: pd.Index,
    require_full_coverage: bool = True,
) -> AnnotationValidationReport:
    report = validate_annotation_table(
        table, adata_obs_names, require_full_coverage=require_full_coverage
    )
    if not report.ok:
        raise AnnotationValidationError("; ".join(report.errors))
    return report
