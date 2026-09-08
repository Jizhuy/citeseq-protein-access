"""Load externally supplied cell-type annotation tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ("cell_id", "cell_type")
OPTIONAL_COLUMNS = (
    "cell_type_coarse",
    "cell_type_fine",
    "annotation_source",
    "annotation_confidence",
)
TEMPLATE_COLUMNS = (
    "cell_id",
    "cell_type_coarse",
    "cell_type_fine",
    "annotation_source",
    "annotation_confidence",
)


class AnnotationLoadError(ValueError):
    """Raised when an annotation table cannot be read as specified."""


def load_annotation_table(path: str | Path) -> pd.DataFrame:
    """Load a CSV with at least cell_id and cell_type.

    Unmatched cells are not dropped here. Validation is a separate step.
    """
    table_path = Path(path)
    if not table_path.exists():
        raise AnnotationLoadError(f"Annotation file does not exist: {table_path}")
    table = pd.read_csv(table_path)
    if table.empty and list(table.columns):
        # Header-only template is valid to load, but not valid to evaluate.
        pass
    if "cell_id" not in table.columns:
        raise AnnotationLoadError(
            f"Annotation table is missing cell_id. Found columns: {list(table.columns)}"
        )
    out = table.copy()
    if "cell_type" not in out.columns:
        if "cell_type_fine" in out.columns:
            out["cell_type"] = out["cell_type_fine"]
        elif "cell_type_coarse" in out.columns:
            out["cell_type"] = out["cell_type_coarse"]
        else:
            raise AnnotationLoadError(
                "Annotation table needs cell_type, or cell_type_fine / cell_type_coarse. "
                f"Found columns: {list(table.columns)}"
            )
    out["cell_id"] = out["cell_id"].astype(str)
    out["cell_type"] = out["cell_type"].astype(str)
    for column in OPTIONAL_COLUMNS:
        if column in out.columns:
            out[column] = out[column].astype(str)
    return out
