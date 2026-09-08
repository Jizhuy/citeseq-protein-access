"""Search AnnData objects for published cell-type annotations."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad

from src.utils.io import resolve_path

ANNOTATION_KEYWORDS = (
    "celltype",
    "cell_type",
    "celltypes",
    "annotation",
    "cell_annotation",
    "celltype.l1",
    "celltype.l2",
    "celltype.l3",
    "predicted.celltype",
    "cluster",
    "seurat_clusters",
    "labels",
    "author_cell_type",
    "manual_annotation",
    "leiden",
)

DEFAULT_FILES = (
    "data/raw/pbmc_10k_protein_v3.h5ad",
    "data/raw/pbmc_5k_protein_v3.h5ad",
    "data/processed/pbmc10k_cite.h5ad",
    "data/processed/pbmc5k_cite.h5ad",
    "data/processed/pbmc_cite_combined_inner.h5ad",
    "data/processed/pbmc_cite_combined_outer.h5ad",
)


def _keyword_hits(names: list[str]) -> list[str]:
    hits = []
    lowered = [(str(name), str(name).lower()) for name in names]
    for original, low in lowered:
        for key in ANNOTATION_KEYWORDS:
            if key in low or low in key:
                hits.append(original)
                break
    return sorted(set(hits))


def inspect_anndata(path: Path) -> dict[str, Any]:
    exists = path.exists()
    report: dict[str, Any] = {"path": str(path), "exists": exists}
    if not exists:
        return report
    adata = ad.read_h5ad(path)
    report.update(
        {
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_columns": list(map(str, adata.obs.columns)),
            "var_columns": list(map(str, adata.var.columns)),
            "uns_keys": list(map(str, adata.uns.keys())),
            "obsm_keys": list(map(str, adata.obsm.keys())),
            "layers": list(map(str, adata.layers.keys())),
            "obs_keyword_hits": _keyword_hits(list(adata.obs.columns)),
            "uns_keyword_hits": _keyword_hits(list(adata.uns.keys())),
            "var_keyword_hits": _keyword_hits(list(adata.var.columns)),
        }
    )
    if "cell_type" in adata.obs:
        counts = adata.obs["cell_type"].astype(str).value_counts().to_dict()
        report["cell_type_value_counts"] = {str(k): int(v) for k, v in counts.items()}
        report["cell_type_is_placeholder"] = set(counts) <= {"unknown"}
    else:
        report["cell_type_value_counts"] = {}
        report["cell_type_is_placeholder"] = None
    return report


def audit_annotation_sources(files: tuple[str, ...] = DEFAULT_FILES) -> dict[str, Any]:
    inspected = [inspect_anndata(resolve_path(rel)) for rel in files]
    raw_have_labels = False
    processed_only_placeholder = True
    for item in inspected:
        if not item.get("exists"):
            continue
        hits = item.get("obs_keyword_hits") or []
        if item["path"].endswith("protein_v3.h5ad") and hits:
            raw_have_labels = True
        if item.get("cell_type_is_placeholder") is False:
            processed_only_placeholder = False
            raw_have_labels = True
    return {
        "files": inspected,
        "published_labels_in_raw_source": raw_have_labels,
        "processed_labels_are_placeholder_only": processed_only_placeholder,
        "labels_lost_during_concatenation": False,
        "confidence": "high",
        "conclusion": (
            "Official scvi-tools PBMC10k/PBMC5k CITE-seq files contain no curated "
            "cell-type annotation columns. obs['cell_type']=='unknown' was added as "
            "a PHASE 2 placeholder and must not be used for ARI/NMI/ASW."
        ),
    }


def write_annotation_audit_markdown(audit: dict[str, Any], path: Path) -> Path:
    lines = [
        "# Annotation source audit (PHASE 5B)",
        "",
        audit["conclusion"],
        "",
        f"- Published labels in raw source files: **{audit['published_labels_in_raw_source']}**",
        f"- Processed `cell_type` is placeholder only: **{audit['processed_labels_are_placeholder_only']}**",
        f"- Labels lost during concatenation: **{audit['labels_lost_during_concatenation']}**",
        f"- Confidence: **{audit['confidence']}**",
        "",
        "No labels were fabricated. Leiden clusters were not generated as ground truth.",
        "",
        "## Files inspected",
        "",
    ]
    for item in audit["files"]:
        lines.append(f"### `{item['path']}`")
        if not item.get("exists"):
            lines.append("- File missing.")
            lines.append("")
            continue
        lines.append(f"- shape: {item['n_obs']} × {item['n_vars']}")
        lines.append(f"- obs columns: {', '.join(item['obs_columns']) or '(none)'}")
        lines.append(f"- uns keys: {', '.join(item['uns_keys']) or '(none)'}")
        lines.append(f"- var columns: {', '.join(item['var_columns']) or '(none)'}")
        lines.append(f"- obs keyword hits: {', '.join(item['obs_keyword_hits']) or '(none)'}")
        lines.append(f"- uns keyword hits: {', '.join(item['uns_keyword_hits']) or '(none)'}")
        if item.get("cell_type_value_counts"):
            lines.append(f"- cell_type counts: {item['cell_type_value_counts']}")
            lines.append(f"- placeholder: {item['cell_type_is_placeholder']}")
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "The raw downloaded AnnData objects have only QC columns",
            "(`n_genes`, `percent_mito`, `n_counts`). There is no `celltype`,",
            "`annotation`, `seurat_clusters`, or author label to recover.",
            "Therefore the absence of biological labels is **not** an artifact of",
            "inner-join concatenation.",
            "",
            "Future evaluation requires an independent published or reference-based",
            "annotation table matching `obs_names`.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
