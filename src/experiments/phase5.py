"""PHASE 5A/5B post-training analysis. Does not train MOFA+ or invent labels."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from src.annotation.audit import audit_annotation_sources, write_annotation_audit_markdown
from src.evaluation.biological_metrics import (
    CELLTYPE_UNAVAILABLE_REASON,
    write_biological_metrics_schema,
)
from src.evaluation.latent_geometry import compare_latent_geometries
from src.evaluation.protein_marker_sanity import write_protein_marker_sanity
from src.utils.io import load_config, resolve_path, save_json

CELLTYPE_NA_REASON = "Validated ground-truth cell-type labels unavailable"


def _load_json(path: Path) -> dict[str, Any]:
    import json

    return json.loads(path.read_text(encoding="utf-8"))


def _require_cell_order(cells_csv: Path, reference: pd.Index, name: str) -> None:
    cells = pd.read_csv(cells_csv)
    if list(cells["cell_id"].astype(str)) != list(reference.astype(str)):
        raise ValueError(f"{name} cell order does not match the source AnnData.")


def write_three_model_table(cfg: dict[str, Any] | None = None, seed: int = 0) -> Path:
    cfg = cfg or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    scvi = _load_json(logs_dir / f"scvi_seed{seed}_manifest.json")
    matched = _load_json(logs_dir / f"scvi_matched_seed{seed}_manifest.json")
    totalvi = _load_json(logs_dir / f"totalvi_seed{seed}_manifest.json")
    rows = [
        {
            "model": "scVI_default",
            "modalities": "RNA",
            "gene_likelihood": scvi["hyperparameters"]["gene_likelihood"],
            "n_hidden": scvi["hyperparameters"]["n_hidden"],
            "n_layers": scvi["hyperparameters"]["n_layers"],
            "dropout": scvi["hyperparameters"]["dropout_rate"],
            "latent_dim": scvi["n_latent"],
            "n_cells": scvi["n_cells"],
            "n_genes": scvi["n_genes"],
            "n_proteins": 0,
            "epochs": scvi["epochs_completed"],
            "runtime_seconds": scvi["runtime_seconds"],
            "peak_rss_mb": scvi["peak_rss_mb"],
            "device": scvi["actual_device"],
            "batch_silhouette": scvi["batch_silhouette"],
            "ARI": np.nan,
            "NMI": np.nan,
            "celltype_ASW": np.nan,
            "celltype_metrics_reason": CELLTYPE_NA_REASON,
        },
        {
            "model": "scVI_matched",
            "modalities": "RNA",
            "gene_likelihood": matched["gene_likelihood"],
            "n_hidden": matched["n_hidden"],
            "n_layers": matched["n_layers"],
            "dropout": matched["dropout_rate"],
            "latent_dim": matched["n_latent"],
            "n_cells": matched["n_cells"],
            "n_genes": matched["n_genes"],
            "n_proteins": 0,
            "epochs": matched["epochs_completed"],
            "runtime_seconds": matched["runtime_seconds"],
            "peak_rss_mb": matched["peak_rss_mb"],
            "device": matched["actual_device"],
            "batch_silhouette": matched["batch_silhouette"],
            "ARI": np.nan,
            "NMI": np.nan,
            "celltype_ASW": np.nan,
            "celltype_metrics_reason": CELLTYPE_NA_REASON,
        },
        {
            "model": "totalVI",
            "modalities": "RNA+Protein",
            "gene_likelihood": totalvi["effective_model_hyperparameters"]["gene_likelihood"],
            "n_hidden": totalvi["effective_model_hyperparameters"]["n_hidden"],
            "n_layers": (
                f"encoder={totalvi['effective_model_hyperparameters']['n_layers_encoder']},"
                f"decoder={totalvi['effective_model_hyperparameters']['n_layers_decoder']}"
            ),
            "dropout": totalvi["effective_model_hyperparameters"]["dropout_rate_encoder"],
            "latent_dim": totalvi["n_latent"],
            "n_cells": totalvi["n_cells"],
            "n_genes": totalvi["n_genes"],
            "n_proteins": totalvi["n_proteins"],
            "epochs": totalvi["epochs_completed"],
            "runtime_seconds": totalvi["runtime_seconds"],
            "peak_rss_mb": totalvi["peak_rss_mb"],
            "device": totalvi["actual_device"],
            "batch_silhouette": totalvi["batch_silhouette"],
            "ARI": np.nan,
            "NMI": np.nan,
            "celltype_ASW": np.nan,
            "celltype_metrics_reason": CELLTYPE_NA_REASON,
        },
    ]
    path = tables_dir / f"scvi_default_vs_matched_vs_totalvi_seed{seed}.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def write_geometry_comparison(cfg: dict[str, Any] | None = None, seed: int = 0) -> Path:
    cfg = cfg or load_config()
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    processed = resolve_path(cfg["paths"]["processed_dir"]) / "pbmc_cite_combined_inner.h5ad"
    adata = ad.read_h5ad(processed)
    names = {
        "X_scVI": "scvi",
        "X_scVI_matched": "scvi_matched",
        "X_totalVI": "totalvi",
    }
    latents = {}
    for key, prefix in names.items():
        npy = embeddings_dir / f"{prefix}_seed{seed}_latent.npy"
        cells = embeddings_dir / f"{prefix}_seed{seed}_cells.csv"
        if not npy.exists() or not cells.exists():
            raise FileNotFoundError(f"Missing {npy} or {cells}")
        _require_cell_order(cells, adata.obs_names, prefix)
        latents[key] = np.load(npy)
        if latents[key].shape[0] != adata.n_obs:
            raise ValueError(f"{key} has {latents[key].shape[0]} rows, expected {adata.n_obs}.")
    table = compare_latent_geometries(latents, seed=seed, knn_ks=(15, 30, 50))
    path = tables_dir / f"latent_geometry_comparison_seed{seed}.csv"
    table.to_csv(path, index=False)
    return path


def write_fairness_note(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """# Fairness-control interpretation (PHASE 5A)

The comparison **scVI_default vs totalVI** mixes several changes at once:

- protein modality (absent vs present)
- RNA likelihood (ZINB vs NB)
- network capacity (n_hidden 128 / n_layers 1 / dropout 0.1 vs 256 / encoder 2 / 0.2)

A difference between those two models therefore cannot be attributed cleanly to
the protein measurements.

The comparison **scVI_matched vs totalVI** is the preferred comparison for
estimating the incremental effect of adding protein information. It is a
**more closely matched control**, not a perfectly identical architecture.

Matched where SCVI 1.3.3 permits:

- latent dimension = 20
- RNA likelihood = NB
- hidden width = 256
- encoder depth = 2
- dropout = 0.2
- same cells, genes, batch variable, seed, batch size, and 200-epoch budget

Still structurally different:

- totalVI has protein-specific decoder / background components
- totalVI is a two-modality generative model
- SCVI `n_layers=2` applies to **both** encoder and decoder; totalVI uses
  encoder=2 and decoder=1
- SCVI training uses the SCVI plan (typical lr 1e-3, no `reduce_lr_on_plateau`);
  totalVI uses lr 4e-3 and reduce-on-plateau
- ELBO values are not the same objective and must not be compared numerically

Do not rank models biologically from batch silhouette or latent-geometry
similarity. Validated cell-type labels are still unavailable.
""",
        encoding="utf-8",
    )
    return path


def run_phase5_analysis(config: dict[str, Any] | None = None, seed: int = 0) -> dict[str, Any]:
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    processed = resolve_path(cfg["paths"]["processed_dir"]) / "pbmc_cite_combined_inner.h5ad"

    comparison = write_three_model_table(cfg, seed=seed)
    geometry = write_geometry_comparison(cfg, seed=seed)
    fairness = write_fairness_note(logs_dir / "fairness_control_interpretation.md")

    audit = audit_annotation_sources()
    audit_md = write_annotation_audit_markdown(audit, logs_dir / "annotation_source_audit.md")
    save_json(audit, logs_dir / "annotation_source_audit.json")

    schema = tables_dir / "biological_representation_metrics.csv"
    write_biological_metrics_schema(schema)

    adata = ad.read_h5ad(processed)
    exploratory = write_protein_marker_sanity(
        adata, resolve_path("results/exploratory")
    )

    return {
        "three_model_table": str(comparison),
        "geometry_table": str(geometry),
        "fairness_note": str(fairness),
        "annotation_audit": str(audit_md),
        "biological_metrics_schema": str(schema),
        "biological_metrics_status": "Pending validated annotations",
        "published_labels_recovered": bool(audit["published_labels_in_raw_source"]),
        "exploratory_protein_summary": str(exploratory),
        "celltype_metrics_reason": CELLTYPE_UNAVAILABLE_REASON,
    }
