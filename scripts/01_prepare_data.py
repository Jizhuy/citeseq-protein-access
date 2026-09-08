#!/usr/bin/env python
"""PHASE 1-2: inspect the environment, load PBMC CITE-seq, and validate counts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from src.data.load_citeseq import (
    load_citeseq_bundle,
    protein_panel_table,
    summarize_matrix,
)
from src.data.preprocess import (
    dataset_summary_row,
    inventory_metadata,
    make_train_test_split,
    prepare_processed_adata,
)
from src.data.qc_plots import plot_dataset_overview
from src.utils.device import detect_compute_environment
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "config.yaml"),
        help="Path to config.yaml",
    )
    return parser.parse_args()


def _write_adata(adata, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    logger = get_logger("prepare_data", logs_dir / "01_prepare_data.log")

    seed = int(cfg["reproducibility"]["default_seed"])
    set_global_seed(seed)
    env = detect_compute_environment()
    save_json(env, logs_dir / "environment.json")
    logger.info("Python %s | torch %s | scvi %s", env["python"], env["torch_version"], env["scvi_version"])
    logger.info("Accelerator: %s (%s)", env["accelerator"], env["accelerator_detail"])
    if env["import_errors"]:
        logger.warning("Import issues: %s", env["import_errors"])

    bundle = load_citeseq_bundle(cfg)
    store_log1p = bool(cfg["preprocessing"]["store_log1p_layer"])
    processed = {
        key: prepare_processed_adata(bundle[key], store_log1p_layer=store_log1p)
        for key in ("pbmc10k", "pbmc5k", "combined_inner", "combined_outer")
    }

    summary_rows = []
    metadata_frames = []
    matrix_rows = []
    issues = []

    for key, adata in processed.items():
        summary_rows.append(dataset_summary_row(adata, key))
        metadata_frames.append(inventory_metadata(adata, key))
        matrix_rows.append(summarize_matrix(adata.layers["counts"], f"{key}:rna"))
        protein = adata.obsm["protein_counts"]
        observed = adata.obsm["protein_observed_mask"]
        # Unobserved outer-join entries are NaN. They are counted separately in
        # dataset_summary and are filled with 0 only for this numeric dump.
        matrix_rows.append(summarize_matrix(protein.fillna(0).to_numpy(), f"{key}:protein_incl_nan0"))
        if not bool(adata.uns.get("raw_counts_preserved", False)):
            issues.append(f"{key}: raw_counts_preserved flag is false")
        if not summary_rows[-1]["rna_approx_integer"]:
            issues.append(f"{key}: RNA matrix is not approximately integer-valued")
        if summary_rows[-1]["rna_has_negative"] or summary_rows[-1]["protein_has_negative"]:
            issues.append(f"{key}: negative values detected")
        n_unobserved = int((~observed.to_numpy().astype(bool)).sum())
        if key != "combined_outer" and n_unobserved:
            issues.append(f"{key}: unexpected unobserved protein entries ({n_unobserved})")
        if key == "combined_outer" and n_unobserved == 0:
            issues.append("combined_outer: expected unobserved proteins from panel mismatch, found none")
        logger.info(
            "%s: %s cells x %s genes x %s proteins | RNA sparsity=%.4f | protein sparsity=%.4f",
            key,
            adata.n_obs,
            adata.n_vars,
            protein.shape[1],
            summary_rows[-1]["rna_sparsity"],
            summary_rows[-1]["protein_sparsity"],
        )
        if not summary_rows[-1]["cell_type_labels_available"]:
            issues.append(
                f"{key}: no usable cell-type labels in obs "
                f"(detected column={adata.uns.get('celltype_column_detected')})"
            )

    summary = pd.DataFrame(summary_rows)
    metadata = pd.concat(metadata_frames, ignore_index=True)
    protein_panel = protein_panel_table(processed["pbmc10k"], processed["pbmc5k"])
    matrix_report = pd.DataFrame(matrix_rows)

    split = make_train_test_split(
        processed["combined_inner"],
        test_fraction=float(cfg["splits"]["test_fraction"]),
        seed=seed,
        stratify_by=cfg["splits"]["stratify_by"],
    )
    processed["combined_inner"].obs["split"] = split.set_index("cell_id").loc[
        processed["combined_inner"].obs_names, "split"
    ]
    processed["combined_outer"].obs["split"] = processed["combined_inner"].obs["split"].to_numpy()

    tables_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(tables_dir / "dataset_summary.csv", index=False)
    metadata.to_csv(tables_dir / "metadata_inventory.csv", index=False)
    protein_panel.to_csv(tables_dir / "protein_panel_overlap.csv", index=False)
    matrix_report.to_csv(tables_dir / "matrix_numeric_report.csv", index=False)
    split.to_csv(tables_dir / "combined_inner_train_test_split.csv", index=False)

    _write_adata(processed["pbmc10k"], processed_dir / "pbmc10k_cite.h5ad")
    _write_adata(processed["pbmc5k"], processed_dir / "pbmc5k_cite.h5ad")
    _write_adata(processed["combined_inner"], processed_dir / "pbmc_cite_combined_inner.h5ad")
    _write_adata(processed["combined_outer"], processed_dir / "pbmc_cite_combined_outer.h5ad")

    plot_dataset_overview(
        processed["pbmc10k"],
        processed["pbmc5k"],
        processed["combined_inner"],
        protein_panel,
        summary,
        figures_dir / "dataset_overview_phase2",
    )

    unique_issues = sorted(set(issues))
    report = {
        "phase": 2,
        "seed": seed,
        "environment": env,
        "files_written": {
            "processed": [
                "data/processed/pbmc10k_cite.h5ad",
                "data/processed/pbmc5k_cite.h5ad",
                "data/processed/pbmc_cite_combined_inner.h5ad",
                "data/processed/pbmc_cite_combined_outer.h5ad",
            ],
            "tables": [
                "results/tables/dataset_summary.csv",
                "results/tables/metadata_inventory.csv",
                "results/tables/protein_panel_overlap.csv",
                "results/tables/matrix_numeric_report.csv",
                "results/tables/combined_inner_train_test_split.csv",
            ],
        },
        "summary": summary.to_dict(orient="records"),
        "protein_panel": {
            "n_proteins_pbmc10k": int(protein_panel["in_pbmc10k"].sum()),
            "n_proteins_pbmc5k": int(protein_panel["in_pbmc5k"].sum()),
            "n_shared": int(protein_panel["shared"].sum()),
            "pbmc10k_only": protein_panel.loc[
                protein_panel["in_pbmc10k"] & ~protein_panel["in_pbmc5k"], "protein"
            ].tolist(),
            "pbmc5k_only": protein_panel.loc[
                protein_panel["in_pbmc5k"] & ~protein_panel["in_pbmc10k"], "protein"
            ].tolist(),
            "shared": protein_panel.loc[protein_panel["shared"], "protein"].tolist(),
        },
        "issues": unique_issues,
        "notes": [
            "X and layers['counts'] store raw RNA counts; they are not overwritten with normalized values.",
            "Protein raw counts are in obsm['protein_counts'] and copied to obsm['protein_expression'].",
            "combined_inner uses shared genes and shared proteins only.",
            "combined_outer keeps the union of proteins; unmeasured proteins are NaN, not zero.",
            "No model training is performed in this script.",
        ],
    }
    save_json(report, logs_dir / "phase2_data_report.json")
    pd.Series(unique_issues, name="issue").to_csv(tables_dir / "data_issues.csv", index=False)

    logger.info("Wrote processed AnnData objects to %s", processed_dir)
    logger.info("Issues recorded: %s", unique_issues or ["none"])
    logger.info("PHASE 2 complete. Do not train scVI/totalVI yet unless starting PHASE 3.")


if __name__ == "__main__":
    main()
