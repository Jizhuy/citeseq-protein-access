"""PHASE 6 orchestrator: independent annotation + biological evaluation."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scvi
from sklearn.metrics import f1_score

from src.annotation.annotation_qc import (
    celltype_by_batch,
    independent_rna_umap,
    protein_marker_validation,
)
from src.annotation.seurat_v4 import (
    N_HVG,
    assign_confidence_tier,
    choose_confidence_tiers,
    compare_gene_identifiers,
    download_seurat_v4_reference,
    inspect_reference,
    l2_to_l1_map,
    label_distribution,
    select_reference_hvgs,
    train_reference_scanvi,
    transfer_labels,
)
from src.evaluation.phase6 import (
    LEIDEN_RESOLUTIONS,
    MODELS,
    N_BOOT,
    PRIMARY_RESOLUTION,
    RARE_TEST_N,
    celltype_specific_purity,
    classification_bootstrap_f1,
    confusion_frame,
    evaluate_subset,
    heldout_with_predictions,
    load_latents,
    neighbor_overlap_per_cell,
    neighborhood_purity_per_cell,
    paired_bootstrap,
)
from src.evaluation.phase6_plots import (
    plot_celltype_delta,
    plot_metric_bars,
    plot_query_umap_by_label,
    plot_reference_label_distribution,
    plot_restructuring_vs_delta,
)
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

ANNOTATION_SOURCE = "Seurat v4 PBMC CITE-seq (Stuart et al., Cell 2021) via scvi.data.pbmc_seurat_v4_cite_seq; SCANVI/scArches label transfer"


def _choose_batch_key(ref: ad.AnnData) -> str:
    for key in ("orig.ident", "donor", "Donor", "lane", "batch"):
        if key in ref.obs and ref.obs[key].astype(str).nunique() >= 2:
            return key
    ref.obs["annotation_batch"] = "seurat_v4"
    return "annotation_batch"


def run_phase6(
    config: dict[str, Any] | None = None,
    seed: int = 0,
    skip_train: bool = False,
    skip_transfer: bool = False,
) -> dict[str, Any]:
    cfg = config or load_config()
    set_global_seed(seed)
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    annot_dir = resolve_path("results/annotation")
    annot_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    logger = get_logger("phase6", logs_dir / "phase6.log")

    query_path = processed_dir / "pbmc_cite_combined_inner.h5ad"
    query = ad.read_h5ad(query_path)
    original_cell_type = query.obs["cell_type"].astype(str).copy()
    logger.info("Query %s cells x %s genes; original object will not be overwritten.", query.n_obs, query.n_vars)

    logger.info("Downloading / loading Seurat v4 reference.")
    ref = download_seurat_v4_reference(resolve_path("data/raw"))
    ref_info = inspect_reference(ref)
    save_json(ref_info, logs_dir / "seurat_v4_reference_inspection.json")
    logger.info(
        "Reference: %s cells, %s genes, proteins=%s, l1=%s, l2=%s, raw_counts=%s",
        ref_info["n_cells"],
        ref_info["n_genes"],
        ref_info["n_proteins"],
        ref_info.get("has_celltype.l1"),
        ref_info.get("has_celltype.l2"),
        ref_info.get("raw_rna_counts_available"),
    )
    if not ref_info.get("has_celltype.l1") or not ref_info.get("has_celltype.l2"):
        raise RuntimeError("Reference is missing celltype.l1 / celltype.l2.")
    if not ref_info.get("compatible", True) and not ref_info.get("raw_rna_counts_available"):
        logger.warning("Reference RNA may not be raw integer counts; inspect before interpreting.")

    for key, name in (("celltype.l1", "l1"), ("celltype.l2", "l2")):
        dist = label_distribution(ref, key)
        dist.to_csv(annot_dir / f"reference_label_distribution_{name}.csv", index=False)
        plot_reference_label_distribution(
            dist,
            figures_dir / f"phase6_figureA_reference_label_distribution_{name}",
            f"Seurat v4 reference {key} counts",
        )

    gene_report = compare_gene_identifiers(query, ref)
    save_json(gene_report, logs_dir / "gene_identifier_comparison.json")
    logger.info(
        "Gene overlap: intersection=%s query_shared=%.3f ref_shared=%.3f kinds=%s/%s",
        gene_report["n_exact_intersection"],
        gene_report["fraction_query_shared"],
        gene_report["fraction_reference_shared"],
        gene_report["query_gene_id_kind"],
        gene_report["reference_gene_id_kind"],
    )
    if not gene_report["compatible"]:
        raise RuntimeError(
            "Gene identifiers are incompatible or overlap is too small. Label transfer aborted."
        )
    if gene_report["n_query_duplicated"] or gene_report["n_reference_duplicated"]:
        logger.warning(
            "Duplicated gene names: query=%s reference=%s. First occurrence is kept for HVGs.",
            gene_report["n_query_duplicated"],
            gene_report["n_reference_duplicated"],
        )

    shared = pd.Index(query.var_names.astype(str)).intersection(pd.Index(ref.var_names.astype(str)))
    hvgs = select_reference_hvgs(ref, shared, n_top=N_HVG, seed=seed)
    Path(annot_dir / "annotation_hvgs.txt").write_text("\n".join(hvgs) + "\n", encoding="utf-8")
    batch_key = _choose_batch_key(ref)
    logger.info("Reference batch_key=%s; annotation HVGs=%s", batch_key, len(hvgs))

    model_dir = resolve_path("results/models/annotation_reference")
    mapping = l2_to_l1_map(ref)
    mapping.rename_axis("celltype.l2").reset_index(name="celltype.l1").to_csv(
        annot_dir / "l2_to_l1_mapping.csv", index=False
    )

    if skip_transfer and (resolve_path("data/annotations/pbmc_reference_annotations.csv")).exists():
        logger.info("Reusing existing transferred annotations.")
        annot = pd.read_csv(resolve_path("data/annotations/pbmc_reference_annotations.csv"))
        train_info = {"skipped": True}
    else:
        if skip_train and (model_dir / "model.pt").exists():
            logger.info("Loading existing annotation-only SCANVI model.")
            scanvi = scvi.model.SCANVI.load(str(model_dir))
            train_info = {"loaded_existing": True}
        else:
            scanvi, train_info = train_reference_scanvi(
                ref,
                hvgs,
                model_dir,
                labels_key="celltype.l2",
                batch_key=batch_key,
                seed=seed,
            )
        annot = transfer_labels(
            query,
            scanvi,
            hvgs,
            mapping,
            batch_key=batch_key,
            query_surgery_epochs=50,
            seed=seed,
        )
        annot.to_csv(resolve_path("data/annotations/pbmc_reference_annotations.csv"), index=False)

    if list(annot["cell_id"].astype(str)) != list(query.obs_names.astype(str)):
        annot = annot.set_index("cell_id").reindex(query.obs_names.astype(str)).reset_index()
        if annot["predicted_celltype_l1"].isna().any():
            raise RuntimeError("Annotation table does not cover every query cell.")

    tiers = choose_confidence_tiers(annot["prediction_confidence_l2"].to_numpy())
    annot["confidence_tier_l2"] = [
        assign_confidence_tier(s, tiers) for s in annot["prediction_confidence_l2"]
    ]
    annot["confidence_tier_l1"] = [
        assign_confidence_tier(s, tiers) for s in annot["prediction_confidence_l1"]
    ]
    annot.to_csv(resolve_path("data/annotations/pbmc_reference_annotations.csv"), index=False)
    save_json(tiers, annot_dir / "confidence_thresholds.json")
    annot["prediction_confidence_l2"].to_frame().describe().to_csv(
        annot_dir / "confidence_distribution_l2.csv"
    )

    for level, col in (("l1", "predicted_celltype_l1"), ("l2", "predicted_celltype_l2")):
        celltype_by_batch(annot[col], annot["batch"], level).to_csv(
            annot_dir / f"celltype_by_batch_counts_{level}.csv", index=False
        )
        protein_marker_validation(query, annot[col], level).to_csv(
            annot_dir / f"protein_marker_validation_{level}.csv", index=False
        )
    pd.concat(
        [
            pd.read_csv(annot_dir / "protein_marker_validation_l1.csv"),
            pd.read_csv(annot_dir / "protein_marker_validation_l2.csv"),
        ],
        ignore_index=True,
    ).to_csv(annot_dir / "protein_marker_validation.csv", index=False)
    pd.concat(
        [
            pd.read_csv(annot_dir / "celltype_by_batch_counts_l1.csv"),
            pd.read_csv(annot_dir / "celltype_by_batch_counts_l2.csv"),
        ],
        ignore_index=True,
    ).to_csv(annot_dir / "celltype_by_batch_counts.csv", index=False)

    annotated = query.copy()
    annotated.obs["cell_type_l1"] = annot["predicted_celltype_l1"].to_numpy()
    annotated.obs["cell_type_l2"] = annot["predicted_celltype_l2"].to_numpy()
    annotated.obs["annotation_confidence_l1"] = annot["prediction_confidence_l1"].to_numpy()
    annotated.obs["annotation_confidence_l2"] = annot["prediction_confidence_l2"].to_numpy()
    annotated.obs["annotation_tier_l2"] = annot["confidence_tier_l2"].to_numpy()
    annotated.obs["annotation_method"] = annot["annotation_method"].to_numpy()
    annotated.obs["annotation_source"] = annot["reference_source"].to_numpy()
    if not annotated.obs["cell_type"].astype(str).equals(original_cell_type):
        raise RuntimeError("Placeholder obs['cell_type'] was altered.")
    annotated_path = processed_dir / "pbmc_cite_combined_inner_annotated.h5ad"
    annotated.write_h5ad(annotated_path, compression="gzip")
    logger.info("Wrote annotated copy to %s; original inner h5ad unchanged.", annotated_path)

    umap = independent_rna_umap(query, seed=seed)
    np.save(embeddings_dir / "query_independent_rna_umap.npy", umap)
    plot_query_umap_by_label(
        umap,
        annot["predicted_celltype_l1"],
        figures_dir / "phase6_figureB_query_rna_umap_l1",
        "Independent RNA UMAP colored by transferred celltype.l1",
    )

    latents = load_latents(embeddings_dir, query.obs_names, seed=seed)
    high_mask = annot["confidence_tier_l2"].to_numpy() == "high"
    subsets = {
        "all_mapped": np.ones(query.n_obs, dtype=bool),
        "high_confidence": high_mask,
    }
    logger.info(
        "High-confidence cells: %s / %s (threshold=%.3f)",
        int(high_mask.sum()),
        query.n_obs,
        tiers["high_threshold"],
    )

    metric_rows = []
    leiden_rows = []
    confusion_rows = []
    specific_rows = []
    primary_rows = []
    extras_store: dict[str, dict[str, Any]] = {}

    for subset_name, mask in subsets.items():
        for level, col in (("l1", "predicted_celltype_l1"), ("l2", "predicted_celltype_l2")):
            labels = pd.Series(annot[col].to_numpy(), index=query.obs_names, name=level)
            valid = mask.copy()
            latents_sub = {k: v[valid] for k, v in latents.items()}
            labels_valid = labels.loc[valid]
            split_valid = query.obs.loc[valid, "split"].astype(str)
            batch_valid = query.obs.loc[valid, "batch"].astype(str)
            metrics, leiden, extras = evaluate_subset(
                {k: np.asarray(v) for k, v in latents_sub.items()},
                labels_valid.reset_index(drop=True),
                split_valid.reset_index(drop=True),
                batch_valid.reset_index(drop=True),
                label_level=level,
                annotation_subset=subset_name,
                annotation_source=ANNOTATION_SOURCE,
                confidence_rule=tiers["rule"],
                seed=seed,
            )
            metric_rows.append(metrics)
            leiden_rows.append(leiden)
            extras_store[f"{subset_name}:{level}"] = extras
            for model, _k, _m, _p in MODELS:
                pred = extras["predictions"][model]
                confusion_rows.append(
                    confusion_frame(
                        pred["y_test"], pred["knn_pred"], model, level, subset_name, "knn"
                    )
                )
                confusion_rows.append(
                    confusion_frame(
                        pred["y_test"], pred["logreg_pred"], model, level, subset_name, "logreg"
                    )
                )
            spec = celltype_specific_purity(
                {k: np.asarray(v) for k, v in latents_sub.items()},
                labels_valid.reset_index(drop=True),
                np.arange(labels_valid.shape[0]),
                predictions=extras["predictions"],
            )
            spec["label_level"] = level
            spec["annotation_subset"] = subset_name
            specific_rows.append(spec)

            boot = paired_bootstrap(extras["per_cell"], seed=seed, n_boot=N_BOOT)
            f1_knn = classification_bootstrap_f1(
                extras["predictions"]["scVI_matched"]["y_test"],
                extras["predictions"]["scVI_matched"]["knn_pred"],
                extras["predictions"]["totalVI"]["knn_pred"],
                seed=seed,
                metric_name="knn_macro_f1",
            )
            f1_logreg = classification_bootstrap_f1(
                extras["predictions"]["scVI_matched"]["y_test"],
                extras["predictions"]["scVI_matched"]["logreg_pred"],
                extras["predictions"]["totalVI"]["logreg_pred"],
                seed=seed,
                metric_name="logreg_macro_f1",
            )
            point = metrics.set_index("model")
            for _, brow in pd.concat(
                [boot, pd.DataFrame([f1_knn, f1_logreg])], ignore_index=True
            ).iterrows():
                metric = brow["metric"]
                primary_rows.append(
                    {
                        "metric": metric,
                        "label_level": level,
                        "annotation_subset": subset_name,
                        "scVI_matched": float(point.loc["scVI_matched", metric])
                        if metric in point.columns
                        else np.nan,
                        "totalVI": float(point.loc["totalVI", metric])
                        if metric in point.columns
                        else np.nan,
                        "delta_totalVI_minus_scVI": float(brow["delta_mean"]),
                        "bootstrap_ci_low": float(brow["bootstrap_ci_low"]),
                        "bootstrap_ci_high": float(brow["bootstrap_ci_high"]),
                        "n_boot": int(brow["n_boot"]),
                    }
                )

    metrics_all = pd.concat(metric_rows, ignore_index=True)
    metrics_all.to_csv(tables_dir / "biological_representation_metrics.csv", index=False)
    pd.concat(leiden_rows, ignore_index=True).to_csv(
        tables_dir / "leiden_resolution_grid_ari_nmi.csv", index=False
    )
    pd.concat(confusion_rows, ignore_index=True).to_csv(
        tables_dir / "celltype_confusion_matrices.csv", index=False
    )
    pd.concat(specific_rows, ignore_index=True).to_csv(
        tables_dir / "celltype_specific_representation_metrics.csv", index=False
    )
    primary = pd.DataFrame(primary_rows)
    primary.to_csv(tables_dir / "primary_scvi_matched_vs_totalvi.csv", index=False)

    rare_rows = []
    for key, extras in extras_store.items():
        subset_name, level = key.split(":")
        for model, _k, _m, _p in MODELS:
            pred = extras["predictions"][model]
            for cell_type, n in pred["test_class_counts"].items():
                rare_rows.append(
                    {
                        "model": model,
                        "label_level": level,
                        "annotation_subset": subset_name,
                        "cell_type": cell_type,
                        "n_test": int(n),
                        "too_few_test_examples": bool(n < RARE_TEST_N),
                    }
                )
    pd.DataFrame(rare_rows).to_csv(tables_dir / "rare_test_classes.csv", index=False)

    high_l1 = metrics_all[
        (metrics_all["annotation_subset"] == "high_confidence")
        & (metrics_all["label_level"] == "l1")
    ]
    plot_metric_bars(
        high_l1,
        "celltype_ASW",
        figures_dir / "phase6_figureC_asw_highconf",
        "Cell-type ASW (high-confidence)",
    )
    plot_metric_bars(
        high_l1,
        "knn_purity_k15",
        figures_dir / "phase6_figureC_purity_k15_highconf",
        "kNN purity k=15 (high-confidence)",
    )
    plot_metric_bars(
        high_l1,
        "logreg_macro_f1",
        figures_dir / "phase6_figureD_logreg_macro_f1_highconf",
        "Held-out logreg macro-F1 (high-confidence)",
    )

    spec_high_l1 = pd.concat(specific_rows, ignore_index=True)
    spec_high_l1 = spec_high_l1[
        (spec_high_l1["annotation_subset"] == "high_confidence")
        & (spec_high_l1["label_level"] == "l1")
    ]
    plot_celltype_delta(
        spec_high_l1,
        figures_dir / "phase6_figureE_delta_purity_by_celltype_l1",
        "Cell-type-specific Δ k=15 purity (high-confidence l1)",
    )

    high_idx = np.where(high_mask)[0]
    overlap = neighbor_overlap_per_cell(
        latents["scVI_matched"][high_idx], latents["totalVI"][high_idx], k=15
    )
    y_high = annot.loc[high_mask, "predicted_celltype_l1"].astype(str).to_numpy()
    pur_m = neighborhood_purity_per_cell(latents["scVI_matched"][high_idx], y_high, 15)
    pur_t = neighborhood_purity_per_cell(latents["totalVI"][high_idx], y_high, 15)
    delta_p = pur_t - pur_m
    pd.DataFrame(
        {
            "cell_id": query.obs_names[high_idx].astype(str),
            "neighbor_overlap_k15": overlap,
            "purity_scVI_matched": pur_m,
            "purity_totalVI": pur_t,
            "delta_purity": delta_p,
            "cell_type_l1": y_high,
        }
    ).to_csv(tables_dir / "neighborhood_restructuring_vs_delta_purity.csv", index=False)
    plot_restructuring_vs_delta(
        overlap,
        delta_p,
        figures_dir / "phase6_figureF_restructuring_vs_delta_purity",
        "Neighborhood restructuring vs Δ biological purity",
    )

    reloaded = ad.read_h5ad(query_path)
    if "cell_type_l1" in reloaded.obs:
        raise RuntimeError("Original inner h5ad was modified.")
    if not reloaded.obs["cell_type"].astype(str).equals(original_cell_type):
        raise RuntimeError("Original placeholder cell_type changed.")

    provenance = _write_provenance(
        logs_dir / "annotation_provenance.md",
        ref_info=ref_info,
        gene_report=gene_report,
        hvgs=hvgs,
        batch_key=batch_key,
        tiers=tiers,
        n_high=int(high_mask.sum()),
        train_info=train_info,
    )
    logger.info("PHASE 6 complete. Provenance: %s", provenance)
    return {
        "n_reference_cells": ref_info["n_cells"],
        "n_shared_genes": gene_report["n_exact_intersection"],
        "n_hvgs": len(hvgs),
        "n_query": int(query.n_obs),
        "n_high_confidence": int(high_mask.sum()),
        "confidence_rule": tiers["rule"],
        "annotated_h5ad": str(annotated_path),
    }


def _write_provenance(path: Path, **kwargs) -> Path:
    ref_info = kwargs["ref_info"]
    gene_report = kwargs["gene_report"]
    hvgs = kwargs["hvgs"]
    tiers = kwargs["tiers"]
    train_info = kwargs["train_info"]
    text = f"""# Annotation provenance (PHASE 6)

- Date (UTC): {datetime.now(timezone.utc).isoformat()}
- Reference: Seurat v4 multimodal PBMC CITE-seq (Stuart, Butler, et al., Cell 2021; doi:10.1016/j.cell.2021.04.048)
- Loader: `scvi.data.pbmc_seurat_v4_cite_seq(apply_filters=True, aggregate_proteins=True)`
- scvi-tools notes that this object uses GEO GSE164378 UMI counts, not the scTransform assay from the Seurat tutorial page
- Reference cells after published-style filters: {ref_info["n_cells"]}
- Reference genes: {ref_info["n_genes"]}
- Reference proteins: {ref_info["n_proteins"]}
- Annotation levels used: celltype.l1 (coarse), celltype.l2 (fine). celltype.l3 was not transferred.
- Raw RNA counts available on reference X: {ref_info.get("raw_rna_counts_available")}
- Query/reference gene ID kinds: {gene_report["query_gene_id_kind"]} / {gene_report["reference_gene_id_kind"]}
- Exact shared genes: {gene_report["n_exact_intersection"]}
- Fraction of query genes shared: {gene_report["fraction_query_shared"]:.4f}
- Fraction of reference genes shared: {gene_report["fraction_reference_shared"]:.4f}
- Annotation model genes: {len(hvgs)} highly variable genes selected **on the reference only** among shared genes
- Mapping method: new SCVI pretrain + SCANVI on the reference (`celltype.l2`), then scArches query surgery
- The PHASE 3/5 scVI, scVI_matched, and totalVI models were not loaded or used
- Query surgery uses query RNA counts only; classifier weights stay frozen
- Confidence: SCANVI softmax max-class probability for l2; l1 confidence is the summed probability of l2 classes mapping to the predicted l1
- {tiers["rule"]}
- {tiers["justification"]}
- High-confidence query cells: {kwargs["n_high"]}
- Reference batch key: {kwargs["batch_key"]}
- Software: scvi-tools {scvi.__version__}
- Annotation model notes: {train_info}
"""
    path.write_text(text, encoding="utf-8")
    return path
