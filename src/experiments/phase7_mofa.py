"""PHASE 7: MOFA+ baseline + biological comparison. Does not run sparsity tests."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from src.evaluation.phase6 import (
    MODELS,
    N_BOOT,
    celltype_specific_purity,
    classification_bootstrap_f1,
    confusion_frame,
    evaluate_subset,
    load_latents,
    paired_bootstrap,
)
from src.evaluation.phase6_plots import plot_metric_bars, plot_query_umap_by_label
from src.models.run_mofa import run_mofa
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

PHASE7_MODELS = MODELS + (
    ("PCA_RNA", "X_PCA_RNA", "RNA", "pca_rna"),
    ("MOFA+", "X_MOFA", "RNA+Protein", "mofa"),
)
ANNOTATION_SOURCE = (
    "Seurat v4 PBMC CITE-seq (Stuart et al., Cell 2021) via PHASE 6 SCANVI/scArches; "
    "labels independent of MOFA+/PCA/VAE latents"
)
COMPARE_PAIRS = (
    ("scVI_matched", "totalVI"),
    ("PCA_RNA", "MOFA+"),
    ("MOFA+", "totalVI"),
    ("scVI_matched", "MOFA+"),
)


def run_phase7(config: dict[str, Any] | None = None, seed: int = 0, skip_train: bool = False) -> dict[str, Any]:
    cfg = config or load_config()
    set_global_seed(seed)
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    logger = get_logger("phase7", logs_dir / "phase7.log")

    if skip_train and (embeddings_dir / f"mofa_seed{seed}_latent.npy").exists():
        logger.info("Reusing existing MOFA+ / PCA embeddings.")
        train_info = {"skipped": True}
    else:
        train_info = run_mofa(config=cfg, seed=seed)

    query = ad.read_h5ad(processed_dir / "pbmc_cite_combined_inner.h5ad")
    annotated = ad.read_h5ad(processed_dir / "pbmc_cite_combined_inner_annotated.h5ad")
    if list(annotated.obs_names.astype(str)) != list(query.obs_names.astype(str)):
        raise RuntimeError("Annotated object cell order does not match the inner-join query.")
    if "cell_type_l1" in query.obs:
        raise RuntimeError("Original inner h5ad contains annotation columns.")

    latents = load_latents(embeddings_dir, query.obs_names, seed=seed, model_specs=PHASE7_MODELS)
    high_mask = annotated.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    confidence_rule = "PHASE 6 high-confidence rule: l2_confidence >= 0.85"
    logger.info("High-confidence cells: %s / %s", int(high_mask.sum()), query.n_obs)

    umap = np.load(embeddings_dir / f"mofa_seed{seed}_umap.npy")
    plot_query_umap_by_label(
        umap,
        annotated.obs["cell_type_l1"],
        figures_dir / "phase7_mofa_umap_l1",
        "MOFA+ UMAP colored by independent transferred celltype.l1",
    )

    subsets = {"all_mapped": np.ones(query.n_obs, dtype=bool), "high_confidence": high_mask}
    metric_rows = []
    leiden_rows = []
    confusion_rows = []
    specific_rows = []
    primary_rows = []

    for subset_name, mask in subsets.items():
        for level, col in (("l1", "cell_type_l1"), ("l2", "cell_type_l2")):
            labels = annotated.obs[col].astype(str)
            labels.name = level
            latents_sub = {k: v[mask] for k, v in latents.items()}
            labels_valid = labels.loc[mask]
            split_valid = query.obs.loc[mask, "split"].astype(str)
            batch_valid = query.obs.loc[mask, "batch"].astype(str)
            metrics, leiden, extras = evaluate_subset(
                {k: np.asarray(v) for k, v in latents_sub.items()},
                labels_valid.reset_index(drop=True),
                split_valid.reset_index(drop=True),
                batch_valid.reset_index(drop=True),
                label_level=level,
                annotation_subset=subset_name,
                annotation_source=ANNOTATION_SOURCE,
                confidence_rule=confidence_rule,
                seed=seed,
                model_specs=PHASE7_MODELS,
            )
            metric_rows.append(metrics)
            leiden_rows.append(leiden)
            for model, _k, _m, _p in PHASE7_MODELS:
                pred = extras["predictions"][model]
                confusion_rows.append(
                    confusion_frame(pred["y_test"], pred["knn_pred"], model, level, subset_name, "knn")
                )
                confusion_rows.append(
                    confusion_frame(pred["y_test"], pred["logreg_pred"], model, level, subset_name, "logreg")
                )
            spec = celltype_specific_purity(
                {k: np.asarray(v) for k, v in latents_sub.items()},
                labels_valid.reset_index(drop=True),
                np.arange(labels_valid.shape[0]),
                predictions=extras["predictions"],
                model_specs=PHASE7_MODELS,
            )
            spec["label_level"] = level
            spec["annotation_subset"] = subset_name
            specific_rows.append(spec)

            point = metrics.set_index("model")
            for a, b in COMPARE_PAIRS:
                boot = paired_bootstrap(
                    extras["per_cell"], seed=seed, n_boot=N_BOOT, model_a=a, model_b=b
                )
                f1_knn = classification_bootstrap_f1(
                    extras["predictions"][a]["y_test"],
                    extras["predictions"][a]["knn_pred"],
                    extras["predictions"][b]["knn_pred"],
                    seed=seed,
                    metric_name="knn_macro_f1",
                )
                f1_logreg = classification_bootstrap_f1(
                    extras["predictions"][a]["y_test"],
                    extras["predictions"][a]["logreg_pred"],
                    extras["predictions"][b]["logreg_pred"],
                    seed=seed,
                    metric_name="logreg_macro_f1",
                )
                for _, brow in pd.concat(
                    [boot, pd.DataFrame([f1_knn, f1_logreg])], ignore_index=True
                ).iterrows():
                    metric = brow["metric"]
                    primary_rows.append(
                        {
                            "comparison": f"{b} - {a}",
                            "model_a": a,
                            "model_b": b,
                            "metric": metric,
                            "label_level": level,
                            "annotation_subset": subset_name,
                            "model_a_value": float(point.loc[a, metric]) if metric in point.columns else np.nan,
                            "model_b_value": float(point.loc[b, metric]) if metric in point.columns else np.nan,
                            "delta_b_minus_a": float(brow["delta_mean"]),
                            "bootstrap_ci_low": float(brow["bootstrap_ci_low"]),
                            "bootstrap_ci_high": float(brow["bootstrap_ci_high"]),
                            "n_boot": int(brow["n_boot"]),
                        }
                    )

    metrics_all = pd.concat(metric_rows, ignore_index=True)
    metrics_all.to_csv(tables_dir / "mofa_biological_representation_metrics.csv", index=False)
    pd.concat(leiden_rows, ignore_index=True).to_csv(
        tables_dir / "mofa_leiden_resolution_grid_ari_nmi.csv", index=False
    )
    pd.concat(confusion_rows, ignore_index=True).to_csv(
        tables_dir / "mofa_celltype_confusion_matrices.csv", index=False
    )
    pd.concat(specific_rows, ignore_index=True).to_csv(
        tables_dir / "mofa_celltype_specific_representation_metrics.csv", index=False
    )
    primary = pd.DataFrame(primary_rows)
    primary.to_csv(tables_dir / "primary_mofa_vs_totalvi.csv", index=False)

    high = metrics_all[metrics_all["annotation_subset"] == "high_confidence"]
    five = ["PCA_RNA", "scVI_default", "scVI_matched", "MOFA+", "totalVI"]
    plot_metric_bars(
        high, "celltype_ASW", figures_dir / "phase7_figureC_asw_highconf",
        "Cell-type ASW (high-confidence)", models=five,
    )
    plot_metric_bars(
        high, "knn_purity_k15", figures_dir / "phase7_figureC_purity_k15_highconf",
        "kNN purity k=15 (high-confidence)", models=five,
    )
    plot_metric_bars(
        high, "logreg_macro_f1", figures_dir / "phase7_figureD_logreg_macro_f1_highconf",
        "Held-out logreg macro-F1 (high-confidence)", models=five,
    )

    provenance = logs_dir / "mofa_provenance.md"
    provenance.write_text(
        f"""# MOFA+ provenance (PHASE 7)

- Date (UTC): {datetime.now(timezone.utc).isoformat()}
- Model: MOFA+ (Argelaguet et al., Genome Biology 2020) via mofapy2
- This is a linear Bayesian factor model on processed RNA + protein views
- It is **not** a likelihood-matched or architecture-matched control for totalVI
- RNA view: 2000 seurat_v3 HVGs, library-size normalized, log1p, scaled
- Protein view: CLR (pseudocount 1) then per-protein z-score on original observed counts
- Groups: PBMC10k, PBMC5k (batch)
- Factors: 20 (matched to VAE n_latent)
- Gaussian likelihoods; scale_views=True; weight_views=True
- RNA PCA control uses the identical 2000-gene scaled RNA matrix and 20 PCs
- PHASE 3/5 VAE models were not loaded or retrained
- Original inner-join h5ad was not modified
- Labels: PHASE 6 independent Seurat v4 SCANVI transfer
- Sparsity stress tests and missing-modality experiments were not started
- Training notes: {train_info}
""",
        encoding="utf-8",
    )
    logger.info("PHASE 7 complete. Provenance: %s", provenance)
    return {
        "n_query": int(query.n_obs),
        "n_high_confidence": int(high_mask.sum()),
        "n_models": len(PHASE7_MODELS),
        "mofa_latent": str(embeddings_dir / f"mofa_seed{seed}_latent.npy"),
        "pca_latent": str(embeddings_dir / f"pca_rna_seed{seed}_latent.npy"),
    }
