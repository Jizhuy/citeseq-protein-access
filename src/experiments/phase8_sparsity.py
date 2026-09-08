"""PHASE 8: protein sparsity / corruption stress test.

Does not study true missing modality. Does not retrain scVI_matched.
model_seed is fixed at 0; perturbation_seed varies the corruption mask.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import precision_recall_fscore_support

from src.data.perturbations import (
    CORRUPTION_LEVELS,
    MODEL_SEED,
    PERTURBATION_SEEDS,
    SMOKE_CONDITIONS,
    apply_mask,
    corruption_code,
    create_or_load_mask,
    mask_path,
    summarize_corruption,
)
from src.evaluation.phase6 import (
    N_BOOT,
    PRIMARY_RESOLUTION,
    classification_bootstrap_f1,
    evaluate_subset,
    heldout_with_predictions,
    neighborhood_purity_per_cell,
)
from src.evaluation.phase8_plots import (
    plot_asw_purity,
    plot_celltype_degradation,
    plot_delta_f1,
    plot_masked_recovery,
    plot_primary_f1,
    plot_protein_recovery,
    plot_sparsity_vs_corruption,
)
from src.models.run_totalvi_corrupted import train_corrupted_totalvi
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

ANNOTATION_SOURCE = (
    "Seurat v4 PBMC CITE-seq (Stuart et al., Cell 2021) via PHASE 6 SCANVI/scArches; "
    "labels independent of corrupted totalVI latents"
)
PRIMARY_SUBSET = "high_confidence"
PRIMARY_LEVEL = "l2"
PRIMARY_METRIC = "logreg_macro_f1"
PROTEIN_KEY = "protein_expression"
SCVI_L2_REFERENCE_NOTE = 0.722


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if x.size < 3 or float(np.std(x)) < 1e-12 or float(np.std(y)) < 1e-12:
        return float("nan")
    fn = pearsonr if method == "pearson" else spearmanr
    return float(fn(x, y)[0])


def _rmse(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() == 0:
        return float("nan")
    return float(np.sqrt(np.mean((x[keep] - y[keep]) ** 2)))


def _tag(fraction: float, seed: int) -> str:
    return f"corruption_{corruption_code(fraction)}_seed{seed}"


def _condition_dir(models_dir: Path, fraction: float, seed: int) -> Path:
    return models_dir / "stress_sparsity" / f"corruption_{corruption_code(fraction)}" / f"seed{seed}"


def _read_protein_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = [str(c) for c in frame.columns]
    return frame


def evaluate_latent(
    latent: np.ndarray,
    labels: pd.Series,
    split: pd.Series,
    batch: pd.Series,
    *,
    model_name: str,
    label_level: str,
    subset_name: str,
    seed: int,
) -> tuple[pd.Series, pd.DataFrame, dict[str, Any]]:
    specs = ((model_name, "X", "RNA+Protein" if "totalVI" in model_name or "MOFA" in model_name else "RNA", "x"),)
    metrics, leiden, extras = evaluate_subset(
        {model_name: np.asarray(latent)},
        labels.reset_index(drop=True),
        split.reset_index(drop=True),
        batch.reset_index(drop=True),
        label_level=label_level,
        annotation_subset=subset_name,
        annotation_source=ANNOTATION_SOURCE,
        confidence_rule="PHASE 6 high-confidence rule: l2_confidence >= 0.85",
        seed=seed,
        model_specs=specs,
    )
    return metrics.iloc[0], leiden, extras


def protein_recovery_table(
    original: np.ndarray,
    mask: np.ndarray,
    estimate: pd.DataFrame,
    foreground: pd.DataFrame | None,
    protein_names: list[str],
    fraction: float,
    seed: int,
) -> pd.DataFrame:
    estimate = estimate.reindex(columns=protein_names)
    est = estimate.to_numpy(dtype=float)
    rows = []
    for j, name in enumerate(protein_names):
        masked = mask[:, j]
        all_entries = np.ones(original.shape[0], dtype=bool)
        orig_pos = original[:, j] > 0
        orig_zero = original[:, j] <= 0
        row = {
            "corruption_fraction": float(fraction),
            "seed": int(seed),
            "protein": name,
            "n_masked_entries": int(masked.sum()),
            "pearson_masked": _safe_corr(est[masked, j], original[masked, j], "pearson") if masked.any() else np.nan,
            "spearman_masked": _safe_corr(est[masked, j], original[masked, j], "spearman") if masked.any() else np.nan,
            "rmse_masked": _rmse(est[masked, j], original[masked, j]) if masked.any() else np.nan,
            "pearson_all_measured": _safe_corr(est[all_entries, j], original[all_entries, j], "pearson"),
            "spearman_all_measured": _safe_corr(est[all_entries, j], original[all_entries, j], "spearman"),
            "rmse_all_measured": _rmse(est[all_entries, j], original[all_entries, j]),
            "foreground_mean_masked": np.nan,
            "foreground_mean_original_positive_unmasked": np.nan,
            "foreground_mean_original_zero": np.nan,
            "foreground_recovery_metric_if_available": np.nan,
        }
        if foreground is not None and name in foreground.columns:
            fg = foreground[name].to_numpy(dtype=float)
            unmasked_pos = orig_pos & ~masked
            if masked.any():
                row["foreground_mean_masked"] = float(np.nanmean(fg[masked]))
                row["foreground_recovery_metric_if_available"] = float(np.nanmean(fg[masked]))
            if unmasked_pos.any():
                row["foreground_mean_original_positive_unmasked"] = float(np.nanmean(fg[unmasked_pos]))
            if orig_zero.any():
                row["foreground_mean_original_zero"] = float(np.nanmean(fg[orig_zero]))
        rows.append(row)
    return pd.DataFrame(rows)


def celltype_table(
    y_test: np.ndarray,
    pred: np.ndarray,
    purity: np.ndarray,
    y_all: np.ndarray,
    fraction: float,
    seed: int,
    label_level: str,
) -> pd.DataFrame:
    labels = sorted(set(y_all))
    precision, recall, f1, support = precision_recall_fscore_support(
        y_test, pred, labels=labels, zero_division=0
    )
    rows = []
    for i, cell_type in enumerate(labels):
        mask = y_all == cell_type
        rows.append(
            {
                "corruption_fraction": float(fraction),
                "seed": int(seed),
                "label_level": label_level,
                "cell_type": cell_type,
                "support": int((y_test == cell_type).sum()),
                "n_cells": int(mask.sum()),
                "precision": float(precision[i]),
                "recall": float(recall[i]),
                "f1": float(f1[i]),
                "knn_purity": float(purity[mask].mean()) if mask.any() else np.nan,
                "too_few_test": bool((y_test == cell_type).sum() < 10),
            }
        )
    return pd.DataFrame(rows)


def reuse_phase4_totalvi(embeddings_dir: Path, obs_names: pd.Index) -> dict[str, Any]:
    latent = np.load(embeddings_dir / "totalvi_seed0_latent.npy")
    cells = pd.read_csv(embeddings_dir / "totalvi_seed0_cells.csv")
    if list(cells["cell_id"].astype(str)) != list(obs_names.astype(str)):
        raise RuntimeError("PHASE 4 totalVI cell order does not match the query.")
    if latent.shape != (len(obs_names), 20) or not np.isfinite(latent).all():
        raise RuntimeError("PHASE 4 totalVI latent is unusable.")
    protein = _read_protein_csv(embeddings_dir / "totalvi_seed0_normalized_protein.csv")
    foreground = _read_protein_csv(embeddings_dir / "totalvi_seed0_protein_foreground_probability.csv")
    protein = protein.reindex(obs_names.astype(str))
    foreground = foreground.reindex(obs_names.astype(str))
    manifest_path = resolve_path("results/logs/totalvi_seed0_manifest.json")
    import json

    with manifest_path.open("r", encoding="utf-8") as handle:
        saved = json.load(handle)
    return {
        "model": "totalVI",
        "reused_phase4_model": True,
        "runtime_seconds": saved.get("runtime_seconds"),
        "epochs": saved.get("epochs_completed"),
        "device": saved.get("actual_device"),
        "peak_rss_mb": saved.get("peak_rss_mb"),
        "early_stopped": saved.get("early_stopped"),
        "model_seed": 0,
        "perturbation_seed": 0,
        "corruption_fraction": 0.0,
        "latent": latent,
        "normalized_protein": protein,
        "foreground": foreground,
        "manifest_path": str(manifest_path),
        "note": "0% seed 0 reuses the PHASE 4 seed-0 totalVI model; there is no corruption randomness.",
    }


def apply_corruption_to_copy(
    query: ad.AnnData,
    original_protein: np.ndarray,
    mask: np.ndarray,
    protein_names: list[str],
) -> ad.AnnData:
    adata = query.copy()
    corrupted = apply_mask(original_protein, mask)
    frame = pd.DataFrame(corrupted, index=query.obs_names, columns=protein_names)
    adata.obsm["protein_expression"] = frame.copy()
    adata.obsm["protein_counts"] = frame.copy()
    return adata


def estimate_breakpoint(summary: pd.DataFrame) -> dict[str, Any]:
    ordered = summary.sort_values("corruption_fraction")
    levels = ordered["corruption_fraction"].to_numpy()
    deltas = ordered["delta_mean_across_masks"].to_numpy()
    crossing = None
    for i in range(len(levels) - 1):
        if deltas[i] > 0 and deltas[i + 1] <= 0:
            crossing = (float(levels[i]), float(levels[i + 1]))
            break
    ci_includes_zero_from = None
    if "ci_includes_zero_seed0" in ordered.columns:
        for level, flag in zip(levels, ordered["ci_includes_zero_seed0"]):
            if bool(flag):
                ci_includes_zero_from = float(level)
                break
    return {
        "mean_delta_crossing": crossing,
        "seed0_ci_includes_zero_from": ci_includes_zero_from,
        "all_mean_deltas_positive": bool(np.all(deltas > 0)),
        "all_mean_deltas_nonpositive": bool(np.all(deltas <= 0)),
        "language": (
            f"the advantage was no longer detectable between {crossing[0]:.0%} and {crossing[1]:.0%} corruption"
            if crossing
            else (
                "mean Delta F1 remained positive at all evaluated corruption levels"
                if np.all(deltas > 0)
                else "mean Delta F1 was not positive at the evaluated levels"
            )
        ),
    }


def run_phase8(
    config: dict[str, Any] | None = None,
    smoke: bool = True,
    skip_train: bool = False,
) -> dict[str, Any]:
    cfg = config or load_config()
    set_global_seed(MODEL_SEED)
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    models_dir = resolve_path(cfg["paths"].get("models_dir", "results/models"))
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    perturbations_dir = resolve_path("results/perturbations")
    stress_embed = embeddings_dir / "stress_sparsity"
    logger = get_logger("phase8", logs_dir / "phase8.log")
    perturbations_dir.mkdir(parents=True, exist_ok=True)
    stress_embed.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    source_path = processed_dir / "pbmc_cite_combined_inner.h5ad"
    query = ad.read_h5ad(source_path)
    annotated = ad.read_h5ad(processed_dir / "pbmc_cite_combined_inner_annotated.h5ad")
    if list(annotated.obs_names.astype(str)) != list(query.obs_names.astype(str)):
        raise RuntimeError("Annotated object cell order does not match the inner-join query.")
    if "cell_type_l1" in query.obs:
        raise RuntimeError("Original inner h5ad contains annotation columns.")
    protein = query.obsm[PROTEIN_KEY]
    if not isinstance(protein, pd.DataFrame):
        raise TypeError("protein_expression must be a DataFrame of original observed counts.")
    protein_names = [str(c) for c in protein.columns]
    original_protein = protein.to_numpy(dtype=np.float64)
    rna_checksum = float(np.asarray(query.layers["counts"]).sum())
    protein_checksum = float(original_protein.sum())
    high_mask = annotated.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"

    scvi_latent = np.load(embeddings_dir / "scvi_matched_seed0_latent.npy")
    scvi_cells = pd.read_csv(embeddings_dir / "scvi_matched_seed0_cells.csv")
    if list(scvi_cells["cell_id"].astype(str)) != list(query.obs_names.astype(str)):
        raise RuntimeError("scVI_matched cell order mismatch.")
    pca_path = embeddings_dir / "pca_rna_seed0_latent.npy"
    pca_latent = np.load(pca_path) if pca_path.exists() else None

    conditions = list(SMOKE_CONDITIONS) if smoke else _full_conditions()
    logger.info(
        "PHASE 8 %s: %s conditions; model_seed=%s; RNA unchanged; not missing-modality",
        "smoke" if smoke else "full grid",
        len(conditions),
        MODEL_SEED,
    )

    metric_rows = []
    leiden_rows = []
    delta_rows = []
    celltype_rows = []
    recovery_rows = []
    sparsity_rows = []
    failed = []
    smoke_checks = []
    ref_cache: dict[str, Any] = {}

    for level, label_col in (("l1", "cell_type_l1"), ("l2", "cell_type_l2")):
        labels = annotated.obs.loc[high_mask, label_col].astype(str)
        split = query.obs.loc[high_mask, "split"].astype(str)
        batch = query.obs.loc[high_mask, "batch"].astype(str)
        scvi_row, _, scvi_extra = evaluate_latent(
            scvi_latent[high_mask], labels, split, batch,
            model_name="scVI_matched", label_level=level, subset_name=PRIMARY_SUBSET, seed=MODEL_SEED,
        )
        ref_cache[level] = {"row": scvi_row, "extra": scvi_extra, "labels": labels, "split": split, "batch": batch}
        if pca_latent is not None:
            pca_row, _, _ = evaluate_latent(
                pca_latent[high_mask], labels, split, batch,
                model_name="PCA_RNA", label_level=level, subset_name=PRIMARY_SUBSET, seed=MODEL_SEED,
            )
            ref_cache[f"pca_{level}"] = pca_row

    for fraction, p_seed in conditions:
        tag = _tag(fraction, p_seed)
        try:
            payload = create_or_load_mask(
                original_protein, fraction, p_seed,
                mask_path(perturbations_dir, fraction, p_seed),
                query.obs_names, protein_names,
            )
            mask = payload["mask"]
            corrupted = apply_mask(original_protein, mask)
            summary = summarize_corruption(
                original_protein, corrupted, mask, protein_names, fraction, p_seed
            )
            sparsity_rows.append(
                {
                    **{k: v for k, v in summary.items() if not isinstance(v, (dict, list))},
                    "corruption_fraction": float(fraction),
                }
            )
            per_prot = []
            for name in protein_names:
                orig = summary["original_per_protein"][name]
                fin = summary["final_per_protein"][name]
                per_prot.append(
                    {
                        "corruption_fraction": fraction,
                        "seed": p_seed,
                        "protein": name,
                        "original_detection_fraction": orig["detection_fraction"],
                        "corrupted_detection_fraction": fin["detection_fraction"],
                        "original_mean_count": orig["mean_count"],
                        "original_median_count": orig["median_count"],
                        "corrupted_mean_count": fin["mean_count"],
                        "corrupted_median_count": fin["median_count"],
                    }
                )
            pd.DataFrame(per_prot).to_csv(
                tables_dir / f"protein_sparsity_per_protein_{tag}.csv", index=False
            )

            cond_dir = _condition_dir(models_dir, fraction, p_seed)
            latent_path = stress_embed / f"{tag}_latent.npy"
            protein_path = stress_embed / f"{tag}_normalized_protein.csv"
            fg_path = stress_embed / f"{tag}_protein_foreground_probability.csv"
            if fraction == 0.0 and p_seed == 0:
                reused = reuse_phase4_totalvi(embeddings_dir, query.obs_names)
                latent = reused["latent"]
                protein_hat = reused["normalized_protein"]
                foreground = reused["foreground"]
                train_info = reused
                np.save(latent_path, latent)
                pd.DataFrame(
                    {"cell_id": query.obs_names.astype(str), "batch": query.obs["batch"].astype(str),
                     "split": query.obs["split"].astype(str)}
                ).to_csv(stress_embed / f"{tag}_cells.csv", index=False)
            elif latent_path.exists() and protein_path.exists():
                logger.info("Reusing existing embedding %s", latent_path)
                latent = np.load(latent_path)
                protein_hat = _read_protein_csv(protein_path)
                foreground = _read_protein_csv(fg_path) if fg_path.exists() else None
                manifest = cond_dir / "manifest.json"
                train_info = {"reused_existing_embedding": True, "runtime_seconds": np.nan, "epochs": np.nan, "device": None, "peak_rss_mb": np.nan, "early_stopped": np.nan}
                if manifest.exists():
                    import json
                    with manifest.open("r", encoding="utf-8") as handle:
                        saved = json.load(handle)
                    train_info = {**train_info, **{k: saved.get(k) for k in ("runtime_seconds", "epochs", "device", "peak_rss_mb", "early_stopped")}}
            elif skip_train:
                raise FileNotFoundError(f"Missing embedding for {tag} and --skip-train was set.")
            else:
                adata_c = apply_corruption_to_copy(query, original_protein, mask, protein_names)
                train_info = train_corrupted_totalvi(
                    adata_c,
                    condition_dir=cond_dir,
                    embeddings_dir=stress_embed,
                    fraction=fraction,
                    perturbation_seed=p_seed,
                    model_seed=MODEL_SEED,
                    config=cfg,
                    logger=logger,
                )
                latent = np.load(latent_path)
                protein_hat = _read_protein_csv(stress_embed / f"{tag}_normalized_protein.csv")
                foreground = _read_protein_csv(stress_embed / f"{tag}_protein_foreground_probability.csv")

            if latent.shape[0] != query.n_obs or not np.isfinite(latent).all():
                raise RuntimeError(f"{tag} latent is missing or non-finite.")
            source_check = ad.read_h5ad(source_path)
            if abs(float(np.asarray(source_check.layers["counts"]).sum()) - rna_checksum) > 1e-3:
                raise RuntimeError("Original RNA counts were modified.")
            if abs(float(source_check.obsm[PROTEIN_KEY].to_numpy().sum()) - protein_checksum) > 1e-3:
                raise RuntimeError("Original protein counts were modified.")
            if "X_totalVI" in source_check.obsm or "cell_type_l1" in source_check.obs:
                raise RuntimeError("Original inner h5ad was modified.")

            rec = protein_recovery_table(
                original_protein, mask, protein_hat, foreground, protein_names, fraction, p_seed
            )
            recovery_rows.append(rec)

            for level in ("l1", "l2"):
                labels = ref_cache[level]["labels"]
                split = ref_cache[level]["split"]
                batch = ref_cache[level]["batch"]
                row, leiden, extras = evaluate_latent(
                    latent[high_mask], labels, split, batch,
                    model_name="totalVI_corrupted", label_level=level,
                    subset_name=PRIMARY_SUBSET, seed=MODEL_SEED,
                )
                row = row.to_dict()
                row.update(
                    {
                        "model": "totalVI",
                        "corruption_fraction": float(fraction),
                        "final_protein_sparsity": summary["final_protein_sparsity"],
                        "perturbation_seed": int(p_seed),
                        "model_seed": int(MODEL_SEED),
                        "runtime_seconds": train_info.get("runtime_seconds"),
                        "epochs": train_info.get("epochs"),
                        "device": train_info.get("device"),
                        "peak_rss_mb": train_info.get("peak_rss_mb"),
                        "early_stopped": train_info.get("early_stopped"),
                        "reused_phase4_model": bool(train_info.get("reused_phase4_model", False)),
                        "ASW": row["celltype_ASW"],
                        "knn_purity_15": row["knn_purity_k15"],
                        "knn_purity_30": row["knn_purity_k30"],
                        "knn_purity_50": row["knn_purity_k50"],
                    }
                )
                metric_rows.append(row)
                leiden["corruption_fraction"] = fraction
                leiden["perturbation_seed"] = p_seed
                leiden_rows.append(leiden)

                pred = extras["predictions"]["totalVI_corrupted"]
                scvi_pred = ref_cache[level]["extra"]["predictions"]["scVI_matched"]
                boot = classification_bootstrap_f1(
                    pred["y_test"], scvi_pred["logreg_pred"], pred["logreg_pred"],
                    seed=MODEL_SEED, n_boot=N_BOOT, metric_name="logreg_macro_f1",
                )
                if level == PRIMARY_LEVEL:
                    delta_rows.append(
                        {
                            "corruption_fraction": float(fraction),
                            "perturbation_seed": int(p_seed),
                            "model_seed": int(MODEL_SEED),
                            "totalVI_l2_macro_f1": float(row["logreg_macro_f1"]),
                            "scVI_matched_l2_macro_f1": float(ref_cache[level]["row"]["logreg_macro_f1"]),
                            "delta": float(boot["delta_mean"]),
                            "bootstrap_ci_low": float(boot["bootstrap_ci_low"]),
                            "bootstrap_ci_high": float(boot["bootstrap_ci_high"]),
                            "n_boot": int(boot["n_boot"]),
                            "point_delta": float(row["logreg_macro_f1"] - ref_cache[level]["row"]["logreg_macro_f1"]),
                        }
                    )
                purity = extras["per_cell"]["totalVI_corrupted"]["purity"][15]
                celltype_rows.append(
                    celltype_table(
                        pred["y_test"], pred["logreg_pred"], purity,
                        labels.reset_index(drop=True).to_numpy(),
                        fraction, p_seed, level,
                    )
                )

            smoke_checks.append(
                {
                    "condition": tag,
                    "sparsity_ok": True,
                    "latent_finite": True,
                    "metrics_ok": True,
                    "original_unchanged": True,
                    "masked_recovery_ok": bool(fraction == 0.0 or rec["pearson_masked"].notna().any()),
                    "realized_mask_fraction": summary["realized_nonzero_masked_fraction"],
                    "final_sparsity": summary["final_protein_sparsity"],
                    "l2_macro_f1": next(
                        r["logreg_macro_f1"] for r in metric_rows
                        if r["corruption_fraction"] == fraction and r["perturbation_seed"] == p_seed and r["label_level"] == "l2"
                    ),
                }
            )
            logger.info("Finished %s l2_macro_F1=%.4f", tag, smoke_checks[-1]["l2_macro_f1"])
        except Exception as exc:  # noqa: BLE001
            logger.exception("Failed %s", tag)
            failed.append({"condition": tag, "error": f"{type(exc).__name__}: {exc}"})
            smoke_checks.append({"condition": tag, "error": str(exc), "metrics_ok": False})

    metrics_all = pd.DataFrame(metric_rows)
    if len(metrics_all):
        metrics_all.to_csv(tables_dir / "protein_sparsity_robustness.csv", index=False)
    if leiden_rows:
        pd.concat(leiden_rows, ignore_index=True).to_csv(
            tables_dir / "protein_sparsity_leiden_grid.csv", index=False
        )
    delta = pd.DataFrame(delta_rows)
    if len(delta):
        delta.to_csv(tables_dir / "protein_sparsity_primary_delta.csv", index=False)
        grouped = delta.groupby("corruption_fraction").agg(
            n_masks=("perturbation_seed", "nunique"),
            totalVI_l2_macro_f1_mean=("totalVI_l2_macro_f1", "mean"),
            totalVI_l2_macro_f1_sd=("totalVI_l2_macro_f1", "std"),
            delta_mean_across_masks=("point_delta", "mean"),
            delta_sd_across_masks=("point_delta", "std"),
        ).reset_index()
        seed0 = delta[delta["perturbation_seed"] == 0][
            ["corruption_fraction", "bootstrap_ci_low", "bootstrap_ci_high", "delta"]
        ].rename(columns={"delta": "seed0_bootstrap_delta"})
        grouped = grouped.merge(seed0, on="corruption_fraction", how="left")
        grouped["ci_includes_zero_seed0"] = (
            (grouped["bootstrap_ci_low"] <= 0) & (grouped["bootstrap_ci_high"] >= 0)
        )
        grouped["scVI_matched_l2_macro_f1"] = float(ref_cache["l2"]["row"]["logreg_macro_f1"])
        grouped.to_csv(tables_dir / "protein_sparsity_primary_delta_summary.csv", index=False)
    else:
        grouped = pd.DataFrame()
    if celltype_rows:
        pd.concat(celltype_rows, ignore_index=True).to_csv(
            tables_dir / "protein_sparsity_celltype_specific.csv", index=False
        )
    if recovery_rows:
        rec_all = pd.concat(recovery_rows, ignore_index=True)
        rec_all.to_csv(tables_dir / "protein_corruption_recovery.csv", index=False)
    else:
        rec_all = pd.DataFrame()
    sparsity_table = pd.DataFrame(sparsity_rows)
    if len(sparsity_table):
        sparsity_table.to_csv(tables_dir / "protein_sparsity_realized_levels.csv", index=False)

    scvi_l2 = ref_cache["l2"]["row"]
    pca_l2 = float(ref_cache["pca_l2"]["logreg_macro_f1"]) if "pca_l2" in ref_cache else None
    l2_metrics = metrics_all[metrics_all["label_level"] == "l2"] if len(metrics_all) else metrics_all
    if len(sparsity_table):
        plot_sparsity_vs_corruption(sparsity_table, figures_dir / "phase8_figure1_sparsity_vs_corruption")
    if len(l2_metrics):
        plot_primary_f1(l2_metrics, float(scvi_l2["logreg_macro_f1"]), pca_l2, figures_dir / "phase8_figure2_l2_macro_f1")
        plot_asw_purity(l2_metrics, scvi_l2, figures_dir / "phase8_figure4_asw_purity")
    if len(delta):
        plot_delta_f1(delta, figures_dir / "phase8_figure3_delta_macro_f1")
    if len(rec_all):
        plot_masked_recovery(rec_all, figures_dir / "phase8_figure5_masked_protein_recovery")
        plot_protein_recovery(rec_all, figures_dir / "phase8_figure7_per_protein_recovery")
    if celltype_rows:
        ct = pd.concat(celltype_rows, ignore_index=True)
        plot_celltype_degradation(
            ct[ct["label_level"] == "l2"], figures_dir / "phase8_figure6_celltype_f1"
        )

    breakpoint = estimate_breakpoint(grouped) if len(grouped) else {}
    provenance = logs_dir / "phase8_provenance.md"
    provenance.write_text(
        f"""# PHASE 8 protein sparsity / corruption stress test

- Date (UTC): {datetime.now(timezone.utc).isoformat()}
- This is a **protein sparsity / corruption** experiment, NOT a missing-protein-modality experiment
- Masked nonzero protein counts are set to 0 and still presented to totalVI as observations
- RNA, labels, batch, and the PHASE 2 train/test split were not perturbed
- Primary reference: fixed scVI_matched latent (not retrained)
- model_seed = {MODEL_SEED} for all totalVI runs; perturbation_seed varies the corruption mask
- 0% seed 0 reuses the PHASE 4 totalVI model (no fake perturbation replicates)
- Primary metric: held-out logreg macro-F1 on celltype.l2, high-confidence cells
- PCA_RNA is a fixed secondary RNA reference; MOFA+ was not rerun in this PHASE
- Failed conditions: {failed}
""",
        encoding="utf-8",
    )
    result = {
        "mode": "smoke" if smoke else "full",
        "n_conditions_attempted": len(conditions),
        "n_successful": int(len(conditions) - len(failed)),
        "n_failed": len(failed),
        "failed": failed,
        "smoke_checks": smoke_checks,
        "scVI_matched_l2_macro_f1": float(scvi_l2["logreg_macro_f1"]),
        "breakpoint": breakpoint,
        "original_rna_checksum": rna_checksum,
        "original_protein_checksum": protein_checksum,
        "model_seed": MODEL_SEED,
        "primary_metric": PRIMARY_METRIC,
        "primary_label_level": PRIMARY_LEVEL,
        "primary_subset": PRIMARY_SUBSET,
    }
    save_json(result, logs_dir / ("phase8_smoke_summary.json" if smoke else "phase8_full_summary.json"))
    logger.info("PHASE 8 %s complete. failed=%s", "smoke" if smoke else "full", failed)
    return result


def _full_conditions() -> list[tuple[float, int]]:
    out = [(0.0, 0)]
    for fraction in CORRUPTION_LEVELS:
        if fraction == 0.0:
            continue
        for seed in PERTURBATION_SEEDS:
            out.append((float(fraction), int(seed)))
    return out
