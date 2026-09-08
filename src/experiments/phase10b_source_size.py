#!/usr/bin/env python
"""PHASE 10B: source-size sensitivity control for PHASE 10 Direction A.

Subsamples the PBMC10k source from 6,855 to 3,994 cells (the PBMC5k source
size) and repeats the PBMC10k -> PBMC5k transfer. Reuses the PHASE 10 training,
scArches query and evaluation code unchanged so the two experiments stay
comparable. Direction B is not retrained; its PHASE 10 embeddings are reused.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

from src.evaluation.phase10b_plots import write_phase10b_figures
from src.evaluation.phase6 import neighborhood_purity_per_cell
from src.experiments.phase10_cross_dataset import (
    DEFAULT_HYPERPARAMETERS,
    KNN_KS,
    MATCHED_HYPERPARAMETERS,
    MIN_SOURCE_SUPPORT,
    MIN_TARGET_SUPPORT,
    N_BOOT,
    _cleanup,
    _latent_drift,
    _ram,
    _require_cuda,
    _save_latent,
    class_plan,
    classify_transfer,
    ensure_csr,
    load_source_target,
    neighbor_agreement,
    paired_bootstrap_delta,
    per_class_table,
    query_adapt,
    train_scvi_source,
    train_totalvi_source,
)
from src.utils.device import detect_compute_environment
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

SOURCE_BATCH = "PBMC10k"
TARGET_BATCH = "PBMC5k"
DIRECTION_A = "pbmc10k_to_pbmc5k"
DIRECTION_B = "pbmc5k_to_pbmc10k"
SUBSAMPLE_N = 3994
SUBSAMPLE_SEEDS = (0, 1, 2, 3, 4)
MODEL_SEED = 0
PHASE10_SEEDS = (0, 1, 2, 3, 4)
LEVELS = (("l2", "cell_type_l2"), ("l1", "cell_type_l1"))


def _root() -> Path:
    return resolve_path("results/phase10b_source_size")


def _paths() -> dict[str, Path]:
    root = _root()
    return {
        name: root / name
        for name in ("models", "embeddings", "tables", "figures", "logs", "subsamples")
    }


# --------------------------------------------------------------------------
# subsampling
# --------------------------------------------------------------------------


def subsample_index(n_source: int, n_keep: int, subsample_seed: int) -> np.ndarray:
    """Uniform random source subsample without replacement.

    The stream is keyed only on ``subsample_seed`` and is independent of the
    model seed, which is held at 0 throughout PHASE 10B. No target data and no
    labels of any kind enter this draw.
    """
    rng = np.random.default_rng(subsample_seed)
    return np.sort(rng.choice(n_source, size=n_keep, replace=False))


def composition_rows(
    full_labels: np.ndarray, sub_labels: np.ndarray, level: str, subsample_seed: int
) -> tuple[pd.DataFrame, float]:
    full_counts = pd.Series(full_labels).value_counts()
    sub_counts = pd.Series(sub_labels).value_counts()
    cats = sorted(set(full_counts.index) | set(sub_counts.index))
    fc = full_counts.reindex(cats, fill_value=0)
    sc = sub_counts.reindex(cats, fill_value=0)
    ff = fc / max(len(full_labels), 1)
    sf = sc / max(len(sub_labels), 1)
    frame = pd.DataFrame(
        {
            "subsample_seed": subsample_seed,
            "label_level": level,
            "cell_type": cats,
            "full_source_count": fc.to_numpy(),
            "subsample_count": sc.to_numpy(),
            "full_source_fraction": ff.to_numpy(),
            "subsample_fraction": sf.to_numpy(),
            "absolute_fraction_difference": np.abs(ff.to_numpy() - sf.to_numpy()),
        }
    )
    tvd = float(0.5 * np.abs(ff.to_numpy() - sf.to_numpy()).sum())
    return frame, tvd


# --------------------------------------------------------------------------
# fixed-class evaluation helpers
# --------------------------------------------------------------------------


def eval_index(adata: ad.AnnData, label_col: str, classes: list[str]) -> tuple[np.ndarray, int]:
    """High-confidence cells whose label is in ``classes`` (PHASE 6 rule)."""
    high = adata.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    y = adata.obs[label_col].astype(str).to_numpy()
    keep = high & np.isin(y, list(classes))
    return np.where(keep)[0], int(high.sum())


def _metrics_row(clf: dict[str, Any], prefix: str) -> dict[str, float]:
    return {
        "accuracy": clf[f"{prefix}_accuracy"],
        "balanced_accuracy": clf[f"{prefix}_balanced_accuracy"],
        "macro_f1": clf[f"{prefix}_macro_f1"],
        "weighted_f1": clf[f"{prefix}_weighted_f1"],
    }


def _agreement_cols(z_s, y_s, z_t, y_t) -> dict[str, float]:
    out = {}
    for k in KNN_KS:
        out[f"neighbor_agreement_k{k}"] = neighbor_agreement(z_s, y_s, z_t, y_t, k)["mean_agreement"]
    return out


def _target_only_quality(z_t_all: np.ndarray, target: ad.AnnData, col: str) -> dict[str, float]:
    high = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    y = target.obs[col].astype(str).to_numpy()
    keep = high & ~pd.isna(target.obs[col]).to_numpy()
    if keep.sum() <= 20 or pd.Series(y[keep]).nunique() <= 1:
        return {}
    return {
        "target_asw": float(silhouette_score(z_t_all[keep], y[keep], metric="euclidean")),
        "target_knn_purity_15": float(neighborhood_purity_per_cell(z_t_all[keep], y[keep], 15).mean()),
    }


# --------------------------------------------------------------------------
# reference conditions reused from PHASE 10 (no retraining)
# --------------------------------------------------------------------------


def _phase10_emb(direction: str) -> Path:
    return resolve_path("results/phase10_cross_dataset/embeddings") / direction


def load_phase10_latents(
    direction: str, model: str, seed: int, src_adata: ad.AnnData, tgt_adata: ad.AnnData
) -> tuple[np.ndarray, np.ndarray]:
    emb = _phase10_emb(direction)
    tag = f"{model}_{direction}_seed{seed}_cuda"
    z_s = np.load(emb / f"{tag}_source_latent.npy")
    z_t = np.load(emb / f"{tag}_target_latent.npy")
    src_cells = pd.read_csv(emb / f"{tag}_source_cells.csv")["cell_id"].astype(str).tolist()
    tgt_cells = pd.read_csv(emb / f"{tag}_target_cells.csv")["cell_id"].astype(str).tolist()
    if src_cells != list(src_adata.obs_names.astype(str)):
        raise RuntimeError(f"{tag}: PHASE 10 source cell order does not match reloaded AnnData.")
    if tgt_cells != list(tgt_adata.obs_names.astype(str)):
        raise RuntimeError(f"{tag}: PHASE 10 target cell order does not match reloaded AnnData.")
    return z_s, z_t


def reference_condition(
    direction: str,
    condition: str,
    src_adata: ad.AnnData,
    tgt_adata: ad.AnnData,
    class_sets: dict[str, list[str]],
    logger,
    bootstrap_levels: tuple[str, ...] = ("l2",),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Re-score stored PHASE 10 embeddings on the common class set."""
    prim_rows, delta_rows, cell_rows = [], [], []
    for level, col in LEVELS:
        classes = class_sets[level]
        si, _ = eval_index(src_adata, col, classes)
        ti, n_high_t = eval_index(tgt_adata, col, classes)
        y_s = src_adata.obs[col].astype(str).to_numpy()[si]
        y_t = tgt_adata.obs[col].astype(str).to_numpy()[ti]
        src_counts = pd.Series(y_s).value_counts()
        tgt_counts = pd.Series(y_t).value_counts()
        for seed in PHASE10_SEEDS:
            preds = {}
            for model in ("scvi_matched", "totalvi"):
                z_s, z_t = load_phase10_latents(direction, model, seed, src_adata, tgt_adata)
                clf = classify_transfer(z_s[si], y_s, z_t[ti], y_t, MODEL_SEED)
                preds[model] = clf
                for classifier, prefix in (("logreg", "logreg"), ("knn15", "knn15"), ("knn30", "knn30")):
                    row = {
                        "condition": condition,
                        "direction": direction,
                        "subsample_seed": np.nan,
                        "model_seed": seed,
                        "model": model,
                        "source_n": int(src_adata.n_obs),
                        "target_n": int(tgt_adata.n_obs),
                        "label_level": level,
                        "evaluation_class_set": "common",
                        "n_classes": len(classes),
                        "target_coverage": float(len(ti) / n_high_t) if n_high_t else np.nan,
                        "classifier": classifier,
                        **_metrics_row(clf, prefix),
                    }
                    if classifier == "logreg":
                        row.update(_agreement_cols(z_s[si], y_s, z_t[ti], y_t))
                    prim_rows.append(row)
                if level == "l2":
                    spec = per_class_table(y_t, clf["logreg_pred"], classes, src_counts, tgt_counts)
                    spec["condition"] = condition
                    spec["direction"] = direction
                    spec["model"] = model
                    spec["model_seed"] = seed
                    spec["subsample_seed"] = np.nan
                    cell_rows.append(spec)
                del z_s, z_t
            if level in bootstrap_levels:
                boot = paired_bootstrap_delta(
                    y_t, preds["scvi_matched"]["logreg_pred"], preds["totalvi"]["logreg_pred"],
                    seed=seed, n_boot=N_BOOT,
                )
            else:
                boot = {"bootstrap_ci_low": np.nan, "bootstrap_ci_high": np.nan}
            delta_rows.append(
                {
                    "condition": condition,
                    "direction": direction,
                    "subsample_seed": np.nan,
                    "model_seed": seed,
                    "label_level": level,
                    "evaluation_class_set": "common",
                    "scvi_macro_f1": preds["scvi_matched"]["logreg_macro_f1"],
                    "totalvi_macro_f1": preds["totalvi"]["logreg_macro_f1"],
                    "delta": float(
                        preds["totalvi"]["logreg_macro_f1"] - preds["scvi_matched"]["logreg_macro_f1"]
                    ),
                    "bootstrap_ci_low": boot["bootstrap_ci_low"],
                    "bootstrap_ci_high": boot["bootstrap_ci_high"],
                    "target_n_eval": int(len(ti)),
                    "n_classes": len(classes),
                }
            )
        logger.info("%s %s re-scored on %s common classes", condition, level, len(classes))
    cell = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    return pd.DataFrame(prim_rows), pd.DataFrame(delta_rows), cell


# --------------------------------------------------------------------------
# size-matched training
# --------------------------------------------------------------------------


def train_subsample(
    subsample_seed: int,
    source_sub: ad.AnnData,
    target: ad.AnnData,
    cfg: dict[str, Any],
    logger,
) -> dict[str, dict[str, Any]]:
    from scvi.model import SCVI, TOTALVI

    paths = _paths()
    max_epochs = int(cfg["models"]["max_epochs"])
    batch_size = int(cfg["models"]["batch_size"])
    query_epochs = int(cfg.get("phase10", {}).get("query_max_epochs", 200))
    fitted: dict[str, dict[str, Any]] = {}

    for model_name in ("scvi_matched", "totalvi"):
        tag = (
            f"{model_name}_pbmc10k_sub{SUBSAMPLE_N}"
            f"_subsample_seed{subsample_seed}_modelseed{MODEL_SEED}_cuda"
        )
        model_dir = (
            paths["models"] / f"subsample_seed{subsample_seed}" / model_name / f"modelseed{MODEL_SEED}_cuda"
        )
        latent_src_path = paths["embeddings"] / f"{tag}_source_latent.npy"
        latent_tgt_path = paths["embeddings"] / f"{tag}_target_latent.npy"
        manifest_path = paths["logs"] / f"{tag}_manifest.json"

        if latent_src_path.exists() and latent_tgt_path.exists() and manifest_path.exists():
            logger.info("Reusing completed PHASE 10B artifacts for %s", tag)
            import json

            fitted[model_name] = {
                "source_latent": np.load(latent_src_path),
                "target_latent": np.load(latent_tgt_path),
                "manifest": json.loads(manifest_path.read_text(encoding="utf-8")),
                "tag": tag,
            }
            continue

        set_global_seed(MODEL_SEED)
        if model_name == "scvi_matched":
            hparams = dict(MATCHED_HYPERPARAMETERS)
            hparams["max_epochs"] = max_epochs
            hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
            hparams["batch_size"] = batch_size
            src = train_scvi_source(source_sub, hparams, model_dir, logger)
            q = query_adapt(SCVI, target, source_sub, model_dir, query_epochs, batch_size, None, logger)
        else:
            hparams = dict(DEFAULT_HYPERPARAMETERS)
            hparams["max_epochs"] = max_epochs
            hparams["batch_size"] = batch_size
            src = train_totalvi_source(source_sub, hparams, model_dir, logger)
            q = query_adapt(
                TOTALVI,
                target,
                source_sub,
                model_dir,
                query_epochs,
                batch_size,
                {"lr": DEFAULT_HYPERPARAMETERS["lr"], "reduce_lr_on_plateau": True},
                logger,
            )

        source_latent = q["source_latent_adapted"]
        if source_latent.shape[1] != q["latent"].shape[1]:
            raise RuntimeError("Adapted source and query latent dimensions differ.")
        drift = _latent_drift(src["latent"], source_latent)
        _save_latent(latent_src_path, source_latent, source_sub)
        _save_latent(latent_tgt_path, q["latent"], target)
        np.save(paths["embeddings"] / f"{tag}_source_latent_source_model.npy", src["latent"])

        manifest = {
            "phase": "10B",
            "model": model_name,
            "direction": DIRECTION_A,
            "condition": "size_matched_A",
            "subsample_seed": subsample_seed,
            "model_seed": MODEL_SEED,
            "source_n": int(source_sub.n_obs),
            "target_n": int(target.n_obs),
            "source_training_epochs": src["epochs"],
            "query_adaptation_epochs": q["epochs"],
            "source_training_runtime": src["runtime_seconds"],
            "query_runtime": q["runtime_seconds"],
            "source_batch_size": src["batch_size"],
            "query_batch_size": q["batch_size"],
            "source_device": src["device"],
            "query_device": q["device"],
            "source_gpu": src["gpu"],
            "query_gpu": q["gpu"],
            "source_ram_before": src["ram_before"],
            "source_ram_after": src["ram_after"],
            "query_ram_before": q["ram_before"],
            "query_ram_after": q["ram_after"],
            "n_trainable_query_params": q["n_trainable_params"],
            "trainable_query_params": q["trainable_param_names"],
            "source_embedding_used_for_metrics": "adapted_model",
            "adapted_source_drift_vs_source_model": drift,
            "query_api": "prepare_query_anndata + load_query_data",
            "zero_shot": False,
            "target_labels_used_in_fitting": False,
            "phase10_artifact_overwritten": False,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        save_json(manifest, manifest_path)
        logger.info(
            "%s drift=%.4f src_ep=%s q_ep=%s src_s=%.1f q_s=%.1f",
            tag,
            drift["mean_euclidean_shift"],
            src["epochs"],
            q["epochs"],
            src["runtime_seconds"],
            q["runtime_seconds"],
        )
        fitted[model_name] = {
            "source_latent": source_latent,
            "target_latent": q["latent"],
            "manifest": manifest,
            "tag": tag,
        }
        del src, q
        _cleanup()
    return fitted


def evaluate_subsample(
    subsample_seed: int,
    source_sub: ad.AnnData,
    target: ad.AnnData,
    fitted: dict[str, dict[str, Any]],
    common_classes: dict[str, list[str]],
    run_plans: dict[str, dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prim_rows, delta_rows, cell_rows = [], [], []
    for level, col in LEVELS:
        for class_set_name in ("common", "per_run_eligible"):
            if class_set_name == "common":
                classes = common_classes[level]
            else:
                classes = run_plans[level]["eligible_shared_classes"]
            si, _ = eval_index(source_sub, col, classes)
            ti, n_high_t = eval_index(target, col, classes)
            y_s = source_sub.obs[col].astype(str).to_numpy()[si]
            y_t = target.obs[col].astype(str).to_numpy()[ti]
            src_counts = pd.Series(y_s).value_counts()
            tgt_counts = pd.Series(y_t).value_counts()
            preds = {}
            for model_name, pack in fitted.items():
                z_s = pack["source_latent"][si]
                z_t = pack["target_latent"][ti]
                clf = classify_transfer(z_s, y_s, z_t, y_t, MODEL_SEED)
                preds[model_name] = clf
                man = pack["manifest"]
                for classifier, prefix in (("logreg", "logreg"), ("knn15", "knn15"), ("knn30", "knn30")):
                    row = {
                        "condition": "size_matched_A",
                        "direction": DIRECTION_A,
                        "subsample_seed": subsample_seed,
                        "model_seed": MODEL_SEED,
                        "model": model_name,
                        "source_n": int(source_sub.n_obs),
                        "target_n": int(target.n_obs),
                        "label_level": level,
                        "evaluation_class_set": class_set_name,
                        "n_classes": len(classes),
                        "target_coverage": float(len(ti) / n_high_t) if n_high_t else np.nan,
                        "classifier": classifier,
                        **_metrics_row(clf, prefix),
                        "source_training_epochs": man["source_training_epochs"],
                        "query_adaptation_epochs": man["query_adaptation_epochs"],
                        "source_runtime": man["source_training_runtime"],
                        "query_runtime": man["query_runtime"],
                        "gpu_peak_allocated_mb": man["source_gpu"].get("peak_cuda_allocated_mb"),
                        "gpu_peak_reserved_mb": man["source_gpu"].get("peak_cuda_reserved_mb"),
                        "process_peak_rss_mb": man["query_ram_after"].get("process_rss_mb"),
                        "source_latent_drift": man["adapted_source_drift_vs_source_model"][
                            "mean_euclidean_shift"
                        ],
                    }
                    if classifier == "logreg":
                        row.update(_agreement_cols(z_s, y_s, z_t, y_t))
                        if class_set_name == "common":
                            row.update(_target_only_quality(pack["target_latent"], target, col))
                    prim_rows.append(row)
                if level == "l2" and class_set_name == "common":
                    spec = per_class_table(y_t, clf["logreg_pred"], classes, src_counts, tgt_counts)
                    spec["condition"] = "size_matched_A"
                    spec["direction"] = DIRECTION_A
                    spec["model"] = model_name
                    spec["subsample_seed"] = subsample_seed
                    spec["model_seed"] = MODEL_SEED
                    cell_rows.append(spec)
            boot = paired_bootstrap_delta(
                y_t, preds["scvi_matched"]["logreg_pred"], preds["totalvi"]["logreg_pred"],
                seed=subsample_seed, n_boot=N_BOOT,
            )
            delta_rows.append(
                {
                    "condition": "size_matched_A",
                    "direction": DIRECTION_A,
                    "subsample_seed": subsample_seed,
                    "model_seed": MODEL_SEED,
                    "label_level": level,
                    "evaluation_class_set": class_set_name,
                    "scvi_macro_f1": preds["scvi_matched"]["logreg_macro_f1"],
                    "totalvi_macro_f1": preds["totalvi"]["logreg_macro_f1"],
                    "delta": float(
                        preds["totalvi"]["logreg_macro_f1"] - preds["scvi_matched"]["logreg_macro_f1"]
                    ),
                    "bootstrap_ci_low": boot["bootstrap_ci_low"],
                    "bootstrap_ci_high": boot["bootstrap_ci_high"],
                    "target_n_eval": int(len(ti)),
                    "n_classes": len(classes),
                }
            )
    cell = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    return pd.DataFrame(prim_rows), pd.DataFrame(delta_rows), cell


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def run_phase10b(
    *, subsample_seeds: tuple[int, ...] = SUBSAMPLE_SEEDS, config: dict[str, Any] | None = None
) -> dict[str, Any]:
    tmp = Path.home() / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmp)
    os.environ.setdefault("MPLBACKEND", "Agg")
    cfg = config or load_config()
    paths = _paths()
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    logger = get_logger("phase10b", paths["logs"] / "phase10b.log")
    env = detect_compute_environment()
    logger.info(
        "PHASE 10B env: %s torch=%s scvi=%s gpu=%s ram=%s",
        env.get("hostname"), env.get("torch_version"), env.get("scvi_version"),
        env.get("gpu_names"), _ram(),
    )
    _require_cuda()
    if str(env.get("scvi_version")) != "1.3.3":
        raise RuntimeError(f"Expected scvi-tools 1.3.3, found {env.get('scvi_version')}")

    source_full, target, audit = load_source_target(SOURCE_BATCH, TARGET_BATCH)
    if source_full.n_obs <= SUBSAMPLE_N:
        raise RuntimeError("Source is not larger than the requested subsample size.")

    # ---- 1. subsets + class plans (labels only; no models, no target labels) ----
    subsets: dict[int, dict[str, Any]] = {}
    comp_frames = []
    for s in subsample_seeds:
        idx = subsample_index(source_full.n_obs, SUBSAMPLE_N, s)
        if len(idx) != SUBSAMPLE_N or len(set(idx.tolist())) != SUBSAMPLE_N:
            raise RuntimeError(f"Subsample seed {s} did not yield {SUBSAMPLE_N} unique cells.")
        view = source_full[idx]
        plans = {lvl: class_plan(view, target, col) for lvl, col in LEVELS}
        tvds = {}
        for lvl, col in LEVELS:
            frame, tvd = composition_rows(
                source_full.obs[col].astype(str).to_numpy(),
                source_full.obs[col].astype(str).to_numpy()[idx],
                lvl,
                s,
            )
            comp_frames.append(frame)
            tvds[lvl] = tvd
        cells = pd.DataFrame(
            {
                "cell_id": source_full.obs_names.astype(str).to_numpy()[idx],
                "source_row_index": idx,
                "cell_type_l1": source_full.obs["cell_type_l1"].astype(str).to_numpy()[idx],
                "cell_type_l2": source_full.obs["cell_type_l2"].astype(str).to_numpy()[idx],
            }
        )
        cells.to_csv(paths["subsamples"] / f"pbmc10k_n{SUBSAMPLE_N}_subsample_seed{s}.csv", index=False)
        subsets[s] = {"index": idx, "plans": plans, "tvd": tvds}
        logger.info(
            "subsample seed %s: n=%s eligible_l2=%s TVD_l2=%.4f TVD_l1=%.4f",
            s, len(idx), plans["l2"]["n_shared_classes"], tvds["l2"], tvds["l1"],
        )

    # ---- 2. common evaluation class set ----
    orig_a_plans = {lvl: class_plan(source_full, target, col) for lvl, col in LEVELS}
    orig_b_plans = {lvl: class_plan(target, source_full, col) for lvl, col in LEVELS}
    common_classes: dict[str, list[str]] = {}
    class_set_audit: dict[str, Any] = {}
    for lvl, _ in LEVELS:
        sets = [set(subsets[s]["plans"][lvl]["eligible_shared_classes"]) for s in subsample_seeds]
        inter = set.intersection(*sets) if sets else set()
        inter &= set(orig_a_plans[lvl]["eligible_shared_classes"])
        inter &= set(orig_b_plans[lvl]["eligible_shared_classes"])
        common_classes[lvl] = sorted(inter)
        union_all = set().union(*sets) | set(orig_a_plans[lvl]["eligible_shared_classes"]) | set(
            orig_b_plans[lvl]["eligible_shared_classes"]
        )
        class_set_audit[lvl] = {
            "common_classes": common_classes[lvl],
            "n_common": len(common_classes[lvl]),
            "original_A_eligible": orig_a_plans[lvl]["eligible_shared_classes"],
            "original_B_eligible": orig_b_plans[lvl]["eligible_shared_classes"],
            "per_subsample_eligible": {
                str(s): subsets[s]["plans"][lvl]["eligible_shared_classes"] for s in subsample_seeds
            },
            "dropped_relative_to_union": sorted(union_all - inter),
        }
        logger.info("common %s classes: %s -> %s", lvl, len(common_classes[lvl]), common_classes[lvl])
    if not common_classes["l2"]:
        raise RuntimeError("Common l2 evaluation class set is empty. STOP.")

    save_json(
        {
            "direction": DIRECTION_A,
            "subsample_n": SUBSAMPLE_N,
            "subsample_seeds": list(subsample_seeds),
            "model_seed": MODEL_SEED,
            "source_full_n": int(source_full.n_obs),
            "target_n": int(target.n_obs),
            "min_source_support": MIN_SOURCE_SUPPORT,
            "min_target_support": MIN_TARGET_SUPPORT,
            "dataset_audit": audit,
            "class_sets": class_set_audit,
            "total_variation_distance": {str(s): subsets[s]["tvd"] for s in subsample_seeds},
            "subsampling": "uniform without replacement, numpy default_rng(subsample_seed), labels not used",
        },
        paths["logs"] / "phase10b_design_audit.json",
    )
    pd.concat(comp_frames, ignore_index=True).to_csv(
        paths["tables"] / "source_subsample_composition.csv", index=False
    )

    # ---- 3. reference conditions from stored PHASE 10 embeddings ----
    ref_prim_a, ref_delta_a, ref_cell_a = reference_condition(
        DIRECTION_A, "original_A", source_full, target, common_classes, logger
    )
    ref_prim_b, ref_delta_b, ref_cell_b = reference_condition(
        DIRECTION_B, "original_B", target, source_full, common_classes, logger
    )

    # ---- 4. size-matched runs ----
    prim_parts = [ref_prim_a, ref_prim_b]
    delta_parts = [ref_delta_a, ref_delta_b]
    cell_parts = [ref_cell_a, ref_cell_b]
    for s in subsample_seeds:
        logger.info("=== size-matched subsample seed %s ===", s)
        source_sub = ensure_csr(source_full[subsets[s]["index"]].copy())
        if source_sub.n_obs != SUBSAMPLE_N:
            raise RuntimeError("Subsampled source has the wrong number of cells.")
        fitted = train_subsample(s, source_sub, target, cfg, logger)
        p, d, c = evaluate_subsample(s, source_sub, target, fitted, common_classes, subsets[s]["plans"])
        prim_parts.append(p)
        delta_parts.append(d)
        cell_parts.append(c)
        del source_sub, fitted
        _cleanup()

    primary = pd.concat(prim_parts, ignore_index=True)
    delta = pd.concat(delta_parts, ignore_index=True)
    celltype = pd.concat([c for c in cell_parts if len(c)], ignore_index=True)
    primary.to_csv(paths["tables"] / "source_size_primary.csv", index=False)
    delta.to_csv(paths["tables"] / "source_size_delta.csv", index=False)
    celltype.to_csv(paths["tables"] / "source_size_celltype_specific.csv", index=False)

    summary = build_summary(delta, primary, celltype, common_classes)
    summary.to_csv(paths["tables"] / "source_size_summary.csv", index=False)
    rare = rare_class_stability(celltype, source_full, target, subsets, subsample_seeds, common_classes)
    rare.to_csv(paths["tables"] / "source_size_rare_class_stability.csv", index=False)

    write_phase10b_figures(paths["tables"], paths["figures"])

    del source_full, target
    _cleanup()
    out = {
        "subsample_seeds": list(subsample_seeds),
        "model_seed": MODEL_SEED,
        "subsample_n": SUBSAMPLE_N,
        "n_common_l2_classes": len(common_classes["l2"]),
        "n_common_l1_classes": len(common_classes["l1"]),
        "tables": sorted(p.name for p in paths["tables"].glob("*.csv")),
        "ram": _ram(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    save_json(out, paths["logs"] / "phase10b_summary.json")
    return out


def build_summary(
    delta: pd.DataFrame, primary: pd.DataFrame, celltype: pd.DataFrame, common_classes: dict[str, list[str]]
) -> pd.DataFrame:
    rows = []
    for level in ("l2", "l1"):
        d = delta[(delta["label_level"] == level) & (delta["evaluation_class_set"] == "common")]
        a = d[d["condition"] == "original_A"]
        b = d[d["condition"] == "original_B"]
        m = d[d["condition"] == "size_matched_A"]
        a0 = a[a["model_seed"] == MODEL_SEED]
        row = {
            "label_level": level,
            "n_common_classes": len(common_classes[level]),
            "original_A_scVI": a["scvi_macro_f1"].mean(),
            "original_A_totalVI": a["totalvi_macro_f1"].mean(),
            "original_A_delta": a["delta"].mean(),
            "original_A_delta_sd": a["delta"].std(ddof=1),
            "original_A_modelseed0_scVI": a0["scvi_macro_f1"].mean(),
            "original_A_modelseed0_totalVI": a0["totalvi_macro_f1"].mean(),
            "original_A_modelseed0_delta": a0["delta"].mean(),
            "size_matched_A_scVI_mean": m["scvi_macro_f1"].mean(),
            "size_matched_A_scVI_SD": m["scvi_macro_f1"].std(ddof=1),
            "size_matched_A_totalVI_mean": m["totalvi_macro_f1"].mean(),
            "size_matched_A_totalVI_SD": m["totalvi_macro_f1"].std(ddof=1),
            "size_matched_A_delta_mean": m["delta"].mean(),
            "size_matched_A_delta_SD": m["delta"].std(ddof=1),
            "original_B_scVI": b["scvi_macro_f1"].mean(),
            "original_B_totalVI": b["totalvi_macro_f1"].mean(),
            "original_B_delta": b["delta"].mean(),
            "original_B_delta_sd": b["delta"].std(ddof=1),
        }
        row["delta_size_effect"] = row["size_matched_A_delta_mean"] - row["original_A_delta"]
        row["delta_size_effect_vs_modelseed0"] = (
            row["size_matched_A_delta_mean"] - row["original_A_modelseed0_delta"]
        )
        row["scvi_size_effect"] = row["size_matched_A_scVI_mean"] - row["original_A_scVI"]
        row["totalvi_size_effect"] = row["size_matched_A_totalVI_mean"] - row["original_A_totalVI"]
        row["scvi_degrades_more_than_totalvi"] = bool(
            row["scvi_size_effect"] < row["totalvi_size_effect"]
        )
        rows.append(row)
    return pd.DataFrame(rows)


def rare_class_stability(
    celltype: pd.DataFrame,
    source_full: ad.AnnData,
    target: ad.AnnData,
    subsets: dict[int, dict[str, Any]],
    subsample_seeds: tuple[int, ...],
    common_classes: dict[str, list[str]],
) -> pd.DataFrame:
    high_s = source_full.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    high_t = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    ys = source_full.obs["cell_type_l2"].astype(str).to_numpy()
    yt = target.obs["cell_type_l2"].astype(str).to_numpy()
    full_counts = pd.Series(ys[high_s]).value_counts()
    tgt_counts = pd.Series(yt[high_t]).value_counts()
    sub_counts = {}
    for s in subsample_seeds:
        idx = subsets[s]["index"]
        mask = high_s[idx]
        sub_counts[s] = pd.Series(ys[idx][mask]).value_counts()
    sm = celltype[celltype["condition"] == "size_matched_A"]
    rows = []
    for ct in common_classes["l2"]:
        counts = [int(sub_counts[s].get(ct, 0)) for s in subsample_seeds]
        row = {
            "cell_type": ct,
            "full_source_count": int(full_counts.get(ct, 0)),
            "mean_subsample_source_count": float(np.mean(counts)),
            "min_subsample_source_count": int(np.min(counts)),
            "target_count": int(tgt_counts.get(ct, 0)),
        }
        for model in ("scvi_matched", "totalvi"):
            vals = sm[(sm["cell_type"] == ct) & (sm["model"] == model)]["f1"]
            row[f"{model}_f1_mean"] = float(vals.mean()) if len(vals) else np.nan
            row[f"{model}_f1_sd_across_subsets"] = float(vals.std(ddof=1)) if len(vals) > 1 else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    out["delta_f1_mean"] = out["totalvi_f1_mean"] - out["scvi_matched_f1_mean"]
    return out.sort_values("full_source_count").reset_index(drop=True)
