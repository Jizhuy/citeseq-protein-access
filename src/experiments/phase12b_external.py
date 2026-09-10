"""PHASE 12B: independent Lawlor Baseline donor-held-out external validation.

Does not modify PBMC10k/PBMC5k historical artifacts. Does not use GSE164378.
"""

from __future__ import annotations

import gc
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import torch
from scipy import sparse
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    log_loss,
)

from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.experiments.phase10_cross_dataset import (
    _apply_run_seed,
    _train_kwargs,
    ensure_csr,
    train_scvi_source,
    train_totalvi_source,
)
from src.experiments.phase11_uncertainty import (
    ECE_BINS,
    calibration_metrics,
    classification_metrics,
    coverage_table,
    error_detection,
    fit_transfer_classifier,
    predictive_uncertainty,
    risk_coverage,
)
from src.experiments.phase11b_robustness import query_adapt_and_save
from src.experiments.phase11b_stats import assert_post_adaptation, latent_drift
from src.models.run_scvi import epochs_completed_from_history, history_to_frame
from src.models.run_scvi_matched import MATCHED_HYPERPARAMETERS
from src.models.run_totalvi import DEFAULT_HYPERPARAMETERS, PROTEIN_OBSM_KEY
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

MODELS = ("scvi_matched", "totalvi")
SEEDS = tuple(range(10))
FOLDS = tuple(range(5))
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
LABEL_COL = "author_cell_type"
INTEGRITY_TOL = 1e-5


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def phase12b_root() -> Path:
    return resolve_path("results/phase12b_external_validation/formal_validation")


def phase12b_paths() -> dict[str, Path]:
    root = phase12b_root()
    return {
        "root": root,
        "processed": root / "processed",
        "models": root / "models",
        "embeddings": root / "embeddings",
        "posterior": root / "posterior",
        "probabilities": root / "probabilities",
        "tables": root / "tables",
        "figures": root / "figures",
        "logs": root / "logs",
        "done": root / "done",
    }


def ensure_dirs() -> dict[str, Path]:
    paths = phase12b_paths()
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def experiment_tag(fold: int, model: str, seed: int) -> str:
    return f"fold{fold}__{model}__seed{seed:02d}"


def cleanup() -> None:
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    try:
        import scvi

        if hasattr(scvi.settings, "manager_store"):
            scvi.settings.manager_store.clear()
    except Exception:
        pass
    gc.collect()


def ram_snapshot() -> dict[str, float]:
    vm = psutil.virtual_memory()
    return {
        "available_mb": float(vm.available) / 1024**2,
        "percent": float(vm.percent),
        "process_rss_mb": float(peak_rss_mb()),
    }


def load_prepared() -> ad.AnnData:
    path = phase12b_paths()["processed"] / "lawlor_baseline_sng_shared.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Missing prepared object: {path}")
    adata = ad.read_h5ad(path)
    if not sparse.issparse(adata.X):
        raise RuntimeError("X must remain sparse")
    adata.X = sparse.csr_matrix(adata.X)
    if "counts" in adata.layers:
        adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])
    if PROTEIN_OBSM_KEY not in adata.obsm:
        raise RuntimeError(f"Missing obsm['{PROTEIN_OBSM_KEY}']")
    if (adata.obs["condition"].astype(str) != "Baseline").any():
        raise RuntimeError("Non-Baseline cells present")
    return adata


def load_folds() -> pd.DataFrame:
    path = phase12b_paths()["tables"] / "phase12b_donor_folds.csv"
    return pd.read_csv(path)


def load_primary_classes() -> list[str]:
    path = phase12b_paths()["logs"] / "primary_class_set.json"
    return list(json.loads(path.read_text())["primary_eligible_classes"])


def split_fold(adata: ad.AnnData, fold: int, folds: pd.DataFrame) -> tuple[ad.AnnData, ad.AnnData, dict]:
    tgt_donors = folds.loc[
        (folds["fold"] == fold) & (folds["role_when_fold_target"] == "target"), "donor"
    ].astype(str).tolist()
    src_donors = folds.loc[
        (folds["fold"] == fold) & (folds["role_when_fold_target"] == "source"), "donor"
    ].astype(str).tolist()
    src_mask = adata.obs["donor"].astype(str).isin(src_donors).to_numpy()
    tgt_mask = adata.obs["donor"].astype(str).isin(tgt_donors).to_numpy()
    if src_mask.sum() == 0 or tgt_mask.sum() == 0:
        raise RuntimeError(f"Empty split for fold {fold}")
    source = ensure_csr(adata[src_mask].copy())
    target = ensure_csr(adata[tgt_mask].copy())
    # categorical batch for scvi
    source.obs["batch"] = source.obs["batch"].astype("category")
    target.obs["batch"] = target.obs["batch"].astype("category")
    meta = {
        "fold": fold,
        "source_donors": src_donors,
        "target_donors": tgt_donors,
        "n_source_all": int(source.n_obs),
        "n_target_all": int(target.n_obs),
        "n_source_labeled": int(source.obs["has_author_label"].sum()),
        "n_target_labeled": int(target.obs["has_author_label"].sum()),
    }
    return source, target, meta


def labeled_index(adata: ad.AnnData, classes: list[str]) -> np.ndarray:
    lab = adata.obs[LABEL_COL].astype(str)
    mask = adata.obs["has_author_label"].to_numpy() & lab.isin(classes).to_numpy()
    return np.where(mask)[0]


def multiclass_brier(y_true: np.ndarray, probs: np.ndarray, classes: np.ndarray) -> float:
    class_to_i = {c: i for i, c in enumerate(classes)}
    Y = np.zeros_like(probs)
    for i, y in enumerate(y_true):
        Y[i, class_to_i[y]] = 1.0
    return float(np.mean(np.sum((probs - Y) ** 2, axis=1)))


def evaluate_run(
    *,
    z_src: np.ndarray,
    z_tgt: np.ndarray,
    u_latent_tgt: np.ndarray,
    source: ad.AnnData,
    target: ad.AnnData,
    classes: list[str],
    seed: int,
) -> dict[str, Any]:
    si = labeled_index(source, classes)
    ti = labeled_index(target, classes)
    if si.size == 0 or ti.size == 0:
        raise RuntimeError("No labeled cells for evaluation")
    y_src = source.obs[LABEL_COL].astype(str).to_numpy()[si]
    y_tgt = target.obs[LABEL_COL].astype(str).to_numpy()[ti]
    # Restrict to classes present in BOTH source and target labeled sets
    present = sorted(set(y_src) & set(y_tgt) & set(classes))
    if len(present) < 2:
        raise RuntimeError(f"Too few shared labeled classes: {present}")
    src_keep = np.isin(y_src, present)
    tgt_keep = np.isin(y_tgt, present)
    si, ti = si[src_keep], ti[tgt_keep]
    y_src, y_tgt = y_src[src_keep], y_tgt[tgt_keep]
    z_s, z_t = z_src[si], z_tgt[ti]
    u_lat = u_latent_tgt[ti]

    fit = fit_transfer_classifier(z_s, y_src, z_t, seed=seed)
    # Align probability columns to a fixed class order (present)
    class_order = list(fit["classes"])
    probs = fit["probs"]
    pred = fit["pred"].astype(str)
    metrics = classification_metrics(y_tgt, pred)
    error = (pred != y_tgt).astype(int)
    u = predictive_uncertainty(probs)
    u_pred = u["u_pred"]
    pred_det = error_detection(u_pred, error)
    lat_det = error_detection(u_lat, error)

    # NLL on true class
    class_to_i = {c: i for i, c in enumerate(class_order)}
    true_idx = np.array([class_to_i[y] for y in y_tgt], dtype=int)
    p_true = probs[np.arange(probs.shape[0]), true_idx]
    nll = float(-np.mean(np.log(np.clip(p_true, 1e-12, 1.0))))
    brier = multiclass_brier(y_tgt, probs, np.asarray(class_order, dtype=object))
    cal = calibration_metrics(probs, y_tgt, np.asarray(class_order, dtype=object), pred)
    rc = risk_coverage(u_pred, (error == 0).astype(float))
    cov_rows = coverage_table(u_pred, y_tgt, pred)

    # Spearman correlations
    from scipy.stats import spearmanr

    spe_pred_lat = spearmanr(u_pred, u_lat)
    spe_lat_nll = spearmanr(u_lat, -np.log(np.clip(p_true, 1e-12, 1.0)))

    # per-class retention under abstention
    order = np.argsort(u_pred, kind="mergesort")
    n = y_tgt.size
    celltype_rows = []
    for ct in present:
        ct_mask = y_tgt == ct
        n_ct = int(ct_mask.sum())
        row = {
            "cell_type": ct,
            "source_support": int((y_src == ct).sum()),
            "target_support": n_ct,
            "accuracy": float(accuracy_score(y_tgt[ct_mask], pred[ct_mask])) if n_ct else float("nan"),
            "f1": float(f1_score(y_tgt == ct, pred == ct, zero_division=0)),
            "median_predictive_uncertainty": float(np.median(u_pred[ct_mask])) if n_ct else float("nan"),
            "median_latent_uncertainty": float(np.median(u_lat[ct_mask])) if n_ct else float("nan"),
        }
        for cov in COVERAGE_GRID:
            if cov >= 1.0:
                row[f"retention_{int(cov*100)}"] = 1.0
                continue
            keep = max(1, int(round(cov * n)))
            kept = order[:keep]
            kept_ct = ct_mask[kept].sum()
            row[f"retention_{int(cov*100)}"] = float(kept_ct / n_ct) if n_ct else float("nan")
        celltype_rows.append(row)

    return {
        "classes": present,
        "n_source_labeled_eval": int(si.size),
        "n_target_labeled_eval": int(ti.size),
        "y_true": y_tgt,
        "y_pred": pred,
        "probs": probs,
        "class_order": class_order,
        "u_pred": u_pred,
        "u_latent": u_lat,
        "error": error,
        "metrics": metrics,
        "error_prevalence": float(error.mean()),
        "predictive_error_auroc": pred_det["auroc"],
        "predictive_error_auprc": pred_det["auprc"],
        "latent_error_auroc": lat_det["auroc"],
        "latent_error_auprc": lat_det["auprc"],
        "nll": nll,
        "brier": brier,
        "ece_equal_frequency": cal["ece_equal_frequency"],
        "ece_equal_width": cal["ece_equal_width"],
        "aurc": rc["aurc"],
        "coverage_rows": cov_rows,
        "celltype_rows": celltype_rows,
        "spearman_u_pred_u_latent": float(spe_pred_lat.correlation) if spe_pred_lat.correlation is not None else float("nan"),
        "spearman_u_latent_neglog_p_true": float(spe_lat_nll.correlation) if spe_lat_nll.correlation is not None else float("nan"),
        "target_eval_index": ti,
        "source_eval_index": si,
    }


def run_experiment(
    *,
    fold: int,
    model_name: str,
    seed: int,
    adata: ad.AnnData | None = None,
    cfg: dict[str, Any] | None = None,
    logger=None,
) -> dict[str, Any] | None:
    paths = ensure_dirs()
    logger = logger or get_logger("phase12b", paths["logs"] / "phase12b_train.log")
    cfg = cfg or load_config()
    tag = experiment_tag(fold, model_name, seed)
    done_marker = paths["done"] / f"{tag}.json"
    if done_marker.exists():
        try:
            if json.loads(done_marker.read_text()).get("complete"):
                logger.info("skip (done): %s", tag)
                return None
        except json.JSONDecodeError:
            logger.warning("corrupt done marker for %s; rerunning", tag)

    if not torch.cuda.is_available():
        raise RuntimeError("PHASE 12B requires CUDA")

    adata = adata if adata is not None else load_prepared()
    folds = load_folds()
    classes = load_primary_classes()
    source, target, split_meta = split_fold(adata, fold, folds)

    # Target labels must not influence training/adaptation
    target_labels_backup = target.obs[LABEL_COL].astype(str).copy()
    # Do not remove the column (needed later); never pass it into models.

    mdir = paths["models"] / f"fold{fold}" / model_name / f"seed{seed:02d}"
    source_dir = mdir / "source_model"
    adapted_dir = mdir / "adapted_query_model"
    max_epochs = int(cfg["models"]["max_epochs"])
    batch_size = int(cfg["models"]["batch_size"])
    query_epochs = int(cfg.get("phase10", {}).get("query_max_epochs", 200))

    logger.info(
        "=== %s source=%s target=%s genes=%s proteins=%s ===",
        tag, source.n_obs, target.n_obs, source.n_vars,
        source.obsm[PROTEIN_OBSM_KEY].shape[1],
    )
    set_global_seed(seed)
    ram_before = ram_snapshot()
    reset_cuda_peak_stats()
    timer = Timer()
    timer.start()

    from scvi.model import SCVI, TOTALVI

    if model_name == "scvi_matched":
        hparams = dict(MATCHED_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
        hparams["batch_size"] = batch_size
        src = train_scvi_source(source, hparams, source_dir, logger, run_seed=None)
        extra_train = None
        model_cls = SCVI
    elif model_name == "totalvi":
        hparams = dict(DEFAULT_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["batch_size"] = batch_size
        src = train_totalvi_source(source, hparams, source_dir, logger, run_seed=None)
        extra_train = {
            "lr": DEFAULT_HYPERPARAMETERS["lr"],
            "reduce_lr_on_plateau": DEFAULT_HYPERPARAMETERS["reduce_lr_on_plateau"],
        }
        model_cls = TOTALVI
    else:
        raise ValueError(model_name)

    q = query_adapt_and_save(
        model_cls,
        target,
        source,
        source_dir,
        adapted_dir,
        query_epochs,
        batch_size,
        extra_train,
        tag,
        logger,
        run_seed=None,
    )
    source_latent = q["source_latent_post"]
    assert_post_adaptation(source_latent, src["latent"], q["source_latent_post"], tag)
    drift = latent_drift(src["latent"], source_latent)

    emb = paths["embeddings"] / f"fold{fold}"
    emb.mkdir(parents=True, exist_ok=True)
    np.save(emb / f"{tag}_source_latent_post.npy", source_latent.astype(np.float32))
    np.save(emb / f"{tag}_target_latent_post.npy", q["target_latent"].astype(np.float32))
    pd.DataFrame({"cell_id": source.obs_names.astype(str)}).to_csv(emb / f"{tag}_source_cells.csv", index=False)
    pd.DataFrame({"cell_id": target.obs_names.astype(str)}).to_csv(emb / f"{tag}_target_cells.csv", index=False)

    post = paths["posterior"] / f"fold{fold}"
    post.mkdir(parents=True, exist_ok=True)
    u_latent_target = np.asarray(q["target_var"]).mean(axis=1)
    np.savez_compressed(
        post / f"{tag}_posterior.npz",
        target_mu=np.asarray(q["target_mu"]).astype(np.float32),
        target_var=np.asarray(q["target_var"]).astype(np.float32),
        u_latent_target=u_latent_target.astype(np.float32),
    )

    # Restore labels for evaluation (they were never used above)
    target.obs[LABEL_COL] = target_labels_backup
    ev = evaluate_run(
        z_src=source_latent,
        z_tgt=q["target_latent"],
        u_latent_tgt=u_latent_target,
        source=source,
        target=target,
        classes=classes,
        seed=seed,
    )

    prob_dir = paths["probabilities"] / f"fold{fold}"
    prob_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        prob_dir / f"{tag}_probs.npz",
        probs=ev["probs"].astype(np.float32),
        y_true=ev["y_true"].astype(object),
        y_pred=ev["y_pred"].astype(object),
        classes=np.asarray(ev["class_order"], dtype=object),
        u_pred=ev["u_pred"].astype(np.float32),
        u_latent=ev["u_latent"].astype(np.float32),
        target_eval_index=ev["target_eval_index"],
    )

    runtime = float(timer.stop())
    gpu = cuda_memory_stats()
    ram_after = ram_snapshot()

    primary = {
        "fold": fold,
        "source_donors": ";".join(split_meta["source_donors"]),
        "target_donors": ";".join(split_meta["target_donors"]),
        "model": model_name,
        "seed": seed,
        "n_source_all": split_meta["n_source_all"],
        "n_target_all": split_meta["n_target_all"],
        "n_source_labeled": split_meta["n_source_labeled"],
        "n_target_labeled": split_meta["n_target_labeled"],
        "n_source_labeled_eval": ev["n_source_labeled_eval"],
        "n_target_labeled_eval": ev["n_target_labeled_eval"],
        "n_genes": int(source.n_vars),
        "n_proteins": int(source.obsm[PROTEIN_OBSM_KEY].shape[1]),
        "macro_f1": ev["metrics"]["macro_f1"],
        "balanced_accuracy": ev["metrics"]["balanced_accuracy"],
        "accuracy": ev["metrics"]["accuracy"],
        "weighted_f1": ev["metrics"]["weighted_f1"],
        "error_prevalence": ev["error_prevalence"],
        "predictive_error_auroc": ev["predictive_error_auroc"],
        "predictive_error_auprc": ev["predictive_error_auprc"],
        "latent_error_auroc": ev["latent_error_auroc"],
        "latent_error_auprc": ev["latent_error_auprc"],
        "nll": ev["nll"],
        "brier": ev["brier"],
        "ece_equal_frequency": ev["ece_equal_frequency"],
        "ece_equal_width": ev["ece_equal_width"],
        "aurc": ev["aurc"],
        "spearman_u_pred_u_latent": ev["spearman_u_pred_u_latent"],
        "spearman_u_latent_neglog_p_true": ev["spearman_u_latent_neglog_p_true"],
        "runtime_sec": runtime,
        "peak_gpu_mb": float(gpu.get("peak_cuda_allocated_mb") or np.nan),
        "peak_rss_mb": ram_after["process_rss_mb"],
        "integrity_ok": bool(q.get("integrity_ok", False)),
        "latent_drift_mean_cell": float(drift.get("mean_cell_displacement", np.nan)),
        "batch_key": "batch(=run_identifier)",
        "timestamp_utc": utc_now(),
    }
    save_json(primary, paths["logs"] / f"{tag}_primary.json")
    pd.DataFrame(ev["coverage_rows"]).assign(fold=fold, model=model_name, seed=seed).to_csv(
        paths["tables"] / f"{tag}_coverage.csv", index=False
    )
    pd.DataFrame(ev["celltype_rows"]).assign(fold=fold, model=model_name, seed=seed).to_csv(
        paths["tables"] / f"{tag}_celltype.csv", index=False
    )

    save_json(
        {
            "complete": True,
            "tag": tag,
            "timestamp_utc": utc_now(),
            "primary": primary,
            "ram_before": ram_before,
            "ram_after": ram_after,
        },
        done_marker,
    )
    del source, target, src, q, ev
    cleanup()
    logger.info("done %s macro_f1=%.4f pred_auroc=%.4f lat_auroc=%.4f",
                tag, primary["macro_f1"], primary["predictive_error_auroc"], primary["latent_error_auroc"])
    return primary


def run_gate(jobs: list[tuple[int, str, int]], logger=None) -> list[dict[str, Any]]:
    paths = ensure_dirs()
    logger = logger or get_logger("phase12b", paths["logs"] / "phase12b_train.log")
    cfg = load_config()
    adata = load_prepared()
    results = []
    for fold, model, seed in jobs:
        try:
            out = run_experiment(fold=fold, model_name=model, seed=seed, adata=adata, cfg=cfg, logger=logger)
            if out is not None:
                results.append(out)
        except Exception as exc:
            logger.exception("FAILED %s: %s", experiment_tag(fold, model, seed), exc)
            fail_path = paths["logs"] / f"FAIL__{experiment_tag(fold, model, seed)}.json"
            save_json({"tag": experiment_tag(fold, model, seed), "error": str(exc), "timestamp_utc": utc_now()}, fail_path)
            raise
        cleanup()
    return results
