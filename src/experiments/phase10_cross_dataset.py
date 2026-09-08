#!/usr/bin/env python
"""PHASE 10: cross-dataset transfer / query generalization on CUDA.

Trains source-only scVI_matched and totalVI, then maps the held-out dataset
with scvi-tools 1.3.3 scArches (prepare_query_anndata + load_query_data).
Does not overwrite historical Mac/MPS artifacts.
Does not reopen PHASE 9.
"""

from __future__ import annotations

import gc
import json
import os
import warnings
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import torch
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    normalized_mutual_info_score,
    precision_recall_fscore_support,
    silhouette_score,
)
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

from src.evaluation.phase10_plots import write_phase10_figures
from src.evaluation.phase6 import PRIMARY_RESOLUTION, leiden_clusters, neighborhood_purity_per_cell
from src.evaluation.representation_metrics import batch_silhouette
from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.models.run_scvi import (
    _final_metric,
    epochs_completed_from_history,
    history_to_frame,
)
from src.models.run_scvi_matched import MATCHED_HYPERPARAMETERS
from src.models.run_totalvi import DEFAULT_HYPERPARAMETERS, PROTEIN_OBSM_KEY
from src.utils.device import detect_compute_environment
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

DIRECTIONS = (
    ("pbmc10k_to_pbmc5k", "PBMC10k", "PBMC5k"),
    ("pbmc5k_to_pbmc10k", "PBMC5k", "PBMC10k"),
)
MIN_SOURCE_SUPPORT = 20
MIN_TARGET_SUPPORT = 10
N_BOOT = 500
KNN_KS = (15, 30, 50)
CONFIDENCE_RULE = "PHASE 6 high-confidence rule: l2_confidence >= 0.85 (annotation_tier_l2 == high)"
WITHIN_L2_HIGH = {
    "scVI_matched": 0.7218643946226526,
    "totalVI": 0.7792121076915633,
}
WITHIN_L1_HIGH = {
    "scVI_matched": 0.8959883913628512,
    "totalVI": 0.8967101331990678,
}
WITHIN_NOTE = (
    "PHASE 6 within-combined-dataset held-out logreg macro-F1 on high-confidence "
    "cells. Design differs from PHASE 10 source→target transfer; gap is descriptive."
)


def _ram() -> dict[str, float]:
    vm = psutil.virtual_memory()
    return {
        "total_mb": float(vm.total) / (1024.0 * 1024.0),
        "available_mb": float(vm.available) / (1024.0 * 1024.0),
        "percent": float(vm.percent),
        "process_rss_mb": float(peak_rss_mb()),
    }


def _cleanup() -> None:
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def _require_cuda() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("PHASE 10 requires CUDA. torch.cuda.is_available() is False.")


def ensure_csr(adata: ad.AnnData) -> ad.AnnData:
    """Convert counts to CSR without a dense round-trip when already sparse."""
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        adata.layers["counts"] = counts.tocsr()
    else:
        adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])
    adata.X = adata.layers["counts"]
    return adata


def _phase10_root() -> Path:
    return resolve_path("results/phase10_cross_dataset")


def _dir_paths(direction: str) -> dict[str, Path]:
    root = _phase10_root()
    return {
        "models": root / "models" / direction,
        "embeddings": root / "embeddings" / direction,
        "tables": root / "tables",
        "figures": root / "figures",
        "logs": root / "logs",
    }


def load_source_target(source_batch: str, target_batch: str) -> tuple[ad.AnnData, ad.AnnData, dict[str, Any]]:
    path = resolve_path("data/processed/pbmc_cite_combined_inner_annotated.h5ad")
    full = ad.read_h5ad(path)
    ensure_csr(full)
    batch = full.obs["batch"].astype(str)
    source = ensure_csr(full[batch == source_batch].copy())
    target = ensure_csr(full[batch == target_batch].copy())
    audit = {
        "source_h5ad": str(path),
        "source_batch": source_batch,
        "target_batch": target_batch,
        "source_n": int(source.n_obs),
        "target_n": int(target.n_obs),
        "n_genes": int(source.n_vars),
        "n_proteins": int(np.asarray(source.obsm[PROTEIN_OBSM_KEY]).shape[1]),
        "gene_order_identical": list(source.var_names) == list(target.var_names),
        "protein_order_identical": list(source.obsm[PROTEIN_OBSM_KEY].columns)
        == list(target.obsm[PROTEIN_OBSM_KEY].columns),
        "disjoint_cells": len(set(source.obs_names) & set(target.obs_names)) == 0,
        "rna_sum_source": float(source.layers["counts"].sum()),
        "rna_sum_target": float(target.layers["counts"].sum()),
        "protein_sum_source": float(np.nansum(np.asarray(source.obsm[PROTEIN_OBSM_KEY]))),
        "protein_sum_target": float(np.nansum(np.asarray(target.obsm[PROTEIN_OBSM_KEY]))),
        "source_has_labels": all(c in source.obs for c in ("cell_type_l1", "cell_type_l2", "annotation_tier_l2")),
        "target_has_labels": all(c in target.obs for c in ("cell_type_l1", "cell_type_l2", "annotation_tier_l2")),
        "counts_csr": True,
    }
    if not audit["gene_order_identical"] or not audit["protein_order_identical"]:
        raise RuntimeError("Gene or protein order differs between source and target.")
    if not audit["disjoint_cells"]:
        raise RuntimeError("Source and target cells are not disjoint.")
    if source.n_vars != 15792 or audit["n_proteins"] != 14:
        raise RuntimeError("Unexpected gene/protein dimensions.")
    del full
    _cleanup()
    return source, target, audit


def class_plan(source: ad.AnnData, target: ad.AnnData, label_col: str) -> dict[str, Any]:
    """Pre-specified shared-class rule. Computed from labels only, before models."""
    high_s = source.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    high_t = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    ys = source.obs[label_col].astype(str).to_numpy()
    yt = target.obs[label_col].astype(str).to_numpy()
    src_counts = pd.Series(ys[high_s]).value_counts()
    tgt_counts = pd.Series(yt[high_t]).value_counts()
    source_classes = set(src_counts.index)
    target_classes = set(tgt_counts.index)
    intersection = sorted(source_classes & target_classes)
    eligible = [
        c
        for c in intersection
        if int(src_counts.get(c, 0)) >= MIN_SOURCE_SUPPORT and int(tgt_counts.get(c, 0)) >= MIN_TARGET_SUPPORT
    ]
    excluded = {
        c: {"source": int(src_counts.get(c, 0)), "target": int(tgt_counts.get(c, 0))}
        for c in sorted(source_classes | target_classes)
        if c not in eligible
    }
    src_keep = high_s & np.isin(ys, eligible)
    tgt_keep = high_t & np.isin(yt, eligible)
    n_high_t = int(high_t.sum())
    return {
        "label_col": label_col,
        "min_source_support": MIN_SOURCE_SUPPORT,
        "min_target_support": MIN_TARGET_SUPPORT,
        "source_classes": sorted(source_classes),
        "target_classes": sorted(target_classes),
        "intersection": intersection,
        "source_only": sorted(source_classes - target_classes),
        "target_only": sorted(target_classes - source_classes),
        "eligible_shared_classes": eligible,
        "excluded_classes": excluded,
        "n_shared_classes": len(eligible),
        "n_high_target": n_high_t,
        "n_target_in_eval": int(tgt_keep.sum()),
        "target_coverage": float(tgt_keep.sum() / n_high_t) if n_high_t else float("nan"),
        "n_source_in_eval": int(src_keep.sum()),
        "source_eval_index": np.where(src_keep)[0],
        "target_eval_index": np.where(tgt_keep)[0],
        "confidence_rule": CONFIDENCE_RULE,
    }


def _train_kwargs(max_epochs: int, batch_size: int, early_stopping: bool, extra: dict | None = None) -> dict[str, Any]:
    kwargs = {
        "max_epochs": max_epochs,
        "accelerator": "gpu",
        "devices": [0],
        "batch_size": batch_size,
        "train_size": 0.9,
        "early_stopping": early_stopping,
        "early_stopping_monitor": "elbo_validation",
        "early_stopping_patience": 45,
        "early_stopping_mode": "min",
        "check_val_every_n_epoch": 1,
    }
    if extra:
        kwargs.update(extra)
    return kwargs


def _fit_with_oom_retry(train_fn, batch_size: int, logger) -> tuple[Any, int, bool]:
    try:
        return train_fn(batch_size), batch_size, False
    except RuntimeError as exc:
        if "out of memory" not in str(exc).lower() and "cuda" not in str(exc).lower():
            raise
        if batch_size <= 128:
            raise
        logger.warning("CUDA OOM at batch_size=%s. Retrying with 128. Exact error: %s", batch_size, exc)
        torch.cuda.empty_cache()
        return train_fn(128), 128, True


def train_scvi_source(adata: ad.AnnData, hparams: dict[str, Any], out_dir: Path, logger) -> dict[str, Any]:
    import scvi
    from scvi.model import SCVI

    _require_cuda()
    reset_cuda_peak_stats()
    ram_before = _ram()
    adata = ensure_csr(adata.copy())
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

    def _go(bs: int):
        model = SCVI(
            adata,
            n_hidden=hparams["n_hidden"],
            n_latent=hparams["n_latent"],
            n_layers=hparams["n_layers"],
            dropout_rate=hparams["dropout_rate"],
            dispersion=hparams["dispersion"],
            gene_likelihood=hparams["gene_likelihood"],
            use_observed_lib_size=hparams["use_observed_lib_size"],
            latent_distribution=hparams["latent_distribution"],
        )
        timer = Timer()
        timer.start()
        model.train(**_train_kwargs(hparams["max_epochs"], bs, hparams["early_stopping"]))
        runtime = timer.stop()
        return model, runtime

    (model, runtime), used_bs, oom_retry = _fit_with_oom_retry(_go, hparams["batch_size"], logger)
    device = str(next(model.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"scVI source training did not use CUDA; device={device}")
    out_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(out_dir), overwrite=True, save_anndata=False)
    latent = np.asarray(model.get_latent_representation())
    if not np.isfinite(latent).all():
        raise RuntimeError("scVI source latent contains NaN/Inf")
    history = history_to_frame(model.history)
    epochs = epochs_completed_from_history(history, hparams["max_epochs"])
    gpu = cuda_memory_stats()
    del model
    _cleanup()
    return {
        "latent": latent,
        "runtime_seconds": float(runtime),
        "epochs": int(epochs),
        "batch_size": int(used_bs),
        "oom_retry_batch_size_128": bool(oom_retry),
        "device": device,
        "gpu": gpu,
        "ram_before": ram_before,
        "ram_after": _ram(),
        "elbo_train_final": _final_metric(history, "elbo_train"),
        "elbo_validation_final": _final_metric(history, "elbo_validation"),
        "history": history,
    }


def train_totalvi_source(adata: ad.AnnData, hparams: dict[str, Any], out_dir: Path, logger) -> dict[str, Any]:
    from scvi.model import TOTALVI

    _require_cuda()
    reset_cuda_peak_stats()
    ram_before = _ram()
    adata = ensure_csr(adata.copy())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        TOTALVI.setup_anndata(
            adata,
            protein_expression_obsm_key=PROTEIN_OBSM_KEY,
            batch_key="batch",
            layer="counts",
        )

    def _go(bs: int):
        model = TOTALVI(
            adata,
            n_latent=hparams["n_latent"],
            gene_dispersion=hparams["gene_dispersion"],
            protein_dispersion=hparams["protein_dispersion"],
            gene_likelihood=hparams["gene_likelihood"],
            latent_distribution=hparams["latent_distribution"],
            empirical_protein_background_prior=hparams["empirical_protein_background_prior"],
            override_missing_proteins=hparams["override_missing_proteins"],
            n_hidden=hparams["n_hidden"],
            n_layers_encoder=hparams["n_layers_encoder"],
            n_layers_decoder=hparams["n_layers_decoder"],
            dropout_rate_encoder=hparams["dropout_rate_encoder"],
            dropout_rate_decoder=hparams["dropout_rate_decoder"],
        )
        timer = Timer()
        timer.start()
        extra = {"lr": hparams["lr"], "reduce_lr_on_plateau": hparams["reduce_lr_on_plateau"]}
        model.train(**_train_kwargs(hparams["max_epochs"], bs, True, extra))
        runtime = timer.stop()
        return model, runtime

    (model, runtime), used_bs, oom_retry = _fit_with_oom_retry(_go, hparams["batch_size"], logger)
    device = str(next(model.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"totalVI source training did not use CUDA; device={device}")
    out_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(out_dir), overwrite=True, save_anndata=False)
    latent = np.asarray(model.get_latent_representation())
    if not np.isfinite(latent).all():
        raise RuntimeError("totalVI source latent contains NaN/Inf")
    history = history_to_frame(model.history)
    epochs = epochs_completed_from_history(history, hparams["max_epochs"])
    gpu = cuda_memory_stats()
    del model
    _cleanup()
    return {
        "latent": latent,
        "runtime_seconds": float(runtime),
        "epochs": int(epochs),
        "batch_size": int(used_bs),
        "oom_retry_batch_size_128": bool(oom_retry),
        "device": device,
        "gpu": gpu,
        "ram_before": ram_before,
        "ram_after": _ram(),
        "elbo_train_final": _final_metric(history, "elbo_train"),
        "elbo_validation_final": _final_metric(history, "elbo_validation"),
        "history": history,
    }


def query_adapt(
    model_cls,
    target: ad.AnnData,
    source: ad.AnnData,
    reference_dir: Path,
    max_epochs: int,
    batch_size: int,
    extra_train: dict | None,
    logger,
) -> dict[str, Any]:
    """scArches query adaptation.

    After adaptation both source and target cells are re-encoded with the
    *adapted* model. scArches leaves ``l_encoder`` trainable and lets gradients
    reach every batch one-hot column of the first encoder layer, including the
    source batch column, so the source embedding shifts during adaptation.
    Encoding both sets with the adapted model is the documented scvi-tools
    scArches procedure and is what puts them in one coordinate system.
    """
    _require_cuda()
    reset_cuda_peak_stats()
    ram_before = _ram()
    query = target.copy()
    model_cls.prepare_query_anndata(query, str(reference_dir))
    qmodel = model_cls.load_query_data(query, str(reference_dir), accelerator="gpu", device="auto")
    trainable = [n for n, p in qmodel.module.named_parameters() if p.requires_grad]
    plan = {"weight_decay": 0.0}
    extra = {"plan_kwargs": plan}
    if extra_train:
        extra.update(extra_train)

    def _go(bs: int):
        timer = Timer()
        timer.start()
        qmodel.train(**_train_kwargs(max_epochs, bs, True, extra))
        return timer.stop()

    try:
        runtime, used_bs, oom_retry = _go(batch_size), batch_size, False
    except RuntimeError as exc:
        if "out of memory" not in str(exc).lower() or batch_size <= 128:
            raise
        logger.warning("Query CUDA OOM at batch_size=%s. Retrying 128.", batch_size)
        torch.cuda.empty_cache()
        runtime, used_bs, oom_retry = _go(128), 128, True
    device = str(next(qmodel.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"Query adaptation did not use CUDA; device={device}")
    latent = np.asarray(qmodel.get_latent_representation())
    if not np.isfinite(latent).all():
        raise RuntimeError("Query latent contains NaN/Inf")
    # Re-encode source cells with the adapted model so both sets share one space.
    source_ref = ensure_csr(source.copy())
    source_latent_adapted = np.asarray(qmodel.get_latent_representation(source_ref))
    if not np.isfinite(source_latent_adapted).all():
        raise RuntimeError("Adapted-model source latent contains NaN/Inf")
    if source_latent_adapted.shape[0] != source.n_obs:
        raise RuntimeError("Adapted-model source latent has the wrong number of cells.")
    del source_ref
    history = history_to_frame(qmodel.history)
    epochs = epochs_completed_from_history(history, max_epochs)
    gpu = cuda_memory_stats()
    del qmodel, query
    _cleanup()
    return {
        "latent": latent,
        "source_latent_adapted": source_latent_adapted,
        "runtime_seconds": float(runtime),
        "epochs": int(epochs),
        "batch_size": int(used_bs),
        "oom_retry_batch_size_128": bool(oom_retry),
        "device": device,
        "gpu": gpu,
        "ram_before": ram_before,
        "ram_after": _ram(),
        "n_trainable_params": len(trainable),
        "trainable_param_names": trainable,
        "procedure": "scArches query adaptation: prepare_query_anndata + load_query_data + train",
        "zero_shot": False,
        "history": history,
    }


def _clf_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "macro_f1": float(f1_score(y, pred, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y, pred, average="weighted", zero_division=0)),
    }


def classify_transfer(z_src: np.ndarray, y_src: np.ndarray, z_tgt: np.ndarray, y_tgt: np.ndarray, seed: int) -> dict[str, Any]:
    logreg = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed, class_weight=None)
    logreg.fit(z_src, y_src)
    logreg_pred = logreg.predict(z_tgt)
    out = {"logreg_pred": logreg_pred, "classes": list(logreg.classes_)}
    out.update({f"logreg_{k}": v for k, v in _clf_metrics(y_tgt, logreg_pred).items()})
    for k in (15, 30):
        knn = KNeighborsClassifier(n_neighbors=k)
        knn.fit(z_src, y_src)
        pred = knn.predict(z_tgt)
        out[f"knn{k}_pred"] = pred
        out.update({f"knn{k}_{name}": val for name, val in _clf_metrics(y_tgt, pred).items()})
    return out


def neighbor_agreement(z_src: np.ndarray, y_src: np.ndarray, z_tgt: np.ndarray, y_tgt: np.ndarray, k: int) -> dict[str, float]:
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(z_src)
    idx = nn.kneighbors(z_tgt, return_distance=False)
    agree = np.array([float(np.mean(y_src[neigh] == y_tgt[i])) for i, neigh in enumerate(idx)])
    return {"mean_agreement": float(agree.mean()), "median_agreement": float(np.median(agree))}


def paired_bootstrap_delta(y: np.ndarray, pred_scvi: np.ndarray, pred_tot: np.ndarray, seed: int = 0, n_boot: int = N_BOOT) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    draws = []
    n = len(y)
    for _ in range(n_boot):
        boot = rng.choice(n, size=n, replace=True)
        draws.append(
            f1_score(y[boot], pred_tot[boot], average="macro", zero_division=0)
            - f1_score(y[boot], pred_scvi[boot], average="macro", zero_division=0)
        )
    arr = np.asarray(draws)
    return {
        "delta_mean": float(arr.mean()),
        "bootstrap_ci_low": float(np.quantile(arr, 0.025)),
        "bootstrap_ci_high": float(np.quantile(arr, 0.975)),
        "n_boot": n_boot,
    }


def per_class_table(y: np.ndarray, pred: np.ndarray, classes: list[str], src_counts: pd.Series, tgt_counts: pd.Series) -> pd.DataFrame:
    p, r, f, _ = precision_recall_fscore_support(y, pred, labels=classes, zero_division=0)
    return pd.DataFrame(
        {
            "cell_type": classes,
            "source_support": [int(src_counts.get(c, 0)) for c in classes],
            "target_support": [int(tgt_counts.get(c, 0)) for c in classes],
            "precision": p,
            "recall": r,
            "f1": f,
        }
    )


def _latent_drift(z_source_model: np.ndarray, z_adapted_model: np.ndarray) -> dict[str, float]:
    """How far source cells move when re-encoded by the adapted model.

    Reported for transparency about what scArches adaptation does to the
    reference embedding. Large drift is not an error, but it is the reason the
    adapted-model source embedding is the one used for cross-space metrics.
    """
    shift = np.linalg.norm(z_adapted_model - z_source_model, axis=1)
    scale = float(np.linalg.norm(z_source_model, axis=1).mean())
    return {
        "mean_euclidean_shift": float(shift.mean()),
        "median_euclidean_shift": float(np.median(shift)),
        "max_euclidean_shift": float(shift.max()),
        "mean_source_latent_norm": scale,
        "relative_shift": float(shift.mean() / scale) if scale else float("nan"),
    }


def _save_latent(path: Path, latent: np.ndarray, adata: ad.AnnData) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, latent)
    cells_path = path.with_name(path.name.replace("_latent.npy", "_cells.csv"))
    pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
        }
    ).to_csv(cells_path, index=False)


def _load_completed_run(
    latent_src: Path,
    latent_tgt: Path,
    model_dir: Path,
    manifest_path: Path,
) -> dict[str, Any] | None:
    """Resume from a finished CUDA run. Never loads historical Mac/MPS artifacts."""
    if not (
        latent_src.exists()
        and latent_tgt.exists()
        and manifest_path.exists()
        and (model_dir / "model.pt").exists()
    ):
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    z_s = np.load(latent_src)
    z_t = np.load(latent_tgt)
    if z_s.ndim != 2 or z_t.ndim != 2 or z_s.shape[1] != z_t.shape[1]:
        return None
    if not np.isfinite(z_s).all() or not np.isfinite(z_t).all():
        return None
    return {
        "source_latent": z_s,
        "target_latent": z_t,
        "manifest": manifest,
        "tag": latent_src.name.replace("_source_latent.npy", ""),
        "skipped_retrain": True,
    }


def run_one_model(
    *,
    model_name: str,
    source: ad.AnnData,
    target: ad.AnnData,
    direction: str,
    seed: int,
    cfg: dict[str, Any],
    logger,
) -> dict[str, Any]:
    from scvi.model import SCVI, TOTALVI

    paths = _dir_paths(direction)
    tag = f"{model_name}_{direction}_seed{seed}_cuda"
    model_dir = paths["models"] / model_name / f"seed{seed}_cuda"
    emb_dir = paths["embeddings"]
    latent_src_path = emb_dir / f"{tag}_source_latent.npy"
    latent_tgt_path = emb_dir / f"{tag}_target_latent.npy"
    max_epochs = int(cfg["models"]["max_epochs"])
    batch_size = int(cfg["models"]["batch_size"])
    query_epochs = int(cfg.get("phase10", {}).get("query_max_epochs", 200))

    completed = _load_completed_run(latent_src_path, latent_tgt_path, model_dir, paths["logs"] / f"{tag}_manifest.json")
    if completed is not None:
        logger.info("Reusing completed PHASE 10 artifacts for %s", tag)
        return completed

    set_global_seed(seed)
    if model_name == "scvi_matched":
        hparams = deepcopy(MATCHED_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
        hparams["batch_size"] = batch_size
        src = train_scvi_source(source, hparams, model_dir, logger)
        q = query_adapt(SCVI, target, source, model_dir, query_epochs, batch_size, None, logger)
        extra: dict[str, Any] = {}
    elif model_name == "totalvi":
        hparams = deepcopy(DEFAULT_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["batch_size"] = batch_size
        src = train_totalvi_source(source, hparams, model_dir, logger)
        q = query_adapt(
            TOTALVI,
            target,
            source,
            model_dir,
            query_epochs,
            batch_size,
            {"lr": DEFAULT_HYPERPARAMETERS["lr"], "reduce_lr_on_plateau": True},
            logger,
        )
        extra = {}
    else:
        raise ValueError(model_name)

    if src["latent"].shape[1] != q["latent"].shape[1]:
        raise RuntimeError("Source and query latent dimensions differ; spaces are not comparable.")
    # Primary source embedding is the adapted-model one: same coordinate system
    # as the target embedding. The source-model embedding is kept for the drift
    # diagnostic only and is never used in a cross-space metric.
    source_latent = q["source_latent_adapted"]
    if source_latent.shape[1] != q["latent"].shape[1]:
        raise RuntimeError("Adapted source and query latent dimensions differ.")
    drift = _latent_drift(src["latent"], source_latent)
    logger.info(
        "%s adapted-source drift: mean_shift=%.4f relative=%.4f",
        tag,
        drift["mean_euclidean_shift"],
        drift["relative_shift"],
    )
    _save_latent(latent_src_path, source_latent, source)
    _save_latent(latent_tgt_path, q["latent"], target)
    np.save(emb_dir / f"{tag}_source_latent_source_model.npy", src["latent"])
    if "history" in src and len(src["history"]):
        src["history"].to_csv(model_dir / "source_training_history.csv", index=False)
    if "history" in q and len(q["history"]):
        q["history"].to_csv(model_dir / "query_training_history.csv", index=False)
    manifest = {
        "model": model_name,
        "direction": direction,
        "model_seed": seed,
        "cuda_name_suffix": "_cuda",
        "historical_artifact_overwritten": False,
        "query_api": "prepare_query_anndata + load_query_data",
        "procedure": "cross-dataset transfer / query generalization (scArches)",
        "zero_shot": False,
        "target_labels_used_in_fitting": False,
        "source_n": int(source.n_obs),
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
        "latent_alignment": (
            "Source and target cells are both encoded by the post-adaptation "
            "model, so they share one coordinate system (documented scvi-tools "
            "scArches procedure)."
        ),
        "adapted_source_drift_vs_source_model": drift,
        "source_latent_source_model_path": str(emb_dir / f"{tag}_source_latent_source_model.npy"),
        "source_latent_path": str(latent_src_path),
        "target_latent_path": str(latent_tgt_path),
        "model_dir": str(model_dir),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        **extra,
    }
    save_json(manifest, paths["logs"] / f"{tag}_manifest.json")
    logger.info(
        "%s %s seed%s CUDA source_epochs=%s query_epochs=%s source_s=%.1f query_s=%.1f peak_alloc=%.1f",
        model_name,
        direction,
        seed,
        src["epochs"],
        q["epochs"],
        src["runtime_seconds"],
        q["runtime_seconds"],
        q["gpu"].get("peak_cuda_allocated_mb") or float("nan"),
    )
    return {
        "source_latent": src["latent"],
        "target_latent": q["latent"],
        "manifest": manifest,
        "tag": tag,
    }


def evaluate_direction(
    direction: str,
    source: ad.AnnData,
    target: ad.AnnData,
    fitted: dict[str, dict[str, Any]],
    plans: dict[str, dict[str, Any]],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary_rows = []
    delta_rows = []
    cell_rows = []
    neigh_rows = []
    gap_rows = []
    confusion_frames = []
    for label_level, col in (("l2", "cell_type_l2"), ("l1", "cell_type_l1")):
        plan = plans[label_level]
        classes = plan["eligible_shared_classes"]
        si = plan["source_eval_index"]
        ti = plan["target_eval_index"]
        y_src = source.obs[col].astype(str).to_numpy()[si]
        y_tgt = target.obs[col].astype(str).to_numpy()[ti]
        src_counts = pd.Series(y_src).value_counts()
        tgt_counts = pd.Series(y_tgt).value_counts()
        preds = {}
        for model_name, pack in fitted.items():
            z_s = pack["source_latent"][si]
            z_t = pack["target_latent"][ti]
            clf = classify_transfer(z_s, y_src, z_t, y_tgt, seed)
            preds[model_name] = clf
            man = pack["manifest"]
            row_base_n = len(primary_rows)
            for classifier, prefix in (("logreg", "logreg"), ("knn15", "knn15"), ("knn30", "knn30")):
                primary_rows.append(
                    {
                        "direction": direction,
                        "model": model_name,
                        "model_seed": seed,
                        "source_n": int(source.n_obs),
                        "target_n": int(target.n_obs),
                        "label_level": label_level,
                        "n_shared_classes": plan["n_shared_classes"],
                        "target_coverage": plan["target_coverage"],
                        "classifier": classifier,
                        "accuracy": clf[f"{prefix}_accuracy"],
                        "balanced_accuracy": clf[f"{prefix}_balanced_accuracy"],
                        "macro_f1": clf[f"{prefix}_macro_f1"],
                        "weighted_f1": clf[f"{prefix}_weighted_f1"],
                        "query_adaptation_epochs": man["query_adaptation_epochs"],
                        "source_training_epochs": man["source_training_epochs"],
                        "source_training_runtime": man["source_training_runtime"],
                        "query_runtime": man["query_runtime"],
                        "gpu_peak_allocated_mb": man["source_gpu"].get("peak_cuda_allocated_mb"),
                        "gpu_peak_reserved_mb": man["source_gpu"].get("peak_cuda_reserved_mb"),
                        "query_gpu_peak_allocated_mb": man["query_gpu"].get("peak_cuda_allocated_mb"),
                        "process_peak_rss_mb": man["query_ram_after"].get("process_rss_mb"),
                        "n_eval_source": plan["n_source_in_eval"],
                        "n_eval_target": plan["n_target_in_eval"],
                    }
                )
            spec = per_class_table(y_tgt, clf["logreg_pred"], classes, src_counts, tgt_counts)
            spec["direction"] = direction
            spec["model"] = model_name
            spec["seed"] = seed
            spec["label_level"] = label_level
            spec["subset"] = "high_confidence_shared"
            cell_rows.append(spec)
            labels = classes
            mat = confusion_matrix(y_tgt, clf["logreg_pred"], labels=labels, normalize="true")
            cm = pd.DataFrame(mat, index=labels, columns=labels)
            long = cm.reset_index().melt(id_vars="index", var_name="predicted", value_name="normalized")
            long = long.rename(columns={"index": "true"})
            long["direction"] = direction
            long["model"] = model_name
            long["seed"] = seed
            long["label_level"] = label_level
            confusion_frames.append(long)
            for k in KNN_KS:
                agr = neighbor_agreement(z_s, y_src, z_t, y_tgt, k)
                neigh_rows.append(
                    {
                        "direction": direction,
                        "model": model_name,
                        "seed": seed,
                        "label_level": label_level,
                        "k": k,
                        "mean_agreement": agr["mean_agreement"],
                        "median_agreement": agr["median_agreement"],
                    }
                )
            z_t_all_high = pack["target_latent"]
            high_t = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
            y_all = target.obs[col].astype(str).to_numpy()
            keep = high_t & ~pd.isna(target.obs[col])
            if keep.sum() > 20 and pd.Series(y_all[keep]).nunique() > 1:
                asw = float(silhouette_score(z_t_all_high[keep], y_all[keep], metric="euclidean"))
                pur = neighborhood_purity_per_cell(z_t_all_high[keep], y_all[keep], 15).mean()
                clusters = leiden_clusters(z_t_all_high[keep], PRIMARY_RESOLUTION, seed=seed)
                primary_rows[row_base_n]["target_asw"] = asw
                primary_rows[row_base_n]["target_knn_purity_15"] = float(pur)
                primary_rows[row_base_n]["target_leiden_ari"] = float(adjusted_rand_score(y_all[keep], clusters))
                primary_rows[row_base_n]["target_leiden_nmi"] = float(normalized_mutual_info_score(y_all[keep], clusters))
            z_joint = np.vstack([pack["source_latent"], pack["target_latent"]])
            dataset_ids = np.array(
                ["source"] * pack["source_latent"].shape[0] + ["target"] * pack["target_latent"].shape[0]
            )
            primary_rows[row_base_n]["dataset_asw"] = float(batch_silhouette(z_joint, dataset_ids))
            primary_rows[row_base_n]["latent_spaces_aligned_by"] = (
                "scArches: source encoded by the source model; query encoded by the adapted model"
            )
            within_map = WITHIN_L2_HIGH if label_level == "l2" else WITHIN_L1_HIGH
            key = "scVI_matched" if model_name == "scvi_matched" else "totalVI"
            within = within_map[key]
            gap_rows.append(
                {
                    "model": model_name,
                    "direction": direction,
                    "label_level": label_level,
                    "within_dataset_macro_f1": within,
                    "cross_dataset_macro_f1": clf["logreg_macro_f1"],
                    "generalization_gap": float(clf["logreg_macro_f1"] - within),
                    "within_dataset_reference": WITHIN_NOTE,
                    "model_seed": seed,
                }
            )
        boot = paired_bootstrap_delta(
            y_tgt,
            preds["scvi_matched"]["logreg_pred"],
            preds["totalvi"]["logreg_pred"],
            seed=seed,
        )
        delta_rows.append(
            {
                "direction": direction,
                "model_seed": seed,
                "label_level": label_level,
                "scvi_matched_macro_f1": preds["scvi_matched"]["logreg_macro_f1"],
                "totalvi_macro_f1": preds["totalvi"]["logreg_macro_f1"],
                "delta_totalvi_minus_scvi": float(
                    preds["totalvi"]["logreg_macro_f1"] - preds["scvi_matched"]["logreg_macro_f1"]
                ),
                "bootstrap_ci_low": boot["bootstrap_ci_low"],
                "bootstrap_ci_high": boot["bootstrap_ci_high"],
                "bootstrap_delta_mean": boot["delta_mean"],
                "target_n": plan["n_target_in_eval"],
                "shared_class_n": plan["n_shared_classes"],
                "target_coverage": plan["target_coverage"],
                "classifier": "logreg",
            }
        )
    primary = pd.DataFrame(primary_rows)
    delta = pd.DataFrame(delta_rows)
    cell = pd.concat(cell_rows, ignore_index=True)
    neigh = pd.DataFrame(neigh_rows)
    gap = pd.DataFrame(gap_rows)
    confusion = pd.concat(confusion_frames, ignore_index=True)
    return primary, delta, cell, neigh, gap, confusion


def run_direction(direction: str, source_batch: str, target_batch: str, seed: int, cfg: dict[str, Any], logger) -> dict[str, Any]:
    source, target, audit = load_source_target(source_batch, target_batch)
    paths = _dir_paths(direction)
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    plans = {
        "l2": class_plan(source, target, "cell_type_l2"),
        "l1": class_plan(source, target, "cell_type_l1"),
    }
    audit_out = {
        **audit,
        "direction": direction,
        "model_seed": seed,
        "l2_plan": {k: v for k, v in plans["l2"].items() if k not in {"source_eval_index", "target_eval_index"}},
        "l1_plan": {k: v for k, v in plans["l1"].items() if k not in {"source_eval_index", "target_eval_index"}},
        "min_source_support_pre_specified": MIN_SOURCE_SUPPORT,
        "min_target_support_pre_specified": MIN_TARGET_SUPPORT,
        "phase9_note": (
            "PHASE 9 was not empirically executed because native totalVI 1.3.3 does not "
            "support arbitrary cell-level missing-protein masking without conflating "
            "missingness with observed zeros."
        ),
    }
    save_json(audit_out, paths["logs"] / f"{direction}_seed{seed}_dataset_audit.json")
    logger.info(
        "%s audit source=%s target=%s l2 shared classes=%s coverage=%.3f",
        direction,
        audit["source_n"],
        audit["target_n"],
        plans["l2"]["n_shared_classes"],
        plans["l2"]["target_coverage"],
    )
    fitted = {}
    for model_name in ("scvi_matched", "totalvi"):
        fitted[model_name] = run_one_model(
            model_name=model_name,
            source=source,
            target=target,
            direction=direction,
            seed=seed,
            cfg=cfg,
            logger=logger,
        )
        _cleanup()
    primary, delta, cell, neigh, gap, confusion = evaluate_direction(
        direction, source, target, fitted, plans, seed
    )
    del source, target, fitted
    _cleanup()
    return {
        "primary": primary,
        "delta": delta,
        "cell": cell,
        "neigh": neigh,
        "gap": gap,
        "confusion": confusion,
        "audit": audit_out,
    }


def _append(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        old = pd.read_csv(path)
        frame = pd.concat([old, frame], ignore_index=True)
        frame = frame.drop_duplicates(
            subset=[c for c in ("direction", "model", "model_seed", "label_level", "classifier", "cell_type", "k", "true", "predicted") if c in frame.columns],
            keep="last",
        )
    frame.to_csv(path, index=False)


TABLE_KEYS = {
    "cross_dataset_primary.csv": ["direction", "model", "model_seed", "label_level", "classifier"],
    "cross_dataset_primary_delta.csv": ["direction", "model_seed", "label_level", "classifier"],
    "cross_dataset_celltype_specific.csv": ["direction", "model", "seed", "label_level", "cell_type"],
    "cross_dataset_neighbor_agreement.csv": ["direction", "model", "seed", "label_level", "k"],
    "cross_dataset_generalization_gap.csv": ["model", "direction", "label_level", "model_seed"],
    "cross_dataset_confusion.csv": ["direction", "model", "seed", "label_level", "true", "predicted"],
}


def _merge_existing(path: Path, frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if path.exists():
        old = pd.read_csv(path)
        frame = pd.concat([old, frame], ignore_index=True)
        present = [c for c in keys if c in frame.columns]
        if present:
            frame = frame.drop_duplicates(subset=present, keep="last")
    return frame


def write_tables(pieces: list[dict[str, Any]]) -> dict[str, Path]:
    tables = _phase10_root() / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    mapping = {
        "cross_dataset_primary.csv": pd.concat([p["primary"] for p in pieces], ignore_index=True),
        "cross_dataset_primary_delta.csv": pd.concat([p["delta"] for p in pieces], ignore_index=True),
        "cross_dataset_celltype_specific.csv": pd.concat([p["cell"] for p in pieces], ignore_index=True),
        "cross_dataset_neighbor_agreement.csv": pd.concat([p["neigh"] for p in pieces], ignore_index=True),
        "cross_dataset_generalization_gap.csv": pd.concat([p["gap"] for p in pieces], ignore_index=True),
        "cross_dataset_confusion.csv": pd.concat([p["confusion"] for p in pieces], ignore_index=True),
    }
    paths = {}
    for name, frame in mapping.items():
        out = tables / name
        frame = _merge_existing(out, frame, TABLE_KEYS[name])
        frame.to_csv(out, index=False)
        paths[name] = out
    delta = pd.read_csv(tables / "cross_dataset_primary_delta.csv")
    if delta["model_seed"].nunique() > 1:
        summary = (
            delta.groupby(["direction", "label_level"], as_index=False)
            .agg(
                n_seeds=("model_seed", "nunique"),
                scvi_matched_macro_f1_mean=("scvi_matched_macro_f1", "mean"),
                scvi_matched_macro_f1_sd=("scvi_matched_macro_f1", "std"),
                totalvi_macro_f1_mean=("totalvi_macro_f1", "mean"),
                totalvi_macro_f1_sd=("totalvi_macro_f1", "std"),
                delta_mean=("delta_totalvi_minus_scvi", "mean"),
                delta_sd=("delta_totalvi_minus_scvi", "std"),
            )
        )
        summary.to_csv(tables / "cross_dataset_primary_delta_summary.csv", index=False)
    return paths


def run_phase10(*, seeds: list[int], directions: tuple = DIRECTIONS, config: dict[str, Any] | None = None) -> dict[str, Any]:
    tmp = Path.home() / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmp)
    cfg = config or load_config()
    logs = _phase10_root() / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    logger = get_logger("phase10", logs / "phase10.log")
    env = detect_compute_environment()
    logger.info("PHASE 10 CUDA env: %s torch=%s scvi=%s gpu=%s ram=%s", env.get("hostname"), env.get("torch_version"), env.get("scvi_version"), env.get("gpu_names"), _ram())
    _require_cuda()
    if str(env.get("scvi_version")) != "1.3.3":
        raise RuntimeError(f"Expected scvi-tools 1.3.3, found {env.get('scvi_version')}")
    pieces = []
    for seed in seeds:
        set_global_seed(seed)
        for direction, src, tgt in directions:
            logger.info("=== %s seed %s ===", direction, seed)
            pieces.append(run_direction(direction, src, tgt, seed, cfg, logger))
    paths = write_tables(pieces)
    write_phase10_figures(_phase10_root() / "tables", _phase10_root() / "figures")
    summary = {
        "seeds": seeds,
        "directions": [d[0] for d in directions],
        "tables": {k: str(v) for k, v in paths.items()},
        "cuda": True,
        "query_procedure": "scArches prepare_query_anndata + load_query_data",
        "zero_shot": False,
        "phase9_not_run": True,
        "historical_mps_untouched": True,
        "ram": _ram(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    save_json(summary, logs / f"phase10_summary_seeds{'_'.join(map(str, seeds))}.json")
    return summary
