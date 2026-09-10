#!/usr/bin/env python
"""PHASE 11B: uncertainty robustness, posterior recovery, variance decomposition.

Prospective replication that fixes the PHASE 10 artifact-design problem. Every
run persists BOTH the source reference model and the scArches-adapted query
model, so target-cell posterior uncertainty is recoverable afterwards.

Three hard rules enforced in code:
  1. Downstream evaluation always uses the POST-adaptation source embedding.
     The PHASE 10 bug returned the pre-adaptation array; `assert_post_adaptation`
     is the regression guard against it.
  2. Target labels never touch training, adaptation, early stopping, or the
     classifier fit.
  3. Nothing is written outside results/phase11b_uncertainty_robustness/.
"""

from __future__ import annotations

import gc
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import psutil
import torch

from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.experiments.phase10_cross_dataset import (  # noqa: F401
    _apply_run_seed,
    _train_kwargs,
    ensure_csr,
    train_scvi_source,
    train_totalvi_source,
)
from src.models.run_scvi import epochs_completed_from_history, history_to_frame
from src.models.run_totalvi import DEFAULT_HYPERPARAMETERS, PROTEIN_OBSM_KEY
from src.experiments.phase11b_stats import assert_post_adaptation, latent_drift
from src.models.run_scvi_matched import MATCHED_HYPERPARAMETERS
from src.utils.io import resolve_path, save_json

__all__ = ["assert_post_adaptation", "latent_drift"]

DIRECTIONS = (
    ("direction_A", "pbmc10k_to_pbmc5k", "PBMC10k", "PBMC5k"),
    ("direction_B", "pbmc5k_to_pbmc10k", "PBMC5k", "PBMC10k"),
)
MODELS = ("scvi_matched", "totalvi")
N_SEEDS = 20
SEEDS = tuple(range(N_SEEDS))
CROSSED_SUBSET_SEEDS = (0, 1, 2, 3, 4)
CROSSED_MODEL_SEEDS = (0, 1, 2, 3, 4)
MC_CHECK_SEEDS = (0, 5, 10, 15)
MC_SAMPLES = 50
INTEGRITY_SUBSET = 256
INTEGRITY_TOL = 1e-5
NEIGHBOR_K = 15
ECE_BINS = 15
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
N_BOOT_CELLS = 500
N_BOOT_SEEDS = 10_000
BOOTSTRAP_SEED = 11_000


def phase11b_root() -> Path:
    return resolve_path("results/phase11b_uncertainty_robustness")


def phase11b_paths() -> dict[str, Path]:
    root = phase11b_root()
    return {
        "root": root,
        "models": root / "models",
        "embeddings": root / "embeddings",
        "probabilities": root / "probabilities",
        "posterior": root / "posterior",
        "tables": root / "tables",
        "figures": root / "figures",
        "logs": root / "logs",
        "done": root / "logs" / "done",
    }


def ensure_dirs() -> dict[str, Path]:
    paths = phase11b_paths()
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def experiment_tag(direction_key: str, model: str, model_seed: int, subset_seed: int | None,
                   replicate: int | None = None) -> str:
    if subset_seed is None:
        tag = f"{direction_key}__{model}__seed{model_seed:02d}"
    else:
        tag = f"{direction_key}__{model}__subset{subset_seed}__seed{model_seed:02d}"
    return tag if replicate is None else f"{tag}__rep{replicate:02d}"


def model_dir_for(paths: dict[str, Path], direction_key: str, model: str,
                  model_seed: int, subset_seed: int | None,
                  replicate: int | None = None) -> Path:
    base = paths["models"] / direction_key / model
    if subset_seed is None:
        base = base / f"seed{model_seed:02d}"
    else:
        base = base / f"subset{subset_seed}_seed{model_seed:02d}"
    return base if replicate is None else base / f"replicate{replicate:02d}"


def ram_snapshot() -> dict[str, float]:
    vm = psutil.virtual_memory()
    return {
        "available_mb": float(vm.available) / 1024**2,
        "percent": float(vm.percent),
        "process_rss_mb": float(peak_rss_mb()),
    }


STORE_ATTRS = ("_setup_adata_manager_store", "_per_instance_manager_store")


def purge_scvi_manager_store() -> int:
    """Drop scvi-tools' class-level AnnDataManager registries.

    `setup_anndata`, `load_query_data`, `load`, and every `_validate_anndata`
    call on a fresh object registers an AnnDataManager in two class-level dicts.
    Each manager holds its AnnData, so deleting the model frees nothing: one run
    of this pipeline strands three full copies. Across a 130-run grid on a
    7.6 GB machine that exhausts RAM and pushes the process into swap.

    Every model subclass carries its OWN pair of dicts rather than sharing the
    ones on BaseModelClass, so the whole subclass tree has to be walked; clearing
    only BaseModelClass is a no-op.

    Safe only when no model is live, which is why this runs between runs and
    never inside one. Returns the number of entries dropped.
    """
    try:
        from scvi.model.base import BaseModelClass
    except ImportError:
        return 0

    seen: set[int] = set()
    classes: list[type] = []

    def walk(cls: type) -> None:
        if id(cls) in seen:
            return
        seen.add(id(cls))
        classes.append(cls)
        for sub in cls.__subclasses__():
            walk(sub)

    walk(BaseModelClass)
    dropped = 0
    for cls in classes:
        for name in STORE_ATTRS:
            store = cls.__dict__.get(name)
            if isinstance(store, dict) and store:
                dropped += len(store)
                store.clear()
    return dropped


def cleanup(purge_manager_store: bool = False) -> None:
    """Free memory between runs. Pass purge_manager_store only when no model is live."""
    if purge_manager_store:
        purge_scvi_manager_store()
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ---------------------------------------------------------------------------
# Query adaptation that PERSISTS the adapted model
# ---------------------------------------------------------------------------
def query_adapt_and_save(
    model_cls,
    target,
    source,
    reference_dir: Path,
    adapted_dir: Path,
    max_epochs: int,
    batch_size: int,
    extra_train: dict | None,
    tag: str,
    logger,
    run_seed: int | None = None,
) -> dict[str, Any]:
    """scArches adaptation; saves the adapted model and verifies it reloads.

    ``run_seed`` (PHASE 12A) redraws the run-level randomness of the adaptation
    step -- the new batch-specific weight initialisation plus the split,
    shuffling and dropout of the adaptation optimiser -- while leaving the
    reference model on disk untouched. ``None`` is the PHASE 11B path.
    """
    reset_cuda_peak_stats()
    query = target.copy()
    model_cls.prepare_query_anndata(query, str(reference_dir))
    _apply_run_seed(run_seed)
    qmodel = model_cls.load_query_data(query, str(reference_dir), accelerator="gpu", device="auto")
    trainable = [n for n, p in qmodel.module.named_parameters() if p.requires_grad]

    extra: dict[str, Any] = {"plan_kwargs": {"weight_decay": 0.0}}
    if extra_train:
        extra.update(extra_train)
    timer = Timer()
    timer.start()
    qmodel.train(**_train_kwargs(max_epochs, batch_size, True, extra))
    runtime = float(timer.stop())

    device = str(next(qmodel.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"{tag}: query adaptation did not use CUDA; device={device}")

    # Encode BOTH sets with the adapted model: this is the shared coordinate system.
    target_latent = np.asarray(qmodel.get_latent_representation())
    source_ref = ensure_csr(source.copy())
    source_latent_post = np.asarray(qmodel.get_latent_representation(source_ref))
    for name, arr in (("target", target_latent), ("source", source_latent_post)):
        if not np.isfinite(arr).all():
            raise RuntimeError(f"{tag}: adapted {name} latent contains NaN/Inf")

    # Posterior q(z|x) for the target under the ADAPTED model.
    tgt_mu, tgt_var = qmodel.get_latent_representation(return_dist=True)
    src_mu, src_var = qmodel.get_latent_representation(source_ref, return_dist=True)
    if tgt_var.shape != (target.n_obs, 20):
        raise RuntimeError(f"{tag}: target posterior variance shape {tgt_var.shape}, expected "
                           f"({target.n_obs}, 20)")
    if not np.all(tgt_var > 0):
        raise RuntimeError(f"{tag}: non-positive target posterior variance")

    history = history_to_frame(qmodel.history)
    epochs = epochs_completed_from_history(history, max_epochs)
    gpu = cuda_memory_stats()

    # ---- persist the adapted model, then prove it reloads faithfully (§8) ----
    if adapted_dir.exists():
        shutil.rmtree(adapted_dir)
    adapted_dir.mkdir(parents=True, exist_ok=True)
    qmodel.save(str(adapted_dir), overwrite=True, save_anndata=False)
    n_check = min(INTEGRITY_SUBSET, query.n_obs)
    check_idx = np.arange(n_check)
    pre_save_latent = np.asarray(qmodel.get_latent_representation(indices=check_idx))
    pre_save_mu, pre_save_var = qmodel.get_latent_representation(indices=check_idx, return_dist=True)

    del qmodel
    cleanup()

    reloaded = model_cls.load(str(adapted_dir), adata=query, accelerator="gpu", device="auto")
    post_load_latent = np.asarray(reloaded.get_latent_representation(indices=check_idx))
    post_load_mu, post_load_var = reloaded.get_latent_representation(indices=check_idx, return_dist=True)
    max_latent_diff = float(np.abs(post_load_latent - pre_save_latent).max())
    max_var_diff = float(np.abs(post_load_var - pre_save_var).max())
    integrity_ok = max_latent_diff <= INTEGRITY_TOL and max_var_diff <= INTEGRITY_TOL
    if not integrity_ok:
        raise RuntimeError(
            f"{tag}: adapted model reload changed outputs "
            f"(latent {max_latent_diff:.3e}, var {max_var_diff:.3e} > {INTEGRITY_TOL})"
        )
    del reloaded, source_ref, query
    cleanup()

    return {
        "target_latent": target_latent,
        "source_latent_post": source_latent_post,
        "target_mu": tgt_mu,
        "target_var": tgt_var,
        "source_mu": src_mu,
        "source_var": src_var,
        "runtime_seconds": runtime,
        "epochs": int(epochs),
        "device": device,
        "gpu": gpu,
        "n_trainable_params": len(trainable),
        "trainable_param_names": trainable,
        "integrity_max_latent_diff": max_latent_diff,
        "integrity_max_variance_diff": max_var_diff,
        "integrity_n_cells_checked": int(n_check),
        "integrity_ok": bool(integrity_ok),
        "history": history,
    }


# ---------------------------------------------------------------------------
# One full experiment
# ---------------------------------------------------------------------------
def run_experiment(
    *,
    direction_key: str,
    direction: str,
    model_name: str,
    model_seed: int,
    source,
    target,
    cfg: dict[str, Any],
    paths: dict[str, Path],
    logger,
    subset_seed: int | None = None,
    replicate: int | None = None,
    run_replicate_seed: int | None = None,
) -> dict[str, Any] | None:
    """Train, adapt, persist, verify, and extract. Returns None if already done.

    PHASE 12A adds ``replicate`` (an index that only widens the tag and output
    directory) and ``run_replicate_seed``. ``model_seed`` always fixes the
    source-model weight initialisation; ``run_replicate_seed``, when given,
    redraws every downstream stochastic choice. Both default to None, which
    reproduces PHASE 11B byte for byte.
    """
    from scvi.model import SCVI, TOTALVI
    from src.utils.seed import set_global_seed

    tag = experiment_tag(direction_key, model_name, model_seed, subset_seed, replicate)
    done_marker = paths["done"] / f"{tag}.json"
    if done_marker.exists():
        try:
            marker = json.loads(done_marker.read_text(encoding="utf-8"))
            if marker.get("complete"):
                logger.info("skip (done): %s", tag)
                return None
        except json.JSONDecodeError:
            logger.warning("corrupt done marker for %s; rerunning", tag)

    if not torch.cuda.is_available():
        raise RuntimeError("PHASE 11B requires CUDA.")

    mdir = model_dir_for(paths, direction_key, model_name, model_seed, subset_seed, replicate)
    source_dir = mdir / "source_model"
    adapted_dir = mdir / "adapted_query_model"
    max_epochs = int(cfg["models"]["max_epochs"])
    batch_size = int(cfg["models"]["batch_size"])
    query_epochs = int(cfg.get("phase10", {}).get("query_max_epochs", 200))

    logger.info("=== %s (source n=%s target n=%s) ===", tag, source.n_obs, target.n_obs)
    set_global_seed(model_seed)
    ram_before = ram_snapshot()

    if model_name == "scvi_matched":
        hparams = dict(MATCHED_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
        hparams["batch_size"] = batch_size
        src = train_scvi_source(source, hparams, source_dir, logger, run_seed=run_replicate_seed)
        extra_train = None
        model_cls = SCVI
    elif model_name == "totalvi":
        hparams = dict(DEFAULT_HYPERPARAMETERS)
        hparams["max_epochs"] = max_epochs
        hparams["batch_size"] = batch_size
        src = train_totalvi_source(source, hparams, source_dir, logger, run_seed=run_replicate_seed)
        extra_train = {"lr": DEFAULT_HYPERPARAMETERS["lr"],
                       "reduce_lr_on_plateau": DEFAULT_HYPERPARAMETERS["reduce_lr_on_plateau"]}
        model_cls = TOTALVI
    else:
        raise ValueError(model_name)

    q = query_adapt_and_save(
        model_cls, target, source, source_dir, adapted_dir,
        query_epochs, batch_size, extra_train, tag, logger,
        run_seed=None if run_replicate_seed is None else run_replicate_seed + 1,
    )

    # The evaluator must receive the post-adaptation source embedding (§5).
    source_latent = q["source_latent_post"]
    assert_post_adaptation(source_latent, src["latent"], q["source_latent_post"], tag)
    drift = latent_drift(src["latent"], source_latent)

    emb = paths["embeddings"] / direction_key
    emb.mkdir(parents=True, exist_ok=True)
    np.save(emb / f"{tag}_source_latent_post.npy", source_latent.astype(np.float32))
    np.save(emb / f"{tag}_target_latent_post.npy", q["target_latent"].astype(np.float32))
    np.save(emb / f"{tag}_source_latent_pre.npy", src["latent"].astype(np.float32))
    pd.DataFrame({"cell_id": source.obs_names.astype(str)}).to_csv(
        emb / f"{tag}_source_cells.csv", index=False)
    pd.DataFrame({"cell_id": target.obs_names.astype(str)}).to_csv(
        emb / f"{tag}_target_cells.csv", index=False)

    post = paths["posterior"] / direction_key
    post.mkdir(parents=True, exist_ok=True)
    u_latent_target = q["target_var"].mean(axis=1)
    np.savez_compressed(
        post / f"{tag}_posterior.npz",
        target_mu=q["target_mu"].astype(np.float32),
        target_var=q["target_var"].astype(np.float32),
        source_var=q["source_var"].astype(np.float32),
        u_latent_target=u_latent_target.astype(np.float32),
        target_cell_ids=target.obs_names.astype(str).to_numpy(),
        source_cell_ids=source.obs_names.astype(str).to_numpy(),
    )

    manifest = {
        "phase": "11B",
        "tag": tag,
        "direction_key": direction_key,
        "direction": direction,
        "model": model_name,
        "model_seed": model_seed,
        "subset_seed": subset_seed,
        "replicate": replicate,
        "run_replicate_seed": run_replicate_seed,
        "source_n": int(source.n_obs),
        "target_n": int(target.n_obs),
        "source_training_epochs": src["epochs"],
        "query_epochs": q["epochs"],
        "source_runtime_seconds": src["runtime_seconds"],
        "query_runtime_seconds": q["runtime_seconds"],
        "total_runtime_seconds": src["runtime_seconds"] + q["runtime_seconds"],
        "source_device": src["device"],
        "query_device": q["device"],
        "gpu_peak_allocated_mb": q["gpu"].get("peak_cuda_allocated_mb"),
        "gpu_peak_reserved_mb": q["gpu"].get("peak_cuda_reserved_mb"),
        "ram_before": ram_before,
        "ram_after": ram_snapshot(),
        "n_trainable_query_params": q["n_trainable_params"],
        "trainable_query_params": q["trainable_param_names"],
        "adapted_query_model_saved": True,
        "adapted_query_model_dir": str(adapted_dir),
        "source_model_dir": str(source_dir),
        "artifact_integrity": {
            "n_cells_checked": q["integrity_n_cells_checked"],
            "max_latent_diff": q["integrity_max_latent_diff"],
            "max_variance_diff": q["integrity_max_variance_diff"],
            "tolerance": INTEGRITY_TOL,
            "passed": q["integrity_ok"],
        },
        "source_latent_drift_pre_to_post_adaptation": drift,
        "source_embedding_used_for_metrics": "post_adaptation",
        "target_labels_used_in_fitting": False,
        "posterior_quantity": "get_latent_representation(return_dist=True) -> (qz.loc, qz.scale**2)",
        "median_u_latent_target": float(np.median(u_latent_target)),
        "median_target_posterior_variance": float(np.median(q["target_var"])),
        "max_target_posterior_variance": float(q["target_var"].max()),
        "timestamp_utc": utc_now(),
    }
    save_json(manifest, paths["logs"] / "manifests" / f"{tag}_manifest.json")
    save_json({"complete": True, "tag": tag, "timestamp_utc": utc_now()}, done_marker)
    del src, q, source_latent
    # No model is live at this point, so the manager registries can be purged.
    cleanup(purge_manager_store=True)
    ram_after = ram_snapshot()
    logger.info(
        "%s done: src_ep=%s q_ep=%s runtime=%.0fs drift_mean=%.6f integrity=%.2e "
        "u_latent=%.5f ram_avail=%.0fMB",
        tag, manifest["source_training_epochs"], manifest["query_epochs"],
        manifest["total_runtime_seconds"],
        drift["mean_cell_displacement"], manifest["artifact_integrity"]["max_latent_diff"],
        manifest["median_u_latent_target"], ram_after["available_mb"],
    )
    return manifest
