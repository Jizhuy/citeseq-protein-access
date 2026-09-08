"""Train totalVI on a protein-corrupted copy. Never writes the original inner h5ad."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb
from src.models.run_scvi import (
    counts_to_csr,
    epochs_completed_from_history,
    history_step_frame,
    history_to_frame,
    requested_accelerator,
    _final_metric,
)
from src.models.run_totalvi import (
    DEFAULT_HYPERPARAMETERS,
    PROTEIN_OBSM_KEY,
    train_totalvi_with_device_fallback,
)
from src.utils.device import detect_compute_environment, lightning_devices
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed


def train_corrupted_totalvi(
    adata_corrupted: ad.AnnData,
    *,
    condition_dir: Path,
    embeddings_dir: Path,
    fraction: float,
    perturbation_seed: int,
    model_seed: int = 0,
    config: dict[str, Any] | None = None,
    logger=None,
) -> dict[str, Any]:
    """Train PHASE-4 totalVI on corrupted proteins. RNA and labels are unchanged."""
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    logger = logger or get_logger(
        "phase8_totalvi",
        logs_dir / f"phase8_totalvi_{int(round(fraction * 100)):03d}_seed{perturbation_seed}.log",
    )
    set_global_seed(model_seed)
    env = detect_compute_environment()
    adata, _csr = counts_to_csr(adata_corrupted)
    protein = adata.obsm[PROTEIN_OBSM_KEY]
    if not isinstance(protein, pd.DataFrame):
        raise TypeError("corrupted protein_expression must remain a DataFrame")
    protein_names = [str(c) for c in protein.columns]
    original_order = adata.obs_names.copy()
    original_split = adata.obs["split"].astype(str).copy()
    rna_checksum = float(
        adata.layers["counts"].sum()
        if hasattr(adata.layers["counts"], "sum")
        else np.asarray(adata.layers["counts"]).sum()
    )

    hparams = deepcopy(DEFAULT_HYPERPARAMETERS)
    hparams["n_latent"] = int(cfg["models"]["latent_dim"])
    hparams["max_epochs"] = int(cfg["models"]["max_epochs"])
    hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
    hparams["batch_size"] = int(cfg["models"]["batch_size"])

    requested_device, requested_detail = requested_accelerator(
        env,
        prefer_cuda=bool(cfg["device"].get("prefer_cuda", True)),
        prefer_mps=bool(cfg["device"].get("prefer_mps", True)),
        allow_cpu_fallback=bool(cfg["device"].get("allow_cpu_fallback", cfg["device"].get("fallback_cpu", True))),
    )
    condition_dir.mkdir(parents=True, exist_ok=True)
    embeddings_dir.mkdir(parents=True, exist_ok=True)
    timer = Timer()
    timer.start()
    model, actual_device, train_warnings, fallback_error = train_totalvi_with_device_fallback(
        adata, hparams, requested_device, logger, devices=lightning_devices(cfg, requested_device)
    )
    runtime = timer.stop()
    gpu_mem = cuda_memory_stats()
    model.save(str(condition_dir), overwrite=True, save_anndata=False)

    latent = np.asarray(model.get_latent_representation())
    if latent.shape != (adata.n_obs, hparams["n_latent"]):
        raise RuntimeError(f"Unexpected latent shape {latent.shape}")
    if not np.isfinite(latent).all():
        raise RuntimeError("totalVI latent contains NaN or Inf.")
    if not np.array_equal(adata.obs_names, original_order):
        raise RuntimeError("Cell order changed during corrupted totalVI training.")
    if not adata.obs["split"].astype(str).equals(original_split):
        raise RuntimeError("Train/test split changed during corrupted totalVI training.")

    current_rna = float(adata.layers["counts"].sum())
    if abs(current_rna - rna_checksum) > 1e-3:
        raise RuntimeError("RNA counts changed during corrupted totalVI training.")

    tag = f"corruption_{int(round(fraction * 100)):03d}_seed{perturbation_seed}"
    np.save(embeddings_dir / f"{tag}_latent.npy", latent)
    pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "split": adata.obs["split"].astype(str).to_numpy(),
        }
    ).to_csv(embeddings_dir / f"{tag}_cells.csv", index=False)

    rna_norm, protein_norm = model.get_normalized_expression(
        return_mean=True, include_protein_background=False, scale_protein=False
    )
    del rna_norm
    protein_norm.to_csv(embeddings_dir / f"{tag}_normalized_protein.csv", index_label="cell_id")
    foreground = model.get_protein_foreground_probability(return_mean=True)
    foreground.to_csv(
        embeddings_dir / f"{tag}_protein_foreground_probability.csv", index_label="cell_id"
    )

    history = history_to_frame(model.history)
    step_history = history_step_frame(model.history)
    epochs_completed = epochs_completed_from_history(history, hparams["max_epochs"])
    history.to_csv(condition_dir / "training_history.csv", index=False)
    if len(step_history):
        step_history.to_csv(condition_dir / "training_history_step.csv", index=False)

    final_metrics = {
        column: _final_metric(history, column)
        for column in history.columns
        if column not in {"epoch", "step"}
    }
    info = {
        "model": "totalVI_corrupted",
        "experiment": "protein_sparsity_stress_test",
        "not_missing_modality": True,
        "corruption_fraction": float(fraction),
        "perturbation_seed": int(perturbation_seed),
        "model_seed": int(model_seed),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_proteins": int(len(protein_names)),
        "protein_names": protein_names,
        "runtime_seconds": float(runtime),
        "epochs": int(epochs_completed),
        "max_epochs": int(hparams["max_epochs"]),
        "early_stopped": bool(epochs_completed < hparams["max_epochs"]),
        "device": actual_device,
        "requested_device": requested_device,
        "requested_device_detail": requested_detail,
        "device_fallback_error": fallback_error,
        "gpu_name": (env.get("gpu_names") or [None])[0],
        "cuda_version": env.get("cuda_version"),
        "driver_version": env.get("driver_version"),
        "peak_rss_mb": float(peak_rss_mb()),
        "peak_cuda_allocated_mb": gpu_mem.get("peak_cuda_allocated_mb"),
        "peak_cuda_reserved_mb": gpu_mem.get("peak_cuda_reserved_mb"),
        "final_logged_metrics": {k: v for k, v in final_metrics.items() if v is not None},
        "warnings": train_warnings,
        "reused_phase4_model": False,
        "hparams": hparams,
        "scvi_tools_version": env["scvi_version"],
        "torch_version": env["torch_version"],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model_path": str(condition_dir),
        "latent_path": str(embeddings_dir / f"{tag}_latent.npy"),
    }
    save_json(info, condition_dir / "manifest.json")
    return info
