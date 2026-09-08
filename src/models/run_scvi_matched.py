"""PHASE 5A: RNA-only scVI control matched to totalVI where SCVI allows."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from src.evaluation.plots import plot_training_history, plot_umap_by_batch
from src.evaluation.representation_metrics import compute_representation_metrics
from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb
from src.models.run_scvi import (
    EXPECTED_BATCHES,
    EXPECTED_N_OBS,
    ScviDataError,
    _final_metric,
    _n_indices,
    counts_to_csr,
    epochs_completed_from_history,
    history_to_frame,
    requested_accelerator,
    train_scvi_with_device_fallback,
    validate_rna_adata,
)
from src.utils.device import detect_compute_environment, lightning_devices
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

MATCHED_HYPERPARAMETERS = {
    "n_latent": 20,
    "n_hidden": 256,
    "n_layers": 2,
    "dropout_rate": 0.2,
    "dispersion": "gene",
    "gene_likelihood": "nb",
    "use_observed_lib_size": True,
    "latent_distribution": "normal",
    "train_size": 0.9,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_mode": "min",
}


def run_scvi_matched(config: dict[str, Any] | None = None, seed: int = 0) -> dict[str, Any]:
    """Train scVI_matched. Does not overwrite scVI_default or totalVI artifacts."""
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    models_dir = resolve_path(cfg["paths"].get("models_dir", "results/models"))
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    logger = get_logger("run_scvi_matched", logs_dir / f"scvi_matched_seed{seed}.log")

    seed_info = set_global_seed(seed)
    env = detect_compute_environment()
    logger.info(
        "Python %s | torch %s | scvi-tools %s | accelerator auto-detect=%s",
        env["python"],
        env["torch_version"],
        env["scvi_version"],
        env["accelerator"],
    )
    if str(env["scvi_version"]) != "1.3.3":
        logger.warning("Expected scvi-tools 1.3.3, found %s", env["scvi_version"])

    source_path = processed_dir / "pbmc_cite_combined_inner.h5ad"
    adata_original = ad.read_h5ad(source_path)
    validate_rna_adata(adata_original)
    original_order = adata_original.obs_names.copy()
    original_split = adata_original.obs["split"].astype(str).copy()
    if adata_original.obs["cell_type"].astype(str).nunique() == 1:
        logger.info("obs['cell_type'] is a placeholder and will not be used for metrics.")

    adata, csr_report = counts_to_csr(adata_original)
    sparse_path = processed_dir / "pbmc_cite_combined_inner_sparse.h5ad"
    if not sparse_path.exists():
        adata.write_h5ad(sparse_path, compression="gzip")
        logger.info("Wrote CSR training cache to %s", sparse_path)
    else:
        logger.info("Reusing existing CSR cache without rewriting it.")

    hparams = deepcopy(MATCHED_HYPERPARAMETERS)
    hparams["max_epochs"] = int(cfg["models"]["max_epochs"])
    hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
    hparams["batch_size"] = int(cfg["models"]["batch_size"])
    hparams["gene_selection"] = (
        "all_15792_shared_genes; config n_top_genes=4000 was not applied"
    )
    # Do not read config n_hidden/n_layers; those belong to scVI_default.

    requested_device, requested_detail = requested_accelerator(
        env,
        prefer_cuda=bool(cfg["device"].get("prefer_cuda", True)),
        prefer_mps=bool(cfg["device"].get("prefer_mps", True)),
        allow_cpu_fallback=bool(cfg["device"].get("allow_cpu_fallback", cfg["device"].get("fallback_cpu", True))),
    )
    logger.info("Requested device: %s (%s)", requested_device, requested_detail)
    models_dir.mkdir(parents=True, exist_ok=True)
    model_dir = models_dir / f"scvi_matched_seed{seed}"

    timer = Timer()
    timer.start()
    model, actual_device, train_warnings, fallback_error = train_scvi_with_device_fallback(
        adata, hparams, requested_device, logger, devices=lightning_devices(cfg, requested_device)
    )
    runtime = timer.stop()
    gpu_mem = cuda_memory_stats()
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    logger.info("Saved scVI_matched model to %s immediately after training.", model_dir)

    latent = np.asarray(model.get_latent_representation())
    if latent.shape != (EXPECTED_N_OBS, hparams["n_latent"]):
        raise ScviDataError(
            f"Latent shape {latent.shape} != ({EXPECTED_N_OBS}, {hparams['n_latent']})"
        )
    if not np.array_equal(adata.obs_names, original_order):
        raise ScviDataError("Cell order changed relative to the original AnnData.")
    if not adata.obs["split"].astype(str).equals(original_split):
        raise ScviDataError("obs['split'] was altered during scVI_matched training.")
    if "X_scVI" in adata.obsm or "X_totalVI" in adata.obsm:
        raise ScviDataError("Refusing to write onto an object that already holds other embeddings.")

    adata.obsm["X_scVI_matched"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI_matched", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm["X_scVI_matched_umap"] = np.asarray(adata.obsm["X_umap"])

    metrics = compute_representation_metrics(
        latent,
        adata.obs["batch"],
        cell_type_labels=adata.obs["cell_type"],
    )
    history = history_to_frame(model.history)
    epochs_completed = epochs_completed_from_history(history, hparams["max_epochs"])
    early_stopped = bool(epochs_completed < hparams["max_epochs"])
    train_n = _n_indices(getattr(model, "train_indices_", None))
    val_n = _n_indices(getattr(model, "validation_indices_", None))
    test_n = _n_indices(getattr(model, "test_indices_", None))

    embeddings_dir.mkdir(parents=True, exist_ok=True)
    np.save(embeddings_dir / f"scvi_matched_seed{seed}_latent.npy", latent)
    cells = pd.DataFrame(
        {
            "row_index": np.arange(adata.n_obs),
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "dataset": adata.obs["dataset"].astype(str).to_numpy(),
            "split": adata.obs["split"].astype(str).to_numpy(),
        }
    )
    cells.to_csv(embeddings_dir / f"scvi_matched_seed{seed}_cells.csv", index=False)

    embed_adata = ad.AnnData(X=latent.copy())
    embed_adata.obs = adata.obs.copy()
    embed_adata.obs_names = adata.obs_names
    embed_adata.obsm["X_scVI_matched"] = latent.copy()
    embed_adata.obsm["X_scVI_matched_umap"] = np.asarray(adata.obsm["X_scVI_matched_umap"])
    embed_adata.uns["scvi_phase"] = "5A"
    embed_adata.uns["x_is_raw_counts"] = False
    embed_adata.uns["X_contents"] = "scVI_matched latent representation (not RNA counts)"
    embed_adata.write_h5ad(embeddings_dir / f"scvi_matched_seed{seed}.h5ad", compression="gzip")

    tables_dir.mkdir(parents=True, exist_ok=True)
    history.to_csv(tables_dir / f"scvi_matched_seed{seed}_training_history.csv", index=False)
    plot_training_history(
        history,
        figures_dir / f"scvi_matched_seed{seed}_training_curve",
        title="scVI_matched training history (scvi-tools logged quantities)",
    )
    plot_umap_by_batch(
        np.asarray(adata.obsm["X_scVI_matched_umap"]),
        adata.obs["batch"],
        figures_dir / f"scvi_matched_seed{seed}_umap_batch",
        title="scVI_matched latent UMAP colored by batch",
    )

    final_elbo_train = _final_metric(history, "elbo_train")
    final_elbo_val = _final_metric(history, "elbo_validation")
    _upsert_baseline_row(
        tables_dir / "baseline_results.csv",
        {
            "model": "scVI_matched",
            "modalities": "RNA",
            "dataset": "PBMC10k+PBMC5k_inner",
            "seed": seed,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_proteins_used": 0,
            "latent_dim": hparams["n_latent"],
            "batch_silhouette": metrics["batch_silhouette"],
            "runtime_seconds": runtime,
            "epochs": epochs_completed,
            "device": actual_device,
            "ARI": metrics["ARI"],
            "NMI": metrics["NMI"],
            "celltype_silhouette": metrics["celltype_silhouette"],
            "celltype_metrics_reason": metrics["celltype_metrics_reason"],
            "elbo_train_final": final_elbo_train,
            "elbo_validation_final": final_elbo_val,
            "peak_rss_mb": peak_rss_mb(),
        },
    )

    source_reloaded = ad.read_h5ad(source_path)
    if "X_scVI" in source_reloaded.obsm or "X_scVI_matched" in source_reloaded.obsm:
        raise ScviDataError("Original processed h5ad unexpectedly gained embeddings.")

    manifest = {
        "dataset": "PBMC10k+PBMC5k_inner",
        "source_h5ad": str(source_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "seed": seed,
        "n_latent": hparams["n_latent"],
        "n_hidden": hparams["n_hidden"],
        "n_layers": hparams["n_layers"],
        "dropout_rate": hparams["dropout_rate"],
        "gene_likelihood": hparams["gene_likelihood"],
        "dispersion": hparams["dispersion"],
        "latent_distribution": hparams["latent_distribution"],
        "batch_size": hparams["batch_size"],
        "train_size": hparams["train_size"],
        "max_epochs": hparams["max_epochs"],
        "epochs_completed": epochs_completed,
        "early_stopping": hparams["early_stopping"],
        "early_stopped": early_stopped,
        "runtime_seconds": runtime,
        "peak_rss_mb": peak_rss_mb(),
        "peak_cuda_allocated_mb": gpu_mem.get("peak_cuda_allocated_mb"),
        "peak_cuda_reserved_mb": gpu_mem.get("peak_cuda_reserved_mb"),
        "requested_device": requested_device,
        "requested_device_detail": requested_detail,
        "actual_device": actual_device,
        "device_fallback_error": fallback_error,
        "gpu_name": (env.get("gpu_names") or [None])[0],
        "cuda_version": env.get("cuda_version"),
        "driver_version": env.get("driver_version"),
        "scvi_tools_version": env["scvi_version"],
        "torch_version": env["torch_version"],
        "python_version": env["python"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "model_class": "scvi.model.SCVI",
        "model_variant": "scvi_matched",
        "seed_backends": seed_info,
        "counts_layer": "counts",
        "batch_key": "batch",
        "counts_layout": "sparse CSR",
        "gene_selection": hparams["gene_selection"],
        "hyperparameters": hparams,
        "unavoidable_differences_vs_totalvi": {
            "scvi_n_layers_applies_to": "both encoder and decoder",
            "totalvi_encoder_layers": 2,
            "totalvi_decoder_layers": 1,
            "scvi_has_protein_decoder": False,
            "totalvi_models_protein": True,
            "scvi_training_lr_default": "SCVI TrainingPlan default (typically 1e-3); not forced to totalVI 4e-3",
            "scvi_reduce_lr_on_plateau": False,
            "totalvi_reduce_lr_on_plateau": True,
        },
        "internal_unsupervised_split": {
            "description": (
                "scvi-tools DataSplitter used for ELBO/early-stopping only. "
                "Independent of obs['split'] from PHASE 2."
            ),
            "train_size_argument": hparams["train_size"],
            "n_train_indices": train_n,
            "n_validation_indices": val_n,
            "n_test_indices": test_n,
            "phase2_obs_split_counts": original_split.value_counts().to_dict(),
        },
        "final_logged_metrics": {
            "elbo_train": final_elbo_train,
            "elbo_validation": final_elbo_val,
            "history_columns": list(history.columns),
        },
        "batch_silhouette": metrics["batch_silhouette"],
        "latent_diagnostics": {
            "dim_mean": metrics["dim_mean"],
            "dim_sd": metrics["dim_sd"],
            "global_mean": metrics["global_mean"],
            "global_sd": metrics["global_sd"],
            "n_nan": metrics["n_nan"],
            "n_inf": metrics["n_inf"],
        },
        "celltype_metrics": "NA",
        "celltype_metrics_reason": metrics["celltype_metrics_reason"],
        "csr_conversion": csr_report,
        "warnings": train_warnings,
        "original_processed_h5ad_modified": False,
        "scvi_default_overwritten": False,
        "totalvi_overwritten": False,
    }
    save_json(manifest, logs_dir / f"scvi_matched_seed{seed}_manifest.json")
    logger.info(
        "PHASE 5A scVI_matched complete: device=%s epochs=%s runtime=%.1fs batch_silhouette=%.4f",
        actual_device,
        epochs_completed,
        runtime,
        metrics["batch_silhouette"],
    )
    return manifest


def _upsert_baseline_row(path, row: dict[str, Any]) -> None:
    new_row = pd.DataFrame([row])
    if path.exists():
        existing = pd.read_csv(path)
        keep = ~(
            (existing["model"].astype(str) == str(row["model"]))
            & (existing["seed"].astype(int) == int(row["seed"]))
        )
        out = pd.concat([existing.loc[keep], new_row], ignore_index=True)
    else:
        out = new_row
    out.to_csv(path, index=False)
