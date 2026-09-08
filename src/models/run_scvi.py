"""PHASE 3: train RNA-only scVI on the combined inner-join CITE-seq object."""

from __future__ import annotations

import traceback
import warnings
from copy import deepcopy
from datetime import datetime, timezone
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse

from src.evaluation.plots import plot_training_history, plot_umap_by_batch
from src.evaluation.representation_metrics import compute_representation_metrics
from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.utils.device import detect_compute_environment, lightning_devices, requested_accelerator
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

EXPECTED_N_OBS = 10849
EXPECTED_N_VARS = 15792
EXPECTED_BATCHES = {"PBMC10k", "PBMC5k"}
CELLTYPE_PLACEHOLDER = "unknown"

DEFAULT_HYPERPARAMETERS = {
    "n_latent": 20,
    "n_hidden": 128,
    "n_layers": 1,
    "dropout_rate": 0.1,
    "dispersion": "gene",
    "gene_likelihood": "zinb",
    "use_observed_lib_size": True,
    "latent_distribution": "normal",
    "train_size": 0.9,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_mode": "min",
}


class ScviDataError(ValueError):
    """Raised when the processed AnnData is unsafe to train on."""


def validate_rna_adata(adata: ad.AnnData) -> None:
    if adata.n_obs != EXPECTED_N_OBS or adata.n_vars != EXPECTED_N_VARS:
        raise ScviDataError(
            f"Unexpected dimensions {adata.n_obs} x {adata.n_vars}; "
            f"expected {EXPECTED_N_OBS} x {EXPECTED_N_VARS}."
        )
    if "counts" not in adata.layers:
        raise ScviDataError("layers['counts'] is missing.")
    if "batch" not in adata.obs:
        raise ScviDataError("obs['batch'] is missing.")
    if "split" not in adata.obs:
        raise ScviDataError("obs['split'] is missing.")
    batches = set(adata.obs["batch"].astype(str).unique())
    if batches != EXPECTED_BATCHES:
        raise ScviDataError(f"obs['batch'] categories are {batches}, expected {EXPECTED_BATCHES}.")
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        data = counts.data
        dense_for_nan = None
    else:
        data = np.asarray(counts)
        dense_for_nan = data
    if dense_for_nan is not None:
        if np.isnan(dense_for_nan).any():
            raise ScviDataError("layers['counts'] contains NaN.")
        if np.any(dense_for_nan < 0):
            raise ScviDataError("layers['counts'] contains negative values.")
        if not np.all(np.abs(dense_for_nan - np.round(dense_for_nan)) < 1e-6):
            raise ScviDataError("layers['counts'] is not integer-valued raw counts.")
    else:
        if np.isnan(data).any():
            raise ScviDataError("layers['counts'] contains NaN.")
        if np.any(data < 0):
            raise ScviDataError("layers['counts'] contains negative values.")
        if not np.all(np.abs(data - np.round(data)) < 1e-6):
            raise ScviDataError("layers['counts'] is not integer-valued raw counts.")


def counts_to_csr(adata: ad.AnnData) -> tuple[ad.AnnData, dict[str, Any]]:
    """Convert counts (and matching X) to CSR without changing values.

    Operates in memory. The caller must not write this object back onto the
    original processed inner-join h5ad path.
    """
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        dense = counts.toarray()
        already_sparse = True
    else:
        dense = np.asarray(counts)
        already_sparse = False
    csr = sparse.csr_matrix(dense)
    identity = bool(np.array_equal(csr.toarray(), dense))
    if not identity:
        raise ScviDataError("CSR conversion changed count values.")
    if sparse.issparse(adata.X):
        x_dense = adata.X.toarray()
    else:
        x_dense = np.asarray(adata.X)
    x_matches_counts = bool(np.array_equal(x_dense, dense))
    del x_dense, dense
    adata.layers["counts"] = csr
    if x_matches_counts:
        adata.X = csr
    report = {
        "counts_layout_before": "csr/sparse" if already_sparse else "dense",
        "counts_layout_after": "csr",
        "csr_nnz": int(csr.nnz),
        "csr_identity_verified": identity,
        "x_also_csr": x_matches_counts,
        "x_matches_counts": x_matches_counts,
    }
    return adata, report


def _n_indices(value: Any) -> int:
    if value is None:
        return 0
    return int(len(value))


def history_to_frame(history: dict[str, pd.DataFrame] | None) -> pd.DataFrame:
    """Align scvi-tools history on epoch when metrics share a common length.

    totalVI logs some quantities on optimizer steps. Concatenating on the raw
    index therefore mixes step numbers with epoch numbers. Metrics that share
    the ELBO/epoch length are reset to a positional epoch column.
    """
    if not history:
        return pd.DataFrame()
    frames: dict[str, pd.DataFrame] = {}
    for key, frame in history.items():
        part = frame.copy()
        if part.shape[1] == 1 and str(part.columns[0]) != key:
            part.columns = [key]
        frames[key] = part
    preferred = next(
        (name for name in ("elbo_train", "elbo_validation", "train_loss_epoch") if name in frames),
        next(iter(frames)),
    )
    epoch_len = len(frames[preferred])
    parts = []
    for key, part in frames.items():
        if len(part) != epoch_len:
            continue
        aligned = part.reset_index(drop=True)
        aligned.index.name = "epoch"
        parts.append(aligned)
    out = pd.concat(parts, axis=1)
    out = out.loc[:, ~out.columns.duplicated()]
    return out.reset_index()


def history_step_frame(history: dict[str, pd.DataFrame] | None) -> pd.DataFrame:
    """Return metrics whose length differs from the epoch-level ELBO series."""
    if not history:
        return pd.DataFrame()
    frames = {key: frame.copy() for key, frame in history.items()}
    preferred = next(
        (name for name in ("elbo_train", "elbo_validation", "train_loss_epoch") if name in frames),
        None,
    )
    if preferred is None:
        return pd.DataFrame()
    epoch_len = len(frames[preferred])
    parts = []
    for key, part in frames.items():
        if len(part) == epoch_len:
            continue
        aligned = part.copy()
        if aligned.shape[1] == 1 and str(aligned.columns[0]) != key:
            aligned.columns = [key]
        aligned.index.name = "step"
        parts.append(aligned)
    if not parts:
        return pd.DataFrame()
    out = pd.concat(parts, axis=1)
    out = out.loc[:, ~out.columns.duplicated()]
    return out.reset_index()


def epochs_completed_from_history(history: pd.DataFrame, max_epochs: int) -> int:
    if "elbo_train" in history.columns and history["elbo_train"].notna().any():
        return int(history["elbo_train"].notna().sum())
    if "epoch" in history.columns and len(history):
        return int(history["epoch"].max()) + 1
    if len(history):
        return min(int(len(history)), int(max_epochs))
    return 0


def _final_metric(history: pd.DataFrame, column: str) -> float | None:
    if column not in history.columns or history[column].dropna().empty:
        return None
    return float(history[column].dropna().iloc[-1])


def train_scvi_with_device_fallback(
    adata: ad.AnnData,
    hparams: dict[str, Any],
    requested_device: str,
    logger,
    devices: int | str | list[int] = "auto",
) -> tuple[scvi.model.SCVI, str, list[str], str | None]:
    """Train scVI, falling back to CPU only after logging the error."""
    warnings_recorded: list[str] = []
    fallback_error = None
    devices_to_try = [requested_device]
    if requested_device != "cpu":
        devices_to_try.append("cpu")

    last_error: Exception | None = None
    for device in devices_to_try:
        if device != requested_device:
            logger.warning(
                "Retrying scVI on CPU after %s failure. Exact error follows.",
                requested_device,
            )
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
        model = scvi.model.SCVI(
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
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                reset_cuda_peak_stats()
                train_devices = devices if device == "gpu" else "auto"
                model.train(
                    max_epochs=hparams["max_epochs"],
                    accelerator=device,
                    devices=train_devices,
                    train_size=hparams["train_size"],
                    batch_size=hparams["batch_size"],
                    early_stopping=hparams["early_stopping"],
                    early_stopping_monitor=hparams["early_stopping_monitor"],
                    early_stopping_patience=hparams["early_stopping_patience"],
                    early_stopping_mode=hparams["early_stopping_mode"],
                    check_val_every_n_epoch=1,
                )
                warnings_recorded.extend(
                    f"{w.category.__name__}: {w.message}" for w in caught
                )
            logger.info("scVI training finished on requested accelerator=%s", device)
            if device != requested_device:
                logger.warning(
                    "Training used fallback device %s rather than requested %s",
                    device,
                    requested_device,
                )
            return model, device, warnings_recorded, fallback_error
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            tb = traceback.format_exc()
            logger.error("scVI training failed on accelerator=%s\n%s", device, tb)
            fallback_error = f"{type(exc).__name__}: {exc}"
            warnings_recorded.append(f"TRAIN_FAIL[{device}]: {fallback_error}")
            continue
    raise RuntimeError(f"scVI training failed on all devices. Last error: {last_error}") from last_error


def run_scvi(config: dict[str, Any] | None = None, seed: int = 0) -> dict[str, Any]:
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    models_dir = resolve_path(cfg["paths"].get("models_dir", "results/models"))
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    logger = get_logger("run_scvi", logs_dir / f"scvi_seed{seed}.log")

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
    adata.write_h5ad(sparse_path, compression="gzip")
    logger.info("Wrote CSR training cache to %s (original inner h5ad was not replaced).", sparse_path)

    hparams = deepcopy(DEFAULT_HYPERPARAMETERS)
    hparams["n_latent"] = int(cfg["models"]["latent_dim"])
    hparams["n_hidden"] = int(cfg["models"]["n_hidden"])
    hparams["n_layers"] = int(cfg["models"]["n_layers"])
    hparams["max_epochs"] = int(cfg["models"]["max_epochs"])
    hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
    hparams["batch_size"] = int(cfg["models"]["batch_size"])
    hparams["gene_selection"] = (
        "all_15792_shared_genes; config n_top_genes=4000 was not applied in PHASE 3"
    )

    requested_device, requested_detail = requested_accelerator(
        env,
        prefer_cuda=bool(cfg["device"].get("prefer_cuda", True)),
        prefer_mps=bool(cfg["device"].get("prefer_mps", True)),
        allow_cpu_fallback=bool(cfg["device"].get("allow_cpu_fallback", cfg["device"].get("fallback_cpu", True))),
    )
    logger.info("Requested device: %s (%s)", requested_device, requested_detail)
    models_dir.mkdir(parents=True, exist_ok=True)
    model_dir = models_dir / f"scvi_seed{seed}"

    timer = Timer()
    timer.start()
    model, actual_device, train_warnings, fallback_error = train_scvi_with_device_fallback(
        adata, hparams, requested_device, logger, devices=lightning_devices(cfg, requested_device)
    )
    runtime = timer.stop()
    gpu_mem = cuda_memory_stats()
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    logger.info("Saved scVI model to %s immediately after training.", model_dir)

    latent = model.get_latent_representation()
    if latent.shape != (EXPECTED_N_OBS, hparams["n_latent"]):
        raise ScviDataError(f"Latent shape {latent.shape} != ({EXPECTED_N_OBS}, {hparams['n_latent']})")
    if not np.array_equal(adata.obs_names, original_order):
        raise ScviDataError("Cell order changed relative to the original AnnData.")
    if not adata.obs["split"].astype(str).equals(original_split):
        raise ScviDataError("obs['split'] was altered during scVI training.")

    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm["X_scVI_umap"] = np.asarray(adata.obsm["X_umap"])

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
    np.save(embeddings_dir / f"scvi_seed{seed}_latent.npy", latent)
    cells = pd.DataFrame(
        {
            "row_index": np.arange(adata.n_obs),
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "dataset": adata.obs["dataset"].astype(str).to_numpy(),
            "split": adata.obs["split"].astype(str).to_numpy(),
        }
    )
    cells.to_csv(embeddings_dir / f"scvi_seed{seed}_cells.csv", index=False)

    embed_adata = ad.AnnData(X=latent.copy())
    embed_adata.obs = adata.obs.copy()
    embed_adata.obs_names = adata.obs_names
    embed_adata.obsm["X_scVI"] = latent.copy()
    embed_adata.obsm["X_scVI_umap"] = np.asarray(adata.obsm["X_scVI_umap"])
    embed_adata.uns["scvi_phase"] = 3
    embed_adata.uns["x_is_raw_counts"] = False
    embed_adata.uns["X_contents"] = "scVI latent representation (not RNA counts)"
    embed_adata.write_h5ad(embeddings_dir / f"scvi_seed{seed}.h5ad", compression="gzip")

    tables_dir.mkdir(parents=True, exist_ok=True)
    history_path = tables_dir / f"scvi_seed{seed}_training_history.csv"
    history.to_csv(history_path, index=False)
    plot_training_history(
        history,
        figures_dir / f"scvi_seed{seed}_training_curve",
        title="scVI training history (scvi-tools logged quantities)",
    )
    plot_umap_by_batch(
        np.asarray(adata.obsm["X_scVI_umap"]),
        adata.obs["batch"],
        figures_dir / f"scvi_seed{seed}_umap_batch",
        title="scVI latent UMAP colored by batch",
    )

    final_elbo_train = _final_metric(history, "elbo_train")
    final_elbo_val = _final_metric(history, "elbo_validation")
    baseline_row = {
        "model": "scVI",
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
    }
    baseline = pd.DataFrame([baseline_row])
    baseline.to_csv(tables_dir / "baseline_results.csv", index=False)

    manifest = {
        "dataset": "PBMC10k+PBMC5k_inner",
        "source_h5ad": str(source_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "batch_categories": sorted(EXPECTED_BATCHES),
        "scvi_tools_version": env["scvi_version"],
        "torch_version": env["torch_version"],
        "python_version": env["python"],
        "lightning_auto_would_select": env["accelerator"],
        "requested_device": requested_device,
        "requested_device_detail": requested_detail,
        "actual_device": actual_device,
        "device_fallback_error": fallback_error,
        "gpu_name": (env.get("gpu_names") or [None])[0],
        "cuda_version": env.get("cuda_version"),
        "driver_version": env.get("driver_version"),
        "peak_cuda_allocated_mb": gpu_mem.get("peak_cuda_allocated_mb"),
        "peak_cuda_reserved_mb": gpu_mem.get("peak_cuda_reserved_mb"),
        "seed": seed,
        "seed_backends": seed_info,
        "n_latent": hparams["n_latent"],
        "max_epochs": hparams["max_epochs"],
        "epochs_completed": epochs_completed,
        "early_stopping": hparams["early_stopping"],
        "early_stopped": early_stopped,
        "runtime_seconds": runtime,
        "counts_layer": "counts",
        "counts_layout": "sparse CSR",
        "gene_selection": hparams["gene_selection"],
        "model_class": "scvi.model.SCVI",
        "hyperparameters": hparams,
        "internal_unsupervised_split": {
            "description": (
                "scvi-tools DataSplitter used for ELBO/early-stopping only. "
                "This is independent of obs['split'] from PHASE 2, which was preserved "
                "for later protein-prediction work and was not used as a supervised holdout."
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
        "leiden": "not generated; leidenalg/igraph are not installed",
        "csr_conversion": csr_report,
        "warnings": train_warnings,
        "peak_rss_mb": peak_rss_mb(),
        "peak_cuda_allocated_mb": gpu_mem.get("peak_cuda_allocated_mb"),
        "peak_cuda_reserved_mb": gpu_mem.get("peak_cuda_reserved_mb"),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "original_processed_h5ad_modified": False,
    }
    save_json(manifest, logs_dir / f"scvi_seed{seed}_manifest.json")
    logger.info(
        "PHASE 3 complete: device=%s epochs=%s runtime=%.1fs batch_silhouette=%.4f",
        actual_device,
        epochs_completed,
        runtime,
        metrics["batch_silhouette"],
    )
    return manifest
