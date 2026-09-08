"""PHASE 4: train RNA+protein totalVI on the combined inner-join CITE-seq object."""

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
from scvi.model import TOTALVI

from src.data.load_citeseq import clean_protein_name
from src.evaluation.plots import (
    plot_side_by_side_umap_by_batch,
    plot_training_history,
    plot_umap_by_batch,
)
from src.evaluation.representation_metrics import (
    compute_representation_metrics,
    exploratory_latent_comparison,
)
from src.evaluation.resource_metrics import Timer, cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.models.run_scvi import (
    EXPECTED_BATCHES,
    EXPECTED_N_OBS,
    EXPECTED_N_VARS,
    ScviDataError,
    _final_metric,
    _n_indices,
    counts_to_csr,
    epochs_completed_from_history,
    history_step_frame,
    history_to_frame,
    requested_accelerator,
    validate_rna_adata,
)
from src.utils.device import detect_compute_environment, lightning_devices
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

EXPECTED_N_PROTEINS = 14
PROTEIN_OBSM_KEY = "protein_expression"
CELLTYPE_METRICS_REASON = "Ground-truth cell-type labels unavailable"

# Architecture defaults are TOTALVI/TOTALVAE 1.3.3 defaults, not scVI's.
# Aligned quantities: n_latent, cells, genes, batch, seed, max_epochs, batch_size.
DEFAULT_HYPERPARAMETERS = {
    "n_latent": 20,
    "gene_dispersion": "gene",
    "protein_dispersion": "protein",
    "gene_likelihood": "nb",
    "latent_distribution": "normal",
    "empirical_protein_background_prior": None,
    "override_missing_proteins": False,
    "n_hidden": 256,
    "n_layers_encoder": 2,
    "n_layers_decoder": 1,
    "dropout_rate_encoder": 0.2,
    "dropout_rate_decoder": 0.2,
    "train_size": 0.9,
    "lr": 4e-3,
    "reduce_lr_on_plateau": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_mode": "min",
}


class TotalviDataError(ValueError):
    """Raised when the processed AnnData is unsafe for totalVI."""


def _protein_frame(adata: ad.AnnData) -> pd.DataFrame:
    if PROTEIN_OBSM_KEY not in adata.obsm:
        raise TotalviDataError("obsm['protein_expression'] is missing.")
    protein = adata.obsm[PROTEIN_OBSM_KEY]
    if not isinstance(protein, pd.DataFrame):
        raise TotalviDataError(
            "obsm['protein_expression'] must be a DataFrame so protein names "
            "and column order can be recovered."
        )
    return protein


def validate_totalvi_adata(adata: ad.AnnData) -> pd.DataFrame:
    try:
        validate_rna_adata(adata)
    except ScviDataError as exc:
        raise TotalviDataError(str(exc)) from exc
    protein = _protein_frame(adata)
    if protein.shape != (EXPECTED_N_OBS, EXPECTED_N_PROTEINS):
        raise TotalviDataError(
            f"protein_expression shape {protein.shape} != "
            f"({EXPECTED_N_OBS}, {EXPECTED_N_PROTEINS})."
        )
    if not protein.index.equals(adata.obs_names):
        raise TotalviDataError(
            "protein_expression row order does not match adata.obs_names."
        )
    values = protein.to_numpy()
    if np.isnan(values).any():
        raise TotalviDataError("protein_expression contains NaN.")
    if np.any(values < 0):
        raise TotalviDataError("protein_expression contains negative values.")
    if not np.all(np.abs(values - np.round(values)) < 1e-6):
        raise TotalviDataError("protein_expression is not integer-valued raw counts.")
    if "protein_names" in adata.uns:
        uns_names = [str(x) for x in np.asarray(adata.uns["protein_names"])]
        frame_names = [str(c) for c in protein.columns]
        if uns_names != frame_names:
            raise TotalviDataError(
                "uns['protein_names'] does not match protein_expression columns. "
                f"uns={uns_names} columns={frame_names}"
            )
    if "protein_counts" in adata.obsm:
        counts = adata.obsm["protein_counts"]
        if isinstance(counts, pd.DataFrame) and not np.array_equal(
            counts.to_numpy(), values
        ):
            raise TotalviDataError(
                "obsm['protein_counts'] and obsm['protein_expression'] differ."
            )
    return protein


def protein_name_table(protein: pd.DataFrame) -> pd.DataFrame:
    originals = [str(c) for c in protein.columns]
    return pd.DataFrame(
        {
            "protein_index": np.arange(len(originals), dtype=int),
            "protein_name_original": originals,
            "protein_name_cleaned": [clean_protein_name(name) for name in originals],
        }
    )


def train_totalvi_with_device_fallback(
    adata: ad.AnnData,
    hparams: dict[str, Any],
    requested_device: str,
    logger,
    devices: int | str | list[int] = "auto",
) -> tuple[TOTALVI, str, list[str], str | None]:
    """Train totalVI, falling back to CPU only after logging the error."""
    warnings_recorded: list[str] = []
    fallback_error = None
    devices_to_try = [requested_device]
    if requested_device != "cpu":
        devices_to_try.append("cpu")

    last_error: Exception | None = None
    for device in devices_to_try:
        if device != requested_device:
            logger.warning(
                "Retrying totalVI on CPU after %s failure. Exact error follows.",
                requested_device,
            )
        # scvi-tools 1.3.3: protein_expression_obsm_key is required and first
        # after adata. setup_anndata still works; setup_mudata is recommended
        # but documented as not changing model performance.
        with warnings.catch_warnings(record=True) as setup_caught:
            warnings.simplefilter("always")
            TOTALVI.setup_anndata(
                adata,
                protein_expression_obsm_key=PROTEIN_OBSM_KEY,
                protein_names_uns_key=None,
                batch_key="batch",
                layer="counts",
            )
            warnings_recorded.extend(
                f"{w.category.__name__}: {w.message}" for w in setup_caught
            )
        model = TOTALVI(
            adata,
            n_latent=hparams["n_latent"],
            gene_dispersion=hparams["gene_dispersion"],
            protein_dispersion=hparams["protein_dispersion"],
            gene_likelihood=hparams["gene_likelihood"],
            latent_distribution=hparams["latent_distribution"],
            empirical_protein_background_prior=hparams[
                "empirical_protein_background_prior"
            ],
            override_missing_proteins=hparams["override_missing_proteins"],
            n_hidden=hparams["n_hidden"],
            n_layers_encoder=hparams["n_layers_encoder"],
            n_layers_decoder=hparams["n_layers_decoder"],
            dropout_rate_encoder=hparams["dropout_rate_encoder"],
            dropout_rate_decoder=hparams["dropout_rate_decoder"],
        )
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                reset_cuda_peak_stats()
                train_devices = devices if device == "gpu" else "auto"
                model.train(
                    max_epochs=hparams["max_epochs"],
                    lr=hparams["lr"],
                    accelerator=device,
                    devices=train_devices,
                    train_size=hparams["train_size"],
                    batch_size=hparams["batch_size"],
                    early_stopping=hparams["early_stopping"],
                    early_stopping_monitor=hparams["early_stopping_monitor"],
                    early_stopping_patience=hparams["early_stopping_patience"],
                    early_stopping_mode=hparams["early_stopping_mode"],
                    reduce_lr_on_plateau=hparams["reduce_lr_on_plateau"],
                    check_val_every_n_epoch=1,
                )
                warnings_recorded.extend(
                    f"{w.category.__name__}: {w.message}" for w in caught
                )
            logger.info("totalVI training finished on requested accelerator=%s", device)
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
            logger.error("totalVI training failed on accelerator=%s\n%s", device, tb)
            fallback_error = f"{type(exc).__name__}: {exc}"
            warnings_recorded.append(f"TRAIN_FAIL[{device}]: {fallback_error}")
            continue
    raise RuntimeError(
        f"totalVI training failed on all devices. Last error: {last_error}"
    ) from last_error


def _save_protein_matrix(
    frame: pd.DataFrame,
    expected_index: pd.Index,
    expected_columns: list[str],
    path,
    label: str,
) -> None:
    if not frame.index.equals(expected_index):
        raise TotalviDataError(f"{label} cell order does not match source AnnData.")
    got = [str(c) for c in frame.columns]
    if got != expected_columns:
        raise TotalviDataError(
            f"{label} protein columns {got} != expected order {expected_columns}."
        )
    if not np.isfinite(frame.to_numpy()).all():
        raise TotalviDataError(f"{label} contains NaN or Inf.")
    frame.to_csv(path, index_label="cell_id")


def _upsert_baseline_row(path, row: dict[str, Any]) -> pd.DataFrame:
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
    return out


def run_totalvi(config: dict[str, Any] | None = None, seed: int = 0) -> dict[str, Any]:
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    models_dir = resolve_path(cfg["paths"].get("models_dir", "results/models"))
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    logger = get_logger("run_totalvi", logs_dir / f"totalvi_seed{seed}.log")

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
    protein = validate_totalvi_adata(adata_original)
    original_order = adata_original.obs_names.copy()
    original_split = adata_original.obs["split"].astype(str).copy()
    original_counts_checksum = float(np.asarray(adata_original.layers["counts"]).sum())
    original_protein_checksum = float(protein.to_numpy().sum())
    if adata_original.obs["cell_type"].astype(str).nunique() == 1:
        logger.info("obs['cell_type'] is a placeholder and will not be used.")

    protein_table = protein_name_table(protein)
    logs_dir.mkdir(parents=True, exist_ok=True)
    protein_table.to_csv(logs_dir / f"totalvi_seed{seed}_proteins.csv", index=False)
    protein_names = protein_table["protein_name_original"].tolist()
    logger.info("totalVI proteins in column order: %s", protein_names)
    if adata_original.uns.get("protein_names_clean") is None:
        logger.warning(
            "PHASE 2 combined object has uns['protein_names_clean'] = None. "
            "Cleaned names are computed for the protein table only; the original "
            "h5ad is not modified."
        )

    adata, csr_report = counts_to_csr(adata_original)
    # Keep protein DataFrame columns/order unchanged after RNA CSR conversion.
    if not isinstance(adata.obsm[PROTEIN_OBSM_KEY], pd.DataFrame):
        raise TotalviDataError("protein_expression lost DataFrame type after CSR conversion.")
    if [str(c) for c in adata.obsm[PROTEIN_OBSM_KEY].columns] != protein_names:
        raise TotalviDataError("protein column order changed during CSR conversion.")

    sparse_path = processed_dir / "pbmc_cite_combined_inner_sparse.h5ad"
    if not sparse_path.exists():
        adata.write_h5ad(sparse_path, compression="gzip")
        logger.info("Wrote CSR training cache to %s", sparse_path)
    else:
        logger.info("Reusing existing CSR cache %s without rewriting it.", sparse_path)

    hparams = deepcopy(DEFAULT_HYPERPARAMETERS)
    hparams["n_latent"] = int(cfg["models"]["latent_dim"])
    hparams["max_epochs"] = int(cfg["models"]["max_epochs"])
    hparams["early_stopping"] = bool(cfg["models"]["early_stopping"])
    hparams["batch_size"] = int(cfg["models"]["batch_size"])
    hparams["gene_selection"] = (
        "all_15792_shared_genes; config n_top_genes=4000 was not applied in PHASE 4"
    )
    hparams["architecture_note"] = (
        "n_hidden/n_layers/dropout/gene_likelihood use TOTALVI 1.3.3 defaults "
        "(256 / encoder=2,decoder=1 / 0.2 / nb), not the scVI config values "
        "(128 / 1 / 0.1 / zinb)."
    )

    requested_device, requested_detail = requested_accelerator(
        env,
        prefer_cuda=bool(cfg["device"].get("prefer_cuda", True)),
        prefer_mps=bool(cfg["device"].get("prefer_mps", True)),
        allow_cpu_fallback=bool(cfg["device"].get("allow_cpu_fallback", cfg["device"].get("fallback_cpu", True))),
    )
    logger.info("Requested device: %s (%s)", requested_device, requested_detail)
    models_dir.mkdir(parents=True, exist_ok=True)
    model_dir = models_dir / f"totalvi_seed{seed}"

    timer = Timer()
    timer.start()
    model, actual_device, train_warnings, fallback_error = train_totalvi_with_device_fallback(
        adata, hparams, requested_device, logger, devices=lightning_devices(cfg, requested_device)
    )
    runtime = timer.stop()
    gpu_mem = cuda_memory_stats()
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    logger.info("Saved totalVI model to %s immediately after training.", model_dir)

    latent = np.asarray(model.get_latent_representation())
    if latent.shape != (EXPECTED_N_OBS, hparams["n_latent"]):
        raise TotalviDataError(
            f"Latent shape {latent.shape} != ({EXPECTED_N_OBS}, {hparams['n_latent']})"
        )
    if not np.isfinite(latent).all():
        raise TotalviDataError("totalVI latent representation contains NaN or Inf.")
    if not np.array_equal(adata.obs_names, original_order):
        raise TotalviDataError("Cell order changed relative to the original AnnData.")
    if not adata.obs["split"].astype(str).equals(original_split):
        raise TotalviDataError("obs['split'] was altered during totalVI training.")
    if "X_scVI" in adata.obsm:
        raise TotalviDataError("Unexpected X_scVI on the training object; refuse to overwrite.")

    current_counts_sum = float(
        adata.layers["counts"].toarray().sum()
        if hasattr(adata.layers["counts"], "toarray")
        else np.asarray(adata.layers["counts"]).sum()
    )
    current_protein_sum = float(adata.obsm[PROTEIN_OBSM_KEY].to_numpy().sum())
    if abs(current_counts_sum - original_counts_checksum) > 1e-3:
        raise TotalviDataError("Raw RNA counts changed during totalVI training.")
    if abs(current_protein_sum - original_protein_checksum) > 1e-3:
        raise TotalviDataError("Raw protein counts changed during totalVI training.")

    adata.obsm["X_totalVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_totalVI", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.obsm["X_totalVI_umap"] = np.asarray(adata.obsm["X_umap"])

    metrics = compute_representation_metrics(
        latent,
        adata.obs["batch"],
        cell_type_labels=adata.obs["cell_type"],
    )
    history = history_to_frame(model.history)
    step_history = history_step_frame(model.history)
    epochs_completed = epochs_completed_from_history(history, hparams["max_epochs"])
    early_stopped = bool(epochs_completed < hparams["max_epochs"])

    train_n = _n_indices(getattr(model, "train_indices_", None))
    val_n = _n_indices(getattr(model, "validation_indices_", None))
    test_n = _n_indices(getattr(model, "test_indices_", None))

    embeddings_dir.mkdir(parents=True, exist_ok=True)
    np.save(embeddings_dir / f"totalvi_seed{seed}_latent.npy", latent)
    cells = pd.DataFrame(
        {
            "row_index": np.arange(adata.n_obs),
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "dataset": adata.obs["dataset"].astype(str).to_numpy(),
            "split": adata.obs["split"].astype(str).to_numpy(),
        }
    )
    cells.to_csv(embeddings_dir / f"totalvi_seed{seed}_cells.csv", index=False)

    protein_outputs: dict[str, str] = {}
    rna_norm, protein_norm = model.get_normalized_expression(
        return_mean=True,
        include_protein_background=False,
        scale_protein=False,
    )
    _save_protein_matrix(
        protein_norm,
        adata.obs_names,
        protein_names,
        embeddings_dir / f"totalvi_seed{seed}_normalized_protein.csv",
        "get_normalized_expression(include_protein_background=False)",
    )
    protein_outputs["normalized_protein"] = (
        "TOTALVI.get_normalized_expression protein return; "
        "include_protein_background=False (API default foreground component)"
    )
    del rna_norm

    _, protein_with_bg = model.get_normalized_expression(
        return_mean=True,
        include_protein_background=True,
        scale_protein=False,
    )
    _save_protein_matrix(
        protein_with_bg,
        adata.obs_names,
        protein_names,
        embeddings_dir / f"totalvi_seed{seed}_normalized_protein_with_background.csv",
        "get_normalized_expression(include_protein_background=True)",
    )
    protein_outputs["normalized_protein_with_background"] = (
        "TOTALVI.get_normalized_expression protein return; include_protein_background=True"
    )

    foreground = model.get_protein_foreground_probability(return_mean=True)
    _save_protein_matrix(
        foreground,
        adata.obs_names,
        protein_names,
        embeddings_dir / f"totalvi_seed{seed}_protein_foreground_probability.csv",
        "get_protein_foreground_probability",
    )
    protein_outputs["protein_foreground_probability"] = (
        "TOTALVI.get_protein_foreground_probability; (1 - pi_nt) in the totalVI paper"
    )

    background_mean_status = "saved"
    try:
        background = model.get_protein_background_mean(
            None, np.arange(adata.n_obs), hparams["batch_size"]
        )
        background = np.asarray(background)
        if background.shape != (EXPECTED_N_OBS, EXPECTED_N_PROTEINS):
            raise TotalviDataError(
                f"get_protein_background_mean shape {background.shape} unexpected."
            )
        if not np.isfinite(background).all():
            raise TotalviDataError("get_protein_background_mean contains NaN or Inf.")
        bg_frame = pd.DataFrame(
            background, index=adata.obs_names, columns=protein_names
        )
        _save_protein_matrix(
            bg_frame,
            adata.obs_names,
            protein_names,
            embeddings_dir / f"totalvi_seed{seed}_protein_background_mean.csv",
            "get_protein_background_mean",
        )
        protein_outputs["protein_background_mean"] = (
            "TOTALVI.get_protein_background_mean (py_['rate_back'])"
        )
    except Exception as exc:  # noqa: BLE001
        background_mean_status = f"{type(exc).__name__}: {exc}"
        logger.warning("get_protein_background_mean was not saved: %s", background_mean_status)
        protein_outputs["protein_background_mean"] = (
            f"API exists but was not saved: {background_mean_status}"
        )

    embed_adata = ad.AnnData(X=latent.copy())
    embed_adata.obs = adata.obs.copy()
    embed_adata.obs_names = adata.obs_names
    embed_adata.obsm["X_totalVI"] = latent.copy()
    embed_adata.obsm["X_totalVI_umap"] = np.asarray(adata.obsm["X_totalVI_umap"])
    embed_adata.uns["totalvi_phase"] = 4
    embed_adata.uns["x_is_raw_counts"] = False
    embed_adata.uns["X_contents"] = "totalVI joint latent representation (not RNA counts)"
    embed_adata.uns["protein_names"] = np.array(protein_names, dtype=object)
    embed_adata.write_h5ad(embeddings_dir / f"totalvi_seed{seed}.h5ad", compression="gzip")

    scvi_latent_path = embeddings_dir / f"scvi_seed{seed}_latent.npy"
    scvi_embed_path = embeddings_dir / f"scvi_seed{seed}.h5ad"
    scvi_cells_path = embeddings_dir / f"scvi_seed{seed}_cells.csv"
    comparison_status = "created"
    exploratory = None
    if scvi_latent_path.exists() and scvi_embed_path.exists() and scvi_cells_path.exists():
        scvi_latent = np.load(scvi_latent_path)
        scvi_cells = pd.read_csv(scvi_cells_path)
        if list(scvi_cells["cell_id"].astype(str)) != list(adata.obs_names.astype(str)):
            raise TotalviDataError("scVI cell order does not match totalVI/source AnnData.")
        if scvi_latent.shape != latent.shape:
            raise TotalviDataError(
                f"scVI latent shape {scvi_latent.shape} != totalVI {latent.shape}."
            )
        scvi_embed = ad.read_h5ad(scvi_embed_path)
        if "X_scVI" not in scvi_embed.obsm:
            raise TotalviDataError("scVI embedding AnnData is missing obsm['X_scVI'].")
        if not np.array_equal(scvi_embed.obs_names, adata.obs_names):
            raise TotalviDataError("scVI embedding cell order does not match totalVI.")
        comparison = ad.AnnData(X=latent.copy())
        comparison.obs = adata.obs.copy()
        comparison.obs_names = adata.obs_names
        comparison.obsm["X_scVI"] = np.asarray(scvi_embed.obsm["X_scVI"])
        comparison.obsm["X_totalVI"] = latent.copy()
        comparison.obsm["X_scVI_umap"] = np.asarray(scvi_embed.obsm["X_scVI_umap"])
        comparison.obsm["X_totalVI_umap"] = np.asarray(adata.obsm["X_totalVI_umap"])
        comparison.uns["comparison_phase"] = 4
        comparison.uns["X_contents"] = "totalVI latent; paired embeddings are in obsm"
        comparison.uns["note"] = (
            "Exploratory paired embeddings only. Identical cells and row order. "
            "Not a biological ranking."
        )
        comparison.write_h5ad(
            embeddings_dir / f"scvi_totalvi_seed{seed}_comparison.h5ad",
            compression="gzip",
        )
        exploratory = exploratory_latent_comparison(scvi_latent, latent, seed=seed)
        exploratory_path = tables_dir / f"scvi_vs_totalvi_seed{seed}_exploratory_latent_comparison.csv"
        tables_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([exploratory]).to_csv(exploratory_path, index=False)
        plot_side_by_side_umap_by_batch(
            np.asarray(scvi_embed.obsm["X_scVI_umap"]),
            np.asarray(adata.obsm["X_totalVI_umap"]),
            adata.obs["batch"],
            figures_dir / f"scvi_vs_totalvi_seed{seed}_batch",
            title_left="scVI latent UMAP (RNA)",
            title_right="totalVI latent UMAP (RNA+protein)",
            suptitle=(
                "Exploratory batch-colored UMAPs. Not a biological ranking; "
                "validated cell-type labels are unavailable."
            ),
        )
    else:
        comparison_status = "skipped; PHASE 3 scVI embeddings were not found"
        logger.warning(comparison_status)

    tables_dir.mkdir(parents=True, exist_ok=True)
    history_path = tables_dir / f"totalvi_seed{seed}_training_history.csv"
    history.to_csv(history_path, index=False)
    if len(step_history):
        step_history.to_csv(
            tables_dir / f"totalvi_seed{seed}_training_history_step.csv", index=False
        )
    plot_training_history(
        history,
        figures_dir / f"totalvi_seed{seed}_training_curve",
        title="totalVI training history (scvi-tools logged quantities)",
    )
    plot_umap_by_batch(
        np.asarray(adata.obsm["X_totalVI_umap"]),
        adata.obs["batch"],
        figures_dir / f"totalvi_seed{seed}_umap_batch",
        title="totalVI latent UMAP colored by batch",
    )

    final_metrics = {
        column: _final_metric(history, column)
        for column in history.columns
        if column not in {"epoch", "step"}
    }
    baseline_row = {
        "model": "totalVI",
        "modalities": "RNA+Protein",
        "dataset": "PBMC10k+PBMC5k_inner",
        "seed": seed,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_proteins_used": EXPECTED_N_PROTEINS,
        "latent_dim": hparams["n_latent"],
        "batch_silhouette": metrics["batch_silhouette"],
        "runtime_seconds": runtime,
        "epochs": epochs_completed,
        "device": actual_device,
        "ARI": metrics["ARI"],
        "NMI": metrics["NMI"],
        "celltype_silhouette": metrics["celltype_silhouette"],
        "celltype_metrics_reason": metrics["celltype_metrics_reason"],
        "elbo_train_final": final_metrics.get("elbo_train"),
        "elbo_validation_final": final_metrics.get("elbo_validation"),
        "peak_rss_mb": peak_rss_mb(),
    }
    baseline = _upsert_baseline_row(tables_dir / "baseline_results.csv", baseline_row)

    comparison_rows = []
    scvi_rows = baseline[baseline["model"].astype(str) == "scVI"]
    if len(scvi_rows):
        scvi_row = scvi_rows.iloc[-1]
        comparison_rows.append(
            {
                "model": "scVI",
                "modalities": "RNA",
                "n_cells": int(scvi_row["n_cells"]),
                "n_genes": int(scvi_row["n_genes"]),
                "n_proteins": int(scvi_row.get("n_proteins_used", 0)),
                "latent_dim": int(scvi_row["latent_dim"]),
                "epochs": int(scvi_row["epochs"]),
                "runtime_seconds": float(scvi_row["runtime_seconds"]),
                "peak_rss_mb": float(scvi_row["peak_rss_mb"]),
                "device": str(scvi_row["device"]),
                "batch_silhouette": float(scvi_row["batch_silhouette"]),
                "ARI": np.nan,
                "NMI": np.nan,
                "celltype_ASW": np.nan,
                "celltype_metrics_reason": CELLTYPE_METRICS_REASON,
            }
        )
    comparison_rows.append(
        {
            "model": "totalVI",
            "modalities": "RNA+Protein",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_proteins": EXPECTED_N_PROTEINS,
            "latent_dim": hparams["n_latent"],
            "epochs": epochs_completed,
            "runtime_seconds": runtime,
            "peak_rss_mb": peak_rss_mb(),
            "device": actual_device,
            "batch_silhouette": metrics["batch_silhouette"],
            "ARI": np.nan,
            "NMI": np.nan,
            "celltype_ASW": np.nan,
            "celltype_metrics_reason": CELLTYPE_METRICS_REASON,
        }
    )
    pd.DataFrame(comparison_rows).to_csv(
        tables_dir / f"scvi_vs_totalvi_seed{seed}.csv", index=False
    )

    source_reloaded = ad.read_h5ad(source_path)
    source_counts_sum = float(np.asarray(source_reloaded.layers["counts"]).sum())
    source_protein_sum = float(source_reloaded.obsm[PROTEIN_OBSM_KEY].to_numpy().sum())
    if abs(source_counts_sum - original_counts_checksum) > 1e-3:
        raise TotalviDataError("Original processed h5ad RNA counts were modified.")
    if abs(source_protein_sum - original_protein_checksum) > 1e-3:
        raise TotalviDataError("Original processed h5ad protein counts were modified.")
    if "X_scVI" in source_reloaded.obsm or "X_totalVI" in source_reloaded.obsm:
        raise TotalviDataError("Original processed h5ad unexpectedly gained embeddings.")

    empirical_prior_used = bool(model.summary_stats.n_proteins > 10) if hparams[
        "empirical_protein_background_prior"
    ] is None else bool(hparams["empirical_protein_background_prior"])

    manifest = {
        "dataset": "PBMC10k+PBMC5k_inner",
        "source_h5ad": str(source_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_proteins": EXPECTED_N_PROTEINS,
        "protein_names": protein_names,
        "protein_names_cleaned": protein_table["protein_name_cleaned"].tolist(),
        "protein_obsm_key": PROTEIN_OBSM_KEY,
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
        "seed": seed,
        "seed_backends": seed_info,
        "n_latent": hparams["n_latent"],
        "batch_size": hparams["batch_size"],
        "max_epochs": hparams["max_epochs"],
        "epochs_completed": epochs_completed,
        "early_stopping": hparams["early_stopping"],
        "early_stopped": early_stopped,
        "runtime_seconds": runtime,
        "peak_rss_mb": peak_rss_mb(),
        "peak_cuda_allocated_mb": gpu_mem.get("peak_cuda_allocated_mb"),
        "peak_cuda_reserved_mb": gpu_mem.get("peak_cuda_reserved_mb"),
        "counts_layer": "counts",
        "counts_layout": "sparse CSR",
        "gene_selection": hparams["gene_selection"],
        "model_class": "scvi.model.TOTALVI",
        "setup_api": (
            "TOTALVI.setup_anndata(adata, protein_expression_obsm_key="
            "'protein_expression', batch_key='batch', layer='counts'). "
            "scvi-tools 1.3.3 emits a DeprecationWarning recommending setup_mudata; "
            "the warning states this does not influence model performance."
        ),
        "effective_model_hyperparameters": hparams,
        "empirical_protein_background_prior_resolved": empirical_prior_used,
        "use_adversarial_classifier": bool(model._use_adversarial_classifier),
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
            **{k: v for k, v in final_metrics.items() if v is not None},
            "history_columns": list(history.columns),
        },
        "batch_silhouette": metrics["batch_silhouette"],
        "batch_silhouette_interpretation": (
            "Descriptive only. PBMC10k and PBMC5k are separate datasets and may "
            "differ in composition. A lower batch silhouette is not automatically "
            "a better biological representation."
        ),
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
        "protein_outputs": protein_outputs,
        "normalized_rna_saved": False,
        "normalized_rna_reason": (
            "TOTALVI.get_normalized_expression also returns an RNA matrix of shape "
            f"{EXPECTED_N_OBS} x {EXPECTED_N_VARS}. It was computed then discarded; "
            "not written because a dense float matrix is not storage-efficient. "
            "Raw layers['counts'] were left unchanged."
        ),
        "protein_background_mean_status": background_mean_status,
        "comparison_anndata": comparison_status,
        "exploratory_latent_comparison": exploratory,
        "csr_conversion": csr_report,
        "warnings": train_warnings,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "original_processed_h5ad_modified": False,
        "scvi_baseline_retrained": False,
    }
    save_json(manifest, logs_dir / f"totalvi_seed{seed}_manifest.json")
    logger.info(
        "PHASE 4 complete: device=%s epochs=%s runtime=%.1fs batch_silhouette=%.4f",
        actual_device,
        epochs_completed,
        runtime,
        metrics["batch_silhouette"],
    )
    return manifest
