#!/usr/bin/env python
"""Non-scientific CUDA smoke tests for scVI and totalVI.

Trains tiny models on a COPY of ~400 cells for a few epochs.
Writes only under results/logs/server_smoke/. Never touches scientific tables.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch
from scvi.model import SCVI, TOTALVI

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.evaluation.resource_metrics import cuda_memory_stats, peak_rss_mb, reset_cuda_peak_stats
from src.utils.device import detect_compute_environment
from src.utils.io import resolve_path, save_json
from src.utils.seed import set_global_seed

N_CELLS = 400
MAX_EPOCHS = 3
BATCH_SIZE = 64
N_LATENT = 8
N_HIDDEN = 64


def _subset(adata: ad.AnnData) -> ad.AnnData:
    rng = np.random.default_rng(0)
    idx = rng.choice(adata.n_obs, size=min(N_CELLS, adata.n_obs), replace=False)
    idx = np.sort(idx)
    return adata[idx].copy()


def _require_cuda() -> None:
    if not torch.cuda.is_available():
        raise SystemExit("torch.cuda.is_available() is False. Stop. Fix CUDA/PyTorch.")
    name = torch.cuda.get_device_name(0)
    print(f"CUDA device: {name}")
    print(f"torch: {torch.__version__} cuda: {torch.version.cuda}")
    print(f"scvi-tools: {scvi.__version__}")


def _scvi_smoke(adata: ad.AnnData) -> dict:
    reset_cuda_peak_stats()
    SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    model = SCVI(
        adata,
        n_latent=N_LATENT,
        n_hidden=N_HIDDEN,
        n_layers=1,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=MAX_EPOCHS,
        accelerator="gpu",
        devices=[0],
        batch_size=BATCH_SIZE,
        early_stopping=False,
        check_val_every_n_epoch=1,
    )
    latent = np.asarray(model.get_latent_representation())
    gpu = cuda_memory_stats()
    if not np.isfinite(latent).all():
        raise RuntimeError("scVI smoke latent contains NaN/Inf")
    device = str(next(model.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"scVI parameters are on {device}, expected CUDA")
    return {
        "model": "scVI_smoke",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "latent_shape": list(latent.shape),
        "latent_finite": True,
        "parameter_device": device,
        "epochs": MAX_EPOCHS,
        "gpu_memory": gpu,
        "peak_rss_mb": peak_rss_mb(),
    }


def _totalvi_smoke(adata: ad.AnnData) -> dict:
    reset_cuda_peak_stats()
    protein_key = "protein_expression" if "protein_expression" in adata.obsm else "protein_counts"
    TOTALVI.setup_anndata(
        adata,
        protein_expression_obsm_key=protein_key,
        batch_key="batch",
        layer="counts",
    )
    model = TOTALVI(
        adata,
        n_latent=N_LATENT,
        n_hidden=N_HIDDEN,
        n_layers_encoder=1,
        n_layers_decoder=1,
        gene_likelihood="nb",
        empirical_protein_background_prior=False,
    )
    model.train(
        max_epochs=MAX_EPOCHS,
        accelerator="gpu",
        devices=[0],
        batch_size=BATCH_SIZE,
        early_stopping=False,
        check_val_every_n_epoch=1,
    )
    latent = np.asarray(model.get_latent_representation())
    rna_norm, protein_norm = model.get_normalized_expression(
        return_mean=True, include_protein_background=False, scale_protein=False
    )
    gpu = cuda_memory_stats()
    if not np.isfinite(latent).all():
        raise RuntimeError("totalVI smoke latent contains NaN/Inf")
    protein_arr = np.asarray(protein_norm)
    if not np.isfinite(protein_arr).all():
        raise RuntimeError("totalVI normalized protein contains NaN/Inf")
    device = str(next(model.module.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"totalVI parameters are on {device}, expected CUDA")
    return {
        "model": "totalVI_smoke",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_proteins": int(np.asarray(adata.obsm[protein_key]).shape[1]),
        "latent_shape": list(latent.shape),
        "latent_finite": True,
        "normalized_protein_finite": True,
        "normalized_protein_shape": list(protein_arr.shape),
        "parameter_device": device,
        "epochs": MAX_EPOCHS,
        "gpu_memory": gpu,
        "peak_rss_mb": peak_rss_mb(),
        "rna_norm_ignored": True,
        "rna_norm_shape": list(np.asarray(rna_norm).shape),
    }


def main() -> None:
    set_global_seed(0)
    _require_cuda()
    env = detect_compute_environment()
    source = resolve_path("data/processed/pbmc_cite_combined_inner.h5ad")
    full = ad.read_h5ad(source)
    subset = _subset(full)
    del full
    out_dir = resolve_path("results/logs/server_smoke")
    out_dir.mkdir(parents=True, exist_ok=True)

    scvi_result = _scvi_smoke(subset.copy())
    totalvi_result = _totalvi_smoke(subset.copy())
    report = {
        "scientific_result": False,
        "purpose": "CUDA / scvi-tools / DataLoader smoke test only",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": env.get("hostname"),
        "gpu_name": (env.get("gpu_names") or [None])[0],
        "torch_version": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "scvi_tools": scvi.__version__,
        "source_h5ad": str(source),
        "source_unmodified": True,
        "scvi": scvi_result,
        "totalvi": totalvi_result,
    }
    save_json(report, out_dir / "cuda_smoke.json")
    print(json.dumps(report, indent=2, default=str))
    print(f"Smoke outputs isolated in {out_dir}")


if __name__ == "__main__":
    main()
