"""Inspect scvi-tools 1.3.3 SCVI/TOTALVI query (scArches) support.

Not a scientific result. Uses a small copy of real cells.
"""

from __future__ import annotations

import inspect
from datetime import datetime, timezone
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch
from scipy import sparse
from scvi.model import SCVI, TOTALVI

from src.utils.io import resolve_path, save_json


def _csr_copy(adata: ad.AnnData) -> ad.AnnData:
    out = adata.copy()
    counts = out.layers["counts"]
    if sparse.issparse(counts):
        out.layers["counts"] = counts.tocsr()
    else:
        out.layers["counts"] = sparse.csr_matrix(np.asarray(counts))
    out.X = out.layers["counts"]
    return out


def _inspect_class(cls) -> dict[str, Any]:
    names = {m for m in dir(cls) if not m.startswith("_")}
    return {
        "class": cls.__name__,
        "mro": [c.__name__ for c in cls.__mro__],
        "has_prepare_query_anndata": "prepare_query_anndata" in names,
        "has_load_query_data": "load_query_data" in names,
        "has_prepare_query_mudata": "prepare_query_mudata" in names,
        "load_query_data_signature": str(inspect.signature(cls.load_query_data)),
        "prepare_query_anndata_signature": str(inspect.signature(cls.prepare_query_anndata)),
        "arches_mixin": "ArchesMixin" in [c.__name__ for c in cls.__mro__],
    }


def _trainable_param_names(model) -> list[str]:
    return [name for name, par in model.module.named_parameters() if par.requires_grad]


def run_api_validation() -> dict[str, Any]:
    if not torch.cuda.is_available():
        raise RuntimeError("PHASE 10 API validation requires CUDA.")
    source_path = resolve_path("data/processed/pbmc_cite_combined_inner_annotated.h5ad")
    full = ad.read_h5ad(source_path)
    rng = np.random.default_rng(0)
    b10 = np.where(full.obs["batch"].astype(str).to_numpy() == "PBMC10k")[0]
    b5 = np.where(full.obs["batch"].astype(str).to_numpy() == "PBMC5k")[0]
    src_idx = np.sort(rng.choice(b10, size=250, replace=False))
    tgt_idx = np.sort(rng.choice(b5, size=200, replace=False))
    source = _csr_copy(full[src_idx])
    target_scvi = _csr_copy(full[tgt_idx])
    target_total = _csr_copy(full[tgt_idx])
    del full

    tests: dict[str, Any] = {
        "scvi_class": _inspect_class(SCVI),
        "totalvi_class": _inspect_class(TOTALVI),
    }

    SCVI.setup_anndata(source, layer="counts", batch_key="batch")
    scvi_ref = SCVI(source, n_latent=4, n_hidden=32, n_layers=1, gene_likelihood="nb")
    scvi_ref.train(max_epochs=2, accelerator="gpu", devices=[0], batch_size=64, early_stopping=False)
    z_src = np.asarray(scvi_ref.get_latent_representation())
    SCVI.prepare_query_anndata(target_scvi, scvi_ref)
    scvi_q = SCVI.load_query_data(target_scvi, scvi_ref, accelerator="gpu", device="auto")
    trainable_scvi = _trainable_param_names(scvi_q)
    scvi_q.train(
        max_epochs=2,
        accelerator="gpu",
        devices=[0],
        batch_size=64,
        early_stopping=False,
        plan_kwargs={"weight_decay": 0.0},
    )
    z_tgt = np.asarray(scvi_q.get_latent_representation())
    z_src_adapted = np.asarray(scvi_q.get_latent_representation(source.copy()))
    tests["scvi_query"] = {
        "ok": True,
        "adapted_model_can_encode_source": bool(
            z_src_adapted.shape == z_src.shape and np.isfinite(z_src_adapted).all()
        ),
        "source_drift_under_adaptation": float(
            np.linalg.norm(z_src_adapted - z_src, axis=1).mean()
        ),
        "source_latent_shape": list(z_src.shape),
        "target_latent_shape": list(z_tgt.shape),
        "source_finite": bool(np.isfinite(z_src).all()),
        "target_finite": bool(np.isfinite(z_tgt).all()),
        "same_latent_dim": z_src.shape[1] == z_tgt.shape[1],
        "trainable_parameter_names": trainable_scvi,
        "n_trainable": len(trainable_scvi),
        "query_parameter_device": str(next(scvi_q.module.parameters()).device),
        "new_batch_registered": "PBMC5k" in set(target_scvi.obs["batch"].astype(str)),
        "zero_shot": False,
        "procedure": "scArches query adaptation via prepare_query_anndata + load_query_data + train",
    }
    del scvi_ref, scvi_q

    TOTALVI.setup_anndata(
        source,
        protein_expression_obsm_key="protein_expression",
        batch_key="batch",
        layer="counts",
    )
    tot_ref = TOTALVI(
        source,
        n_latent=4,
        n_hidden=32,
        n_layers_encoder=1,
        n_layers_decoder=1,
        gene_likelihood="nb",
        empirical_protein_background_prior=False,
    )
    tot_ref.train(max_epochs=2, accelerator="gpu", devices=[0], batch_size=64, early_stopping=False)
    z_src_t = np.asarray(tot_ref.get_latent_representation())
    TOTALVI.prepare_query_anndata(target_total, tot_ref)
    tot_q = TOTALVI.load_query_data(target_total, tot_ref, accelerator="gpu", device="auto")
    trainable_tot = _trainable_param_names(tot_q)
    tot_q.train(
        max_epochs=2,
        accelerator="gpu",
        devices=[0],
        batch_size=64,
        early_stopping=False,
        plan_kwargs={"weight_decay": 0.0},
    )
    z_tgt_t = np.asarray(tot_q.get_latent_representation())
    z_src_t_adapted = np.asarray(tot_q.get_latent_representation(source.copy()))
    rna_n, prot_n = tot_q.get_normalized_expression(return_mean=True, n_samples=1)
    tests["totalvi_query"] = {
        "ok": True,
        "adapted_model_can_encode_source": bool(
            z_src_t_adapted.shape == z_src_t.shape and np.isfinite(z_src_t_adapted).all()
        ),
        "source_drift_under_adaptation": float(
            np.linalg.norm(z_src_t_adapted - z_src_t, axis=1).mean()
        ),
        "source_latent_shape": list(z_src_t.shape),
        "target_latent_shape": list(z_tgt_t.shape),
        "source_finite": bool(np.isfinite(z_src_t).all()),
        "target_finite": bool(np.isfinite(z_tgt_t).all()),
        "protein_norm_finite": bool(np.isfinite(np.asarray(prot_n)).all()),
        "same_latent_dim": z_src_t.shape[1] == z_tgt_t.shape[1],
        "trainable_parameter_names": trainable_tot,
        "n_trainable": len(trainable_tot),
        "query_parameter_device": str(next(tot_q.module.parameters()).device),
        "zero_shot": False,
        "procedure": "scArches query adaptation via prepare_query_anndata + load_query_data + train",
        "query_uses_protein": True,
    }
    del tot_ref, tot_q, rna_n, prot_n

    tests["comparable"] = bool(
        tests["scvi_class"]["arches_mixin"]
        and tests["totalvi_class"]["arches_mixin"]
        and tests["scvi_query"]["ok"]
        and tests["totalvi_query"]["ok"]
        and tests["scvi_query"]["same_latent_dim"]
        and tests["totalvi_query"]["same_latent_dim"]
        and tests["scvi_query"]["adapted_model_can_encode_source"]
        and tests["totalvi_query"]["adapted_model_can_encode_source"]
    )
    tests["terminology"] = (
        "cross-dataset transfer / query generalization (scArches). "
        "Not strict zero-shot: unlabeled target RNA (and protein for totalVI) "
        "update new-batch / unfrozen parameters."
    )
    tests["zero_shot_supported"] = False
    tests["zero_shot_reason"] = (
        "load_query_data expands batch embeddings for unseen categories with "
        "random padding. Embedding query cells without adaptation would use "
        "untrained batch parameters. Mapping query cells as the source batch "
        "would ignore dataset identity. Neither is a documented zero-shot mode."
    )
    return {
        "scvi_tools": scvi.__version__,
        "torch": torch.__version__,
        "cuda": torch.cuda.is_available(),
        "gpu": torch.cuda.get_device_name(0),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "tests": tests,
        "phase9_note": (
            "PHASE 9 was not empirically executed: totalVI 1.3.3 has no honest "
            "cell-level missing-protein mask. That is an API limitation, not a "
            "failed scientific result."
        ),
    }


def main() -> None:
    result = run_api_validation()
    out = resolve_path("results/phase10_cross_dataset/logs")
    save_json(result, out / "phase10_query_api.json")
    md = out / "phase10_query_api.md"
    tests = result["tests"]
    md.write_text(
        f"""# PHASE 10 query API validation (scvi-tools {result['scvi_tools']})

Date (UTC): {result['timestamp_utc']}
GPU: {result['gpu']}

## SCVI
- ArchesMixin: {tests['scvi_class']['arches_mixin']}
- prepare_query_anndata / load_query_data: present
- tiny query test: {tests['scvi_query']['ok']}
- trainable params after freeze: {tests['scvi_query']['n_trainable']}
- adapted model can re-encode source cells: {tests['scvi_query']['adapted_model_can_encode_source']}
- mean source latent drift under adaptation: {tests['scvi_query']['source_drift_under_adaptation']:.4f}
- device: {tests['scvi_query']['query_parameter_device']}

## TOTALVI
- ArchesMixin: {tests['totalvi_class']['arches_mixin']}
- prepare_query_anndata / load_query_data: present
- tiny query test: {tests['totalvi_query']['ok']}
- trainable params after freeze: {tests['totalvi_query']['n_trainable']}
- query uses protein: {tests['totalvi_query']['query_uses_protein']}
- adapted model can re-encode source cells: {tests['totalvi_query']['adapted_model_can_encode_source']}
- mean source latent drift under adaptation: {tests['totalvi_query']['source_drift_under_adaptation']:.4f}
- device: {tests['totalvi_query']['query_parameter_device']}

## Terminology
{tests['terminology']}

Zero-shot secondary experiment: **not added**. {tests['zero_shot_reason']}

Comparable SCVI vs TOTALVI query workflow: **{tests['comparable']}**

{result['phase9_note']}
""",
        encoding="utf-8",
    )
    print(md.read_text(encoding="utf-8"))
    if not tests["comparable"]:
        raise SystemExit("SCVI and TOTALVI query workflows are not comparable. STOP.")


if __name__ == "__main__":
    main()
