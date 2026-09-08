"""Independent Seurat v4 PBMC reference mapping. Does not use evaluated latents."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse

from src.models.run_scvi import requested_accelerator
from src.utils.device import detect_compute_environment
from src.utils.io import resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

UNLABELED = "Unknown"
REFERENCE_LABELS = ("celltype.l1", "celltype.l2", "celltype.l3")
N_HVG = 2000


def download_seurat_v4_reference(save_dir) -> ad.AnnData:
    """Load the official scvi-tools Seurat v4 CITE-seq object with published filters."""
    save_dir = resolve_path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    return scvi.data.pbmc_seurat_v4_cite_seq(
        save_path=str(save_dir),
        apply_filters=True,
        aggregate_proteins=True,
        mask_protein_batches=0,
    )


def inspect_reference(ref: ad.AnnData) -> dict[str, Any]:
    obs = [str(c) for c in ref.obs.columns]
    report: dict[str, Any] = {
        "n_cells": int(ref.n_obs),
        "n_genes": int(ref.n_vars),
        "n_proteins": int(ref.obsm["protein_counts"].shape[1]) if "protein_counts" in ref.obsm else 0,
        "obs_columns": obs,
        "obsm_keys": list(map(str, ref.obsm.keys())),
        "uns_keys": list(map(str, ref.uns.keys())),
        "layers": list(map(str, ref.layers.keys())),
        "var_name_examples": [str(v) for v in ref.var_names[:8]],
        "gene_id_kind": _gene_id_kind(ref.var_names),
        "x_is_sparse": bool(sparse.issparse(ref.X)),
    }
    if sparse.issparse(ref.X):
        data = np.asarray(ref.X.data)
    else:
        data = np.asarray(ref.X)
    report["x_min"] = float(np.min(data)) if data.size else None
    report["x_max"] = float(np.max(data)) if data.size else None
    report["x_approx_integer"] = bool(np.all(np.abs(data - np.round(data)) < 1e-6)) if data.size else False
    report["raw_rna_counts_available"] = bool(report["x_approx_integer"] and (report["x_min"] or 0) >= 0)
    for key in REFERENCE_LABELS:
        report[f"has_{key}"] = key in ref.obs
        if key in ref.obs:
            counts = ref.obs[key].astype(str).value_counts()
            report[f"{key}_n_classes"] = int(counts.shape[0])
            report[f"{key}_classes"] = counts.to_dict()
    for col in ("donor", "Donor", "orig.ident", "lane", "Phase", "time", "batch"):
        if col in ref.obs:
            report[f"{col}_n_unique"] = int(ref.obs[col].astype(str).nunique())
            report[f"{col}_examples"] = ref.obs[col].astype(str).value_counts().head(8).to_dict()
    return report


def _gene_id_kind(names: pd.Index) -> str:
    sample = [str(x) for x in names[:50]]
    ensembl = sum(s.startswith("ENSG") for s in sample)
    symbolish = sum(any(c.isalpha() for c in s) and not s.startswith("ENSG") for s in sample)
    if ensembl > 40:
        return "ensembl"
    if symbolish > 30:
        return "gene_symbol"
    return "unknown"


def compare_gene_identifiers(query: ad.AnnData, reference: ad.AnnData) -> dict[str, Any]:
    q = pd.Index(query.var_names.astype(str))
    r = pd.Index(reference.var_names.astype(str))
    q_dups = q[q.duplicated()].unique().tolist()
    r_dups = r[r.duplicated()].unique().tolist()
    shared = q.intersection(r)
    report = {
        "query_n_genes": int(len(q)),
        "reference_n_genes": int(len(r)),
        "query_gene_id_kind": _gene_id_kind(q),
        "reference_gene_id_kind": _gene_id_kind(r),
        "query_examples": q[:8].tolist(),
        "reference_examples": r[:8].tolist(),
        "n_exact_intersection": int(len(shared)),
        "fraction_query_shared": float(len(shared) / len(q)) if len(q) else 0.0,
        "fraction_reference_shared": float(len(shared) / len(r)) if len(r) else 0.0,
        "query_duplicated_genes": q_dups[:20],
        "reference_duplicated_genes": r_dups[:20],
        "n_query_duplicated": int(len(q_dups)),
        "n_reference_duplicated": int(len(r_dups)),
        "compatible": (
            _gene_id_kind(q) == _gene_id_kind(r)
            and _gene_id_kind(q) != "unknown"
            and len(shared) >= 1000
        ),
    }
    return report


def label_distribution(ref: ad.AnnData, key: str) -> pd.DataFrame:
    counts = ref.obs[key].astype(str).value_counts()
    out = counts.rename_axis("cell_type").reset_index(name="n_cells")
    out["proportion"] = out["n_cells"] / out["n_cells"].sum()
    out["rare"] = out["n_cells"] < 100
    return out


def l2_to_l1_map(ref: ad.AnnData) -> pd.Series:
    table = (
        ref.obs[["celltype.l2", "celltype.l1"]]
        .astype(str)
        .groupby("celltype.l2")["celltype.l1"]
        .agg(lambda s: s.value_counts().idxmax())
    )
    return table


def _counts_matrix(adata: ad.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def select_reference_hvgs(
    reference: ad.AnnData,
    shared_genes: pd.Index,
    n_top: int = N_HVG,
    seed: int = 0,
) -> list[str]:
    """HVGs from the reference only, restricted to genes shared with the query."""
    keep = reference.var_names.isin(shared_genes)
    tmp = reference[:, keep].copy()
    tmp.X = _counts_matrix(tmp)
    if tmp.var_names.duplicated().any():
        tmp = tmp[:, ~tmp.var_names.duplicated()].copy()
    try:
        sc.pp.highly_variable_genes(tmp, n_top_genes=n_top, flavor="seurat_v3", span=0.3)
    except ImportError:
        # seurat_v3 needs scikit-misc; fall back to the documented Seurat flavor on log-normalized counts
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(tmp, n_top_genes=n_top, flavor="seurat")
    hvgs = tmp.var_names[tmp.var["highly_variable"]].astype(str).tolist()
    if len(hvgs) < 500:
        raise ValueError(f"Only {len(hvgs)} HVGs selected; gene overlap is insufficient.")
    return hvgs


def _subset_for_scanvi(adata: ad.AnnData, genes: list[str], *, is_query: bool) -> ad.AnnData:
    out = adata.copy()
    missing = [g for g in genes if g not in out.var_names]
    if missing:
        raise ValueError(f"{len(missing)} annotation genes missing from {'query' if is_query else 'reference'}.")
    out = out[:, genes].copy()
    counts = _counts_matrix(out)
    if sparse.issparse(counts):
        out.X = counts.tocsr()
    else:
        out.X = sparse.csr_matrix(np.asarray(counts))
    if "counts" in out.layers:
        del out.layers["counts"]
    return out


def train_reference_scanvi(
    reference: ad.AnnData,
    genes: list[str],
    model_dir,
    *,
    labels_key: str = "celltype.l2",
    batch_key: str,
    seed: int = 0,
    max_epochs_scvi: int = 80,
    max_epochs_scanvi: int = 20,
) -> tuple[scvi.model.SCANVI, dict[str, Any]]:
    """Train a NEW reference SCANVI model. Never loads PHASE 3/5 weights."""
    set_global_seed(seed)
    logger = get_logger("annotation_reference", resolve_path("results/logs") / "annotation_reference.log")
    env = detect_compute_environment()
    requested, detail = requested_accelerator(env)
    ref = _subset_for_scanvi(reference, genes, is_query=False)
    if labels_key not in ref.obs:
        raise KeyError(f"Reference is missing {labels_key}")
    ref.obs[labels_key] = ref.obs[labels_key].astype(str)
    ref.obs[batch_key] = ref.obs[batch_key].astype(str)
    logger.info(
        "Training annotation-only SCVI/SCANVI on reference: %s cells x %s genes, labels=%s, batch=%s, device=%s",
        ref.n_obs,
        ref.n_vars,
        labels_key,
        batch_key,
        requested,
    )
    scvi.model.SCVI.setup_anndata(ref, batch_key=batch_key)
    scvi_model = scvi.model.SCVI(
        ref,
        n_hidden=128,
        n_latent=30,
        n_layers=2,
        dropout_rate=0.1,
        gene_likelihood="nb",
        dispersion="gene",
    )
    scvi_model.train(
        max_epochs=max_epochs_scvi,
        accelerator=requested,
        devices="auto",
        batch_size=256,
        train_size=0.9,
        early_stopping=True,
        early_stopping_patience=15,
        check_val_every_n_epoch=1,
    )
    scanvi = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=UNLABELED,
        labels_key=labels_key,
    )
    scanvi.train(
        max_epochs=max_epochs_scanvi,
        n_samples_per_label=100,
        accelerator=requested,
        devices="auto",
        batch_size=256,
        train_size=0.9,
        check_val_every_n_epoch=1,
    )
    model_dir = resolve_path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    scanvi.save(str(model_dir), overwrite=True, save_anndata=False)
    info = {
        "model_dir": str(model_dir),
        "n_cells": int(ref.n_obs),
        "n_genes": int(ref.n_vars),
        "labels_key": labels_key,
        "batch_key": batch_key,
        "requested_device": requested,
        "requested_device_detail": detail,
        "n_latent": 30,
        "n_layers": 2,
        "gene_likelihood": "nb",
        "max_epochs_scvi": max_epochs_scvi,
        "max_epochs_scanvi": max_epochs_scanvi,
        "initialized_from_query_models": False,
        "uses_evaluated_latents": False,
        "scvi_tools_version": env["scvi_version"],
        "torch_version": env["torch_version"],
    }
    save_json(info, model_dir / "annotation_model_info.json")
    return scanvi, info


def transfer_labels(
    query: ad.AnnData,
    scanvi: scvi.model.SCANVI,
    genes: list[str],
    l2_to_l1: pd.Series,
    *,
    batch_key: str,
    query_surgery_epochs: int = 50,
    seed: int = 0,
) -> pd.DataFrame:
    """Map query RNA into the reference SCANVI model. Surgery uses query RNA only."""
    set_global_seed(seed)
    logger = get_logger("annotation_query", resolve_path("results/logs") / "annotation_query.log")
    env = detect_compute_environment()
    requested, _ = requested_accelerator(env)
    q = _subset_for_scanvi(query, genes, is_query=True)
    q.obs["celltype.l2"] = UNLABELED
    q.obs[batch_key] = query.obs.loc[q.obs_names, "batch"].astype(str)
    scvi.model.SCANVI.prepare_query_anndata(q, scanvi, inplace=True)
    q_model = scvi.model.SCANVI.load_query_data(q, scanvi, freeze_classifier=True)
    logger.info("Running scArches query surgery on %s query cells; classifier frozen.", q.n_obs)
    q_model.train(
        max_epochs=query_surgery_epochs,
        accelerator=requested,
        devices="auto",
        batch_size=256,
        plan_kwargs={"weight_decay": 0.0},
    )
    soft = q_model.predict(soft=True)
    if isinstance(soft, tuple):
        soft = soft[0]
    soft = pd.DataFrame(soft)
    soft.index = q.obs_names
    pred_l2 = soft.idxmax(axis=1).astype(str)
    conf_l2 = soft.max(axis=1).astype(float)
    pred_l1 = pred_l2.map(l2_to_l1).astype(str)
    conf_l1 = []
    for cell_id, l1 in pred_l1.items():
        members = [c for c in soft.columns if l2_to_l1.get(str(c), None) == l1]
        conf_l1.append(float(soft.loc[cell_id, members].sum()) if members else float(conf_l2.loc[cell_id]))
    table = pd.DataFrame(
        {
            "cell_id": q.obs_names.astype(str),
            "batch": query.obs.loc[q.obs_names, "batch"].astype(str).to_numpy(),
            "predicted_celltype_l1": pred_l1.to_numpy(),
            "predicted_celltype_l2": pred_l2.to_numpy(),
            "prediction_confidence_l1": np.asarray(conf_l1),
            "prediction_confidence_l2": conf_l2.to_numpy(),
            "annotation_method": "SCANVI_scArches_SeuratV4",
            "reference_source": "scvi.data.pbmc_seurat_v4_cite_seq; Stuart et al. Cell 2021",
        }
    )
    return table


def choose_confidence_tiers(scores: np.ndarray) -> dict[str, Any]:
    """Choose tiers from the observed score distribution. Do not invent probabilities."""
    scores = np.asarray(scores, dtype=float)
    q50, q75, q90 = np.quantile(scores, [0.5, 0.75, 0.9])
    # SCANVI max-class probabilities are typically high; 0.5 is not used.
    high_cut = 0.85 if q75 >= 0.85 else float(q75)
    if (scores >= high_cut).mean() < 0.40:
        high_cut = float(q75)
    medium_cut = min(float(q50), high_cut)
    if medium_cut >= high_cut:
        medium_cut = max(high_cut - 0.10, float(np.quantile(scores, 0.25)))
    return {
        "high_threshold": float(high_cut),
        "medium_threshold": float(medium_cut),
        "quantiles": {
            "q25": float(np.quantile(scores, 0.25)),
            "q50": float(q50),
            "q75": float(q75),
            "q90": float(q90),
            "min": float(scores.min()),
            "max": float(scores.max()),
            "mean": float(scores.mean()),
        },
        "rule": (
            f"high: l2_confidence >= {high_cut:.3f}; "
            f"medium: [{medium_cut:.3f}, {high_cut:.3f}); "
            f"low: < {medium_cut:.3f}. "
            "Thresholds taken from the observed SCANVI max-class probability distribution."
        ),
        "justification": (
            "0.5 was not used because the observed score mass is far above 0.5. "
            "High confidence uses 0.85 when the 75th percentile supports it, otherwise the 75th percentile."
        ),
    }


def assign_confidence_tier(score: float, tiers: dict[str, Any]) -> str:
    if score >= tiers["high_threshold"]:
        return "high"
    if score >= tiers["medium_threshold"]:
        return "medium"
    return "low"
