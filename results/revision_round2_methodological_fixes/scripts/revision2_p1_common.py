#!/usr/bin/env python3
"""Shared helpers for Revision Round 2 P1 controls."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
)
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]
CTRL = ROOT / "results" / "revision_round2_methodological_fixes" / "controls"
TABLES = CTRL / "tables"
LOGS = CTRL / "logs"
FIGS = CTRL / "figures_preview"

N_HVG = 2000
PRIMARY_RNA_PCS = 10
PRIMARY_PROT_PCS = 10
SECONDARY_RNA_PCS = 20
PROTEIN_OBSM = "protein_expression"
SEED = 0


def clr_rows(matrix: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    values = np.clip(values, 0.0, None) + float(pseudocount)
    logx = np.log(values)
    return logx - logx.mean(axis=1, keepdims=True)


def _to_dense(X) -> np.ndarray:
    if sparse.issparse(X):
        return np.asarray(X.toarray(), dtype=np.float64)
    return np.asarray(X, dtype=np.float64)


def counts_matrix(adata: ad.AnnData) -> np.ndarray:
    if "counts" in adata.layers:
        return _to_dense(adata.layers["counts"])
    return _to_dense(adata.X)


def protein_matrix(adata: ad.AnnData) -> tuple[np.ndarray, list[str]]:
    prot = adata.obsm[PROTEIN_OBSM]
    if isinstance(prot, pd.DataFrame):
        return prot.to_numpy(dtype=np.float64), [str(c) for c in prot.columns]
    return np.asarray(prot, dtype=np.float64), [f"p{i}" for i in range(prot.shape[1])]


def fit_rna_pca(
    adata_fit: ad.AnnData,
    adata_transform: ad.AnnData,
    n_pcs: int,
    n_hvg: int = N_HVG,
) -> dict[str, Any]:
    """Fit HVG+normalize+log1p+scale+PCA on fit cells; transform transform cells."""
    # Build fit RNA AnnData from counts
    Xf = counts_matrix(adata_fit)
    rna_fit = ad.AnnData(X=Xf)
    rna_fit.obs_names = adata_fit.obs_names
    rna_fit.var_names = adata_fit.var_names
    sc.pp.highly_variable_genes(rna_fit, n_top_genes=n_hvg, flavor="seurat_v3")
    hvgs = rna_fit.var_names[rna_fit.var["highly_variable"]].astype(str).tolist()
    rna_fit = rna_fit[:, hvgs].copy()
    sc.pp.normalize_total(rna_fit, target_sum=1e4)
    sc.pp.log1p(rna_fit)
    sc.pp.scale(rna_fit, max_value=10)
    X_fit = _to_dense(rna_fit.X)
    n_pcs_eff = min(n_pcs, X_fit.shape[0] - 1, X_fit.shape[1])
    pca = PCA(n_components=n_pcs_eff, random_state=SEED)
    Z_fit = pca.fit_transform(X_fit)

    # Transform
    Xt = counts_matrix(adata_transform)
    rna_t = ad.AnnData(X=Xt)
    rna_t.obs_names = adata_transform.obs_names
    rna_t.var_names = adata_transform.var_names
    missing = [g for g in hvgs if g not in set(rna_t.var_names.astype(str))]
    if missing:
        raise ValueError(f"{len(missing)} HVGs missing in transform object")
    rna_t = rna_t[:, hvgs].copy()
    sc.pp.normalize_total(rna_t, target_sum=1e4)
    sc.pp.log1p(rna_t)
    # Apply fit scaling: use stored means/vars from scanpy scale on fit
    # scanpy scale stores mean/var in var; recompute from fit dense matrix
    means = X_fit.mean(axis=0)  # after log1p before scale? We scaled in place.
    # Better: transform using pca on similarly processed matrix with fit parameters.
    # Re-extract unscaled log-normalized then apply fit mean/std from pre-scale state.
    # Rebuild fit log-normalized without scale for parameter capture:
    rna_fit2 = ad.AnnData(X=counts_matrix(adata_fit))
    rna_fit2.obs_names = adata_fit.obs_names
    rna_fit2.var_names = adata_fit.var_names
    rna_fit2 = rna_fit2[:, hvgs].copy()
    sc.pp.normalize_total(rna_fit2, target_sum=1e4)
    sc.pp.log1p(rna_fit2)
    L_fit = _to_dense(rna_fit2.X)
    mu = L_fit.mean(axis=0)
    sd = L_fit.std(axis=0, ddof=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    L_fit_s = np.clip((L_fit - mu) / sd, -10, 10)
    pca = PCA(n_components=n_pcs_eff, random_state=SEED)
    Z_fit = pca.fit_transform(L_fit_s)

    rna_t2 = ad.AnnData(X=counts_matrix(adata_transform))
    rna_t2.obs_names = adata_transform.obs_names
    rna_t2.var_names = adata_transform.var_names
    rna_t2 = rna_t2[:, hvgs].copy()
    sc.pp.normalize_total(rna_t2, target_sum=1e4)
    sc.pp.log1p(rna_t2)
    L_t = _to_dense(rna_t2.X)
    L_t_s = np.clip((L_t - mu) / sd, -10, 10)
    Z_t = pca.transform(L_t_s)
    return {
        "Z_fit": Z_fit,
        "Z_transform": Z_t,
        "hvgs": hvgs,
        "n_pcs": int(n_pcs_eff),
        "pca": pca,
        "rna_mu": mu,
        "rna_sd": sd,
    }


def fit_protein_pca(
    adata_fit: ad.AnnData,
    adata_transform: ad.AnnData,
    n_pcs: int,
) -> dict[str, Any]:
    Pf, names = protein_matrix(adata_fit)
    Pt, _ = protein_matrix(adata_transform)
    Cf = clr_rows(Pf)
    mu = Cf.mean(axis=0)
    sd = Cf.std(axis=0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    Zf_in = (Cf - mu) / sd
    Ct = clr_rows(Pt)
    Zt_in = (Ct - mu) / sd
    n_pcs_eff = min(n_pcs, Zf_in.shape[0] - 1, Zf_in.shape[1])
    pca = PCA(n_components=n_pcs_eff, random_state=SEED)
    Z_fit = pca.fit_transform(Zf_in)
    Z_t = pca.transform(Zt_in)
    return {
        "Z_fit": Z_fit,
        "Z_transform": Z_t,
        "protein_names": names,
        "n_pcs": int(n_pcs_eff),
        "prot_mu": mu,
        "prot_sd": sd,
        "pca": pca,
    }


def concat_and_scale(
    rna_fit: np.ndarray,
    prot_fit: np.ndarray,
    rna_t: np.ndarray,
    prot_t: np.ndarray,
) -> dict[str, Any]:
    X_fit = np.hstack([rna_fit, prot_fit])
    X_t = np.hstack([rna_t, prot_t])
    scaler = StandardScaler()
    Z_fit = scaler.fit_transform(X_fit)
    Z_t = scaler.transform(X_t)
    return {"Z_fit": Z_fit, "Z_transform": Z_t, "scaler": scaler, "dim": int(Z_fit.shape[1])}


def logreg_metrics(z_train, y_train, z_test, y_test, seed: int = SEED) -> dict[str, Any]:
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    clf.fit(z_train, y_train)
    pred = clf.predict(z_test)
    labels = sorted(set(y_test) | set(pred))
    per_class = {}
    for c in labels:
        per_class[c] = float(f1_score(y_test == c, pred == c, zero_division=0))
    return {
        "accuracy": float(accuracy_score(y_test, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, pred)),
        "macro_f1": float(f1_score(y_test, pred, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_test, pred, average="weighted", zero_division=0)),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "n_classes_test": int(len(set(y_test))),
        "per_class_f1": per_class,
        "y_test": y_test,
        "y_pred": pred,
        "clf": clf,
    }


def write_done(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def done_ok(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        data = json.loads(path.read_text())
        return bool(data.get("ok"))
    except Exception:
        return False


def align_latent(cells_csv: Path, latent_npy: Path, obs_names: pd.Index) -> np.ndarray:
    cells = pd.read_csv(cells_csv)["cell_id"].astype(str).tolist()
    Z = np.load(latent_npy)
    if len(cells) != Z.shape[0]:
        raise ValueError(f"cell/latent mismatch {len(cells)} vs {Z.shape}")
    idx = pd.Index(cells)
    order = [idx.get_loc(c) for c in obs_names.astype(str)]
    return Z[np.asarray(order)]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
