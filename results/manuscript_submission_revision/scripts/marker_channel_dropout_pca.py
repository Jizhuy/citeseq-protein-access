#!/usr/bin/env python3
"""Marker-channel dropout for protein/concat PCA (RNA fit once; totalVI not retrained)."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data/processed/pbmc_cite_combined_inner_annotated.h5ad"
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 0
N_HVG = 2000
N_RNA = 10
N_PROT = 10
SINGLE = ["CD3", "CD4", "CD8a", "CD14", "CD19", "CD56"]
PANELS = {
    "T_lineage": ["CD3", "CD4", "CD8a", "CD45RA", "CD45RO"],
    "myeloid": ["CD14", "CD16", "CD15"],
    "B_NK": ["CD19", "CD56"],
}


def clr_rows(X: np.ndarray, eps: float = 1.0) -> np.ndarray:
    X = np.clip(X, 0, None) + eps
    logX = np.log(X)
    return logX - logX.mean(axis=1, keepdims=True)


def protein_arr(adata) -> tuple[np.ndarray, list[str]]:
    p = adata.obsm["protein_expression"]
    if isinstance(p, pd.DataFrame):
        return p.to_numpy(dtype=np.float64), [str(c) for c in p.columns]
    return np.asarray(p, dtype=np.float64), [f"p{i}" for i in range(p.shape[1])]


def clean_name(n: str) -> str:
    s = str(n)
    for suf in ("_TotalSeqB", "_TotalSeqBE", "_TotalSeqBF", "_TotalSeqA"):
        s = s.replace(suf, "")
    return s


def resolve(names, wanted):
    cleaned = [clean_name(n) for n in names]
    idxs, labs, miss = [], [], []
    for w in wanted:
        hits = [i for i, c in enumerate(cleaned) if c.lower() == w.lower()]
        if not hits:
            hits = [i for i, c in enumerate(cleaned) if w.lower() in c.lower()]
        if hits:
            idxs.append(hits[0])
            labs.append(names[hits[0]])
        else:
            miss.append(w)
    return idxs, labs, miss


def fit_rna_pca(adata_fit, adata_all):
    """Train-only HVG/norm/scale/PCA; transform all cells."""
    af = adata_fit.copy()
    sc.pp.highly_variable_genes(af, n_top_genes=N_HVG, flavor="seurat_v3", subset=False)
    hvgs = af.var_names[af.var["highly_variable"]].tolist()
    af = af[:, hvgs].copy()
    sc.pp.normalize_total(af, target_sum=1e4)
    sc.pp.log1p(af)
    sc.pp.scale(af, max_value=10)
    Xf = np.asarray(af.X)
    n = min(N_RNA, Xf.shape[0] - 1, Xf.shape[1])
    pca = PCA(n_components=n, random_state=SEED)
    Z_fit = pca.fit_transform(Xf)

    aa = adata_all[:, hvgs].copy()
    sc.pp.normalize_total(aa, target_sum=1e4)
    sc.pp.log1p(aa)
    # apply fit scale params from af
    means = np.asarray(af.var["mean"]) if "mean" in af.var else Xf.mean(axis=0)
    # scanpy stores mean/std in var after scale
    if "mean" in af.var.columns and "std" in af.var.columns:
        mu = af.var["mean"].to_numpy()
        sd = af.var["std"].to_numpy()
        sd = np.where(sd < 1e-8, 1.0, sd)
        Xt = np.asarray(aa.X)
        Xt = (Xt - mu) / sd
        Xt = np.clip(Xt, -10, 10)
    else:
        sc.pp.scale(aa, max_value=10)
        Xt = np.asarray(aa.X)
    Z_all = pca.transform(Xt)
    return Z_fit, Z_all


def protein_pca(P_fit, P_all):
    Cf = clr_rows(P_fit)
    mu = Cf.mean(axis=0)
    sd = Cf.std(axis=0, ddof=1)
    sd = np.where(sd < 1e-8, 1.0, sd)
    Zf = (Cf - mu) / sd
    Zt = (clr_rows(P_all) - mu) / sd
    n = min(N_PROT, Zf.shape[0] - 1, Zf.shape[1])
    pca = PCA(n_components=n, random_state=SEED)
    return pca.fit_transform(Zf), pca.transform(Zt)


def f1_logreg(Ztr, ytr, Zte, yte):
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=SEED)
    clf.fit(Ztr, ytr)
    return float(f1_score(yte, clf.predict(Zte), average="macro", zero_division=0))


def f1_concat(Zr_tr, Zp_tr, Zr_te, Zp_te, ytr, yte):
    scx = StandardScaler()
    Xtr = scx.fit_transform(np.hstack([Zr_tr, Zp_tr]))
    Xte = scx.transform(np.hstack([Zr_te, Zp_te]))
    return f1_logreg(Xtr, ytr, Xte, yte)


def main():
    print("Loading", H5AD, flush=True)
    adata = ad.read_h5ad(H5AD)
    high = adata.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    hc = adata[high].copy()
    del adata
    split = hc.obs["split"].astype(str).to_numpy()
    tr = split == "train"
    te = split == "test"
    y = hc.obs["cell_type_l2"].astype(str).to_numpy()
    train = hc[tr].copy()
    print(f"HC={hc.n_obs} train={tr.sum()} test={te.sum()}", flush=True)

    P_all, names = protein_arr(hc)
    P_tr = P_all[tr]
    print("proteins:", names, flush=True)

    print("RNA PCA...", flush=True)
    Zrna_tr, Zrna_all = fit_rna_pca(train, hc)
    f1_rna = f1_logreg(Zrna_all[tr], y[tr], Zrna_all[te], y[te])

    Zp_tr0, Zp_all0 = protein_pca(P_tr, P_all)
    f1_p0 = f1_logreg(Zp_tr0, y[tr], Zp_all0[te], y[te])
    f1_c0 = f1_concat(Zrna_all[tr], Zp_tr0, Zrna_all[te], Zp_all0[te], y[tr], y[te])
    print(f"clean rna={f1_rna:.4f} prot={f1_p0:.4f} concat={f1_c0:.4f}", flush=True)

    rows = [
        {
            "setting": "clean",
            "dropped_markers": "",
            "n_channels_dropped": 0,
            "F1_RNA_PCA": f1_rna,
            "F1_protein_PCA": f1_p0,
            "F1_concat": f1_c0,
            "missing_requested": "",
        }
    ]
    jobs = [(f"drop_{m}", [m]) for m in SINGLE] + [(f"panel_{k}", v) for k, v in PANELS.items()]
    for setting, wanted in jobs:
        idxs, labs, miss = resolve(names, wanted)
        print(f"{setting}: {labs} miss={miss}", flush=True)
        Ptr = P_tr.copy()
        Pall = P_all.copy()
        for i in idxs:
            Ptr[:, i] = 0.0
            Pall[:, i] = 0.0
        Zp_tr, Zp_all = protein_pca(Ptr, Pall)
        fp = f1_logreg(Zp_tr, y[tr], Zp_all[te], y[te])
        fc = f1_concat(Zrna_all[tr], Zp_tr, Zrna_all[te], Zp_all[te], y[tr], y[te])
        rows.append(
            {
                "setting": setting,
                "dropped_markers": ";".join(labs),
                "n_channels_dropped": len(idxs),
                "F1_RNA_PCA": f1_rna,
                "F1_protein_PCA": fp,
                "F1_concat": fc,
                "missing_requested": ";".join(miss),
            }
        )

    out = pd.DataFrame(rows)
    out["delta_concat_vs_clean"] = out["F1_concat"] - f1_c0
    out["delta_protein_vs_clean"] = out["F1_protein_PCA"] - f1_p0
    out.to_csv(OUT / "marker_channel_dropout_pca_development.csv", index=False)
    (OUT / "marker_channel_dropout_meta.json").write_text(
        json.dumps(
            {
                "design": "zero protein channels before CLR/z/PCA; RNA PCA fixed from clean data",
                "scope": "development HC; protein/concat PCA only; totalVI NOT retrained",
                "protein_columns": names,
                "seed": SEED,
            },
            indent=2,
        )
    )
    print(out.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
