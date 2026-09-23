#!/usr/bin/env python3
"""Strict train-only internal development evaluation (HC 7573 / 1921).

Trains scVI_matched and totalVI ONLY on high-confidence training cells,
encodes held-out HC test cells after fitting (no query adaptation training),
fits PCA/concat on train-only, then classifier on train labels only.

Preserves verified hyperparameters from MATCHED_HYPERPARAMETERS / totalVI defaults.
"""
from __future__ import annotations

import json
import sys
import warnings
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler
from scvi.model import SCVI, TOTALVI

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "results/revision_round2_methodological_fixes/scripts"))
from revision2_p1_common import (  # noqa: E402
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    SEED,
    concat_and_scale,
    fit_protein_pca,
    fit_rna_pca,
)

OUT = Path(__file__).resolve().parents[1] / "strict_train_only"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "models").mkdir(exist_ok=True)
(OUT / "embeddings").mkdir(exist_ok=True)
(OUT / "tables").mkdir(exist_ok=True)
(OUT / "logs").mkdir(exist_ok=True)

N_BOOT = 10_000
BOOT_SEED = 20260923
STRAT_SEED = 20260924
CLF_SEED = 0

MATCHED = {
    "n_latent": 20,
    "n_hidden": 256,
    "n_layers": 2,
    "dropout_rate": 0.2,
    "dispersion": "gene",
    "gene_likelihood": "nb",
    "use_observed_lib_size": True,
    "latent_distribution": "normal",
    "train_size": 0.9,
    "max_epochs": 200,
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_mode": "min",
    "batch_size": 256,
}

TOTALVI_HP = {
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
    "max_epochs": 200,
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_mode": "min",
    "batch_size": 256,
}


def log(msg: str) -> None:
    line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
    print(line, flush=True)
    with open(OUT / "logs" / "strict_train_only.log", "a") as f:
        f.write(line + "\n")


def macro_f1(y_true, y_pred) -> float:
    return float(f1_score(y_true, y_pred, average="macro", zero_division=0))


def ensure_counts_csr(adata: ad.AnnData) -> ad.AnnData:
    a = adata.copy()
    if "counts" not in a.layers:
        # X is raw counts in this project object
        from scipy import sparse

        X = a.X
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        else:
            X = X.tocsr()
        a.layers["counts"] = X
    else:
        from scipy import sparse

        if not sparse.issparse(a.layers["counts"]):
            a.layers["counts"] = sparse.csr_matrix(a.layers["counts"])
        else:
            a.layers["counts"] = a.layers["counts"].tocsr()
    # protein must be DataFrame
    if "protein_expression" in a.obsm and not isinstance(a.obsm["protein_expression"], pd.DataFrame):
        raise TypeError("protein_expression must be DataFrame")
    return a


def accelerator() -> str:
    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "gpu"
    return "cpu"


def fit_predict(Z_tr, y_tr, Z_te, seed: int = CLF_SEED):
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    clf.fit(Z_tr, y_tr)
    return clf.predict(Z_te), clf


def encode_heldout(model, adata_query: ad.AnnData) -> np.ndarray:
    """Encode held-out cells with frozen trained weights (no query training).

    scvi-tools transfers AnnData setup automatically when adata= is passed to
    get_latent_representation on a trained model. This does not update weights.
    """
    return np.asarray(model.get_latent_representation(adata=adata_query))


def train_scvi_matched(adata_train: ad.AnnData, device: str) -> tuple[SCVI, Path]:
    model_dir = OUT / "models" / "scvi_matched_trainonly_seed0"
    if (model_dir / "model.pt").exists() or (model_dir / "attr.pkl").exists() or list(model_dir.glob("*")):
        # try load if complete
        try:
            log("Loading existing scVI_matched train-only model")
            model = SCVI.load(str(model_dir), adata=adata_train)
            return model, model_dir
        except Exception as e:
            log(f"Reload failed ({e}); retraining scVI_matched")

    SCVI.setup_anndata(adata_train, layer="counts", batch_key="batch")
    hp = MATCHED
    model = SCVI(
        adata_train,
        n_latent=hp["n_latent"],
        n_hidden=hp["n_hidden"],
        n_layers=hp["n_layers"],
        dropout_rate=hp["dropout_rate"],
        dispersion=hp["dispersion"],
        gene_likelihood=hp["gene_likelihood"],
        latent_distribution=hp["latent_distribution"],
        use_observed_lib_size=hp["use_observed_lib_size"],
    )
    log(f"Training scVI_matched on n={adata_train.n_obs} device={device}")
    try:
        model.train(
            max_epochs=hp["max_epochs"],
            accelerator=device,
            devices=1 if device != "cpu" else "auto",
            train_size=hp["train_size"],
            batch_size=hp["batch_size"],
            early_stopping=hp["early_stopping"],
            early_stopping_monitor=hp["early_stopping_monitor"],
            early_stopping_patience=hp["early_stopping_patience"],
            early_stopping_mode=hp["early_stopping_mode"],
            check_val_every_n_epoch=1,
        )
    except Exception as e:
        log(f"scVI_matched on {device} failed: {e}; retrying CPU")
        model.train(
            max_epochs=hp["max_epochs"],
            accelerator="cpu",
            train_size=hp["train_size"],
            batch_size=hp["batch_size"],
            early_stopping=hp["early_stopping"],
            early_stopping_monitor=hp["early_stopping_monitor"],
            early_stopping_patience=hp["early_stopping_patience"],
            early_stopping_mode=hp["early_stopping_mode"],
            check_val_every_n_epoch=1,
        )
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    log(f"Saved scVI_matched to {model_dir}")
    return model, model_dir


def train_totalvi(adata_train: ad.AnnData, device: str) -> tuple[TOTALVI, Path]:
    model_dir = OUT / "models" / "totalvi_trainonly_seed0"
    if list(model_dir.glob("*")):
        try:
            log("Loading existing totalVI train-only model")
            model = TOTALVI.load(str(model_dir), adata=adata_train)
            return model, model_dir
        except Exception as e:
            log(f"Reload failed ({e}); retraining totalVI")

    TOTALVI.setup_anndata(
        adata_train,
        protein_expression_obsm_key="protein_expression",
        batch_key="batch",
        layer="counts",
    )
    hp = TOTALVI_HP
    model = TOTALVI(
        adata_train,
        n_latent=hp["n_latent"],
        gene_dispersion=hp["gene_dispersion"],
        protein_dispersion=hp["protein_dispersion"],
        gene_likelihood=hp["gene_likelihood"],
        latent_distribution=hp["latent_distribution"],
        empirical_protein_background_prior=hp["empirical_protein_background_prior"],
        override_missing_proteins=hp["override_missing_proteins"],
        n_hidden=hp["n_hidden"],
        n_layers_encoder=hp["n_layers_encoder"],
        n_layers_decoder=hp["n_layers_decoder"],
        dropout_rate_encoder=hp["dropout_rate_encoder"],
        dropout_rate_decoder=hp["dropout_rate_decoder"],
    )
    log(f"Training totalVI on n={adata_train.n_obs} device={device}")
    try:
        model.train(
            max_epochs=hp["max_epochs"],
            lr=hp["lr"],
            accelerator=device,
            devices=1 if device != "cpu" else "auto",
            train_size=hp["train_size"],
            batch_size=hp["batch_size"],
            early_stopping=hp["early_stopping"],
            early_stopping_monitor=hp["early_stopping_monitor"],
            early_stopping_patience=hp["early_stopping_patience"],
            early_stopping_mode=hp["early_stopping_mode"],
            reduce_lr_on_plateau=hp["reduce_lr_on_plateau"],
            check_val_every_n_epoch=1,
        )
    except Exception as e:
        log(f"totalVI on {device} failed: {e}; retrying CPU")
        model.train(
            max_epochs=hp["max_epochs"],
            lr=hp["lr"],
            accelerator="cpu",
            train_size=hp["train_size"],
            batch_size=hp["batch_size"],
            early_stopping=hp["early_stopping"],
            early_stopping_monitor=hp["early_stopping_monitor"],
            early_stopping_patience=hp["early_stopping_patience"],
            early_stopping_mode=hp["early_stopping_mode"],
            reduce_lr_on_plateau=hp["reduce_lr_on_plateau"],
            check_val_every_n_epoch=1,
        )
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    log(f"Saved totalVI to {model_dir}")
    return model, model_dir


def bootstrap_paired(y, pred_rna, pred_cat, pred_tot, seed: int, stratified: bool = False):
    rng = np.random.default_rng(seed)
    n = len(y)
    rows = []
    classes = np.unique(y)
    for b in range(N_BOOT):
        if stratified:
            parts = []
            for lab in classes:
                idx = np.flatnonzero(y == lab)
                parts.append(rng.choice(idx, size=len(idx), replace=True))
            idx = np.concatenate(parts)
        else:
            idx = rng.choice(n, size=n, replace=True)
        f_rna = macro_f1(y[idx], pred_rna[idx])
        f_cat = macro_f1(y[idx], pred_cat[idx])
        f_tot = macro_f1(y[idx], pred_tot[idx])
        g_access = f_cat - f_rna
        g_assoc = f_tot - f_cat
        gap = f_tot - f_rna
        r = g_access / gap if abs(gap) > 1e-12 else np.nan
        rows.append(
            {
                "F1_RNA_PCA": f_rna,
                "F1_concat": f_cat,
                "F1_totalVI": f_tot,
                "G_access": g_access,
                "G_assoc": g_assoc,
                "Gap_total": gap,
                "R_access": r,
            }
        )
    return pd.DataFrame(rows)


def summarize_boot(df: pd.DataFrame, point: dict, method: str) -> pd.DataFrame:
    rows = []
    for key, pe in point.items():
        s = df[key]
        rows.append(
            {
                "method": method,
                "quantity": key,
                "point_estimate": pe,
                "boot_mean": float(s.mean()),
                "boot_median": float(s.median()),
                "boot_sd": float(s.std(ddof=1)),
                "ci_2.5": float(np.nanpercentile(s, 2.5)),
                "ci_97.5": float(np.nanpercentile(s, 97.5)),
                "prop_gt_0": float(np.mean(s > 0))
                if key in ("G_access", "G_assoc", "Gap_total")
                else np.nan,
                "n_boot": int(len(s)),
                "n_nan": int(s.isna().sum()),
            }
        )
    return pd.DataFrame(rows)


def main():
    scvi.settings.seed = SEED
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    path = ROOT / "data/processed/pbmc_cite_combined_inner_annotated.h5ad"
    log(f"Loading {path}")
    adata = ad.read_h5ad(path)
    high = adata.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    adata_hc = adata[high].copy()
    split = adata_hc.obs["split"].astype(str).to_numpy()
    tr = split == "train"
    te = split == "test"
    assert int(tr.sum()) == 7573, tr.sum()
    assert int(te.sum()) == 1921, te.sum()

    adata_train = ensure_counts_csr(adata_hc[tr])
    adata_test = ensure_counts_csr(adata_hc[te])
    y_tr = adata_train.obs["cell_type_l2"].astype(str).to_numpy()
    y_te = adata_test.obs["cell_type_l2"].astype(str).to_numpy()
    names_tr = adata_train.obs_names.astype(str).to_numpy()
    names_te = adata_test.obs_names.astype(str).to_numpy()

    meta = {
        "n_train": int(adata_train.n_obs),
        "n_test": int(adata_test.n_obs),
        "n_genes": int(adata_train.n_vars),
        "n_proteins": int(adata_train.obsm["protein_expression"].shape[1]),
        "scvi_tools": scvi.__version__,
        "torch": torch.__version__,
        "device": accelerator(),
        "seed": SEED,
        "protocol": "strict_train_only_HC_7573_1921",
    }
    (OUT / "tables" / "run_meta.json").write_text(json.dumps(meta, indent=2))
    log(json.dumps(meta))

    # --- PCA path (train-only unsupervised fit; transform test) ---
    log("Fitting RNA/protein/concat PCA on train only")
    rna = fit_rna_pca(adata_train, adata_test, PRIMARY_RNA_PCS)
    prot = fit_protein_pca(adata_train, adata_test, PRIMARY_PROT_PCS)
    cat = concat_and_scale(
        rna["Z_fit"], prot["Z_fit"], rna["Z_transform"], prot["Z_transform"]
    )
    Z_rna_tr, Z_rna_te = rna["Z_fit"], rna["Z_transform"]
    Z_prot_tr, Z_prot_te = prot["Z_fit"], prot["Z_transform"]
    Z_cat_tr, Z_cat_te = cat["Z_fit"], cat["Z_transform"]
    assert Z_rna_tr.shape[0] == 7573 and Z_rna_te.shape[0] == 1921
    assert Z_cat_tr.shape[0] == 7573 and Z_cat_te.shape[0] == 1921

    pred_rna, _ = fit_predict(Z_rna_tr, y_tr, Z_rna_te)
    pred_prot, _ = fit_predict(Z_prot_tr, y_tr, Z_prot_te)
    pred_cat, _ = fit_predict(Z_cat_tr, y_tr, Z_cat_te)
    f1_rna = macro_f1(y_te, pred_rna)
    f1_prot = macro_f1(y_te, pred_prot)
    f1_cat = macro_f1(y_te, pred_cat)
    log(f"PCA F1 RNA={f1_rna:.6f} protein={f1_prot:.6f} concat={f1_cat:.6f}")

    # --- Deep models ---
    device = accelerator()
    scvi_model, scvi_dir = train_scvi_matched(adata_train, device)
    Z_scvi_tr = np.asarray(scvi_model.get_latent_representation())
    log("Encoding scVI_matched test (frozen weights; no query train)")
    Z_scvi_te = encode_heldout(scvi_model, adata_test)
    assert Z_scvi_tr.shape[0] == 7573 and Z_scvi_te.shape[0] == 1921

    tot_model, tot_dir = train_totalvi(adata_train, device)
    Z_tot_tr = np.asarray(tot_model.get_latent_representation())
    log("Encoding totalVI test (frozen weights; no query train)")
    Z_tot_te = encode_heldout(tot_model, adata_test)
    assert Z_tot_tr.shape[0] == 7573 and Z_tot_te.shape[0] == 1921

    np.save(OUT / "embeddings" / "scvi_matched_train_latent.npy", Z_scvi_tr)
    np.save(OUT / "embeddings" / "scvi_matched_test_latent.npy", Z_scvi_te)
    np.save(OUT / "embeddings" / "totalvi_train_latent.npy", Z_tot_tr)
    np.save(OUT / "embeddings" / "totalvi_test_latent.npy", Z_tot_te)
    pd.DataFrame({"cell_id": names_tr}).to_csv(OUT / "embeddings" / "train_cells.csv", index=False)
    pd.DataFrame({"cell_id": names_te, "y": y_te}).to_csv(
        OUT / "embeddings" / "test_cells.csv", index=False
    )

    pred_scvi, _ = fit_predict(Z_scvi_tr, y_tr, Z_scvi_te)
    pred_tot, _ = fit_predict(Z_tot_tr, y_tr, Z_tot_te)
    f1_scvi = macro_f1(y_te, pred_scvi)
    f1_tot = macro_f1(y_te, pred_tot)
    log(f"Deep F1 scVI_matched={f1_scvi:.6f} totalVI={f1_tot:.6f}")

    g_access = f1_cat - f1_rna
    g_assoc = f1_tot - f1_cat
    gap = f1_tot - f1_rna
    r_access = g_access / gap if abs(gap) > 1e-12 else float("nan")

    point = {
        "F1_RNA_PCA": f1_rna,
        "F1_protein_PCA": f1_prot,
        "F1_concat": f1_cat,
        "F1_scVI_matched": f1_scvi,
        "F1_totalVI": f1_tot,
        "G_access": g_access,
        "G_assoc": g_assoc,
        "Gap_total": gap,
        "R_access": r_access,
    }
    pd.DataFrame([point]).to_csv(OUT / "tables" / "strict_train_only_point_estimates.csv", index=False)
    log(f"Contrasts: {point}")

    # predictions cache for bootstrap (RNA, concat, totalVI primary)
    np.savez_compressed(
        OUT / "tables" / "strict_train_only_heldout_predictions.npz",
        y=y_te,
        pred_rna=pred_rna,
        pred_protein=pred_prot,
        pred_concat=pred_cat,
        pred_scvi_matched=pred_scvi,
        pred_totalvi=pred_tot,
        cell_id=names_te,
    )

    log("Running global paired bootstrap B=10000")
    boot = bootstrap_paired(y_te, pred_rna, pred_cat, pred_tot, BOOT_SEED, stratified=False)
    boot.to_csv(OUT / "tables" / "strict_train_only_bootstrap_draws.csv", index=False)
    boot_point = {
        "F1_RNA_PCA": f1_rna,
        "F1_concat": f1_cat,
        "F1_totalVI": f1_tot,
        "G_access": g_access,
        "G_assoc": g_assoc,
        "Gap_total": gap,
        "R_access": r_access,
    }
    summ = summarize_boot(boot, boot_point, "global_paired_strict_train_only")
    summ.to_csv(OUT / "tables" / "strict_train_only_bootstrap_summary.csv", index=False)

    log("Running class-stratified paired bootstrap B=10000")
    boot_s = bootstrap_paired(y_te, pred_rna, pred_cat, pred_tot, STRAT_SEED, stratified=True)
    summ_s = summarize_boot(boot_s, boot_point, "class_stratified_paired_strict_train_only")
    summ_s.to_csv(OUT / "tables" / "strict_train_only_bootstrap_class_stratified.csv", index=False)
    compare = pd.concat([summ, summ_s], ignore_index=True)
    compare.to_csv(OUT / "tables" / "strict_train_only_bootstrap_global_vs_stratified.csv", index=False)

    log("DONE")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
