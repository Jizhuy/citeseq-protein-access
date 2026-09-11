#!/usr/bin/env python3
"""PART A — Simple RNA+protein concatenated PCA baseline (P1).

Primary concat dimensionality is PREDECLARED: 10 RNA PCs + 10 protein PCs = 20.
All PCA/scaling fit on TRAIN or SOURCE only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from revision2_p1_common import (  # noqa: E402
    CTRL,
    FIGS,
    LOGS,
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    ROOT,
    SECONDARY_RNA_PCS,
    TABLES,
    align_latent,
    concat_and_scale,
    done_ok,
    fit_protein_pca,
    fit_rna_pca,
    logreg_metrics,
    write_done,
)

OUT = CTRL / "simple_multimodal_baseline"
DONE = OUT / "done" / "SIMPLE_MULTIMODAL_DONE.json"
EMB = OUT / "embeddings"


def load_annotated() -> ad.AnnData:
    path = ROOT / "data/processed/pbmc_cite_combined_inner_annotated.h5ad"
    return ad.read_h5ad(path)


def development_benchmark(adata: ad.AnnData, *, secondary: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """High-confidence development held-out benchmark with train-only PCA fits."""
    high = adata.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    split = adata.obs["split"].astype(str).to_numpy()
    train_mask = high & (split == "train")
    test_mask = high & (split == "test")
    adata_train = adata[train_mask].copy()
    adata_all_hc = adata[high].copy()  # transform all HC then index train/test

    # Fit on train HC only; transform all HC cells
    rna = fit_rna_pca(adata_train, adata_all_hc, PRIMARY_RNA_PCS)
    prot = fit_protein_pca(adata_train, adata_all_hc, PRIMARY_PROT_PCS)
    cat = concat_and_scale(rna["Z_fit"], prot["Z_fit"], rna["Z_transform"], prot["Z_transform"])

    # Map back: Z_transform rows correspond to adata_all_hc order
    hc_names = adata_all_hc.obs_names.astype(str)
    split_hc = adata_all_hc.obs["split"].astype(str).to_numpy()
    y_l2 = adata_all_hc.obs["cell_type_l2"].astype(str).to_numpy()
    tr = split_hc == "train"
    te = split_hc == "test"

    # Frozen deep / MOFA embeddings aligned to full adata then subset HC
    emb_dir = ROOT / "results/embeddings"
    latents = {
        "RNA_PCA_10": rna["Z_transform"],
        "protein_PCA_10": prot["Z_transform"],
        "concat_PCA_10_10": cat["Z_transform"],
        "MOFA+": align_latent(emb_dir / "mofa_seed0_cells.csv", emb_dir / "mofa_seed0_latent.npy", hc_names),
        "scVI_matched": align_latent(
            emb_dir / "scvi_matched_seed0_cells.csv", emb_dir / "scvi_matched_seed0_latent.npy", hc_names
        ),
        "totalVI": align_latent(emb_dir / "totalvi_seed0_cells.csv", emb_dir / "totalvi_seed0_latent.npy", hc_names),
    }

    # Secondary higher-capacity concat (optional)
    if secondary:
        n_prot = prot["Z_fit"].shape[1]
        # Use all usable protein PCs (rank-limited) + 20 RNA
        rna20 = fit_rna_pca(adata_train, adata_all_hc, SECONDARY_RNA_PCS)
        n_prot_sec = min(13, protein_n_usable(adata_train))
        prot_sec = fit_protein_pca(adata_train, adata_all_hc, n_prot_sec)
        cat_sec = concat_and_scale(
            rna20["Z_fit"], prot_sec["Z_fit"], rna20["Z_transform"], prot_sec["Z_transform"]
        )
        latents["concat_PCA_20_plus_prot"] = cat_sec["Z_transform"]

    rows = []
    per_class_rows = []
    for name, Z in latents.items():
        m = logreg_metrics(Z[tr], y_l2[tr], Z[te], y_l2[te])
        rows.append(
            {
                "setting": "development_HC",
                "representation": name,
                "label_level": "l2",
                "fit_cells": "train_HC_only",
                "macro_f1": m["macro_f1"],
                "accuracy": m["accuracy"],
                "balanced_accuracy": m["balanced_accuracy"],
                "weighted_f1": m["weighted_f1"],
                "n_train": m["n_train"],
                "n_test": m["n_test"],
                "n_classes_test": m["n_classes_test"],
                "dim": int(Z.shape[1]),
                "is_primary_concat": name == "concat_PCA_10_10",
            }
        )
        for ct, f1 in m["per_class_f1"].items():
            per_class_rows.append(
                {
                    "setting": "development_HC",
                    "representation": name,
                    "cell_type": ct,
                    "f1": f1,
                    "n_test": int((m["y_test"] == ct).sum()),
                }
            )

    # Save primary concat embedding for HC cells
    EMB.mkdir(parents=True, exist_ok=True)
    np.save(EMB / "development_concat_PCA_10_10.npy", cat["Z_transform"])
    pd.DataFrame({"cell_id": hc_names}).to_csv(EMB / "development_HC_cells.csv", index=False)
    meta = {
        "primary": "10_RNA_PC + 10_protein_PC",
        "rna_pcs": rna["n_pcs"],
        "protein_pcs": prot["n_pcs"],
        "concat_dim": cat["dim"],
        "n_hvg": len(rna["hvgs"]),
        "fit": "train high-confidence only; transform all HC",
    }
    (EMB / "development_concat_meta.json").write_text(json.dumps(meta, indent=2))
    return pd.DataFrame(rows), pd.DataFrame(per_class_rows)


def protein_n_usable(adata: ad.AnnData) -> int:
    from revision2_p1_common import protein_matrix, clr_rows

    P, _ = protein_matrix(adata)
    C = clr_rows(P)
    # rank after standardization
    C = (C - C.mean(0)) / (C.std(0, ddof=1) + 1e-8)
    return int(min(C.shape[1], np.linalg.matrix_rank(C)))


def transfer_class_plan(source, target, label_col="cell_type_l2", min_src=20, min_tgt=10):
    """Local copy of PHASE10 class_plan without importing scvi/torch deps."""
    high_s = source.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    high_t = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    ys = source.obs[label_col].astype(str).to_numpy()
    yt = target.obs[label_col].astype(str).to_numpy()
    src_counts = pd.Series(ys[high_s]).value_counts()
    tgt_counts = pd.Series(yt[high_t]).value_counts()
    intersection = sorted(set(src_counts.index) & set(tgt_counts.index))
    eligible = [
        c
        for c in intersection
        if int(src_counts.get(c, 0)) >= min_src and int(tgt_counts.get(c, 0)) >= min_tgt
    ]
    src_keep = high_s & np.isin(ys, eligible)
    tgt_keep = high_t & np.isin(yt, eligible)
    return {
        "eligible_shared_classes": eligible,
        "source_eval_index": np.where(src_keep)[0],
        "target_eval_index": np.where(tgt_keep)[0],
    }


def transfer_benchmark(adata: ad.AnnData) -> pd.DataFrame:
    """Inductive concat PCA on Direction A/B; source-only fit."""
    from sklearn.preprocessing import StandardScaler

    rows = []
    directions = {
        "direction_A": ("PBMC10k", "PBMC5k"),
        "direction_B": ("PBMC5k", "PBMC10k"),
    }
    p11 = ROOT / "results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv"
    frozen = {}
    if p11.exists():
        df = pd.read_csv(p11)
        for (dname, model), sub in df.groupby(["direction", "model"]):
            frozen[(dname, model)] = {
                "macro_f1_mean": float(sub.macro_f1.mean()),
                "macro_f1_sd": float(sub.macro_f1.std(ddof=1)),
                "n_seeds": int(len(sub)),
            }

    for dname, (src_batch, tgt_batch) in directions.items():
        source = adata[adata.obs["batch"].astype(str) == src_batch].copy()
        target = adata[adata.obs["batch"].astype(str) == tgt_batch].copy()
        plan = transfer_class_plan(source, target, "cell_type_l2")
        classes = plan["eligible_shared_classes"]
        src_idx = np.asarray(plan["source_eval_index"])
        tgt_idx = np.asarray(plan["target_eval_index"])

        rna_pack = _fit_rna_once(source, [source, target], PRIMARY_RNA_PCS)
        prot_pack = _fit_prot_once(source, [source, target], PRIMARY_PROT_PCS)
        X_s = np.hstack([rna_pack["Zs"][0], prot_pack["Zs"][0]])
        X_t = np.hstack([rna_pack["Zs"][1], prot_pack["Zs"][1]])
        scaler = StandardScaler().fit(X_s)
        Z_s = scaler.transform(X_s)
        Z_t = scaler.transform(X_t)

        y_s = source.obs["cell_type_l2"].astype(str).to_numpy()[src_idx]
        y_t = target.obs["cell_type_l2"].astype(str).to_numpy()[tgt_idx]
        keep_s = np.isin(y_s, classes)
        keep_t = np.isin(y_t, classes)
        m = logreg_metrics(Z_s[src_idx][keep_s], y_s[keep_s], Z_t[tgt_idx][keep_t], y_t[keep_t])
        rows.append(
            {
                "direction": dname,
                "representation": "concat_PCA_10_10",
                "adaptation": "none_inductive_source_fit_only",
                "macro_f1": m["macro_f1"],
                "accuracy": m["accuracy"],
                "balanced_accuracy": m["balanced_accuracy"],
                "n_train": m["n_train"],
                "n_test": m["n_test"],
                "n_classes": m["n_classes_test"],
                "dim": Z_s.shape[1],
            }
        )
        for rep_name, Zs, Zt in [
            ("RNA_PCA_10", rna_pack["Zs"][0], rna_pack["Zs"][1]),
            ("protein_PCA_10", prot_pack["Zs"][0], prot_pack["Zs"][1]),
        ]:
            mm = logreg_metrics(Zs[src_idx][keep_s], y_s[keep_s], Zt[tgt_idx][keep_t], y_t[keep_t])
            rows.append(
                {
                    "direction": dname,
                    "representation": rep_name,
                    "adaptation": "none_inductive_source_fit_only",
                    "macro_f1": mm["macro_f1"],
                    "accuracy": mm["accuracy"],
                    "balanced_accuracy": mm["balanced_accuracy"],
                    "n_train": mm["n_train"],
                    "n_test": mm["n_test"],
                    "n_classes": mm["n_classes_test"],
                    "dim": Zs.shape[1],
                }
            )
        for model_key, label in [("scvi_matched", "scVI_matched"), ("totalvi", "totalVI")]:
            key = (dname, model_key)
            if key in frozen:
                rows.append(
                    {
                        "direction": dname,
                        "representation": label,
                        "adaptation": "scArches_unsupervised_transductive",
                        "macro_f1": frozen[key]["macro_f1_mean"],
                        "macro_f1_sd": frozen[key]["macro_f1_sd"],
                        "n_seeds": frozen[key]["n_seeds"],
                        "accuracy": np.nan,
                        "balanced_accuracy": np.nan,
                        "n_train": np.nan,
                        "n_test": np.nan,
                        "n_classes": np.nan,
                        "dim": 20,
                    }
                )
    return pd.DataFrame(rows)


def _fit_rna_once(fit_adata, transform_list, n_pcs):
    # Fit on fit_adata; transform each
    first = fit_rna_pca(fit_adata, transform_list[0], n_pcs)
    Zs = [first["Z_transform"]]
    # Reuse hvgs/mu/sd/pca manually for remaining
    from revision2_p1_common import counts_matrix, _to_dense
    import scanpy as sc

    hvgs = first["hvgs"]
    mu, sd, pca = first["rna_mu"], first["rna_sd"], first["pca"]
    for i, a in enumerate(transform_list):
        if i == 0:
            continue
        rna = ad.AnnData(X=counts_matrix(a))
        rna.obs_names = a.obs_names
        rna.var_names = a.var_names
        rna = rna[:, hvgs].copy()
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
        L = _to_dense(rna.X)
        Ls = np.clip((L - mu) / sd, -10, 10)
        Zs.append(pca.transform(Ls))
    return {"Zs": Zs, "hvgs": hvgs, "n_pcs": first["n_pcs"]}


def _fit_prot_once(fit_adata, transform_list, n_pcs):
    first = fit_protein_pca(fit_adata, transform_list[0], n_pcs)
    Zs = [first["Z_transform"]]
    from revision2_p1_common import protein_matrix, clr_rows

    mu, sd, pca = first["prot_mu"], first["prot_sd"], first["pca"]
    for i, a in enumerate(transform_list):
        if i == 0:
            continue
        P, _ = protein_matrix(a)
        C = clr_rows(P)
        Zin = (C - mu) / sd
        Zs.append(pca.transform(Zin))
    return {"Zs": Zs, "n_pcs": first["n_pcs"]}


def lawlor_benchmark() -> pd.DataFrame:
    path = ROOT / "results/phase12b_external_validation/formal_validation/processed/lawlor_baseline_sng_shared.h5ad"
    if not path.exists():
        raise FileNotFoundError(path)
    adata = ad.read_h5ad(path)
    folds = pd.read_csv(
        ROOT / "results/phase12b_external_validation/formal_validation/tables/phase12b_donor_folds.csv"
    )
    # primary classes
    import json as _json

    cls_path = ROOT / "results/phase12b_external_validation/formal_validation/logs/primary_class_set.json"
    if cls_path.exists():
        classes = _json.loads(cls_path.read_text()).get("primary_eligible_classes")
    else:
        classes = None
    if not classes:
        classes = ["B", "CD14_Mono", "CD4T_Mem", "CD4T_Naive", "CD8T_Mem", "CD8T_Naive", "NK"]

    label_col = "author_cell_type" if "author_cell_type" in adata.obs else "cell_type"
    rows = []
    donor_rows = []
    # Map chip->Donor labels already in obs?
    if "donor" not in adata.obs.columns:
        # map via annotations
        ann = pd.read_csv(ROOT / "results/phase12b_external_validation/raw_audit/baseline_sng_annotations.csv")
        counts = pd.read_csv(ROOT / "results/phase12b_external_validation/tables/external_donor_condition_counts.csv")
        chip = counts[["donor", "Donor_of_Origin"]].drop_duplicates().set_index("Donor_of_Origin")["donor"]
        ann_map = ann.set_index("barcode")["Donor_of_Origin"].map(chip)
        adata.obs["donor"] = adata.obs_names.astype(str).map(ann_map.to_dict())

    for fold in range(5):
        tgt_donors = folds.loc[
            (folds.fold == fold) & (folds.role_when_fold_target == "target"), "donor"
        ].astype(str).tolist()
        src_donors = folds.loc[
            (folds.fold == fold) & (folds.role_when_fold_target == "source"), "donor"
        ].astype(str).tolist()
        source = adata[adata.obs["donor"].astype(str).isin(src_donors)].copy()
        target = adata[adata.obs["donor"].astype(str).isin(tgt_donors)].copy()
        rna_pack = _fit_rna_once(source, [source, target], PRIMARY_RNA_PCS)
        prot_pack = _fit_prot_once(source, [source, target], PRIMARY_PROT_PCS)
        from sklearn.preprocessing import StandardScaler

        X_s = np.hstack([rna_pack["Zs"][0], prot_pack["Zs"][0]])
        X_t = np.hstack([rna_pack["Zs"][1], prot_pack["Zs"][1]])
        scaler = StandardScaler().fit(X_s)
        Z_s, Z_t = scaler.transform(X_s), scaler.transform(X_t)

        y_s_all = source.obs[label_col].astype(str).to_numpy()
        y_t_all = target.obs[label_col].astype(str).to_numpy()
        # labeled eligible
        src_lab = np.isin(y_s_all, classes) & (y_s_all != "nan") & (~pd.isna(source.obs[label_col]).to_numpy())
        tgt_lab = np.isin(y_t_all, classes) & (y_t_all != "nan") & (~pd.isna(target.obs[label_col]).to_numpy())
        # Also require non-null label string not 'nan'
        src_lab = src_lab & np.array([c in classes for c in y_s_all])
        tgt_lab = tgt_lab & np.array([c in classes for c in y_t_all])

        m = logreg_metrics(Z_s[src_lab], y_s_all[src_lab], Z_t[tgt_lab], y_t_all[tgt_lab])
        rows.append(
            {
                "fold": fold,
                "target_donors": ";".join(tgt_donors),
                "representation": "concat_PCA_10_10",
                "macro_f1": m["macro_f1"],
                "accuracy": m["accuracy"],
                "balanced_accuracy": m["balanced_accuracy"],
                "n_train": m["n_train"],
                "n_test": m["n_test"],
            }
        )
        # per-donor target metrics
        tgt_donor_arr = target.obs["donor"].astype(str).to_numpy()
        for d in tgt_donors:
            mask = tgt_lab & (tgt_donor_arr == d)
            if mask.sum() < 5:
                continue
            # Need predictions: reuse classifier
            clf = m["clf"]
            pred = clf.predict(Z_t[mask])
            yt = y_t_all[mask]
            donor_rows.append(
                {
                    "fold": fold,
                    "donor": d,
                    "representation": "concat_PCA_10_10",
                    "macro_f1": float(
                        __import__("sklearn.metrics", fromlist=["f1_score"]).f1_score(
                            yt, pred, average="macro", zero_division=0
                        )
                    ),
                    "accuracy": float(accuracy_score_safe(yt, pred)),
                    "n_cells": int(mask.sum()),
                }
            )
        # RNA / protein only fold-level
        for rep, Zs, Zt in [
            ("RNA_PCA_10", rna_pack["Zs"][0], rna_pack["Zs"][1]),
            ("protein_PCA_10", prot_pack["Zs"][0], prot_pack["Zs"][1]),
        ]:
            mm = logreg_metrics(Zs[src_lab], y_s_all[src_lab], Zt[tgt_lab], y_t_all[tgt_lab])
            rows.append(
                {
                    "fold": fold,
                    "target_donors": ";".join(tgt_donors),
                    "representation": rep,
                    "macro_f1": mm["macro_f1"],
                    "accuracy": mm["accuracy"],
                    "balanced_accuracy": mm["balanced_accuracy"],
                    "n_train": mm["n_train"],
                    "n_test": mm["n_test"],
                }
            )

    # Attach Lawlor deep model donor-level means from revision P0 if present
    donor_delta = ROOT / "results/revision_round2_methodological_fixes/statistics/lawlor_donor_seedmean.csv"
    if donor_delta.exists():
        dd = pd.read_csv(donor_delta)
        for model in ("scvi_matched", "totalvi"):
            sub = dd[dd.model == model]
            for _, r in sub.iterrows():
                donor_rows.append(
                    {
                        "fold": np.nan,
                        "donor": r["donor"],
                        "representation": "scVI_matched" if model == "scvi_matched" else "totalVI",
                        "macro_f1": r.get("present_class_macro_f1", np.nan),
                        "accuracy": r.get("accuracy", np.nan),
                        "n_cells": np.nan,
                        "source": "P0_seedmean",
                    }
                )

    pd.DataFrame(donor_rows).to_csv(TABLES / "simple_multimodal_lawlor_per_donor.csv", index=False)
    return pd.DataFrame(rows)


def accuracy_score_safe(y, p):
    from sklearn.metrics import accuracy_score

    return accuracy_score(y, p)


def write_log(dev: pd.DataFrame, transfer: pd.DataFrame, lawlor: pd.DataFrame) -> None:
    def get(df, rep, col="macro_f1"):
        sub = df[df.representation == rep]
        return float(sub[col].iloc[0]) if len(sub) else float("nan")

    scvi = get(dev, "scVI_matched")
    tot = get(dev, "totalVI")
    concat = get(dev, "concat_PCA_10_10")
    md = f"""# Simple multimodal baseline (PART A)

## Predeclared primary representation

**10 RNA PCs + 10 protein PCs = 20-D concatenated representation**
(matched to scVI/totalVI latent dimensionality).

Secondary (optional): 20 RNA PCs + rank-limited protein PCs — not used for primary claims.

## Preprocessing

- RNA: counts → HVG(2000, seurat_v3) → normalize 1e4 → log1p → scale (clip 10) → PCA
- Protein: CLR → train/source z-score → PCA
- Concat: StandardScaler on concatenated PCs
- **All fits on TRAIN (development) or SOURCE (transfer/Lawlor) only**

## Development high-confidence l2 macro-F1 (logreg)

| Representation | macro-F1 |
|---|---:|
| RNA PCA (10) | {get(dev,'RNA_PCA_10'):.4f} |
| protein PCA (10) | {get(dev,'protein_PCA_10'):.4f} |
| **concat PCA 10+10** | **{concat:.4f}** |
| MOFA+ | {get(dev,'MOFA+'):.4f} |
| scVI_matched | {scvi:.4f} |
| totalVI | {tot:.4f} |

Delta(totalVI − scVI) = {tot - scvi:.4f}
Delta(totalVI − concat) = {tot - concat:.4f}
Delta(concat − scVI) = {concat - scvi:.4f}

Fraction of totalVI−scVI gain reproduced by concat:
{((concat - scvi) / (tot - scvi)) if (tot - scvi) != 0 else float('nan'):.3f}

## Transfer (inductive PCA vs transductive VAEs)

See `tables/simple_multimodal_transfer.csv`.

## Lawlor (source-fit PCA per fold)

See `tables/simple_multimodal_lawlor.csv` (deterministic; one value per fold).
"""
    (LOGS / "simple_multimodal_baseline.md").write_text(md)
    (OUT / "logs" / "simple_multimodal_baseline.md").write_text(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-secondary", action="store_true")
    ap.add_argument("--skip-lawlor", action="store_true")
    ap.add_argument("--skip-transfer", action="store_true")
    args = ap.parse_args()

    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    if args.status:
        print("DONE" if done_ok(DONE) else "PENDING", DONE)
        return

    if args.resume and done_ok(DONE):
        print("Already complete:", DONE)
        return

    if args.dry_run:
        print("dry-run OK: would run development/transfer/Lawlor concat PCA baselines")
        return

    print("Loading annotated PBMC…")
    adata = load_annotated()
    print("Development benchmark…")
    dev, per_class = development_benchmark(adata, secondary=not args.skip_secondary)
    # GATE 1 print
    print("\n=== GATE 1: development macro-F1 ===")
    print(dev[["representation", "macro_f1", "dim"]].to_string(index=False))
    dev.to_csv(TABLES / "simple_multimodal_development.csv", index=False)
    per_class.to_csv(TABLES / "simple_multimodal_per_class.csv", index=False)

    if args.skip_transfer:
        transfer = pd.DataFrame()
    else:
        print("Transfer benchmark…")
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        transfer = transfer_benchmark(adata)
        transfer.to_csv(TABLES / "simple_multimodal_transfer.csv", index=False)
        print(transfer[["direction", "representation", "macro_f1"]].to_string(index=False))

    if args.skip_lawlor:
        lawlor = pd.DataFrame()
    else:
        print("Lawlor benchmark…")
        lawlor = lawlor_benchmark()
        lawlor.to_csv(TABLES / "simple_multimodal_lawlor.csv", index=False)
        print(lawlor.groupby("representation")["macro_f1"].mean())

    write_log(dev, transfer, lawlor)
    write_done(
        DONE,
        {
            "ok": True,
            "outputs": [
                str(TABLES / "simple_multimodal_development.csv"),
                str(TABLES / "simple_multimodal_transfer.csv"),
                str(TABLES / "simple_multimodal_lawlor.csv"),
                str(TABLES / "simple_multimodal_per_class.csv"),
            ],
            "primary_concat": "10+10",
        },
    )
    print("Wrote DONE", DONE)


if __name__ == "__main__":
    main()
