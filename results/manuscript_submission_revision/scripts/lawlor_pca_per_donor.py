#!/usr/bin/env python3
"""Extend Lawlor PCA baselines to per-donor present-class macro-F1 (n=10).

Uses the same source-donor-only PCA fits as the Lawlor simple-multimodal control;
evaluates each target donor separately. Does NOT retrain VAEs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "results/revision_round2_methodological_fixes/scripts"))
from revision2_p1_common import (  # noqa: E402
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    SEED,
    fit_protein_pca,
    fit_rna_pca,
    logreg_metrics,
)

OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)

CLASSES = ["B", "CD14_Mono", "CD4T_Mem", "CD4T_Naive", "CD8T_Mem", "CD8T_Naive", "NK"]


def main():
    path = ROOT / "results/phase12b_external_validation/formal_validation/processed/lawlor_baseline_sng_shared.h5ad"
    adata = ad.read_h5ad(path)
    folds = pd.read_csv(
        ROOT / "results/phase12b_external_validation/formal_validation/tables/phase12b_donor_folds.csv"
    )
    label_col = "author_cell_type" if "author_cell_type" in adata.obs else "cell_type"
    if "donor" not in adata.obs.columns:
        ann = pd.read_csv(ROOT / "results/phase12b_external_validation/raw_audit/baseline_sng_annotations.csv")
        counts = pd.read_csv(ROOT / "results/phase12b_external_validation/tables/external_donor_condition_counts.csv")
        chip = counts[["donor", "Donor_of_Origin"]].drop_duplicates().set_index("Donor_of_Origin")["donor"]
        ann_map = ann.set_index("barcode")["Donor_of_Origin"].map(chip)
        adata.obs["donor"] = adata.obs_names.astype(str).map(ann_map.to_dict())

    rows = []
    for fold in range(5):
        tgt_donors = folds.loc[
            (folds.fold == fold) & (folds.role_when_fold_target == "target"), "donor"
        ].astype(str).tolist()
        src_donors = folds.loc[
            (folds.fold == fold) & (folds.role_when_fold_target == "source"), "donor"
        ].astype(str).tolist()
        source = adata[adata.obs["donor"].astype(str).isin(src_donors)].copy()
        target_all = adata[adata.obs["donor"].astype(str).isin(tgt_donors)].copy()

        rna = fit_rna_pca(source, target_all, PRIMARY_RNA_PCS)
        rna_s = fit_rna_pca(source, source, PRIMARY_RNA_PCS)
        prot = fit_protein_pca(source, target_all, PRIMARY_PROT_PCS)
        prot_s = fit_protein_pca(source, source, PRIMARY_PROT_PCS)

        Zs_rna = rna_s["Z_fit"]
        Zt_rna = rna["Z_transform"]
        Zs_prot = prot_s["Z_fit"]
        Zt_prot = prot["Z_transform"]

        y_s = source.obs[label_col].astype(str).to_numpy()
        src_lab = np.isin(y_s, CLASSES)
        donors_t = target_all.obs["donor"].astype(str).to_numpy()
        y_t = target_all.obs[label_col].astype(str).to_numpy()

        for name, Zs, Zt in [
            ("RNA_PCA_10", Zs_rna, Zt_rna),
            ("protein_PCA_10", Zs_prot, Zt_prot),
            (
                "concat_PCA_10_10",
                np.hstack([Zs_rna, Zs_prot]),
                np.hstack([Zt_rna, Zt_prot]),
            ),
        ]:
            sc = StandardScaler().fit(Zs)
            Zs_sc = sc.transform(Zs)
            Zt_sc = sc.transform(Zt)
            for d in tgt_donors:
                mask_d = donors_t == d
                tgt_lab = np.isin(y_t, CLASSES) & mask_d
                if tgt_lab.sum() < 5:
                    continue
                m = logreg_metrics(
                    Zs_sc[src_lab], y_s[src_lab], Zt_sc[tgt_lab], y_t[tgt_lab], seed=SEED
                )
                present = sorted(set(y_t[tgt_lab].tolist()))
                rows.append(
                    {
                        "fold": fold,
                        "donor": d,
                        "representation": name,
                        "macro_f1": m["macro_f1"],
                        "accuracy": m["accuracy"],
                        "n_cells": int(tgt_lab.sum()),
                        "n_present_classes": len(present),
                        "present_classes": ";".join(present),
                        "metric_note": "sklearn macro-F1 on present labeled classes in donor",
                        "source": "pca_per_donor_extension",
                    }
                )

    pca_df = pd.DataFrame(rows)
    pca_df.to_csv(OUT / "lawlor_pca_per_donor.csv", index=False)

    # Prefer previously audited concat donor values if present (sanity).
    existing = ROOT / (
        "results/revision_round2_methodological_fixes/controls/tables/"
        "simple_multimodal_lawlor_per_donor.csv"
    )
    if existing.exists():
        old = pd.read_csv(existing)
        old_concat = old[old["representation"] == "concat_PCA_10_10"][
            ["donor", "representation", "macro_f1"]
        ].copy()
        new_concat = pca_df[pca_df["representation"] == "concat_PCA_10_10"][
            ["donor", "representation", "macro_f1"]
        ].copy()
        merged = old_concat.merge(new_concat, on=["donor", "representation"], suffixes=("_old", "_new"))
        max_abs = float((merged["macro_f1_old"] - merged["macro_f1_new"]).abs().max())
        print(f"concat per-donor max abs vs prior artifact: {max_abs:.6g}")

    vae = pd.read_csv(
        ROOT / "results/revision_round2_methodological_fixes/statistics/lawlor_donor_seedmean.csv"
    )
    vae = vae[["donor", "model", "present_class_macro_f1"]].rename(
        columns={"model": "representation", "present_class_macro_f1": "macro_f1"}
    )
    vae["representation"] = vae["representation"].map(
        {"scvi_matched": "scVI_matched", "totalvi": "totalVI"}
    )

    combined = pd.concat(
        [
            pca_df[["donor", "representation", "macro_f1"]],
            vae[["donor", "representation", "macro_f1"]],
        ],
        ignore_index=True,
    )
    summary = (
        combined.groupby("representation")["macro_f1"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )
    iqr = combined.groupby("representation")["macro_f1"].quantile([0.25, 0.75]).unstack()
    summary["q25"] = summary["representation"].map(iqr[0.25])
    summary["q75"] = summary["representation"].map(iqr[0.75])
    summary["iqr"] = summary["q75"] - summary["q25"]
    summary.to_csv(OUT / "lawlor_all_methods_donor_level_summary.csv", index=False)
    combined.to_csv(OUT / "lawlor_all_methods_per_donor.csv", index=False)
    print(summary.to_string(index=False))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
