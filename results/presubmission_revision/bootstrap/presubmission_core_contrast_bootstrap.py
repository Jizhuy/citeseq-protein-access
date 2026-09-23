#!/usr/bin/env python3
"""Presubmission paired cell-level bootstrap for core protein-access contrasts.

Reconstructs held-out predictions from frozen development embeddings / train-only
PCA fits matching revision2 simple multimodal development benchmark.
Does NOT retrain VAEs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler

# Repo root: results/presubmission_revision/bootstrap/ -> ../../../
PROJ = Path(__file__).resolve().parents[3]
P1 = PROJ / "results/revision_round2_methodological_fixes/scripts"
sys.path.insert(0, str(P1))

from revision2_p1_common import (  # noqa: E402
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    SEED,
    align_latent,
    concat_and_scale,
    fit_protein_pca,
    fit_rna_pca,
)

OUT = Path(__file__).resolve().parent
N_BOOT = 10_000
BOOT_SEED = 20260923
NAVY = "#1B3A5F"
BLUE = "#2F6FAE"


def macro_f1(y_true, y_pred) -> float:
    return float(f1_score(y_true, y_pred, average="macro", zero_division=0))


def fit_predict(Z_tr, y_tr, Z_te, seed: int = SEED):
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    clf.fit(Z_tr, y_tr)
    return clf.predict(Z_te)


def load_development_predictions():
    adata = ad.read_h5ad(PROJ / "data/processed/pbmc_cite_combined_inner_annotated.h5ad")
    high = adata.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    # Match development_benchmark: healthy-control high-confidence if column exists
    if "condition" in adata.obs.columns:
        hc_mask = high & (adata.obs["condition"].astype(str).to_numpy() == "HC")
        # Some datasets use healthy / control naming
        if hc_mask.sum() < 1000:
            cond = adata.obs["condition"].astype(str).str.lower()
            hc_mask = high & cond.isin(["hc", "healthy", "control", "healthy_control"])
    else:
        hc_mask = high
    # Prefer exact n=9494 if available via prior filter
    adata_hc = adata[hc_mask].copy()
    if "split" not in adata_hc.obs.columns:
        raise RuntimeError("split column missing; cannot reproduce held-out cells")

    split = adata_hc.obs["split"].astype(str).to_numpy()
    tr = split == "train"
    te = split == "test"
    y = adata_hc.obs["cell_type_l2"].astype(str).to_numpy()
    names = adata_hc.obs_names

    adata_train = adata_hc[tr].copy()
    # RNA / protein / concat train-only fits (same as revision2)
    rna = fit_rna_pca(adata_train, adata_hc, PRIMARY_RNA_PCS)
    prot = fit_protein_pca(adata_train, adata_hc, PRIMARY_PROT_PCS)
    cat = concat_and_scale(rna["Z_fit"], prot["Z_fit"], rna["Z_transform"], prot["Z_transform"])

    emb = PROJ / "results/embeddings"
    Z_total = align_latent(emb / "totalvi_seed0_cells.csv", emb / "totalvi_seed0_latent.npy", names)

    # Verify frozen concat if present
    frozen_cat = PROJ / "results/revision_round2_methodological_fixes/controls/simple_multimodal_baseline/embeddings/development_concat_PCA_10_10.npy"
    frozen_cells = PROJ / "results/revision_round2_methodological_fixes/controls/simple_multimodal_baseline/embeddings/development_HC_cells.csv"
    if frozen_cat.exists() and frozen_cells.exists():
        Z_cat_frozen = np.load(frozen_cat)
        cells = pd.read_csv(frozen_cells)["cell_id"].astype(str).tolist()
        # Align frozen to current HC order
        idx = pd.Index(cells)
        if list(names.astype(str)) == cells:
            # use frozen concat
            Z_concat = Z_cat_frozen
        else:
            # map
            pos = {c: i for i, c in enumerate(cells)}
            order = [pos[c] for c in names.astype(str)]
            Z_concat = Z_cat_frozen[order]
            # compare to recomputed
            if not np.allclose(Z_concat, cat["Z_transform"], atol=1e-4, rtol=1e-4):
                # prefer recomputed matching current code path; log mismatch
                Z_concat = cat["Z_transform"]
                concat_source = "recomputed_train_only_fit"
            else:
                concat_source = "frozen_npy_aligned"
        if list(names.astype(str)) == cells and np.allclose(Z_cat_frozen, cat["Z_transform"], atol=1e-4, rtol=1e-4):
            Z_concat = Z_cat_frozen
            concat_source = "frozen_npy"
        elif list(names.astype(str)) == cells:
            Z_concat = Z_cat_frozen
            concat_source = "frozen_npy_prefer_over_recompute"
        else:
            Z_concat = cat["Z_transform"]
            concat_source = "recomputed_order_mismatch"
    else:
        Z_concat = cat["Z_transform"]
        concat_source = "recomputed_no_frozen"

    Z_rna = rna["Z_transform"]
    pred_rna = fit_predict(Z_rna[tr], y[tr], Z_rna[te])
    pred_cat = fit_predict(Z_concat[tr], y[tr], Z_concat[te])
    pred_tot = fit_predict(Z_total[tr], y[tr], Z_total[te])
    y_te = y[te]

    f1_rna = macro_f1(y_te, pred_rna)
    f1_cat = macro_f1(y_te, pred_cat)
    f1_tot = macro_f1(y_te, pred_tot)

    audit = {
        "n_hc": int(adata_hc.n_obs),
        "n_train": int(tr.sum()),
        "n_test": int(te.sum()),
        "n_classes_test": int(len(np.unique(y_te))),
        "f1_rna": f1_rna,
        "f1_concat": f1_cat,
        "f1_totalVI": f1_tot,
        "target_rna": 0.666,
        "target_concat": 0.745,
        "target_totalVI": 0.779,
        "match_rna": abs(f1_rna - 0.666) < 0.002,
        "match_concat": abs(f1_cat - 0.745) < 0.002,
        "match_totalVI": abs(f1_tot - 0.779) < 0.002,
        "concat_source": concat_source,
        "classifier": "LogisticRegression max_iter=2000 lbfgs",
        "seed": SEED,
        "note": "Predictions reconstructed from frozen/train-only latents; VAEs not retrained.",
    }
    return y_te, pred_rna, pred_cat, pred_tot, audit


def run_bootstrap(y, p_rna, p_cat, p_tot, n_boot=N_BOOT, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    n = len(y)
    rows = []
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt = y[idx]
        fr = macro_f1(yt, p_rna[idx])
        fc = macro_f1(yt, p_cat[idx])
        ft = macro_f1(yt, p_tot[idx])
        g_access = fc - fr
        g_assoc = ft - fc
        gap = ft - fr
        r_access = g_access / gap if gap != 0 else np.nan
        rows.append(
            {
                "replicate": b,
                "F1_RNA_PCA": fr,
                "F1_concat": fc,
                "F1_totalVI": ft,
                "G_access": g_access,
                "G_assoc": g_assoc,
                "Gap_total": gap,
                "R_access": r_access,
            }
        )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, point: dict) -> pd.DataFrame:
    out = []
    for key, pe in point.items():
        s = df[key]
        out.append(
            {
                "quantity": key,
                "point_estimate": pe,
                "boot_mean": float(s.mean()),
                "boot_median": float(s.median()),
                "boot_sd": float(s.std(ddof=1)),
                "ci_2.5": float(np.nanpercentile(s, 2.5)),
                "ci_97.5": float(np.nanpercentile(s, 97.5)),
                "prop_gt_0": float(np.mean(s > 0)) if key in ("G_access", "G_assoc", "Gap_total") else np.nan,
                "n_boot": int(len(s)),
                "n_nan": int(s.isna().sum()),
            }
        )
    # denominator stability for R_access
    gap = df["Gap_total"]
    out.append(
        {
            "quantity": "Gap_total_min",
            "point_estimate": float(gap.min()),
            "boot_mean": float(gap.mean()),
            "boot_median": float(gap.median()),
            "boot_sd": float(gap.std(ddof=1)),
            "ci_2.5": float(np.nanpercentile(gap, 2.5)),
            "ci_97.5": float(np.nanpercentile(gap, 97.5)),
            "prop_gt_0": float(np.mean(gap > 0)),
            "n_boot": int(len(gap)),
            "n_nan": int((gap <= 0).sum()),
        }
    )
    return pd.DataFrame(out)


def plot_summary(df, summary, point):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6))
    for ax, key, title, color in [
        (axes[0], "G_access", "Protein-access gain", BLUE),
        (axes[1], "G_assoc", "Residual model-associated gain", NAVY),
        (axes[2], "R_access", "Descriptive recovery ratio", BLUE),
    ]:
        ax.hist(df[key].dropna(), bins=40, color=color, alpha=0.85, edgecolor="white")
        pe = point[key]
        lo = float(summary.loc[summary.quantity == key, "ci_2.5"].iloc[0])
        hi = float(summary.loc[summary.quantity == key, "ci_97.5"].iloc[0])
        ax.axvline(pe, color="black", lw=1.5, label=f"point {pe:.3f}")
        ax.axvline(lo, color="gray", ls="--", lw=1)
        ax.axvline(hi, color="gray", ls="--", lw=1, label=f"95% [{lo:.3f},{hi:.3f}]")
        ax.set_title(title, fontweight="bold", color=NAVY)
        ax.set_xlabel(key)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(
        "Cell-level conditional paired-bootstrap (n=10,000) — not biological CI",
        color=NAVY,
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(OUT / "core_contrast_bootstrap.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "core_contrast_bootstrap.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    y, p_rna, p_cat, p_tot, audit = load_development_predictions()
    # Prefer frozen manuscript point estimates for reporting contrasts
    # but bootstrap uses reconstructed predictions; require close match
    pe = {
        "F1_RNA_PCA": audit["f1_rna"],
        "F1_concat": audit["f1_concat"],
        "F1_totalVI": audit["f1_totalVI"],
        "G_access": audit["f1_concat"] - audit["f1_rna"],
        "G_assoc": audit["f1_totalVI"] - audit["f1_concat"],
        "Gap_total": audit["f1_totalVI"] - audit["f1_rna"],
        "R_access": (audit["f1_concat"] - audit["f1_rna"]) / (audit["f1_totalVI"] - audit["f1_rna"]),
    }
    # Manuscript rounded points
    pe_ms = {
        "F1_RNA_PCA": 0.666,
        "F1_concat": 0.745,
        "F1_totalVI": 0.779,
        "G_access": 0.079,
        "G_assoc": 0.034,
        "Gap_total": 0.113,
        "R_access": 0.079 / 0.113,
    }

    df = run_bootstrap(y, p_rna, p_cat, p_tot)
    df.to_csv(OUT / "core_contrast_bootstrap_all.csv", index=False)
    summary = summarize(df, pe)
    # also attach manuscript-rounded points
    summary["manuscript_rounded_point"] = summary["quantity"].map(
        {**pe_ms, "Gap_total_min": np.nan}
    )
    summary.to_csv(OUT / "core_contrast_bootstrap_summary.csv", index=False)
    plot_summary(df, summary, pe)

    gap = df["Gap_total"]
    audit_md = [
        "# Core contrast bootstrap audit",
        "",
        f"- n_boot: {N_BOOT}",
        f"- bootstrap_seed: {BOOT_SEED}",
        f"- classifier_seed: {SEED}",
        f"- n_test_cells: {audit['n_test']}",
        f"- n_train_cells: {audit['n_train']}",
        f"- concat_source: {audit['concat_source']}",
        "",
        "## Reconstructed point estimates vs manuscript rounding",
        "",
        f"- RNA PCA F1: {audit['f1_rna']:.6f} (target 0.666; match={audit['match_rna']})",
        f"- concat F1: {audit['f1_concat']:.6f} (target 0.745; match={audit['match_concat']})",
        f"- totalVI F1: {audit['f1_totalVI']:.6f} (target 0.779; match={audit['match_totalVI']})",
        f"- G_access: {pe['G_access']:.6f}",
        f"- G_assoc: {pe['G_assoc']:.6f}",
        f"- Gap_total: {pe['Gap_total']:.6f}",
        f"- R_access: {pe['R_access']:.6f}",
        "",
        "## Ratio denominator stability",
        "",
        f"- min Gap_total^(b): {float(gap.min()):.6e}",
        f"- n Gap_total <= 0: {int((gap <= 0).sum())}",
        f"- n Gap_total < 0.01: {int((gap < 0.01).sum())}",
        f"- n Gap_total < 0.05: {int((gap < 0.05).sum())}",
        "",
        "## Interval interpretation",
        "",
        "Cell-level conditional paired-bootstrap 95% intervals.",
        "NOT biological-generalization confidence intervals.",
        "",
        "## Summary table",
        "",
        summary.to_string(index=False),
        "",
        "## Full audit JSON",
        "",
        "```json",
        json.dumps(audit, indent=2),
        "```",
    ]
    (OUT / "core_contrast_bootstrap_audit.md").write_text("\n".join(audit_md) + "\n")
    print(json.dumps({"audit": audit, "point": pe, "summary_path": str(OUT / "core_contrast_bootstrap_summary.csv")}, indent=2))


if __name__ == "__main__":
    main()
