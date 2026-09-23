#!/usr/bin/env python3
"""Figures for information-vs-model reframed manuscript (frozen numbers + matched transfer)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
REV2 = ROOT / "results/revision_round2_methodological_fixes"
REV3 = ROOT / "results/revision_round3_information_vs_model"
OUT = ROOT / "results/final_information_vs_model_manuscript/figures"
OUT.mkdir(parents=True, exist_ok=True)

NAVY, BLUE, LBLUE, GRAY, LGRAY, WHITE = "#1B3A5F", "#2F6FAE", "#9DC3E6", "#5A5A5A", "#D9DEE7", "#FFFFFF"
TEAL = "#3D7C8A"


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": WHITE,
            "axes.facecolor": WHITE,
        }
    )


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)
    print("wrote", name)


def fig1():
    style()
    fig, ax = plt.subplots(figsize=(11.2, 6.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("Figure 1. Central logic: information gain vs model-specific advantage", loc="left", fontweight="bold", color=NAVY)

    def box(x, y, w, h, text, fc=LBLUE, tc=NAVY):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", lw=1.2, fc=fc, ec=NAVY))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.5, color=tc, fontweight="bold" if fc == NAVY else "normal")

    box(0.4, 5.8, 2.8, 1.5, "RNA only\nRNA PCA\nF1 = 0.666", LGRAY)
    box(4.0, 5.8, 3.2, 1.5, "Simple multimodal\nRNA+protein concat PCA\nF1 = 0.745", LBLUE)
    box(8.2, 5.8, 3.4, 1.5, "Complex multimodal model\ntotalVI\nF1 = 0.779", NAVY, tc=WHITE)
    ax.annotate("", xy=(4.0, 6.55), xytext=(3.2, 6.55), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
    ax.annotate("", xy=(8.2, 6.55), xytext=(7.2, 6.55), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
    ax.text(3.6, 7.55, "G_info = +0.079\n(~70% of RNA→totalVI gap)", ha="center", fontsize=8, color=BLUE)
    ax.text(7.7, 7.55, "G_model = +0.034\n(residual / context-dependent)", ha="center", fontsize=8, color=NAVY)

    box(0.4, 3.2, 11.2, 1.8,
        "Stress-test the residual model-specific advantage under:\n"
        "dataset transfer (inductive + target-data-access-matched)  ·  protein degradation  ·  training stochasticity  ·  uncertainty diagnostics\n"
        "External donor evidence (Lawlor): multimodal information value across 10 donors — not universal architecture superiority",
        "#EEF3F9")
    box(0.4, 0.5, 11.2, 2.0,
        "Descriptive contrasts only (not causal variance components):\n"
        "G_info = F1_concat − F1_RNA_PCA    ·    G_model = F1_totalVI − F1_concat\n"
        "Central question: how much multimodal gain is protein access, and when does a complex model add reliable value beyond simple multimodal features?",
        WHITE)
    save(fig, "Figure1_central_logic")


def fig2():
    style()
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), gridspec_kw={"width_ratios": [1.4, 0.9, 0.85]})
    # A primary trio emphasized
    ax = axes[0]
    primary = [("RNA PCA", 0.666, LGRAY), ("concat PCA", 0.745, BLUE), ("totalVI", 0.779, NAVY)]
    secondary = [("protein PCA", 0.512), ("MOFA+", 0.549), ("scVI", 0.722)]
    x = np.arange(3)
    ax.bar(x, [v for _, v, _ in primary], color=[c for *_, c in primary], edgecolor=NAVY, width=0.7)
    for i, (lab, v, _) in enumerate(primary):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=8, color=NAVY, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([p[0] for p in primary])
    ax.set_ylim(0.55, 0.85)
    ax.set_ylabel("macro-F1")
    ax.set_title("A. Primary information→model sequence", fontweight="bold")
    ax.text(0.02, 0.08, "Secondary: protein 0.512 · MOFA+ 0.549 · scVI 0.722", transform=ax.transAxes, fontsize=7, color=GRAY)

    ax = axes[1]
    ax.bar([0, 1], [0.079, 0.034], color=[BLUE, NAVY], edgecolor=NAVY, width=0.65)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["G_info\nconcat−RNA", "G_model\ntotalVI−concat"])
    ax.set_ylabel("Δ macro-F1")
    ax.set_title("B. Descriptive increments", fontweight="bold")
    ax.text(0, 0.082, "+0.079", ha="center", fontsize=8)
    ax.text(1, 0.037, "+0.034", ha="center", fontsize=8)
    ax.text(0.5, -0.22, "≈70% of RNA→totalVI gap\nrecovered by concat", ha="center", transform=ax.transAxes, fontsize=7.5, color=BLUE)

    ax = axes[2]
    ax.axis("off")
    ax.set_title("C. RNA-only label audit", fontweight="bold")
    ax.text(0.05, 0.85,
            "Development labels already\nRNA-only (SCANVI/scArches).\n\n"
            "n = 9,494 high-confidence cells\nagreement = 1.0\n\n"
            "Rules out simple label-mismatch\nas the explanation of the\ninternal multimodal gain.\n"
            "Not independent biological\nvalidation of labels.",
            va="top", fontsize=8.5, color=NAVY,
            bbox=dict(boxstyle="round,pad=0.4", fc="#EEF3F9", ec=BLUE))
    fig.tight_layout()
    save(fig, "Figure2_information_vs_model")


def fig3():
    style()
    mt = pd.read_csv(REV3 / "transfer_control/tables/matched_transfer_control.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2))

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A. Different target-data access", fontweight="bold")
    ax.text(
        0.05,
        0.88,
        "SOURCE-ONLY inductive concat\n  · fit on SOURCE only\n  · transform TARGET\n\n"
        "scVI / totalVI (scArches)\n  · unlabeled TARGET features\n    used in query adaptation\n\n"
        "Target-data-access-matched\ntransductive concat (new)\n  · joint unsupervised fit on\n    SOURCE + unlabeled TARGET\n  · FEATURE ACCESS matched only\n  · NOT algorithmically identical",
        va="top",
        fontsize=8.2,
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE),
    )

    def grouped(ax, direction, title, case_note):
        order = [
            ("RNA_PCA_10", "none_inductive_source_fit_only", "RNA\nind.", LGRAY),
            ("concat_PCA_10_10", "none_inductive_source_fit_only", "concat\nind.", LBLUE),
            ("concat_PCA_10_10", "transductive_joint_unsupervised_feature_fit", "concat\ntrans.", BLUE),
            ("scVI_matched", "scArches_unsupervised_transductive", "scVI\nscArches", LBLUE),
            ("totalVI", "scArches_unsupervised_transductive", "totalVI\nscArches", NAVY),
        ]
        vals, labs, cols, errs = [], [], [], []
        for rep, adapt, lab, col in order:
            r = mt[(mt.direction == direction) & (mt.representation == rep) & (mt.adaptation == adapt)].iloc[0]
            vals.append(float(r.macro_f1))
            labs.append(lab)
            cols.append(col)
            errs.append(float(r.macro_f1_sd) if pd.notna(r.macro_f1_sd) else 0.0)
        x = np.arange(len(vals))
        ax.bar(x, vals, yerr=errs, color=cols, edgecolor=NAVY, capsize=3, width=0.72)
        ax.set_xticks(x)
        ax.set_xticklabels(labs, fontsize=7.5)
        ax.set_ylabel("macro-F1")
        ax.set_title(title, fontweight="bold")
        # Full range avoids exaggerating small Direction A gaps.
        ax.set_ylim(0.0, 0.90)
        for i, v in enumerate(vals):
            ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=7, color=NAVY)
        ax.text(
            0.98,
            0.05,
            case_note,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.5,
            color=NAVY,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=BLUE, alpha=0.95),
        )

    grouped(
        axes[0, 1],
        "direction_A",
        "B. Direction A",
        "Case C: totalVI − concat_trans ≈ +0.016",
    )
    grouped(
        axes[1, 0],
        "direction_B",
        "C. Direction B",
        "Case C: totalVI − concat_trans ≈ −0.113",
    )

    ax = axes[1, 1]
    labs = ["Dir A\ntotalVI−scVI", "Dir A\ntotalVI−concat_tr", "Dir B\ntotalVI−scVI", "Dir B\ntotalVI−concat_tr"]
    vals = [-0.0032, 0.0161, 0.0525, -0.1131]
    cols = [LGRAY, BLUE, NAVY, "#C0392B"]
    bars = ax.bar(np.arange(4), vals, color=cols, edgecolor=NAVY)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("Δ macro-F1")
    ax.set_title("D. Case C contrasts (direction-dependent)", fontweight="bold")
    ax.set_ylim(-0.14, 0.08)
    for i, v in enumerate(vals):
        y = v + (0.008 if v >= 0 else -0.012)
        ax.text(i, y, f"{v:+.3f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8, color=NAVY, fontweight="bold")
    ax.text(
        0.5,
        -0.22,
        "Matched target-feature access does not yield a stable totalVI advantage",
        transform=ax.transAxes,
        ha="center",
        fontsize=8,
        color=NAVY,
    )
    fig.tight_layout()
    save(fig, "Figure3_transfer_fairness")


def fig4():
    style()
    sm = pd.read_csv(REV2 / "controls/test_time_corruption/tables/test_time_corruption_summary.csv")
    train = pd.read_csv(ROOT / "results/tables/protein_sparsity_primary_delta_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0))

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A. Training-time degradation", fontweight="bold")
    ax.text(0.05, 0.7, "Corrupt protein → retrain totalVI\n(scVI fixed RNA-only reference)\nWithin-dataset HC development",
            fontsize=9, color=NAVY, va="center", bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    ax = axes[0, 1]
    ax.axis("off")
    ax.set_title("B. Post-adaptation target-only", fontweight="bold")
    ax.text(0.05, 0.7, "Clean train/adapt → freeze\n→ corrupt TARGET protein only\nPrimary: ΔF1 from clean",
            fontsize=9, color=NAVY, va="center", bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    ax = axes[1, 0]
    # training-time F1
    if "corruption_fraction" in train.columns:
        fr = train["corruption_fraction"] if "corruption_fraction" in train.columns else train.iloc[:, 0]
        # try common column names
        f1col = [c for c in train.columns if "macro" in c.lower() or "f1" in c.lower()]
        f1c = f1col[0] if f1col else train.columns[1]
        ax.plot(train.iloc[:, 0] if "corruption" not in train.columns[0].lower() else train["corruption_fraction"],
                train[f1c] if f1c in train.columns else train.iloc[:, 1], "o-", color=NAVY, label="totalVI")
    # use known values if needed
    frs = [0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85]
    tr_f1 = [0.779, 0.771, 0.758, 0.760, 0.750, 0.736, 0.755]
    ax.clear()
    ax.plot(frs, tr_f1, "o-", color=NAVY, label="totalVI (retrain)")
    ax.axhline(0.722, color=LBLUE, ls="--", label="scVI fixed")
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("macro-F1")
    ax.set_title("C. Training-time F1", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1, 1]
    for d, c in [("direction_A", NAVY), ("direction_B", BLUE)]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.macro_f1_mean, "o-", color=c, label=d.replace("direction_", "Dir "))
    ax.axhline(0.601, color=LBLUE, ls=":", label="Dir B scVI clean")
    ax.axhline(0.656, color=LGRAY, ls=":", label="Dir A scVI clean")
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("macro-F1")
    ax.set_title("D. Post-adaptation: advantage narrows", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "Figure4_protein_degradation")


def fig5():
    style()
    var = pd.read_csv(REV2 / "statistics/variance_components_corrected.csv")
    dvar = pd.read_csv(REV2 / "statistics/delta_variance_components_corrected.csv")
    law = pd.read_csv(REV2 / "calibration/lawlor_reliability_summary_by_model.csv")
    pb = pd.read_csv(REV2 / "calibration/pbmc_20seed_reliability_summary.csv")
    sm = pd.read_csv(REV2 / "controls/test_time_corruption/tables/test_time_corruption_summary.csv")

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 6.8))

    # A predictive vs latent AUROC transfer
    ax = axes[0, 0]
    x = np.arange(2)
    for off, m, col, lab in [(-0.18, "scvi_matched", LBLUE, "scVI"), (0.18, "totalvi", NAVY, "totalVI")]:
        pred, lat = [], []
        for d in ["PBMC_direction_A", "PBMC_direction_B"]:
            r = pb[(pb.experiment == d) & (pb.model == m)].iloc[0]
            pred.append(r.error_auroc_mean)
            lat.append(r.latent_error_auroc_mean)
        ax.plot([0, 1], pred, "o-", color=col, label=f"{lab} pred")
        ax.plot([0, 1], lat, "s--", color=col, alpha=0.7, label=f"{lab} latent")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_ylabel("error AUROC")
    ax.set_title("A. Predictive vs latent AUROC", fontweight="bold")
    ax.legend(frameon=False, fontsize=6, ncol=2)

    # B E-AURC under corruption
    ax = axes[0, 1]
    for d, c in [("direction_A", NAVY), ("direction_B", BLUE)]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.e_aurc_mean, "o-", color=c, label=d.replace("direction_", "Dir "))
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("E-AURC")
    ax.set_title("B. Post-adaptation E-AURC", fontweight="bold")
    ax.legend(frameon=False)

    # C paired delta variance
    ax = axes[1, 0]
    r = dvar[(dvar.model == "delta_totalvi_minus_scvi") & (dvar.metric == "macro_f1")].iloc[0]
    comps = ["source_subset_fraction", "model_seed_fraction", "interaction_fraction", "residual_fraction"]
    names = ["source", "seed", "interaction", "residual"]
    colors = [LBLUE, BLUE, NAVY, LGRAY]
    bottom = 0
    for c, n, col in zip(comps, names, colors):
        v = 100 * float(r[c])
        ax.bar(0, v, bottom=bottom, color=col, edgecolor=NAVY, label=n)
        bottom += v
    ax.set_xticks([0])
    ax.set_xticklabels(["Δ totalVI−scVI\nmacro-F1"])
    ax.set_ylabel("% variance")
    ax.set_ylim(0, 100)
    ax.set_title("C. Paired Δ variance components", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)

    # D Lawlor risk coverage quick
    ax = axes[1, 1]
    coverages = [100, 90, 80, 70, 60, 50]
    for m, col, lab in [("scvi_matched", LBLUE, "scVI"), ("totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        risks = [r[f"risk_{c}_mean"] for c in coverages]
        ax.plot(coverages, risks, "o-", color=col, label=lab)
    ax.set_xlabel("Coverage (%)")
    ax.set_ylabel("Risk")
    ax.set_title("D. Lawlor risk–coverage", fontweight="bold")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Figure5_uncertainty_stochasticity")


def fig6():
    style()
    delta = pd.read_csv(REV2 / "statistics/lawlor_donor_delta.csv")
    law_base = pd.read_csv(REV2 / "controls/tables/simple_multimodal_lawlor.csv")
    law = pd.read_csv(REV2 / "calibration/lawlor_reliability_summary_by_model.csv")
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8))

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A. Lawlor design", fontweight="bold")
    ax.text(0.05, 0.85, "10 biological donors\n5 folds · 10 seeds/fold\nPrimary unit = DONOR\n\nAuthor labels are\nprotein-gated.",
            va="top", fontsize=9, color=NAVY, bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    dlt = delta[(delta.metric == "present_class_macro_f1") & (delta.donor != "SUMMARY")]
    ax = axes[0, 1]
    ax.scatter(dlt.scvi_mean_over_seeds, dlt.totalvi_mean_over_seeds, c=NAVY, s=40)
    lims = [0.75, 0.92]
    ax.plot(lims, lims, "--", color=GRAY)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("scVI F1")
    ax.set_ylabel("totalVI F1")
    ax.set_title("B. Per-donor paired F1 (10/10)", fontweight="bold")

    ax = axes[0, 2]
    ax.barh(dlt.donor, dlt.delta, color=BLUE, edgecolor=NAVY)
    ax.axvline(0.090, color=NAVY, ls="--", lw=1)
    ax.set_xlabel("ΔF1")
    ax.set_title("C. Donor ΔF1 (mean≈+0.090)", fontweight="bold")

    ax = axes[1, 0]
    means = {
        "RNA": law_base[law_base.representation == "RNA_PCA_10"].macro_f1.mean(),
        "protein": law_base[law_base.representation == "protein_PCA_10"].macro_f1.mean(),
        "concat": law_base[law_base.representation == "concat_PCA_10_10"].macro_f1.mean(),
        "scVI": 0.800,
        "totalVI": 0.888,
    }
    cols = [LGRAY, BLUE, NAVY, LBLUE, NAVY]
    ax.bar(range(5), list(means.values()), color=cols, edgecolor=NAVY)
    ax.set_xticks(range(5))
    ax.set_xticklabels(list(means.keys()), rotation=20, ha="right")
    ax.set_ylim(0.6, 1.0)
    ax.set_ylabel("macro-F1")
    ax.set_title("D. Methods (protein baselines highest)", fontweight="bold")
    for i, v in enumerate(means.values()):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=7)

    ax = axes[1, 1]
    x = np.arange(2)
    for off, m, col, lab in [(-0.18, "scvi_matched", LBLUE, "scVI"), (0.18, "totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        ax.bar(x + off, [r.error_auroc_mean, r.latent_error_auroc_mean], width=0.34, color=col, edgecolor=NAVY, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(["predictive", "latent"])
    ax.set_title("E. Uncertainty AUROC", fontweight="bold")
    ax.legend(frameon=False)

    ax = axes[1, 2]
    for off, m, col, lab in [(-0.18, "scvi_matched", LBLUE, "scVI"), (0.18, "totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        ax.bar(np.arange(3) + off, [r.full_aurc_mean, r["paurc_0.5_1.0_mean"], r.e_aurc_mean], width=0.34, color=col, edgecolor=NAVY, label=lab)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["full AURC", "pAURC", "E-AURC"], rotation=15, ha="right")
    ax.set_title("F. Selective areas", fontweight="bold")
    ax.legend(frameon=False)

    fig.text(0.5, 0.01, "Author labels are protein-gated: strong protein PCA / concat performance indicates label–modality alignment, not universal totalVI superiority.",
             ha="center", fontsize=8, color=GRAY, style="italic")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save(fig, "Figure6_lawlor_external")


def main():
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    print("All figures ->", OUT)


if __name__ == "__main__":
    main()
