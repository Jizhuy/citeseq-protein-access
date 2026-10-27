#!/usr/bin/env python3
"""Redraw of main Figures 2-6 from frozen result tables (plotting only).

Every plotted value is read from a frozen CSV under results/ or copied from a hard-coded constant
documented in docs/RESULT_MANIFEST.md. No analysis is rerun. Writes vector PDFs (plus preview PNGs)
to figures/source/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "figures" / "source"

NAVY, BLUE, LBLUE, GRAY, LGRAY, WHITE = "#1B3A5F", "#2F6FAE", "#9DC3E6", "#5A5A5A", "#D9DEE7", "#FFFFFF"
RED = "#C0392B"
BOX = dict(boxstyle="round,pad=0.4", fc="#EEF3F9", ec=BLUE, lw=0.8)
W = 7.2
MINUS = "\u2212"


def style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.titlesize": 8,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.unicode_minus": True,
            "figure.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def fmt(v: float, signed: bool = False) -> str:
    s = f"{v:+.3f}" if signed else f"{v:.3f}"
    return s.replace("-", MINUS)


def panel(ax, letter: str, title: str, dx: float = -30):
    """Bold panel letter at the upper-left corner; title left-aligned on its own line below it."""
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold", color=NAVY, pad=4)
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(dx, 16), textcoords="offset points",
                fontsize=11, fontweight="bold", ha="left", va="bottom", color="black")


def save(fig, name: str):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.04, facecolor=WHITE)
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight", pad_inches=0.04, facecolor=WHITE)
    plt.close(fig)
    print("wrote", name)


def fig2():
    style()
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.9), gridspec_kw={"width_ratios": [1.25, 1.0, 0.95]})
    ax = axes[0]
    primary = [("RNA PCA", 0.667, LGRAY), ("concat PCA", 0.745, BLUE), ("totalVI", 0.741, NAVY)]
    x = np.arange(3)
    ax.bar(x, [v for _, v, _ in primary], color=[c for *_, c in primary], edgecolor=NAVY, width=0.7, lw=0.7)
    for i, (_, v, _) in enumerate(primary):
        ax.text(i, v + 0.006, fmt(v), ha="center", va="bottom", fontsize=7, color=NAVY, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([p[0] for p in primary])
    ax.set_ylim(0.55, 0.80)
    ax.set_ylabel("Macro-F1")
    panel(ax, "A", "Primary train-only sequence")
    ax.text(0.0, -0.20, "Secondary: protein PCA 0.512 · scVI_matched 0.709 ·\nMOFA+ 0.549 (historical)",
            transform=ax.transAxes, fontsize=6.5, color=GRAY, va="top")

    ax = axes[1]
    vals, lo, hi = [0.078, -0.004], [0.042, -0.060], [0.131, 0.030]
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.bar([0, 1], vals, color=[BLUE, NAVY], edgecolor=NAVY, width=0.6, lw=0.7)
    ax.errorbar([0, 1], vals, yerr=[[v - l for v, l in zip(vals, lo)], [h - v for v, h in zip(vals, hi)]],
                fmt="none", ecolor=GRAY, capsize=3, lw=0.9)
    for i, (v, h) in enumerate(zip(vals, hi)):
        ax.text(i, h + 0.006, fmt(v, signed=True), ha="center", va="bottom", fontsize=7, color=NAVY, fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"G_access\nconcat{MINUS}RNA", f"G_assoc\ntotalVI{MINUS}concat"])
    ax.set_ylim(-0.075, 0.16)
    ax.set_ylabel("Δ macro-F1")
    panel(ax, "B", "Descriptive increments (±95% CI)")
    ax.text(0.5, -0.30, "G_assoc near null\n(CI crosses 0; prop>0 = 0.355)", ha="center", va="top",
            transform=ax.transAxes, fontsize=6.5, color=BLUE)

    ax = axes[2]
    ax.axis("off")
    panel(ax, "C", "Protocol note", dx=-6)
    ax.text(0.0, 0.92,
            "Strict train-only (PRIMARY):\n"
            "fit PCA / scVI_matched / totalVI\n"
            "on 7,573 HC train cells only;\n"
            "encode 1,921 held-out cells\nafter fit.\n\n"
            "Historical feature-transductive\n"
            "deep-model result moved to\n"
            "Supplementary sensitivity.",
            va="top", fontsize=7, color=NAVY, bbox=BOX, linespacing=1.3)
    fig.tight_layout(w_pad=2.0)
    save(fig, "Figure2_protein_access_vs_model")


def fig3():
    style()
    mt = pd.read_csv(RES / "transfer/matched_concat_control/tables/matched_transfer_control.csv")
    fig, axes = plt.subplots(2, 2, figsize=(W, 5.9))
    ax = axes[0, 0]
    ax.axis("off")
    panel(ax, "A", "Different target-data access", dx=-6)
    ax.text(0.0, 0.94,
            "SOURCE-ONLY inductive concat\n  · fit on SOURCE only\n  · transform TARGET\n\n"
            "scVI_matched / totalVI (scArches)\n  · unlabeled TARGET features\n    used in query adaptation\n\n"
            "Feature-access-matched\ntransductive concat\n  · joint unsupervised fit on\n    SOURCE + unlabeled TARGET\n"
            "  · FEATURE ACCESS matched only\n  · NOT algorithmically identical",
            va="top", fontsize=7, color=NAVY, bbox=BOX, linespacing=1.25)

    def grouped(ax, direction, letter, title, note):
        order = [
            ("RNA_PCA_10", "none_inductive_source_fit_only", "RNA\nind.", LGRAY),
            ("concat_PCA_10_10", "none_inductive_source_fit_only", "concat\nind.", LBLUE),
            ("concat_PCA_10_10", "transductive_joint_unsupervised_feature_fit", "concat\ntrans.", BLUE),
            ("scVI_matched", "scArches_unsupervised_transductive", "scVI_matched\nscArches", LBLUE),
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
        ax.bar(x, vals, yerr=errs, color=cols, edgecolor=NAVY, capsize=2.5, width=0.72, lw=0.7,
               error_kw=dict(lw=0.8))
        ax.set_xticks(x)
        ax.set_xticklabels(labs, fontsize=6.6)
        ax.set_ylabel("Macro-F1")
        panel(ax, letter, title)
        ax.set_ylim(0.0, 1.05)
        ax.set_yticks(np.arange(0, 0.91, 0.1))
        for i, (v, e) in enumerate(zip(vals, errs)):
            ax.text(i, v + e + 0.015, fmt(v), ha="center", va="bottom", fontsize=6.8, color=NAVY)
        ax.text(0.98, 0.985, note, transform=ax.transAxes, ha="right", va="top", fontsize=7, color=NAVY,
                fontweight="bold", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=BLUE, lw=0.8))

    grouped(axes[0, 1], "direction_A", "B", "Direction A", f"totalVI {MINUS} concat_trans ≈ +0.016")
    grouped(axes[1, 0], "direction_B", "C", "Direction B", f"totalVI {MINUS} concat_trans ≈ {MINUS}0.113")

    ax = axes[1, 1]
    labs = [f"Dir A\ntotalVI {MINUS}\nscVI_matched", f"Dir A\ntotalVI {MINUS}\nconcat_tr",
            f"Dir B\ntotalVI {MINUS}\nscVI_matched", f"Dir B\ntotalVI {MINUS}\nconcat_tr"]
    vals = [-0.0032, 0.0161, 0.0525, -0.1131]
    cols = [LGRAY, BLUE, NAVY, RED]
    ax.bar(np.arange(4), vals, color=cols, edgecolor=NAVY, lw=0.7)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(labs, fontsize=6.8)
    ax.set_ylabel("Δ macro-F1")
    panel(ax, "D", "Feature-access contrasts (direction-dependent)")
    ax.set_ylim(-0.14, 0.08)
    for i, v in enumerate(vals):
        y = v + (0.006 if v >= 0 else -0.006)
        ax.text(i, y, fmt(v, signed=True), ha="center", va="bottom" if v >= 0 else "top", fontsize=7,
                color=NAVY, fontweight="bold")
    fig.tight_layout(h_pad=2.6, w_pad=2.0)
    fig.text(0.75, -0.005, "Feature-access matching does not yield a stable totalVI advantage",
             ha="center", va="top", fontsize=7, color=NAVY)
    save(fig, "Figure3_transfer_fairness")


def fig4():
    style()
    summ = pd.read_csv(RES / "degradation/training_time/training_time_degradation_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(W, 4.3), gridspec_kw={"height_ratios": [0.32, 1.0]})
    ax = axes[0, 0]
    ax.axis("off")
    panel(ax, "A", "Training-time degradation", dx=-6)
    ax.text(0.0, 0.88,
            "Corrupt protein → retrain totalVI\n(scVI_matched fixed RNA-only reference)\n"
            "5 mask/retraining replicates\n(nonzero levels; model_seed fixed)",
            va="top", fontsize=7.2, color=NAVY, bbox=BOX, linespacing=1.3)
    ax = axes[0, 1]
    ax.axis("off")
    panel(ax, "B", "Post-adaptation target-only", dx=-6)
    ax.text(0.0, 0.88, "Clean train/adapt → freeze\n→ corrupt TARGET protein only\nPrimary: ΔF1 from clean",
            va="top", fontsize=7.2, color=NAVY, bbox=BOX, linespacing=1.3)

    ax = axes[1, 0]
    fr = summ["corruption_fraction"].to_numpy()
    mu = summ["mean_macro_f1"].to_numpy()
    sd = summ["sd_macro_f1"].fillna(0).to_numpy()
    ax.errorbar(fr, mu, yerr=sd, fmt="o-", color=NAVY, capsize=2.5, ms=4, lw=1.1, label="totalVI mean ± SD")
    ax.axhline(0.722, color=LBLUE, ls="--", lw=1.1, label="scVI_matched fixed")
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("Macro-F1")
    panel(ax, "C", "Training-time F1 (mask/retraining replicates)")
    ax.legend(frameon=False, loc="upper right")
    ax.text(0.03, 0.04, "Non-monotonic 0.70→0.85 within replicate variability", transform=ax.transAxes,
            fontsize=6.5, color=GRAY)

    ax = axes[1, 1]
    frs = [0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85]
    a = [0.652, 0.650, 0.647, 0.644, 0.641, 0.637, 0.635]
    b = [0.654, 0.652, 0.647, 0.643, 0.637, 0.629, 0.621]
    ax.plot(frs, a, "o-", color=NAVY, ms=4, lw=1.1, label="Dir A totalVI")
    ax.plot(frs, b, "o-", color=BLUE, ms=4, lw=1.1, label="Dir B totalVI")
    ax.axhline(0.656, color=GRAY, ls=":", lw=1.1, label="Dir A scVI_matched clean")
    ax.axhline(0.601, color=LBLUE, ls=":", lw=1.3, label="Dir B scVI_matched clean")
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("Macro-F1")
    panel(ax, "D", "Post-adaptation: advantage narrows")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.42))
    fig.tight_layout(h_pad=1.6, w_pad=2.4)
    save(fig, "Figure4_protein_degradation")


def fig5():
    style()
    dvar = pd.read_csv(RES / "stochasticity/variance_decomposition/statistics/delta_variance_components_corrected.csv")
    law = pd.read_csv(RES / "uncertainty/calibration/lawlor_reliability_summary_by_model.csv")
    pb = pd.read_csv(RES / "uncertainty/calibration/pbmc_20seed_reliability_summary.csv")
    sm = pd.read_csv(RES / "degradation/test_time_corruption/tables/test_time_corruption_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(W, 5.0))

    ax = axes[0, 0]
    for m, col, lab in [("scvi_matched", LBLUE, "scVI_matched"), ("totalvi", NAVY, "totalVI")]:
        pred, lat = [], []
        for d in ["PBMC_direction_A", "PBMC_direction_B"]:
            r = pb[(pb.experiment == d) & (pb.model == m)].iloc[0]
            pred.append(r.error_auroc_mean)
            lat.append(r.latent_error_auroc_mean)
        ax.plot([0, 1], pred, "o-", color=col, ms=4, lw=1.1, label=f"{lab} predictive")
        ax.plot([0, 1], lat, "s--", color=col, alpha=0.75, ms=4, lw=1.1, label=f"{lab} latent")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(0.30, 0.90)
    ax.set_ylabel("Error AUROC")
    panel(ax, "A", "Predictive vs latent AUROC")
    ax.legend(frameon=False, ncol=2, loc="lower left", fontsize=6.5, columnspacing=1.0)

    ax = axes[0, 1]
    for d, c in [("direction_A", NAVY), ("direction_B", BLUE)]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.e_aurc_mean, "o-", color=c, ms=4, lw=1.1,
                label=d.replace("direction_", "Dir ") + " totalVI")
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("E-AURC")
    panel(ax, "B", "Post-adaptation E-AURC")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1, 0]
    r = dvar[(dvar.model == "delta_totalvi_minus_scvi") & (dvar.metric == "macro_f1")].iloc[0]
    comps = ["source_subset_fraction", "model_seed_fraction", "interaction_fraction", "residual_fraction"]
    names = ["source", "seed", "interaction", "residual"]
    colors = [LBLUE, BLUE, NAVY, LGRAY]
    bottom = 0.0
    for c, n, col in zip(comps, names, colors):
        v = 100 * float(r[c])
        ax.bar(0, v, bottom=bottom, color=col, edgecolor=NAVY, lw=0.7, width=0.6, label=n)
        bottom += v
    ax.set_xticks([0])
    ax.set_xticklabels([f"Δ totalVI {MINUS} scVI_matched\nmacro-F1"])
    ax.set_xlim(-0.6, 1.3)
    ax.set_ylabel("% variance")
    ax.set_ylim(0, 100)
    panel(ax, "C", "Paired Δ variance components")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], frameon=False, loc="center right")

    ax = axes[1, 1]
    coverages = [100, 90, 80, 70, 60, 50]
    for m, col, lab in [("scvi_matched", LBLUE, "scVI_matched"), ("totalvi", NAVY, "totalVI")]:
        rr = law[law.model == m].iloc[0]
        ax.plot(coverages, [rr[f"risk_{c}_mean"] for c in coverages], "o-", color=col, ms=4, lw=1.1, label=lab)
    ax.set_xlabel("Coverage (%)")
    ax.set_ylabel("Risk")
    panel(ax, "D", "Lawlor risk–coverage")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout(h_pad=2.2, w_pad=2.4)
    save(fig, "Figure5_uncertainty_stochasticity")


def fig6():
    style()
    summary = pd.read_csv(RES / "donor_validation/summary/lawlor_all_methods_donor_level_summary.csv")
    per = pd.read_csv(RES / "donor_validation/summary/lawlor_all_methods_per_donor.csv")
    order = ["RNA_PCA_10", "protein_PCA_10", "concat_PCA_10_10", "scVI_matched", "totalVI"]
    labels = ["RNA PCA", "protein PCA", "concat PCA", "scVI_matched", "totalVI"]
    colors = [LGRAY, "#A8C5D8", BLUE, LBLUE, NAVY]
    means = [float(summary.loc[summary.representation == r, "mean"].iloc[0]) for r in order]
    stds = [float(summary.loc[summary.representation == r, "std"].iloc[0]) for r in order]

    fig, axes = plt.subplots(1, 2, figsize=(W, 3.1), gridspec_kw={"width_ratios": [1.2, 1.0]})
    ax = axes[0]
    x = np.arange(len(order))
    ax.bar(x, means, yerr=stds, color=colors, edgecolor=NAVY, capsize=2.5, width=0.72, lw=0.7, error_kw=dict(lw=0.8))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Present-class macro-F1")
    ax.set_ylim(0.55, 1.0)
    panel(ax, "A", "Donor-level means (n = 10 biological donors)")
    for i, v in enumerate(means):
        ax.text(i, v + stds[i] + 0.008, fmt(v), ha="center", va="bottom", fontsize=7, color=NAVY)
    ax.text(0.0, -0.33, "Error bars: SD across donors. Aggregation unit: biological donor.",
            transform=ax.transAxes, fontsize=6.5, color=GRAY, va="top")

    ax = axes[1]
    wide = per.pivot_table(index="donor", columns="representation", values="macro_f1")
    delta = (wide["totalVI"] - wide["scVI_matched"]).sort_values()
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.barh(np.arange(len(delta)), delta.values, color=NAVY, edgecolor=NAVY, height=0.7)
    ax.set_yticks(np.arange(len(delta)))
    ax.set_yticklabels(delta.index.astype(str), fontsize=7)
    ax.set_xlim(0, 0.16)
    ax.set_xlabel(f"Δ present-class macro-F1 (totalVI {MINUS} scVI_matched)")
    panel(ax, "B", "Paired donor differences (10/10 > 0)", dx=-46)
    ax.text(0.98, 0.04, f"mean Δ ≈ +0.090\nsign test P = 0.00195", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=7, color=NAVY, bbox=dict(boxstyle="round,pad=0.3", fc="#EEF3F9", ec=BLUE, lw=0.8))
    fig.tight_layout(w_pad=2.4)
    save(fig, "Figure6_lawlor_external")


if __name__ == "__main__":
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
