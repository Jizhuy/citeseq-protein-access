#!/usr/bin/env python3
"""Regenerate Figure 1 and Figure 6 with final terminology and donor-level Lawlor aggregation."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
REV = Path(__file__).resolve().parents[1]
OUTS = [
    REV / "figures",
    ROOT / "results/final_information_vs_model_manuscript/figures",
]
# Parent-directory DOCX figure pack (may be outside sandbox write scope).
_PARENT_FIGS = ROOT.parent / "presubmission_revision" / "figures"
if _PARENT_FIGS.exists() or True:
    try:
        _PARENT_FIGS.mkdir(parents=True, exist_ok=True)
        OUTS.append(_PARENT_FIGS)
    except OSError:
        pass
for o in OUTS:
    o.mkdir(parents=True, exist_ok=True)

NAVY, BLUE, LBLUE, GRAY, LGRAY, WHITE = "#1B3A5F", "#2F6FAE", "#9DC3E6", "#5A5A5A", "#D9DEE7", "#FFFFFF"


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


def save(fig, name: str):
    for out in OUTS:
        for ext in ("pdf", "png"):
            path = out / f"{name}.{ext}"
            try:
                fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=WHITE)
            except OSError as exc:
                print(f"skip {path}: {exc}")
    plt.close(fig)
    print("wrote", name)


def fig1():
    style()
    fig, ax = plt.subplots(figsize=(11.2, 6.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title(
        "Figure 1. Protein-access gain vs residual totalVI-associated contrast",
        loc="left",
        fontweight="bold",
        color=NAVY,
    )

    def box(x, y, w, h, text, fc=LBLUE, tc=NAVY):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                lw=1.2,
                fc=fc,
                ec=NAVY,
            )
        )
        ax.text(
            x + w / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=8.5,
            color=tc,
            fontweight="bold" if fc == NAVY else "normal",
        )

    box(0.4, 5.8, 2.8, 1.5, "RNA only\nRNA PCA\nF1 = 0.667", LGRAY)
    box(4.0, 5.8, 3.2, 1.5, "Simple multimodal\nRNA+protein concat PCA\nF1 = 0.745", LBLUE)
    box(8.2, 5.8, 3.4, 1.5, "Complex multimodal\ntotalVI\nF1 = 0.741", NAVY, tc=WHITE)
    ax.annotate("", xy=(4.0, 6.55), xytext=(3.2, 6.55), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
    ax.annotate("", xy=(8.2, 6.55), xytext=(7.2, 6.55), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
    ax.text(3.6, 7.55, "Protein-access gain\nG_access = +0.078\n(strict train-only)", ha="center", fontsize=8, color=BLUE)
    ax.text(
        7.7,
        7.55,
        "Residual totalVI-associated contrast\nG_assoc = −0.004\n(near null; CI crosses 0)",
        ha="center",
        fontsize=8,
        color=NAVY,
    )

    box(
        0.4,
        3.2,
        11.2,
        1.8,
        "Stress-test the residual totalVI-associated contrast under:\n"
        "dataset transfer (inductive + feature-access-matched transductive)  ·  protein degradation  ·  training stochasticity  ·  uncertainty diagnostics\n"
        "External donor evidence (Lawlor): multimodal measurement value across 10 donors — not universal architecture superiority",
        "#EEF3F9",
    )
    box(
        0.4,
        0.5,
        11.2,
        2.0,
        "Descriptive contrasts only (not causal architecture effects; primary = strict train-only 7,573/1,921):\n"
        "G_access = F1_concat − F1_RNA_PCA    ·    G_assoc = F1_totalVI − F1_concat\n"
        "Central question: How much of the observed multimodal gain is associated with protein access,\n"
        "and is there a stable residual totalVI-associated contrast beyond a simple RNA–protein representation?\n"
        "(R_access is a secondary descriptive ratio only; not a recovery fraction.)",
        WHITE,
    )
    save(fig, "Figure1_central_logic")


def fig2():
    """Primary internal hierarchy + contrasts under strict train-only protocol."""
    style()
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), gridspec_kw={"width_ratios": [1.4, 0.95, 0.85]})
    ax = axes[0]
    primary = [("RNA PCA", 0.667, LGRAY), ("concat PCA", 0.745, BLUE), ("totalVI", 0.741, NAVY)]
    x = np.arange(3)
    ax.bar(x, [v for _, v, _ in primary], color=[c for *_, c in primary], edgecolor=NAVY, width=0.7)
    for i, (lab, v, _) in enumerate(primary):
        ax.text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=8, color=NAVY, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([p[0] for p in primary])
    ax.set_ylim(0.55, 0.82)
    ax.set_ylabel("macro-F1")
    ax.set_title("A. Primary train-only sequence", fontweight="bold")
    ax.text(
        0.02,
        0.06,
        "Secondary: protein 0.512 · scVI_matched 0.709\nMOFA+ 0.549 (historical)",
        transform=ax.transAxes,
        fontsize=7,
        color=GRAY,
    )

    ax = axes[1]
    vals = [0.078, -0.004]
    colors = [BLUE, NAVY]
    ax.axhline(0, color=GRAY, lw=0.9)
    ax.bar([0, 1], vals, color=colors, edgecolor=NAVY, width=0.65)
    ax.errorbar([0, 1], vals, yerr=[[0.078 - 0.042, -0.004 - (-0.060)], [0.131 - 0.078, 0.030 - (-0.004)]], fmt="none", ecolor=GRAY, capsize=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["G_access\nconcat−RNA", "G_assoc\ntotalVI−concat"])
    ax.set_ylabel("Δ macro-F1")
    ax.set_title("B. Descriptive increments (±95% CI)", fontweight="bold")
    ax.text(0, 0.088, "+0.078", ha="center", fontsize=8, color=BLUE)
    ax.text(1, 0.012, "−0.004", ha="center", fontsize=8, color=NAVY)
    ax.text(
        0.5,
        -0.24,
        "G_assoc near null\n(CI crosses 0; prop>0=0.355)",
        ha="center",
        transform=ax.transAxes,
        fontsize=7.5,
        color=BLUE,
    )

    ax = axes[2]
    ax.axis("off")
    ax.set_title("C. Protocol note", fontweight="bold")
    ax.text(
        0.05,
        0.9,
        "Strict train-only (PRIMARY):\n"
        "fit PCA / scVI_matched / totalVI\n"
        "on 7,573 HC train cells only;\n"
        "encode 1,921 held-out cells after fit.\n\n"
        "Historical feature-transductive\n"
        "deep-model result moved to\n"
        "Supplementary sensitivity.",
        va="top",
        fontsize=8.5,
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.4", fc="#EEF3F9", ec=BLUE),
    )
    fig.tight_layout()
    save(fig, "Figure2_protein_access_vs_model")
    # Keep legacy filename used by older figure packs.
    for out in OUTS:
        src = out / "Figure2_protein_access_vs_model.png"
        dst = out / "Figure2_information_vs_model.png"
        if src.exists():
            try:
                dst.write_bytes(src.read_bytes())
            except OSError:
                pass
        srcp = out / "Figure2_protein_access_vs_model.pdf"
        dstp = out / "Figure2_information_vs_model.pdf"
        if srcp.exists():
            try:
                dstp.write_bytes(srcp.read_bytes())
            except OSError:
                pass


def fig6():
    style()
    summary = pd.read_csv(REV / "tables/lawlor_all_methods_donor_level_summary.csv")
    per = pd.read_csv(REV / "tables/lawlor_all_methods_per_donor.csv")
    order = ["RNA_PCA_10", "protein_PCA_10", "concat_PCA_10_10", "scVI_matched", "totalVI"]
    labels = ["RNA PCA", "protein PCA", "concat PCA", "scVI_matched", "totalVI"]
    colors = [LGRAY, "#A8C5D8", BLUE, LBLUE, NAVY]
    means = [float(summary.loc[summary.representation == r, "mean"].iloc[0]) for r in order]
    stds = [float(summary.loc[summary.representation == r, "std"].iloc[0]) for r in order]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax = axes[0]
    x = np.arange(len(order))
    ax.bar(x, means, yerr=stds, color=colors, edgecolor=NAVY, capsize=3, width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("present-class macro-F1")
    ax.set_ylim(0.55, 1.0)
    ax.set_title("A. Donor-level means (n=10 biological donors)", fontweight="bold")
    for i, v in enumerate(means):
        ax.text(i, v + stds[i] + 0.01, f"{v:.3f}", ha="center", fontsize=8, color=NAVY)
    ax.text(
        0.02,
        0.02,
        "Error bars: SD across donors\nAggregation unit: biological donor",
        transform=ax.transAxes,
        fontsize=7.5,
        color=GRAY,
        va="bottom",
    )

    ax = axes[1]
    # paired donor deltas totalVI - scVI_matched
    wide = per.pivot_table(index="donor", columns="representation", values="macro_f1")
    delta = (wide["totalVI"] - wide["scVI_matched"]).sort_values()
    ax.axvline(0, color=GRAY, lw=1)
    ax.barh(np.arange(len(delta)), delta.values, color=NAVY, edgecolor=NAVY, height=0.7)
    ax.set_yticks(np.arange(len(delta)))
    ax.set_yticklabels(delta.index.astype(str), fontsize=8)
    ax.set_xlabel("Δ present-class macro-F1 (totalVI − scVI_matched)")
    ax.set_title("B. Paired donor differences (10/10 > 0)", fontweight="bold")
    ax.text(
        0.98,
        0.05,
        "mean Δ ≈ +0.090\nsign test P = 0.00195",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.3", fc="#EEF3F9", ec=BLUE),
    )
    fig.tight_layout()
    save(fig, "Figure6_lawlor_external")


if __name__ == "__main__":
    fig1()
    fig2()
    fig6()
