#!/usr/bin/env python3
"""Generate main Figures 1–6 from frozen revision_round2 tables only."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "results" / "revision_round2_methodological_fixes"
OUT = ROOT / "results" / "final_application_manuscript_v2" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Palette: blue / navy / gray / white
NAVY = "#1B3A5F"
BLUE = "#2F6FAE"
LBLUE = "#9DC3E6"
GRAY = "#5A5A5A"
LGRAY = "#D9DEE7"
WHITE = "#FFFFFF"
TEAL = "#3A7CA5"  # still blue family


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.facecolor": WHITE,
            "figure.facecolor": WHITE,
            "axes.edgecolor": GRAY,
            "text.color": NAVY,
            "axes.labelcolor": NAVY,
            "xtick.color": GRAY,
            "ytick.color": GRAY,
            "grid.color": LGRAY,
            "grid.linewidth": 0.5,
        }
    )


def save(fig, name: str):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)
    print("wrote", name)


def fig1_study_design():
    style()
    fig, ax = plt.subplots(figsize=(11.0, 7.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("Figure 1. Study design and evaluation units", loc="left", fontweight="bold", color=NAVY)

    def box(x, y, w, h, text, fc=LBLUE, ec=NAVY):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", linewidth=1.2, facecolor=fc, edgecolor=ec)
        ax.add_patch(p)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5, color=NAVY, wrap=True)

    box(0.3, 6.4, 2.4, 1.2, "Development\nPBMC10k / PBMC5k\nRNA + protein\nn=10,849 cells")
    box(3.0, 6.4, 2.6, 1.2, "Representations\nRNA PCA · protein PCA\nconcat PCA · MOFA+\nscVI · totalVI")
    box(5.9, 6.4, 2.6, 1.2, "Internal HC benchmark\nn=9,494 cells\nRNA-only labels\nagreement=1.0")
    box(8.8, 6.4, 2.8, 1.2, "External Lawlor CITE-seq\n10 donors · 5 folds\ndonor-level biology\nprotein-gated labels")

    box(0.3, 4.4, 3.5, 1.4, "Reciprocal transfer\nA: PBMC10k→PBMC5k\nB: PBMC5k→PBMC10k\nscArches: unsupervised\ntransductive adaptation")
    box(4.1, 4.4, 3.6, 1.4, "Protein degradation\n(1) training-time retraining\n(2) post-adaptation\ntarget-only corruption\nfractions 0–0.85")
    box(8.0, 4.4, 3.6, 1.4, "Reliability diagnostics\nclassifier predictive U\nlatent posterior variance\nNLL · Brier · ECE\nAURC / pAURC / E-AURC")

    box(0.3, 2.0, 5.5, 1.8, "Computational replication\n• model seeds (e.g., 20 transfer seeds)\n• fold×seed pairs (Lawlor: 50 paired comparisons)\n• run replicates in variance design\nNOT biological replicates", fc="#EEF3F9")
    box(6.2, 2.0, 5.4, 1.8, "Biological replication\n• Lawlor: n=10 donors (primary unit)\n• between-donor SD of ΔF1 ≈ 0.010\n• within-donor seed SD smaller\nDonor-level inference preferred", fc="#EEF3F9")

    box(0.3, 0.3, 11.3, 1.3, "Scope note: selected CITE-seq representations (PCA baselines, MOFA+, scVI, totalVI).\nNot a universal multi-omics method benchmark. Architecture-specific claims are softened (STATE 2):\npaired protein information is primary; residual totalVI advantages are context-dependent.", fc=WHITE, ec=BLUE)

    for x0, x1 in [(2.7, 3.0), (5.6, 5.9), (8.5, 8.8)]:
        ax.annotate("", xy=(x1, 7.0), xytext=(x0, 7.0), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))
    save(fig, "Figure1_study_design")


def fig2_information_vs_representation():
    style()
    dev = pd.read_csv(REV / "controls/tables/simple_multimodal_development.csv")
    order = ["RNA_PCA_10", "protein_PCA_10", "concat_PCA_10_10", "MOFA+", "scVI_matched", "totalVI"]
    labels = ["RNA PCA", "protein PCA", "concat PCA", "MOFA+", "scVI", "totalVI"]
    vals = [float(dev.loc[dev.representation == r, "macro_f1"].iloc[0]) for r in order]
    colors = [LGRAY, LGRAY, BLUE, LGRAY, LBLUE, NAVY]

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), gridspec_kw={"width_ratios": [1.35, 0.9, 0.9]})

    ax = axes[0]
    x = np.arange(len(vals))
    ax.bar(x, vals, color=colors, edgecolor=NAVY, linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("macro-F1")
    ax.set_ylim(0.45, 0.82)
    ax.set_title("A. Development HC representations", fontweight="bold")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=7, color=NAVY)
    ax.axhline(0.745, color=BLUE, ls="--", lw=0.8, alpha=0.7)

    ax = axes[1]
    deltas = [0.779 - 0.722, 0.779 - 0.745]
    ax.bar([0, 1], deltas, color=[NAVY, BLUE], edgecolor=NAVY, width=0.65)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["totalVI−scVI", "totalVI−concat"], rotation=20, ha="right")
    ax.set_ylabel("Δ macro-F1")
    ax.set_title("B. Residual differences", fontweight="bold")
    ax.axhline(0, color=GRAY, lw=0.8)
    for i, v in enumerate(deltas):
        ax.text(i, v + 0.002, f"+{v:.3f}", ha="center", fontsize=8)

    ax = axes[2]
    ax.axis("off")
    ax.set_title("C. RNA-only label audit", fontweight="bold")
    txt = (
        "Development annotation provenance\n"
        "was already RNA-only\n"
        "(SCANVI / scArches RNA query).\n\n"
        "High-confidence cells: n = 9,494\n"
        "Agreement with historical labels:\n"
        "l1 = 1.0, l2 = 1.0\n\n"
        "Internal totalVI advantage is not\n"
        "explained by protein-informed\n"
        "label circularity."
    )
    ax.text(0.05, 0.92, txt, va="top", ha="left", fontsize=8.5, color=NAVY, family="DejaVu Serif",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#EEF3F9", edgecolor=BLUE))

    fig.tight_layout()
    save(fig, "Figure2_information_vs_representation")


def fig3_protein_degradation():
    style()
    train = pd.read_csv(ROOT / "results/tables/protein_sparsity_primary_delta_summary.csv")
    tt = pd.read_csv(REV / "controls/test_time_corruption/tables/training_vs_testtime_corruption.csv")
    sm = pd.read_csv(REV / "controls/test_time_corruption/tables/test_time_corruption_summary.csv")
    pc = pd.read_csv(REV / "controls/test_time_corruption/tables/test_time_corruption_per_class.csv")

    fig = plt.figure(figsize=(11.2, 8.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.85, 1.1, 1.15], hspace=0.45, wspace=0.28)

    # A design training-time
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.set_title("A. Training-time degradation", fontweight="bold", loc="left")
    ax.text(0.02, 0.55, "Corrupt training protein → retrain totalVI → evaluate\n(scVI fixed reference)\nWithin-dataset HC development setting",
            fontsize=8.5, color=NAVY, va="center",
            bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    # B design post-adapt
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("B. Post-adaptation target-only degradation", fontweight="bold", loc="left")
    ax.text(0.02, 0.55, "Clean source train → clean scArches adapt → freeze\n→ corrupt TARGET protein only → re-infer\nPrimary endpoint: Δ from clean within seed",
            fontsize=8.5, color=NAVY, va="center",
            bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    # C Dir A F1
    ax = fig.add_subplot(gs[1, 0])
    for d, c, ls in [("direction_A", NAVY, "-"), ("direction_B", BLUE, "-")]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.macro_f1_mean, "o-", color=c, label=d.replace("direction_", "Dir "), lw=1.6)
        # scVI clean ref
        ref = tt.loc[tt.direction == d, "scVI_reference_clean_transfer"].iloc[0]
        ax.axhline(ref, color=c, ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("macro-F1")
    ax.set_title("C–D. Post-adaptation totalVI F1 (dashed = clean scVI)", fontweight="bold")
    ax.legend(frameon=False)

    # E delta
    ax = fig.add_subplot(gs[1, 1])
    for d, c in [("direction_A", NAVY), ("direction_B", BLUE)]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.delta_macro_f1_mean, "o-", color=c, label=d.replace("direction_", "Dir "), lw=1.6)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("ΔF1 from clean")
    ax.set_title("E. ΔF1 from clean (primary)", fontweight="bold")
    ax.legend(frameon=False)

    # F per-class Dir B at 0.85
    ax = fig.add_subplot(gs[2, :])
    tvi = pc[(pc.model == "totalvi") & (pc.direction == "direction_B")].copy()
    clean = tvi[tvi.corruption_fraction == 0].groupby("cell_type")["f1"].mean()
    high = tvi[tvi.corruption_fraction == 0.85].groupby("cell_type").agg(f1=("f1", "mean"), support=("support", "median"))
    df = high.join(clean.rename("f1_clean"), how="inner")
    df["delta"] = df["f1"] - df["f1_clean"]
    df = df[df.support.fillna(0) >= 50].sort_values("delta")
    focus = df.head(10)
    ax.barh(focus.index.astype(str), focus.delta, color=BLUE, edgecolor=NAVY)
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("ΔF1 at p=0.85 vs clean (Dir B; support≥50)")
    ax.set_title("F. Most degradation-sensitive cell types (Direction B)", fontweight="bold")

    save(fig, "Figure4_protein_degradation")


def fig4_transfer_uncertainty():
    style()
    tr = pd.read_csv(REV / "controls/tables/simple_multimodal_transfer.csv")
    pb = pd.read_csv(REV / "calibration/pbmc_20seed_reliability_summary.csv")

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.6))
    # A schematic
    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A. Transfer designs", fontweight="bold")
    ax.text(0.05, 0.75, "Direction A: PBMC10k → PBMC5k\nDirection B: PBMC5k → PBMC10k\n\nDeep models: unsupervised\ntransductive scArches adaptation\n\nConcat PCA: inductive\nsource-fitted transform only",
            fontsize=8, color=NAVY, va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE))

    # B 20-seed F1
    ax = axes[0, 1]
    rows = []
    for d in ["direction_A", "direction_B"]:
        for m, col in [("scVI", LBLUE), ("totalVI", NAVY)]:
            key = "scvi_matched" if m == "scVI" else "totalvi"
            r = pb[(pb.experiment == f"PBMC_{d}") & (pb.model == key)].iloc[0]
            rows.append((d, m, r.macro_f1_mean, r.macro_f1_sd, col))
    xpos = [0, 1, 3, 4]
    for i, (d, m, mu, sd, col) in enumerate(rows):
        ax.bar(xpos[i], mu, yerr=sd, color=col, edgecolor=NAVY, width=0.8, capsize=3)
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_ylabel("macro-F1")
    ax.set_title("B. 20-seed transfer (mean±SD)", fontweight="bold")
    ax.legend(handles=[mpatches.Patch(color=LBLUE, label="scVI"), mpatches.Patch(color=NAVY, label="totalVI")], frameon=False)

    # C concat
    ax = axes[0, 2]
    for i, d in enumerate(["direction_A", "direction_B"]):
        vals = []
        labs = ["RNA", "protein", "concat", "scVI", "totalVI"]
        mapr = {
            "RNA": "RNA_PCA_10",
            "protein": "protein_PCA_10",
            "concat": "concat_PCA_10_10",
            "scVI": "scVI_matched",
            "totalVI": "totalVI",
        }
        g = tr[tr.direction == d]
        ys = [float(g.loc[g.representation == mapr[l], "macro_f1"].iloc[0]) for l in labs]
        ax.plot(np.arange(len(labs)), ys, "o-", color=[NAVY, BLUE][i], label=d.replace("direction_", "Dir "))
    ax.set_xticks(range(5))
    ax.set_xticklabels(labs, rotation=25, ha="right")
    ax.set_ylabel("macro-F1")
    ax.set_title("C. Inductive concat vs adapted VAEs", fontweight="bold")
    ax.legend(frameon=False)

    # D predictive AUROC
    ax = axes[1, 0]
    for i, d in enumerate(["PBMC_direction_A", "PBMC_direction_B"]):
        for m, col, off in [("scvi_matched", LBLUE, -0.15), ("totalvi", NAVY, 0.15)]:
            r = pb[(pb.experiment == d) & (pb.model == m)].iloc[0]
            ax.bar(i + off, r.error_auroc_mean, width=0.28, color=col, edgecolor=NAVY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_ylabel("error AUROC")
    ax.set_ylim(0.8, 0.9)
    ax.set_title("D. Predictive error AUROC", fontweight="bold")

    # E latent AUROC
    ax = axes[1, 1]
    for i, d in enumerate(["PBMC_direction_A", "PBMC_direction_B"]):
        for m, col, off in [("scvi_matched", LBLUE, -0.15), ("totalvi", NAVY, 0.15)]:
            r = pb[(pb.experiment == d) & (pb.model == m)].iloc[0]
            ax.bar(i + off, r.latent_error_auroc_mean, width=0.28, color=col, edgecolor=NAVY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_ylabel("error AUROC")
    ax.set_ylim(0.4, 0.75)
    ax.set_title("E. Latent posterior error AUROC", fontweight="bold")

    # F note
    ax = axes[1, 2]
    ax.axis("off")
    ax.set_title("F. Uncertainty definitions", fontweight="bold")
    ax.text(
        0.02,
        0.9,
        "Classifier predictive uncertainty:\n"
        "U_pred = 1 − max_k p(y=k | z)\n\n"
        "Latent posterior variance:\n"
        "U_latent = (1/d) Σ_j Var[z_j | x]\n"
        "d = 20\n\n"
        "These are distinct quantities;\n"
        "do not equate latent variance with\n"
        "epistemic uncertainty.",
        va="top",
        fontsize=8,
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE),
    )
    fig.tight_layout()
    save(fig, "Figure3_transfer_and_uncertainty")


def fig5_stochasticity_selective():
    style()
    var = pd.read_csv(REV / "statistics/variance_components_corrected.csv")
    dvar = pd.read_csv(REV / "statistics/delta_variance_components_corrected.csv")
    law = pd.read_csv(REV / "calibration/lawlor_reliability_summary_by_model.csv")
    sm = pd.read_csv(REV / "controls/test_time_corruption/tables/test_time_corruption_summary.csv")
    ret = pd.read_csv(ROOT / "results/phase11b_uncertainty_robustness/tables/phase11b_retention_imbalance.csv")

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8))

    def stack_components(ax, title, rows):
        comps = ["source_subset_fraction", "model_seed_fraction", "interaction_fraction", "residual_fraction"]
        names = ["source subset", "model seed", "interaction", "residual"]
        colors = [LBLUE, BLUE, NAVY, LGRAY]
        bottom = np.zeros(len(rows))
        x = np.arange(len(rows))
        for c, name, col in zip(comps, names, colors):
            vals = [100 * float(r[c]) for r in rows]
            ax.bar(x, vals, bottom=bottom, color=col, edgecolor=NAVY, linewidth=0.4, label=name)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels([r["_lab"] for r in rows], rotation=15, ha="right")
        ax.set_ylabel("% variance")
        ax.set_ylim(0, 100)
        ax.set_title(title, fontweight="bold")

    rows = []
    for m, lab in [("scvi_matched", "scVI"), ("totalvi", "totalVI")]:
        r = var[(var.model == m) & (var.metric == "macro_f1")].iloc[0]
        rows.append({**r.to_dict(), "_lab": lab})
    stack_components(axes[0, 0], "A. Within-model variance (macro-F1)", rows)
    axes[0, 0].legend(fontsize=6, frameon=False, loc="upper right")

    r = dvar[(dvar.model == "delta_totalvi_minus_scvi") & (dvar.metric == "macro_f1")].iloc[0]
    stack_components(axes[0, 1], "B. Paired Δ variance (macro-F1)", [{**r.to_dict(), "_lab": "Δ totalVI−scVI"}])

    # C risk-coverage Lawlor
    ax = axes[0, 2]
    coverages = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    for m, col, lab in [("scvi_matched", LBLUE, "scVI"), ("totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        risks = [r[f"risk_{int(c*100)}_mean"] for c in coverages]
        ax.plot([c * 100 for c in coverages], risks, "o-", color=col, label=lab)
    ax.set_xlabel("Coverage (%)")
    ax.set_ylabel("Risk (error rate)")
    ax.set_title("C. Lawlor risk–coverage", fontweight="bold")
    ax.legend(frameon=False)

    # D AURC trio
    ax = axes[1, 0]
    metrics = ["full_aurc_mean", "paurc_0.5_1.0_mean", "e_aurc_mean"]
    labs = ["full AURC", "pAURC[0.5,1]", "E-AURC"]
    x = np.arange(len(metrics))
    for off, m, col in [(-0.18, "scvi_matched", LBLUE), (0.18, "totalvi", NAVY)]:
        r = law[law.model == m].iloc[0]
        ax.bar(x + off, [r[k] for k in metrics], width=0.34, color=col, edgecolor=NAVY, label="scVI" if m.startswith("scvi") else "totalVI")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=15, ha="right")
    ax.set_ylabel("Area (lower better)")
    ax.set_title("D. Lawlor selective areas", fontweight="bold")
    ax.legend(frameon=False)

    # E retention CV
    ax = axes[1, 1]
    sub = ret[ret.coverage == 0.5]
    for i, d in enumerate(["direction_A", "direction_B"]):
        for m, col, off in [("scvi_matched", LBLUE, -0.15), ("totalvi", NAVY, 0.15)]:
            v = float(sub[(sub.direction == d) & (sub.model == m)].retention_cv.iloc[0])
            ax.bar(i + off, v, width=0.28, color=col, edgecolor=NAVY)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Dir A", "Dir B"])
    ax.set_ylabel("Retention CV")
    ax.set_title("E. Cell-type retention CV @50%", fontweight="bold")

    # F Part C E-AURC
    ax = axes[1, 2]
    for d, c in [("direction_A", NAVY), ("direction_B", BLUE)]:
        g = sm[(sm.model == "totalvi") & (sm.direction == d)].sort_values("corruption_fraction")
        ax.plot(g.corruption_fraction, g.e_aurc_mean, "o-", color=c, label=d.replace("direction_", "Dir "))
    ax.set_xlabel("Corruption fraction")
    ax.set_ylabel("E-AURC")
    ax.set_title("F. Post-adaptation E-AURC", fontweight="bold")
    ax.legend(frameon=False)

    fig.tight_layout()
    save(fig, "Figure5_stochasticity_and_selective")


def fig6_lawlor():
    style()
    delta = pd.read_csv(REV / "statistics/lawlor_donor_delta.csv")
    seedm = pd.read_csv(REV / "statistics/lawlor_donor_seedmean.csv")
    law_base = pd.read_csv(REV / "controls/tables/simple_multimodal_lawlor.csv")
    law = pd.read_csv(REV / "calibration/lawlor_reliability_summary_by_model.csv")

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8))

    ax = axes[0, 0]
    ax.axis("off")
    ax.set_title("A. Lawlor design", fontweight="bold")
    ax.text(
        0.05,
        0.85,
        "10 biological donors\n5 donor folds\n10 seeds / fold\n100 deep-model runs\n50 matched fold×seed pairs\n\nPrimary biological unit: DONOR\nFold×seed pairs ≠ biological n",
        va="top",
        fontsize=8.5,
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.35", fc="#EEF3F9", ec=BLUE),
    )

    # B paired F1
    ax = axes[0, 1]
    dlt = delta[(delta.metric == "present_class_macro_f1") & (delta.donor != "SUMMARY")].copy()
    donors = dlt.donor.tolist()
    ax.scatter(dlt.scvi_mean_over_seeds, dlt.totalvi_mean_over_seeds, c=NAVY, s=40, zorder=3)
    lims = [0.75, 0.92]
    ax.plot(lims, lims, "--", color=GRAY, lw=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("scVI present-class macro-F1")
    ax.set_ylabel("totalVI present-class macro-F1")
    ax.set_title("B. Per-donor paired F1 (10/10)", fontweight="bold")

    # C delta
    ax = axes[0, 2]
    ax.barh(donors, dlt.delta, color=BLUE, edgecolor=NAVY)
    ax.axvline(0.0902, color=NAVY, ls="--", lw=1, label="mean Δ=+0.090")
    ax.set_xlabel("Δ macro-F1 (totalVI−scVI)")
    ax.set_title("C. Per-donor ΔF1", fontweight="bold")
    ax.legend(frameon=False, fontsize=7)

    # D method comparison means
    ax = axes[1, 0]
    means = {
        "RNA PCA": law_base[law_base.representation == "RNA_PCA_10"].macro_f1.mean(),
        "protein PCA": law_base[law_base.representation == "protein_PCA_10"].macro_f1.mean(),
        "concat PCA": law_base[law_base.representation == "concat_PCA_10_10"].macro_f1.mean(),
        "scVI": 0.799820,
        "totalVI": 0.888304,
    }
    cols = [LGRAY, BLUE, NAVY, LBLUE, NAVY]
    ax.bar(range(len(means)), list(means.values()), color=cols, edgecolor=NAVY)
    ax.set_xticks(range(len(means)))
    ax.set_xticklabels(list(means.keys()), rotation=25, ha="right")
    ax.set_ylabel("macro-F1")
    ax.set_title("D. Methods incl. protein baselines", fontweight="bold")
    ax.set_ylim(0.6, 1.0)
    for i, v in enumerate(means.values()):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=7)

    # E AUROCs
    ax = axes[1, 1]
    x = np.arange(2)
    for off, m, col, lab in [(-0.18, "scvi_matched", LBLUE, "scVI"), (0.18, "totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        ax.bar(x + off, [r.error_auroc_mean, r.latent_error_auroc_mean], width=0.34, color=col, edgecolor=NAVY, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(["predictive", "latent"])
    ax.set_ylabel("error AUROC")
    ax.set_title("E. Uncertainty discrimination", fontweight="bold")
    ax.legend(frameon=False)

    # F AURC
    ax = axes[1, 2]
    for off, m, col, lab in [(-0.18, "scvi_matched", LBLUE, "scVI"), (0.18, "totalvi", NAVY, "totalVI")]:
        r = law[law.model == m].iloc[0]
        ax.bar(np.arange(3) + off, [r.full_aurc_mean, r["paurc_0.5_1.0_mean"], r.e_aurc_mean], width=0.34, color=col, edgecolor=NAVY, label=lab)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["full AURC", "pAURC", "E-AURC"], rotation=15, ha="right")
    ax.set_ylabel("Area (lower better)")
    ax.set_title("F. Selective prediction areas", fontweight="bold")
    ax.legend(frameon=False)

    # annotation about label modality
    fig.text(
        0.5,
        0.01,
        "Note: Lawlor author labels are protein-gated; strong protein PCA / concat performance indicates label–modality alignment, not universal architecture superiority.",
        ha="center",
        fontsize=8,
        color=GRAY,
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save(fig, "Figure6_lawlor_donor_validation")


def main():
    fig1_study_design()
    fig2_information_vs_representation()
    fig3_protein_degradation()
    fig4_transfer_uncertainty()
    fig5_stochasticity_selective()
    fig6_lawlor()
    print("All figures written to", OUT)


if __name__ == "__main__":
    main()
