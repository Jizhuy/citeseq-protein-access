"""PHASE 10B figures: source-size sensitivity control."""

from __future__ import annotations

from pathlib import Path

import matplotlib

# The CUDA server is headless WSL2; the interactive Tk backend is unavailable.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save  # noqa: E402

CONDITION_ORDER = ("original_A", "size_matched_A", "original_B")
CONDITION_LABEL = {
    "original_A": "Original A\nPBMC10k→PBMC5k\nsource n=6855",
    "size_matched_A": "Size-matched A\nPBMC10k→PBMC5k\nsource n=3994",
    "original_B": "Original B\nPBMC5k→PBMC10k\nsource n=3994",
}


def write_phase10b_figures(tables_dir: Path, figures_dir: Path) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(tables_dir / "source_size_primary.csv")
    delta = pd.read_csv(tables_dir / "source_size_delta.csv")
    celltype = pd.read_csv(tables_dir / "source_size_celltype_specific.csv")
    summary = pd.read_csv(tables_dir / "source_size_summary.csv")
    rare = pd.read_csv(tables_dir / "source_size_rare_class_stability.csv")

    _fig_design(figures_dir / "phase10b_figure1_design")
    _fig_f1_conditions(delta, figures_dir / "phase10b_figure2_macro_f1")
    _fig_delta(delta, figures_dir / "phase10b_figure3_delta")
    _fig_size_effect(summary, figures_dir / "phase10b_figure4_size_effect")
    _fig_celltype(celltype, figures_dir / "phase10b_figure5_celltype_delta")
    _fig_rare(rare, figures_dir / "phase10b_figure6_rare_class_stability")
    _fig_neighbor(primary, figures_dir / "phase10b_figure7_neighbor_agreement")


def _common_l2(frame: pd.DataFrame) -> pd.DataFrame:
    return frame[(frame["label_level"] == "l2") & (frame["evaluation_class_set"] == "common")]


def _fig_design(prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 2.9))
    ax.set_axis_off()
    boxes = [
        (0.02, "PBMC10k\nsource\nn=6855"),
        (0.23, "random subsample\nn=3994\n(5 seeds, no labels)"),
        (0.46, "scVI_matched\nor totalVI\nmodel_seed=0"),
        (0.68, "scArches\nquery mapping"),
        (0.88, "PBMC5k\ntarget\nn=3994"),
    ]
    for x, text in boxes:
        ax.add_patch(
            plt.Rectangle((x, 0.26), 0.155, 0.52, facecolor="#d6e6f2", edgecolor=BLUE, linewidth=1.2)
        )
        ax.text(x + 0.0775, 0.52, text, ha="center", va="center", fontsize=8, color=GRAY)
    for x in (0.185, 0.415, 0.635, 0.845):
        ax.annotate("", xy=(x + 0.04, 0.52), xytext=(x, 0.52), arrowprops=dict(arrowstyle="->", color=BLUE))
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1)
    ax.set_title("PHASE 10B: source-size sensitivity control (Direction A only)")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_f1_conditions(delta: pd.DataFrame, prefix: Path) -> None:
    sub = _common_l2(delta)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    x = np.arange(len(CONDITION_ORDER))
    width = 0.35
    for i, (col, label, color) in enumerate(
        (("scvi_macro_f1", "scVI_matched", LIGHT_BLUE), ("totalvi_macro_f1", "totalVI", BLUE))
    ):
        means, sds = [], []
        for cond in CONDITION_ORDER:
            vals = sub[sub["condition"] == cond][col]
            means.append(float(vals.mean()))
            sds.append(float(vals.std(ddof=1)) if len(vals) > 1 else 0.0)
        ax.bar(x + (i - 0.5) * width, means, width, yerr=sds, color=color, label=label, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([CONDITION_LABEL[c] for c in CONDITION_ORDER], fontsize=8)
    ax.set_ylabel("Cross-dataset l2 logreg macro-F1 (common classes)")
    ax.set_title("Does size-matched Direction A move toward Direction B?")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_delta(delta: pd.DataFrame, prefix: Path) -> None:
    sub = _common_l2(delta)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for i, cond in enumerate(CONDITION_ORDER):
        part = sub[sub["condition"] == cond]
        mean = float(part["delta"].mean())
        ax.scatter(np.full(len(part), i), part["delta"], color=LIGHT_BLUE, zorder=3, s=34)
        lo = float(part["bootstrap_ci_low"].mean())
        hi = float(part["bootstrap_ci_high"].mean())
        if np.isfinite(lo) and np.isfinite(hi):
            ax.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="o", color=BLUE, capsize=4, zorder=4)
        else:
            ax.scatter([i], [mean], color=BLUE, zorder=4, s=60)
    ax.axhline(0, color=GRAY, linewidth=0.8)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels([CONDITION_LABEL[c] for c in CONDITION_ORDER], fontsize=8)
    ax.set_ylabel("Δ macro-F1 (totalVI − scVI_matched)")
    ax.set_title("Individual runs (light) and mean with paired bootstrap 95% CI (dark)")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_size_effect(summary: pd.DataFrame, prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    levels = list(summary["label_level"])
    x = np.arange(len(levels))
    width = 0.35
    for i, (col, label, color) in enumerate(
        (("scvi_size_effect", "scVI_matched", LIGHT_BLUE), ("totalvi_size_effect", "totalVI", BLUE))
    ):
        ax.bar(x + (i - 0.5) * width, summary[col], width, color=color, label=label)
    ax.axhline(0, color=GRAY, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(levels)
    ax.set_ylabel("macro-F1 change after downsampling source\n(3994 − 6855)")
    ax.set_title("Model-specific cost of removing source cells")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_celltype(celltype: pd.DataFrame, prefix: Path) -> None:
    piv = {}
    for cond in CONDITION_ORDER:
        part = celltype[celltype["condition"] == cond]
        if not len(part):
            continue
        m = part.groupby(["cell_type", "model"])["f1"].mean().unstack("model")
        if "totalvi" in m and "scvi_matched" in m:
            piv[cond] = m["totalvi"] - m["scvi_matched"]
    if not piv:
        return
    mat = pd.DataFrame(piv).sort_values("size_matched_A" if "size_matched_A" in piv else list(piv)[0])
    fig, ax = plt.subplots(figsize=(6.4, 0.30 * len(mat) + 1.9))
    vmax = float(np.nanmax(np.abs(mat.to_numpy()))) or 0.1
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([c.replace("_", "\n") for c in mat.columns], fontsize=8)
    ax.set_title("Δ F1 (totalVI − scVI_matched), common l2 classes")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Δ F1")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_rare(rare: pd.DataFrame, prefix: Path) -> None:
    cols = ("scvi_matched_f1_sd_across_subsets", "totalvi_f1_sd_across_subsets")
    # A single subsample gives no between-subset SD, so there is nothing to plot.
    if not any(rare[c].notna().any() for c in cols if c in rare):
        return
    x = rare["mean_subsample_source_count"]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for col, label, color in (
        (cols[0], "scVI_matched", LIGHT_BLUE),
        (cols[1], "totalVI", BLUE),
    ):
        ax.scatter(x, rare[col], color=color, label=label, s=36)
    if (x > 0).any():
        ax.set_xscale("log")
    ax.set_xlabel("Mean high-confidence source cells per class in the 3994-cell subsets (log)")
    ax.set_ylabel("F1 SD across the five source subsets")
    ax.set_title("Rare source classes drive transfer instability")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_neighbor(primary: pd.DataFrame, prefix: Path) -> None:
    sub = primary[
        (primary["classifier"] == "logreg") & (primary["evaluation_class_set"] == "common")
    ]
    ks = [15, 30, 50]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), sharey=True)
    for ax, level in zip(axes, ("l2", "l1")):
        part = sub[sub["label_level"] == level]
        x = np.arange(len(CONDITION_ORDER))
        for j, k in enumerate(ks):
            col = f"neighbor_agreement_k{k}"
            if col not in part:
                continue
            deltas = []
            for cond in CONDITION_ORDER:
                q = part[part["condition"] == cond]
                t = q[q["model"] == "totalvi"][col].mean()
                s = q[q["model"] == "scvi_matched"][col].mean()
                deltas.append(t - s)
            ax.bar(x + (j - 1) * 0.27, deltas, 0.27, label=f"k={k}", color=[LIGHT_BLUE, BLUE, GRAY][j])
        ax.axhline(0, color=GRAY, linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("_", "\n") for c in CONDITION_ORDER], fontsize=8)
        ax.set_title(f"{level} neighbor label agreement")
    axes[0].set_ylabel("Δ mean agreement (totalVI − scVI_matched)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    _save(fig, prefix)
