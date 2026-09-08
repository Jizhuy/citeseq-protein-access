"""PHASE 10 figures: cross-dataset transfer / query generalization."""

from __future__ import annotations

from pathlib import Path

import matplotlib

# The CUDA server is headless WSL2; the interactive Tk backend is unavailable.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save


def write_phase10_figures(tables_dir: Path, figures_dir: Path) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(tables_dir / "cross_dataset_primary.csv")
    delta = pd.read_csv(tables_dir / "cross_dataset_primary_delta.csv")
    cell = pd.read_csv(tables_dir / "cross_dataset_celltype_specific.csv")
    neigh = pd.read_csv(tables_dir / "cross_dataset_neighbor_agreement.csv")
    gap = pd.read_csv(tables_dir / "cross_dataset_generalization_gap.csv")
    confusion = pd.read_csv(tables_dir / "cross_dataset_confusion.csv")
    _fig_design(figures_dir / "phase10_figure1_design")
    _fig_f1(primary, figures_dir / "phase10_figure2_l2_macro_f1")
    _fig_delta(delta, figures_dir / "phase10_figure3_delta")
    _fig_l1_l2(primary, figures_dir / "phase10_figure4_l1_vs_l2")
    _fig_celltype(cell, figures_dir / "phase10_figure5_celltype_delta")
    _fig_confusion(confusion, figures_dir / "phase10_figure6_confusion")
    _fig_neighbor(neigh, figures_dir / "phase10_figure7_neighbor_agreement")
    _fig_gap(gap, figures_dir / "phase10_figure8_generalization_gap")


def _fig_design(prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 2.8))
    ax.set_axis_off()
    boxes = [
        (0.04, "SOURCE\ntrain only"),
        (0.28, "scVI_matched\nor totalVI"),
        (0.52, "scArches\nquery mapping"),
        (0.76, "TARGET\nevaluation only"),
    ]
    for x, text in boxes:
        ax.add_patch(plt.Rectangle((x, 0.28), 0.18, 0.5, fill=True, facecolor="#d6e6f2", edgecolor=BLUE, linewidth=1.2))
        ax.text(x + 0.09, 0.53, text, ha="center", va="center", fontsize=9, color=GRAY)
    for x in (0.22, 0.46, 0.70):
        ax.annotate("", xy=(x + 0.06, 0.53), xytext=(x, 0.53), arrowprops=dict(arrowstyle="->", color=BLUE))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("PHASE 10: cross-dataset transfer / query generalization")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_f1(primary: pd.DataFrame, prefix: Path) -> None:
    sub = primary[(primary["label_level"] == "l2") & (primary["classifier"] == "logreg")]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    dirs = list(sub["direction"].unique())
    x = np.arange(len(dirs))
    width = 0.35
    for i, (model, color) in enumerate((("scvi_matched", LIGHT_BLUE), ("totalvi", BLUE))):
        means, sds = [], []
        for d in dirs:
            vals = sub[(sub["direction"] == d) & (sub["model"] == model)]["macro_f1"]
            means.append(float(vals.mean()))
            sds.append(float(vals.std(ddof=1)) if len(vals) > 1 else 0.0)
        ax.bar(x + (i - 0.5) * width, means, width, yerr=sds, color=color, label=model, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([d.replace("_", "\n") for d in dirs])
    ax.set_ylabel("Cross-dataset l2 logreg macro-F1")
    ax.set_title("High-confidence l2 transfer performance")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_delta(delta: pd.DataFrame, prefix: Path) -> None:
    sub = delta[delta["label_level"] == "l2"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    dirs = list(sub["direction"].unique())
    for i, d in enumerate(dirs):
        part = sub[sub["direction"] == d]
        y = float(part["delta_totalvi_minus_scvi"].mean())
        lo = float(part["bootstrap_ci_low"].mean())
        hi = float(part["bootstrap_ci_high"].mean())
        ax.errorbar(i, y, yerr=[[y - lo], [hi - y]], fmt="o", color=BLUE, capsize=4)
        ax.scatter(np.full(len(part), i), part["delta_totalvi_minus_scvi"], color=LIGHT_BLUE, zorder=3)
    ax.axhline(0, color=GRAY, linewidth=0.8)
    ax.set_xticks(range(len(dirs)))
    ax.set_xticklabels([d.replace("_", "\n") for d in dirs])
    ax.set_ylabel("Δ macro-F1 (totalVI − scVI_matched)")
    ax.set_title("Paired bootstrap 95% CI on target cells")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_l1_l2(primary: pd.DataFrame, prefix: Path) -> None:
    sub = primary[primary["classifier"] == "logreg"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    styles = {"pbmc10k_to_pbmc5k": "-", "pbmc5k_to_pbmc10k": "--"}
    seen = set()
    for model, color, marker in (("scvi_matched", LIGHT_BLUE, "o"), ("totalvi", BLUE, "s")):
        for d in sub["direction"].unique():
            part = sub[(sub["model"] == model) & (sub["direction"] == d)]
            l1 = part[part["label_level"] == "l1"]["macro_f1"].mean()
            l2 = part[part["label_level"] == "l2"]["macro_f1"].mean()
            label = f"{model} {d.replace('_', ' ')}"
            if label in seen:
                continue
            seen.add(label)
            ax.plot([1, 2], [l1, l2], color=color, marker=marker, linestyle=styles.get(d, "-"), alpha=0.9, label=label)
    ax.legend(frameon=False, fontsize=8)
    ax.set_xticks([1, 2], ["l1", "l2"])
    ax.set_ylabel("Cross-dataset logreg macro-F1")
    ax.set_title("Coarse vs fine transfer")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_celltype(cell: pd.DataFrame, prefix: Path) -> None:
    l2 = cell[cell["label_level"] == "l2"]
    rows = []
    for (direction, seed, ctype), g in l2.groupby(["direction", "seed", "cell_type"]):
        scvi = g[g["model"] == "scvi_matched"]["f1"]
        tot = g[g["model"] == "totalvi"]["f1"]
        if len(scvi) and len(tot):
            rows.append({"direction": direction, "seed": seed, "cell_type": ctype, "delta": float(tot.mean() - scvi.mean())})
    if not rows:
        return
    frame = pd.DataFrame(rows)
    pivot = frame.groupby(["cell_type", "direction"])["delta"].mean().unstack(fill_value=np.nan)
    fig, ax = plt.subplots(figsize=(7.2, max(3.5, 0.28 * len(pivot))))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-0.3, vmax=0.3)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(list(pivot.columns), rotation=20, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(list(pivot.index), fontsize=7)
    fig.colorbar(im, ax=ax, label="Δ F1 (totalVI − scVI_matched)")
    ax.set_title("Cell-type-specific transfer Δ F1")
    fig.tight_layout()
    _save(fig, prefix)


def _fig_confusion(confusion: pd.DataFrame, prefix: Path) -> None:
    l2 = confusion[confusion["label_level"] == "l2"]
    dirs = list(l2["direction"].unique())
    models = ["scvi_matched", "totalvi"]
    if not dirs:
        return
    fig, axes = plt.subplots(len(dirs), 2, figsize=(9.5, 4.2 * len(dirs)))
    axes = np.atleast_2d(axes)
    for i, d in enumerate(dirs):
        for j, m in enumerate(models):
            ax = axes[i, j]
            sub = l2[(l2["direction"] == d) & (l2["model"] == m)]
            if sub.empty:
                continue
            mat = sub.pivot_table(index="true", columns="predicted", values="normalized", aggfunc="mean")
            im = ax.imshow(mat.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
            ax.set_title(f"{d}\n{m}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_neighbor(neigh: pd.DataFrame, prefix: Path) -> None:
    sub = neigh[neigh["label_level"] == "l2"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for model, color in (("scvi_matched", LIGHT_BLUE), ("totalvi", BLUE)):
        part = sub[sub["model"] == model].groupby("k", as_index=False)["mean_agreement"].mean()
        ax.plot(part["k"], part["mean_agreement"], marker="o", color=color, label=model)
    ax.set_xlabel("k SOURCE neighbors")
    ax.set_ylabel("Mean same-label fraction")
    ax.set_title("Cross-dataset neighborhood label agreement (l2)")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)


def _fig_gap(gap: pd.DataFrame, prefix: Path) -> None:
    sub = gap[gap["label_level"] == "l2"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    dirs = list(sub["direction"].unique())
    x = np.arange(len(dirs))
    width = 0.35
    for i, (model, color) in enumerate((("scvi_matched", LIGHT_BLUE), ("totalvi", BLUE))):
        vals = [sub[(sub["direction"] == d) & (sub["model"] == model)]["generalization_gap"].mean() for d in dirs]
        ax.bar(x + (i - 0.5) * width, vals, width, color=color, label=model)
    ax.axhline(0, color=GRAY, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([d.replace("_", "\n") for d in dirs])
    ax.set_ylabel("Cross F1 − PHASE 6 within F1")
    ax.set_title("Descriptive generalization gap (designs differ)")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, prefix)
