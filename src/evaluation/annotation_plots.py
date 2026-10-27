"""PHASE 6 figures. Neutral RNA UMAP is not an scVI/totalVI embedding."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save


def plot_reference_label_distribution(table: pd.DataFrame, output_prefix: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    order = table.sort_values("n_cells", ascending=False)
    ax.bar(np.arange(len(order)), order["n_cells"], color=BLUE)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order["cell_type"], rotation=60, ha="right")
    ax.set_ylabel("Reference cells")
    ax.set_title(title)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_query_umap_by_label(umap: np.ndarray, labels: pd.Series, output_prefix: Path, title: str) -> None:
    labels = labels.astype(str)
    types = sorted(labels.unique())
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for i, name in enumerate(types):
        mask = labels.to_numpy() == name
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=4,
            alpha=0.7,
            c=[cmap(i % 20)],
            label=name,
            linewidths=0,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=3, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_metric_bars(
    table: pd.DataFrame,
    metric: str,
    output_prefix: Path,
    title: str,
    models: list[str] | None = None,
) -> None:
    models = models or ["scVI_default", "scVI_matched", "totalVI"]
    palette = {
        "scVI_default": GRAY,
        "scVI_matched": LIGHT_BLUE,
        "PCA_RNA": "#9E9E9E",
        "MOFA+": "#5B8FA8",
        "totalVI": BLUE,
    }
    levels = [x for x in ("l1", "l2") if x in set(table["label_level"])]
    fig, axes = plt.subplots(1, len(levels), figsize=(1.35 * len(models) * max(len(levels), 1), 4.6), sharey=True)
    if len(levels) == 1:
        axes = [axes]
    for ax, level in zip(axes, levels):
        sub = table[(table["label_level"] == level) & (table["model"].isin(models))]
        for i, model in enumerate(models):
            row = sub[sub["model"] == model]
            if row.empty:
                continue
            ax.bar(i, float(row[metric].iloc[0]), color=palette.get(model, GRAY), label=model)
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels(models, rotation=20, ha="right")
        ax.set_title(f"{title} ({level})")
        ax.set_ylabel(metric)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_celltype_delta(table: pd.DataFrame, output_prefix: Path, title: str) -> None:
    order = table.sort_values("delta_totalVI_minus_scVI_matched")
    colors = [BLUE if v >= 0 else GRAY for v in order["delta_totalVI_minus_scVI_matched"]]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.bar(np.arange(len(order)), order["delta_totalVI_minus_scVI_matched"], color=colors)
    ax.axhline(0, color=GRAY, linewidth=1)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order["cell_type"], rotation=60, ha="right")
    ax.set_ylabel("Δ purity (totalVI − scVI_matched)")
    ax.set_title(title)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_restructuring_vs_delta(
    overlap: np.ndarray,
    delta_purity: np.ndarray,
    output_prefix: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.hexbin(overlap, delta_purity, gridsize=40, cmap="Blues", mincnt=1)
    ax.axhline(0, color=GRAY, linewidth=1)
    corr = float(np.corrcoef(overlap, delta_purity)[0, 1]) if len(overlap) > 2 else np.nan
    ax.set_xlabel("k=15 neighbor overlap (scVI_matched vs totalVI)")
    ax.set_ylabel("Δ purity (totalVI − scVI_matched)")
    ax.set_title(f"{title}\nPearson r = {corr:.3f}")
    fig.tight_layout()
    _save(fig, output_prefix)
