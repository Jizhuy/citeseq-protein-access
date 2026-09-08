"""Matplotlib figures for scVI / totalVI baselines. No seaborn."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save

BATCH_COLORS = {"PBMC10k": BLUE, "PBMC5k": LIGHT_BLUE}


def plot_training_history(
    history: pd.DataFrame,
    output_prefix: Path,
    title: str = "Training history (scvi-tools logged quantities)",
) -> None:
    """Plot whichever metric columns scvi-tools actually logged."""
    exclude = {"epoch", "step"}
    metric_cols = [c for c in history.columns if c not in exclude]
    if not metric_cols:
        raise ValueError("Training history has no metric columns to plot.")
    x_col = "epoch" if "epoch" in history.columns else history.columns[0]
    n = len(metric_cols)
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 2.4 * n), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, col in zip(axes, metric_cols):
        ax.plot(history[x_col], history[col], color=BLUE, linewidth=1.5)
        ax.set_ylabel(col)
        ax.set_title(col)
    axes[-1].set_xlabel(x_col)
    fig.suptitle(title, color=BLUE, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, output_prefix)


def plot_umap_by_batch(
    umap: np.ndarray,
    batch: pd.Series,
    output_prefix: Path,
    title: str = "Latent UMAP colored by batch",
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    _scatter_umap_by_batch(ax, umap, batch)
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=3)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_side_by_side_umap_by_batch(
    umap_left: np.ndarray,
    umap_right: np.ndarray,
    batch: pd.Series,
    output_prefix: Path,
    title_left: str,
    title_right: str,
    suptitle: str,
) -> None:
    """Identical visual conventions; exploratory only."""
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), sharex=False, sharey=False)
    for ax, coords, title in (
        (axes[0], umap_left, title_left),
        (axes[1], umap_right, title_right),
    ):
        _scatter_umap_by_batch(ax, coords, batch)
        ax.set_title(title)
    axes[1].legend(frameon=False, markerscale=3)
    fig.suptitle(suptitle, color=GRAY, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, output_prefix)


def _scatter_umap_by_batch(ax, umap: np.ndarray, batch: pd.Series) -> None:
    batch = batch.astype(str)
    for name in sorted(batch.unique()):
        mask = batch.to_numpy() == name
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=4,
            alpha=0.55,
            c=BATCH_COLORS.get(name, GRAY),
            label=name,
            linewidths=0,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
