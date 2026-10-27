"""PHASE 8 publication figures: protein sparsity / corruption stress test."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save


def _errorbar_by_level(ax, table: pd.DataFrame, y: str, yerr: str | None, color, label):
    levels = np.sort(table["corruption_fraction"].unique())
    means, lows, highs = [], [], []
    for level in levels:
        sub = table[table["corruption_fraction"] == level]
        means.append(float(sub[y].mean()))
        if yerr == "sd":
            sd = float(sub[y].std(ddof=1)) if len(sub) > 1 else 0.0
            lows.append(means[-1] - sd)
            highs.append(means[-1] + sd)
        elif yerr in sub.columns:
            lows.append(float(sub[yerr].mean()) if "ci_low" not in yerr else float(sub["bootstrap_ci_low"].iloc[0]))
            highs.append(float(sub["bootstrap_ci_high"].iloc[0]))
        else:
            lows.append(means[-1])
            highs.append(means[-1])
    ax.plot(levels, means, color=color, marker="o", linewidth=1.6, label=label)
    ax.fill_between(levels, lows, highs, color=color, alpha=0.18)


def plot_sparsity_vs_corruption(table: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    xcol = "corruption_fraction" if "corruption_fraction" in table.columns else "target_corruption_fraction"
    ax.plot(
        table[xcol],
        table["final_protein_sparsity"],
        color=BLUE,
        marker="o",
        linestyle="none",
        alpha=0.7,
    )
    grouped = table.groupby(xcol, as_index=False)["final_protein_sparsity"].mean()
    ax.plot(grouped[xcol], grouped["final_protein_sparsity"], color=BLUE, linewidth=1.6)
    ax.set_xlabel("Imposed nonzero-entry corruption fraction")
    ax.set_ylabel("Final protein matrix sparsity")
    ax.set_title("Protein sparsity vs corruption level")
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_primary_f1(table: pd.DataFrame, scvi_f1: float, pca_f1: float | None, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for seed, sub in table.groupby("perturbation_seed"):
        ax.plot(
            sub["corruption_fraction"],
            sub["logreg_macro_f1"],
            color=LIGHT_BLUE,
            marker="o",
            linewidth=1,
            alpha=0.55,
            label="_nolegend_",
        )
    grouped = table.groupby("corruption_fraction", as_index=False)["logreg_macro_f1"].agg(["mean", "std"]).reset_index()
    grouped["std"] = grouped["std"].fillna(0.0)
    ax.plot(grouped["corruption_fraction"], grouped["mean"], color=BLUE, marker="o", linewidth=1.8, label="totalVI (mean across masks)")
    ax.fill_between(
        grouped["corruption_fraction"],
        grouped["mean"] - grouped["std"],
        grouped["mean"] + grouped["std"],
        color=BLUE,
        alpha=0.15,
    )
    ax.axhline(scvi_f1, color=GRAY, linestyle="--", linewidth=1.4, label=f"scVI_matched ({scvi_f1:.3f})")
    if pca_f1 is not None:
        ax.axhline(pca_f1, color="#9E9E9E", linestyle=":", linewidth=1.4, label=f"PCA_RNA ({pca_f1:.3f})")
    ax.set_xlabel("Protein corruption fraction")
    ax.set_ylabel("Held-out logreg macro-F1 (celltype.l2)")
    ax.set_title("Primary robustness curve")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_delta_f1(delta: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    seed0 = delta[delta["perturbation_seed"] == 0].sort_values("corruption_fraction")
    if len(seed0):
        ax.plot(seed0["corruption_fraction"], seed0["delta"], color=BLUE, marker="o", linewidth=1.8, label="seed 0")
        ax.fill_between(
            seed0["corruption_fraction"],
            seed0["bootstrap_ci_low"],
            seed0["bootstrap_ci_high"],
            color=BLUE,
            alpha=0.18,
            label="seed 0 paired bootstrap 95% CI",
        )
    others = delta[delta["perturbation_seed"] != 0]
    if len(others):
        ax.scatter(
            others["corruption_fraction"],
            others["delta"],
            color=LIGHT_BLUE,
            s=22,
            alpha=0.7,
            label="other perturbation masks",
        )
    ax.axhline(0.0, color=GRAY, linewidth=1)
    ax.set_xlabel("Protein corruption fraction")
    ax.set_ylabel("Δ logreg macro-F1 (totalVI − scVI_matched)")
    ax.set_title("Multimodal advantage under protein corruption")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_asw_purity(table: pd.DataFrame, scvi: pd.Series, output_prefix: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), sharex=True)
    metrics = (("celltype_ASW", "Cell-type ASW"), ("knn_purity_k15", "kNN purity k=15"))
    for ax, (col, title) in zip(axes, metrics):
        grouped = table.groupby("corruption_fraction")[col].agg(["mean", "std"]).reset_index()
        grouped["std"] = grouped["std"].fillna(0.0)
        ax.plot(grouped["corruption_fraction"], grouped["mean"], color=BLUE, marker="o", linewidth=1.6, label="totalVI")
        ax.fill_between(
            grouped["corruption_fraction"],
            grouped["mean"] - grouped["std"],
            grouped["mean"] + grouped["std"],
            color=BLUE,
            alpha=0.15,
        )
        if col in scvi.index:
            ax.axhline(float(scvi[col]), color=GRAY, linestyle="--", label="scVI_matched")
        ax.set_title(title)
        ax.set_xlabel("Protein corruption fraction")
        ax.set_ylabel(col)
        ax.legend(frameon=False)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_masked_recovery(recovery: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    sub = recovery.dropna(subset=["pearson_masked"])
    if sub.empty:
        ax.text(0.5, 0.5, "No masked entries (0% corruption)", ha="center", va="center")
    else:
        grouped = sub.groupby("corruption_fraction")["pearson_masked"].mean().reset_index()
        ax.plot(grouped["corruption_fraction"], grouped["pearson_masked"], color=BLUE, marker="o", linewidth=1.6)
        ax.set_ylabel("Mean protein-wise Pearson (masked entries)")
    ax.set_xlabel("Protein corruption fraction")
    ax.set_title("Masked-entry protein recovery")
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_celltype_degradation(celltype: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    counts = celltype.groupby("cell_type")["support"].max()
    keep = counts[counts >= 10].index
    sub = celltype[celltype["cell_type"].isin(keep)]
    for name, part in sub.groupby("cell_type"):
        g = part.groupby("corruption_fraction")["f1"].mean()
        ax.plot(g.index, g.values, linewidth=1.2, marker="o", markersize=3, label=str(name))
    ax.set_xlabel("Protein corruption fraction")
    ax.set_ylabel("Held-out logreg F1")
    ax.set_title("Per-cell-type F1 vs protein corruption (l2, support ≥ 10)")
    ax.legend(frameon=False, fontsize=7, ncol=3)
    fig.tight_layout()
    _save(fig, output_prefix)


def plot_protein_recovery(recovery: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    sub = recovery.dropna(subset=["pearson_masked"])
    if sub.empty:
        ax.text(0.5, 0.5, "No masked entries", ha="center", va="center")
    else:
        for name, part in sub.groupby("protein"):
            g = part.groupby("corruption_fraction")["pearson_masked"].mean()
            ax.plot(g.index, g.values, linewidth=1.2, marker="o", markersize=3, label=str(name))
        ax.set_ylabel("Pearson (masked entries)")
        ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.set_xlabel("Protein corruption fraction")
    ax.set_title("Per-protein masked recovery")
    fig.tight_layout()
    _save(fig, output_prefix)
