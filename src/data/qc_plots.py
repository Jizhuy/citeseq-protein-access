"""Publication-style matplotlib figures for data validation."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE = "#1f4e79"
LIGHT_BLUE = "#7ba3c9"
PALE_BLUE = "#d6e6f2"
GRAY = "#4a4a4a"

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": GRAY,
        "axes.labelcolor": GRAY,
        "xtick.color": GRAY,
        "ytick.color": GRAY,
        "text.color": GRAY,
        "axes.grid": False,
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_dataset_overview(
    adata_10k: ad.AnnData,
    adata_5k: ad.AnnData,
    combined: ad.AnnData,
    protein_panel: pd.DataFrame,
    summary: pd.DataFrame,
    output_prefix: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))

    counts = [
        adata_10k.n_obs,
        adata_5k.n_obs,
        combined.n_obs,
    ]
    labels = ["PBMC10k", "PBMC5k", "Combined"]
    axes[0, 0].bar(labels, counts, color=[BLUE, LIGHT_BLUE, PALE_BLUE], edgecolor=BLUE, width=0.65)
    axes[0, 0].set_ylabel("Number of cells")
    axes[0, 0].set_title("Cell counts")
    for x, y in zip(labels, counts):
        axes[0, 0].text(x, y, f"{y:,}", ha="center", va="bottom", fontsize=9)

    rna_s = summary.set_index("dataset")["rna_sparsity"]
    prot_s = summary.set_index("dataset")["protein_sparsity"]
    overview_datasets = [d for d in ["pbmc10k", "pbmc5k", "combined_inner"] if d in rna_s.index]
    x = np.arange(len(overview_datasets))
    width = 0.35
    axes[0, 1].bar(
        x - width / 2,
        [rna_s.loc[d] for d in overview_datasets],
        width=width,
        color=BLUE,
        label="RNA",
    )
    axes[0, 1].bar(
        x + width / 2,
        [prot_s.loc[d] for d in overview_datasets],
        width=width,
        color=LIGHT_BLUE,
        label="Protein",
    )
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(overview_datasets)
    axes[0, 1].set_ylabel("Fraction of zeros")
    axes[0, 1].set_ylim(0, 1.05)
    axes[0, 1].set_title("Modality sparsity")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].hist(
        np.log1p(adata_10k.obs["rna_counts"]),
        bins=40,
        histtype="stepfilled",
        color=BLUE,
        alpha=0.55,
        label="PBMC10k",
    )
    axes[1, 0].hist(
        np.log1p(adata_5k.obs["rna_counts"]),
        bins=40,
        histtype="stepfilled",
        color=LIGHT_BLUE,
        alpha=0.55,
        label="PBMC5k",
    )
    axes[1, 0].set_xlabel("log1p RNA library size")
    axes[1, 0].set_ylabel("Cells")
    axes[1, 0].set_title("RNA library size")
    axes[1, 0].legend(frameon=False)

    shared = int(protein_panel["shared"].sum())
    only_10k = int((protein_panel["in_pbmc10k"] & ~protein_panel["in_pbmc5k"]).sum())
    only_5k = int((protein_panel["in_pbmc5k"] & ~protein_panel["in_pbmc10k"]).sum())
    axes[1, 1].bar(
        ["Shared", "PBMC10k only", "PBMC5k only"],
        [shared, only_10k, only_5k],
        color=[BLUE, LIGHT_BLUE, PALE_BLUE],
        edgecolor=BLUE,
        width=0.65,
    )
    axes[1, 1].set_ylabel("Proteins")
    axes[1, 1].set_title("Protein panel overlap")
    for x_lab, y in zip(["Shared", "PBMC10k only", "PBMC5k only"], [shared, only_10k, only_5k]):
        axes[1, 1].text(x_lab, y, str(y), ha="center", va="bottom", fontsize=9)

    fig.suptitle("PBMC CITE-seq dataset overview", color=BLUE, fontsize=13, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, output_prefix)
