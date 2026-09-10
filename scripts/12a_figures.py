#!/usr/bin/env python
"""PHASE 12A figures 1-5. Saves PNG and PDF.

Each figure is skipped with a notice if its input table is absent, so the
script can be run repeatedly as stages complete.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase12a_replication import ensure_dirs  # noqa: E402

MODEL_LABEL = {"scvi_matched": "scVI_matched", "totalvi": "totalVI",
               "delta_totalvi_minus_scvi": "$\\Delta$ totalVI $-$ scVI"}
COMPONENTS = ("source_subset", "model_seed", "interaction", "residual")
COMP_LABEL = {"source_subset": "source subset", "model_seed": "model seed",
              "interaction": "subset $\\times$ seed interaction",
              "residual": "pure residual (run-level)"}
COMP_COLOR = {"source_subset": "#E1892B", "model_seed": "#3B6FB6",
              "interaction": "#5C9E5C", "residual": "#9A9A9A"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
PATHS = ensure_dirs()
TAB, FIG, LOGS = PATHS["tables"], PATHS["figures"], PATHS["logs"]
METRIC_ORDER = ["macro_f1", "error_auroc", "brier", "ece", "aurc"]
METRIC_LABEL = {"macro_f1": "macro-F1", "error_auroc": "error AUROC",
                "error_auprc": "error AUPRC", "nll": "NLL", "brier": "Brier",
                "ece": "ECE", "aurc": "AURC"}


def save(fig, name: str) -> None:
    fig.savefig(FIG / f"{name}.png", dpi=300)
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}", flush=True)


def read(name: str) -> pd.DataFrame | None:
    path = TAB / name
    if not path.exists():
        print(f"  skip: {name} not found", flush=True)
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
def fig1_experimental_design() -> None:
    """Schematic of the replicated 5 x 5 x 2 crossed design."""
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    boxes = [
        (0.4, 4.2, 2.2, 1.2, "5 source subsets\n(PHASE 10B,\nn=3994 each)", "#E1892B"),
        (3.2, 4.2, 2.2, 1.2, "5 model seeds\n(initialisation\nheld fixed)", "#3B6FB6"),
        (6.0, 4.2, 2.2, 1.2, "2 run replicates\nrep0 = PHASE 11B\nrep1 = run_seed", "#5C9E5C"),
        (8.2, 4.2, 1.5, 1.2, "2 models\nscVI\ntotalVI", "#9A9A9A"),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, alpha=0.25,
                                   edgecolor=color, lw=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)

    ax.annotate("", xy=(3.1, 4.8), xytext=(2.7, 4.8),
                arrowprops=dict(arrowstyle="->", lw=1.5, color="black"))
    ax.annotate("", xy=(5.9, 4.8), xytext=(5.5, 4.8),
                arrowprops=dict(arrowstyle="->", lw=1.5, color="black"))
    ax.annotate("", xy=(8.1, 4.8), xytext=(8.3 - 0.1, 4.8),
                arrowprops=dict(arrowstyle="->", lw=1.5, color="black"))
    ax.text(5.0, 5.7, "Direction A only: PBMC10k $\\rightarrow$ PBMC5k",
            ha="center", fontsize=11, fontweight="bold")
    ax.text(5.0, 3.6,
            "100 observations = 5 $\\times$ 5 $\\times$ 2 $\\times$ 2\n"
            "With n=2 replicates, subset$\\times$seed interaction is separated from pure run-level residual.",
            ha="center", fontsize=9.5)

    ax.add_patch(plt.Rectangle((0.6, 0.5), 4.0, 2.6, facecolor="#F7F7F7",
                               edgecolor="#555555", lw=1.2))
    ax.text(2.6, 2.8, "Within each grid cell", ha="center", fontsize=9, fontweight="bold")
    ax.text(2.6, 1.7,
            "Held fixed:\n• exact source subset\n• model_seed (init)\n"
            "• architecture / hparams\n• query direction / labels\n\n"
            "Redrawn via run_replicate_seed:\n• train/val split, shuffle\n"
            "• dropout / adapter path",
            ha="center", va="center", fontsize=8)

    ax.add_patch(plt.Rectangle((5.2, 0.5), 4.2, 2.6, facecolor="#F7F7F7",
                               edgecolor="#555555", lw=1.2))
    ax.text(7.3, 2.8, "Identifiability", ha="center", fontsize=9, fontweight="bold")
    ax.text(7.3, 1.7,
            "PHASE 11B (n=1):\ninteraction and residual\nconfounded (~80% lumped)\n\n"
            "PHASE 12A (n=2):\n$\\hat\\sigma^2_{AB}$ and $\\hat\\sigma^2_E$\nseparately estimated\n"
            "(balanced MoM / EMS)",
            ha="center", va="center", fontsize=8)

    fig.suptitle("Figure 1. Replicated 5$\\times$5$\\times$2 experimental design",
                 fontweight="bold")
    save(fig, "figure01_replicated_experimental_design")


def fig1b_old_vs_new() -> None:
    """The headline: PHASE 11B's single lumped term split into two."""
    v = read("replicated_variance_components.csv")
    if v is None:
        return
    v = v[v.metric.isin(METRIC_ORDER)]
    models = ["scvi_matched", "totalvi"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, model in zip(axes, models):
        sub = v[v.model == model].set_index("metric").reindex(METRIC_ORDER)
        x = np.arange(len(METRIC_ORDER))
        w = 0.38
        # left bar: PHASE 11B (three terms, interaction and residual fused)
        b = np.zeros(len(x))
        for comp, col in (("phase11b_source_fraction", COMP_COLOR["source_subset"]),
                          ("phase11b_seed_fraction", COMP_COLOR["model_seed"])):
            ax.bar(x - w / 2, sub[comp], w, bottom=b, color=col, edgecolor="white", linewidth=0.4)
            b += sub[comp].to_numpy()
        ax.bar(x - w / 2, sub["phase11b_interaction_residual_fraction"], w, bottom=b,
               color="#5C9E5C", edgecolor="white", linewidth=0.4, hatch="///")
        # right bar: PHASE 12A (four terms, interaction and residual separated)
        b = np.zeros(len(x))
        for comp in COMPONENTS:
            ax.bar(x + w / 2, sub[f"{comp}_fraction"], w, bottom=b,
                   color=COMP_COLOR[comp], edgecolor="white", linewidth=0.4)
            b += sub[f"{comp}_fraction"].to_numpy()
        ax.set_xticks(x)
        ax.set_xticklabels([METRIC_LABEL[m] for m in METRIC_ORDER], rotation=20, ha="right")
        ax.set_title(f"{MODEL_LABEL[model]}\nleft: PHASE 11B (n=1)   right: PHASE 12A (n=2)")
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("fraction of total variance")
    handles = [plt.Rectangle((0, 0), 1, 1, color=COMP_COLOR[c]) for c in COMPONENTS]
    handles.append(plt.Rectangle((0, 0), 1, 1, color="#5C9E5C", hatch="///"))
    labels = [COMP_LABEL[c] for c in COMPONENTS] + ["PHASE 11B: interaction + residual (confounded)"]
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.16))
    fig.suptitle("Figure 1b. Replication splits the PHASE 11B interaction/residual term",
                 fontweight="bold")
    save(fig, "figure01b_old_vs_new_components")


def _component_panel(ax, sub, title):
    x = np.arange(len(METRIC_ORDER))
    b = np.zeros(len(x))
    for comp in COMPONENTS:
        vals = sub[f"{comp}_fraction"].to_numpy()
        ax.bar(x, vals, 0.62, bottom=b, color=COMP_COLOR[comp],
               edgecolor="white", linewidth=0.4, label=COMP_LABEL[comp])
        b += vals
    # mark components the F-test cannot distinguish from zero
    for i, m in enumerate(METRIC_ORDER):
        p = sub.loc[m, "p_interaction"] if "p_interaction" in sub else np.nan
        if np.isfinite(p) and p < 0.05:
            ax.text(i, 1.02, "*", ha="center", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([METRIC_LABEL[m] for m in METRIC_ORDER], rotation=20, ha="right")
    ax.set_ylim(0, 1.10)
    ax.set_title(title)


def fig2_by_model() -> None:
    v = read("replicated_variance_components.csv")
    if v is None:
        return
    v = v[v.metric.isin(METRIC_ORDER)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, model in zip(axes, ["scvi_matched", "totalvi"]):
        sub = v[v.model == model].set_index("metric").reindex(METRIC_ORDER)
        _component_panel(ax, sub, MODEL_LABEL[model])
    axes[0].set_ylabel("fraction of total variance")
    axes[1].legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    fig.suptitle("Figure 2. Replicated variance decomposition by model "
                 "(* interaction F-test p < 0.05)", fontweight="bold")
    save(fig, "figure02_variance_by_model")


def fig3_delta() -> None:
    d = read("replicated_delta_variance_components.csv")
    if d is None:
        return
    d = d[d.metric.isin(METRIC_ORDER)].set_index("metric").reindex(METRIC_ORDER)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    _component_panel(axes[0], d, "$\\Delta$ = totalVI $-$ scVI: variance components")
    axes[0].set_ylabel("fraction of total variance")
    axes[0].legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    ax = axes[1]
    x = np.arange(len(METRIC_ORDER))
    means = d["grand_mean"].to_numpy()
    sd = np.sqrt(d["total_variance_nonnegative"].to_numpy())
    ax.errorbar(x, means, yerr=sd, fmt="o", color="#D1495B", capsize=4, lw=1.4)
    ax.axhline(0, color="black", ls="--", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([METRIC_LABEL[m] for m in METRIC_ORDER], rotation=20, ha="right")
    ax.set_ylabel("$\\Delta$ (mean $\\pm$ total SD over 50 runs)")
    ax.set_title("$\\Delta$ magnitude vs its own variability")
    fig.suptitle("Figure 3. Decomposition of the totalVI $-$ scVI difference",
                 fontweight="bold")
    save(fig, "figure03_delta_variance_decomposition")


def fig4_celltype() -> None:
    c = read("celltype_replicated_variance.csv")
    if c is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
    for ax, model in zip(axes, ["scvi_matched", "totalvi"]):
        sub = c[c.model == model].sort_values("median_interaction_fraction", ascending=True)
        y = np.arange(len(sub))
        left = np.zeros(len(sub))
        for comp, key in (("source_subset", "median_source_fraction"),
                          ("model_seed", "median_seed_fraction"),
                          ("interaction", "median_interaction_fraction"),
                          ("residual", "median_residual_fraction")):
            vals = sub[key].to_numpy()
            ax.barh(y, vals, 0.72, left=left, color=COMP_COLOR[comp],
                    edgecolor="white", linewidth=0.4, label=COMP_LABEL[comp])
            left += vals
        ax.set_yticks(y)
        ax.set_yticklabels([f"{t} (n={n})" for t, n in
                            zip(sub.cell_type, sub.n_cells)], fontsize=7.5)
        ax.set_xlabel("median fraction of per-cell variance")
        ax.set_title(MODEL_LABEL[model])
        ax.set_xlim(0, 1)
    axes[1].legend(loc="lower right", fontsize=7.5, framealpha=0.95)
    fig.suptitle("Figure 4. Per-cell-type variance decomposition "
                 "(true-class probability, 5$\\times$5$\\times$2 design)", fontweight="bold")
    save(fig, "figure04_celltype_variance_decomposition")


def fig5_literature() -> None:
    """Conceptual comparison built from the audited method matrix."""
    path = TAB / "method_novelty_matrix.csv"
    if not path.exists():
        print("  skip: method_novelty_matrix.csv not found", flush=True)
        return
    m = pd.read_csv(path)
    prop_cols = [c for c in m.columns if c not in ("method", "citation", "year", "fusion_rule")]
    labels = {c: c.replace("_", " ") for c in prop_cols}
    mat = m[prop_cols].to_numpy()
    fig, ax = plt.subplots(figsize=(1.05 * len(prop_cols) + 4.2, 0.46 * len(m) + 2.4))
    # 1 = yes, 0 = no, -1 = not verified
    cmap = matplotlib.colors.ListedColormap(["#BEBEBE", "#EDEDED", "#3B6FB6"])
    norm = matplotlib.colors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(np.arange(len(prop_cols)))
    ax.set_xticklabels([labels[c] for c in prop_cols], rotation=42, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(m)))
    ax.set_yticklabels([f"{r.method} ({r.year})" for r in m.itertuples()], fontsize=8.5)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, {1: "yes", 0: "no", -1: "?"}.get(v, "?"),
                    ha="center", va="center", fontsize=7,
                    color="white" if v == 1 else "black")
    ax.set_xticks(np.arange(-0.5, len(prop_cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(m), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)
    fig.suptitle("Figure 5. Method novelty audit: what existing multimodal "
                 "integration methods already do", fontweight="bold")
    save(fig, "figure05_literature_method_comparison")


def main() -> int:
    for fn in (fig1_experimental_design, fig1b_old_vs_new, fig2_by_model,
               fig3_delta, fig4_celltype, fig5_literature):
        print(f"[{fn.__name__}]", flush=True)
        fn()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
