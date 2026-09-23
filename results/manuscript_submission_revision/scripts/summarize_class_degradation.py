#!/usr/bin/env python3
"""Summarize cell-type-specific post-adaptation protein degradation (Direction B).

Reads frozen per-class F1 table; does not retrain models.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = (
    ROOT
    / "results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_corruption_per_class.csv"
)
OUT = Path(__file__).resolve().parents[1] / "tables"
FIG = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

FOCUS = ["CD8 Naive", "CD8 TCM", "CD8 TEM", "CD4 TCM", "CD4 TEM", "CD4 Naive"]


def main():
    df = pd.read_csv(SRC)
    tv = df[(df["direction"] == "direction_B") & (df["model"] == "totalvi")].copy()
    # mean over seeds
    g = (
        tv.groupby(["corruption_fraction", "cell_type"], as_index=False)
        .agg(mean_f1=("f1", "mean"), sd_f1=("f1", "std"), n=("f1", "count"), mean_support=("support", "mean"))
    )
    clean = g[g.corruption_fraction == 0.0][["cell_type", "mean_f1"]].rename(
        columns={"mean_f1": "f1_clean"}
    )
    hi = g[g.corruption_fraction == 0.85][["cell_type", "mean_f1", "sd_f1", "mean_support"]].rename(
        columns={"mean_f1": "f1_p085", "sd_f1": "sd_p085"}
    )
    delta = clean.merge(hi, on="cell_type")
    delta["delta_f1"] = delta["f1_p085"] - delta["f1_clean"]
    delta = delta.sort_values("delta_f1")
    delta.to_csv(OUT / "class_degradation_dirB_totalVI_p085_delta.csv", index=False)

    focus = delta[delta.cell_type.isin(FOCUS)].copy()
    focus.to_csv(OUT / "class_degradation_dirB_focus_subsets.csv", index=False)

    meta = {
        "source": str(SRC.relative_to(ROOT)),
        "direction": "direction_B",
        "model": "totalvi",
        "comparison": "mean F1 over seeds at corruption 0.85 vs 0.00",
        "cd8_naive_delta": float(delta.loc[delta.cell_type == "CD8 Naive", "delta_f1"].iloc[0])
        if (delta.cell_type == "CD8 Naive").any()
        else None,
        "note": "Computational seed means; not biological replicates.",
    }
    (OUT / "class_degradation_meta.json").write_text(json.dumps(meta, indent=2))

    # Figure: focus subsets
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    plot_df = focus.sort_values("delta_f1")
    ax.barh(plot_df.cell_type, plot_df.delta_f1, color="#2F6FAE")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Δ macro-F1 (corruption 0.85 − clean)")
    ax.set_title("Direction B totalVI: cell-type sensitivity to target protein corruption")
    fig.tight_layout()
    fig.savefig(FIG / "class_degradation_dirB_focus.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Full trajectory for focus types
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for ct in FOCUS:
        sub = g[g.cell_type == ct].sort_values("corruption_fraction")
        if sub.empty:
            continue
        ax.plot(sub.corruption_fraction, sub.mean_f1, marker="o", label=ct)
    ax.set_xlabel("Protein corruption fraction")
    ax.set_ylabel("Mean per-class F1 (over seeds)")
    ax.set_title("Direction B totalVI per-class F1 under post-adaptation corruption")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "class_degradation_dirB_trajectories.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(focus.to_string(index=False))
    print("CD8 Naive delta:", meta["cd8_naive_delta"])


if __name__ == "__main__":
    main()
