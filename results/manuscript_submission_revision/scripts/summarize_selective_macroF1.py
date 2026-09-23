#!/usr/bin/env python3
"""Summarize selective-prediction macro-F1 / balanced accuracy vs coverage.

Primary selective metric remains risk = 1 - accuracy (AURC family).
This exports class-sensitive companions already present in phase11b tables.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results/phase11b_uncertainty_robustness/tables/phase11b_selective_prediction.csv"
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    df = pd.read_csv(SRC)
    # Aggregate over seeds: mean ± sd of macro_f1 and risk at each coverage
    g = (
        df.groupby(["direction", "model", "coverage"], as_index=False)
        .agg(
            n_seeds=("seed", "nunique"),
            mean_risk=("risk", "mean"),
            sd_risk=("risk", "std"),
            mean_accuracy=("accuracy", "mean"),
            mean_macro_f1=("macro_f1", "mean"),
            sd_macro_f1=("macro_f1", "std"),
            mean_balanced_accuracy=("balanced_accuracy", "mean"),
            sd_balanced_accuracy=("balanced_accuracy", "std"),
        )
    )
    g.to_csv(OUT / "selective_prediction_macroF1_by_coverage.csv", index=False)

    # Focus: coverage 1.0 vs 0.5 for totalVI Dir B
    focus = g[
        (g.direction == "direction_B")
        & (g.model == "totalvi")
        & (g.coverage.isin([0.5, 1.0]))
    ]
    focus.to_csv(OUT / "selective_prediction_dirB_totalVI_cov50_vs_100.csv", index=False)
    print(focus.to_string(index=False))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
