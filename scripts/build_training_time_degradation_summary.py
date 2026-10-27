#!/usr/bin/env python
"""Summarise frozen training-time protein-degradation runs by corruption fraction.

Input : results/degradation/training_time/protein_sparsity_primary_delta.csv (frozen per-run table)
Output: results/degradation/training_time/training_time_degradation_summary.csv

Pure aggregation of existing rows (count, unique seeds, mean, SD with ddof=1, min, max);
no model, bootstrap or random process is run. `manuscript_reported_mean` is the mean
rounded to three decimals, as reported in the manuscript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEG = ROOT / "results/degradation/training_time"


def build(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for frac, sub in df.groupby("corruption_fraction", sort=True):
        f1 = sub["totalVI_l2_macro_f1"]
        rows.append({
            "corruption_fraction": frac,
            "n_mask_replicates": int(len(sub)),
            "n_unique_perturbation_seeds": int(sub["perturbation_seed"].nunique()),
            "n_unique_model_seeds": int(sub["model_seed"].nunique()),
            "model_seed_values": ",".join(str(s) for s in sorted(sub["model_seed"].unique())),
            "mean_macro_f1": f1.mean(),
            "sd_macro_f1": f1.std(ddof=1),
            "min_macro_f1": f1.min(),
            "max_macro_f1": f1.max(),
            "scVI_fixed_reference": sub["scVI_matched_l2_macro_f1"].iloc[0],
            "manuscript_reported_mean": round(f1.mean(), 3),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEG / "protein_sparsity_primary_delta.csv")
    ap.add_argument("--output", type=Path, default=DEG / "training_time_degradation_summary.csv")
    args = ap.parse_args()
    df = pd.read_csv(args.input)
    assert df["scVI_matched_l2_macro_f1"].nunique() == 1, "scVI_matched reference is expected to be fixed"
    build(df).to_csv(args.output, index=False)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
