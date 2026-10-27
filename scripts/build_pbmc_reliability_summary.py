#!/usr/bin/env python
"""Summarise frozen 20-seed PBMC transfer reliability metrics by direction and model.

Input : results/transfer/replicate_runs/tables/transfer_replicates_20seed_primary.csv (frozen per-seed table)
Output: results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv

Pure aggregation of existing per-seed rows (seed means; SD with ddof=1 for macro-F1 and
full AURC); no model, bootstrap or random process is run. The historical `aurc` column is
the full-coverage AURC.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results/transfer/replicate_runs/tables/transfer_replicates_20seed_primary.csv"
OUT = ROOT / "results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv"
NOTE = "paurc/E-AURC need cell-level probs (not present locally for 11B)"

MEANS = [
    ("accuracy_mean", "accuracy"),
    ("balanced_accuracy_mean", "balanced_accuracy"),
    ("error_prevalence_mean", "error_prevalence"),
    ("error_auroc_mean", "error_auroc_predictive"),
    ("error_auprc_mean", "error_auprc_predictive"),
    ("latent_error_auroc_mean", "error_auroc_latent"),
    ("nll_mean", "nll"),
    ("brier_mean", "brier"),
    ("ece_equal_frequency_mean", "ece_equal_frequency"),
    ("ece_equal_width_mean", "ece_equal_width"),
]


def build(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (direction, model), sub in df.groupby(["direction", "model"], sort=True):
        rec = {
            "table_block": "pbmc_20seed_summary",
            "experiment": f"PBMC_{direction}",
            "model": model,
            "n_seeds": int(len(sub)),
            "macro_f1_mean": sub["macro_f1"].mean(),
            "macro_f1_sd": sub["macro_f1"].std(ddof=1),
        }
        for new, col in MEANS:
            rec[new] = sub[col].mean()
        rec["full_aurc_mean"] = sub["aurc"].mean()
        rec["full_aurc_sd"] = sub["aurc"].std(ddof=1)
        rec["note"] = NOTE
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=IN)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    build(pd.read_csv(args.input)).to_csv(args.output, index=False)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
