#!/usr/bin/env python
"""Add across-seed variability to the scVI vs totalVI comparison table.

The paired target-cell bootstrap is conditional on a single model seed. For
several metrics the between-seed spread is far larger than that bootstrap, so
reporting the bootstrap CI alone would overstate the evidence. Both are kept,
side by side and separately labelled, as PHASE 11 §38 requires.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import ensure_dirs  # noqa: E402

METRIC_TO_COLUMN = {
    "error_auroc": "error_detection_auroc_pred",
    "error_auprc": "error_detection_auprc_pred",
    "ece": "ece_equal_frequency",
    "brier": "brier",
    "aurc": "aurc",
}


def main() -> int:
    paths = ensure_dirs()
    tables = paths["tables"]
    delta = pd.read_csv(tables / "uncertainty_scvi_vs_totalvi_delta.csv")
    primary = pd.read_csv(tables / "uncertainty_primary.csv")

    extra = []
    for _, row in delta.iterrows():
        col = METRIC_TO_COLUMN[row["metric"]]
        sub = primary[primary["direction"] == row["direction"]]
        piv = sub.pivot_table(index="model_seed", columns="model", values=col)
        d = piv["totalvi"] - piv["scvi_matched"]
        extra.append(
            {
                "across_seed_delta_mean": float(d.mean()),
                "across_seed_delta_sd": float(d.std(ddof=1)),
                "across_seed_delta_min": float(d.min()),
                "across_seed_delta_max": float(d.max()),
                "n_seeds": int(d.size),
                "seed0_delta": float(d.loc[0]),
                "bootstrap_conditional_on_seed": 0,
                "bootstrap_excludes_zero": bool(
                    (row["bootstrap_ci_low"] > 0) or (row["bootstrap_ci_high"] < 0)
                ),
                "across_seed_consistent_sign": bool((d > 0).all() or (d < 0).all()),
            }
        )
    out = pd.concat([delta, pd.DataFrame(extra)], axis=1)
    out["evidence"] = [
        "robust: bootstrap and all 5 seeds agree in sign"
        if r["across_seed_consistent_sign"] and r["bootstrap_excludes_zero"]
        else "seed-dependent: bootstrap CI excludes zero but seeds disagree in sign"
        if r["bootstrap_excludes_zero"]
        else "no reliable difference"
        for _, r in out.iterrows()
    ]
    out.to_csv(tables / "uncertainty_scvi_vs_totalvi_delta.csv", index=False)
    cols = [
        "direction", "metric", "delta_totalvi_minus_scvi", "bootstrap_ci_low",
        "bootstrap_ci_high", "across_seed_delta_mean", "across_seed_delta_sd", "evidence",
    ]
    print(out[cols].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
