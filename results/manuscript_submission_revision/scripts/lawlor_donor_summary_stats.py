#!/usr/bin/env python3
"""Donor-level Lawlor summary statistics (biological unit = donor, n=10).

Does not treat fold×seed as biological n.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[3]
DELTA = ROOT / "results/revision_round2_methodological_fixes/statistics/lawlor_donor_delta.csv"
LAWLOR_PCA = (
    ROOT / "results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_lawlor.csv"
)
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    d = pd.read_csv(DELTA)
    per = d[(d.metric == "present_class_macro_f1") & (d.donor != "SUMMARY")].copy()
    deltas = per["delta"].to_numpy(dtype=float)
    n_pos = int((deltas > 0).sum())
    n = len(deltas)
    # Exact two-sided sign test under H0: P(Δ>0)=0.5 (no ties here)
    sign = binomtest(n_pos, n=n, p=0.5, alternative="two-sided")

    summary = {
        "metric": "present_class_macro_f1",
        "n_donors": n,
        "mean_paired_delta": float(np.mean(deltas)),
        "median_paired_delta": float(np.median(deltas)),
        "sd_paired_delta": float(np.std(deltas, ddof=1)),
        "iqr_paired_delta": float(np.subtract(*np.percentile(deltas, [75, 25]))),
        "min_paired_delta": float(np.min(deltas)),
        "max_paired_delta": float(np.max(deltas)),
        "n_donors_favor_totalVI": n_pos,
        "n_donors_favor_scVI": int((deltas < 0).sum()),
        "exact_sign_test_pvalue": float(sign.pvalue),
        "scvi_mean": float(per["scvi_mean_over_seeds"].mean()),
        "totalvi_mean": float(per["totalvi_mean_over_seeds"].mean()),
        "note": "Biological unit = donor. Seed averaging is within-donor computational aggregation.",
    }
    pd.DataFrame([summary]).to_csv(OUT / "lawlor_donor_summary_stats.csv", index=False)
    per.to_csv(OUT / "lawlor_donor_paired_deltas_present_class.csv", index=False)

    pca = pd.read_csv(LAWLOR_PCA)
    fold_means = (
        pca.groupby("representation")["macro_f1"].agg(["mean", "std", "min", "max"]).reset_index()
    )
    fold_means.to_csv(OUT / "lawlor_simple_baseline_fold_means.csv", index=False)

    meta = {
        "sign_test": "exact binomial two-sided on sign of donor paired deltas",
        "summary": summary,
        "pca_fold_means": fold_means.to_dict(orient="records"),
    }
    (OUT / "lawlor_donor_summary_meta.json").write_text(json.dumps(meta, indent=2, default=float))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
