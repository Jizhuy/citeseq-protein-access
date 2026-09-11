#!/usr/bin/env python3
"""PART 2: Lawlor pairwise consistency from canonical phase12b_primary_runs.csv."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "revision_round2_methodological_fixes"
STATS = OUT / "statistics"
LOGS = OUT / "logs"

PRIMARY = (
    ROOT
    / "results"
    / "phase12b_external_validation"
    / "formal_validation"
    / "tables"
    / "phase12b_primary_runs.csv"
)

# Higher-is-better vs lower-is-better
HIGHER_BETTER = {
    "macro_f1",
    "accuracy",
    "balanced_accuracy",
    "predictive_error_auroc",
    "latent_error_auroc",
    "predictive_error_auprc",
    "latent_error_auprc",
}
LOWER_BETTER = {"nll", "brier", "ece_equal_frequency", "ece_equal_width", "aurc"}

METRICS = [
    "macro_f1",
    "accuracy",
    "balanced_accuracy",
    "predictive_error_auroc",
    "latent_error_auroc",
    "nll",
    "brier",
    "ece_equal_frequency",
    "ece_equal_width",
    "aurc",
]


def pairwise_consistency(df: pd.DataFrame) -> pd.DataFrame:
    scvi = df[df["model"] == "scvi_matched"].set_index(["fold", "seed"])
    tot = df[df["model"] == "totalvi"].set_index(["fold", "seed"])
    common = scvi.index.intersection(tot.index)
    rows = []
    for metric in METRICS:
        a = scvi.loc[common, metric].astype(float)
        b = tot.loc[common, metric].astype(float)
        delta = b - a
        if metric in LOWER_BETTER:
            n_tot_better = int((delta < 0).sum())
            n_scvi_better = int((delta > 0).sum())
            # for reporting columns: totalVI > scVI means better performance
            n_tot_gt = n_tot_better
            n_tot_lt = n_scvi_better
            interpretation = "lower_better"
        else:
            n_tot_gt = int((delta > 0).sum())
            n_tot_lt = int((delta < 0).sum())
            interpretation = "higher_better"
        n_eq = int((delta == 0).sum())
        # Also raw numeric comparisons regardless of better direction
        n_raw_gt = int((delta > 0).sum())
        n_raw_lt = int((delta < 0).sum())
        rows.append(
            {
                "metric": metric,
                "n_pairs": int(len(common)),
                "interpretation": interpretation,
                "n_totalVI_better": n_tot_gt,
                "n_scVI_better": n_tot_lt,
                "n_equal": n_eq,
                "n_delta_positive_raw": n_raw_gt,
                "n_delta_negative_raw": n_raw_lt,
                "mean_delta_totalVI_minus_scVI": float(delta.mean()),
                "median_delta": float(delta.median()),
                "mean_scVI": float(a.mean()),
                "mean_totalVI": float(b.mean()),
            }
        )
    return pd.DataFrame(rows)


def audit_19_20_origin() -> str:
    """Document where 19/20 appears and that it is not Lawlor."""
    return """# Lawlor inconsistency audit: 50 vs 19/20

## Design reminder

Lawlor formal validation = 5 donor folds × 10 seeds × 2 models = **100 runs**.
Matched comparisons = **50** fold×seed pairs.

## Canonical counts (from `phase12b_primary_runs.csv`)

Computed in `statistics/lawlor_pairwise_consistency_corrected.csv`.

Headline historical report values (phase12b_report.md) already used fold×seed pairs:

| Metric | Historical report | Nature |
|---|---|---|
| macro-F1 | **50/50** totalVI higher | Lawlor |
| predictive AUROC | **45/50** | Lawlor |
| latent AUROC | **38/50** | Lawlor |
| NLL / Brier / AURC | **50/50** totalVI better | Lawlor |
| ECE | mixed (~29/50) | Lawlor |

## Where did **19/20** come from?

**Not from Lawlor.** It is an **internal PBMC transfer** sign-consistency count:

1. **Direction B macro-F1 (PBMC5k→PBMC10k), 20-seed grid:** totalVI favored in **19/20** seeds  
   — see `results/phase13_manuscript/manuscript/manuscript_v2_claim_checked.md`,  
   `results/phase14_submission_prep/manuscript/manuscript_genome_biology_v1.md`,  
   figure legends for cross-dataset transfer.

2. **Latent error AUROC in PBMC transfer:** also written as **19/20 seeds per direction**  
   — same manuscript family (PHASE 11B transfer, not Lawlor).

3. Lawlor latent sign consistency was correctly reported elsewhere as **38/50**, not 19/20.

## Root cause of the apparent inconsistency

The **19/20** figure is an **internal-transfer seed-consistency** number that coexists in the same manuscripts as Lawlor **50/50** / **38/50**.  
If prose or a table ever juxtaposed 19/20 with Lawlor without naming the experiment, that is a **cross-experiment copy/context error**, not a miscount of the Lawlor 50 pairs.

No Lawlor retraining is required for this correction.

## Action

- Keep Lawlor pairwise denominators as **/50**.
- Keep PBMC 20-seed denominators as **/20**.
- Never mix them in a single sentence without naming the experiment.
"""


def main() -> None:
    STATS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PRIMARY)
    assert len(df) == 100, len(df)
    assert set(df["model"]) == {"scvi_matched", "totalvi"}
    out = pairwise_consistency(df)
    out_path = STATS / "lawlor_pairwise_consistency_corrected.csv"
    out.to_csv(out_path, index=False)
    audit_path = LOGS / "lawlor_inconsistency_audit.md"
    audit_path.write_text(audit_19_20_origin())
    print(out.to_string(index=False))
    print("Wrote", out_path)
    print("Wrote", audit_path)


if __name__ == "__main__":
    main()
