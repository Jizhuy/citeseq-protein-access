# Final PHASE 10 integrity check

**Date:** 2026-09-10  
**Scope:** supervisor manuscript + packaged tables/legends/figures

## Search targets
Obsolete mixed-coordinate markers: `-0.0102`, `0.6504`, uncorrected `tables/cross_dataset_primary*.csv` (non-`_corrected`), “sign flip” as primary interpretation.

## Results
| Check | Result |
|---|---|
| `-0.0102` in supervisor manuscript | **Absent** |
| `0.6504` in supervisor manuscript | **Absent** |
| Transfer narrative uses 20-seed near-null Dir A + positive Dir B | **Yes** |
| Mentions post-adaptation coordinate rule | **Yes** (Methods) |
| Packaged Table 2 uses PHASE 11B 20-seed + Lawlor primaries | **Yes** |
| Figure 3 assembly uses design + 11B/10B panels (not obsolete primary delta plot as sole source) | **Yes** (see figure QC) |

## Canonical sources to use (unchanged)
- `results/phase10_cross_dataset/tables_corrected/*`
- `results/phase10_cross_dataset/logs/phase10_erratum_latent_alignment.md`
- `results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv`

## Conclusion
**No PHASE 10 mixed-coordinate contamination found** in the supervisor package.
