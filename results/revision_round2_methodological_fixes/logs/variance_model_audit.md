# Variance model audit (PART 4)

## A. Was the original variance decomposition mathematically wrong?

**No.** Code in `src/experiments/phase12a_stats.py` and `scripts/12a_analyze.py`
decomposes **scVI_matched and totalVI separately**, and additionally decomposes
the paired **Delta = totalVI − scVI** cube. Absence of a model index inside a
per-model fit is correct.

If the manuscript prose described “the variance decomposition” without saying
“within each model / for the paired difference,” that is an **ambiguity of
description**, not an invalid estimator.

## B. Corrected within-model variance model

Y_ijr^(a) = mu_a + S_i^(a) + G_j^(a) + (SG)_ij^(a) + eps_ijr^(a)

for each model a in {scVI, totalVI}, with i=source subset (5),
j=model seed (5), r=run replicate (2).

## C. Corrected Delta variance model

Delta_ijr = Y_totalVI,ijr - Y_scVI,ijr
         = mu_delta + S_i_delta + G_j_delta + (SG)_ij_delta + eps_ijr_delta

This is the preferred model-comparison variance analysis.

## D. REML vs MoM

Primary estimator remains **method-of-moments from expected mean squares**.

For this **balanced** design with equal replication, MoM EMS estimators are
**identical to REML** under the usual normal random-effects model (Searle,
Casella & McCulloch). A separate MixedLM REML fit with only 5x5 levels is
exploratory and often singular; we do not elevate unstable REML fits over MoM.

Raw negative MoM components are preserved; nonnegative truncated summaries are
also reported (as in the historical tables).

## Simulation check (MoM recovery)

```json
{
  "source_subset": {
    "true": 0.02,
    "mean_raw": 0.01922679792736951,
    "sd_raw": 0.01651878921514566,
    "mean_nonneg": 0.01929932122287509
  },
  "model_seed": {
    "true": 0.01,
    "mean_raw": 0.010436623920246104,
    "sd_raw": 0.009685165167597496,
    "mean_nonneg": 0.010698787683047621
  },
  "interaction": {
    "true": 0.005,
    "mean_raw": 0.005670415114096678,
    "sd_raw": 0.00826129010693891,
    "mean_nonneg": 0.0068393967720673165
  },
  "residual": {
    "true": 0.03,
    "mean_raw": 0.02932510504413358,
    "sd_raw": 0.00823920961400382,
    "mean_nonneg": 0.02932510504413358
  }
}
```

## Artifacts

- `statistics/variance_components_corrected.csv`
- `statistics/delta_variance_components_corrected.csv`
- Historical files under `results/phase12a_variance_replication/` are untouched.

## Exploratory caveat

With only 5 subset and 5 seed levels, variance components are
**exploratory/descriptive**. Do not overstate precision.
