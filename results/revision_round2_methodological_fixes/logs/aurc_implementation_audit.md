# AURC implementation audit (PART 8)

## Historical implementation

`src/experiments/phase11_uncertainty.py::risk_coverage`:

- rank cells by increasing uncertainty
- risk_k = 1 - accuracy of the first k cells
- **AURC = mean(risk_k) for k=1..n**

This is the **full risk-coverage AURC over all prefixes** (coverage grid is NOT used
for the `aurc` scalar).

`COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)` is used only for discrete
selective tables (accuracy/risk/macro-F1 at those coverages).

## Correction required in manuscript language

Do **not** rename historical `aurc` to partial AURC — that would be incorrect.

Instead:

- Keep calling the stored scalar **full AURC**.
- Additionally report **pAURC_[0.5,1.0]** = mean risk over prefixes with coverage ≥ 0.5.
- Report **E-AURC** = AURC − AURC* where AURC* is the AURC of a perfect uncertainty
  ranking (all correct before all incorrect), following the selective-classification
  excess-AURC idea (Geifman & El-Yaniv, 2017; standard excess over optimal ranking).

## Artifacts

- `calibration/selective_prediction_corrected.csv` (Lawlor run-level)
- `calibration/full_reliability_metrics.csv`
