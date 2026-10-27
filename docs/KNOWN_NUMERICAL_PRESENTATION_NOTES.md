# Known numerical presentation notes

These notes record differences in how values are *presented*. None of them is an analysis discrepancy, and no value was changed to remove them.

## 1. Direction A totalVI − scVI_matched: −0.003 (Figure 3) vs −0.004 (Table 3)

| Quantity | Value | Source |
| --- | --- | --- |
| scVI_matched, Direction A, 20-seed mean macro-F1 | 0.655623 | `results/transfer/matched_concat_control/tables/matched_transfer_control.csv` |
| totalVI, Direction A, 20-seed mean macro-F1 | 0.652383 | same |
| Exact difference (totalVI − scVI_matched) | −0.00324 | computed from the two frozen means |

- **Figure 3** rounds the exact difference: −0.00324 → **−0.003**.
- **Table 3** reports the rounded means, 0.652 and 0.656, and their difference 0.652 − 0.656 = **−0.004**.

The two presentations differ only in the order of rounding: rounding the difference versus differencing the rounded means. This is a **rounding-order difference, not an analysis discrepancy**. Both figures come from the same frozen means, and neither was changed.

## 2. PBMC 50% retention CV range

The manuscript reports the PBMC class-retention coefficient of variation at 50% coverage as approximately 0.88–1.02. The frozen values in `results/transfer/replicate_runs/tables/transfer_replicates_retention_imbalance.csv` at coverage 0.5 are 0.880 (Direction A scVI_matched), 0.947 (Direction A totalVI), 1.021 (Direction B scVI_matched) and 0.807 (Direction B totalVI). This note records the frozen values next to the reported range. Neither the range nor the values were changed.

## 3. Lawlor macro-F1 definition

Lawlor donor-level scores are **present-class** macro-F1: the average over classes with labelled support > 0 in that donor. This definition was declared before the donor deltas were computed, and fixed-ontology macro-F1 is the sensitivity analysis (`results/donor_validation/logs/lawlor_macro_f1_primary_definition.md`). Fold-level run scores in `results/donor_validation/lawlor_runs/tables/lawlor_primary_runs.csv` are scored per held-out fold, not per donor. They are therefore not numerically identical to the donor-level summary in `results/donor_validation/summary/`.
