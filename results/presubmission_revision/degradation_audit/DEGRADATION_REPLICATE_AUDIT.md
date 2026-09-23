# DEGRADATION REPLICATE AUDIT

Source: `results/tables/protein_sparsity_primary_delta.csv` (historical; not rerun).

## Design verified from artifacts

- Corruption fractions: 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85
- At p=0.00: 1 row (perturbation_seed=0, model_seed=0)
- At each nonzero fraction: **5 corruption-mask replicates** (perturbation_seed 0–4)
- **model_seed is fixed at 0** for all rows (mask/retraining varies; model/training seed does not additionally vary in this table)
- scVI_matched reference is fixed across rows (~0.722)
- Each nonzero mask appears associated with its own totalVI retraining outcome (distinct macro-F1 per perturbation_seed)

## Mean ± SD (computational mask/retraining replicates)

 corruption_fraction  n_mask_replicates  n_unique_perturbation_seeds  n_unique_model_seeds model_seed_values  mean_macro_f1  sd_macro_f1  min_macro_f1  max_macro_f1  scVI_fixed_reference  manuscript_reported_mean
                0.00                  1                            1                     1                 0       0.779212          NaN      0.779212      0.779212              0.721864                     0.779
                0.10                  5                            5                     1                 0       0.771091     0.010491      0.756591      0.780642              0.721864                     0.771
                0.25                  5                            5                     1                 0       0.758287     0.011822      0.745848      0.777580              0.721864                     0.758
                0.40                  5                            5                     1                 0       0.759680     0.011684      0.745393      0.772723              0.721864                     0.760
                0.55                  5                            5                     1                 0       0.750314     0.016274      0.732916      0.774993              0.721864                     0.750
                0.70                  5                            5                     1                 0       0.735911     0.019897      0.715795      0.758571              0.721864                     0.736
                0.85                  5                            5                     1                 0       0.755452     0.010766      0.740903      0.766311              0.721864                     0.755

## Manuscript mean check

Reported manuscript means 0.779, 0.771, 0.758, 0.760, 0.750, 0.736, 0.755 match the historical summary means to rounding.

## Rebound 0.70 → 0.85

- p=0.70: mean=0.7359, SD=0.0199, min=0.7158, max=0.7586, values=[0.7158, 0.7238, 0.7253, 0.7561, 0.7586]
- p=0.85: mean=0.7555, SD=0.0108, min=0.7409, max=0.7663, values=[0.7409, 0.7511, 0.753, 0.7659, 0.7663]
- Range overlap between p=0.70 and p=0.85 replicate sets: True

Interpretation supported by overlap:

> The non-monotonic high-corruption means occurred within mask/retraining variability and should not be interpreted as evidence of improved performance under more severe corruption.

## Language for manuscript

- Call these **corruption-mask computational replicates** / **mask/retraining replicates**.
- Do **not** call them biological replicates.
- Primary figure: mean ± SD across the 5 masks at nonzero levels.
