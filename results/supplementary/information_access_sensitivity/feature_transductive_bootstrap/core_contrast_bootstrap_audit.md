# Core contrast bootstrap audit

- n_boot: 10000
- bootstrap_seed: 20260923
- classifier_seed: 0
- n_test_cells: 1921
- n_train_cells: 7573
- concat_source: frozen_npy

## Reconstructed point estimates vs manuscript rounding

- RNA PCA F1: 0.666495 (target 0.666; match=True)
- concat F1: 0.744897 (target 0.745; match=True)
- totalVI F1: 0.778906 (target 0.779; match=True)
- G_access: 0.078402
- G_assoc: 0.034009
- Gap_total: 0.112411
- R_access: 0.697456

## Ratio denominator stability

- min Gap_total^(b): 5.639160e-02
- n Gap_total <= 0: 0
- n Gap_total < 0.01: 0
- n Gap_total < 0.05: 0

## Interval interpretation

Cell-level conditional paired-bootstrap 95% intervals.
NOT biological-generalization confidence intervals.

## Summary table

     quantity  point_estimate  boot_mean  boot_median  boot_sd    ci_2.5  ci_97.5  prop_gt_0  n_boot  n_nan  manuscript_rounded_point
   F1_RNA_PCA        0.666495   0.652613     0.653372 0.023607  0.604878 0.696657        NaN   10000      0                  0.666000
    F1_concat        0.744897   0.738889     0.737394 0.020307  0.702120 0.782772        NaN   10000      0                  0.745000
   F1_totalVI        0.778906   0.769836     0.770386 0.020816  0.728609 0.811766        NaN   10000      0                  0.779000
     G_access        0.078402   0.086276     0.085089 0.023828  0.042944 0.131534     1.0000   10000      0                  0.079000
      G_assoc        0.034009   0.030948     0.031864 0.022701 -0.011778 0.073462     0.9045   10000      0                  0.034000
    Gap_total        0.112411   0.117224     0.116556 0.018333  0.083896 0.153974     1.0000   10000      0                  0.113000
     R_access        0.697456   0.739975     0.722833 0.189399  0.404649 1.114481        NaN   10000      0                  0.699115
Gap_total_min        0.056392   0.117224     0.116556 0.018333  0.083896 0.153974     1.0000   10000      0                       NaN

## Full audit JSON

```json
{
  "n_hc": 9494,
  "n_train": 7573,
  "n_test": 1921,
  "n_classes_test": 30,
  "f1_rna": 0.6664948061673877,
  "f1_concat": 0.7448967230318437,
  "f1_totalVI": 0.7789060134068557,
  "target_rna": 0.666,
  "target_concat": 0.745,
  "target_totalVI": 0.779,
  "match_rna": true,
  "match_concat": true,
  "match_totalVI": true,
  "concat_source": "frozen_npy",
  "classifier": "LogisticRegression max_iter=2000 lbfgs",
  "seed": 0,
  "note": "Predictions reconstructed from frozen/train-only latents; VAEs not retrained."
}
```
