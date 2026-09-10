# PHASE 12A interpretation

Generated 2026-09-10T03:52:15.007655+00:00

## The question

PHASE 11B reported that roughly 80% of the variance in the crossed design fell
into a single "interaction/residual" term. That term was unidentified: with one
observation per grid cell, a real subset x seed interaction and pure run-level
noise produce the same number. This phase added a second independent replicate
per cell so the two can be told apart, and asks which of them it actually was.

## Metric-level components

```
       model      metric  source_subset_fraction  model_seed_fraction  interaction_fraction  residual_fraction  p_interaction  p_source_subset  p_model_seed  phase11b_interaction_residual_fraction
scvi_matched    macro_f1                  0.0231               0.0000                0.3545             0.6223         0.0428           0.3593        0.4790                                  1.0000
scvi_matched error_auroc                  0.0000               0.0845                0.2070             0.7085         0.1469           0.6804        0.1879                                  0.5889
scvi_matched error_auprc                  0.0114               0.0000                0.2839             0.7047         0.0899           0.3949        0.8346                                  0.9569
scvi_matched         nll                  0.1023               0.1681                0.2691             0.4605         0.0401           0.1393        0.0693                                  0.6809
scvi_matched       brier                  0.0703               0.0878                0.2727             0.5691         0.0640           0.2151        0.1806                                  0.6928
scvi_matched         ece                  0.0243               0.0000                0.0937             0.8820         0.3240           0.3384        0.4426                                  0.7326
scvi_matched        aurc                  0.0613               0.0977                0.3330             0.5079         0.0293           0.2428        0.1719                                  0.6396
     totalvi    macro_f1                  0.0000               0.0562                0.0324             0.9113         0.4271           0.6693        0.2287                                  0.7907
     totalvi error_auroc                  0.0000               0.2004                0.3017             0.4979         0.0365           0.4450        0.0603                                  0.7683
     totalvi error_auprc                  0.0020               0.0310                0.3878             0.5792         0.0276           0.4290        0.3378                                  0.9386
     totalvi         nll                  0.0000               0.0413                0.0000             0.9587         0.8300           0.4816        0.2027                                  0.8170
     totalvi       brier                  0.0000               0.0279                0.1149             0.8573         0.2894           0.7030        0.3275                                  0.9213
     totalvi         ece                  0.0246               0.0326                0.0000             0.9429         0.6313           0.3086        0.2757                                  0.7276
     totalvi        aurc                  0.0000               0.0402                0.1695             0.7902         0.2060           0.9331        0.2927                                  1.0000
```

## How the PHASE 11B lumped term split

```
       model      metric  phase11b_interaction_residual_fraction  interaction_fraction  residual_fraction  residual_share_of_lumped
scvi_matched    macro_f1                                  1.0000                0.3545             0.6223                    0.6371
scvi_matched error_auroc                                  0.5889                0.2070             0.7085                    0.7739
scvi_matched error_auprc                                  0.9569                0.2839             0.7047                    0.7128
scvi_matched         nll                                  0.6809                0.2691             0.4605                    0.6312
scvi_matched       brier                                  0.6928                0.2727             0.5691                    0.6760
scvi_matched         ece                                  0.7326                0.0937             0.8820                    0.9039
scvi_matched        aurc                                  0.6396                0.3330             0.5079                    0.6040
     totalvi    macro_f1                                  0.7907                0.0324             0.9113                    0.9656
     totalvi error_auroc                                  0.7683                0.3017             0.4979                    0.6227
     totalvi error_auprc                                  0.9386                0.3878             0.5792                    0.5990
     totalvi         nll                                  0.8170                0.0000             0.9587                    1.0000
     totalvi       brier                                  0.9213                0.1149             0.8573                    0.8818
     totalvi         ece                                  0.7276                0.0000             0.9429                    1.0000
     totalvi        aurc                                  1.0000                0.1695             0.7902                    0.8233
```

The interaction F-test (MS_interaction / MS_residual, which only exists once
n > 1) is significant at the 0.05 level in 5 of 14 model-by-metric
combinations. The median share of the old lumped term attributable to pure
run-level residual is 0.743.

## Delta decomposition

```
                   model      metric  source_subset_fraction  model_seed_fraction  interaction_fraction  residual_fraction  p_interaction  p_source_subset  p_model_seed  phase11b_interaction_residual_fraction
delta_totalvi_minus_scvi    macro_f1                  0.0509               0.1204                0.1411             0.6875         0.2144           0.2420        0.1101                                  1.0000
delta_totalvi_minus_scvi error_auroc                  0.0000               0.1310                0.1289             0.7400         0.2446           0.6157        0.1021                                  0.6803
delta_totalvi_minus_scvi error_auprc                  0.1671               0.0100                0.1742             0.6487         0.1629           0.0699        0.3899                                  0.9233
delta_totalvi_minus_scvi         nll                  0.0686               0.1959                0.0064             0.7292         0.4717           0.1554        0.0272                                  0.6772
delta_totalvi_minus_scvi       brier                  0.0222               0.1135                0.2627             0.6016         0.0773           0.3499        0.1419                                  0.7505
delta_totalvi_minus_scvi         ece                  0.0304               0.0303                0.1209             0.8184         0.2733           0.3165        0.3168                                  0.7723
delta_totalvi_minus_scvi        aurc                  0.0000               0.0560                0.3581             0.5860         0.0356           0.4518        0.2695                                  0.8197
```

## Per-cell decomposition (median over target cells)

```
       model  source_subset_fraction  model_seed_fraction  interaction_fraction  residual_fraction
scvi_matched                  0.1011               0.0107                0.0000             0.7763
     totalvi                  0.0545               0.0204                0.0100             0.8078
```

## Reading rule

Component fractions are descriptive and, after truncation at zero, upward
biased for null factors: on pure noise the three null components collectively
absorb roughly 15% of the truncated total. The F-tests, not the fractions,
decide whether a component is real. Method-of-moments component estimates are unbiased on the raw scale but fluctuate around zero for a null factor. After truncation at zero a null factor therefore receives a strictly positive fraction, and with a single 5x5x2 table that upward bias is large: on pure noise the residual typically retains only about half of the truncated total. Fractions are descriptive; the F-tests, not the fractions, decide whether a component is real.

## Variance model

Balanced two-factor random-effects ANOVA with n replicates per cell (Searle, Casella & McCulloch, Variance Components, ch. 4). Method-of-moments estimators from expected mean squares. Interaction and pure residual are separately identified only because n > 1; with n = 1 (the PHASE 11B design) MS_E is undefined and the two collapse into one term.
