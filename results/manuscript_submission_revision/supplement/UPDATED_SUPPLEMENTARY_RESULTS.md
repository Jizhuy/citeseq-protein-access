# UPDATED Supplementary Results

Companion to: *Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift*

Status legend for main-text Supplementary references: **COMPLETE** / **PARTIAL** / **MISSING**.

---

## SR1. Class-stratified paired bootstrap (sensitivity) — COMPLETE

Primary intervals remain the **global** cell-level paired bootstrap under the **strict train-only** protocol (*B*=10,000, seed 20260923).

| Quantity | Global 95% CI | Class-stratified 95% CI | Prop>0 (strat.) |
| --- | --- | --- | ---: |
| \(G_{\mathrm{access}}\) | [0.042, 0.131] | [0.046, 0.120] | 1.000 |
| \(G_{\mathrm{assoc}}\) | [−0.060, 0.030] | [−0.035, 0.023] | 0.358 |
| \(Gap_{\mathrm{total}}\) | [0.029, 0.130] | [0.043, 0.117] | 1.000 |
| \(R_{\mathrm{access}}\) | [0.680, 2.304] | [0.768, 1.968] | — |

**Interpretation:** Stratification modestly narrows intervals but does **not** change conclusions: protein-access contrast remains entirely above zero; residual totalVI-associated contrast remains near null and still crosses zero. Files: `strict_train_only/tables/strict_train_only_bootstrap_*.csv`.

### SR1b. Effect of held-out feature access on the internal deep-model comparison (not primary)

| Metric | Transductive internal (historical) | Strict train-only internal (**PRIMARY**) | Notes |
| --- | ---: | ---: | --- |
| RNA PCA | 0.666 | 0.667 | PCA was train-only in both |
| concat PCA | 0.745 | 0.745 | PCA was train-only in both |
| scVI_matched | 0.722 | 0.709 | Deep model access changed |
| totalVI | 0.779 | 0.741 | Deep model access changed |
| \(G_{\mathrm{access}}\) | +0.079 [0.043, 0.132] | +0.078 [0.042, 0.131] | Essentially unchanged |
| \(G_{\mathrm{assoc}}\) | +0.034 [−0.012, 0.073] | −0.004 [−0.060, 0.030] | Apparently attenuated |
| \(Gap_{\mathrm{total}}\) | +0.113 [0.084, 0.154] | +0.074 [0.029, 0.130] | Reduced with train-only deep models |
| \(R_{\mathrm{access}}\) | ≈0.70 [0.405, 1.114] | ≈1.05 [0.680, 2.304] | Secondary descriptive ratio only |

**Interpretation.** In the historical analysis, scVI_matched and totalVI were fit using features from all 10,849 cells while PCA remained train-only. Under the primary strict train-only protocol, all representations were fit on the same 7,573 HC training cells and evaluated on the same 1,921 held-out cells. The apparent totalVI advantage over concat was substantially attenuated when held-out-cell feature access was removed. Deep-model performance, particularly totalVI, was sensitive to the representation-learning information-access regime. This does **not** claim that transductive access alone causally explains the entire numerical difference as a controlled causal effect; it documents that information-access regime can materially affect apparent representation-level gains. Files: `STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md`; `strict_train_only/tables/`.


---

## SR2. Cell-type–specific post-adaptation degradation (Direction B, totalVI) — COMPLETE

Mean per-class F1 (over seeds) at clean vs protein corruption 0.85:

| Cell type | F1 clean | F1 @0.85 | ΔF1 |
| --- | ---: | ---: | ---: |
| CD8 Naive | 0.784 | 0.519 | **−0.264** |
| CD8 TCM | 0.261 | 0.189 | −0.072 |
| CD8 TEM | 0.826 | 0.755 | −0.071 |
| CD4 Naive | 0.853 | 0.797 | −0.056 |
| CD4 TCM | 0.712 | 0.668 | −0.043 |
| CD4 TEM | 0.575 | 0.533 | −0.042 |

CD8 Naive shows the largest drop among focus subsets (≈ −0.26). Descriptive only. Files: `tables/class_degradation_dirB_focus_subsets.csv`, `figures/class_degradation_dirB_*.png`.

---

## SR3. Selective prediction: accuracy vs macro-F1 / balanced accuracy — COMPLETE

Primary selective risk uses accuracy. At 50% vs 100% coverage (Direction B, totalVI; mean over 20 seeds):

| Coverage | Mean accuracy | Mean macro-F1 | Mean balanced accuracy |
| ---: | ---: | ---: | ---: |
| 1.00 | 0.787 | 0.654 | 0.659 |
| 0.50 | 0.966 | 0.703 | 0.698 |

Retention remains cell-type imbalanced (main text). Files: `tables/selective_prediction_dirB_totalVI_cov50_vs_100.csv`, `tables/selective_prediction_macroF1_by_coverage.csv`.

---

## SR4. Lawlor donor-level summaries — COMPLETE (updated to matched aggregation)

### Paired totalVI − scVI_matched (seed-averaged within donor; *n*=10)

| Statistic | Value |
| --- | ---: |
| Mean paired Δ | +0.090 |
| Median paired Δ | +0.089 |
| SD | 0.010 |
| IQR | 0.013 |
| Min / Max | +0.077 / +0.110 |
| Donors favoring totalVI | 10/10 |
| Exact two-sided sign test *P* | 0.00195 |

### Donor-mean present-class macro-F1 for all five representations (*n*=10)

| Representation | Mean | Median | SD | IQR | Min | Max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RNA PCA | 0.690 | 0.698 | 0.031 | 0.016 | 0.626 | 0.730 |
| Protein PCA | 0.923 | 0.920 | 0.013 | 0.017 | 0.898 | 0.940 |
| Concat PCA | 0.942 | 0.943 | 0.010 | 0.013 | 0.927 | 0.954 |
| scVI_matched | 0.794 | 0.787 | 0.016 | 0.030 | 0.772 | 0.815 |
| totalVI | 0.884 | 0.883 | 0.013 | 0.017 | 0.857 | 0.897 |

Files: `tables/lawlor_all_methods_donor_level_summary.csv`, `tables/lawlor_all_methods_per_donor.csv`, `tables/lawlor_donor_summary_stats.csv`.

### Historical fold-mean PCA aggregates (continuity only; not preferred for Figure 6)

Protein PCA fold mean ≈ 0.927; concat fold mean ≈ 0.945 (`tables/lawlor_simple_baseline_fold_means.csv`).

Computational consistency (latent AUROC): 38/50 fold×seed — **not** biological *n*.

Fixed-ontology Lawlor sensitivity: seed-mean tables include `fixed_ontology_macro_f1` identical to present-class for these donors in the audited seedmean file (no separate conflicting sensitivity table required for the reported classes).

---

## SR5. Marker-channel dropout (PCA baselines; totalVI not retrained) — COMPLETE

Clean concat macro-F1 matched the frozen baseline (0.745).

| Setting | Protein PCA | Concat | Δ concat vs clean |
| --- | ---: | ---: | ---: |
| Clean | 0.512 | 0.745 | 0.000 |
| Drop CD3 | 0.457 | 0.770 | +0.025 |
| Drop CD4 | 0.448 | 0.762 | +0.017 |
| Drop CD8a | 0.485 | 0.755 | +0.010 |
| Drop CD14 | 0.469 | 0.742 | −0.003 |
| Drop CD19 | 0.463 | 0.754 | +0.009 |
| Drop CD56 | 0.465 | 0.738 | −0.007 |
| Panel T-lineage (5 markers) | 0.328 | 0.720 | **−0.025** |
| Panel myeloid (3) | 0.421 | 0.750 | +0.005 |
| Panel B/NK (2) | 0.467 | 0.741 | −0.004 |

File: `tables/marker_channel_dropout_pca_development.csv`.

---

## SR6. Blocked analyses (no fabricated results)

1. **Strict train-only scVI_matched / totalVI retrain on HC train partition** — COMPLETE (PRIMARY); see `strict_train_only/` and `STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md`. Historical feature-transductive CASE B retained as SR1b sensitivity.
2. **totalVI retrain under marker-channel dropout** — not run.
3. **Additional multimodal comparator (WNN/MultiVI/Cobolt)** — not newly trained.
4. **Second external cohort with non–protein-gated labels** — not available locally.

---

## Supplementary completeness checklist vs main-text promises

| Claim | Status |
| --- | --- |
| Class-stratified paired-bootstrap sensitivity | COMPLETE |
| Cell-type-specific Direction B degradation | COMPLETE |
| PCA-only marker-channel dropout | COMPLETE |
| macro-F1 vs selective coverage | COMPLETE |
| balanced-accuracy vs selective coverage | COMPLETE |
| Lawlor donor median/SD/IQR/min/max + exact sign test | COMPLETE |
| Fixed-ontology Lawlor sensitivity | PARTIAL (column present; no separate divergent analysis) |
| Strict train-only deep internal retrain | MISSING (characterized as transductive; not claimed as completed) |
