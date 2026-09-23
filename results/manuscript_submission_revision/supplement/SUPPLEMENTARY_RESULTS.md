# Supplementary Results

## SR1. Class-stratified paired bootstrap (sensitivity)

Primary intervals remain the **global** cell-level paired bootstrap (*B*=10,000).

| Quantity | Global 95% CI | Class-stratified 95% CI | Prop>0 (strat.) |
| --- | --- | --- | --- |
| \(G_{\mathrm{access}}\) | [0.043, 0.132] | [0.045, 0.119] | 1.000 |
| \(G_{\mathrm{assoc}}\) | [−0.012, 0.073] | [−0.008, 0.065] | 0.936 |
| \(Gap_{\mathrm{total}}\) | [0.084, 0.154] | [0.081, 0.142] | 1.000 |
| \(R_{\mathrm{access}}\) | [0.405, 1.114] | [0.448, 1.078] | — |

**Interpretation:** Stratification modestly narrows intervals but does **not** change conclusions: protein-access contrast remains entirely above zero; residual totalVI-associated contrast still crosses zero. Files: `tables/bootstrap_global_vs_class_stratified.csv`, `figures/bootstrap_global_vs_class_stratified.png`.

---

## SR2. Cell-type–specific post-adaptation degradation (Direction B, totalVI)

Mean per-class F1 (over seeds) at clean vs protein corruption 0.85:

| Cell type | F1 clean | F1 @0.85 | ΔF1 |
| --- | ---: | ---: | ---: |
| CD8 Naive | 0.784 | 0.519 | **−0.264** |
| CD8 TCM | 0.261 | 0.189 | −0.072 |
| CD8 TEM | 0.826 | 0.755 | −0.071 |
| CD4 Naive | 0.853 | 0.797 | −0.056 |
| CD4 TCM | 0.712 | 0.668 | −0.043 |
| CD4 TEM | 0.575 | 0.533 | −0.042 |

CD8 Naive shows the largest drop among focus subsets (≈ −0.26), consistent with surface-marker dependence for fine T-cell states. This is descriptive; no causal mechanism is claimed. Files: `tables/class_degradation_dirB_focus_subsets.csv`, `figures/class_degradation_dirB_*.png`.

---

## SR3. Selective prediction: accuracy vs macro-F1 (Direction B, totalVI)

Primary selective risk uses accuracy. At 50% vs 100% coverage (mean over 20 seeds):

| Coverage | Mean accuracy | Mean macro-F1 | Mean balanced accuracy |
| ---: | ---: | ---: | ---: |
| 1.00 | 0.787 | 0.654 | 0.659 |
| 0.50 | 0.966 | 0.703 | 0.698 |

Retaining low-uncertainty cells improves accuracy and also improves macro-F1 / balanced accuracy in this setting, but retention remains cell-type imbalanced (main text). File: `tables/selective_prediction_dirB_totalVI_cov50_vs_100.csv`.

---

## SR4. Lawlor donor-level summary (biological *n*=10)

Present-class macro-F1, totalVI − scVI_matched (seed-averaged within donor):

| Statistic | Value |
| --- | ---: |
| Mean paired Δ | +0.090 |
| Median paired Δ | +0.089 |
| SD | 0.010 |
| IQR | 0.013 |
| Min / Max | +0.077 / +0.110 |
| Donors favoring totalVI | 10/10 |
| Exact two-sided sign test *P* | 0.00195 |
| scVI_matched donor mean | 0.794 |
| totalVI donor mean | 0.884 |

Fold-mean simple baselines: protein PCA ≈ 0.927; concat ≈ 0.945.  
Computational consistency (latent AUROC): 38/50 fold×seed pairs — **not** biological *n*.

Files: `tables/lawlor_donor_summary_stats.csv`, `tables/lawlor_simple_baseline_fold_means.csv`.

---

## SR5. Marker-channel dropout (PCA baselines; totalVI not retrained)

Entire protein channels were zeroed before CLR/z/PCA on the development HC held-out set; RNA PCA was fit once on clean data. **totalVI was not retrained.** Clean concat macro-F1 matched the frozen baseline (0.745).

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

**Interpretation (cautious):** Protein-only PCA is sensitive to channel loss, especially T-lineage panel dropout. Concatenated RNA–protein PCA is more buffered by RNA PCs; single-channel drops can even slightly raise concat F1 via reweighting of remaining protein PCs. This sensitivity does **not** replace entry-wise training-time totalVI degradation and does not claim totalVI robustness under panel failure. File: `tables/marker_channel_dropout_pca_development.csv`.

---

## SR6. Blocked analyses (no fabricated results)

1. **totalVI retrain under marker-channel dropout** — GPU/time; TODO.  
2. **Additional multimodal comparator (WNN/MultiVI/Cobolt)** — not implemented in repo.  
3. **Second external cohort with non–protein-gated labels** — not available locally.
