# Result manifest

For each principal result this file gives the frozen value, the file that holds it, the script that generated it, and where it appears in the manuscript. Paths are relative to the repository root. Values are reported as in the manuscript (three decimals unless noted). Full-precision values are in the cited files. No value was recomputed during migration.

Abbreviations:

- Dir A = PBMC10k → PBMC5k.
- Dir B = PBMC5k → PBMC10k.
- `P` = `results/primary/strict_train_only/tables/`.

## 1. Internal analysis — PRIMARY (strict train-only)

| Result | Frozen value | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- |
| RNA PCA macro-F1 | 0.667 | `P/strict_train_only_point_estimates.csv` (`F1_RNA_PCA`) | `scripts/run_strict_train_only_internal.py` | Table 1; Figure 2 |
| Protein PCA macro-F1 | 0.512 | same (`F1_protein_PCA`) | same | Table 1 |
| Concat RNA–protein PCA macro-F1 | 0.745 | same (`F1_concat`) | same | Table 1; Figure 2 |
| scVI_matched macro-F1 | 0.709 | same (`F1_scVI_matched`) | same | Table 1; Figure 2 |
| totalVI macro-F1 | 0.741 | same (`F1_totalVI`) | same | Table 1; Figure 2 |
| MOFA+ macro-F1 (historical secondary baseline) | 0.549 | `results/supplementary/controls/tables/simple_multimodal_development.csv` (`MOFA+`) | historical run; not re-fit train-only | Table 1 (historical secondary baseline) |
| G_access | +0.078 [0.042, 0.131]; prop>0 = 1.0 | `P/strict_train_only_bootstrap_summary.csv` (`G_access`) | `scripts/run_strict_train_only_internal.py` (B = 10,000; seed 20260923) | Table 2; Figure 2; Abstract |
| G_assoc | −0.004 [−0.060, 0.030]; prop>0 = 0.355 | same (`G_assoc`) | same | Table 2; Figure 2; Abstract |
| Gap_total | +0.074 [0.029, 0.130] | same (`Gap_total`) | same | Table 2; Figure 2 |
| R_access (secondary descriptive ratio) | ≈ 1.05 [0.680, 2.304] | same (`R_access`) | same | Table 2 (secondary) |
| Global vs class-stratified bootstrap | — | `P/strict_train_only_bootstrap_global_vs_stratified.csv`, `P/strict_train_only_bootstrap_class_stratified.csv` | same (seed 20260924) | Supplementary |
| Cell counts | 7,573 train / 1,921 held-out | `P/run_meta.json` | same | Methods |

## 2. Internal analysis — SECONDARY information-access sensitivity

| Result | Source file | Generating script | Manuscript |
| --- | --- | --- | --- |
| Feature-transductive contrasts and bootstrap | `results/supplementary/information_access_sensitivity/feature_transductive_bootstrap/core_contrast_bootstrap_summary.csv` | `scripts/run_feature_transductive_bootstrap.py` (inputs from `scripts/run_feature_transductive_models.py`) | Supplementary (information-access sensitivity) |
| Class-stratified bootstrap (feature-transductive) | `results/supplementary/information_access_sensitivity/class_stratified_bootstrap/` | `scripts/cache_feature_transductive_predictions.py`, `scripts/run_class_stratified_bootstrap.py` | Supplementary |
| Strict vs feature-transductive comparison | `results/supplementary/information_access_sensitivity/information_access_sensitivity_comparison.csv` | FROZEN DERIVED ARTIFACT — NO STANDALONE GENERATOR RETAINED (§9) | Supplementary |

## 3. Bidirectional transfer and feature-access-matched control

Deep-model values are 20-seed means. They are frozen from the transfer replicate runs and restated in the control table.

| Result | Dir A | Dir B | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- | --- |
| Inductive concat PCA | 0.735 | 0.734 | `results/transfer/matched_concat_control/tables/matched_transfer_control.csv` | `scripts/run_matched_concat_transfer.py` | Table 3; Figure 3 |
| Feature-access-matched transductive concat PCA | 0.636 | 0.767 | same | same | Table 3; Figure 3 |
| scVI_matched (scArches, 20 seeds) | 0.656 | 0.601 | same; `results/transfer/replicate_runs/tables/transfer_replicates_20seed_summary.csv` | `scripts/run_bidirectional_transfer.py`, `scripts/analyze_bidirectional_transfer.py` | Table 3; Figure 3 |
| totalVI (scArches, 20 seeds) | 0.652 | 0.654 | same | same | Table 3; Figure 3 |
| totalVI − scVI_matched | −0.004 (Table 3) / −0.003 (Figure 3) | +0.052 | derived from the rows above; per-seed deltas in `transfer_replicates_seed_delta_summary.csv` | same | Table 3; Figure 3. See `KNOWN_NUMERICAL_PRESENTATION_NOTES.md`. |
| totalVI − transductive concat | +0.016 | −0.113 | derived from `matched_transfer_control.csv` | — | Table 3; Figure 3 |

## 4. Protein degradation

| Result | Frozen value | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- |
| Training-time totalVI macro-F1 by masked fraction 0 / 0.10 / 0.25 / 0.40 / 0.55 / 0.70 / 0.85 | 0.779 (n = 1); 0.771 ± 0.010; 0.758 ± 0.012; 0.760 ± 0.012; 0.750 ± 0.016; 0.736 ± 0.020; 0.755 ± 0.011 (mean ± SD, 5 mask replicates) | `results/degradation/training_time/training_time_degradation_summary.csv` (replicates: `training_time_degradation_replicates.csv`; per-run: `protein_sparsity_primary_delta.csv`) | `scripts/run_protein_degradation.py` (per-run table); `scripts/build_training_time_degradation_summary.py` (summary; EXACT GENERATOR ADDED, §9) | Table 4; Figure 4 |
| Post-adaptation Dir B totalVI macro-F1, clean → 85% corruption | 0.654 → 0.621 | `results/degradation/test_time_corruption/tables/test_time_corruption_summary.csv` | `scripts/run_test_time_protein_corruption.py` | Figure 4; Results text |
| Training-time vs test-time comparison | — | `results/degradation/test_time_corruption/tables/training_vs_testtime_corruption.csv` | same | Supplementary |
| Class-level degradation | — | `results/degradation/class_level/` | `scripts/summarize_class_degradation.py` | Supplementary |
| Marker-channel dropout (PCA only) | — | `results/supplementary/marker_channel_dropout/marker_channel_dropout_pca_development.csv` | `scripts/run_marker_channel_dropout.py` | Supplementary |

## 5. Uncertainty and selective prediction

| Result | Frozen value | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- |
| Predictive error AUROC (PBMC, 20 seeds) | ≈ 0.84–0.86 | `results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv` (`error_auroc_mean`) | `scripts/build_pbmc_reliability_summary.py` (EXACT GENERATOR ADDED, §9) | Table 5; Figure 5 |
| Latent error AUROC, Dir A scVI_matched / totalVI | 0.544 ± 0.054 / 0.669 ± 0.021 | means: same file (`latent_error_auroc_mean`); SDs: `results/transfer/replicate_runs/tables/transfer_replicates_20seed_summary.csv` (`error_auroc_latent`) | `scripts/analyze_bidirectional_transfer.py` | Table 5; Figure 5 |
| Latent error AUROC, Dir B scVI_matched / totalVI | 0.469 ± 0.075 / 0.610 ± 0.038 | same | same | Table 5; Figure 5 |
| Dir B totalVI predictive AUROC, clean → 85% corruption | 0.856 → 0.841 | `results/degradation/test_time_corruption/tables/test_time_corruption_summary.csv` (`error_auroc_mean`) | `scripts/run_test_time_protein_corruption.py` | Table 5 |
| Dir B totalVI full AURC | 0.060 → 0.079 | same (`full_aurc_mean`) | same | Table 5 |
| Dir B totalVI pAURC (coverage 0.5–1.0) | 0.112 → 0.144 | same (`paurc_0.5_1.0_mean`) | same | Table 5 |
| Dir B totalVI E-AURC | 0.0350 → 0.0448 | same (`e_aurc_mean`) | same | Table 5 |
| PBMC class-retention CV at 50% coverage | ≈ 0.88–1.02 | `results/transfer/replicate_runs/tables/transfer_replicates_retention_imbalance.csv` (coverage 0.5) | `scripts/analyze_bidirectional_transfer.py --stage retention` | Table 5 |
| Lawlor class-retention CV at 50% coverage | ≈ 0.57–0.59 | `results/donor_validation/lawlor_runs/tables/lawlor_replication_summary.csv` (selective-prediction row: "mean CV=0.578") | `scripts/analyze_lawlor_donor_validation.py` | Table 5 |
| Selective macro-F1 by coverage | — | `results/uncertainty/selective_prediction/` | `scripts/summarize_selective_macro_f1.py` | Supplementary |
| Corrected reliability / selective metrics; statistical-unit ledger | — | `results/uncertainty/calibration/full_reliability_metrics.csv`, `selective_prediction_corrected.csv`, `results/uncertainty/statistics/statistical_unit_ledger.csv` | `scripts/run_uncertainty_selective_prediction.py` | Supplementary |

## 6. Training stochasticity

Variance fractions for macro-F1, in the order source subset / model seed / subset × seed interaction / residual:

| Result | Frozen value | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- |
| scVI_matched | 2.31 / 0.00 / 35.45 / 62.23 % | `results/stochasticity/variance_decomposition/statistics/variance_components_corrected.csv` (`macro_f1`) | `scripts/run_variance_decomposition.py` (inputs: `scripts/run_training_stochasticity.py`, `scripts/analyze_training_stochasticity.py`) | Table 5; Figure 5 |
| totalVI | 0.00 / 5.62 / 3.24 / 91.13 % | same | same | Table 5; Figure 5 |
| Delta (totalVI − scVI_matched) | 5.09 / 12.04 / 14.11 / 68.75 % | `results/stochasticity/variance_decomposition/statistics/delta_variance_components_corrected.csv` (`macro_f1`) | same | Table 5; Figure 5 |
| Replicated crossed design (raw per-run scores) | — | `results/stochasticity/replicate_runs/tables/replicated_crossed_design.csv` | `scripts/run_training_stochasticity.py` | Supplementary |

## 7. Lawlor donor-level external validation (biological unit: 10 donors)

| Result | Frozen value | Source file | Generating script | Manuscript |
| --- | --- | --- | --- | --- |
| Donor-mean present-class macro-F1: RNA PCA / protein PCA / concat PCA / scVI_matched / totalVI | 0.690 / 0.923 / 0.942 / 0.794 / 0.884 | `results/donor_validation/summary/lawlor_all_methods_donor_level_summary.csv` | `scripts/run_lawlor_pca_per_donor.py`, `scripts/summarize_lawlor_donor_statistics.py` (deep-model runs: `scripts/run_lawlor_donor_validation.py`) | Table 6; Figure 6 |
| totalVI − scVI_matched mean donor delta | ≈ +0.090 | `results/donor_validation/summary/lawlor_donor_summary_stats.csv` (`mean_paired_delta`) | `scripts/summarize_lawlor_donor_statistics.py` | Table 6; Figure 6 |
| Donors favouring totalVI | 10/10 | same (`n_donors_favor_totalVI`) | same | Table 6; Figure 6 |
| Exact two-sided sign test | P = 0.00195 | same (`exact_sign_test_pvalue` = 0.001953125) | same | Table 6; Figure 6 |
| Latent-error-signal consistency (computational) | 38/50 fold × seed pairs | `results/donor_validation/statistics/lawlor_pairwise_consistency_corrected.csv` (`latent_error_auroc`) | `scripts/run_lawlor_pairwise_consistency.py` | Table 6 |
| Per-donor metrics and deltas | — | `results/donor_validation/statistics/lawlor_donor_delta.csv`, `lawlor_per_donor_metrics.csv` | `scripts/run_lawlor_donor_level.py` | Supplementary |
| Simple baselines on Lawlor | — | `results/supplementary/controls/tables/simple_multimodal_lawlor.csv`, `simple_multimodal_lawlor_per_donor.csv` | `scripts/run_simple_multimodal_baselines.py` | Supplementary |

## 8. Figures

| Manuscript figure | Final file | Vector source | Script | Data inputs |
| --- | --- | --- | --- | --- |
| Figure 1 | `figures/Figure_1.tiff` | `figures/source/Figure1_representation_framework.*` | `scripts/generate_figure1_representation_framework.py` | none (schematic) |
| Figure 2 | `figures/Figure_2.tiff` | `figures/source/Figure2_protein_access_vs_model.*` | `scripts/generate_publication_figures.py` | section 1 |
| Figure 3 | `figures/Figure_3.tiff` | `figures/source/Figure3_transfer_fairness.*` | same | section 3 |
| Figure 4 | `figures/Figure_4.tiff` | `figures/source/Figure4_protein_degradation.*` | same | section 4 |
| Figure 5 | `figures/Figure_5.tiff` | `figures/source/Figure5_uncertainty_stochasticity.*` | same | sections 5–6 |
| Figure 6 | `figures/Figure_6.tiff` | `figures/source/Figure6_lawlor_external.*` | same | section 7 |
| Supplementary Figure S1 | `figures/Supplementary_Figure_S1.*` | same | `scripts/generate_supplementary_figure_s1.py` | none (design schematic; illustrative icons) |

## 9. Provenance of derived artifacts without an original generator

Four outputs had no retained standalone generating script after migration. Each was audited against the frozen inputs and the existing project code. A helper was added only where the output could be reproduced exactly from frozen tables, without assumptions, models, bootstraps, random processes or downloads.

| Artifact | Outcome | Source input(s) | Generating script | Validation | Result |
| --- | --- | --- | --- | --- | --- |
| `results/degradation/training_time/training_time_degradation_summary.csv` | **EXACT GENERATOR ADDED** (exactly recoverable) | `results/degradation/training_time/protein_sparsity_primary_delta.csv` | `scripts/build_training_time_degradation_summary.py`: group by corruption fraction; count, unique seeds, mean, SD (ddof = 1), min, max; fixed scVI_matched reference; mean rounded to 3 decimals | Run in a temporary directory; `cmp` against the frozen file | **Byte-identical** |
| `results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv` | **EXACT GENERATOR ADDED** (exactly recoverable) | `results/transfer/replicate_runs/tables/transfer_replicates_20seed_primary.csv` | `scripts/build_pbmc_reliability_summary.py`: group by direction × model; seed means, SD (ddof = 1) for macro-F1 and full AURC | Run in a temporary directory; `cmp` against the frozen file | **Byte-identical** |
| `results/supplementary/information_access_sensitivity/information_access_sensitivity_comparison.csv` | **FROZEN DERIVED ARTIFACT — NO STANDALONE GENERATOR RETAINED** (partially recoverable) | Numbers trace to `results/primary/strict_train_only/tables/strict_train_only_point_estimates.csv`, `strict_train_only_bootstrap_summary.csv`, `.../feature_transductive_bootstrap/core_contrast_bootstrap_summary.csv` and `results/supplementary/controls/tables/simple_multimodal_development.csv` | none | Every number checked against those files | Not regenerated: approximate entries (for example "~0", "~-0.038"), interval formatting and role text were hand-formatted, so a generator would need invented formatting rules |
| Lawlor Matrix-Market → `lawlor_baseline_sng_shared.h5ad` conversion | **FROZEN DERIVED ARTIFACT — NO STANDALONE GENERATOR RETAINED** (partially recoverable) | HCA contributor matrices → `scripts/prepare_lawlor_baseline_matrices.R` (retained, original) → `baseline_mm/` | none for the Python step | Documented steps: Baseline + Demuxlet SNG subset; Ensembl GRCh37.87 → symbol mapping with sum aggregation of 55 duplicated symbols; 12 776 genes shared with the PBMC universe; 12-protein panel; lane batch key (`results/donor_validation/lawlor_runs/logs/data_prep_provenance.md`, `data_prep_summary.json`, `results/donor_validation/lawlor_dataset_audit/gene_mapping_external.csv`) | Not reconstructed: symbol tie-breaking, protein-name canonicalisation, `obs` field construction and cell ordering are not fully recorded, and the output could only be validated by re-processing raw data |

The two added helpers are reconstructions validated by exact byte-level agreement with the frozen outputs. They are not the original ad hoc commands. No frozen artifact was modified.
