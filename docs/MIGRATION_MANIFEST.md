# Migration manifest

Source: development repository `Jizhuy/citeseq-protein-access` (local folder `multiomics_robustness`).
Target: publication repository `Jizhuy/citeseq-multimodal-reliability`.

Actions: COPIED UNCHANGED (byte-identical), RENAMED ONLY (byte-identical, new name), PATH NORMALIZED (code: paths, imports, file names and module names only; no algorithm, parameter or seed changed; frozen logs/JSON: development-repository path references in provenance text only, no value changed), SYNTAX FIXED (configuration file made parseable; no configured value changed), CONSOLIDATED (content merged into a new document), EXCLUDED, SUPERSEDED.

## Migrated files

| SOURCE | TARGET | ACTION | REASON |
| --- | --- | --- | --- |
| `scripts/01_prepare_data.py` | `scripts/prepare_data.py` | RENAMED ONLY | [A] Data preparation (PBMC10k/PBMC5k combined inner-join object). |
| `scripts/06_run_phase6.py` | `scripts/run_reference_annotation.py` | PATH NORMALIZED | [A] Reference-based (Seurat v4 / scArches) cell-type labels used by every analysis. |
| `scripts/02_run_baseline.py` | `scripts/run_feature_transductive_models.py` | PATH NORMALIZED | [A] Historical feature-transductive scVI / scVI_matched / totalVI fits on all development cells: reference model for training-time degradation and inputs to the secondary information-access sensitivity. |
| `results/manuscript_submission_revision/scripts/run_strict_train_only_internal.py` | `scripts/run_strict_train_only_internal.py` | PATH NORMALIZED | [A] PRIMARY strict train-only internal analysis and its bootstraps. |
| `results/revision_round2_methodological_fixes/scripts/revision2_simple_multimodal_baseline.py` | `scripts/run_simple_multimodal_baselines.py` | PATH NORMALIZED | [A] RNA/protein/concat PCA baselines (development, transfer, Lawlor). |
| `scripts/11b_run_grid.py` | `scripts/run_bidirectional_transfer.py` | PATH NORMALIZED | [A] 20-seed bidirectional transfer runs (scVI_matched, totalVI) and crossed design replicate 0. |
| `scripts/11b_analyze.py` | `scripts/analyze_bidirectional_transfer.py` | PATH NORMALIZED | [A] Scores transfer runs; 20-seed primary and selective-prediction tables. |
| `results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py` | `scripts/run_matched_concat_transfer.py` | PATH NORMALIZED | [A] Feature-access-matched concat PCA transfer control (Table 3 / Figure 3). |
| `scripts/08_run_sparsity_stress.py` | `scripts/run_protein_degradation.py` | PATH NORMALIZED | [A] Training-time protein degradation (Table 4 / Figure 4). |
| `results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py` | `scripts/run_test_time_protein_corruption.py` | PATH NORMALIZED | [A] Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/scripts/07_09_units_aurc_reliability.py` | `scripts/run_uncertainty_selective_prediction.py` | PATH NORMALIZED | [A] Statistical-unit ledger, selective prediction, reliability metrics (Figure 5). |
| `scripts/11b_crossed.py` | `scripts/analyze_crossed_design.py` | PATH NORMALIZED | [A] Scores the crossed source-subset x model-seed design (replicate 0). |
| `scripts/12a_run_replicates.py` | `scripts/run_training_stochasticity.py` | PATH NORMALIZED | [A] Trains replicate 1 of the crossed stochasticity design. |
| `scripts/12a_analyze.py` | `scripts/analyze_training_stochasticity.py` | PATH NORMALIZED | [A] Replicated variance components (Table 5). |
| `results/revision_round2_methodological_fixes/scripts/04_variance_decomposition_corrected.py` | `scripts/run_variance_decomposition.py` | PATH NORMALIZED | [A] Corrected variance decomposition used in Figure 5 / Table 5. |
| `scripts/12b_run_gate.py` | `scripts/run_lawlor_donor_validation.py` | PATH NORMALIZED | [A] Lawlor donor-held-out training runs (scVI_matched, totalVI). |
| `scripts/12b_analyze.py` | `scripts/analyze_lawlor_donor_validation.py` | PATH NORMALIZED | [A] Aggregates Lawlor runs. |
| `results/revision_round2_methodological_fixes/scripts/03_lawlor_donor_level.py` | `scripts/run_lawlor_donor_level.py` | PATH NORMALIZED | [A] Donor-level Lawlor metrics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/scripts/02_lawlor_pairwise_consistency.py` | `scripts/run_lawlor_pairwise_consistency.py` | PATH NORMALIZED | [A] Fold-by-seed computational consistency (38/50). |
| `results/manuscript_submission_revision/scripts/lawlor_pca_per_donor.py` | `scripts/run_lawlor_pca_per_donor.py` | PATH NORMALIZED | [A] Per-donor PCA baselines; all-method donor table (Figure 6). |
| `results/manuscript_submission_revision/scripts/lawlor_donor_summary_stats.py` | `scripts/summarize_lawlor_donor_statistics.py` | PATH NORMALIZED | [A] Donor paired-difference summary and exact sign test (Table 6). |
| `results/GPB_submission_FINAL_v2/scripts/generate_gpb_figures.py` | `scripts/generate_publication_figures.py` | PATH NORMALIZED | [A] Figures 2-6 from frozen tables (plotting only). |
| `scripts/generate_figure1_representation_framework.py` | `scripts/generate_figure1_representation_framework.py` | PATH NORMALIZED | [A] Figure 1 schematic (no data). |
| `scripts/print_environment.py` | `scripts/print_environment.py` | PATH NORMALIZED | [A] Records runtime package versions. |
| `results/manuscript_submission_revision/scripts/run_class_stratified_bootstrap.py` | `scripts/run_class_stratified_bootstrap.py` | PATH NORMALIZED | [B] Class-stratified bootstrap of the historical feature-transductive analysis (SR1b). |
| `results/presubmission_revision/bootstrap/presubmission_core_contrast_bootstrap.py` | `scripts/run_feature_transductive_bootstrap.py` | PATH NORMALIZED | [B] Historical feature-transductive core-contrast bootstrap (SR1b). |
| `results/manuscript_submission_revision/scripts/cache_development_predictions.py` | `scripts/cache_feature_transductive_predictions.py` | PATH NORMALIZED | [B] Caches feature-transductive held-out predictions for SR1b. |
| `results/manuscript_submission_revision/scripts/marker_channel_dropout_pca.py` | `scripts/run_marker_channel_dropout.py` | PATH NORMALIZED | [B] PCA-only marker-channel dropout (SR5). |
| `results/manuscript_submission_revision/scripts/summarize_class_degradation.py` | `scripts/summarize_class_degradation.py` | PATH NORMALIZED | [B] Per-class Direction B degradation (SR2). |
| `results/manuscript_submission_revision/scripts/summarize_selective_macroF1.py` | `scripts/summarize_selective_macro_f1.py` | PATH NORMALIZED | [B] Selective-prediction macro-F1 by coverage (SR3). |
| `scripts/generate_figure1_study_framework.py` | `scripts/generate_supplementary_figure_s1.py` | PATH NORMALIZED | [B] Supplementary Figure S1 study-design schematic (no data). |
| `results/phase12b_external_validation/raw_audit/convert_baseline_sparse.R` | `scripts/prepare_lawlor_baseline_matrices.R` | RENAMED ONLY | [B] Subsets Lawlor HCA matrices to Baseline Demuxlet singlets. |
| `results/revision_round2_methodological_fixes/scripts/revision2_p1_common.py` | `src/utils/analysis_common.py` | PATH NORMALIZED | Shared PCA/classifier helpers imported by the strict train-only, baseline and transfer-control scripts. |
| `src/__init__.py` | `src/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/experiments/__init__.py` | `src/experiments/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/experiments/phase8_sparsity.py` | `src/experiments/protein_sparsity_stress.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/utils/__init__.py` | `src/utils/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/utils/io.py` | `src/utils/io.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/utils/logging_utils.py` | `src/utils/logging_utils.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/utils/device.py` | `src/utils/device.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/utils/seed.py` | `src/utils/seed.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/__init__.py` | `src/evaluation/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/evaluation/phase6.py` | `src/evaluation/annotation_evaluation.py` | RENAMED ONLY | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/phase8_plots.py` | `src/evaluation/protein_sparsity_plots.py` | RENAMED ONLY | Imported (directly or transitively) by a migrated script. |
| `src/models/__init__.py` | `src/models/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/models/run_totalvi_corrupted.py` | `src/models/run_totalvi_corrupted.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/data/__init__.py` | `src/data/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/data/perturbations.py` | `src/data/perturbations.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/data/load_citeseq.py` | `src/data/load_citeseq.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/data/preprocess.py` | `src/data/preprocess.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/models/run_totalvi.py` | `src/models/run_totalvi.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/models/run_scvi.py` | `src/models/run_scvi.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/resource_metrics.py` | `src/evaluation/resource_metrics.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/representation_metrics.py` | `src/evaluation/representation_metrics.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/plots.py` | `src/evaluation/plots.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/data/qc_plots.py` | `src/data/qc_plots.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase12b_external.py` | `src/experiments/lawlor_external.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase10_cross_dataset.py` | `src/experiments/cross_dataset_transfer.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase11_uncertainty.py` | `src/experiments/uncertainty.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/models/run_scvi_matched.py` | `src/models/run_scvi_matched.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase11b_stats.py` | `src/experiments/transfer_replicate_stats.py` | RENAMED ONLY | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase11b_robustness.py` | `src/experiments/transfer_replicates.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/phase10_plots.py` | `src/evaluation/transfer_plots.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase10b_source_size.py` | `src/experiments/source_size_subsampling.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/phase10b_plots.py` | `src/evaluation/source_size_plots.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase6_run.py` | `src/experiments/reference_annotation.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/evaluation/phase6_plots.py` | `src/evaluation/annotation_plots.py` | RENAMED ONLY | Imported (directly or transitively) by a migrated script. |
| `src/annotation/__init__.py` | `src/annotation/__init__.py` | COPIED UNCHANGED | Package initialiser required by migrated modules. |
| `src/annotation/seurat_v4.py` | `src/annotation/seurat_v4.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/annotation/annotation_qc.py` | `src/annotation/annotation_qc.py` | COPIED UNCHANGED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase12a_stats.py` | `src/experiments/stochasticity_stats.py` | RENAMED ONLY | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase11b_analysis.py` | `src/experiments/transfer_replicate_analysis.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `src/experiments/phase12a_replication.py` | `src/experiments/stochasticity_replication.py` | PATH NORMALIZED | Imported (directly or transitively) by a migrated script. |
| `config/protein_gene_map.yaml` | `config/protein_gene_map.yaml` | SYNTAX FIXED | Protein-to-gene map used by the loader. The free-text `notes` item was quoted because its trailing colon made the source file invalid YAML (load_yaml would fail); no mapping entry changed. |
| `environment.yml` | `environment.yml` | COPIED UNCHANGED | Pinned environment (Python 3.11.7, scvi-tools 1.3.3, ...); versions preserved. |
| `requirements.txt` | `requirements.txt` | COPIED UNCHANGED | Full pinned package tree (includes mofapy2 0.7.5, statsmodels); versions preserved. |
| `config/config.yaml` | `config/config.yaml` | PATH NORMALIZED | Removed obsolete project stage/phase keys and the unused missing_modality block; renamed the query-adaptation key; no parameter value changed. |
| `data/annotations/README.md` | `data/annotations/README.md` | COPIED UNCHANGED | Reference-transfer cell labels (derived, required input). |
| `data/annotations/pbmc_reference_annotations.csv` | `data/annotations/pbmc_reference_annotations.csv` | COPIED UNCHANGED | Reference-transfer cell labels (derived, required input). |
| `data/annotations/pbmc_cell_annotations_template.csv` | `data/annotations/pbmc_cell_annotations_template.csv` | COPIED UNCHANGED | Reference-transfer cell labels (derived, required input). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_point_estimates.csv` | `results/primary/strict_train_only/tables/strict_train_only_point_estimates.csv` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_bootstrap_summary.csv` | `results/primary/strict_train_only/tables/strict_train_only_bootstrap_summary.csv` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_bootstrap_global_vs_stratified.csv` | `results/primary/strict_train_only/tables/strict_train_only_bootstrap_global_vs_stratified.csv` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_bootstrap_class_stratified.csv` | `results/primary/strict_train_only/tables/strict_train_only_bootstrap_class_stratified.csv` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_bootstrap_draws.csv` | `results/primary/strict_train_only/tables/strict_train_only_bootstrap_draws.csv` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/strict_train_only_heldout_predictions.npz` | `results/primary/strict_train_only/tables/strict_train_only_heldout_predictions.npz` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/strict_train_only/tables/run_meta.json` | `results/primary/strict_train_only/tables/run_meta.json` | COPIED UNCHANGED | PRIMARY strict train-only output (Tables 1-2, Figure 2). |
| `results/manuscript_submission_revision/STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md` | `results/primary/STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md` | COPIED UNCHANGED | Leakage/protocol audit of the primary analysis. |
| `results/revision_round3_information_vs_model/transfer_control/tables/matched_transfer_control.csv` | `results/transfer/matched_concat_control/tables/matched_transfer_control.csv` | COPIED UNCHANGED | Feature-access-matched transfer control (Table 3, Figure 3). |
| `results/revision_round3_information_vs_model/transfer_control/logs/matched_transfer_integrity.json` | `results/transfer/matched_concat_control/logs/matched_transfer_integrity.json` | COPIED UNCHANGED | Feature-access-matched transfer control (Table 3, Figure 3). |
| `results/revision_round3_information_vs_model/transfer_control/logs/matched_transfer_control_design.md` | `results/transfer/matched_concat_control/logs/matched_transfer_control_design.md` | COPIED UNCHANGED | Feature-access-matched transfer control (Table 3, Figure 3). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_20seed_primary.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_20seed_summary.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_20seed_summary.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_cell_bootstrap.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_cell_bootstrap.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_cell_instability.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_cell_instability.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_celltype_retention.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_celltype_retention.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_celltype_variance.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_celltype_variance.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_crossed_design.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_crossed_design.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_crossed_matrices.npz` | `results/transfer/replicate_runs/tables/transfer_replicates_crossed_matrices.npz` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_ensemble.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_ensemble.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_ensemble_convergence.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_ensemble_convergence.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_ensemble_percell.npz` | `results/transfer/replicate_runs/tables/transfer_replicates_ensemble_percell.npz` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_instability_quantiles.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_instability_quantiles.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_latent_drift.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_latent_drift.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_phase10_replication.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_cross_dataset_replication.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_phase10_replication_summary.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_cross_dataset_replication_summary.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_posterior_mc_check.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_posterior_mc_check.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_retention_imbalance.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_retention_imbalance.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_run_manifest_summary.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_run_manifest_summary.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_seed_delta.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_seed_delta.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_seed_delta_summary.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_seed_delta_summary.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_selective_prediction.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_selective_prediction.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_source_support_association.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_source_support_association.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_uncertainty_instability_assoc.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_uncertainty_instability_assoc.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/phase11b_uncertainty_robustness/tables/phase11b_variance_components.csv` | `results/transfer/replicate_runs/tables/transfer_replicates_variance_components.csv` | RENAMED ONLY | 20-seed transfer replicate tables (deep-model transfer values, uncertainty, crossed design replicate 0). |
| `results/revision_round2_methodological_fixes/controls/tables/rna_only_label_distribution.csv` | `results/supplementary/controls/tables/rna_only_label_distribution.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/rna_only_per_class_metrics.csv` | `results/supplementary/controls/tables/rna_only_per_class_metrics.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/rna_only_representation_benchmark.csv` | `results/supplementary/controls/tables/rna_only_representation_benchmark.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/rna_vs_multimodal_label_agreement.csv` | `results/supplementary/controls/tables/rna_vs_multimodal_label_agreement.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/rna_vs_multimodal_label_confusion.csv` | `results/supplementary/controls/tables/rna_vs_multimodal_label_confusion.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_development.csv` | `results/supplementary/controls/tables/simple_multimodal_development.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_lawlor.csv` | `results/supplementary/controls/tables/simple_multimodal_lawlor.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_lawlor_per_donor.csv` | `results/supplementary/controls/tables/simple_multimodal_lawlor_per_donor.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_per_class.csv` | `results/supplementary/controls/tables/simple_multimodal_per_class.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_transfer.csv` | `results/supplementary/controls/tables/simple_multimodal_transfer.csv` | COPIED UNCHANGED | Simple RNA/protein/concat PCA baselines and RNA-only label control. |
| `results/presubmission_revision/degradation_audit/training_time_degradation_summary.csv` | `results/degradation/training_time/training_time_degradation_summary.csv` | COPIED UNCHANGED | Training-time degradation summary (Table 4, Figure 4). |
| `results/presubmission_revision/degradation_audit/training_time_degradation_replicates.csv` | `results/degradation/training_time/training_time_degradation_replicates.csv` | COPIED UNCHANGED | Training-time degradation summary (Table 4, Figure 4). |
| `results/presubmission_revision/degradation_audit/DEGRADATION_REPLICATE_AUDIT.md` | `results/degradation/training_time/DEGRADATION_REPLICATE_AUDIT.md` | COPIED UNCHANGED | Training-time degradation summary (Table 4, Figure 4). |
| `results/tables/protein_sparsity_primary_delta.csv` | `results/degradation/training_time/protein_sparsity_primary_delta.csv` | COPIED UNCHANGED | Training-time protein-sparsity stress output underlying Table 4. |
| `results/tables/protein_sparsity_primary_delta_summary.csv` | `results/degradation/training_time/protein_sparsity_primary_delta_summary.csv` | COPIED UNCHANGED | Training-time protein-sparsity stress output underlying Table 4. |
| `results/tables/protein_sparsity_realized_levels.csv` | `results/degradation/training_time/protein_sparsity_realized_levels.csv` | COPIED UNCHANGED | Training-time protein-sparsity stress output underlying Table 4. |
| `results/tables/protein_sparsity_robustness.csv` | `results/degradation/training_time/protein_sparsity_robustness.csv` | COPIED UNCHANGED | Training-time protein-sparsity stress output underlying Table 4. |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/clean_baseline_vs_phase10.csv` | `results/degradation/test_time_corruption/tables/clean_baseline_vs_cross_dataset.csv` | RENAMED ONLY | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_corruption_all_runs.csv` | `results/degradation/test_time_corruption/tables/test_time_corruption_all_runs.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_corruption_delta_from_clean.csv` | `results/degradation/test_time_corruption/tables/test_time_corruption_delta_from_clean.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_corruption_per_class.csv` | `results/degradation/test_time_corruption/tables/test_time_corruption_per_class.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_corruption_summary.csv` | `results/degradation/test_time_corruption/tables/test_time_corruption_summary.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/test_time_reliability_metrics.csv` | `results/degradation/test_time_corruption/tables/test_time_reliability_metrics.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption/tables/training_vs_testtime_corruption.csv` | `results/degradation/test_time_corruption/tables/training_vs_testtime_corruption.csv` | COPIED UNCHANGED | Post-adaptation test-time protein corruption (Figure 4). |
| `results/manuscript_submission_revision/tables/class_degradation_dirB_focus_subsets.csv` | `results/degradation/class_level/class_degradation_dirB_focus_subsets.csv` | COPIED UNCHANGED | Per-class degradation (SR2). |
| `results/manuscript_submission_revision/tables/class_degradation_dirB_totalVI_p085_delta.csv` | `results/degradation/class_level/class_degradation_dirB_totalVI_p085_delta.csv` | COPIED UNCHANGED | Per-class degradation (SR2). |
| `results/manuscript_submission_revision/tables/class_degradation_meta.json` | `results/degradation/class_level/class_degradation_meta.json` | PATH NORMALIZED | Per-class degradation (SR2). Development-repository path references (provenance/log text only) replaced by repository-relative paths or a neutral description; no value changed. |
| `results/revision_round2_methodological_fixes/calibration/full_reliability_metrics.csv` | `results/uncertainty/calibration/full_reliability_metrics.csv` | COPIED UNCHANGED | Reliability / selective-prediction metrics (Figure 5). |
| `results/revision_round2_methodological_fixes/calibration/lawlor_reliability_summary_by_model.csv` | `results/uncertainty/calibration/lawlor_reliability_summary_by_model.csv` | COPIED UNCHANGED | Reliability / selective-prediction metrics (Figure 5). |
| `results/revision_round2_methodological_fixes/calibration/pbmc_20seed_reliability_summary.csv` | `results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv` | COPIED UNCHANGED | Reliability / selective-prediction metrics (Figure 5). |
| `results/revision_round2_methodological_fixes/calibration/selective_prediction_corrected.csv` | `results/uncertainty/calibration/selective_prediction_corrected.csv` | COPIED UNCHANGED | Reliability / selective-prediction metrics (Figure 5). |
| `results/revision_round2_methodological_fixes/statistics/statistical_unit_ledger.csv` | `results/uncertainty/statistics/statistical_unit_ledger.csv` | COPIED UNCHANGED | Statistical-unit ledger (computational vs biological units). |
| `results/manuscript_submission_revision/tables/selective_prediction_macroF1_by_coverage.csv` | `results/uncertainty/selective_prediction/selective_prediction_macroF1_by_coverage.csv` | COPIED UNCHANGED | Selective prediction macro-F1 (SR3). |
| `results/manuscript_submission_revision/tables/selective_prediction_dirB_totalVI_cov50_vs_100.csv` | `results/uncertainty/selective_prediction/selective_prediction_dirB_totalVI_cov50_vs_100.csv` | COPIED UNCHANGED | Selective prediction macro-F1 (SR3). |
| `results/revision_round2_methodological_fixes/statistics/delta_variance_components_corrected.csv` | `results/stochasticity/variance_decomposition/statistics/delta_variance_components_corrected.csv` | COPIED UNCHANGED | Variance decomposition (Table 5, Figure 5). |
| `results/revision_round2_methodological_fixes/statistics/variance_components_corrected.csv` | `results/stochasticity/variance_decomposition/statistics/variance_components_corrected.csv` | COPIED UNCHANGED | Variance decomposition (Table 5, Figure 5). |
| `results/revision_round2_methodological_fixes/logs/variance_model_audit.md` | `results/stochasticity/variance_decomposition/logs/variance_model_audit.md` | PATH NORMALIZED | Variance-model audit. Development-repository path references (provenance/log text only) replaced by repository-relative paths or a neutral description; no value changed. |
| `results/revision_round2_methodological_fixes/logs/variance_model_audit.json` | `results/stochasticity/variance_decomposition/logs/variance_model_audit.json` | COPIED UNCHANGED | Variance-model audit. |
| `results/phase12a_variance_replication/tables/celltype_replicated_variance.csv` | `results/stochasticity/replicate_runs/tables/celltype_replicated_variance.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/interaction_power_check.csv` | `results/stochasticity/replicate_runs/tables/interaction_power_check.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/method_novelty_matrix.csv` | `results/stochasticity/replicate_runs/tables/method_novelty_matrix.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/per_cell_replicated_variance.csv` | `results/stochasticity/replicate_runs/tables/per_cell_replicated_variance.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/replicated_crossed_design.csv` | `results/stochasticity/replicate_runs/tables/replicated_crossed_design.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/replicated_cubes.npz` | `results/stochasticity/replicate_runs/tables/replicated_cubes.npz` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/replicated_delta_variance_components.csv` | `results/stochasticity/replicate_runs/tables/replicated_delta_variance_components.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/phase12a_variance_replication/tables/replicated_variance_components.csv` | `results/stochasticity/replicate_runs/tables/replicated_variance_components.csv` | COPIED UNCHANGED | Replicated crossed-design tables (training stochasticity). |
| `results/revision_round2_methodological_fixes/statistics/lawlor_donor_delta.csv` | `results/donor_validation/statistics/lawlor_donor_delta.csv` | COPIED UNCHANGED | Lawlor donor-level statistics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/statistics/lawlor_donor_seedmean.csv` | `results/donor_validation/statistics/lawlor_donor_seedmean.csv` | COPIED UNCHANGED | Lawlor donor-level statistics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/statistics/lawlor_pairwise_consistency_corrected.csv` | `results/donor_validation/statistics/lawlor_pairwise_consistency_corrected.csv` | COPIED UNCHANGED | Lawlor donor-level statistics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/statistics/lawlor_per_donor_metrics.csv` | `results/donor_validation/statistics/lawlor_per_donor_metrics.csv` | COPIED UNCHANGED | Lawlor donor-level statistics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/statistics/lawlor_within_donor_seed_sd.csv` | `results/donor_validation/statistics/lawlor_within_donor_seed_sd.csv` | COPIED UNCHANGED | Lawlor donor-level statistics (biological unit n = 10). |
| `results/revision_round2_methodological_fixes/logs/lawlor_inconsistency_audit.md` | `results/donor_validation/logs/lawlor_inconsistency_audit.md` | PATH NORMALIZED | Lawlor metric-definition / consistency audit. Development-repository path references (provenance/log text only) replaced by repository-relative paths or a neutral description; no value changed. |
| `results/revision_round2_methodological_fixes/logs/lawlor_macro_f1_primary_definition.md` | `results/donor_validation/logs/lawlor_macro_f1_primary_definition.md` | COPIED UNCHANGED | Lawlor metric-definition / consistency audit. |
| `results/manuscript_submission_revision/tables/lawlor_all_methods_donor_level_summary.csv` | `results/donor_validation/summary/lawlor_all_methods_donor_level_summary.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_all_methods_per_donor.csv` | `results/donor_validation/summary/lawlor_all_methods_per_donor.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_donor_paired_deltas_present_class.csv` | `results/donor_validation/summary/lawlor_donor_paired_deltas_present_class.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_donor_summary_meta.json` | `results/donor_validation/summary/lawlor_donor_summary_meta.json` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_donor_summary_stats.csv` | `results/donor_validation/summary/lawlor_donor_summary_stats.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_pca_per_donor.csv` | `results/donor_validation/summary/lawlor_pca_per_donor.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/manuscript_submission_revision/tables/lawlor_simple_baseline_fold_means.csv` | `results/donor_validation/summary/lawlor_simple_baseline_fold_means.csv` | COPIED UNCHANGED | Lawlor donor summaries (Table 6, Figure 6). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_celltype_results.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_celltype_results.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_donor_folds.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_donor_folds.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_download_checksums.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_download_checksums.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_label_support_by_fold.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_label_support_by_fold.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_model_deltas.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_model_deltas.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_primary_protein_panel.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_primary_protein_panel.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_primary_runs.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_primary_runs.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_replication_summary.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_replication_summary.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/tables/phase12b_variability.csv` | `results/donor_validation/lawlor_runs/tables/lawlor_variability.csv` | RENAMED ONLY | Aggregated Lawlor run tables (per-run intermediates not migrated). |
| `results/phase12b_external_validation/formal_validation/logs/primary_class_set.json` | `results/donor_validation/lawlor_runs/logs/primary_class_set.json` | COPIED UNCHANGED | Lawlor preparation provenance, fold rule and class set. |
| `results/phase12b_external_validation/formal_validation/logs/data_prep_provenance.md` | `results/donor_validation/lawlor_runs/logs/data_prep_provenance.md` | COPIED UNCHANGED | Lawlor preparation provenance, fold rule and class set. |
| `results/phase12b_external_validation/formal_validation/logs/donor_fold_rule.md` | `results/donor_validation/lawlor_runs/logs/donor_fold_rule.md` | COPIED UNCHANGED | Lawlor preparation provenance, fold rule and class set. |
| `results/phase12b_external_validation/formal_validation/logs/data_prep_summary.json` | `results/donor_validation/lawlor_runs/logs/data_prep_summary.json` | PATH NORMALIZED | Lawlor preparation provenance, fold rule and class set. Development-repository path references (provenance/log text only) replaced by repository-relative paths or a neutral description; no value changed. |
| `results/phase12b_external_validation/formal_validation/logs/analysis_summary.json` | `results/donor_validation/lawlor_runs/logs/analysis_summary.json` | COPIED UNCHANGED | Lawlor preparation provenance, fold rule and class set. |
| `results/phase12b_external_validation/tables/external_donor_condition_counts.csv` | `results/donor_validation/lawlor_dataset_audit/external_donor_condition_counts.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_file_inventory.csv` | `results/donor_validation/lawlor_dataset_audit/external_file_inventory.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_gene_overlap.csv` | `results/donor_validation/lawlor_dataset_audit/external_gene_overlap.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_label_distribution.csv` | `results/donor_validation/lawlor_dataset_audit/external_label_distribution.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_label_distribution_by_condition.csv` | `results/donor_validation/lawlor_dataset_audit/external_label_distribution_by_condition.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_label_harmonization.csv` | `results/donor_validation/lawlor_dataset_audit/external_label_harmonization.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/external_protein_overlap.csv` | `results/donor_validation/lawlor_dataset_audit/external_protein_overlap.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/tables/gene_mapping_external.csv` | `results/donor_validation/lawlor_dataset_audit/gene_mapping_external.csv` | COPIED UNCHANGED | Lawlor dataset audit (donor/condition counts, panel and gene overlap). |
| `results/phase12b_external_validation/logs/data_provenance.md` | `results/donor_validation/lawlor_dataset_audit/data_provenance.md` | COPIED UNCHANGED | Lawlor provenance and independence audit. |
| `results/phase12b_external_validation/logs/independence_audit.md` | `results/donor_validation/lawlor_dataset_audit/independence_audit.md` | COPIED UNCHANGED | Lawlor provenance and independence audit. |
| `results/manuscript_submission_revision/tables/information_access_sensitivity_comparison.csv` | `results/supplementary/information_access_sensitivity/information_access_sensitivity_comparison.csv` | COPIED UNCHANGED | SECONDARY: historical feature-transductive vs strict train-only comparison (SR1b). |
| `results/presubmission_revision/bootstrap/core_contrast_bootstrap_summary.csv` | `results/supplementary/information_access_sensitivity/feature_transductive_bootstrap/core_contrast_bootstrap_summary.csv` | COPIED UNCHANGED | SECONDARY: historical feature-transductive bootstrap (SR1b). |
| `results/presubmission_revision/bootstrap/core_contrast_bootstrap_audit.md` | `results/supplementary/information_access_sensitivity/feature_transductive_bootstrap/core_contrast_bootstrap_audit.md` | COPIED UNCHANGED | SECONDARY: historical feature-transductive bootstrap (SR1b). |
| `results/manuscript_submission_revision/tables/bootstrap_global_vs_class_stratified.csv` | `results/supplementary/information_access_sensitivity/class_stratified_bootstrap/bootstrap_global_vs_class_stratified.csv` | COPIED UNCHANGED | SECONDARY: historical feature-transductive class-stratified bootstrap (SR1b). |
| `results/manuscript_submission_revision/tables/bootstrap_class_stratified_summary.csv` | `results/supplementary/information_access_sensitivity/class_stratified_bootstrap/bootstrap_class_stratified_summary.csv` | COPIED UNCHANGED | SECONDARY: historical feature-transductive class-stratified bootstrap (SR1b). |
| `results/manuscript_submission_revision/tables/bootstrap_class_stratified_meta.json` | `results/supplementary/information_access_sensitivity/class_stratified_bootstrap/bootstrap_class_stratified_meta.json` | PATH NORMALIZED | SECONDARY: historical feature-transductive class-stratified bootstrap (SR1b). Development-repository path references (provenance/log text only) replaced by repository-relative paths or a neutral description; no value changed. |
| `results/manuscript_submission_revision/tables/development_heldout_predictions_audit.json` | `results/supplementary/information_access_sensitivity/class_stratified_bootstrap/development_heldout_predictions_audit.json` | COPIED UNCHANGED | SECONDARY: historical feature-transductive class-stratified bootstrap (SR1b). |
| `results/manuscript_submission_revision/tables/marker_channel_dropout_pca_development.csv` | `results/supplementary/marker_channel_dropout/marker_channel_dropout_pca_development.csv` | COPIED UNCHANGED | PCA-only marker-channel dropout (SR5). |
| `results/manuscript_submission_revision/tables/marker_channel_dropout_meta.json` | `results/supplementary/marker_channel_dropout/marker_channel_dropout_meta.json` | COPIED UNCHANGED | PCA-only marker-channel dropout (SR5). |
| `results/manuscript_submission_revision/tables/runtime_environment.json` | `results/supplementary/environment/runtime_environment.json` | COPIED UNCHANGED | Recorded runtime package versions. |
| `results/GPB_submission_FINAL_v2/figures/Figure_1.tiff` | `figures/Figure_1.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures/Figure_2.tiff` | `figures/Figure_2.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures/Figure_3.tiff` | `figures/Figure_3.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures/Figure_4.tiff` | `figures/Figure_4.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures/Figure_5.tiff` | `figures/Figure_5.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures/Figure_6.tiff` | `figures/Figure_6.tiff` | COPIED UNCHANGED | Final main figure (600 dpi TIFF). |
| `results/GPB_submission_FINAL_v2/figures_source/Figure2_protein_access_vs_model.pdf` | `figures/source/Figure2_protein_access_vs_model.pdf` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure2_protein_access_vs_model.png` | `figures/source/Figure2_protein_access_vs_model.png` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure3_transfer_fairness.pdf` | `figures/source/Figure3_transfer_fairness.pdf` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure3_transfer_fairness.png` | `figures/source/Figure3_transfer_fairness.png` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure4_protein_degradation.pdf` | `figures/source/Figure4_protein_degradation.pdf` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure4_protein_degradation.png` | `figures/source/Figure4_protein_degradation.png` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure5_uncertainty_stochasticity.pdf` | `figures/source/Figure5_uncertainty_stochasticity.pdf` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure5_uncertainty_stochasticity.png` | `figures/source/Figure5_uncertainty_stochasticity.png` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure6_lawlor_external.pdf` | `figures/source/Figure6_lawlor_external.pdf` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/GPB_submission_FINAL_v2/figures_source/Figure6_lawlor_external.png` | `figures/source/Figure6_lawlor_external.png` | COPIED UNCHANGED | Vector source of main Figures 2-6. |
| `results/model_framework_revision/figures/Figure1_representation_framework.pdf` | `figures/source/Figure1_representation_framework.pdf` | COPIED UNCHANGED | Vector source of main Figure 1. |
| `results/model_framework_revision/supplementary/SupplementaryFigureS1_study_design.pdf` | `figures/Supplementary_Figure_S1.pdf` | RENAMED ONLY | Supplementary Figure S1 (study-design schematic). |
| `results/model_framework_revision/figures/Figure1_representation_framework.png` | `figures/source/Figure1_representation_framework.png` | COPIED UNCHANGED | Vector source of main Figure 1. |
| `results/model_framework_revision/supplementary/SupplementaryFigureS1_study_design.png` | `figures/Supplementary_Figure_S1.png` | RENAMED ONLY | Supplementary Figure S1 (study-design schematic). |
| `results/model_framework_revision/figures/Figure1_representation_framework.svg` | `figures/source/Figure1_representation_framework.svg` | COPIED UNCHANGED | Vector source of main Figure 1. |
| `results/model_framework_revision/supplementary/SupplementaryFigureS1_study_design.svg` | `figures/Supplementary_Figure_S1.svg` | RENAMED ONLY | Supplementary Figure S1 (study-design schematic). |
| `results/GPB_submission_FINAL_v2/supplement/04_GPB_Supplementary_Information.pdf` | `supplementary/Supplementary_Information.pdf` | RENAMED ONLY | Supplementary Information (journal-neutral filename). |
| `results/GPB_submission_FINAL_v2/supplement/05_GPB_Supplementary_Data.zip` | `supplementary/Supplementary_Data.zip` | RENAMED ONLY | Supplementary Data archive (frozen tables). |

## Excluded, superseded and consolidated material (grouped)

| SOURCE | TARGET | ACTION | REASON |
| --- | --- | --- | --- |
| results/GPB_submission_FINAL_v2/02_GPB_Cover_Letter.docx, 03_GPB_Metadata_Checklist.docx, 00_SUBMISSION_README.md | — | EXCLUDED | Journal-submission-only material (cover letter, metadata checklist, upload notes, author-information placeholders, declarations drafts). |
| results/GPB_submission_FINAL_v2/01_GPB_Main_Manuscript.docx, archive/, logs/, scripts/build_gpb_*.py, gpb_text.py, audit_gpb_package.py, export_figures.py | — | EXCLUDED | Manuscript files, manuscript builders, submission audits; generate_gpb_figures.py migrated as scripts/generate_publication_figures.py. |
| Root-level manuscript DOCX/PDF files (CITEseq_Protein_Access_Model_Gains_FINAL*.docx/pdf, Reliability_CITEseq_Manuscript_v2.pdf, *_BACKUP_*, *_v1.pdf) | — | EXCLUDED | Manuscript drafts, backups and historical PDFs. |
| results/final_manuscript_revision/ | — | EXCLUDED | Historical Genome Biology-addressed draft material (manuscript draft and cover letter). |
| results/final_submission_package/, final_application_manuscript_v2/, final_information_vs_model_manuscript/, phase13_manuscript/, phase14_submission_prep/ | — | SUPERSEDED | Earlier manuscript/figure packages superseded by the final figures and frozen tables. |
| results/manuscript_submission_revision/manuscript/, supplement/, audit/, figures/, logs/, *.md reports; scripts/build_final_docx.py, regenerate_fig1_fig6.py | — | SUPERSEDED | Manuscript sources, builders and superseded figure scripts; analysis scripts and frozen tables migrated. |
| results/manuscript_submission_revision/strict_train_only/models/, embeddings/, logs/ | — | EXCLUDED | Model checkpoints (*.pt), latent embeddings and run logs; regenerated by scripts/run_strict_train_only_internal.py. |
| results/manuscript_submission_revision/tables/development_heldout_predictions.npz | — | EXCLUDED | Prediction cache regenerated by scripts/cache_feature_transductive_predictions.py. |
| results/framework_revision/, results/model_framework_revision/ (reports, DOCX/PDF builds) | — | SUPERSEDED | Figure-revision working folders; final Figure 1 and Supplementary Figure S1 files migrated. |
| results/revision_round2_methodological_fixes/ (audit/, baselines/, figures/, figures_preview/, per-run JSON, done markers, progress reports; scripts 00_*, revision2_p1_orchestrate.py, revision2_rna_only_annotation.py) | — | SUPERSEDED | Working audits and orchestration; analysis scripts and frozen tables migrated. |
| results/revision_round3_information_vs_model/transfer_control/done/, figures_preview/ | — | EXCLUDED | Completion markers (contain absolute machine paths) and preview figures. |
| results/phase10_cross_dataset/, phase10b_source_size/, phase11_uncertainty/, phase12a/ logs+figures, phase11b/ logs+figures+models | — | EXCLUDED | Run artifacts (models, embeddings, probabilities, per-run logs, preview figures); regenerated by migrated scripts. |
| results/phase12b_external_validation/raw_audit/ | — | EXCLUDED | Raw/processed third-party Lawlor data and derived matrices (not redistributed); R helper migrated as scripts/prepare_lawlor_baseline_matrices.R. |
| results/phase12b_external_validation/formal_validation/ per-run tables, done/, embeddings/, models/, posterior/, probabilities/, processed/, figures/ | — | EXCLUDED | Per-run intermediates, checkpoints and processed data; aggregated tables migrated. |
| results/annotation/, results/exploratory/, results/tables/ (except protein_sparsity_* tables), results/figures/, results/logs/, results/models/, results/embeddings/, results/perturbations/ | — | EXCLUDED | Pipeline intermediates, exploratory outputs and checkpoints. |
| results/presubmission_revision/bootstrap/core_contrast_bootstrap.png, ../presubmission_revision/ (outside repository) | — | SUPERSEDED | Preview figure; out-of-repository copies are byte-identical duplicates of tracked files. |
| scripts/03_*, 04_*, 05_*, 07_run_mofa.py, 09_*, 10_validate_query_api.py, 10c_rescore_corrected.py, 11_figures.py, 11b_figures.py, 12a_figures.py | — | SUPERSEDED | Historical scVI baseline / MOFA+ / missing-modality wrappers and superseded figure scripts (D). |
| scripts/10_run_cross_dataset.py, 10b_run_source_size.py, 11_run_uncertainty.py, 11_gate.py, 11_subset_latent.py, 11_augment_delta.py, 11b_mc_check.py, 11b_methods_log.py, 12a_determinism_probe.py, 12a_integrity_audit.py, 12a_power_check.py, 12a_write_logs.py, 12b_orchestrate.py, 10_phase10_cuda.sh, run_on_gpu.sh, server_cuda_smoke.py, server_env_check.py | — | EXCLUDED | Development-only gates, probes, server launchers and log writers (C); several contain server paths. |
| src/ modules not imported by migrated scripts (e.g. models/run_mofa.py, models/baseline_pca.py, experiments/phase5.py, phase7_mofa.py, phase9_api_validation.py, phase10_api_validation.py, evaluation/biological_metrics.py, latent_geometry.py, prediction_metrics.py, protein_marker_sanity.py, plots.py, annotation/audit.py) | — | EXCLUDED | Not required by any migrated script (import closure). |
| tests/test_phase11b.py, tests/test_phase12a.py, notebooks/exploratory_analysis.ipynb | — | EXCLUDED | Development tests and exploratory notebook (C). |
| environment-linux-cuda.yml | — | EXCLUDED | Server-specific CUDA environment variant; pinned versions preserved in environment.yml / requirements.txt. |
| README.md (development repository) | see reason | CONSOLIDATED | Rewritten as the publication README; analysis tables and provenance carried over unchanged. |
| results/manuscript_submission_revision/REPRODUCIBILITY_README.md, MANUSCRIPT_OUTPUT_MANIFEST.md | see reason | CONSOLIDATED | Merged into docs/REPRODUCIBILITY.md and docs/RESULT_MANIFEST.md. |
| AI assistant transcripts, prompts, Cursor/ChatGPT conversation exports, AI-disclosure drafts | — | EXCLUDED | Not research outputs; never part of the publication repository. |
| data/raw/, data/processed/ (*.h5ad), model checkpoints (*.pt, *.ckpt), caches | — | EXCLUDED | Raw/processed data and checkpoints are not redistributed. |

## Code normalisation record

Every PATH NORMALIZED code file received only the following kinds of change: repository-root resolution (`parents[n]`), output/input directory strings, module renames in import statements (development stage labels removed from file names), output file names carrying stage labels, usage strings naming renamed scripts, one configuration key name (`phase10` → `query_adaptation`, value unchanged), and an optional `CUDA_LIB_DIR` environment variable replacing a hard-coded server library path. The full edit record is in `docs/migration_edit_record.json`.

### Module renames

| Development module | Publication module |
| --- | --- |
| `src/experiments/phase10_cross_dataset.py` | `src/experiments/cross_dataset_transfer.py` |
| `src/experiments/phase10b_source_size.py` | `src/experiments/source_size_subsampling.py` |
| `src/experiments/phase11_uncertainty.py` | `src/experiments/uncertainty.py` |
| `src/experiments/phase11b_robustness.py` | `src/experiments/transfer_replicates.py` |
| `src/experiments/phase11b_analysis.py` | `src/experiments/transfer_replicate_analysis.py` |
| `src/experiments/phase11b_stats.py` | `src/experiments/transfer_replicate_stats.py` |
| `src/experiments/phase12a_replication.py` | `src/experiments/stochasticity_replication.py` |
| `src/experiments/phase12a_stats.py` | `src/experiments/stochasticity_stats.py` |
| `src/experiments/phase12b_external.py` | `src/experiments/lawlor_external.py` |
| `src/experiments/phase8_sparsity.py` | `src/experiments/protein_sparsity_stress.py` |
| `src/experiments/phase6_run.py` | `src/experiments/reference_annotation.py` |
| `src/evaluation/phase6_plots.py` | `src/evaluation/annotation_plots.py` |
| `src/evaluation/phase6.py` | `src/evaluation/annotation_evaluation.py` |
| `src/evaluation/phase8_plots.py` | `src/evaluation/protein_sparsity_plots.py` |
| `src/evaluation/phase10_plots.py` | `src/evaluation/transfer_plots.py` |
| `src/evaluation/phase10b_plots.py` | `src/evaluation/source_size_plots.py` |

### Directory mapping used by migrated code

| Development path | Publication path |
| --- | --- |
| `results/phase12b_external_validation/formal_validation` | `results/donor_validation/lawlor_runs` |
| `results/phase12b_external_validation/raw_audit` | `data/external/lawlor` |
| `results/phase12b_external_validation/tables` | `results/donor_validation/lawlor_dataset_audit` |
| `results/phase10b_source_size` | `results/transfer/source_size_runs` |
| `results/phase10_cross_dataset` | `results/transfer/cross_dataset_runs` |
| `results/phase11b_uncertainty_robustness` | `results/transfer/replicate_runs` |
| `results/phase11_uncertainty` | `results/uncertainty/single_seed_runs` |
| `results/phase12a_variance_replication` | `results/stochasticity/replicate_runs` |
| `results/revision_round3_information_vs_model/transfer_control` | `results/transfer/matched_concat_control` |
| `results/revision_round2_methodological_fixes/controls/test_time_corruption` | `results/degradation/test_time_corruption` |
| `results/revision_round2_methodological_fixes/controls` | `results/supplementary/controls` |
| `results/revision_round2_methodological_fixes/statistics/lawlor_` | `results/donor_validation/statistics/lawlor_` |
| `results/presubmission_revision/bootstrap` | `results/supplementary/information_access_sensitivity/feature_transductive_bootstrap` |
| `results/tables/protein_sparsity_primary_delta_summary.csv` | `results/degradation/training_time/protein_sparsity_primary_delta_summary.csv` |
