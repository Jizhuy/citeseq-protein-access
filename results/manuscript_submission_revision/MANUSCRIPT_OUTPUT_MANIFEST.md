# Manuscript Output Manifest

| Manuscript item | Source results | Generation script | Final output | Primary numbers | Status |
| --- | --- | --- | --- | --- | --- |
| Development F1 table | `results/revision_round2_methodological_fixes/controls/tables/simple_multimodal_development.csv` | `revision2_simple_multimodal_baseline.py` | Main Table 1 | 0.666 / 0.512 / 0.745 / 0.549 / 0.722 / 0.779 | verified |
| Core contrasts + bootstrap | `results/presubmission_revision/bootstrap/core_contrast_bootstrap_summary.csv` | `presubmission_core_contrast_bootstrap.py` | Main Table 2 / Fig 2 | G_access +0.079; G_assoc +0.034; R≈0.70 | verified |
| Class-stratified bootstrap | `results/manuscript_submission_revision/tables/bootstrap_global_vs_class_stratified.csv` | `run_class_stratified_bootstrap.py` | Suppl. SR1 | G_assoc CI still crosses 0 | verified |
| Transfer + matched concat | `results/revision_round3_information_vs_model/transfer_control/tables/matched_transfer_control.csv` | `revision3_matched_concat_transfer.py` | Main Table 3 / Fig 3 | 0.735/0.734; 0.636/0.767; Δ +0.016/−0.113 | verified |
| Training-time degradation | `results/presubmission_revision/degradation_audit/training_time_degradation_summary.csv` | historical stress + audit | Main Table 4 / Fig 4 | 0.779…0.755 curve | verified |
| Post-adapt Dir B degradation | test-time corruption summary CSVs | `revision2_testtime_corruption.py` | Fig 4 | 0.654→0.621 | verified |
| Class-level Dir B degradation | per-class corruption CSV | `summarize_class_degradation.py` | Suppl. SR2 | CD8 Naive Δ≈−0.264 | verified |
| Selective macro-F1 | `phase11b_selective_prediction.csv` | `summarize_selective_macroF1.py` | Suppl. SR3 | cov50 macro-F1↑ | verified |
| Stochasticity MoM | `phase12a` tables | `12a_analyze.py` | Main Table 5 | residual fractions | verified |
| Lawlor donor Δ | `lawlor_donor_delta.csv` | `03_lawlor_donor_level.py` + `lawlor_donor_summary_stats.py` | Main Table 6 / Fig 6 | Δ+0.090; 10/10; means 0.794/0.884 | verified |
| Lawlor protein/concat | `simple_multimodal_lawlor.csv` | revision2 simple baseline | Table 6 | ≈0.927 / 0.945 | verified |
| Marker-channel dropout | `marker_channel_dropout_pca_development.csv` | `marker_channel_dropout_pca.py` | Suppl. SR5 | T-lineage panel Δ concat −0.025; protein −0.185 | verified |
| Environment | `runtime_environment.json` | `scripts/print_environment.py` | Suppl. S10 | py3.11.7; scvi 1.3.3 | verified |
| Canonical manuscript | this revision folder | manual Methods expansion | `manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md` | protein-access framing | updated |

## Reproduction status legend

- **verified**: frozen numbers match code/tables  
- **running / pending**: script available; output not yet frozen  
- **blocked**: see SUPPLEMENTARY_RESULTS SR6
