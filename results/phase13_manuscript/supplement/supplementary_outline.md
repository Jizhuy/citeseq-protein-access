# Supplementary outline

## Supplementary Methods
1. Independent cell-type annotation protocol (Seurat-v4 reference for PBMC development cohorts only).
2. Matched scVI hyperparameter alignment to totalVI.
3. scArches adaptation and post-adaptation coordinate rule (with PHASE 10 erratum summary).
4. Entry-wise protein corruption mask generation.
5. Uncertainty metric definitions (equations).
6. Selective prediction / risk–coverage / retention statistics.
7. Replicated ANOVA / method-of-moments variance decomposition and power caveats.
8. Lawlor gene symbol mapping (Ensembl GRCh37.87) and 12-protein alias rules.
9. Compute environment and software versions.
10. Why cell-wise missing-protein experiments were not reported (API limitation).

**Populate from:** phase11b_methods.md; phase10 erratum; phase8 provenance; phase12b data_prep / setup_anndata; environment.json.

## Supplementary Results
1. Full metric tables for representation benchmark (l1/l2; ASW; kNN; batch silhouette).
2. Protein-recovery and per-marker corruption analyses.
3. Direction-wise transfer class tables.
4. Source-size subsample-level results.
5. 20-seed distributions and ensemble analyses (incl. ECE trade-offs).
6. Posterior MC variance confirmation.
7. Full variance-component tables and heatmaps.
8. Complete Lawlor per-fold / per-class tables.
9. Replication status table expanding PHASE 12B summary.

## Supplementary Figures (proposed index)
| ID | Content | Artifact source |
|---|---|---|
| S1 | Default vs matched scVI UMAPs / geometry | latent_geometry_comparison_seed0.csv + figures |
| S2 | MOFA variance explained | mofa_seed0_variance_explained.csv |
| S3 | Full corruption curves + CIs | protein_sparsity_primary_delta*.csv |
| S4 | Latent drift / correction diagnostic | phase10 erratum + corrected vs obsolete deltas |
| S5 | Source-size composition diagnostics | phase10b tables |
| S6 | 20-seed metric beeswarms | phase11b_20seed_primary.csv |
| S7 | Ensemble convergence / ECE | phase11b_ensemble*.csv |
| S8 | MC posterior variance check | phase11b_posterior_mc_check.csv |
| S9 | 5×5 variance heatmaps | phase12a tables/figures |
| S10 | Power analysis for interaction | phase12a logs |
| S11 | All-class retention curves (PBMC) | phase11b_celltype_retention.csv |
| S12 | Lawlor gene/protein harmonization | phase12b gene/protein tables |
| S13 | Lawlor per-fold metric panels | phase12b_primary_runs.csv |
| S14 | Runtime/RAM notes | phase12b orchestrator logs |

## Supplementary Tables
| ID | Content |
|---|---|
| ST1 | Dataset inventory |
| ST2 | Full representation metrics |
| ST3 | Corruption delta table |
| ST4 | Corrected transfer seed-level table |
| ST5 | 20-seed uncertainty metrics |
| ST6 | Variance component estimates (raw + nonnegative) |
| ST7 | Lawlor primary 100-run table |
| ST8 | Evidence ledger (internal; optional) |
| ST9 | Software versions |
| ST10 | Novelty audit matrix |
