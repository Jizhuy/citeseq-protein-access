# Supplementary package plan (publication-facing)

Internal phase numbers must not appear in published titles.

| ID | Title | Purpose | Source artifact | Cited from |
|---|---|---|---|---|
| S1 | Matched scVI geometry controls | Show scVI_matched vs default geometry | `results/figures/scvi_matched_seed0_umap_batch.png`, latent geometry tables | Methods / Result 1 |
| S2 | MOFA+ diagnostics | Variance explained / view structure | `mofa_seed0_variance_explained.csv`, UMAP | Result 1 |
| S3 | Additional representation metrics | ASW, kNN purity, l1 metrics | `mofa_biological_representation_metrics.csv` | Result 1 |
| S4 | Per-protein recovery under corruption | Marker-level recovery detail | `phase8_figure7_per_protein_recovery.png` | Result 2 |
| S5 | Transfer coordinate / design diagnostics | Post-adaptation alignment erratum summary | PHASE10 erratum + corrected tables | Result 3 / Methods |
| S6 | Multi-seed transfer distributions | Full 20-seed metric spreads | `phase11b_20seed_primary.csv` / figures | Result 3–5 |
| S7 | Ensemble behavior and ECE | Discrimination vs calibration trade-offs | `phase11b_ensemble*.csv` / fig06 | Result 4 (supp) |
| S8 | Posterior variance Monte Carlo check | Confirm variance semantics | `phase11b_posterior_mc_check.csv` | Methods |
| S9 | Replicated variance heatmaps | Full component maps | phase12a figures 02–03 | Result 6 |
| S10 | Variance power / identifiability | Power caveats for interaction | phase12a logs | Result 6 |
| S11 | External gene/protein harmonization | Mapping + 12-protein panel | phase12b gene/protein tables | Methods / Result 8 |
| S12 | External per-fold metrics | Fold-resolved Lawlor results | `phase12b_primary_runs.csv` / fig09 | Result 8 |
| S13 | Full cell-type retention curves | All coverages / classes | phase11b + phase12b retention tables | Result 7–8 |
| S14 | Runtime and resource notes | Reproducibility environment | orchestrator / methods logs | Methods |

## Supplementary tables
ST1 Dataset inventory · ST2 Full representation metrics · ST3 Corruption deltas · ST4 Corrected transfer seed table · ST5 20-seed uncertainty · ST6 Variance components · ST7 Lawlor 100-run primary · ST8 Software versions · ST9 Novelty/prior-art matrix (optional)

## Excluded from supplement
Internal phase diaries, obsolete mixed-coordinate tables, smoke-test logs, novelty “project planning” notes not needed for reproducibility.
