**Figure 1. Study overview.** Development PBMC10k/PBMC5k CITE-seq cohorts were used for representation benchmarking (PCA, MOFA+, scVI_matched, totalVI), entry-wise protein corruption, reciprocal transfer, uncertainty/stochasticity analyses, and selective prediction. Independent Lawlor Baseline CITE-seq provided donor-held-out external validation with author labels. Not all models are used in every experiment.

**Figure 2. Representation quality and protein-corruption robustness.** (A) Example latent visualizations for PCA, MOFA+, scVI_matched, and totalVI. (B) Fine-grained held-out classification. (C–D) Macro-F1 and totalVI−scVI_matched differences under entry-wise protein corruption. (E) Protein recovery summary. Entry-wise corruption is not missing-modality missingness.

**Figure 3. Cross-dataset transfer and source-context sensitivity.** (A) Reciprocal transfer design with post-adaptation coordinates. (B) Seed-resolved macro-F1 structure from the multi-seed grid. (C) Paired totalVI−scVI deltas. (D) Source-size sensitivity. Direction A is near null; Direction B favors totalVI.

**Figure 4. Predictive and latent uncertainty under dataset shift.** (A) Predictive error-detection AUROC distributions. (B) Latent uncertainty for correct vs incorrect predictions. (C) Predictive vs latent uncertainty association. (D) Paired seed deltas across uncertainty-related metrics. Discrimination, probability quality, and calibration are distinct.

**Figure 5. Stochastic fragility and selective prediction.** (A) Uncertainty versus source-subset instability. (B) Replicated variance decomposition (source, seed, interaction, residual). (C) Risk–coverage. (D) Cell-type retention imbalance under abstention.

**Figure 6. Independent Lawlor donor-held-out validation.** (A) Ten-donor / five-fold design. (B) Macro-F1. (C) Predictive error AUROC. (D) Latent error AUROC. (E) Risk–coverage. (F) Class retention. Author labels only; 12-protein panel; 100 runs.

---
