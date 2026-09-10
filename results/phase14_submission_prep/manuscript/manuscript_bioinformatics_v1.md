# Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction

Bioinformatics Original Paper backup (same claims as Genome Biology version; compressed for ~7-page / ~5,000-word constraint).

---

## Abstract

**Motivation:** Deep generative models such as scVI and totalVI are widely used for CITE-seq integration, yet reliability under protein degradation, dataset/donor shift, training stochasticity, and uncertainty-guided filtering remains incompletely characterized.

**Results:** totalVI improved fine-grained PBMC classification versus matched scVI (macro-F1 0.779 vs 0.722) and retained positive mean advantage under entry-wise protein corruption. Transfer was direction-dependent (near-null one way; Δ≈+0.052 the other; 19/20 seeds). Predictive uncertainty detected errors (AUROC≈0.84–0.86); latent variance tracked error more for totalVI. Residual run variance dominated; abstention cut risk but imbalanced retention. On Lawlor donor-held-out validation (100 runs), macro-F1 rose from 0.800 to 0.888 (Δ=+0.088; 50/50).

**Availability and Implementation:** Code and frozen tables are in the `multiomics_robustness` repository.

**Contact:** [corresponding author email]

**Supplementary Information:** Supplementary Methods, Figures S1–S14, and Tables online.

## Introduction

CITE-seq jointly profiles RNA and surface proteins [1]. Probabilistic models including scVI and totalVI provide latent representations for integration and prediction [2,3], with MOFA+ as a statistical baseline [4] and scArches supporting query adaptation [5]. Benchmarks often emphasize average integration or classification performance [6]. We instead ask how reliable these models are under protein-entry degradation, reciprocal dataset shift, training stochasticity, predictive/latent uncertainty, and selective prediction, and whether findings replicate under independent donor-held-out validation. We do not introduce a new architecture.

## Methods

Development analyses used PBMC10k/PBMC5k CITE-seq with independently derived labels. External validation used Lawlor Baseline CITE-seq with author labels only (16,175 Baseline singlets; 5,207 labeled; 12,776 shared genes; 12 proteins; 5 donor folds × 10 seeds = 100 runs) [9]. Models: PCA, MOFA+, scVI_matched (n_latent=20, n_hidden=256, n_layers=2, dropout=0.2), and totalVI (matched widths). scvi-tools 1.3.3; post-adaptation coordinates for all transfer embeddings. Uncertainty: \(U_{\mathrm{pred}}=1-\max p\); \(U_{\mathrm{latent}}=\mathrm{mean}_j\mathrm{Var}[z_j|x]\). Selective prediction used confidence abstention with per-class retention. Variance decomposition used a replicated source-subset × seed × run design. Full protocols are in Supplementary Methods.

## Results and Discussion

**Representation and corruption.** totalVI exceeded scVI_matched on l2 high-confidence macro-F1 (0.779 vs 0.722; Δ=+0.057, CI [0.020,0.091]); no model dominated every metric. Mean multimodal Δ remained positive under entry-wise protein corruption to 85% of entries (not a missing-modality experiment).

**Transfer.** In 20 seeds, Direction A was near null; Direction B favored totalVI (Δ≈+0.052; 19/20). Source-size matching partially but incompletely explained the asymmetry.

**Uncertainty.** Predictive AUROC≈0.84–0.86. Latent AUROC was higher for totalVI (≈0.61–0.67) than scVI (≈0.47–0.54). NLL/Brier/AURC often favored totalVI; ECE did not uniformly. These are distinct uncertainty facets.

**Stochasticity and abstention.** Replicated decomposition showed residual/run-level variance dominating previously lumped interaction/residual mass. Abstention lowered risk but produced cell-type retention imbalance (PBMC CV≈0.88–1.02 at 50% coverage).

**External validation.** Lawlor donor-held-out: macro-F1 0.800→0.888 (Δ=+0.088; 50/50); predictive AUROC 0.801 vs 0.838; latent AUROC 0.605 vs 0.650; AURC 0.087 vs 0.036. Direction of latent advantage replicated with attenuated magnitude; scVI latent AUROC>0.5 externally.

**Limitations.** PBMC-only cohorts; coarse Lawlor labels; 12/14 protein overlap; no stimulation/39-ADT; no external source-composition repeat; limited variance-factor levels; cell-wise missing-protein modality not honestly supported by totalVI 1.3.3 API; downstream classifier uncertainty.

## Conclusion

Reliability gains of multimodal generative models under shift coexist with multi-dimensional uncertainty behavior and abstention-related biological coverage costs. Independent donor-held-out validation supports these conclusions without architectural novelty.

## Acknowledgements / Funding / Competing interests
[Placeholders as in primary manuscript]

## References
As in Genome Biology version (reformat to OUP author–year if submitting here); see `references/reference_audit_final.csv`.
