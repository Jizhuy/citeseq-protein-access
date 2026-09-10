# Abstract

**Background:** CITE-seq jointly profiles transcriptomes and antibody-derived surface proteins. Probabilistic deep generative models such as scVI and totalVI are widely used for representation learning, yet reliability under protein-entry degradation, dataset or donor shift, training stochasticity, and uncertainty-guided filtering remains incompletely characterized. We evaluated architecture-matched scVI (scVI_matched), totalVI, PCA, and MOFA+ in PBMC CITE-seq cohorts under repeated-seed protocols.

**Results:** In the combined high-confidence PBMC analysis, totalVI improved l2-cell-type F1 over scVI_matched (0.779 vs 0.722; Δ+0.057; target-cell paired bootstrap CI [0.020, 0.091]); MOFA+ reached 0.550. Under Lawlor donor-fold evaluation, scVI achieved F1 0.800 and totalVI 0.888 (Δ+0.088) across 5 folds × 10 seeds, giving 50 paired comparisons within 100 runs. Predictive uncertainty separated correct from incorrect predictions with internal AUROC approximately 0.84-0.86 and Lawlor AUROC 0.801 vs 0.838 for scVI and totalVI. Latent posterior uncertainty was more informative for totalVI than scVI internally (approximately 0.61-0.67 vs 0.47-0.54) and externally (0.650 vs 0.605). Protein-entry corruption, cross-dataset transfer, source-size downsampling, and variance partitioning showed context-dependent totalVI gains, with residual run-level variability largest. Selective prediction risk decreased as coverage was reduced, and totalVI achieved a lower overall AURC, but retained coverage differed across cell types.

**Conclusions:** Across PBMC CITE-seq evaluations, totalVI usually provided stronger representation and uncertainty behavior than scVI_matched, while reliability depended on transfer direction, measured-protein perturbation, and selective-retention imbalance.

## Keywords

CITE-seq; single-cell multi-omics; scVI; totalVI; distribution shift; uncertainty; selective prediction; PBMC
