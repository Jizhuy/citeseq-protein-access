# Main Figure Legends


### Figure 1. Study overview of reliability analyses in PBMC CITE-seq.

Schematic of the evaluated reliability programme. Development analyses used PBMC10k and PBMC5k CITE-seq (combined analysis object: 10,849 cells, 15,792 genes, 14 shared proteins; 9,494 high-confidence cells at annotation confidence ≥0.85). Experiments covered representation quality, entry-wise protein corruption, reciprocal cross-dataset transfer with post-adaptation coordinates, predictive and latent uncertainty, replicated stochastic variability, and selective prediction with cell-type retention. Independent external validation used Lawlor Baseline CITE-seq with author-provided labels only (16,175 Baseline singlets; 5,207 labeled cells; 12,776 shared genes; 12 matched proteins; 5 donor-held-out folds × 10 seeds = 100 runs). Models: PCA, MOFA+, scVI_matched, and totalVI (not all models in every experiment). Abbreviations: PBMC, peripheral blood mononuclear cells.

### Figure 2. Representation quality and entry-wise protein corruption.

(A) Example latent visualizations for PCA, MOFA+, scVI_matched, and totalVI in the development PBMC analysis. (B) Fine-grained (l2) held-out logistic-regression macro-F1 on high-confidence labels (n=9,494 cells). totalVI = 0.779, scVI_matched = 0.722 (Δ = +0.057; target-cell paired bootstrap 95% CI [0.020, 0.091], n=500), MOFA+ = 0.550. (C–D) Macro-F1 and totalVI−scVI_matched differences under entry-wise protein corruption. Randomly selected originally non-zero measured protein entries were set to zero at fractions 0, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85 (five mask seeds for nonzero levels); preexisting zeros were unchanged; totalVI was retrained; scVI_matched remained fixed. (E) Protein recovery summary. Entry-wise corruption is not cell-wise missing modality. Error bars/bands summarize corruption-mask or bootstrap variability as shown in panels.

### Figure 3. Reciprocal dataset transfer under post-adaptation coordinates.

(A) Design for Directions A (PBMC10k→PBMC5k) and B (PBMC5k→PBMC10k) with scArches query adaptation; source and target embeddings used the same post-adaptation model. (B–C) Macro-F1 across 20 model seeds and paired deltas. Direction A was near null (scVI_matched 0.656 ± 0.018; totalVI 0.652 ± 0.015). Direction B favored totalVI (scVI_matched 0.601 ± 0.020; totalVI 0.654 ± 0.025; Δ ≈ +0.052; 19/20 seeds). (D) Source-size sensitivity after downsampling PBMC10k to 3,994 cells (five source subsets). Values are mean ± SD over seeds; SD is over the 20 model seeds unless otherwise noted.

### Figure 4. Predictive confidence versus latent posterior variance.

(A) Predictive error AUROC using U_pred = 1 − max_k p(y=k|z) after multinomial logistic regression on post-adaptation latents (internal AUROC ≈0.84–0.86). (B) Latent uncertainty by prediction correctness, where U_latent is the mean of posterior variances across d=20 latent dimensions from scvi-tools get_latent_representation(..., return_dist=True). (C) Association between predictive and latent uncertainty. (D) Paired seed-level metric deltas. Direction A latent error AUROC: scVI_matched 0.544 ± 0.054 versus totalVI 0.669 ± 0.021; Direction B: 0.469 ± 0.075 versus 0.610 ± 0.038 (19/20 seeds favoring totalVI per direction). Error discrimination, probability quality (NLL/Brier), calibration (ECE), and selective risk (AURC) are reported as distinct properties. Error bars denote seed SD (n=20 seeds).

### Figure 5. Stochastic variability, selective prediction, and retention.

(A) Association between uncertainty and prediction instability under source-subset and seed perturbations. (B) Replicated variance decomposition for Direction A (5 source subsets × 5 model seeds × 2 run replicates × 2 models) under Y = μ + A_s + B_m + (AB)_sm + ε_smr; residual/run-level components were largest and are interpreted as run-level residual under the replication scheme, not irreducible biology. (C) Risk–coverage curves; prediction risk decreased as coverage was reduced from 100% to 50%, and totalVI achieved a lower overall AURC. (D) Cell-type retention imbalance at reduced coverage (PBMC retention CV ≈0.88–1.02 at 50% coverage). Confidence abstention is a diagnostic baseline, not a new selective algorithm.

### Figure 6. Independent Lawlor donor-held-out validation.

(A) Ten-donor / five-fold design (two held-out donors per fold; each donor target once). (B) Macro-F1: scVI_matched 0.800 ± 0.015 versus totalVI 0.888 ± 0.012 (paired Δ = +0.088; 50/50 fold×seed pairs descriptively favor totalVI; pairs are not 50 independent biological replicates). (C) Predictive error AUROC 0.801 versus 0.838. (D) Latent error AUROC 0.605 versus 0.650. (E) Risk–coverage; AURC 0.087 versus 0.036. (F) Class retention (CV ≈0.57–0.59), with CD4/CD8 naive and CD4 memory preferentially rejected relative to monocytes and B cells. Author-provided labels only; 12-protein panel; 12,776 shared genes; 100 runs. Where shown, error bars summarize fold×seed variation as indicated in panels.

