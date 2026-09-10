# Reliability and Uncertainty in Deep Generative Single-Cell Multi-Omics Integration

**Recommended title (from ranked candidates).** Alternate candidates appear at the end of this draft.

---

## Abstract

Single-cell multi-omics assays such as CITE-seq jointly measure RNA and surface proteins, and deep generative models including scVI and totalVI provide probabilistic latent representations for integration and prediction. In practice, however, reliability under protein degradation, dataset and donor shift, stochastic training variability, and uncertainty-guided filtering remains incompletely characterized. Here we systematically evaluate representation quality, robustness, transfer, uncertainty, selective prediction, and external reproducibility for architecture-matched scVI and totalVI, with PCA and MOFA+ as representation baselines.

On an independently annotated peripheral-blood mononuclear cell (PBMC) CITE-seq cohort, totalVI improved fine-grained held-out classification relative to scVI_matched (l2 high-confidence logistic-regression macro-F1 0.779 versus 0.722; bootstrap Δ = +0.057). Mean multimodal advantage remained positive across entry-wise protein corruption fractions up to 0.85, without implying monotonicity or a missing-modality regime. Cross-dataset transfer was direction-dependent: in a 20-seed evaluation, Direction A (PBMC10k→PBMC5k) was effectively near-null for macro-F1, whereas Direction B (PBMC5k→PBMC10k) favored totalVI (mean Δ ≈ +0.052; 19/20 seeds).

Predictive uncertainty strongly detected classification errors under shift (predictive error AUROC typically ≈0.84–0.86). Latent posterior variance was model-dependent: totalVI latent error AUROC exceeded scVI in both transfer directions (≈0.67 vs 0.54 in Direction A; ≈0.61 vs 0.47 in Direction B). Replicated variance decomposition showed that run-level residual stochasticity dominates many metric and Δ variances once source-subset and seed factors are separated. Confidence-based selective prediction reduced retained-set risk but produced severe cell-type retention imbalance.

These reliability patterns were tested on an independent Lawlor Baseline CITE-seq cohort using author labels only (16,175 Baseline singlets; 5,207 labeled; 12,776 shared genes; 12 proteins; 5 donor-held-out folds × 10 seeds = 100 runs). totalVI improved donor-held-out macro-F1 from 0.800 to 0.888 (paired Δ = +0.088; 50/50 pairs), with predictive AUROC 0.801 vs 0.838 and latent AUROC 0.605 vs 0.650. We conclude that deep generative multi-omics models can improve biological prediction under shift, while probabilistic uncertainty provides useful but multi-dimensional reliability information and abstention trades aggregate risk for biological coverage.

---

## Introduction

Single-cell multi-omics technologies enable joint characterization of molecular layers within the same cell. CITE-seq, in particular, couples transcriptomes with antibody-derived surface-protein counts and has become a workhorse for immune-cell phenotyping. Integrating these views is statistically challenging because modalities differ in noise, sparsity, background, and batch structure.

Deep generative models provide a flexible response to this challenge. scVI learns probabilistic RNA representations with explicit library-size and dispersion structure, while totalVI extends the framework to joint RNA–protein likelihoods. Related factor and multimodal models, including MOFA+, offer complementary statistical baselines. Existing benchmarks often emphasize integration scores, clustering agreement, or average predictive accuracy.

Less systematically characterized is whether these models remain reliable when proteins are degraded, when labels and batches shift across datasets or donors, when training stochasticity is non-negligible, and when uncertainty is used to filter predictions. These stresses matter for real analyses: practitioners may abstain on low-confidence cells, potentially altering which biological populations remain visible.

Uncertainty itself is multi-dimensional. Predictive confidence from a downstream classifier, proper scoring rules, calibration error, and latent posterior width answer related but distinct questions. Collapsing them into a single claim that one model has “better uncertainty” obscures when discrimination, probability quality, or calibration diverge.

Here we evaluate reliability of deep generative single-cell multi-omics models across representation quality, protein-entry corruption, cross-dataset transfer, predictive and latent uncertainty, stochastic fragility, selective prediction with cell-type retention, and independent donor-held-out external validation. We do not introduce a new fusion architecture; the contribution is a disciplined reliability characterization with claim-level traceability to frozen computational experiments.

---

## Results

### Result 1. Benchmarking biological representation across statistical and deep generative models

We compared PCA on RNA, MOFA+, architecture-matched scVI (scVI_matched), and totalVI using independently derived cell-type labels on a PBMC CITE-seq query. Metrics were intentionally plural: held-out logistic-regression macro-F1, neighborhood purity, cell-type silhouette (ASW), and batch silhouette.

No model dominated every metric. On fine-grained (l2) high-confidence labels, totalVI achieved the highest logistic-regression macro-F1 (0.779), exceeding scVI_matched (0.722; Δ = +0.057, bootstrap 95% CI [0.020, 0.091]) and MOFA+ (0.550). PCA attained competitive classification in some settings and the highest batch silhouette, which should not be misread as biological failure: batch silhouette and biological discrimination answer different questions. Neighborhood and ASW rankings were not identical to classification rankings, reinforcing that metric choice shapes conclusions about “integration quality.”

### Result 2. Multimodal gains remain directionally robust under severe entry-wise protein corruption

To stress protein information without claiming true missingness, we applied entry-wise corruption to protein counts at fractions from 0 to 0.85 while holding the RNA pathway fixed for the scVI_matched reference comparison. Mean totalVI−scVI_matched l2 macro-F1 differences remained positive at every evaluated corruption level (for example, Δ ≈ +0.057 at 0%, +0.049 at 10%, and +0.034 at 85%), although interval estimates at intermediate levels sometimes included zero and the curve was not monotonic. This experiment is an entry-wise degradation study, not a cell-wise missing-modality study.

### Result 3. Multimodal transfer gains depend on distribution-shift direction and source context

We evaluated bidirectional transfer between PBMC10k and PBMC5k with scArches query adaptation, encoding both source and target with the post-adaptation model. Under a 20-seed primary grid, Direction A (PBMC10k→PBMC5k) macro-F1 was effectively near-null (scVI_matched 0.656 ± 0.018; totalVI 0.652 ± 0.015; paired Δ ≈ −0.003, inconclusive). Direction B (PBMC5k→PBMC10k) showed a reproducible totalVI advantage (0.602 ± 0.020 vs 0.654 ± 0.025; paired Δ ≈ +0.052; 19/20 seeds positive).

Earlier corrected five-seed summaries are consistent with this picture after removal of obsolete mixed-coordinate artifacts, but Direction A should not be summarized as a robust multimodal disadvantage. Source-size sensitivity analyses indicate that matching source cell numbers partially narrows the Direction A/B gap yet does not fully explain it. Transfer conclusions are therefore context-dependent rather than universally multimodal.

### Result 4. Predictive uncertainty identifies unreliable cells under dataset shift

Using predictive uncertainty \(U_{\mathrm{pred}} = 1 - \max_k p(y=k\mid z)\) from multinomial logistic regression on post-adaptation latents, error-detection AUROC was high in both directions and models (approximately 0.84–0.86 across the 20-seed grid). Proper scoring rules often favored totalVI (lower NLL and Brier), and selective-risk summaries (AURC) likewise improved in robust paired comparisons, whereas predictive error-ranking AUROC gains were smaller and more seed-dependent. Expected calibration error (ECE) differences were not uniformly decisive. These results separate error discrimination, probability quality, and calibration.

### Result 5. Latent posterior uncertainty is model-dependent

Latent uncertainty was defined as mean posterior variance across latent dimensions, \(U_{\mathrm{latent}} = \frac{1}{d}\sum_j \mathrm{Var}[z_j\mid x]\), extracted with distributional returns from adapted models. In the 20-seed transfer grid, totalVI latent error AUROC exceeded scVI_matched in Direction A (0.669 ± 0.021 vs 0.544 ± 0.054) and Direction B (0.610 ± 0.038 vs 0.469 ± 0.075), with paired advantages in 19/20 seeds per direction. We interpret this as an association between multimodal training and posterior widths that more consistently track downstream error, not as a demonstrated causal mechanism of “protein supervision causing calibrated posteriors.”

Independent Lawlor donor-held-out validation reproduced the direction of the effect at attenuated magnitude (scVI 0.605; totalVI 0.650; paired Δ ≈ +0.045; 38/50 pairs). Notably, scVI latent AUROC was above chance externally, so scVI posterior variance should not be described as universally uninformative.

### Result 6. Uncertainty predicts stochastic prediction fragility, while run-level residual dominates variance

Full-model uncertainty correlated with prediction instability under source-subset and seed perturbations, indicating that uncertain cells are often fragile cells. Uncertainty does not, however, identify whether fragility arises from source composition, seed, or residual training noise.

A replicated Direction A design (5 source subsets × 5 model seeds × 2 run replicates × 2 models) separated these components. Averaged across key metrics, residual/run-level fractions were large (approximately 0.66 for scVI_matched and 0.80 for totalVI in the metric-level decomposition), whereas source-subset main effects were small. For Δ(totalVI−scVI), residual likewise dominated (mean residual fraction ≈ 0.69), with modest interaction and limited source contribution. We therefore describe the formerly lumped “interaction/residual” mass as largely run-level stochastic variability, not an unidentified biological residual, while noting limited factor-level power for interaction tests.

### Result 7. Uncertainty-guided selective prediction reduces risk but distorts biological coverage

Confidence abstention improved retained-set accuracy and reduced risk as coverage decreased, with lower AURC for totalVI in both transfer and external settings. At 50% coverage, however, cell-type retention was highly imbalanced (PBMC retention CV commonly ≈0.88–1.02). On Lawlor Baseline, retention CV remained high (≈0.57–0.59), with CD4/CD8 naive and CD4 memory preferentially rejected relative to monocytes and B cells. We term this cell-type retention imbalance (biological coverage imbalance), not a social-fairness claim. Selective prediction is used here as a baseline diagnostic of the risk–coverage–biology trade-off, not as a new abstention algorithm.

### Result 8. Independent donor-held-out validation reproduces central reliability findings

We validated key claims on Lawlor Baseline CITE-seq using author labels only. Design: 10 donors, five deterministic two-donor target folds, 16,175 Baseline singlets for unsupervised training, 5,207 labeled cells for classifier evaluation, 12,776 shared genes, 12 overlapping proteins, and 10 seeds (100 runs; 0 failures).

Primary biological performance favored totalVI (macro-F1 0.888 ± 0.012 vs 0.800 ± 0.015; paired Δ = +0.088; 50/50). Predictive uncertainty remained strongly error-informative (AUROC 0.838 vs 0.801). Latent uncertainty directionally favored totalVI as above. AURC improved (0.036 vs 0.087). ECE was essentially tied. For most metrics, between-donor-fold variability exceeded within-fold seed variability. This section demonstrates external reproducibility rather than a new method.

---

## Discussion

Multimodal deep generative modeling can improve fine-grained biological prediction, but advantages are metric- and context-dependent. Cross-dataset transfer in particular is not a universal win for totalVI: one direction is near-null under adequate seeding, while the reverse direction is reproducible.

Uncertainty quality is multi-dimensional. Predictive confidence is a strong error detector; proper scoring rules and selective risk often favor totalVI; calibration (ECE) does not automatically follow. Latent posterior width carries usable reliability information that is stronger for totalVI in our settings, replicates directionally externally, and should not be over-causalized.

Uncertainty foreshadows stochastic fragility without identifying its source. Once run replicates are available, residual stochastic variability—not source composition alone—explains much of the observed metric variance. Evaluation protocols that report a single seed therefore risk overinterpreting differences.

Selective prediction exposes a central practical trade-off: aggregate risk falls while biological coverage becomes uneven. External replication of this imbalance argues it is not a quirk of one cohort.

Finally, architectural novelty is unnecessary to uncover these reliability failures and successes. Precision-weighted Product-of-Experts fusion is prior art and is not positioned as a contribution here. The contribution is systematic reliability characterization under sparsity, shift, stochasticity, uncertainty, abstention, and external donor-held-out validation.

### Limitations

Primary development cohorts are PBMC. External validation is also PBMC, with coarse and partially available author labels and a 12/14 protein overlap. Stimulation shift and full 39-ADT sensitivity were not tested. Source-composition perturbation was not repeated externally. Variance-factor designs have limited levels. Comparisons center on scVI_matched/totalVI with PCA/MOFA+ baselines. Cell-wise missing-protein modality could not be honestly formulated in totalVI 1.3.3 without encoding panel structure. Predictive uncertainty is downstream of latents rather than an end-to-end classification posterior. Abstention improves retained-set reliability at the cost of uneven coverage.

---

## Methods

### Study overview
Frozen computational experiments evaluated representation, entry-wise protein corruption, cross-dataset transfer, source-size sensitivity, uncertainty, replicated stochastic variability, selective prediction, and Lawlor donor-held-out validation. No new models were trained during manuscript consolidation.

### Datasets
**PBMC10k / PBMC5k CITE-seq** development cohorts with independent label transfer for biological evaluation. **Lawlor Baseline CITE-seq** external cohort (Front. Immunol. 2021): Baseline singlets only; author labels only; gene mapping via Ensembl GRCh37.87; 12-protein primary panel (aliases PD-1↔CD279_PD1, CD127↔CD127_IL7Ra); donor-held-out folds frozen before performance inspection; batch key = sequencing lane (`run_identifier`), not donor.

### Models
PCA (RNA); MOFA+; scVI_matched (n_latent=20, n_hidden=256, n_layers=2, dropout=0.2, NB likelihood); totalVI (matched latent/hidden widths; encoder 2 / decoder 1 layers; 12-protein panel externally). Query adaptation used scArches; source and target embeddings always come from the post-adaptation model.

### Protein corruption
Entry-wise masking/corruption of protein counts at predefined fractions; primary endpoint l2 high-confidence macro-F1 delta versus scVI_matched.

### Uncertainty definitions
Predictive: \(U_{\mathrm{pred}}(i)=1-\max_k p_{ik}\).  
Latent: \(U_{\mathrm{latent}}(i)=\frac{1}{d}\sum_{j=1}^{d}\mathrm{Var}[z_{ij}\mid x_i]\).  
Downstream classifier: multinomial logistic regression on latents. Metrics: error AUROC/AUPRC, NLL, multiclass Brier, ECE (equal-frequency/width), AURC, risk–coverage, per-class retention.

### Variance decomposition
Balanced source-subset × model-seed × replicate ANOVA/method-of-moments decomposition for metrics and Δ(totalVI−scVI), reporting raw and nonnegative components.

### Statistical reporting
Effect sizes with seed/fold SDs and sign consistency. Cells are not treated as independent experimental replicates when seed or donor is the relevant factor. No new significance tests were added solely for manuscript decoration.

### Software and compute
Python 3.11.7; scvi-tools 1.3.3; scanpy 1.10.3; torch 2.14.0 (CUDA 12.6 on GPU server). Exact versions in `results/logs/environment.json` and phase methods logs.

### Missing-modality limitation
scvi-tools totalVI 1.3.3 does not support arbitrary cell-wise missing protein modality without panel/batch encoding. Zero-filling would imply observed zeros; no fabricated missing-modality experiment is reported.

---

## Data availability
PBMC public CITE-seq cohorts as used in project manifests. Lawlor data via Human Cell Atlas / ENA accessions recorded in the external-validation audit logs (`PRJEB40376` / HCA dataset UUID). Processed intermediates and tables are in the project `results/` tree.

## Code availability
Analysis code in the project repository (`multiomics_robustness`). Experiment scripts for transfer, uncertainty, variance replication, and external validation are under `scripts/` and `src/experiments/`.

## Acknowledgements
[Placeholder]

## Author contributions
[Placeholder]

## Competing interests
The authors declare no competing interests.

## References
See `references/reference_audit.csv` for citation keys, DOIs, and verification status. Primary method citations: Lopez et al., Nat Methods 2018 (scVI); Gayoso et al., Nat Methods 2021 (totalVI); Stoeckius et al., Nat Methods 2017 (CITE-seq); Argelaguet et al., Genome Biol 2020 (MOFA+); Lotfollahi et al., Nat Biotechnol 2022 (scArches); Hao et al., Cell 2021 (Seurat v4 reference used only for development-cohort annotation); Lawlor et al., Front Immunol 2021 (external cohort); Wu & Goodman, NeurIPS 2018 (PoE multimodal VAE prior art).

---

## Candidate titles (ranked)

1. **Reliability and Uncertainty in Deep Generative Single-Cell Multi-Omics Integration** ← recommended
2. Robustness and Uncertainty of Deep Generative Models for Single-Cell Multi-Omics Integration
3. Evaluating Reliability Under Distribution Shift in Deep Generative Single-Cell Multi-Omics Models
4. Uncertainty-Aware Evaluation of Deep Generative Single-Cell Multi-Omics Integration
5. When Multimodal Single-Cell Models Are Reliable: Transfer, Uncertainty, and Selective Prediction
6. Donor Shift, Stochasticity, and Uncertainty in Deep Generative CITE-seq Models
7. Reliability Trade-offs in Probabilistic Single-Cell Multi-Omics Integration
8. External Validation of Uncertainty and Robustness in Deep Generative Multi-Omics Models
9. Predictive and Latent Uncertainty in Single-Cell Multi-Omics Deep Generative Models
10. Benchmarking Robustness and Selective Prediction in Single-Cell Multi-Omics VAEs
11. From Representation to Reliability: Deep Generative Models for CITE-seq Under Shift
12. Characterizing Failure Modes of Deep Generative Single-Cell Multi-Omics Integration
