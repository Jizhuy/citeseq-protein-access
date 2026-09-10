# Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction

**Final recommended title.**  
Alternate candidates are listed at the end.

---

## Abstract

**Background:** CITE-seq and related assays jointly measure RNA and surface proteins, and deep generative models such as scVI and totalVI provide probabilistic representations for integration and cell-state prediction. Reliability under protein degradation, dataset or donor shift, training stochasticity, and uncertainty-based filtering is less systematically characterized than average predictive accuracy.

**Results:** Across PBMC development cohorts, totalVI improved fine-grained held-out classification relative to architecture-matched scVI (l2 high-confidence macro-F1 0.779 vs 0.722; bootstrap Δ = +0.057), and mean multimodal advantage remained positive under entry-wise protein corruption up to 85% of entries. Cross-dataset transfer was direction-dependent: in a 20-seed evaluation, PBMC10k→PBMC5k macro-F1 differences were near null, whereas PBMC5k→PBMC10k favored totalVI (mean Δ ≈ +0.052; 19/20 seeds). Predictive uncertainty detected errors with AUROC ≈0.84–0.86. Latent posterior variance tracked error more consistently for totalVI (latent AUROC ≈0.61–0.67) than for scVI (≈0.47–0.54). Replicated variance decomposition attributed most previously lumped interaction/residual variability to run-level stochastic residual variance. Confidence abstention reduced risk but produced severe cell-type retention imbalance. On independent Lawlor Baseline CITE-seq (author labels; 5 donor folds × 10 seeds = 100 runs), totalVI raised donor-held-out macro-F1 from 0.800 to 0.888 (Δ = +0.088; 50/50 pairs), with predictive AUROC 0.801 vs 0.838 and latent AUROC 0.605 vs 0.650.

**Conclusions:** Deep generative multi-omics models can improve biological prediction under shift, but uncertainty is multi-dimensional and abstention trades aggregate reliability for biological coverage. These patterns replicate under independent donor-held-out validation without requiring a new architecture.

**Keywords:** single-cell multi-omics; CITE-seq; scVI; totalVI; uncertainty; selective prediction; distribution shift; reliability

---

## Background

Single-cell multi-omics assays such as CITE-seq couple transcriptomes with antibody-derived surface-protein counts, enabling joint characterization of molecular layers within individual cells [1]. Integrating these modalities is challenging because RNA and protein measurements differ in noise, sparsity, background, and batch structure.

Deep generative models address part of this challenge. scVI learns probabilistic RNA representations [2], and totalVI extends the framework to joint RNA–protein likelihoods [3]. Factorization approaches such as MOFA+ provide complementary statistical baselines [4]. Query mapping methods such as scArches support transfer of representations across datasets [5]. Existing evaluations often emphasize integration scores, clustering agreement, or average classification accuracy [6].

Less clear is how reliable these models remain when protein measurements are degraded, when labels and technical covariates shift across datasets or donors, when training stochasticity is non-negligible, and when uncertainty is used to filter predictions. These questions are practical: downstream workflows may abstain on low-confidence cells, potentially changing which biological populations remain visible in a retained analysis set.

Uncertainty is also multi-dimensional. Predictive confidence, proper scoring rules, calibration error, and latent posterior width answer related but distinct questions. Treating them as a single notion of “better uncertainty” can conceal when discrimination, probability quality, and calibration diverge.

Here we systematically evaluate reliability of deep generative single-cell multi-omics models under representation benchmarking, entry-wise protein corruption, cross-dataset transfer, predictive and latent uncertainty, stochastic fragility, selective prediction with cell-type retention reporting, and independent donor-held-out external validation. We do not propose a new fusion architecture; the contribution is a disciplined reliability characterization with claim-level traceability to frozen experiments.

---

## Results

### Multimodal integration improves selected aspects of biological representation

**Question.** Do multimodal deep generative representations improve biological discrimination relative to matched unimodal and classical baselines?

We compared PCA (RNA), MOFA+, architecture-matched scVI (scVI_matched), and totalVI on an independently annotated PBMC CITE-seq query, reporting held-out logistic-regression macro-F1, neighborhood purity, cell-type silhouette, and batch silhouette (Figure 2).

No model dominated every metric. On fine-grained (l2) high-confidence labels, totalVI achieved macro-F1 0.779 versus 0.722 for scVI_matched (Δ = +0.057; bootstrap 95% CI [0.020, 0.091]) and 0.550 for MOFA+. PCA remained competitive for some classification settings and showed the highest batch silhouette, which should not be interpreted as biological failure. Neighborhood and silhouette rankings were not identical to classification rankings, indicating that metric choice materially shapes conclusions about representation quality.

### Multimodal advantages remain robust to severe entry-wise protein corruption

**Question.** Does the multimodal classification advantage collapse when protein entries are heavily corrupted?

We applied entry-wise corruption to protein counts at fractions from 0 to 0.85 while comparing totalVI to a fixed scVI_matched RNA reference (Figure 2). Mean Δ macro-F1 remained positive at every evaluated level (e.g., +0.057 at 0%, +0.049 at 10%, +0.034 at 85%), although some intermediate intervals included zero and the response was not monotonic. This experiment stresses observed protein entries; it is not a cell-wise missing-modality design.

### Cross-dataset transfer benefits depend on shift direction and source context

**Question.** Are multimodal transfer gains stable across reciprocal dataset shifts?

We evaluated PBMC10k↔PBMC5k transfer with scArches adaptation, encoding source and target with the post-adaptation model (Figure 3). In a 20-seed primary grid, Direction A (PBMC10k→PBMC5k) was near null for macro-F1 (scVI_matched 0.656 ± 0.018; totalVI 0.652 ± 0.015; paired Δ ≈ −0.003, inconclusive). Direction B (PBMC5k→PBMC10k) favored totalVI (0.602 ± 0.020 vs 0.654 ± 0.025; Δ ≈ +0.052; 19/20 seeds). Source-size matching partially narrowed the Direction A/B gap but did not fully explain it. Transfer benefits are therefore context-dependent rather than universally multimodal.

### Predictive uncertainty identifies unreliable cells under dataset shift

**Question.** Can predictive confidence detect classification errors after dataset shift?

Using \(U_{\mathrm{pred}}=1-\max_k p(y=k\mid z)\) from multinomial logistic regression on post-adaptation latents, predictive error AUROC was ≈0.84–0.86 across directions and models (Figure 4; Table 2). Proper scoring rules (NLL, Brier) and selective-risk summaries (AURC) often favored totalVI in paired seed comparisons, whereas predictive AUROC gains were smaller and more seed-dependent. ECE differences were not uniformly decisive. Error discrimination, probability quality, and calibration should therefore be reported separately.

### Latent posterior uncertainty is strongly model-dependent

**Question.** Does latent posterior variance carry information about downstream error, and is that property model-dependent?

Defining \(U_{\mathrm{latent}}=\frac{1}{d}\sum_j\mathrm{Var}[z_j\mid x]\), totalVI latent error AUROC exceeded scVI_matched in Direction A (0.669 ± 0.021 vs 0.544 ± 0.054) and Direction B (0.610 ± 0.038 vs 0.469 ± 0.075), with paired advantages in 19/20 seeds per direction (Figure 4). We describe this as posterior variance that more consistently tracked downstream classification error for totalVI, not as proof that protein supervision “calibrates” posteriors. Independent Lawlor validation reproduced the direction at attenuated magnitude (0.650 vs 0.605). Because scVI latent AUROC was above chance externally, scVI posterior variance should not be called uninformative in general.

### Prediction fragility is dominated by run-level stochastic variability

**Question.** When predictions are unstable across perturbations, what variability source dominates?

Full-model uncertainty correlated with instability under source-subset and seed perturbations, but did not identify the causal source of fragility (Figure 5). A replicated Direction A design (5 source subsets × 5 model seeds × 2 run replicates × 2 models) separated components. Residual/run-level fractions were large (mean ≈0.66 for scVI_matched and ≈0.80 for totalVI across key metrics), whereas source-subset main effects were small. For Δ(totalVI−scVI), residual likewise dominated (mean fraction ≈0.69). Most of the previously lumped interaction/residual mass is therefore attributable to run-level stochastic residual variability, with limited power to exclude modest interactions.

### Selective prediction reduces risk but produces cell-type retention imbalance

**Question.** Does confidence abstention improve reliability without distorting biological coverage?

Reducing coverage lowered retained-set risk and AURC, especially for totalVI (Figure 5). At 50% coverage, however, cell-type retention was highly imbalanced (PBMC retention CV ≈0.88–1.02). On Lawlor Baseline, retention CV remained high (≈0.57–0.59), with CD4/CD8 naive and CD4 memory preferentially rejected relative to monocytes and B cells (Figure 6). This is a biological coverage imbalance, not a social-fairness claim, and confidence-based abstention is used here as a diagnostic baseline rather than a new algorithm.

### Independent donor-held-out validation reproduces the central reliability findings

**Question.** Do the central reliability patterns replicate on an independent CITE-seq cohort under donor shift?

We analyzed Lawlor Baseline CITE-seq with author labels only: 16,175 Baseline singlets, 5,207 labeled cells, 12,776 shared genes, 12 matched proteins, five deterministic two-donor target folds, and 10 seeds (100 runs; 0 failures) (Figure 6; Tables 1–2). totalVI improved macro-F1 from 0.800 ± 0.015 to 0.888 ± 0.012 (paired Δ = +0.088; 50/50). Predictive AUROC was 0.801 vs 0.838; latent AUROC 0.605 vs 0.650; AURC 0.087 vs 0.036. ECE was essentially tied. For most metrics, between-fold variability exceeded within-fold seed variability. These results support external reproducibility of the reliability story without introducing a new method.

---

## Discussion

Multimodal deep generative modeling can improve fine-grained biological prediction, but advantages depend on metric and shift context. Reciprocal transfer shows that a near-null direction and a reproducible gain can coexist; adequate seeding is required to avoid overinterpreting small five-seed asymmetries.

Traditional representation metrics are insufficient for reliability assessment. Classification, neighborhood structure, and batch silhouette can disagree, and protein-entry degradation can preserve average multimodal advantages while still stressing protein recovery.

Predictive and latent uncertainty measure distinct properties. Predictive confidence is a strong error detector; probability quality and selective risk often favor totalVI; calibration (ECE) does not automatically follow. Latent posterior width can contain useful error information, but that utility is model-dependent and attenuated—yet directionally preserved—under independent donor-held-out validation.

Training stochasticity is a first-class reliability factor. Once run replicates are available, residual stochastic variability—not source composition alone—explains much of the observed metric and Δ variance. Single-seed evaluations therefore risk fragile conclusions.

Uncertainty-guided abstention creates a reliability–coverage trade-off: aggregate risk falls while cell-type retention becomes uneven. External replication of this imbalance indicates it is not a quirk of one cohort.

Independent Lawlor donor-held-out validation strengthens generalizability for Baseline PBMC CITE-seq. Architectural novelty is unnecessary to expose these reliability successes and failures; precision-weighted Product-of-Experts fusion is prior art and is not claimed as a contribution here [7].

### Limitations

Both development and external cohorts are PBMC. Lawlor labels are coarse and only partially available; the external protein panel overlaps 12 of 14 development proteins. Stimulation shift and full 39-protein sensitivity were not tested. Source-composition perturbation was not repeated externally. Variance-factor designs have limited levels. Arbitrary cell-wise missing protein modality was not honestly testable with scvi-tools totalVI 1.3.3 without panel-level missingness encoding. Predictive uncertainty is downstream of latent representations rather than an end-to-end classification posterior. Abstention can disproportionately exclude difficult populations.

---

## Conclusions

Deep generative single-cell multi-omics models can improve biological prediction under dataset and donor shift, while probabilistic uncertainty provides useful but multi-dimensional reliability information. Uncertainty-guided abstention reduces prediction risk at the cost of biological coverage imbalance. Independent donor-held-out validation reproduces these central findings. Future work may optionally test stimulation shift or broader protein panels, but such experiments are not required to support the present reliability conclusions.

---

## Methods

### Study overview
Frozen computational experiments evaluated representation quality, entry-wise protein corruption, cross-dataset transfer, source-size sensitivity, uncertainty, replicated stochastic variability, selective prediction, and Lawlor donor-held-out validation. No models were retrained during manuscript preparation.

### Datasets
**PBMC10k / PBMC5k CITE-seq** were used for development analyses. Biological labels for these cohorts were derived independently of the evaluated latents (Seurat v4 reference pathway for development annotation only) [8]. **Lawlor Baseline CITE-seq** provided external validation [9]: Baseline singlets only; author labels only; Ensembl GRCh37.87 gene-symbol mapping with duplicate-symbol summation; 12-protein primary panel; donor-held-out folds frozen before performance inspection; batch key = sequencing lane (`run_identifier`), not donor. Accession and HCA identifiers are recorded in project audit logs (ENA `PRJEB40376`).

### Models and adaptation
PCA (RNA); MOFA+; scVI_matched (n_latent=20, n_hidden=256, n_layers=2, dropout=0.2, negative-binomial gene likelihood); totalVI (matched latent/hidden widths; encoder 2 / decoder 1 layers). Query adaptation used scArches [5]; source and target embeddings were always taken from the post-adaptation model. Software: Python 3.11.7; scvi-tools 1.3.3; scanpy 1.10.3; PyTorch 2.14.0 (CUDA 12.6 on GPU runs), as recorded in environment manifests.

### Protein corruption
Entry-wise corruption of protein counts at predefined fractions; primary endpoint l2 high-confidence macro-F1 difference versus scVI_matched.

### Uncertainty and selective prediction
Predictive uncertainty \(U_{\mathrm{pred}}(i)=1-\max_k p_{ik}\). Latent uncertainty \(U_{\mathrm{latent}}(i)=\frac{1}{d}\sum_j\mathrm{Var}[z_{ij}\mid x_i]\). Downstream classifier: multinomial logistic regression. Metrics: error AUROC/AUPRC, NLL, multiclass Brier, ECE, AURC, risk–coverage, per-class retention. Confidence abstention only; no new selective algorithm.

### Variance decomposition
Balanced source-subset × model-seed × run-replicate decomposition for metrics and Δ(totalVI−scVI), reporting raw and nonnegative components.

### External validation design
5 folds × 2 models × 10 seeds = 100 runs; primary eligible Lawlor classes required ≥20 labeled source and ≥10 labeled target cells in every fold (B, CD14_Mono, CD4T_Mem, CD4T_Naive, CD8T_Mem, CD8T_Naive, NK).

### Statistics
Effect sizes with seed or fold standard deviations and sign consistency. Cells are not treated as independent experimental units when seed or donor is the relevant factor. No new significance tests were added solely for manuscript decoration.

---

## Declarations

**Data availability.** Public PBMC CITE-seq cohorts as used in project manifests. Lawlor data via recorded HCA/ENA accessions. Processed tables and figure sources are in the project `results/` tree.

**Code availability.** Analysis code is available in the project repository `multiomics_robustness` (`scripts/`, `src/experiments/`).

**Competing interests.** The authors declare that they have no competing interests.

**Funding.** [Placeholder]

**Authors’ contributions.** [Placeholder]

**Acknowledgements.** [Placeholder]

---

## References

1. Stoeckius M, et al. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods. 2017. doi:10.1038/nmeth.4380.
2. Lopez R, et al. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15:1053–1058. doi:10.1038/s41592-018-0229-2.
3. Gayoso A, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods. 2021;18:272–282. doi:10.1038/s41592-020-01050-x.
4. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol. 2020. doi:10.1186/s13059-020-02015-1.
5. Lotfollahi M, et al. Mapping single-cell data to reference atlases by transfer learning. Nat Biotechnol. 2022. doi:10.1038/s41587-021-01001-7.
6. Luecken MD, et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods. 2022. doi:10.1038/s41592-021-01336-8.
7. Wu M, Goodman N. Multimodal generative models for scalable weakly-supervised learning. NeurIPS. 2018.
8. Hao Y, et al. Integrated analysis of multimodal single-cell data. Cell. 2021. doi:10.1016/j.cell.2021.04.048.
9. Lawlor N, et al. Front Immunol. 2021. doi:10.3389/fimmu.2021.636720.
10. Guo C, et al. On calibration of modern neural networks. ICML. 2017.
11. Geifman Y, El-Yaniv R. Selective classification for deep neural networks. NeurIPS. 2017.
12. Ashuach T, et al. MultiVI: deep generative model for the integration of multimodal data. Nat Methods. 2023. doi:10.1038/s41592-023-01909-9.

*(Full bibliographic audit: `references/reference_audit_final.csv`. Additional prior-art multimodal VAE citations for novelty positioning are retained in the audit file.)*

---

## Figure legends

**Figure 1. Study overview.** Development PBMC10k/PBMC5k CITE-seq cohorts were used for representation benchmarking (PCA, MOFA+, scVI_matched, totalVI), entry-wise protein corruption, reciprocal transfer, uncertainty/stochasticity analyses, and selective prediction. Independent Lawlor Baseline CITE-seq provided donor-held-out external validation with author labels. Not all models are used in every experiment.

**Figure 2. Representation quality and protein-corruption robustness.** (A) Example latent visualizations for PCA, MOFA+, scVI_matched, and totalVI. (B) Fine-grained held-out classification. (C–D) Macro-F1 and totalVI−scVI_matched differences under entry-wise protein corruption. (E) Protein recovery summary. Entry-wise corruption is not missing-modality missingness.

**Figure 3. Cross-dataset transfer and source-context sensitivity.** (A) Reciprocal transfer design with post-adaptation coordinates. (B) Seed-resolved macro-F1 structure from the multi-seed grid. (C) Paired totalVI−scVI deltas. (D) Source-size sensitivity. Direction A is near null; Direction B favors totalVI.

**Figure 4. Predictive and latent uncertainty under dataset shift.** (A) Predictive error-detection AUROC distributions. (B) Latent uncertainty for correct vs incorrect predictions. (C) Predictive vs latent uncertainty association. (D) Paired seed deltas across uncertainty-related metrics. Discrimination, probability quality, and calibration are distinct.

**Figure 5. Stochastic fragility and selective prediction.** (A) Uncertainty versus source-subset instability. (B) Replicated variance decomposition (source, seed, interaction, residual). (C) Risk–coverage. (D) Cell-type retention imbalance under abstention.

**Figure 6. Independent Lawlor donor-held-out validation.** (A) Ten-donor / five-fold design. (B) Macro-F1. (C) Predictive error AUROC. (D) Latent error AUROC. (E) Risk–coverage. (F) Class retention. Author labels only; 12-protein panel; 100 runs.

---

## Tables

See `tables/Table1_datasets.md` and `tables/Table2_primary_results.md`.

---

## Title variants considered

1. **Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction** ← selected  
2. Reliability and Uncertainty in Deep Generative Single-Cell Multi-Omics Integration (prior working title)  
3. Evaluating reliability under dataset and donor shift in deep generative CITE-seq models  
4. Uncertainty-aware reliability assessment of scVI and totalVI for single-cell multi-omics  
5. Distribution shift, stochasticity, and selective prediction in deep generative multi-omics integration  
6. When multimodal single-cell VAEs are reliable: transfer, uncertainty, and biological coverage  
