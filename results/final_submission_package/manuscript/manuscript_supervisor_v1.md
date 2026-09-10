# Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction

**Status:** Frozen scientific draft for supervisor review (Genome Biology — Research).  
**Backup journal:** Bioinformatics — Original Paper.  
**Do not submit without human approval of placeholders below.**

---

## Abstract

**Background:** CITE-seq jointly profiles RNA and surface proteins, and deep generative models such as scVI and totalVI provide probabilistic representations for integration and prediction. How reliable these models remain under protein degradation, dataset or donor shift, training stochasticity, and uncertainty-guided filtering is less well characterized than average accuracy.

**Results:** On independent Lawlor Baseline CITE-seq with author labels (5 donor-held-out folds × 10 seeds = 100 runs), totalVI improved donor-held-out macro-F1 from 0.800 to 0.888 (paired Δ = +0.088; 50/50 pairs), with predictive error AUROC 0.801 versus 0.838 and latent posterior error AUROC 0.605 versus 0.650. In PBMC development analyses, multimodal classification gains were positive on average under entry-wise protein corruption, but reciprocal dataset transfer was direction-dependent (near-null in one direction; Δ ≈ +0.052 with 19/20 seeds favoring totalVI in the other). Predictive uncertainty detected errors with AUROC ≈0.84–0.86. Latent posterior variance tracked downstream error more consistently for totalVI than for matched scVI. Replicated variance decomposition attributed most previously lumped interaction/residual variability to run-level stochastic residual variance. Confidence abstention reduced prediction risk but produced severe cell-type retention imbalance.

**Conclusions:** Deep generative multi-omics models can improve biological prediction under shift, yet uncertainty is multi-dimensional and abstention trades aggregate reliability for biological coverage. These patterns replicate under independent donor-held-out validation without requiring a new architecture.

**Keywords:** single-cell multi-omics; CITE-seq; scVI; totalVI; uncertainty; selective prediction; distribution shift; reliability

---

## Background

Single-cell multi-omics assays such as CITE-seq couple transcriptomes with antibody-derived surface-protein counts [1]. Integrating these modalities is difficult because RNA and protein measurements differ in noise, sparsity, background, and batch structure.

Deep generative models address part of this challenge. scVI learns probabilistic RNA representations [2], and totalVI extends the framework to joint RNA–protein likelihoods [3]. MOFA+ provides a complementary statistical baseline [4], and scArches supports query adaptation across datasets [5]. Existing evaluations often emphasize integration scores, clustering agreement, or average classification accuracy [6].

Less clear is reliability when proteins are degraded, when datasets or donors shift, when training stochasticity is non-negligible, and when uncertainty is used to filter predictions. Abstention on low-confidence cells can change which biological populations remain visible. Uncertainty itself is multi-dimensional: predictive confidence, proper scoring rules, calibration error, and latent posterior width answer related but distinct questions.

Here we evaluate representation quality, entry-wise protein corruption, cross-dataset transfer, predictive and latent uncertainty, stochastic fragility, selective prediction with cell-type retention, and independent donor-held-out external validation. We do not propose a new fusion architecture; the contribution is a disciplined reliability characterization.

---

## Results

### Multimodal representations improve selected aspects of biological prediction

We compared PCA (RNA), MOFA+, architecture-matched scVI (scVI_matched), and totalVI on an independently annotated PBMC CITE-seq query (Figure 2). Metrics included held-out logistic-regression macro-F1, neighborhood purity, cell-type silhouette, and batch silhouette.

No model dominated every metric. On fine-grained (l2) high-confidence labels, totalVI achieved macro-F1 0.779 versus 0.722 for scVI_matched (Δ = +0.057; bootstrap 95% CI [0.020, 0.091]) and 0.550 for MOFA+. PCA remained competitive in some settings and had the highest batch silhouette, which should not be read as biological failure. Neighborhood and silhouette rankings were not identical to classification rankings.

### Multimodal advantages persist under severe entry-wise protein corruption

We corrupted protein count entries at fractions from 0 to 0.85 and compared totalVI to a fixed scVI_matched RNA reference (Figure 2). Mean Δ macro-F1 remained positive at every level (e.g., +0.057 at 0%, +0.049 at 10%, +0.034 at 85%), although some intermediate intervals included zero and the curve was not monotonic. This is entry-wise degradation of observed protein values, not a cell-wise missing-modality experiment.

### Dataset shift reveals that multimodal benefit is context-dependent

We evaluated reciprocal PBMC10k↔PBMC5k transfer with scArches adaptation, encoding source and target with the same post-adaptation model (Figure 3). In a 20-seed grid, Direction A (PBMC10k→PBMC5k) was near null for macro-F1 (0.656 ± 0.018 vs 0.652 ± 0.015; paired Δ ≈ −0.003). Direction B (PBMC5k→PBMC10k) favored totalVI (0.602 ± 0.020 vs 0.654 ± 0.025; Δ ≈ +0.052; 19/20 seeds). Source-size matching partially narrowed the gap but did not fully explain it. Multimodal transfer benefit is therefore not universal across shift directions.

### Predictive uncertainty identifies unreliable cells

Using \(U_{\mathrm{pred}}=1-\max_k p(y=k\mid z)\) from multinomial logistic regression on post-adaptation latents, predictive error AUROC was ≈0.84–0.86 across directions and models (Figure 4; Table 2). NLL, Brier, and AURC often favored totalVI in paired comparisons, whereas predictive AUROC gains were smaller and more seed-dependent. ECE differences were not uniformly decisive. Error discrimination, probability quality, and calibration should be reported separately.

### Latent posterior uncertainty is model-dependent

Defining \(U_{\mathrm{latent}}=\frac{1}{d}\sum_j\mathrm{Var}[z_j\mid x]\), totalVI latent error AUROC exceeded scVI_matched in Direction A (0.669 ± 0.021 vs 0.544 ± 0.054) and Direction B (0.610 ± 0.038 vs 0.469 ± 0.075), with paired advantages in 19/20 seeds per direction (Figure 4). We interpret this as posterior variance that more consistently tracked downstream classification error for totalVI—not as evidence that protein supervision “calibrates” posteriors. Independent Lawlor validation reproduced the direction at attenuated magnitude (0.650 vs 0.605). Because scVI latent AUROC exceeded chance externally, scVI posterior variance should not be described as generally uninformative.

### Reliability varies substantially across stochastic runs

Full-model uncertainty correlated with prediction instability under source-subset and seed perturbations, but did not identify the causal source of fragility (Figure 5). A replicated Direction A design (5 source subsets × 5 model seeds × 2 run replicates × 2 models) separated components. Residual/run-level fractions were large (mean ≈0.66 for scVI_matched and ≈0.80 for totalVI across key metrics), whereas source-subset main effects were small. For Δ(totalVI−scVI), residual likewise dominated (mean fraction ≈0.69). Most of the previously lumped interaction/residual variability is therefore run-level stochastic residual variability, with limited power to exclude modest interactions.

### Abstention lowers risk but creates biological coverage imbalance

Reducing coverage lowered retained-set risk and AURC, especially for totalVI (Figure 5). At 50% coverage, cell-type retention was highly imbalanced (PBMC retention CV ≈0.88–1.02). On Lawlor Baseline, retention CV remained ≈0.57–0.59, with CD4/CD8 naive and CD4 memory preferentially rejected relative to monocytes and B cells (Figure 6). This is cell-type retention imbalance, not a social-fairness claim. Confidence-based abstention is used as a diagnostic baseline, not as a new selective algorithm.

### Independent donor-held-out validation reproduces the central findings

We validated key claims on Lawlor Baseline CITE-seq using author-provided labels only—independent of the Hao/Seurat-v4 reference used solely for development-cohort annotation (Figure 6; Tables 1–2). Design: 16,175 Baseline singlets; 5,207 labeled cells; 12,776 shared genes; 12 matched proteins; five deterministic two-donor target folds; 10 seeds (100 runs; 0 failures).

totalVI improved macro-F1 from 0.800 ± 0.015 to 0.888 ± 0.012 (paired Δ = +0.088; 50/50). Predictive AUROC was 0.801 vs 0.838; latent AUROC 0.605 vs 0.650; AURC 0.087 vs 0.036. ECE was essentially tied. For most metrics, between-fold variability exceeded within-fold seed variability. These results support external reproducibility without introducing a new method.

---

## Discussion

Multimodal deep generative modeling can improve fine-grained biological prediction, but advantages depend on metric and shift context. A near-null transfer direction and a reproducible gain can coexist; adequate seeding is required to avoid overinterpreting small asymmetries.

Average representation scores are insufficient for reliability assessment. Classification, neighborhood structure, and batch silhouette can disagree, and protein-entry degradation can preserve average multimodal advantages while still stressing protein recovery.

Predictive and latent uncertainty measure distinct properties. Predictive confidence is a strong error detector; probability quality and selective risk often favor totalVI; calibration does not automatically follow. Latent posterior width can contain useful error information, but that utility is model-dependent and attenuated—yet directionally preserved—under independent donor-held-out validation.

Training stochasticity is a first-class reliability factor. Once run replicates are available, residual stochastic variability—not source composition alone—explains much of the observed metric and Δ variance.

Uncertainty-guided abstention creates a reliability–coverage trade-off: aggregate risk falls while cell-type retention becomes uneven. External replication indicates this imbalance is not cohort-specific.

Architectural novelty is unnecessary to expose these reliability patterns. Precision-weighted Product-of-Experts fusion is prior art and is not claimed as a contribution [7].

### Limitations

Both development and external cohorts are PBMC. Lawlor labels are coarse and only partially available; the external protein panel overlaps 12 of 14 development proteins. Stimulation shift and full 39-protein sensitivity were not tested. Source-composition perturbation was not repeated externally. Variance-factor designs have limited levels. Arbitrary cell-wise missing protein modality was not honestly testable with scvi-tools totalVI 1.3.3 without panel-level missingness encoding; no fabricated missing-modality experiment is reported. Predictive uncertainty is downstream of latent representations rather than an end-to-end classification posterior. Abstention can disproportionately exclude difficult populations.

---

## Conclusions

Deep generative single-cell multi-omics models can improve biological prediction under dataset and donor shift, while probabilistic uncertainty provides useful but multi-dimensional reliability information. Uncertainty-guided abstention reduces prediction risk at the cost of biological coverage imbalance. Independent donor-held-out validation reproduces these central findings. Optional extensions (stimulation shift, broader protein panels) are not required to support the present conclusions.

---

## Methods

### Study overview
Frozen computational experiments evaluated representation quality, entry-wise protein corruption, cross-dataset transfer, source-size sensitivity, uncertainty, replicated stochastic variability, selective prediction, and Lawlor donor-held-out validation. Models were not retrained for manuscript preparation.

### Datasets and labels
**PBMC10k / PBMC5k CITE-seq** (10x Genomics public CITE-seq) were used for development analyses. Development biological labels were transferred independently of the evaluated latents via a Seurat v4 multimodal reference pathway [8]; those reference labels were **not** applied to the external cohort.

**Lawlor Baseline CITE-seq** provided independent external validation [9]: HCA project `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005`. Analysis used Baseline singlets only (stimulated conditions excluded from the primary analysis). Author protein-gated labels were used exclusively (`obs["author_cell_type"]`). Gene identifiers were mapped with Ensembl GRCh37.87; duplicate gene symbols were aggregated by summing raw UMI counts. Shared gene universe after intersection: **12,776** symbols. Primary protein panel: **12** matched proteins (CD15 and TIGIT absent; high-confidence aliases PD-1↔CD279_PD1 and CD127↔CD127_IL7Ra). Cell counts: **16,175** Baseline singlets; **5,207** author-labeled. Batch covariate for `setup_anndata`: sequencing lane (`run_identifier`); **donor was not used as batch**. Donor-held-out folds were frozen before performance inspection: donors ordered Donor01–Donor10; folds = consecutive pairs (01–02)…(09–10); each donor is target in exactly one fold.

### Models
- **PCA (RNA):** 20 components on the same 2,000 HVG scaled RNA matrix used for MOFA+ RNA view.
- **MOFA+:** mofapy2 **0.7.5**; 20 factors; RNA (2,000 HVGs, library-size normalized, log1p, scaled) + protein (CLR then z-score); Gaussian likelihoods; `scale_views=True`, `weight_views=True` [4].
- **scVI_matched:** n_latent=20, n_hidden=256, n_layers=2, dropout=0.2, gene_likelihood=NB, dispersion=gene, latent_distribution=normal.
- **totalVI:** n_latent=20, n_hidden=256, encoder layers=2, decoder layers=1, dropout=0.2, gene_likelihood=NB; proteins via `obsm["protein_expression"]` (14 in development; 12 externally).

Training (transfer/uncertainty/external grids): batch_size=256, max_epochs=200, train_size=0.9, early stopping on `elbo_validation` (patience=45). scArches query adaptation: `prepare_query_anndata` + `load_query_data`, max_epochs=200, plan `weight_decay=0.0` [5]. **Post-adaptation coordinate rule:** source and target embeddings were always encoded with the adapted query model (never mixed pre-/post-adaptation coordinates).

Software (recorded manifests): Python 3.11.7; scvi-tools 1.3.3; scanpy 1.10.3; anndata 0.11.4; PyTorch 2.14.0 (CUDA 12.6 on GPU runs).

### Protein corruption
Entry-wise corruption of protein count entries at predefined fractions (0–0.85). Primary endpoint: l2 high-confidence logistic-regression macro-F1 difference (totalVI − scVI_matched). This is not a missing-modality design.

### Uncertainty, classifier, and selective prediction
Predictive uncertainty: \(U_{\mathrm{pred}}(i)=1-\max_k p_{ik}\).  
Latent uncertainty: \(U_{\mathrm{latent}}(i)=\frac{1}{d}\sum_{j=1}^{d}\mathrm{Var}[z_{ij}\mid x_i]\) from `get_latent_representation(..., return_dist=True)`.  
Classifier: multinomial logistic regression (`sklearn`, `random_state=seed`) fitted on source post-adaptation latents; evaluated on target labeled cells.  
Metrics: error AUROC/AUPRC, NLL, multiclass Brier, ECE (equal-frequency/width), AURC, risk–coverage at coverages {1.0,0.9,…,0.5}, per-class retention. Confidence abstention only.

### Variance decomposition
Balanced source-subset × model-seed × run-replicate method-of-moments decomposition for metrics and Δ(totalVI−scVI), reporting raw and nonnegative components (Direction A replicated design).

### External validation design
5 folds × 2 models × 10 seeds = 100 runs. Eligible Lawlor classes required ≥20 labeled source and ≥10 labeled target cells in every fold: B, CD14_Mono, CD4T_Mem, CD4T_Naive, CD8T_Mem, CD8T_Naive, NK.

### Statistics
Effect sizes with seed or fold SDs and sign consistency. Cells are not treated as independent experimental units when seed or donor is the relevant factor. No new significance tests were added for manuscript decoration.

### Missing-protein modality limitation
scvi-tools totalVI 1.3.3 does not support arbitrary cell-wise missing protein modality without panel/batch encoding of missingness. Zero-filling would imply observed zeros; therefore no cell-wise missing-modality experiment is reported.

---

## Declarations (human input required where marked)

**Ethics approval and consent to participate.** Not applicable for secondary analysis of publicly deposited human single-cell datasets as used under their original access terms. [Confirm with supervisor/institution if additional statement required.]

**Consent for publication.** Not applicable.

**Availability of data and materials.**  
Public PBMC10k/PBMC5k CITE-seq datasets as specified in project manifests. Lawlor data: HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376`. Processed analysis tables and figure sources are retained in the project `results/` tree. **Public archival DOI (Zenodo):** [NOT YET CREATED — do not claim].

**Code availability.** Analysis code is in the project repository `multiomics_robustness` (`scripts/`, `src/experiments/`). **Public release / archival DOI:** [PENDING HUMAN APPROVAL — repository public status not asserted here].

**Competing interests.** The authors declare that they have no competing interests. [Confirm.]

**Funding.** [PLACEHOLDER — insert grant IDs]

**Authors’ contributions.** [PLACEHOLDER — CRediT roles]

**Acknowledgements.** [PLACEHOLDER]

**Authors and affiliations.** [PLACEHOLDER — names, order, affiliations, corresponding author email]

---

## References

1. Stoeckius M, Hafemeister C, Stephenson W, Houck-Loomis B, Chattopadhyay PK, Swerdlow H, Satija R, Smibert P. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods. 2017;14:865–868. doi:10.1038/nmeth.4380.
2. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15:1053–1058. doi:10.1038/s41592-018-0229-2.
3. Gayoso A, Steier Z, Lopez R, Regier J, Nazor KL, Streets A, Yosef N. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods. 2021;18:272–282. doi:10.1038/s41592-020-01050-x.
4. Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, Stegle O. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol. 2020;21:111. doi:10.1186/s13059-020-02015-1.
5. Lotfollahi M, Naghipourfar M, Luecken MD, Khajavi M, Büttner M, Wagenstetter M, Avsec Ž, Gayoso A, Yosef N, Interlandi M, et al. Mapping single-cell data to reference atlases by transfer learning. Nat Biotechnol. 2022;40:121–130. doi:10.1038/s41587-021-01001-7.
6. Luecken MD, Büttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colomé-Tatché M, Theis FJ. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods. 2022;19:41–50. doi:10.1038/s41592-021-01336-8.
7. Wu M, Goodman N. Multimodal generative models for scalable weakly-supervised learning. Adv Neural Inf Process Syst. 2018;31. https://papers.nips.cc/paper/2018/hash/1102a326d5f7c9e04fc3c29e95bd2ffa-Abstract.html
8. Hao Y, Hao S, Andersen-Nissen E, Mauck WM III, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184:3573–3587.e29. doi:10.1016/j.cell.2021.04.048.
9. Lawlor N, Nehar-Belaid D, Grassmann JDS, Stoeckius M, Smibert P, Stitzel ML, Pascual V, Banchereau J, Williams A, Ucar D. Single Cell Analysis of Blood Mononuclear Cells Stimulated Through Either LPS or Anti-CD3 and Anti-CD28. Front Immunol. 2021;12:636720. doi:10.3389/fimmu.2021.636720.
10. Guo C, Pleiss G, Sun Y, Weinberger KQ. On calibration of modern neural networks. Proc Mach Learn Res. 2017;70:1321–1330.
11. Geifman Y, El-Yaniv R. Selective classification for deep neural networks. Adv Neural Inf Process Syst. 2017;30.
12. Ashuach T, Gabitto MI, Koodli RV, Saldi GA, Jordan MI, Yosef N. MultiVI: deep generative model for the integration of multimodal data. Nat Methods. 2023;20:1222–1231. doi:10.1038/s41592-023-01909-9.

Unresolved bibliographic details (NeurIPS page-range format; expand truncated author lists if required): see `references/unresolved_references.md`.

---

## Figure legends

**Figure 1. Study overview.** What was evaluated? PBMC10k/PBMC5k development analyses (representation, entry-wise protein corruption, reciprocal transfer, uncertainty/stochasticity, selective prediction) and independent Lawlor Baseline donor-held-out validation with author labels. Models: PCA, MOFA+, scVI_matched, totalVI (not all models used in every experiment).

**Figure 2. Representation and protein degradation.** How do representation choice and entry-wise protein corruption affect biological performance? (A) Example latents. (B) Fine-grained held-out classification. (C–D) Macro-F1 and totalVI−scVI_matched differences under entry-wise corruption. (E) Protein recovery. Entry-wise corruption ≠ missing modality.

**Figure 3. Dataset shift.** How does multimodal benefit behave under reciprocal dataset shift? (A) Design with post-adaptation coordinates. (B) Multi-seed macro-F1 structure. (C) Paired deltas. (D) Source-size sensitivity. Direction A near null; Direction B favors totalVI.

**Figure 4. Predictive and latent uncertainty.** What information is contained in predictive confidence versus latent posterior variance? (A) Predictive error AUROC. (B) Latent uncertainty by correctness. (C) Predictive vs latent association. (D) Paired metric deltas. Discrimination ≠ probability quality ≠ calibration.

**Figure 5. Stochastic fragility and selective prediction.** How does stochastic fragility relate to abstention and biological retention? (A) Uncertainty vs instability. (B) Replicated variance decomposition (source, seed, interaction, residual). (C) Risk–coverage. (D) Cell-type retention imbalance.

**Figure 6. Independent Lawlor validation.** Do central findings reproduce under donor-held-out external validation? (A) Ten-donor / five-fold design. (B) Macro-F1. (C) Predictive AUROC. (D) Latent AUROC. (E) Risk–coverage. (F) Class retention. Author labels only; 12-protein panel; 100 runs.

---

## Tables

Table 1 (datasets) and Table 2 (primary reliability estimates) are provided as `tables/Table1_datasets.*` and `tables/Table2_primary_results.*`.

---

## Title decision

**Final recommendation (unchanged):**  
*Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction*

Compared with “Reliability and uncertainty of deep generative models for single-cell multi-omics integration,” the selected title is more specific about distribution shift and uncertainty-guided prediction without implying architectural novelty. Prefer accuracy and clarity over brevity alone.
