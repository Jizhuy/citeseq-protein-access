<!--
Alternative titles (for human review only):
1. Reliability of RNA-only and multimodal CITE-seq representations under cohort shift and measurement degradation
2. Evaluating representation reliability in CITE-seq: protein information, transfer asymmetry, and selective prediction
-->

# Reliability of CITE-seq representations under distribution shift and uncertainty-guided prediction

**Jizhu Yang**  
Research Manuscript  
Unpublished manuscript  
September 2026

## Abstract

**Background:** CITE-seq jointly measures transcriptomes and surface proteins in the same cells, creating an opportunity to test whether RNA-only and multimodal representations remain reliable under dataset shift, protein degradation, stochastic training variation, and uncertainty-guided selective prediction.

**Results:** In a 9,494-cell healthy-control development cohort, paired protein information improved cell-type macro-F1. RNA PCA reached 0.666, concat PCA 0.745, scVI 0.722, and totalVI 0.779; the totalVI–scVI difference was approximately +0.057 (conditional 95% CI 0.020–0.091), and concat PCA recovered about 40% of this gain. Transfer was asymmetric: Direction A was near null, whereas Direction B showed the clearest totalVI advantage (Δ+0.052; 19/20 seeds). Inductive concat PCA remained strong (~0.73). Post-adaptation target-only protein corruption produced gradual declines, stronger in Direction B. Classifier predictive uncertainty and latent posterior variance provided complementary diagnostics. Separately, Lawlor donor-level validation (n=10) showed mean present-class ΔF1 +0.090 (10/10 donors; 0.888 vs 0.800), while protein PCA/concat baselines were higher under protein-gated author labels.

**Conclusions:** Protein information was the strongest source of improved CITE-seq representation quality. Architecture-specific effects, including totalVI advantages, were secondary and context-dependent. Reliable use therefore required stress testing under transfer, degradation, uncertainty, and donor-level validation.


## Background

CITE-seq measures RNA expression and antibody-derived surface proteins in the same single cells, linking broad transcriptional state to targeted phenotypic markers [1]. This paired measurement is especially useful in immune profiling, where cell-type labels often depend on both gene-expression programs and surface-marker patterns. A central question is therefore not only whether protein information improves cell-type prediction, but whether that improvement remains reliable when the data distribution changes or the protein modality is degraded.

Deep generative models provide one way to learn compact representations from single-cell data. scVI models transcript counts with a variational autoencoder and provides a widely used RNA-only latent representation [2]. totalVI extends this framework to paired RNA and protein measurements and can represent CITE-seq data in a joint latent space [3]. MOFA+ and weakly supervised multimodal generative modeling provide related factor-analysis and multimodal latent-variable perspectives [4,8].

The broader single-cell integration literature includes MultiVI, Cobolt, scMaui, atlas-level integration benchmarking, and WNN-style multimodal integration as important related methods and evaluation settings rather than experimental comparators in this manuscript [5,6,7,9,14]. The present analysis focuses on a narrower question: how selected RNA-only and multimodal CITE-seq representations behave under repeated reliability stress tests. This framing separates the value of multimodal information from any claim of universal model ranking.

Reliability is broader than internal accuracy. A representation can perform well in a development cohort but fail under source-target shift, donor differences, protein measurement loss, or stochastic variation in model training. These issues are common in CITE-seq analysis because antibody panels, cell compositions, and annotation conventions can differ across studies. A useful evaluation should therefore include internal comparisons, simple baselines, transfer asymmetry, degradation experiments, repeated runs, and uncertainty diagnostics.

Transfer evaluation is particularly important because source-target direction can change the result. scArches-style adaptation allows a query dataset to be incorporated without labels, but the procedure remains unsupervised transductive adaptation rather than a fully inductive application to unseen cells [10]. Distinguishing these settings matters for interpreting how much information the model used during adaptation and how source-trained classifiers are evaluated on target cells.

Uncertainty was included because representation failures are often not visible from macro-F1 alone. Classifier predictive uncertainty, derived from downstream class probabilities, measures ambiguity at the prediction layer. Latent uncertainty, summarized here as mean posterior variance over latent dimensions, measures instability in the inferred representation. Calibration and uncertainty-aware evaluation are therefore relevant for deciding when automated labels should be trusted or deferred [11,15,16].

Selective prediction offers one operational use of uncertainty: retain low-uncertainty cells and defer high-uncertainty cells for review [12]. This can reduce error among retained cells, but it can also change the retained cell-type composition. In single-cell biology, that compositional change is important because rare, transitional, or marker-ambiguous populations may be preferentially deferred.

**Figure 1.** Study design overview spanning internal development evaluation, reciprocal transfer, protein degradation, uncertainty-guided selective prediction, stochasticity analysis, and Lawlor donor-level validation.

Together, these analyses ask how reliable CITE-seq representations are when multimodal information is useful but imperfect. The governing interpretation is that protein information can strongly improve representation quality, while architecture-specific advantages remain context-dependent and should be judged jointly with simple multimodal baselines, transfer asymmetry, and donor-level validation.

## Results

### Paired protein information improved internal representation quality

The internal development analysis used 9,494 healthy-control cells with matched RNA and protein measurements. Across representation families, paired protein information improved cell-type macro-F1. RNA PCA reached 0.666, protein PCA alone reached 0.512, concatenated RNA-protein PCA reached 0.745, MOFA+ reached 0.549, scVI reached 0.722, and totalVI reached 0.779. The first internal contrast between totalVI and scVI was +0.057 (cell-level conditional paired-bootstrap 95% CI 0.020–0.091), and the totalVI-minus-concat difference was +0.034.

**Figure 2.** Development-cohort comparison of information content and representation strategy, including RNA PCA, protein PCA, concat PCA, MOFA+, scVI, totalVI, totalVI residual differences, and the RNA-only label audit.

The concat PCA baseline was central to interpretation. It substantially exceeded RNA PCA and scVI, recovering approximately 40% of the totalVI-versus-scVI gain without using a specialized deep generative architecture. This result indicates that part of the internal improvement came from the measured protein modality itself, not only from the modeling framework.

Protein PCA alone performed below RNA PCA even though protein information helped strongly when paired with RNA. This pattern is compatible with the structure of CITE-seq panels: the 14 measured proteins were targeted and label-relevant, but too limited to replace genome-wide RNA for the full annotation task. Their strongest value appeared when combined with transcriptomic information.

The RNA-only annotation audit found identical labels across 9,494 cells. This rules out a simple label-mismatch explanation for the RNA-only versus multimodal difference, but it does not independently validate biological correctness. The audit supports the internal comparison as a consistency check while leaving biological label validity to dataset annotation quality and external validation.

The internal development result therefore established a strong protein-information signal and a more cautious model interpretation. totalVI was the best internal representation among the tested methods, but the simple concat baseline showed that paired protein measurements were the primary source of the gain.

### Representation gains were strongly context-dependent under dataset transfer

Bidirectional transfer experiments showed that model-level gains depended strongly on source-target direction. In Direction A, scVI reached 0.656±0.018 and totalVI reached 0.652±0.015, giving a near-null difference. In Direction B, scVI reached 0.601±0.020 and totalVI reached 0.654±0.025, giving Δ+0.052. Direction B therefore provided the clearest evidence of a totalVI advantage under transfer, with 19/20 seeds favoring totalVI.

**Figure 3.** Bidirectional transfer and uncertainty diagnostics, including reciprocal source-target design, unsupervised transductive scArches adaptation, inductive concat PCA, predictive error AUROC, and latent error AUROC.

The inductive concat PCA baseline again helped anchor the interpretation. It reached 0.735 in Direction A and 0.734 in Direction B, above the adapted VAE latent classifiers in both directions. This does not imply that concat PCA is universally preferable; rather, it shows that aligned RNA-protein features can transfer strongly when a source-only feature transform is sufficient for the task.

The Direction A near-null result prevents overgeneralization from the internal cohort. A totalVI advantage observed internally did not automatically translate into every transfer setting. When the source-target relationship was favorable for RNA-only adaptation, or when protein information did not add label-aligned signal after adaptation, scVI and totalVI performed similarly.

Direction B showed the stronger totalVI transfer result. The mean Δ+0.052 was accompanied by a seed-level pattern favoring totalVI in 19/20 repeated runs. This combination of effect size and repeated-run agreement makes Direction B more persuasive than a single split or single seed would be.

These transfer results support a context-dependent reliability view. Protein information was valuable, but the architecture-specific advantage of totalVI depended on source-target direction and adaptation context. Transfer reliability was therefore not reducible to the internal development ranking.

### Protein degradation affected training-time and post-adaptation reliability differently

Protein degradation was evaluated in two distinct designs. In the training-time design, protein entries were corrupted before totalVI training, and the resulting totalVI representations were compared with a fixed scVI RNA-only reference. In the post-adaptation design, a clean source-trained model was adapted with unlabeled query data, source and target cells were encoded with the same adapted model, a classifier was trained on clean source post-adaptation latent coordinates, and target protein values were corrupted only afterward for re-inference.

**Figure 4.** Protein degradation experiments showing training-time corruption, post-adaptation target-only corruption, Direction A/B F1 curves, ΔF1 curves, and Direction B class-level sensitivity.

Training-time totalVI macro-F1 values across corruption fractions 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85 were 0.779, 0.771, 0.758, 0.760, 0.750, 0.736, and 0.755. The scVI reference was 0.722. The non-monotonicity, especially at high corruption, was consistent with additional contribution from training stochasticity and optimization variability.

Post-adaptation degradation gave a different reliability profile. In Direction A, totalVI macro-F1 changed from 0.652 to 0.650, 0.647, 0.644, 0.641, 0.637, and 0.635 as corruption increased. The corresponding ΔF1 values were 0.0000, −0.0025, −0.0058, −0.0087, −0.0113, −0.0153, and −0.0178.

Direction B was more sensitive. totalVI macro-F1 changed from 0.654 to 0.652, 0.647, 0.643, 0.637, 0.629, and 0.621, with ΔF1 values of 0.0000, −0.0023, −0.0069, −0.0108, −0.0173, −0.0253, and −0.0328. The decline was gradual rather than abrupt, but it was stronger than in Direction A.

Class-level effects showed that global macro-F1 did not capture all vulnerability. In Direction B, CD8 Naive cells showed approximately −0.26 ΔF1 at 0.85 corruption, and CD8 TCM, CD8 TEM, CD4 TCM, CD4 TEM, and CD4 Naive populations were also sensitive. These patterns are consistent with T-cell state distinctions depending on measured surface-marker information.

### Predictive and latent uncertainty captured distinct aspects of model failure

Predictive uncertainty and latent uncertainty were intentionally separated. Predictive uncertainty was defined as 1 minus the maximum logistic-regression class probability. Latent uncertainty was defined as mean posterior variance over the 20-dimensional latent representation. Reporting both provides complementary information about classifier-level ambiguity and representation-level instability.

In PBMC transfer analyses, predictive uncertainty detected errors with AUROC values around 0.84–0.86. Under Direction B post-adaptation protein corruption, predictive AUROC declined from 0.856 to 0.841. This indicates that the classifier probability score remained informative but became less effective for ranking errors as protein measurements were degraded.

Latent uncertainty responded differently. Direction A latent error AUROC was 0.544±0.054 for scVI and 0.669±0.021 for totalVI. Direction B latent error AUROC was 0.469±0.075 for scVI and 0.610±0.038 for totalVI. In the Direction B corruption series, totalVI latent AUROC increased from 0.610 to 0.627 while predictive AUROC decreased.

Risk-coverage metrics reinforced this distinction. In Direction A, E-AURC changed from 0.0315 to 0.0339 under degradation. In Direction B, E-AURC increased from 0.0350 to 0.0448. Full AURC in Direction B increased from 0.060 to 0.079, and pAURC_[0.5,1] increased from 0.112 to 0.144. These changes show that protein corruption worsened selective-prediction behavior, especially in the more sensitive transfer direction.

The uncertainty results should be interpreted as diagnostic rather than definitive. Predictive uncertainty was stronger for immediate error ranking, while latent variance reflected representation-level instability that could increase under corrupted protein inputs. Neither score independently establishes biological correctness.

### Training stochasticity contributed materially to reliability variation

Repeated training showed that reliability variation was not explained only by modality or transfer direction. Variance partitioning separated source subset effects, model seed effects, source-by-seed interaction, and residual run-level stochastic variation under the balanced design.

**Figure 5.** Training stochasticity and selective prediction, including variance components, paired Δ components, risk-coverage summaries, E-AURC, and retention composition.

For scVI macro-F1, the source fraction was 2.31%, seed fraction was 0.00%, interaction fraction was 35.45%, and residual fraction was 62.23%. For totalVI, source accounted for 0.00%, seed for 5.62%, interaction for 3.24%, and residual for 91.13%. For the paired totalVI-minus-scVI difference, source accounted for 5.09%, seed for 12.04%, interaction for 14.11%, and residual for 68.75%.

The residual term denotes run-level residual or stochastic variation under the experimental design, not irreducible biology. This distinction matters because the balanced design used computational repetitions to estimate training and optimization variability, whereas biological replication was addressed separately through donor-level Lawlor validation.

The stochasticity results explain why small mean differences require caution. In Direction A, the near-null scVI-totalVI difference could plausibly change ordering across repeated runs. In Direction B, the larger mean difference and 19/20 seed pattern gave stronger evidence for a totalVI advantage under transfer.

### Selective prediction reduced error while changing retained cell composition

Selective prediction retained lower-uncertainty cells and deferred higher-uncertainty cells. This reduced error among retained cells, particularly when classifier predictive uncertainty was used. The benefit was not neutral, however, because retained cells were not compositionally identical to the full target population.

Retention coefficients of variation at 50% coverage were approximately 0.88–1.02 in the internal transfer setting. In Lawlor donor validation, the corresponding values were approximately 0.57–0.59. These values show that selective prediction changed retained cell composition in both settings, with stronger unevenness in the internal transfer analyses.

This compositional shift is important for biological interpretation. Rare, transitional, or marker-ambiguous populations may be deferred at higher rates if they are harder to classify. Selective prediction can therefore support review workflows, but retained-set performance should be reported together with coverage and cell-type retention.

The degradation analyses showed the limits of uncertainty-guided retention. In Direction B, full AURC increased from 0.060 to 0.079 and pAURC_[0.5,1] increased from 0.112 to 0.144 under protein corruption. Selective prediction reduced error on retained cells, but it did not fully offset the reliability loss caused by degraded protein measurements.

### Donor-level external validation confirmed the value of multimodal information

External validation used the Lawlor dataset with 10 donors and a labeled subset for evaluation [13]. Present-class donor mean macro-F1 improved by Δ+0.090 for totalVI versus scVI, with all 10 donors favoring totalVI and between-donor SD approximately 0.010. Deep-representation averages were 0.800 for scVI and 0.888 for totalVI.

**Figure 6.** Lawlor donor-level validation showing present-class donor macro-F1, donor-level ΔF1, uncertainty metrics, selective areas, and protein-aligned PCA baselines.

Across computational fold-seed evaluations, totalVI had higher macro-F1 in 50/50 comparisons. Latent AUROC favored totalVI in 38/50 comparisons. These fold-seed comparisons describe computational repeatability within the evaluation design, not biological replication; the biological units were the 10 donors.

Uncertainty and risk-coverage metrics also favored totalVI in Lawlor. Predictive error AUROC improved from 0.801 for scVI to 0.838 for totalVI, and latent error AUROC improved from 0.605 to 0.650. Full AURC was 0.087 for scVI and 0.036 for totalVI, pAURC_[0.5,1] was 0.157 versus 0.067, and E-AURC was 0.054 versus 0.026.

Simple protein-aligned baselines were very strong in Lawlor. Protein PCA reached approximately 0.927 and concat PCA reached approximately 0.945. These values indicate strong label-modality alignment in the evaluated labeled subset. They also reinforce the main hierarchy of evidence: measured protein information was highly informative, while architecture-specific conclusions remained dependent on the evaluation context.

## Discussion

This manuscript evaluated CITE-seq representation reliability across internal development data, reciprocal transfer, protein degradation, stochastic training variation, selective prediction, and donor-level validation. The strongest result was that paired protein information improved representation quality. That conclusion was supported by internal totalVI performance, strong concat PCA performance, and Lawlor donor-level gains.

The internal development cohort showed that totalVI outperformed scVI, but the concat PCA baseline prevented a purely architecture-centered interpretation. Concat PCA reached 0.745, exceeding RNA PCA and scVI, while totalVI reached 0.779. The most direct explanation is that measured protein markers supplied label-aligned information not fully captured by RNA-only representations.

Transfer analysis made the interpretation more conditional. Direction A showed a near-null totalVI-scVI difference, whereas Direction B showed Δ+0.052 with 19/20 seeds favoring totalVI. This asymmetry is important because it reflects the reality of source-target shift: composition, marker distributions, and label alignment can make one transfer direction much harder than the other.

Protein degradation experiments added another layer. Training-time corruption tested model learning with impaired protein measurements, while post-adaptation corruption tested target-side measurement loss after adaptation and classifier fitting. These designs address different questions and should not be pooled into a single robustness claim. The post-adaptation curves showed gradual declines, stronger in Direction B, and class-level sensitivity for T-cell subsets.

Uncertainty was useful but limited. Predictive uncertainty was the stronger error-ranking signal, with AUROC values around 0.84–0.86 in PBMC transfer, while latent variance captured representation-level instability and responded differently to protein degradation. These results support reporting both scores, but neither should be treated as a calibrated guarantee of biological correctness.

Selective prediction reduced retained-cell error while changing retained cell composition. This is a relevant trade-off for single-cell workflows. Deferred cells may be exactly those that require marker review or manual annotation, but retained-set summaries can become biased if coverage and cell-type retention are not reported.

The Lawlor donor analysis provided external support for the value of multimodal information. TotalVI improved present-class donor macro-F1 over scVI by Δ+0.090 across 10/10 donors, and protein PCA and concat PCA performed even higher under the evaluated protein-aligned labels. This does not establish universal totalVI superiority; it shows that when protein markers closely match label definitions, multimodal information can be highly reliable across donors.

Several limitations remain. The task was cell-type prediction, not trajectory inference, perturbation modeling, or rare-state discovery. The protein panels were limited to available measured markers, and conclusions depend on label-modality alignment. The Lawlor analysis evaluated the labeled subset and present classes per donor. Finally, the bootstrap interval was descriptive and conditional on the sampled target cells, while biological replication was evaluated separately at donor level.

## Conclusions

Protein information was the strongest contributor to improved CITE-seq representation quality in this analysis. It improved internal representation performance, supported stronger donor-level validation, and remained informative under several stress tests.

Architecture was secondary and context-dependent. totalVI showed clear advantages in the development cohort, Direction B transfer, and Lawlor donor validation, but Direction A transfer and strong concat PCA baselines showed that model rankings depended on evaluation setting.

Reliable use of CITE-seq representations therefore requires stress tests. Transfer direction, protein degradation, training stochasticity, uncertainty behavior, selective-prediction composition, and donor-level validation all contributed information that a single internal macro-F1 table would miss.

## Methods

### Datasets and cell annotations

The internal PBMC analysis used two CITE-seq datasets: PBMC10k with 6,855 cells and PBMC5k with 3,994 cells. The combined internal data contained 10,849 cells, 15,792 genes, and 14 measured proteins. A healthy-control development subset of 9,494 cells was used for internal representation comparison.

Cell-type annotations were treated as prediction targets. Macro-F1 was used for internal and transfer analyses to avoid dominance by abundant classes. The RNA-only audit compared labels used for the RNA-only analysis with reference labels in the 9,494-cell healthy-control cohort and found identical labels across 9,494 cells. This rules out a simple label-mismatch explanation but does not independently validate biological correctness.

External validation used the Lawlor dataset with 10 donors [13]. The labeled subset was used for evaluation. Donor-level analyses used present-class macro-F1 so that each donor was evaluated only over classes present in that donor's labeled subset.

### Feature processing and representation methods

RNA counts were processed by selecting 2,000 highly variable genes using the seurat_v3 procedure, normalizing counts to 1e4 per cell, applying log1p transformation, scaling with clipping at 10, and computing 10 RNA principal components. Protein features were centered log-ratio transformed, z-scored using the train or source data, and represented by 10 protein principal components.

The concat baseline concatenated 10 RNA PCs and 10 protein PCs, then applied `StandardScaler` on the train or source data. For transfer, this was an inductive source-only fit applied to target cells without query-label information. RNA PCA, protein PCA, concat PCA, MOFA+, scVI, and totalVI were compared as representation methods.

MOFA+ was included as a frozen historical embedding (`mofa_seed0`). `mofapy2` was used historically, and the frozen MOFA+ embedding with comparable latent dimension was used as a baseline without adding unsupported hyperparameter claims.

The matched scVI configuration used `n_latent=20`, `n_hidden=256`, `n_layers=2`, `dropout=0.2`, `gene_likelihood=nb`, `dispersion=gene`, and `latent_distribution=normal`. The totalVI configuration used `n_latent=20`, `n_hidden=256`, two encoder layers, one decoder layer, `dropout=0.2`, and `gene_likelihood=nb`.

Downstream classifiers used `sklearn.linear_model.LogisticRegression` with `max_iter=2000`, `solver=lbfgs`, and `random_state=seed`. The software environment used Python 3.11.7, scikit-learn 1.5.2, and scvi-tools 1.3.3.

### Internal development-cohort evaluation

The healthy-control development cohort contained 9,494 cells. It was used to compare RNA PCA, protein PCA, concat PCA, MOFA+, scVI, and totalVI using the shared logistic-regression classifier.

Internal macro-F1 values were RNA PCA 0.666, protein PCA 0.512, concat PCA 0.745, MOFA+ 0.549, scVI 0.722, and totalVI 0.779. The totalVI-minus-scVI difference was +0.057, and the totalVI-minus-concat difference was +0.034. The internal totalVI-minus-scVI interval used a target-cell paired bootstrap with n=500 resamples; the interval is descriptive and conditional on the evaluated target cells. Biological replication was addressed separately by Lawlor donor-level validation.

### Bidirectional transfer and scArches adaptation

Bidirectional transfer reversed which internal PBMC dataset served as source and which served as target. Direction A produced scVI macro-F1 0.656±0.018 and totalVI macro-F1 0.652±0.015. Direction B produced scVI macro-F1 0.601±0.020 and totalVI macro-F1 0.654±0.025, for Δ+0.052, with 19/20 seeds favoring totalVI. The inductive concat baseline reached 0.735 in Direction A and 0.734 in Direction B.

Generative-model transfer used unsupervised transductive scArches query adaptation [10]. A clean source model was adapted using unlabeled query data. After adaptation, source and target cells were encoded with the same adapted model. The classifier was trained on clean source post-adaptation latent coordinates and evaluated on target cells. This procedure used unlabeled target data during adaptation and was therefore not fully inductive.

### Protein degradation experiments

Protein corruption was entry-wise. At each nonzero corruption fraction, originally non-zero measured protein entries were randomly selected and set to zero; preexisting zeros were unchanged. Corruption fractions were 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85, with five mask seeds at nonzero fractions.

Training-time degradation retrained totalVI after corrupting protein measurements, while the scVI RNA-only reference remained fixed. Training-time totalVI macro-F1 values were 0.779, 0.771, 0.758, 0.760, 0.750, 0.736, and 0.755 across increasing corruption fractions; the fixed scVI reference was 0.722.

Post-adaptation degradation used a clean source training procedure followed by unsupervised transductive scArches query adaptation. Source and target cells were encoded with the same adapted model, the classifier was trained on clean source post-adaptation latent coordinates and frozen, and then target protein values only were corrupted and re-inferred. Source classifier training remained clean, and the procedure was not fully inductive.

Post-adaptation Direction A macro-F1 values were 0.652, 0.650, 0.647, 0.644, 0.641, 0.637, and 0.635, with ΔF1 values 0.0000, −0.0025, −0.0058, −0.0087, −0.0113, −0.0153, and −0.0178. Direction B macro-F1 values were 0.654, 0.652, 0.647, 0.643, 0.637, 0.629, and 0.621, with ΔF1 values 0.0000, −0.0023, −0.0069, −0.0108, −0.0173, −0.0253, and −0.0328. In Direction B, CD8 Naive cells showed approximately −0.26 ΔF1 at 0.85 corruption.

### Uncertainty and selective prediction

Predictive uncertainty was defined as `U = 1 - max p`, where `p` is the logistic-regression class-probability vector. Latent uncertainty was defined as the mean posterior variance over the 20 latent dimensions. Predictive uncertainty was treated as classifier-level ambiguity, whereas latent uncertainty was treated as representation-level instability.

Predictive error-detection AUROC was approximately 0.84–0.86 in PBMC transfer analyses. Direction A latent AUROC values were 0.544±0.054 for scVI and 0.669±0.021 for totalVI. Direction B latent AUROC values were 0.469±0.075 for scVI and 0.610±0.038 for totalVI. In Direction B post-adaptation degradation, predictive AUROC changed from 0.856 to 0.841, while latent AUROC changed from 0.610 to 0.627.

Risk-coverage analysis evaluated error among retained cells as coverage varied. Historical `aurc` was treated as full AURC. Partial AURC over coverage 0.5 to 1.0 was reported as pAURC_[0.5,1], and excess AURC was reported separately as E-AURC. Under degradation, E-AURC changed from 0.0315 to 0.0339 in Direction A and from 0.0350 to 0.0448 in Direction B. Direction B full AURC changed from 0.060 to 0.079, and pAURC_[0.5,1] changed from 0.112 to 0.144.

Selective prediction retained lower-uncertainty cells and deferred higher-uncertainty cells. Cell-type retention was summarized at 50% coverage. Retention coefficients of variation were approximately 0.88–1.02 internally and approximately 0.57–0.59 in Lawlor donor validation.

### Training stochasticity and variance decomposition

Training stochasticity was evaluated with a balanced design using 5 source subsets, 5 model seeds, and 2 run replicates. Per-model method-of-moments estimates were derived from expected mean squares. A separate paired Delta cube was used for the totalVI-minus-scVI contrast.

Variance components were summarized as source, seed, source-by-seed interaction, and residual. For scVI macro-F1, the fractions were source 2.31%, seed 0.00%, interaction 35.45%, and residual 62.23%. For totalVI, the fractions were source 0.00%, seed 5.62%, interaction 3.24%, and residual 91.13%. For totalVI minus scVI, the fractions were source 5.09%, seed 12.04%, interaction 14.11%, and residual 68.75%.

The residual component denotes run-level residual or stochastic variation under the design, not irreducible biology. Donor-level Lawlor validation was used to address biological replication separately.

### Lawlor donor-level validation

Lawlor donor-level validation used 10 donors and evaluated the labeled subset [13]. The primary donor-level metric was present-class macro-F1. TotalVI improved over scVI by Δ+0.090 on average, with 10/10 donors favoring totalVI and between-donor SD approximately 0.010. Deep-representation averages were 0.800 for scVI and 0.888 for totalVI.

Across computational fold-seed evaluations, totalVI had higher macro-F1 in 50/50 comparisons. Latent AUROC favored totalVI in 38/50 comparisons. These repeated fold-seed evaluations were computational comparisons, not donor-level replication.

Predictive AUROC was 0.801 for scVI and 0.838 for totalVI. Latent AUROC was 0.605 for scVI and 0.650 for totalVI. Full AURC was 0.087 for scVI and 0.036 for totalVI, pAURC_[0.5,1] was 0.157 versus 0.067, and E-AURC was 0.054 versus 0.026. Protein PCA reached approximately 0.927, and concat PCA reached approximately 0.945.

### Software and reproducibility

Analyses were conducted in Python 3.11.7. scVI and totalVI models used scvi-tools 1.3.3. Logistic-regression classifiers used scikit-learn 1.5.2 with `max_iter=2000`, `solver=lbfgs`, and `random_state=seed`.

Analysis code is available upon request.

## Figure Legends

Figure 1. Study design and evaluation units. Schematic of the analysis spanning development PBMC CITE-seq cohorts, representation families, reciprocal transfer, training-time and post-adaptation protein degradation, uncertainty and selective prediction, variance decomposition, and Lawlor 10-donor external validation. Biological donors are distinguished from computational repetitions such as seeds, folds, mask seeds, and run replicates.

Figure 2. Information content versus representation strategy in the development high-confidence cohort. Macro-F1 is shown for RNA PCA (0.666), protein PCA (0.512), concat PCA (0.745), MOFA+ (0.549), scVI (0.722), and totalVI (0.779) on n=9,494 healthy-control cells. Residual differences are totalVI−scVI (+0.057) and totalVI−concat (+0.034). The RNA-only audit showed identical labels across 9,494 cells; the totalVI−scVI interval used a cell-level conditional paired bootstrap (95% CI 0.020–0.091).

Figure 3. Dataset transfer and uncertainty. Reciprocal transfer is shown for Directions A and B, contrasting unsupervised transductive scArches adaptation with inductive source-fitted concat PCA. Panels summarize 20-seed mean±SD macro-F1, inductive concat performance (0.735/0.734), predictive error AUROC, latent error AUROC, and uncertainty definitions. The 19/20 seed pattern applies only to PBMC Direction B.

Figure 4. Protein degradation experiments. Training-time and post-adaptation designs are shown separately. Post-adaptation curves summarize totalVI macro-F1 and ΔF1 from clean for Directions A and B across corruption fractions 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85, with five mask seeds at nonzero fractions. Direction B class-level sensitivity includes CD8 Naive ΔF1 ≈−0.26 at 0.85 corruption.

Figure 5. Training stochasticity and selective prediction. Variance components are shown for within-model macro-F1 and paired totalVI−scVI Δ, with residual denoting run-level residual or stochastic variation under the balanced design. Risk-coverage summaries distinguish full AURC, pAURC_[0.5,1], and E-AURC. Retention coefficients of variation at 50% coverage summarize changes in retained cell composition.

Figure 6. Lawlor donor-level external validation. The Lawlor analysis used 10 biological donors and a labeled subset for evaluation. Panels summarize per-donor present-class macro-F1, donor-level ΔF1 (10/10 donors positive; mean Δ≈+0.090; between-donor SD≈0.010), protein PCA (~0.927), concat PCA (~0.945), predictive and latent AUROC, and selective risk-coverage areas. Fold-seed comparisons are computational repetitions, not donor-level replication.

## Tables

### Table 1. Internal development-cohort representation performance

| Representation | Modalities | Macro-F1 |
|---|---:|---:|
| RNA PCA | RNA | 0.666 |
| Protein PCA | Protein | 0.512 |
| Concatenated PCA | RNA + protein | 0.745 |
| MOFA+ | RNA + protein | 0.549 |
| scVI | RNA | 0.722 |
| totalVI | RNA + protein | 0.779 |

### Table 2. Bidirectional transfer performance

| Direction | scVI macro-F1 | totalVI macro-F1 | Difference | Seed pattern | Inductive concat macro-F1 |
|---|---:|---:|---:|---:|---:|
| A | 0.656±0.018 | 0.652±0.015 | -0.004 | near null | 0.735 |
| B | 0.601±0.020 | 0.654±0.025 | +0.052 | 19/20 favor totalVI | 0.734 |

### Table 3. Post-adaptation protein degradation in totalVI

| Corruption fraction | Direction A macro-F1 | Direction A Δ from clean | Direction B macro-F1 | Direction B Δ from clean |
|---:|---:|---:|---:|---:|
| 0.00 | 0.652 | 0.0000 | 0.654 | 0.0000 |
| 0.10 | 0.650 | -0.0025 | 0.652 | -0.0023 |
| 0.25 | 0.647 | -0.0058 | 0.647 | -0.0069 |
| 0.40 | 0.644 | -0.0087 | 0.643 | -0.0108 |
| 0.55 | 0.641 | -0.0113 | 0.637 | -0.0173 |
| 0.70 | 0.637 | -0.0153 | 0.629 | -0.0253 |
| 0.85 | 0.635 | -0.0178 | 0.621 | -0.0328 |

### Table 4. Training-time protein degradation

| Corruption fraction | totalVI macro-F1 | scVI reference (fixed) |
|---:|---:|---:|
| 0.00 | 0.779 | 0.722 |
| 0.10 | 0.771 | 0.722 |
| 0.25 | 0.758 | 0.722 |
| 0.40 | 0.760 | 0.722 |
| 0.55 | 0.750 | 0.722 |
| 0.70 | 0.736 | 0.722 |
| 0.85 | 0.755 | 0.722 |

### Table 5. Variance partitioning of macro-F1

| Outcome | Source | Seed | Interaction | Residual |
|---|---:|---:|---:|---:|
| scVI | 2.31% | 0.00% | 35.45% | 62.23% |
| totalVI | 0.00% | 5.62% | 3.24% | 91.13% |
| totalVI - scVI | 5.09% | 12.04% | 14.11% | 68.75% |

### Table 6. Lawlor donor-level validation

| Metric | scVI | totalVI | Interpretation |
|---|---:|---:|---|
| Deep average macro-F1 | 0.800 | 0.888 | Multimodal improvement |
| Present-class donor mean delta | NA | +0.090 | 10/10 donors favor totalVI |
| Between-donor SD of delta | NA | ~0.010 | Donor-level consistency |
| Macro-F1 fold-seed wins | NA | 50/50 | Computational comparisons |
| Predictive AUROC | 0.801 | 0.838 | Better error ranking |
| Latent AUROC | 0.605 | 0.650 | Better latent uncertainty signal |
| Latent AUROC wins | NA | 38/50 | Computational comparisons |
| Full AURC | 0.087 | 0.036 | Lower is better |
| pAURC_[0.5,1] | 0.157 | 0.067 | Lower is better |
| E-AURC | 0.054 | 0.026 | Lower is better |
| Protein PCA macro-F1 | ~0.927 | NA | Strong label-protein alignment |
| Concatenated PCA macro-F1 | ~0.945 | NA | Strong label-modality alignment |

## References

1. Stoeckius M, Hafemeister C, Stephenson W, Houck-Loomis B, Chattopadhyay PK, Swerdlow H, et al. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods. 2017;14:865-868. doi:10.1038/nmeth.4380.
2. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15:1053-1058. doi:10.1038/s41592-018-0229-2.
3. Gayoso A, Steier Z, Lopez R, Regier J, Nazor KL, Streets A, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods. 2021;18:272-282. doi:10.1038/s41592-020-01050-x.
4. Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol. 2020;21:111. doi:10.1186/s13059-020-02015-1.
5. Ashuach T, Gabitto MI, Koodli RV, Saldi GA, Jordan MI, Yosef N. MultiVI: deep generative model for the integration of multimodal data. Nat Methods. 2023;20:1222–1231. doi:10.1038/s41592-023-01909-9.
6. Gong B, Zhou Y, Purdom E. Cobolt: integrative analysis of multimodal single-cell sequencing data. Genome Biol. 2021;22:351. doi:10.1186/s13059-021-02556-z.
7. Jeong Y, Ronen J, Kopp W, Lutsik P, Akalin A. scMaui: a widely applicable deep learning framework for single-cell multiomics integration in the presence of batch effects and missing data. BMC Bioinformatics. 2024;25:257. doi:10.1186/s12859-024-05880-w.
8. Wu M, Goodman N. Multimodal generative models for scalable weakly-supervised learning. Adv Neural Inf Process Syst. 2018;31:5580–5590.
9. Luecken MD, Büttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colomé-Tatché M, Theis FJ. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods. 2022;19:41–50. doi:10.1038/s41592-021-01336-8.
10. Lotfollahi M, Naghipourfar M, Luecken MD, Khajavi M, Büttner M, Wagenstetter M, Avsec Ž, Gayoso A, Yosef N, Interlandi M, et al. Mapping single-cell data to reference atlases by transfer learning. Nat Biotechnol. 2022;40:121–130. doi:10.1038/s41587-021-01001-7.
11. Guo C, Pleiss G, Sun Y, Weinberger KQ. On calibration of modern neural networks. Proc Mach Learn Res. 2017;70:1321–1330.
12. Geifman Y, El-Yaniv R. Selective classification for deep neural networks. Adv Neural Inf Process Syst. 2017;30:4878–4887.
13. Lawlor N, Nehar-Belaid D, Grassmann JDS, Stoeckius M, Smibert P, Stitzel ML, Pascual V, Banchereau J, Williams A, Ucar D. Single Cell Analysis of Blood Mononuclear Cells Stimulated Through Either LPS or Anti-CD3 and Anti-CD28. Front Immunol. 2021;12:636720. doi:10.3389/fimmu.2021.636720.
14. Hao Y, Hao S, Andersen-Nissen E, Mauck WM III, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184:3573–3587.e29. doi:10.1016/j.cell.2021.04.048.
15. Naeini MP, Cooper GF, Hauskrecht M. Obtaining well calibrated probabilities using Bayesian binning. AAAI. 2015;29:2901-2907. doi:10.1609/aaai.v29i1.9602.
16. Brier GW. Verification of forecasts expressed in terms of probability. Mon Weather Rev. 1950;78:1-3. doi:10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2.
