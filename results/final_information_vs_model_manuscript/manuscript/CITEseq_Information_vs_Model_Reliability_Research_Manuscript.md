<!--
Alternative titles:
1. Separating multimodal information gain from model-specific advantage in CITE-seq representation learning
2. Protein information and model-specific reliability in CITE-seq representation learning
-->

# Disentangling protein-information gain and model-specific reliability in CITE-seq representation learning

**Jizhu Yang**  
Research Manuscript  
Unpublished manuscript  
September 2026

## Abstract

**Background:** Multimodal CITE-seq models often improve cell-type representations, but gains can confound access to protein measurements with modeling architecture. We asked how much multimodal gain reflects protein information itself, and when a complex model adds reliable value beyond a simple multimodal representation. **Results:** In a 9,494-cell development cohort, RNA PCA, concat PCA, and totalVI reached macro-F1 0.666, 0.745, and 0.779. Descriptive contrasts were \(G_{\mathrm{info}}=+0.079\) and \(G_{\mathrm{model}}=+0.034\); concat recovered about 70% of the RNA-to-totalVI gap without causal attribution. Source-only inductive concat remained strong under transfer (0.735, 0.734). Target-feature-access-matched transductive concat reached 0.636 (Direction A) and 0.767 (Direction B). Relative to that matched control, totalVI was about +0.016 in Direction A and −0.113 in Direction B (Case C), so matched target-data access did not yield a stable totalVI advantage. Direction B totalVI–scVI was +0.052 (19/20 seeds); Direction A was near null. Post-adaptation protein degradation narrowed multimodal/model advantage gradually rather than catastrophically. Lawlor validation across 10 donors supported multimodal information value under protein-gated labels, where protein PCA and concat exceeded totalVI. **Conclusions:** Protein-information gain was strong; residual model-specific gain was smaller and context-dependent.

## Background

CITE-seq jointly measures transcriptomes and antibody-derived surface proteins in the same single cells [1]. This paired measurement is attractive for immune profiling because many cell labels are informed by both transcriptional programs and surface-marker phenotypes. It also creates a central evaluation problem: a model that improves cell-type prediction may be benefiting from access to protein measurements, from a better representation architecture, or from both.

Deep generative models have become common tools for single-cell representation learning. scVI models transcript counts with a variational autoencoder and provides an RNA-only latent representation [2]. totalVI extends this framework to paired RNA and protein measurements and learns a joint latent representation for CITE-seq data [3]. Related approaches include MOFA+ factor analysis [4], MultiVI [5], Cobolt [6], scMaui [7], weakly supervised multimodal generative modeling [8], atlas-level integration benchmarking [9], scArches query mapping [10], and WNN-style multimodal integration [14].

The key gap addressed here is that comparisons between RNA-only and multimodal models often confound modality access with architecture. A contrast such as totalVI versus scVI answers whether a particular multimodal model outperforms a particular RNA-only model, but it does not say how much of the gain could be obtained by a simple representation that merely has access to RNA and protein features. A simple concatenated RNA-protein PCA baseline is therefore not a nuisance control; it is essential for separating protein-information gain from residual model-specific advantage.

This distinction matters for reliability. If most of the observed gain comes from protein information, then the core scientific conclusion is that paired surface markers are informative. If a residual gain remains after comparing against a simple multimodal baseline, then that residual should be stress-tested under transfer, protein degradation, repeated training, uncertainty diagnostics, selective prediction, and donor-level validation. The central question is therefore: how much of multimodal CITE-seq gain comes from access to protein information itself, and when does a complex multimodal model provide reliable additional value beyond a simple multimodal representation?

Transfer evaluation is especially sensitive to target-data access. scArches adaptation uses unlabeled query features during model adaptation [10], whereas an inductive PCA baseline can be fit only on the source dataset and applied directly to target cells. These are both useful practical settings, but they are not matched in their access to unlabeled target features. A matched transfer control should therefore compare against a transductive concat baseline fit on source plus unlabeled target features, while clearly stating that this matches feature access only and is not algorithmically identical to scArches.

Uncertainty and selective prediction add another layer of reliability assessment. Calibration studies show that probability estimates can be miscalibrated even when classification accuracy is high [11,15], and selective classification formalizes the trade-off between retained-set risk and coverage [12]. In single-cell biology, selective prediction must also track retained cell-type composition, because deferring uncertain cells can preferentially remove rare, transitional, or marker-ambiguous populations. Proper scoring rules such as Brier score provide additional probability-quality diagnostics [16]. These diagnostics are used here as stress tests of when the residual model-specific advantage may fail, not as a separate research program.

This study first separates multimodal information gain from model-specific gain, then tests the stability of the remaining model-specific advantage under transfer, degradation, stochastic training variation, and external donor validation. The Lawlor blood CITE-seq dataset provides a 10-donor external setting with protein-gated author labels and paired RNA-protein measurements [13]. Lawlor is used to test whether multimodal information remains valuable across donors under protein-gated labels, not to claim that one architecture is superior in every possible CITE-seq application.

## Results

**Figure 1.** Central logic: information gain versus model-specific advantage.

### Simple multimodal information explained much of the internal gain

The internal development analysis used 9,494 healthy-control cells with matched RNA and protein measurements. The primary contrast was deliberately descriptive. RNA PCA reached macro-F1 0.666, concatenated RNA-protein PCA reached 0.745, and totalVI reached 0.779. This gives:

G_info = F1_concat - F1_RNA_PCA = 0.745 - 0.666 = +0.079  
G_model = F1_totalVI - F1_concat = 0.779 - 0.745 = +0.034

The full RNA-to-totalVI gap was +0.113. The concat baseline therefore recovered approximately 70% of that gap descriptively, without implying that 70% is a mechanistic partition or explanatory fraction. The more limited interpretation is that a large part of the internal multimodal improvement was already obtained by giving a simple representation access to both RNA and protein features.

**Figure 2.** Information-versus-model comparison in the development cohort.

Secondary representations reinforced the same point. Protein PCA alone reached 0.512, MOFA+ reached 0.549, and scVI reached 0.722. Protein PCA did not replace RNA, but protein features were useful when paired with transcriptomic features. Among tested methods, totalVI reached the highest internal macro-F1 (0.779), yet its residual advantage over concat (\(G_{\mathrm{model}}=+0.034\)) was smaller than the descriptive protein-information increment (\(G_{\mathrm{info}}=+0.079\)).

The RNA-only label audit ruled out a simple label-provenance explanation for the internal gain. Development labels were already derived through an RNA-only reference pathway, and agreement across the 9,494 high-confidence cells was 1.0. This audit does not independently validate biological labels, but it supports the internal comparison by showing that the multimodal advantage was not a trivial consequence of protein-informed labels.

### The residual model-specific advantage was not stable across transfer settings

Bidirectional transfer tested whether the residual model-specific advantage persisted under source–target shift. Two simple multimodal baselines were retained so that target-data access could be interpreted separately from architecture.

**Source-fitted inductive concat PCA** (source only): Direction A 0.735; Direction B 0.734.

**Target-data-access-matched transductive concat PCA** (source + unlabeled target features; not algorithmically identical to scArches): Direction A 0.636; Direction B 0.767.

**scArches-adapted VAEs:** scVI reached 0.656±0.018 (A) and 0.601±0.020 (B); totalVI reached 0.652±0.015 (A) and 0.654±0.025 (B). Direction A totalVI–scVI was near null; Direction B gave Δ = +0.052 with 19/20 seeds favoring totalVI.

**Figure 3.** Transfer fairness comparison across inductive concat, transductive concat, scVI, and totalVI.

Relative to the matched transductive concat control (**Case C**):

- Direction A: totalVI − transductive concat ≈ +0.016 (small positive).
- Direction B: totalVI − transductive concat ≈ −0.113 (transductive concat substantially higher).

Thus matching unlabeled-target *feature access* did not produce a uniform totalVI advantage. Unlabeled target-data access itself was also not uniformly beneficial: joint unsupervised concat fitting lowered performance in Direction A (0.735 → 0.636) but raised it in Direction B (0.734 → 0.767). These patterns are direction- and strategy-dependent; they do not imply a general rule that transduction outperforms induction, nor that scArches “failed.”

The Direction B totalVI–scVI contrast shows that protein-aware VAE modeling can help under one transfer design. Case C shows that, once unlabeled target feature access is matched with a simple multimodal representation, residual totalVI advantage is not stable across directions. Transfer therefore supports multimodal information value more strongly than a general architecture ranking.

### Protein degradation progressively narrowed multimodal model advantage

Protein degradation asked what happens to residual multimodal/model advantage when protein signal quality deteriorates. Two designs were used. In training-time degradation, protein entries were corrupted before totalVI training and totalVI was retrained. In post-adaptation degradation, source training and scArches adaptation were clean, the classifier was trained on clean source post-adaptation latents, and target protein values were corrupted only afterward before re-inference.

**Figure 4.** Protein degradation experiments for training-time corruption and post-adaptation target-only corruption.

Training-time totalVI macro-F1 across corruption fractions 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85 was 0.779, 0.771, 0.758, 0.760, 0.750, 0.736, and 0.755, compared with a fixed scVI reference of 0.722. Performance declined gradually rather than catastrophically; the curve remained above the RNA-only scVI reference at the measured levels but was not monotonic, consistent with training and optimization variation in addition to measurement degradation.

Post-adaptation target-only degradation produced clearer progressive loss. In Direction A, totalVI macro-F1 changed from 0.652 to 0.650, 0.647, 0.644, 0.641, 0.637, and 0.635 as corruption increased (Δ from clean to −0.0178 at p=0.85). Direction A totalVI was already near or slightly below clean scVI (0.656) at baseline, so residual model advantage was limited even before corruption. In Direction B, totalVI macro-F1 changed from 0.654 to 0.652, 0.647, 0.643, 0.637, 0.629, and 0.621 (Δ to −0.0328 at p=0.85). Relative to the clean scVI reference (0.601), Direction B totalVI advantage narrowed as protein quality fell but remained positive at p=0.85.

The stronger Direction B decline is important because Direction B was also the direction with the clearest totalVI–scVI transfer advantage. As target protein measurements were degraded, that multimodal model advantage narrowed. At the class level, Direction B CD8 Naive cells showed approximately −0.26 F1 change at 0.85 corruption, and several CD4/CD8 T-cell subsets were also sensitive, consistent with surface-marker dependence for fine immune-state distinctions.

### Uncertainty diagnostics identified distinct forms of reliability loss

Uncertainty diagnostics were subordinated to the central question: when residual multimodal/model advantage becomes less reliable, do uncertainty signals reveal that degradation? Classifier predictive uncertainty and latent posterior variance measured different failure modes. Classifier predictive uncertainty was defined as \(U = 1 - \max_k p_k\), where \(p\) is the logistic-regression class-probability vector. Latent posterior variance was defined as the mean posterior variance over latent dimensions and therefore measured representation-level dispersion, not epistemic uncertainty in a formal Bayesian sense.

**Figure 5.** Uncertainty diagnostics, training stochasticity, and selective prediction.

In PBMC transfer, classifier predictive uncertainty detected classification errors with AUROC values around 0.84–0.86. In Direction B post-adaptation corruption, predictive AUROC declined mildly from 0.856 to 0.841, showing that the classifier confidence score remained informative but became less effective for ranking errors as protein inputs were degraded.

Latent posterior variance responded differently. Direction A latent error AUROC was 0.544±0.054 for scVI and 0.669±0.021 for totalVI. Direction B latent error AUROC was 0.469±0.075 for scVI and 0.610±0.038 for totalVI. Under Direction B corruption, totalVI latent AUROC increased from 0.610 to 0.627 while predictive AUROC decreased. A rise in latent AUROC does not mean that “uncertainty improved”; it indicates that posterior-variance ranking of errors changed under degradation, whereas classifier predictive uncertainty and selective-risk metrics moved in the opposite direction.

Risk-coverage summaries showed the practical consequence. In Direction A, E-AURC changed from 0.0315 to 0.0339 under degradation. In Direction B, E-AURC increased from 0.0350 to 0.0448. Full AURC in Direction B increased from 0.060 to 0.079, and \(\mathrm{pAURC}_{[0.5,1.0]}\) increased from 0.112 to 0.144. Negative log-likelihood and Brier score also worsened under degradation, whereas ECE remained roughly stable. These metrics therefore measure different properties of reliability loss.

### Training stochasticity limited the precision of small model advantages

Training stochasticity was connected directly to the central information-versus-model question: small residual model-specific gains must be interpreted relative to training variability. The balanced stochasticity design separated source subset effects, model seed effects, source-by-seed interaction, and residual run-level variation using method-of-moments estimates from expected mean squares.

For scVI macro-F1, the variance fractions were source 2.31%, seed 0.00%, interaction 35.45%, and residual 62.23%. For totalVI, the fractions were source 0.00%, seed 5.62%, interaction 3.24%, and residual 91.13%. For the paired Δ = totalVI − scVI difference, the fractions were source 5.09%, seed 12.04%, interaction 14.11%, and residual 68.75%.

The residual component denotes run-level residual or stochastic variation under the implemented replication design, not irreducible biological noise. This distinction matters because computational repeats estimate training and optimization variability, whereas biological replication was assessed separately using donor-level external validation.

These results explain why Direction A should be interpreted as near null and why the small internal \(G_{\mathrm{model}}\) estimate should be treated cautiously. The Direction B totalVI–scVI transfer difference was more persuasive because it combined a larger mean difference with 19/20 seed-level agreement, but the matched transfer control still showed that this was not a stable residual model advantage once unlabeled target feature access was matched.

### Selective prediction reduced error but introduced cell-type-dependent retention

Selective prediction retained low-uncertainty cells and deferred high-uncertainty cells. This reduced error among retained cells, especially when predictive uncertainty was used, but it also changed the retained cell-type composition.

At 50% coverage, retention coefficients of variation were approximately 0.88-1.02 in the internal transfer setting and approximately 0.57-0.59 in Lawlor donor validation. These values indicate that retained cells were not a compositionally neutral subset of the evaluated cells.

The biological implication is that retained-set performance is incomplete without coverage and cell-type retention. Rare, transitional, or marker-ambiguous populations may be preferentially deferred. Selective prediction can therefore support review workflows, but it should be reported as a risk-coverage-composition trade-off rather than as a simple improvement in classification.

Protein degradation also limited selective prediction. In Direction B, full AURC increased from 0.060 to 0.079 and \(\mathrm{pAURC}_{[0.5,1.0]}\) increased from 0.112 to 0.144 under target-only protein corruption. Retaining low-uncertainty cells reduced error, but it did not eliminate reliability loss caused by degraded protein measurements.

### External donor validation confirmed multimodal information value under protein-gated labels

The Lawlor blood CITE-seq validation used 10 biological donors as the primary units and author labels in the labeled subset [13]. Present-class donor-level macro-F1 improved by Δ ≈ +0.090 for totalVI versus scVI, with all 10 donors favoring totalVI and between-donor SD approximately 0.010. Deep-representation averages were ≈0.800 for scVI and ≈0.888 for totalVI.

**Figure 6.** Lawlor donor-level external validation.

Classifier predictive and latent posterior-variance diagnostics also differed between models in Lawlor (predictive error AUROC 0.801 vs 0.838; latent error AUROC 0.605 vs 0.650; full AURC 0.087 vs 0.036; \(\mathrm{pAURC}_{[0.5,1.0]}\) 0.157 vs 0.067; E-AURC 0.054 vs 0.026). These are reliability diagnostics under the same protein-gated labels, not independent proof of architecture superiority.

Simple protein-aligned baselines were very strong in the same external setting. Protein PCA reached approximately 0.927 and concat PCA reached approximately 0.945, exceeding totalVI (≈0.888). Because Lawlor labels were defined using protein markers, the strong protein-PCA and concatenated-PCA performance is consistent with substantial label–modality alignment. The cohort therefore provides evidence that protein-aligned multimodal information remains useful across donors, rather than evidence of universal multimodal-architecture superiority.

## Discussion

Across the evaluated PBMC CITE-seq settings, access to paired protein measurements provided the most reproducible multimodal benefit, whereas the additional advantage associated with totalVI relative to a simple RNA–protein representation was smaller and strongly dependent on transfer context. This manuscript therefore separates two quantities that are often conflated: protein-information gain and residual model-specific gain. The development cohort showed a clear multimodal improvement, but concatenated RNA–protein PCA captured much of the RNA-to-totalVI gap, shifting interpretation away from an architecture-only story.

The descriptive contrasts make that point transparent. \(G_{\mathrm{info}}\) was +0.079, while \(G_{\mathrm{model}}\) was +0.034. Concat recovered about 70% of the RNA-to-totalVI gap descriptively, but this fraction should not be read as a causal attribution. The model-specific contrast is descriptive and should not be interpreted as a causal decomposition of architecture effects: representations, objectives, adaptation procedures, and preprocessing assumptions differ across methods.

Transfer analysis showed why this distinction matters. Source-only inductive concat remained strong in both directions (0.735 / 0.734). In the scArches-adapted VAE comparison, Direction A was near null and Direction B favored totalVI over scVI by +0.052 (19/20 seeds). These observations support context-specific totalVI value versus scVI in Direction B, but they do not support a broad model-ranking statement.

The matched transfer control was the most important additional check. Because scArches adaptation uses unlabeled target features, an inductive source-only concat baseline is not matched for target-data access. Transductive concat matched unlabeled-target feature access only, while remaining algorithmically different from scArches. Under Case C, totalVI was only ≈+0.016 above transductive concat in Direction A but ≈−0.113 below it in Direction B. Matching target-feature access therefore changes the central interpretation: residual totalVI advantage is not stable once a simple multimodal representation is given comparable unlabeled-target access.

Protein degradation supported the same reliability framing. Training-time performance declined gradually rather than catastrophically. Post-adaptation target-only degradation showed progressive loss, stronger in Direction B, where totalVI’s advantage over clean scVI narrowed as protein quality fell (0.654 → 0.621 at p=0.85; clean scVI 0.601). Degradation therefore erodes residual multimodal/model advantage as signal quality falls; it does not imply that totalVI is simply “robust to corruption.”

Uncertainty diagnostics were useful because macro-F1 alone hid different forms of failure. Classifier predictive uncertainty ranked errors well and declined mildly under protein corruption, while latent posterior variance tracked representation-level dispersion differently. A rise in latent AUROC under degradation should not be read as uncertainty “improving.” Neither signal guarantees that labels are biologically correct.

Training stochasticity further limited the precision of small advantages. Method-of-moments variance decomposition showed substantial run-level residual/stochastic variability under the implemented replication design, especially for paired contrasts. Small differences such as \(G_{\mathrm{model}}\) therefore require repeated-run context and matched baselines before being interpreted as reliable added value.

Selective prediction illustrates an operational trade-off. Uncertainty-based retention reduced retained-set error, but retention was cell-type dependent. Reporting risk without cell-type retention would overstate practical reliability.

Lawlor donor validation confirmed multimodal information value across 10 biological donors under protein-gated author labels. totalVI improved over scVI by about +0.090 in present-class donor macro-F1, and protein PCA / concat PCA were even stronger (≈0.927 / ≈0.945 vs totalVI ≈0.888). Because labels were defined using protein markers, strong protein-aligned baselines are consistent with label–modality alignment. Donor consistency supports biological reproducibility of multimodal information value; it does not establish universal multimodal-architecture superiority.

Several limitations remain. The analysis is PBMC-focused with limited antibody panels and a cell-type classification focus rather than trajectories or rare-state discovery. Lawlor labels are protein-gated, creating external protein–label alignment. Deep-model transfer uses unsupervised scArches transduction, whereas the matched concat control matches unlabeled-target feature access but not the adaptation algorithm or objective; concat and scArches are not algorithmically identical. Simple multimodal baselines remain competitive. There is no independent external cohort with protein-independent labels, no arbitrary cell-wise missing-protein benchmark, and no universal comparison across multimodal architectures. A useful future validation would involve an independent CITE-seq cohort with labels not directly defined by protein gating; that experiment is not required for the present application manuscript. These limits make the conclusion narrower but more reliable: protein information contributes substantially, while residual model-specific advantage must be demonstrated within each evaluation context.

## Conclusion

Paired protein measurements accounted for a substantial part of the performance improvement observed in CITE-seq representation learning. totalVI provided additional gains in selected settings, but these gains were smaller and dependent on source-target context, protein quality, and computational stochasticity. Reliability claims for multimodal models should therefore be evaluated against strong simple multimodal baselines and interpreted separately from the information value of the additional modality itself.

## Methods

### Datasets and cell annotations

The internal PBMC analysis used public PBMC10k and PBMC5k CITE-seq datasets. After quality control and feature harmonization, the combined internal data contained 10,849 cells, 15,792 shared genes, and 14 shared proteins. The development analysis used 9,494 healthy-control high-confidence cells.

Cell-type annotations were treated as prediction targets. Macro-F1 was the primary metric for internal and transfer analyses to reduce dominance by abundant cell types. An RNA-only label audit compared the development labels with RNA-only reference labels and found agreement of 1.0 across the 9,494 high-confidence cells.

External validation used the Lawlor blood CITE-seq dataset with 10 donors and author labels [13]. Author labels are protein-gated: annotation depends substantially on surface-protein phenotypes. The labeled subset was used for classification evaluation. Donor-level summaries used present-class macro-F1 so that each donor was evaluated only over classes present in that donor's labeled subset. Source-only preprocessing and `run_identifier` as batch were retained. Latent consistency of totalVI over scVI across matched fold×seed computational pairs was 38/50; these pairs are computational comparisons, not biological replicates.

### Feature processing and representation methods

RNA PCA features were built from 2,000 highly variable genes selected with the seurat_v3 procedure, followed by normalization to 1e4 counts per cell, log1p transformation, scaling with clipping at 10, and PCA. Protein features were centered log-ratio transformed, z-scored using training or source data, and represented by protein PCA.

The primary simple multimodal representation concatenated 10 RNA PCs and 10 protein PCs, then applied `StandardScaler` on the train or source data. In inductive transfer, all PCA and scaling fits used source data only and were applied to target cells. In the matched transductive transfer control, the concat representation was fit on source plus unlabeled target features, matching unlabeled-target FEATURE ACCESS only and not attempting to reproduce the scArches algorithm.

MOFA+ was included as a frozen multimodal factor-analysis baseline. scVI used matched capacity settings with `n_latent=20`, `n_hidden=256`, `n_layers=2`, `dropout=0.2`, `gene_likelihood=nb`, `dispersion=gene`, and `latent_distribution=normal`. totalVI used `n_latent=20`, `n_hidden=256`, two encoder layers, one decoder layer, `dropout=0.2`, and `gene_likelihood=nb`.

Downstream classifiers used multinomial logistic regression. The software environment included scikit-learn 1.5.2 and scvi-tools 1.3.3.

### Information-gain and model-specific contrasts

The primary descriptive contrasts were defined as:

\[
G_{\mathrm{info}} = F1_{\mathrm{concat}} - F1_{\mathrm{RNA\,PCA}}
\]

\[
G_{\mathrm{model}} = F1_{\mathrm{totalVI}} - F1_{\mathrm{concat}}
\]

In the development cohort, \(G_{\mathrm{info}} = 0.745 - 0.666 = +0.079\) and \(G_{\mathrm{model}} = 0.779 - 0.745 = +0.034\). The RNA-to-totalVI gap was +0.113, so concat recovered approximately 70% of that gap descriptively. These quantities are descriptive contrasts, not mechanistic attribution claims. The model-specific contrast is descriptive and should not be interpreted as a causal decomposition of architecture effects, because representations, objectives, adaptation procedures, and preprocessing assumptions differ across methods.

### Internal development-cohort evaluation

The development cohort contained 9,494 high-confidence healthy-control cells. RNA PCA, protein PCA, concatenated RNA-protein PCA, MOFA+, scVI, and totalVI were evaluated with the same logistic-regression classifier. Internal macro-F1 values were RNA PCA 0.666, protein PCA 0.512, concat PCA 0.745, MOFA+ 0.549, scVI 0.722, and totalVI 0.779.

The internal totalVI-minus-scVI interval used a cell-level conditional paired bootstrap with n=500 resamples over evaluated cells. This interval is conditional on the observed cell set and does not estimate biological replication.

### Bidirectional transfer, scArches adaptation, and matched concat control

Bidirectional transfer used PBMC10k → PBMC5k as Direction A and PBMC5k → PBMC10k as Direction B, with shared gene/protein ordering retained. Deep-model transfer used unsupervised transductive scArches query adaptation. A source model was trained, unlabeled target features were used for query adaptation, and source and target cells were encoded with the adapted model. A logistic-regression classifier was trained on source post-adaptation latent coordinates and evaluated on target cells.

The inductive concat baseline (10 RNA PCs + 10 protein PCs) was fit on source data only and applied to target cells (0.735 / 0.734). The matched transductive concat control used the same 10+10 dimensionality but fit unsupervised RNA transformation, protein transformation, PCA, and final scaling on source plus unlabeled target features (0.636 / 0.767). Classifier training used source labels only. Target labels were withheld until final scoring and did not enter HVG selection, RNA/protein standardization, PCA, component selection, feature scaling, classifier fitting, or hyperparameter selection. Eligible shared classes for scoring were defined from the evaluation design; they did not supervise feature fitting. This control matches unlabeled-target feature access only and is not algorithmically identical to scArches; differences between joint PCA/scaling and neural query adaptation remain.

Twenty-seed scArches transfer summaries were scVI 0.656±0.018 and totalVI 0.652±0.015 in Direction A, and scVI 0.601±0.020 and totalVI 0.654±0.025 in Direction B. Direction B totalVI–scVI Δ was +0.052 (19/20 seeds). Case C contrasts against transductive concat were approximately +0.016 in Direction A and −0.113 in Direction B.

### Protein degradation experiments

Protein corruption was entry-wise. At each nonzero corruption fraction, originally nonzero measured protein entries were randomly selected and set to zero; preexisting zeros were unchanged. Corruption fractions were 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85.

Training-time degradation corrupted protein entries before totalVI training and retrained totalVI at each condition. The scVI RNA-only reference remained fixed. Post-adaptation degradation used clean source training and clean unsupervised query adaptation, then froze the classifier trained on clean source post-adaptation latents and corrupted target protein values only before re-inference.

Training-time totalVI macro-F1 values were 0.779, 0.771, 0.758, 0.760, 0.750, 0.736, and 0.755 across increasing corruption fractions; the fixed scVI reference was 0.722. Post-adaptation Direction A macro-F1 values were 0.652, 0.650, 0.647, 0.644, 0.641, 0.637, and 0.635. Direction B values were 0.654, 0.652, 0.647, 0.643, 0.637, 0.629, and 0.621.

### Uncertainty and selective prediction

Classifier predictive uncertainty was defined as \(U = 1 - \max_k p_k\), where \(p\) is the logistic-regression class-probability vector. Latent posterior variance was defined as the mean posterior variance over the 20 latent dimensions. Predictive uncertainty was interpreted as classifier-level ambiguity, and latent posterior variance as representation-level dispersion (not formal epistemic uncertainty).

Error-detection AUROC, negative log likelihood, Brier score, expected calibration error, full AURC, \(\mathrm{pAURC}_{[0.5,1.0]}\), and E-AURC were computed where applicable. In PBMC transfer, predictive error AUROC was approximately 0.84–0.86. Direction A latent AUROC was 0.544±0.054 for scVI and 0.669±0.021 for totalVI. Direction B latent AUROC was 0.469±0.075 for scVI and 0.610±0.038 for totalVI.

Selective prediction retained cells with lower uncertainty and deferred higher-uncertainty cells. Risk-coverage summaries were accompanied by cell-type retention at 50% coverage. Retention coefficients of variation were approximately 0.88–1.02 in PBMC transfer and approximately 0.57–0.59 in Lawlor donor validation.

### Training stochasticity and variance decomposition

Training stochasticity used a balanced design with 5 source subsets, 5 model seeds, and 2 run replicates. Variance components were estimated using method-of-moments expected mean squares (MoM EMS). Components were summarized as source, seed, source-by-seed interaction, and residual run-level variation.

For scVI macro-F1, component fractions were source 2.31%, seed 0.00%, interaction 35.45%, and residual 62.23%. For totalVI, fractions were source 0.00%, seed 5.62%, interaction 3.24%, and residual 91.13%. For totalVI minus scVI, fractions were source 5.09%, seed 12.04%, interaction 14.11%, and residual 68.75%.

The residual term denotes run-level residual or stochastic variation under the implemented replication design. It should not be interpreted as irreducible biological variation.

### Lawlor donor-level validation

Lawlor validation used 10 donors and protein-gated author labels in the labeled subset [13]. The primary biological unit was the donor. The primary donor-level metric was present-class macro-F1. Source-only preprocessing and `run_identifier` as batch were retained. Five folds with 8 source / 2 target donors were used for transfer-style folds where applicable. Latent consistency favoring totalVI over scVI was 38/50 among matched fold×seed computational pairs (not biological replicates).

TotalVI improved over scVI by Δ +0.090 on average, with 10/10 donors favoring totalVI and between-donor SD approximately 0.010. Deep-representation averages were 0.800 for scVI and 0.888 for totalVI. Predictive AUROC was 0.801 for scVI and 0.838 for totalVI. Latent AUROC was 0.605 for scVI and 0.650 for totalVI. Full AURC was 0.087 for scVI and 0.036 for totalVI, \(\mathrm{pAURC}_{[0.5,1.0]}\) was 0.157 versus 0.067, and E-AURC was 0.054 versus 0.026. Protein PCA reached approximately 0.927, and concat PCA reached approximately 0.945. Because Lawlor labels were defined using protein markers, the strong protein-PCA and concatenated-PCA performance is consistent with substantial label–modality alignment and should not be interpreted as proof of architecture superiority.

### Software and reproducibility

Analyses used Python 3.11.7, scikit-learn 1.5.2, and scvi-tools 1.3.3. Logistic-regression classifiers used `max_iter=2000`, `solver=lbfgs`, and fixed `random_state` values where repeated seeds were evaluated.

Analysis code available upon request.

## Figure Legends

Figure 1. Central logic: information gain versus model-specific advantage. The schematic separates RNA-only performance, simple multimodal performance, and complex multimodal-model performance. G_info is defined as concat minus RNA PCA, and G_model as totalVI minus concat. These are descriptive contrasts, not mechanistic partitions.

Figure 2. Information-versus-model comparison in the development cohort. Macro-F1 is shown for RNA PCA (0.666), concat PCA (0.745), and totalVI (0.779), with secondary values for protein PCA (0.512), MOFA+ (0.549), and scVI (0.722). The figure highlights G_info=+0.079, G_model=+0.034, and the descriptive recovery of about 70% of the RNA-to-totalVI gap by concat.

Figure 3. Transfer fairness comparison. Panels distinguish source-only inductive concat, target-data-access-matched transductive concat, and scArches-adapted scVI/totalVI for Directions A and B. Case C annotations report totalVI − transductive concat ≈ +0.016 (A) and ≈ −0.113 (B). Transductive concat matches unlabeled target feature access only and is not algorithmically identical to scArches.

Figure 4. Protein degradation experiments. Training-time and post-adaptation target-only designs ask how residual multimodal/model advantage changes as protein signal quality deteriorates. Direction B shows stronger progressive narrowing (0.654 → 0.621) relative to clean scVI (0.601); performance declined gradually rather than catastrophically.

Figure 5. Uncertainty diagnostics, training stochasticity, and selective prediction. Panels summarize predictive and latent error AUROC, risk-coverage behavior, variance-component estimates from the balanced stochasticity design, and cell-type-dependent retention at selective-prediction coverage thresholds.

Figure 6. Lawlor donor-level external validation. The Lawlor analysis uses 10 biological donors and protein-gated author labels. Donor-level present-class macro-F1 favors totalVI over scVI by about +0.090 (10/10 donors), while protein PCA (~0.927) and concat PCA (~0.945) exceed totalVI (~0.888) under protein-gated labels, supporting multimodal information value rather than universal architecture superiority. Annotation: author labels are protein-gated.

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

### Table 2. Information-gain and model-specific descriptive contrasts

| Contrast | Definition | Value | Interpretation |
|---|---|---:|---|
| G_info | F1_concat - F1_RNA_PCA | +0.079 | Simple multimodal information gain |
| G_model | F1_totalVI - F1_concat | +0.034 | Residual model-specific gain in development cohort |
| RNA-to-totalVI gap | F1_totalVI - F1_RNA_PCA | +0.113 | Overall internal multimodal-model gap |
| Descriptive recovery by concat | G_info / RNA-to-totalVI gap | ~70% | Descriptive, not causal |

### Table 3. Matched transfer results and Case C contrasts

| Direction | Inductive concat | Transductive concat | scVI | totalVI | totalVI - scVI | totalVI - transductive concat |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.735 | 0.636 | 0.656±0.018 | 0.652±0.015 | -0.003 | +0.016 |
| B | 0.734 | 0.767 | 0.601±0.020 | 0.654±0.025 | +0.052 | -0.113 |

### Table 4. Protein degradation results

| Corruption fraction | Training-time totalVI | Fixed scVI reference | Post-adapt Direction A totalVI | Direction A Delta | Post-adapt Direction B totalVI | Direction B Delta |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.779 | 0.722 | 0.652 | 0.0000 | 0.654 | 0.0000 |
| 0.10 | 0.771 | 0.722 | 0.650 | -0.0025 | 0.652 | -0.0023 |
| 0.25 | 0.758 | 0.722 | 0.647 | -0.0058 | 0.647 | -0.0069 |
| 0.40 | 0.760 | 0.722 | 0.644 | -0.0087 | 0.643 | -0.0108 |
| 0.55 | 0.750 | 0.722 | 0.641 | -0.0113 | 0.637 | -0.0173 |
| 0.70 | 0.736 | 0.722 | 0.637 | -0.0153 | 0.629 | -0.0253 |
| 0.85 | 0.755 | 0.722 | 0.635 | -0.0178 | 0.621 | -0.0328 |

### Table 5. Uncertainty, selective prediction, and stochasticity diagnostics

| Diagnostic | scVI | totalVI | Notes |
|---|---:|---:|---|
| Direction A latent AUROC | 0.544±0.054 | 0.669±0.021 | Error detection by latent variance |
| Direction B latent AUROC | 0.469±0.075 | 0.610±0.038 | Error detection by latent variance |
| Direction B predictive AUROC under corruption | NA | 0.856 to 0.841 | Clean to highest corruption |
| Direction B full AURC under corruption | NA | 0.060 to 0.079 | Lower is better |
| Direction B pAURC_{[0.5,1.0]} under corruption | NA | 0.112 to 0.144 | Lower is better |
| Direction B E-AURC under corruption | NA | 0.0350 to 0.0448 | Lower is better |
| PBMC 50% retention CV | NA | ~0.88-1.02 | Cell-type-dependent retention |
| Lawlor 50% retention CV | NA | ~0.57-0.59 | Cell-type-dependent retention |
| scVI variance fractions | source 2.31%; seed 0.00%; interaction 35.45%; residual 62.23% | NA | Macro-F1 |
| totalVI variance fractions | NA | source 0.00%; seed 5.62%; interaction 3.24%; residual 91.13% | Macro-F1 |
| Delta variance fractions | source 5.09%; seed 12.04%; interaction 14.11%; residual 68.75% | totalVI - scVI | Paired contrast |

### Table 6. Lawlor donor-level validation

| Metric | scVI | totalVI | Protein PCA | Concat PCA |
|---|---:|---:|---:|---:|
| Deep average macro-F1 | 0.800 | 0.888 | NA | NA |
| Present-class donor mean Delta | NA | +0.090 | NA | NA |
| Donors favoring totalVI over scVI | NA | 10/10 | NA | NA |
| Between-donor SD of Delta | NA | ~0.010 | NA | NA |
| Predictive AUROC | 0.801 | 0.838 | NA | NA |
| Latent AUROC | 0.605 | 0.650 | NA | NA |
| Full AURC | 0.087 | 0.036 | NA | NA |
| pAURC_{[0.5,1.0]} | 0.157 | 0.067 | NA | NA |
| E-AURC | 0.054 | 0.026 | NA | NA |
| Simple baseline macro-F1 | NA | NA | ~0.927 | ~0.945 |

## References

1. Stoeckius M, Hafemeister C, Stephenson W, Houck-Loomis B, Chattopadhyay PK, Swerdlow H, et al. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods. 2017;14:865-868. doi:10.1038/nmeth.4380.
2. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15:1053-1058. doi:10.1038/s41592-018-0229-2.
3. Gayoso A, Steier Z, Lopez R, Regier J, Nazor KL, Streets A, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. Nat Methods. 2021;18:272-282. doi:10.1038/s41592-020-01050-x.
4. Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol. 2020;21:111. doi:10.1186/s13059-020-02015-1.
5. Ashuach T, Gabitto MI, Koodli RV, Saldi GA, Jordan MI, Yosef N. MultiVI: deep generative model for the integration of multimodal data. Nat Methods. 2023;20:1222-1231. doi:10.1038/s41592-023-01909-9.
6. Gong B, Zhou Y, Purdom E. Cobolt: integrative analysis of multimodal single-cell sequencing data. Genome Biol. 2021;22:351. doi:10.1186/s13059-021-02556-z.
7. Jeong Y, Ronen J, Kopp W, Lutsik P, Akalin A. scMaui: a widely applicable deep learning framework for single-cell multiomics integration in the presence of batch effects and missing data. BMC Bioinformatics. 2024;25:257. doi:10.1186/s12859-024-05880-w.
8. Wu M, Goodman N. Multimodal generative models for scalable weakly-supervised learning. Adv Neural Inf Process Syst. 2018;31:5580-5590.
9. Luecken MD, Buttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colome-Tatche M, Theis FJ. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods. 2022;19:41-50. doi:10.1038/s41592-021-01336-8.
10. Lotfollahi M, Naghipourfar M, Luecken MD, Khajavi M, Buttner M, Wagenstetter M, Avsec Z, Gayoso A, Yosef N, Interlandi M, et al. Mapping single-cell data to reference atlases by transfer learning. Nat Biotechnol. 2022;40:121-130. doi:10.1038/s41587-021-01001-7.
11. Guo C, Pleiss G, Sun Y, Weinberger KQ. On calibration of modern neural networks. Proc Mach Learn Res. 2017;70:1321-1330.
12. Geifman Y, El-Yaniv R. Selective classification for deep neural networks. Adv Neural Inf Process Syst. 2017;30:4878-4887.
13. Lawlor N, Nehar-Belaid D, Grassmann JDS, Stoeckius M, Smibert P, Stitzel ML, Pascual V, Banchereau J, Williams A, Ucar D. Single Cell Analysis of Blood Mononuclear Cells Stimulated Through Either LPS or Anti-CD3 and Anti-CD28. Front Immunol. 2021;12:636720. doi:10.3389/fimmu.2021.636720.
14. Hao Y, Hao S, Andersen-Nissen E, Mauck WM III, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184:3573-3587.e29. doi:10.1016/j.cell.2021.04.048.
15. Naeini MP, Cooper GF, Hauskrecht M. Obtaining well calibrated probabilities using Bayesian binning. AAAI. 2015;29:2901-2907. doi:10.1609/aaai.v29i1.9602.
16. Brier GW. Verification of forecasts expressed in terms of probability. Mon Weather Rev. 1950;78:1-3. doi:10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2.
