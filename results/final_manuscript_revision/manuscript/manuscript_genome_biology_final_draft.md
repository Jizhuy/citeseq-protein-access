<!-- Title candidates considered: (1) Reliability of scVI and totalVI for CITE-seq under distribution shift and uncertainty-guided prediction; (2) Stress-testing CITE-seq generative representations under shift and uncertainty; (3) Reliability characterization of variational CITE-seq representations. Selected: candidate 1. -->

# Reliability of scVI and totalVI for CITE-seq under distribution shift and uncertainty-guided prediction

## Title Page

**Title:** Reliability of scVI and totalVI for CITE-seq under distribution shift and uncertainty-guided prediction

**Running title:** Reliability of CITE-seq generative representations

**Authors:** [AUTHOR NAMES AND AFFILIATIONS TO BE COMPLETED]

**Corresponding author:** [CORRESPONDING AUTHOR NAME, ADDRESS, EMAIL]

**Manuscript category:** Research

**Word count:** [TO BE COMPLETED AFTER JOURNAL FORMATTING]

**Figures and tables:** [TO BE COMPLETED]

## Abstract

**Background:** CITE-seq jointly profiles transcriptomes and antibody-derived surface proteins. Probabilistic deep generative models such as scVI and totalVI are widely used for representation learning, yet reliability under protein-entry degradation, dataset or donor shift, training stochasticity, and uncertainty-guided filtering remains incompletely characterized. We evaluated architecture-matched scVI (scVI_matched), totalVI, PCA, and MOFA+ in PBMC CITE-seq cohorts under repeated-seed protocols.

**Results:** In the combined high-confidence PBMC analysis, totalVI improved l2-cell-type F1 over scVI_matched (0.779 vs 0.722; Δ+0.057; target-cell paired bootstrap CI [0.020, 0.091]); MOFA+ reached 0.550. Under Lawlor donor-fold evaluation, scVI achieved F1 0.800 and totalVI 0.888 (Δ+0.088) across 5 folds × 10 seeds, giving 50 paired comparisons within 100 runs. Predictive uncertainty separated correct from incorrect predictions with internal AUROC approximately 0.84-0.86 and Lawlor AUROC 0.801 vs 0.838 for scVI and totalVI. Latent posterior uncertainty was more informative for totalVI than scVI internally (approximately 0.61-0.67 vs 0.47-0.54) and externally (0.650 vs 0.605). Protein-entry corruption, cross-dataset transfer, source-size downsampling, and variance partitioning showed context-dependent totalVI gains, with residual run-level variability largest. Selective prediction risk decreased as coverage was reduced, and totalVI achieved a lower overall AURC, but retained coverage differed across cell types.

**Conclusions:** Across PBMC CITE-seq evaluations, totalVI usually provided stronger representation and uncertainty behavior than scVI_matched, while reliability depended on transfer direction, measured-protein perturbation, and selective-retention imbalance.

## Keywords

CITE-seq; single-cell multi-omics; scVI; totalVI; distribution shift; uncertainty; selective prediction; PBMC

## Background

CITE-seq extends single-cell RNA sequencing by measuring antibody-derived tags alongside transcriptomes in the same cells, enabling joint analyses of gene expression and surface protein state [1]. This paired measurement is especially valuable in immune profiling, where related lymphoid and myeloid populations can share transcriptional programs while differing in canonical surface markers. At the same time, the joint assay creates a reliability challenge: downstream biological conclusions may depend on how computational representations handle sparse counts, unequal modality quality, cohort differences, and missing or uninformative measurements.

Deep generative models are now widely used to represent single-cell data because they combine dimension reduction, batch adjustment, and probabilistic latent-variable modeling. scVI models gene-expression counts with a variational autoencoder and has become a standard reference for transcriptomic representation learning [2]. totalVI extends this framework to paired RNA and protein measurements in CITE-seq by jointly modeling both modalities [3]. These methods are often used upstream of classifiers, label transfer, integration plots, and cell-state comparisons, so their practical value depends not only on mean performance but also on stability and behavior under plausible perturbations.

Alternative multi-omics integration methods provide important reference points. MOFA+ uses factor analysis to decompose multi-modal variation into interpretable latent factors [4]. Later neural frameworks such as MultiVI, Cobolt, and scMaui expand the design space for paired and partially paired assays [5–7], and product-of-experts multimodal variational autoencoders provide related prior art for precision-weighted latent fusion outside the totalVI architecture evaluated here [8]. This manuscript does not attempt to rank all multi-omics methods. Instead, it focuses on a controlled reliability characterization of scVI_matched and totalVI, with PCA and MOFA+ as baselines, in the PBMC CITE-seq cohorts studied here.

Reliability has several distinct meanings in this context. A representation should preserve cell-type information in held-out cells, remain useful when transferred between related cohorts, degrade in understandable ways when measured protein entries are corrupted, and provide uncertainty scores that are informative about errors. Variational models also expose posterior distributions over latent variables, creating the possibility of using latent posterior variance as a signal for ambiguous or shifted observations. Predictive confidence, proper scoring rules, calibration error, and latent posterior width answer related but distinct questions and must be evaluated empirically rather than assumed to be interchangeable [11,15,16].

Single-cell integration benchmarks have emphasized batch mixing, biological conservation, clustering agreement, and neighborhood structure [9]. Those metrics remain essential, but they do not fully characterize reliability under protein degradation, reciprocal dataset transfer, donor-held-out generalization, training stochasticity, or uncertainty-guided abstention. Transfer methods such as scArches support mapping and adaptation across datasets [10], yet reliability still depends on the representation coordinates used for downstream decisions and on whether performance is symmetric across transfer directions.

Prediction uncertainty and selective prediction provide complementary views of reliability. Calibration methods ask whether confidence values match empirical accuracy [11], while selective classification asks whether a model can abstain from low-confidence predictions so that retained predictions have lower error [12]. In biomedical single-cell analysis, abstention is useful only if it is interpreted carefully: lower risk at lower coverage may come with uneven retention of cell types, and those retention patterns can shape which biological conclusions remain well supported.

Here we evaluated scVI_matched and totalVI in PBMC CITE-seq cohorts using frozen annotations, repeated seeds, and explicitly defined perturbations. The study covers internal representation performance, measured-protein corruption, cross-dataset transfer and post-adaptation, source-size downsampling, variance partitioning, predictive and latent uncertainty, selective prediction, and an external donor-fold Lawlor evaluation. The contribution is a systematic reliability characterization for these PBMC CITE-seq settings, rather than a new model architecture.

Throughout, we restrict interpretation to the evaluated PBMC cohorts and the tested pipelines. The results support totalVI as a frequently stronger representation when paired protein measurements are available, but they also identify important boundary conditions: gains are transfer-direction dependent, posterior uncertainty is informative but not a universal detector, and selective prediction can change cell-type composition among retained cells.

## Results

### Study Design And Evaluation Scope

We evaluated scVI_matched and totalVI as the principal generative representations for PBMC CITE-seq analysis, with PCA and MOFA+ used as non-deep and multi-modal baselines. The development analysis combined PBMC10k and PBMC5k after preprocessing and high-confidence annotation. The combined matrix contained 10,849 cells, 15,792 genes, and 14 shared proteins. PBMC10k contributed 6,855 cells, 16,727 genes, and 14 proteins before shared-feature restriction, whereas PBMC5k contributed 3,994 cells, 16,581 genes, and 29 proteins.

Annotations were generated independently of the evaluated latent representations. A Seurat v4 multimodal reference pathway [14] provided fixed labels, and only cells with assignment confidence at or above 0.85 were included in the high-confidence analysis. This retained 9,494 of 10,849 cells. Downstream representation scores therefore evaluated how well model-derived coordinates supported these fixed labels, not whether the evaluated models generated the labels themselves. Those Seurat/Hao-derived labels were used only for development cohorts and were not applied to the external Lawlor analysis.

All classifier-based evaluations used repeated random seeds and a consistent logistic-regression classifier on learned representations. Unless otherwise specified, uncertainty and transfer summaries use seed-level means with standard deviations over seeds. This repeated-seed design separates deterministic conclusions from variation introduced by model initialization, train-test splits, and replicate perturbation masks.

### totalVI Improved High-Confidence PBMC Representation Accuracy

In the high-confidence l2-cell-type task, totalVI outperformed scVI_matched. The F1 score was 0.779 for totalVI and 0.722 for scVI_matched, corresponding to Δ+0.057. A target-cell paired bootstrap with n=500 yielded a confidence interval of [0.020, 0.091], indicating that the observed advantage was consistent at the target-cell level under the specified resampling procedure.

The non-generative and alternative multi-modal baselines provided useful context. MOFA+ reached an F1 score of 0.550, below both scVI_matched and totalVI in this high-confidence PBMC evaluation. PCA was included as a baseline representation throughout the study, but the main verified high-confidence contrast centered on the scVI_matched, totalVI, and MOFA+ values reported above.

These results support the expected advantage of using measured protein information when the downstream task is cell-type recognition in a PBMC CITE-seq setting. They do not imply that totalVI is universally superior in all single-cell multi-omics datasets. Rather, under this annotation, feature intersection, and classifier protocol, totalVI provided a stronger representation of the fixed high-confidence labels than scVI_matched.

### Protein-Entry Corruption Tested Robustness To Missing Measured Signal

We next evaluated robustness by corrupting measured protein entries. For each corruption fraction, randomly selected originally non-zero measured protein entries were set to zero, while preexisting zeros were left unchanged. The tested fractions were 0, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85, with five independent mask seeds at each nonzero condition. totalVI was retrained under each corrupted matrix, whereas scVI_matched remained fixed because it did not use the protein matrix in its representation.

This design targets loss of observed protein signal rather than generic zero inflation. By sampling only entries that were originally non-zero, the perturbation removes measured antibody-derived tag evidence while preserving the original pattern of unmeasured or zero-valued entries. This distinction matters because protein matrices already contain zeros from biological absence, technical dropout, and count sparsity; changing preexisting zeros would not represent removal of measured protein evidence.

Mean Δ macro-F1 (totalVI − scVI_matched) remained positive at every evaluated corruption level (for example, +0.057 at 0%, +0.049 at 10%, and +0.034 at 85%), although some intermediate intervals included zero and the response was not monotonic. Because scVI_matched was fixed across corruption levels, changes in the relative gap reflect the totalVI response to progressively degraded protein input under retraining. This experiment is entry-wise degradation of observed protein values, not a cell-wise missing-modality design.

### Cross-Dataset Transfer Was Direction-Dependent

We assessed cross-dataset transfer between PBMC10k and PBMC5k. In direction A, PBMC10k was used as the source and PBMC5k as the target. Performance was near null between the two evaluated representations: scVI_matched achieved 0.656±0.018 and totalVI achieved 0.652±0.015. The model order is important because the values are nearly identical and do not support a meaningful totalVI advantage in this direction.

In direction B, PBMC5k was used as the source and PBMC10k as the target. Here scVI_matched achieved 0.601±0.020 and totalVI achieved 0.654±0.025, giving Δ≈+0.052. totalVI exceeded scVI_matched in 19 of 20 seeds. Thus, the representation advantage appeared when transferring from the smaller PBMC5k source to the larger PBMC10k target, but not in the reverse direction.

Post-adaptation analyses used the same model coordinates after adaptation. This rule prevents a comparison from mixing pre-adaptation and post-adaptation spaces, and it ensures that downstream classifiers operate on the representation actually produced by the adapted model. Across 20 seeds, seed standard deviations quantify run-to-run variation rather than uncertainty over biological replicates.

The directional asymmetry argues against a simple statement that one model always transfers better. It suggests that cohort size, protein-panel overlap, cell-type composition, and the source-target relationship jointly affect whether totalVI's paired-modality representation improves transfer. For the evaluated PBMC cohorts, totalVI showed a robust directional gain in PBMC5k to PBMC10k transfer but not in PBMC10k to PBMC5k transfer.

### Source-Size Downsampling Probed Training-Set Dependence

To examine whether the transfer pattern could be explained by source size, PBMC10k was downsampled to match the PBMC5k cell count of 3,994 cells. Five downsampled subsets were generated. Models were then evaluated under the same representation and classifier framework used for the main transfer analyses.

The purpose of this experiment was not to create a new biological cohort but to isolate the effect of source-set size within the PBMC10k distribution. If a transfer advantage depended only on the number of source cells, downsampling PBMC10k should move performance toward the PBMC5k-source condition. If it depended on other cohort properties, such as composition or protein-panel structure, downsampling alone would be insufficient to reproduce the full transfer behavior.

The verified results support using source-size analysis as a sensitivity check rather than as a standalone conclusion. Together with the directional transfer results, it shows that model reliability is conditioned by the available training cohort and by how the source relates to the target. This is important for practical deployment, where analysts often train on whichever cohort is available rather than on a source designed to match the target distribution.

### Variance Partitioning Identified Run-Level Residual Variation As Largest

We decomposed performance variability using a 5×5×2×2 replication scheme and the model Y=μ+A_s+B_m+(AB)+ε. The terms represented replicate factors from the design, model effects, their interaction, and residual run-level variation. This analysis was intended to determine whether observed variation was dominated by structured design factors or by run-level residual variation within the replication scheme.

The residual component was the largest. This residual should be interpreted narrowly: it is run-level residual variation under the specified replication design, not irreducible biological variability. It can include effects of random initialization, sampled perturbation masks, classifier stochasticity where applicable, and other unmodeled run-level factors. It should not be reified as an estimate of intrinsic cell-state diversity or unavoidable measurement noise.

The design had limited power for interactions. Consequently, interaction terms were treated as descriptive rather than as evidence for a precise mechanistic interaction between model and replicate factors. The practical implication is that reliability reporting should include replicate seeds and variability summaries, because single-run comparisons may overstate the stability of model differences.

### Predictive Uncertainty Identified Higher-Risk Predictions

Predictive uncertainty was evaluated from classifier outputs on learned representations. The uncertainty score, denoted U_pred, was used to distinguish correct from incorrect predictions. Internally, AUROC values were approximately 0.84-0.86, indicating that prediction-level uncertainty carried useful information about error risk in the development PBMC analyses.

In the Lawlor external donor-fold evaluation, U_pred AUROC was 0.801 for scVI and 0.838 for totalVI. These values show that uncertainty remained informative beyond the development cohorts, and that totalVI provided a stronger error-ranking signal in that external setting. The comparison is descriptive across the fixed donor-fold and seed protocol rather than a claim about all CITE-seq cohorts.

Calibration was also assessed, including expected calibration error and Brier score summaries where relevant. In the Lawlor evaluation, ECE was tied between the two principal models, while selective-prediction behavior differed more clearly. This distinction is important: a model may rank uncertain predictions well without producing globally better calibrated probabilities, and calibration summaries alone may not capture the utility of uncertainty for abstention.

### Latent Posterior Uncertainty Was More Informative For totalVI

Latent posterior uncertainty was defined as the mean posterior variance over d=20 latent dimensions. Using scvi-tools 1.3.3, posterior distributions were obtained with return_dist, extracting qz.loc and qz.scale squared. The scalar latent uncertainty U_latent was then evaluated for its ability to separate correct from incorrect downstream predictions.

In internal transfer direction A, latent AUROC was 0.544±0.054 for scVI and 0.669±0.021 for totalVI. In direction B, the corresponding values were 0.469±0.075 for scVI and 0.610±0.038 for totalVI. Thus, totalVI latent uncertainty was consistently more informative than scVI latent uncertainty in the two internal transfer directions, whereas scVI latent uncertainty was near chance or below chance depending on direction.

The Lawlor external analysis showed the same ordering: latent AUROC was 0.605 for scVI and 0.650 for totalVI, with totalVI higher in 19 of 20 paired fold-seed comparisons. This result supports latent posterior variance as a useful reliability signal for totalVI in the evaluated PBMC CITE-seq settings. It does not imply that latent variance alone is sufficient for error detection, nor that transcriptome-only latent uncertainty is uninformative in every context.

### Lawlor Donor-Fold Evaluation Confirmed External Reliability Patterns

The external Lawlor analysis used public CITE-seq data with author labels only. The dataset contained 16,175 singlets, of which 5,207 were labeled for the evaluated classification task, with 12,776 genes and 12 proteins across 10 donors. Evaluation used five donor folds and ten seeds per fold, producing 100 model runs and 50 paired scVI-totalVI fold-by-seed comparisons. The 50 paired comparisons describe consistency across the fold×seed grid and are not 50 independent biological replicates.

In this setting, scVI achieved F1 0.800 and totalVI achieved 0.888, giving Δ+0.088. The totalVI advantage was consistent across the paired fold×seed comparisons. Hierarchical variability summaries separated between-fold variability from within-fold seed variability, reflecting that donor assignment and seed-level stochasticity contribute different forms of uncertainty.

Uncertainty results in Lawlor aligned with the development analyses. Predictive uncertainty AUROC was 0.801 for scVI and 0.838 for totalVI; latent uncertainty AUROC was 0.605 for scVI and 0.650 for totalVI. For selective prediction, AURC was 0.087 for scVI and 0.036 for totalVI, while ECE was tied. These results indicate that totalVI improved both representation accuracy and uncertainty-guided ranking in the Lawlor donor-fold analysis, while calibration did not separate the models as clearly.

### Selective Prediction Reduced Risk At Lower Coverage With Uneven Retention

Selective prediction evaluated the risk retained after abstaining on the most uncertain predictions. Risk was defined as 1 minus accuracy among retained cells at a given coverage level. Coverage was evaluated from 100% down to 50%. Prediction risk decreased as coverage was reduced, and totalVI achieved a lower overall AURC.

The coverage-risk curves show that uncertainty can support triage: retained predictions become more accurate when lower-confidence cells are rejected. However, the retained set was not compositionally uniform. At 50% coverage, retention coefficient of variation in PBMC analyses was approximately 0.88-1.02, whereas Lawlor showed approximately 0.57-0.59. CD4/CD8 naive and CD4 memory populations were rejected more frequently, while monocyte and B-cell populations were retained more often.

These patterns matter for interpretation. Selective prediction can improve the reliability of retained labels, but it also changes which cell types contribute to downstream summaries. A low-risk retained subset should therefore be reported together with coverage and retention balance by cell type, especially when downstream analyses compare abundance, marker expression, or condition-specific changes across immune populations.

## Discussion

This study provides a reliability-focused evaluation of scVI_matched and totalVI for PBMC CITE-seq representation learning. Across the verified analyses, totalVI often improved classifier performance and uncertainty behavior relative to scVI_matched, especially in high-confidence PBMC labels, PBMC5k to PBMC10k transfer, and the external Lawlor donor-fold evaluation. The results are strongest when stated with their scope: PBMC cohorts, the specified preprocessing and annotation pipeline, and the tested PCA and MOFA+ baselines.

The high-confidence PBMC result is consistent with the value of measured protein information for immune cell-type discrimination. totalVI reached l2 F1 0.779 compared with 0.722 for scVI_matched, with a target-cell paired bootstrap confidence interval for the difference of [0.020, 0.091]. This is a meaningful representation gain under the fixed label and classifier protocol. The lower MOFA+ value of 0.550 further indicates that merely using a multi-modal method is not sufficient; the modeling assumptions and downstream representation geometry matter.

At the same time, the cross-dataset transfer results caution against universalizing the totalVI advantage. PBMC10k to PBMC5k transfer was near null, with scVI_matched 0.656±0.018 and totalVI 0.652±0.015. PBMC5k to PBMC10k transfer favored totalVI, with scVI_matched 0.601±0.020 and totalVI 0.654±0.025. This directional asymmetry is practically important because many applications involve mapping from one available reference to another target cohort, and the reference-target relationship may dominate the observed model difference.

The corruption experiment adds another boundary condition. By setting randomly selected originally non-zero measured protein entries to zero and leaving preexisting zeros unchanged, the analysis directly tested sensitivity to loss of measured protein signal. Since totalVI was retrained under each corrupted matrix and scVI_matched was fixed, the experiment primarily characterizes how a protein-aware representation responds to degraded protein evidence. It should not be interpreted as a claim that one model and another saw identical input information.

Uncertainty analyses were a central part of the reliability characterization. Predictive uncertainty was useful for distinguishing correct from incorrect predictions internally and externally, and totalVI had higher Lawlor U_pred AUROC than scVI. Latent posterior uncertainty was more model-dependent. totalVI U_latent was informative internally and externally, whereas scVI U_latent was weaker and sometimes near chance in internal transfer. These findings support empirical evaluation of uncertainty outputs rather than assuming that any posterior variance estimate is automatically useful for triage.

Selective prediction translated uncertainty into an operational decision rule. As coverage decreased from 100% toward 50%, risk decreased, and totalVI had a lower overall AURC. This is encouraging for workflows where analysts may prefer fewer but more reliable labels. Yet the retained cells were imbalanced by cell type, with naive and memory CD4-related groups rejected more often and monocyte and B-cell groups retained more often. The method is therefore best viewed as uncertainty-guided triage with composition reporting, not as a neutral filtering step.

The variance analysis reinforces the need for replicated evaluation. The largest component was the residual run-level term under the 5×5×2×2 design. Because this residual reflects the replication scheme rather than irreducible biology, it should motivate practical reporting of seeds, folds, and perturbation masks instead of biological overinterpretation. Limited interaction power also means that interaction patterns should be treated as hypothesis-generating unless a larger design is run.

The external Lawlor evaluation strengthens the manuscript because it used donor folds and author labels rather than the development annotation pathway. totalVI improved F1 from 0.800 to 0.888 and also improved predictive and latent uncertainty AUROC. The paired 5 fold × 10 seed grid provided 50 within-grid model comparisons across 100 runs, but those comparisons remain tied to 10 donors and should not be described as 50 independent biological replicates. Separating between-fold from within-fold seed variability is essential for interpreting this external evidence.

Several limitations remain. The study is restricted to PBMC CITE-seq cohorts and does not test tissue atlases, perturbation screens, or disease settings with different protein panels and class imbalance. The missing-modality analysis is a controlled protein-entry corruption experiment rather than a complete model of real missingness mechanisms. The classifier-based evaluation measures usefulness of representations for fixed labels; it does not prove that every biological axis is better represented. These limits do not undermine the main conclusion, but they define where further validation is needed.

## Conclusions

In evaluated PBMC CITE-seq cohorts, totalVI generally provided more reliable representations and uncertainty signals than scVI_matched, while the strength of this advantage depended on data context. totalVI improved high-confidence PBMC cell-type F1, showed a directional transfer gain from PBMC5k to PBMC10k, and performed strongly in the external Lawlor donor-fold evaluation. Predictive uncertainty and latent posterior variance were useful for error ranking, especially for totalVI, and selective prediction reduced retained risk as coverage decreased. However, protein-entry corruption, transfer asymmetry, run-level residual variability, and cell-type retention imbalance show that reliability is conditional. These results support reporting repeated seeds, transfer direction, uncertainty metrics, and retained-cell composition when applying generative CITE-seq representations.

## Methods

### Overview

We evaluated reliability of CITE-seq representations in development PBMC cohorts and in an external Lawlor donor-fold cohort. The principal comparison was scVI_matched versus totalVI. PCA and MOFA+ were included as baselines where applicable. All scientific results reported here were generated before manuscript drafting and were treated as fixed. No additional experiments were performed for this manuscript file.

Reliability was assessed through eight analysis families: high-confidence representation classification, measured-protein entry corruption, cross-dataset transfer, post-adaptation evaluation, source-size downsampling, variance partitioning, predictive and latent uncertainty, and selective prediction. The same broad evaluation logic was used throughout: train or extract a representation, fit a standardized classifier where required, evaluate against fixed labels, and summarize repeated seeds or folds according to the analysis design.

### Development Datasets

The development data comprised two PBMC CITE-seq cohorts referred to as PBMC10k and PBMC5k. PBMC10k contained 6,855 cells, 16,727 genes, and 14 proteins after preprocessing for the evaluated cohort. PBMC5k contained 3,994 cells, 16,581 genes, and 29 proteins. The combined analysis used 10,849 cells, 15,792 genes, and 14 shared proteins after harmonizing features across the two cohorts.

For the high-confidence labeled analysis, 9,494 of 10,849 cells passed the label-confidence threshold of 0.85. These high-confidence cells were used for the main l2-cell-type representation comparison. The shared-protein restriction ensured that totalVI and multi-modal baselines operated on a consistent protein panel in the combined development analysis.

### Preprocessing

RNA and protein matrices were processed consistently within each evaluation branch before model fitting. Gene features were harmonized across cohorts for combined analyses, and protein features were restricted to the shared panel where required. Cells and features failing upstream quality filters were excluded before the reported cell, gene, and protein counts were finalized.

For the external Lawlor analysis, gene identifiers were mapped using Ensembl GRCh37.87, and duplicate gene mappings were summed. Protein features were harmonized to 12 antibody aliases. These steps produced the verified Lawlor analysis matrix of 16,175 singlets, 5,207 labeled cells, 12,776 genes, and 12 proteins.

### Annotation And High-Confidence Labeling

Development-set annotation used a Seurat v4 SCANVI pathway independent of the evaluated latent representations. The evaluated scVI_matched, totalVI, PCA, and MOFA+ coordinates were not used to create the labels against which they were scored. A label-confidence threshold of 0.85 defined the high-confidence subset, retaining 9,494 of 10,849 combined cells.

The annotation strategy was designed to avoid circularity between representation learning and evaluation. Because the labels were fixed before downstream model comparison, classifier performance on evaluated representations reflects alignment with the annotation pathway rather than self-consistency of the evaluated model. Seurat v4 reference-mapping concepts and CITE-seq reference resources motivated the annotation workflow [14].

### PCA Baseline

PCA was used as a baseline low-dimensional representation. The PCA representation provided a non-generative comparator for classifier-based evaluations. It was included to contextualize whether deep variational models improved performance beyond a standard linear dimension-reduction approach.

PCA preprocessing followed the same feature harmonization used for the relevant analysis branch. Classifiers were then trained on the resulting principal-component coordinates under the same classifier protocol used for other representations.

### MOFA+ Baseline

MOFA+ was evaluated as a multi-modal factor-analysis baseline using mofapy2 version 0.7.5. MOFA+ provides latent factors from multiple data views and served as an established multi-omics comparator to totalVI. In the high-confidence PBMC l2-cell-type evaluation, MOFA+ reached F1 0.550.

MOFA+ inputs used the harmonized RNA and protein features appropriate for the combined PBMC analysis. Downstream classification used the same logistic-regression evaluation framework applied to the other latent representations.

### scVI_matched Model

scVI_matched denotes the transcriptome-only scVI representation trained under an architecture-matched evaluation configuration relative to totalVI. Exact configuration: n_latent=20, n_hidden=256, n_layers=2, dropout_rate=0.2, gene_likelihood="nb", dispersion="gene", and latent_distribution="normal". Gene-expression counts were provided as raw UMI counts for scVI. Training used batch_size=256, max_epochs=200, train_size=0.9, and early stopping on elbo_validation with patience=45. Software: scvi-tools 1.3.3.

The scVI_matched representation did not use protein counts. In the protein-entry corruption analysis, scVI_matched was therefore held fixed across protein corruption fractions. scVI and totalVI are not literally identical architectures even under matched latent width and hidden size; the matched comparison controls major capacity settings rather than claiming architectural identity.

### totalVI Model

totalVI was trained on paired RNA and protein matrices with n_latent=20, n_hidden=256, encoder layers=2, decoder layers=1, dropout=0.2, and gene_likelihood=NB. Proteins were supplied via obsm["protein_expression"] (14 proteins in development analyses; 12 in Lawlor). Training used batch_size=256, max_epochs=200, train_size=0.9, early stopping on elbo_validation with patience=45, and learning rate 4e-3 with reduce_lr_on_plateau enabled as recorded in the totalVI runner. Software: scvi-tools 1.3.3.

For protein-entry corruption experiments, totalVI was retrained at each corruption fraction and mask seed. For transfer and Lawlor analyses, representations were extracted under the same downstream classifier and uncertainty protocols as scVI_matched.

### Classifier Protocol

Representation quality was evaluated with scikit-learn LogisticRegression using max_iter=2000, solver=lbfgs, and random_state equal to the evaluation seed. Unspecified keyword arguments, including C, penalty, and multi_class, used scikit-learn defaults for the installed version. Classifiers were fit on the representation coordinates and evaluated against fixed labels.

The same classifier protocol was used across scVI_matched, totalVI, PCA, and MOFA+ where classifier-based metrics were reported. This design isolates differences in representation coordinates as much as possible while keeping the downstream supervised model constant.

### Representation Metrics

The primary representation metric was macro or task-specific F1 for fixed cell-type labels, reported according to each analysis. In the high-confidence combined PBMC l2-cell-type analysis, the verified values were totalVI 0.779, scVI_matched 0.722, and MOFA+ 0.550. The totalVI-scVI_matched difference was Δ+0.057 with a target-cell paired bootstrap confidence interval of [0.020, 0.091] using n=500 bootstrap resamples.

For transfer and external analyses, performance was summarized across seeds or fold-seed combinations. Standard deviations over seeds are reported where specified. No unverified p-values were computed or added.

### Protein-Entry Corruption

Protein-entry corruption was applied by randomly selecting originally non-zero measured protein entries and setting them to zero. Preexisting zeros were unchanged. Corruption fractions were 0, 0.10, 0.25, 0.40, 0.55, 0.70, and 0.85. Five mask seeds were used for each corruption condition.

The design intentionally removed measured protein evidence while preserving the original zero pattern. totalVI was retrained for each corrupted matrix and mask seed. scVI_matched was fixed because its input did not include proteins. The resulting comparisons characterize robustness of the totalVI representation to progressive measured-protein signal loss relative to a stable transcriptome-only baseline.

### Cross-Dataset Transfer And Post-Adaptation

Cross-dataset transfer was evaluated in both directions between PBMC10k and PBMC5k using scArches query adaptation [10]. Direction A used PBMC10k as source and PBMC5k as target. Direction B used PBMC5k as source and PBMC10k as target. Each direction used 20 model seeds, and reported standard deviations are over seeds. Query adaptation used prepare_query_anndata and load_query_data with max_epochs=200 and plan weight_decay=0.0.

Direction A yielded near-null performance differences: scVI_matched 0.656±0.018 and totalVI 0.652±0.015. Direction B favored totalVI: scVI_matched 0.601±0.020 and totalVI 0.654±0.025, with Δ≈+0.052 and totalVI higher in 19 of 20 seeds.

**Post-adaptation coordinate rule:** source and target embeddings were always encoded with the adapted query model. Downstream classifiers and uncertainty metrics were computed in that shared post-adaptation coordinate system, avoiding comparisons that mix pre- and post-adaptation spaces.

### Predictive Uncertainty

Predictive uncertainty was defined as U_pred(i) = 1 − max_k p_ik, where p_ik denotes the multinomial logistic-regression probability for class k on cell i. Incorrect prediction indicators were scored against fixed labels. Error discrimination used AUROC and AUPRC; probability quality used negative log-likelihood and multiclass Brier score [16]; calibration used expected calibration error (ECE) with 15 bins, equal-frequency as primary and equal-width as sensitivity [11,15].

Internal U_pred AUROC values were approximately 0.84–0.86. In the Lawlor donor-fold evaluation, U_pred AUROC was 0.801 for scVI_matched and 0.838 for totalVI. These AUROC values quantify ranking of error risk by uncertainty, not calibration of the numeric probabilities themselves.

### Selective Prediction

Selective prediction ranked cells by ascending predictive uncertainty (equivalently, descending confidence) and retained the most confident predictions at coverages {1.0, 0.9, 0.8, 0.7, 0.6, 0.5}. Risk at coverage c was defined as risk(c) = 1 − accuracy among retained cells. AURC was calculated as the integral of the risk–coverage curve across the evaluated coverage range and is therefore an integrated summary, not a risk value at a single coverage.

Prediction risk decreased as coverage was reduced, and totalVI achieved a lower overall AURC. In Lawlor, AURC was 0.087 for scVI_matched and 0.036 for totalVI. Per-class retention was retained cells of class k divided by original evaluable cells of class k. Retention imbalance was summarized by range and coefficient of variation across classes. At 50% coverage, retention CV was approximately 0.88–1.02 in PBMC analyses and approximately 0.57–0.59 in Lawlor. CD4/CD8 naive and CD4 memory populations were rejected more often, while monocyte and B-cell populations were retained more often.

### Latent Posterior Uncertainty

Latent posterior uncertainty was defined as U_latent(i) = (1/d) Σ_j Var[z_ij | x_i] with d=20. Using scvi-tools version 1.3.3, get_latent_representation(..., return_dist=True) returned qz.loc and qz.scale; posterior variances were taken as qz.scale squared. Monte Carlo confirmation of variance semantics is summarized in Supplementary Methods.

In direction A, U_latent AUROC was 0.544±0.054 for scVI_matched and 0.669±0.021 for totalVI. In direction B, values were 0.469±0.075 for scVI_matched and 0.610±0.038 for totalVI (SD over 20 seeds). In Lawlor, values were 0.605 for scVI_matched and 0.650 for totalVI, with totalVI higher in 19 of 20 fold-seed comparisons.

### Source-Size Downsampling

To assess source-size dependence, PBMC10k was downsampled to 3,994 cells, matching the PBMC5k cell count. Five downsampled subsets were generated. Each subset was evaluated with the same representation and classifier logic used in the cross-dataset transfer analyses. The analysis was a sensitivity check rather than a causal identification of source size as the sole driver of transfer asymmetry.

### Twenty-Seed Robustness And Seed-Pair Bootstrap

Development transfer and uncertainty grids used 20 model seeds with paired scVI_matched/totalVI comparisons at matched seeds. Seed-level paired bootstrap used 10,000 replicates over the 20 seed pairs where reported. Target-cell bootstrap (n=500) was used for the high-confidence representation ΔF1 interval and was not pooled with seed-pair bootstrap.

### Variance Decomposition And Replicate Randomness

Variance partitioning used the balanced design of 5 source subsets × 5 model seeds × 2 run replicates × 2 models under

Y_smr = μ + A_s + B_m + (AB)_sm + ε_smr,

where s indexes source subset, m indexes model seed, and r indexes run replicate. Method-of-moments estimators from expected mean squares were used for components source_subset, model_seed, interaction, and residual, with both raw and non-negative truncated estimates reported. Negative raw estimates can arise under MoM sampling variability and were truncated for nonnegative summaries.

The residual component was largest for most major analyses and is interpreted as run-level residual/stochastic variability under the implemented replication scheme. It is not automatically irreducible biological noise, optimization noise alone, or unexplained biology. Power to exclude modest interactions was limited.

A determinism audit showed that identical full seed state produced identical reruns; therefore the replicated design separated model_seed from run_replicate_seed controlling shuffling/dropout/adaptation stochasticity as implemented.

### Lawlor External Cohort

The external evaluation used the Lawlor public CITE-seq cohort [13], accessed through HCA dataset efea6426-510a-4b60-9a19-277e52bfa815 and ENA project PRJEB40376 / ERP124005. The analysis used author labels only. The processed dataset contained 16,175 singlets, 5,207 labeled cells, 12,776 genes, 12 proteins, and 10 donors.

Gene processing used Ensembl GRCh37.87, with duplicate gene identifiers summed. Protein processing harmonized 12 aliases. Donor folds were used for evaluation, with five folds and ten seeds per fold. The batch covariate was run_identifier rather than donor, preventing donor identity from being used as the batch variable in model fitting.

The Lawlor evaluation comprised 100 runs and 50 paired fold×seed model comparisons. The 50 paired comparisons are descriptive paired consistency checks within the donor-fold and seed grid, not 50 independent biological replicates. Hierarchical variability summaries distinguished between-fold variability from within-fold seed variability.

### Statistics

Bootstrap uncertainty for the high-confidence l2-cell-type difference used a target-cell paired bootstrap with n=500 resamples. Seed-level summaries report means and standard deviations over seeds where specified. Lawlor summaries used donor folds and seed replication as described above.

No unverified p-values were added. Directional consistency counts, such as 19 of 20 seeds or fold-seed comparisons, were reported only where verified. Confidence intervals, standard deviations, and AUROC values are interpreted within the corresponding replication design.

### Missing-Modality Limitation

The protein-entry corruption analysis is a controlled perturbation of measured protein entries, not a full model of missing modality acquisition. It leaves preexisting zeros unchanged and removes a specified fraction of originally non-zero measured entries. Real missingness can arise from antibody-panel design, sample processing, background correction, and cohort-specific technical factors, so the corruption results should be interpreted as a stress test rather than a complete missing-modality benchmark.

## Declarations

### Ethics Approval And Consent To Participate

This study used secondary analysis of publicly available human single-cell data. Confirmation of the applicable ethics status, consent language, and whether additional institutional review was required is pending: [HUMAN/SUPERVISOR CONFIRMATION REQUIRED]. No new human participants were recruited for this manuscript.

### Consent For Publication

Not applicable to newly collected participant data. Final confirmation required: [HUMAN/SUPERVISOR CONFIRMATION REQUIRED].

### Availability Of Data And Materials

The development datasets and processed analysis objects will be made available through [ZENODO OR OTHER REPOSITORY PLACEHOLDER]. The external Lawlor data are available through HCA dataset efea6426-510a-4b60-9a19-277e52bfa815 and ENA project PRJEB40376 / ERP124005. Public repository information for scripts and reproducible analysis workflows will be provided at [PUBLIC REPOSITORY PLACEHOLDER].

### Availability Of Code

Code used for preprocessing, model training, evaluation, and figure generation will be deposited at [PUBLIC REPOSITORY PLACEHOLDER] with release archive [ZENODO PLACEHOLDER].

### Competing Interests

The authors declare the following competing interests: [COMPETING INTERESTS PLACEHOLDER]. If none, replace with: The authors declare that they have no competing interests.

### Funding

Funding information: [FUNDING PLACEHOLDER].

### Authors' Contributions

Author contributions: [AUTHOR CONTRIBUTION PLACEHOLDER].

### Acknowledgements

Acknowledgements: [ACKNOWLEDGEMENTS PLACEHOLDER].

### Use Of Generative AI Tools

Generative AI tools were used to assist with manuscript drafting and editing: [TOOL NAMES AND VERSIONS TO BE CONFIRMED BY HUMAN AUTHORS]. All scientific claims, numerical results, citations, and final text require human author review and confirmation: [HUMAN CONFIRMATION REQUIRED].

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

15. Naeini MP, Cooper GF, Hauskrecht M. Obtaining well calibrated probabilities using Bayesian binning. Proc AAAI Conf Artif Intell. 2015;29:2901–2907.

16. Brier GW. Verification of forecasts expressed in terms of probability. Mon Weather Rev. 1950;78:1–3. doi:10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2.

## Figure Legends

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

## Tables

Table 1 and Table 2 are provided as separate files. In the manuscript, Table 1 should summarize datasets, feature counts, labels, and evaluation splits. Table 2 should summarize model comparisons, transfer directions, uncertainty metrics, and selective-prediction summaries.

## Supplementary Note

This manuscript reports only verified analyses in the evaluated PBMC CITE-seq settings. Interpretation should remain tied to the stated cohorts, labels, seeds, folds, corruption design, and model configurations. The comparison is between scVI_matched, totalVI, and the specified PCA and MOFA+ baselines, not a comprehensive benchmark of all single-cell multi-omics integration methods.

The selective-prediction analyses should be read together with retention summaries. Lower retained risk at lower coverage can be useful for high-confidence annotation workflows, but retained cells are not compositionally neutral. Reports using selective prediction should include both the risk-coverage curve and cell-type retention balance.

The Lawlor fold×seed grid provides paired consistency evidence across the implemented evaluation design. It should not be interpreted as 50 independent biological replicates, because biological donor structure remains the primary external replication unit. Between-fold and within-fold seed variability therefore provide complementary, not interchangeable, information.
