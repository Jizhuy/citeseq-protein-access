<!--
Alternative title: Separating protein-access gains from model-associated gains in CITE-seq representation learning
-->

# Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift

**Jizhu Yang**  
Research Manuscript | Unpublished manuscript | September 2026

## Abstract

**Background:** CITE-seq representation learning can improve cell-type prediction, but gains may reflect access to paired surface proteins rather than a more complex multimodal representation. We quantified protein-access gain and the residual totalVI-associated contrast using simple RNA–protein controls and stress-tested that residual under distribution shift. **Results:** In a 9,494-cell development cohort, RNA PCA, concat PCA, and totalVI reached macro-F1 0.666, 0.745, and 0.779. Protein-access gain was G_access=+0.079 [cell-level conditional paired-bootstrap 95% interval 0.043–0.132]. The residual totalVI-associated contrast was G_assoc=+0.034 [−0.012–0.073]; this interval crossed zero. The RNA-to-totalVI gap was +0.113 [0.084–0.154], and simple concatenation recovered approximately 70% of that gap descriptively [0.40–1.11]. Under transfer, inductive concat remained strong (0.735/0.734); feature-access-matched transductive concat reached 0.636 and 0.767, with totalVI about +0.016 and −0.113 relative to that feature-access control. Protein degradation narrowed the residual advantage, and training stochasticity limited precision for small contrasts. Across 10 Lawlor donors, totalVI exceeded scVI_matched (mean paired Δ=+0.090; 10/10; donor means 0.794 vs 0.884), while protein PCA and concat were stronger under protein-gated labels. **Conclusions:** Protein-access gain was strong and reproducible; residual totalVI-associated advantage was smaller and context-dependent.
## Background

CITE-seq jointly measures transcriptomes and antibody-derived surface proteins in the same single cells [1]. This paired design is valuable for immune profiling because cell identities and activation states often depend on both transcriptional programs and surface-marker phenotypes. It also creates an interpretation problem: an improved representation may benefit from access to protein measurements, from the modeling procedure, or from an interaction between the two.

Deep generative methods are now widely used for single-cell representation learning. scVI provides an RNA-only variational representation [2], while totalVI extends this framework to paired RNA and protein measurements [3]. Other multimodal approaches, including MOFA+, MultiVI, Cobolt, scMaui, weakly supervised multimodal generative modeling, atlas integration, scArches, and WNN-style analysis, have expanded the range of available strategies for paired or partially paired single-cell data [4-10,14].

A central practical question is whether a complex multimodal representation adds value beyond a simple representation that has access to the same broad feature types. Comparing totalVI only with an RNA-only model can answer whether a particular multimodal workflow performs better than a particular RNA-only workflow, but it does not show how much of the improvement is already obtained by adding paired protein features to a simple baseline. A concatenated RNA-protein PCA control therefore provides an important reference point.

This distinction matters under distribution shift. scArches query adaptation uses unlabeled target features during adaptation [10], whereas a source-only PCA baseline can be applied inductively without fitting on target features. A fairer transfer stress test should include a feature-access-matched transductive concat control. This control matches unlabeled target features but not the learning objective or adaptation algorithm of scArches, so it is a feature-access control rather than an algorithmic replica.

Reliability also depends on protein quality, training stochasticity, and uncertainty behavior. Protein degradation can reveal whether residual totalVI-associated advantage depends on intact antibody-derived signal. Repeated training can show whether small contrasts are large enough to interpret confidently. Classifier predictive uncertainty and latent posterior variance can detect different reliability losses, while selective prediction can reduce retained-cell error at the cost of uneven cell-type retention.

External donor validation adds a biological unit of replication that computational resampling cannot replace. The Lawlor blood CITE-seq dataset provides 10 donors with paired RNA-protein measurements and author labels [13]. Because those labels are protein-gated, strong protein-aligned baselines are expected and should be interpreted as evidence of multimodal measurement value under protein-defined annotation rather than universal superiority of any representation class.

This study makes three contributions. First, we quantify the empirical performance improvement associated with access to paired protein measurements using strong simple multimodal controls and distinguish it from the residual gain associated with a more complex multimodal representation. Second, we test whether this residual totalVI-associated contrast persists under source–target shift, feature-access-matched transductive comparison, protein degradation, and training stochasticity. Third, we evaluate reliability using uncertainty-aware selective prediction and donor-level external validation.

## Results

### Simple multimodal protein-access gain explained much of the internal improvement

The internal development analysis used 9,494 high-confidence healthy-control cells with matched RNA and protein measurements. **Figure 1.** summarizes the study logic: RNA-only performance was compared with a simple RNA-protein representation and then with totalVI, separating the protein-access contrast from the residual totalVI-associated contrast.

RNA PCA reached macro-F1 0.666, concatenated RNA-protein PCA reached 0.745, and totalVI reached 0.779. We defined the protein-access gain as \(G_{\mathrm{access}}=F1_{\mathrm{concat}}-F1_{\mathrm{RNA\,PCA}}\), giving +0.079. We defined the residual totalVI-associated contrast as \(G_{\mathrm{assoc}}=F1_{\mathrm{totalVI}}-F1_{\mathrm{concat}}\), giving +0.034. This contrast is descriptive and does not isolate a causal architecture effect. The full RNA-to-totalVI gap was +0.113.

**Figure 2.** shows the internal performance hierarchy and cell-level conditional paired-bootstrap intervals. \(G_{\mathrm{access}}\) was +0.079 [0.043, 0.132], prop>0=1.0. \(G_{\mathrm{assoc}}\) was +0.034 [−0.012, 0.073], prop>0=0.9045, and the 95% interval crossed zero. The total gap was +0.113 [0.084, 0.154]. The descriptive recovery ratio was \(R_{\mathrm{access}}\approx0.70\) [0.405, 1.114]; the denominator was stable, with minimum gap 0.056 and never ≤0. These are cell-level conditional intervals over evaluated cells, not biological confidence intervals.

Secondary internal representations supported the same interpretation. Protein PCA alone reached 0.512, MOFA+ reached 0.549, and scVI reached 0.722. TotalVI was the highest internal method among those tested, but the larger and more certain internal contrast was the gain associated with protein access in a simple multimodal representation.

### Residual totalVI-associated advantage was not stable across transfer settings

Bidirectional transfer tested whether residual totalVI-associated advantage persisted under source-target shift. **Figure 3.** compares source-only inductive concat, the feature-access-matched transductive concat control, scVI, and totalVI across Directions A and B.

The source-only inductive concat baseline remained strong in both transfer directions, reaching 0.735 in Direction A and 0.734 in Direction B. The feature-access-matched transductive concat control reached 0.636 in Direction A and 0.767 in Direction B. This feature-access control used source plus unlabeled target features for unsupervised representation fitting, while source labels alone were used for classifier training.

scVI reached 0.656 in Direction A and 0.601 in Direction B. totalVI reached 0.652 in Direction A and 0.654 in Direction B. Thus, totalVI was near scVI in Direction A, while Direction B showed a totalVI-scVI difference of +0.052 with 19/20 seeds favoring totalVI.

The comparison against the feature-access control changed the interpretation. totalVI minus transductive concat was approximately +0.016 in Direction A but −0.113 in Direction B. Matching unlabeled target feature access therefore did not yield a stable totalVI advantage. The control matches unlabeled target features but not the learning objective or adaptation algorithm of scArches, so this result should be read as a feature-access-matched stress test, not as proof that the two procedures are otherwise equivalent.

### Protein degradation progressively narrowed residual totalVI-associated advantage

Protein degradation asked whether the residual totalVI-associated contrast persisted when antibody-derived signal was corrupted. **Figure 4.** summarizes training-time degradation with retraining and post-adaptation target-only degradation.

Training-time degradation corrupted nonzero protein entries before totalVI training and retrained totalVI at each nonzero corruption level using five mask/retraining replicates with model_seed fixed at 0. The clean point was a single uncorrupted run. Mean macro-F1 values were 0.779 at 0.00 (n=1), 0.771±0.010 at 0.10, 0.758±0.012 at 0.25, 0.760±0.012 at 0.40, 0.750±0.016 at 0.55, 0.736±0.020 at 0.70, and 0.755±0.011 at 0.85. The fixed clean scVI reference was 0.722.

The high-corruption rebound should be interpreted cautiously. Non-monotonic high-corruption means were within mask/retraining variability and are not evidence of improved performance under more severe corruption. The main pattern was a gradual narrowing of the residual totalVI-associated advantage as protein signal quality degraded.

Post-adaptation target-only degradation showed the same direction in the transfer setting where totalVI most clearly exceeded scVI. Direction B totalVI declined from 0.654 to 0.621 at the highest corruption level, while clean scVI was 0.601. Thus, the totalVI advantage over the clean RNA-only reference narrowed but did not disappear in that setting.

### Uncertainty diagnostics identified distinct forms of reliability loss

Uncertainty diagnostics were used to characterize reliability loss rather than to define a new primary performance endpoint. **Figure 5.** shows classifier predictive uncertainty, latent posterior variance, risk-coverage behavior, and training stochasticity summaries.

Classifier predictive uncertainty was \(U=1-\max_k p_k\), where \(p_k\) is the logistic-regression probability for class \(k\). This score ranked classification errors well in PBMC transfer, with predictive error AUROC values around 0.84-0.86. Under Direction B target-only protein corruption, predictive AUROC declined mildly from 0.856 to 0.841, indicating that confidence still ranked errors but became less informative as protein inputs degraded.

Latent posterior variance measured representation-level dispersion, not formal epistemic uncertainty. Direction A latent error AUROC was 0.544±0.054 for scVI and 0.669±0.021 for totalVI. Direction B latent error AUROC was 0.469±0.075 for scVI and 0.610±0.038 for totalVI. Under Direction B corruption, latent AUROC increased from 0.610 to 0.627 while predictive AUROC and risk-coverage summaries worsened. This divergence shows that the two uncertainty signals capture different behavior.

Risk-coverage metrics showed practical reliability loss under degradation. In Direction B, full AURC increased from 0.060 to 0.079, \(\mathrm{pAURC}_{[0.5,1.0]}\) increased from 0.112 to 0.144, and E-AURC increased from 0.0350 to 0.0448. Negative log likelihood and Brier score also worsened, whereas expected calibration error was less sensitive.

### Training stochasticity limited the precision of small residual totalVI-associated contrasts

Training stochasticity was evaluated because small residual contrasts can be difficult to distinguish from run-level variation. A balanced design varied source subsets, model seeds, and run replicates, and variance components were estimated using method-of-moments expected mean squares.

For scVI macro-F1, the variance fractions were source 2.31%, seed 0.00%, source-by-seed interaction 35.45%, and residual run-level variation 62.23%. For totalVI, the corresponding fractions were source 0.00%, seed 5.62%, interaction 3.24%, and residual 91.13%. For the paired totalVI minus scVI contrast, the fractions were source 5.09%, seed 12.04%, interaction 14.11%, and residual 68.75%.

These summaries explain why small residual totalVI-associated contrasts should be interpreted with repeated-run context. Direction A transfer was effectively near null, and the internal \(G_{\mathrm{assoc}}\) bootstrap interval crossed zero. Direction B totalVI versus scVI was more consistent because the mean difference was larger and 19/20 seeds favored totalVI, but the feature-access control still showed that residual totalVI advantage was not stable across transfer directions.

### Selective prediction reduced error but introduced cell-type-dependent retention

Selective prediction retained cells with lower uncertainty and deferred cells with higher uncertainty. This reduced error among retained cells, especially when classifier predictive uncertainty was used, but it also changed the retained cell-type composition.

At 50% coverage, retention coefficients of variation were approximately 0.88-1.02 in PBMC transfer and approximately 0.57-0.59 in Lawlor donor validation. These values indicate that low-uncertainty retained cells were not a compositionally neutral subset of the evaluated cells.

The practical implication is that risk-coverage curves should be reported together with retention by cell type. Selective prediction can support review workflows, but deferring uncertain cells may preferentially remove rare, transitional, or marker-ambiguous populations. Retained-set accuracy alone would overstate reliability if the retained population shifts away from difficult biological classes.

### External donor validation confirmed multimodal information value under protein-gated labels

The Lawlor validation used 10 biological donors and author labels in the labeled subset [13]. **Figure 6.** presents donor-level validation together with simple protein-aligned baselines.

All 10 donors showed a positive totalVI-minus-scVI_matched difference in present-class macro-F1 (exact two-sided sign test, *P*=0.00195; mean Δ≈+0.090, median≈0.089, SD≈0.010). Donor-mean deep-representation scores were 0.794 for scVI_matched and 0.884 for totalVI. Latent consistency favoring totalVI over scVI_matched was 38/50 across matched fold-by-seed computational comparisons; these are computational comparisons, not biological replicates.

Simple protein-aligned baselines were stronger in the same setting. Fold-mean protein PCA reached approximately 0.927, and concatenated RNA-protein PCA reached approximately 0.945, exceeding totalVI. Because the Lawlor labels were protein-gated, this result is consistent with strong label-modality alignment. The donor analysis therefore supports multimodal measurement value across donors under protein-gated labels, while reinforcing the need to compare complex models with simple multimodal baselines.

## Discussion

Access to paired protein measurements accounted for a substantial part of the observed multimodal performance improvement. The remaining advantage associated with totalVI relative to a simple RNA–protein representation was smaller and strongly dependent on evaluation context.

The internal development cohort made this pattern explicit. \(G_{\mathrm{access}}\) was +0.079 with a bootstrap interval entirely above zero, whereas \(G_{\mathrm{assoc}}\) was +0.034 and its 95% interval crossed zero. The descriptive recovery ratio \(R_{\mathrm{access}}\approx0.70\) indicates that a simple RNA-protein representation recovered much of the RNA-to-totalVI gap, but it should not be interpreted as a mechanistic partition. The residual totalVI-associated contrast is descriptive and should not be interpreted as a causal decomposition of architecture effects.

Transfer analysis showed why this caution is necessary. Source-only inductive concat remained strong in both directions. totalVI exceeded scVI consistently in Direction B but not Direction A. More importantly, totalVI did not retain a stable advantage over the feature-access-matched transductive concat control. Because this control matches unlabeled target features but not the learning objective or adaptation algorithm of scArches, it does not prove algorithmic equivalence; it does show that feature access alone can change the apparent residual advantage.

Protein degradation further constrained the interpretation. Training-time corruption gradually narrowed totalVI performance relative to the clean setting, and Direction B post-adaptation degradation reduced totalVI from 0.654 to 0.621. The non-monotonic high-corruption training-time mean was within mask/retraining variability and should not be described as evidence that more severe corruption improves performance.

Uncertainty diagnostics identified different forms of reliability loss. Classifier predictive uncertainty ranked errors and supported selective prediction, while latent posterior variance reflected representation-level dispersion. Because these measures can move differently under degradation, neither should be treated as a complete reliability measure. Selective prediction reduced retained-cell error but introduced cell-type-dependent retention, making coverage composition a necessary companion to risk summaries.

Training stochasticity limited the precision of small residual totalVI-associated contrasts. Residual run-level variation was large in the balanced replication design, especially for totalVI and paired totalVI-scVI contrasts. This does not invalidate larger and direction-consistent effects, but it argues against overinterpreting small differences without matched baselines, repeated runs, and interval estimates.

Lawlor donor validation strengthened the conclusion that paired proteins can be highly valuable across biological donors. The totalVI-scVI donor contrast favored totalVI in all 10 donors, but protein PCA and concat exceeded totalVI under protein-gated labels. This is not a contradiction; it shows that when labels depend heavily on surface markers, direct protein-aligned representations can be extremely strong.

Several limitations remain. The internal and transfer analyses are PBMC-focused and classification-centered; conclusions should not be generalized to trajectory inference, perturbation response, or rare-state discovery without testing. Antibody panels are limited (14 shared internal proteins; 12 matched Lawlor proteins). Only one primary complex multimodal model (totalVI) was stress-tested against concat, with MOFA+ as a secondary factor baseline; additional multimodal architectures were not newly trained in this revision. Lawlor labels are protein-gated, and only 10 biological donors are available. The feature-access-matched concat control does not reproduce scArches’ learning objective. Entry-wise zero-masking approximates protein degradation; an exploratory PCA-only marker-channel dropout sensitivity is reported in Supplementary Results and does not test totalVI under panel failure. Bootstrap intervals are conditional cell-level intervals, not biological confidence intervals. Uncertainty scores do not guarantee biological label correctness, and selective prediction changes retained cell-type composition. These limits make the conclusion narrower but more reliable: protein access is a strong contributor, while residual totalVI-associated gain must be demonstrated in each evaluation context.

## Conclusion

Protein-access gain was strong and reproducible in the evaluated CITE-seq settings. The residual totalVI-associated contrast was smaller, less precisely estimated internally, and context-dependent under transfer, degradation, and training stochasticity. Multimodal representation claims should therefore be evaluated against simple RNA-protein baselines that make feature access explicit.

## Methods

### Datasets and biological evaluation units

**Internal PBMC development and transfer.** Analyses used the public scvi-tools / 10x Genomics CITE-seq example files `pbmc_10k_protein_v3.h5ad` and `pbmc_5k_protein_v3.h5ad`, downloaded from the scverse exampledata URLs with a YosefLab/scVI-data GitHub fallback (project `config/config.yaml`). No GEO/SRA accession was present in local metadata, so none is claimed. Datasets were loaded separately and combined by **inner join** on shared genes and proteins (project loader; the scvi-tools helper was not used as the primary API because it silently intersects features and can fill missing proteins with zero). The combined inner object contained **10,849** cells, **15,792** shared genes, and **14** shared proteins (PBMC10k **6,855** cells; PBMC5k **3,994** cells). Batch labels were dataset identity (`batch` ∈ {PBMC10k, PBMC5k}).

**Quality control.** Project configuration sets `additional_cell_filter: false`. No additional minimum/maximum gene, count, mitochondrial-fraction, or doublet filters were applied beyond those already present in the downloaded processed objects (documented as inherited from totalVI-style reproducibility preprocessing of the source files). Gene and protein panels for modeling used the shared inner-join universe; highly variable gene selection for **PCA baselines** is described below and is distinct from VAE training, which used the full shared gene count matrix.

**Development evaluation set.** Primary internal contrasts used **9,494** high-confidence cells (`annotation_tier_l2 == "high"`, confidence threshold 0.85). Macro-F1 on L2 labels was the primary metric.

**Lawlor external cohort.** External validation used the Lawlor blood CITE-seq dataset [13] (HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005`). The biological unit was the **donor** (*n* = 10). Author labels in the labeled subset are protein-gated. Batch covariate for deep models was sequencing lane (`run_identifier`); donor was **not** supplied as a batch key. Primary donor metric: present-class macro-F1 (classes absent in a donor are excluded for that donor).

### Annotation provenance and label-independence audit

Internal development labels were generated through an **RNA-only** Seurat v4 reference annotation pathway. An RNA-only re-annotation audit agreed with the stored high-confidence labels on all **9,494** evaluated cells (agreement = 1.0). This rules out the simple explanation that the internal protein-access advantage was produced by defining labels from the same protein features used as predictors. It does **not** independently establish biological correctness of those labels.

Lawlor author labels are protein-gated; strong protein PCA / concat performance is therefore expected under label–modality alignment and is not evidence of universal representation superiority.

### Feature preprocessing and PCA representations

Unless noted, unsupervised transforms for internal evaluation were fit on **training** high-confidence cells only and applied to held-out cells. For transfer, inductive transforms were fit on **source** cells only.

**RNA PCA.** Select **2,000** highly variable genes (`scanpy` `flavor="seurat_v3"`); normalize to **1×10⁴** counts per cell; `log1p`; z-score with clipping at **±10**; PCA with **10** components (`random_state=0`). (Project `config.yaml` lists `n_top_genes: 4000` for historical VAE-related settings; the **primary PCA path hard-codes 2000** and does not read that config value.)

**Protein PCA.** Centered log-ratio (CLR) transform on raw protein counts with pseudocount **1.0**; z-score using training/source mean and sample SD (`ddof=1`); PCA with **10** components.

**Concatenated RNA–protein PCA.** Concatenate 10 RNA PCs and 10 protein PCs; apply `StandardScaler` fit on training/source concatenated coordinates only (inductive). In the feature-access-matched transductive control, unsupervised RNA/protein transforms, PCA, and final scaling were fit on **source + unlabeled target** features jointly.

**MOFA+.** Frozen multimodal factor-analysis baseline on paired RNA–protein inputs (secondary comparator).

### Deep representations (scVI_matched and totalVI)

Reported internal and transfer “scVI” results use the capacity-matched setting **`scVI_matched`** (not the project’s PHASE-3 scVI default of `n_hidden=128`, `n_layers=1`, `gene_likelihood=zinb`, `dropout_rate=0.1`).

**scVI_matched:** `SCVI` (scvi-tools **1.3.3**); `n_latent=20`, `n_hidden=256`, `n_layers=2`, `dropout_rate=0.2`, `gene_likelihood="nb"`, `dispersion="gene"`, `latent_distribution="normal"`, `use_observed_lib_size=True`; `setup_anndata` with `layer="counts"`, `batch_key="batch"`. Training: `max_epochs=200`, `early_stopping=True` (monitor `elbo_validation`, patience **45**, mode `min`), `train_size=0.9`, `batch_size=256`. Learning rate was not overridden in code (scvi-tools TrainingPlan default).

**totalVI:** `TOTALVI`; `n_latent=20`, `n_hidden=256`, `n_layers_encoder=2`, `n_layers_decoder=1`, encoder/decoder dropout **0.2**, `gene_likelihood="nb"`, `gene_dispersion="gene"`, `protein_dispersion="protein"`, `latent_distribution="normal"`; protein key `protein_expression` (14 proteins). Training additionally used `lr=4e-3` and `reduce_lr_on_plateau=True`, with the same epoch/early-stopping/batch settings as above.

**scArches transfer.** For each direction and seed: train source model → `prepare_query_anndata` / `load_query_data` → unsupervised query adaptation (`query_max_epochs=200`, `weight_decay=0.0`) → encode **source and target** with the **adapted** model → fit logistic regression on post-adaptation **source** coordinates → score target. Target labels were not used in fitting. This is unsupervised **transductive** adaptation, not fully inductive prediction.

### Downstream classification

All primary representation comparisons used `sklearn.linear_model.LogisticRegression` with `max_iter=2000`, `solver="lbfgs"`, `random_state` equal to the evaluation seed (development seed **0**). `C`, `class_weight`, and `multi_class` were left at scikit-learn **1.5.2** defaults (`C=1.0`, `class_weight=None`, multinomial behavior under `lbfgs`). The same classifier strategy was used across representations so that comparisons focus on representation quality.

### Internal split and leakage controls

A stratified train/test split (`test_fraction=0.20`, stratify by `batch`, seed **0**) was assigned on the combined 10,849-cell object before high-confidence restriction (implied all-cell sizes ≈8,679 / 2,170). Development metrics used the high-confidence subset of that split (**7,573** train / **1,921** test). Chronology: split → train-only unsupervised feature fitting (PCA path) → fixed deep embeddings from trained models → train-only classifier fit → held-out evaluation. Transductive concat is the exception that **intentionally** uses unlabeled target features for unsupervised fitting only; target labels remain withheld until scoring.

### Protein-access and residual totalVI-associated contrasts

\[
G_{\mathrm{access}} = F1_{\mathrm{concat}} - F1_{\mathrm{RNA\,PCA}},\quad
G_{\mathrm{assoc}} = F1_{\mathrm{totalVI}} - F1_{\mathrm{concat}},\quad
Gap_{\mathrm{total}} = F1_{\mathrm{totalVI}} - F1_{\mathrm{RNA\,PCA}},\quad
R_{\mathrm{access}} = G_{\mathrm{access}} / Gap_{\mathrm{total}}.
\]

\(G_{\mathrm{assoc}}\) is a **descriptive residual totalVI-associated contrast**, not a causal architecture effect: totalVI and concatenated PCA differ in objective, optimization, and representation class. Because numerator and denominator of \(R_{\mathrm{access}}\) are recomputed jointly across bootstrap resamples, \(R_{\mathrm{access}}\) is not constrained to [0,1] and is not a causal fraction, variance-explained proportion, or mechanistic decomposition.

Point estimates: RNA PCA 0.666, concat 0.745, totalVI 0.779; \(G_{\mathrm{access}}=+0.079\), \(G_{\mathrm{assoc}}=+0.034\), \(Gap_{\mathrm{total}}=+0.113\), \(R_{\mathrm{access}}\approx0.70\).

### Conditional paired bootstrap

Primary intervals used a **cell-level conditional paired bootstrap** (*B* = **10,000**; seed **20260923**) over held-out cells with fixed prediction vectors for RNA PCA, concat, and totalVI. Each replicate resampled cells with replacement (global, not class-stratified), recomputed macro-F1 and contrasts, and formed **percentile** 2.5%/97.5% intervals. These are conditional on the evaluated cells and prediction vectors—not donor-level, cohort-generalization, or biological confidence intervals. A class-stratified paired-bootstrap sensitivity (within-class resampling) is reported in Supplementary Results; it does not replace the primary analysis.

### Bidirectional transfer and feature-access-matched control

Direction A: PBMC10k→PBMC5k; Direction B: PBMC5k→PBMC10k.

1. **Source-only inductive concat:** fit HVG/RNA transform, RNA PCA, protein CLR/z/PCA, concat scaling, and classifier on source; apply unchanged to target.  
2. **Feature-access-matched transductive concat:** fit unsupervised steps on source + unlabeled target features; classifier on source labels only. Matches **feature access**, not scArches’ objective or neural adaptation.  
3. **scVI_matched / totalVI + scArches** as above (20 seeds for primary summaries).

### Protein degradation

**Training-time:** originally nonzero protein entries randomly set to zero (preexisting zeros unchanged) at fractions 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85; totalVI retrained; five mask/retraining replicates at nonzero levels with `model_seed=0`; clean level *n*=1. Means are reported without forcing monotonicity; high-corruption rebound is interpreted within mask/retraining variability.

**Post-adaptation:** clean source training and clean query adaptation → classifier on clean post-adaptation source coordinates → freeze classifier → corrupt **target** proteins only → re-infer target → evaluate. Cell-type–specific Direction B results are summarized in Supplementary Results.

**Marker-channel dropout (sensitivity):** entire protein channels zeroed before protein PCA / concat on the development cohort (totalVI not retrained in this sensitivity); see Supplement.

### Uncertainty, calibration, and selective prediction

Predictive uncertainty: \(U_{\mathrm{pred},i}=1-\max_k p(y_i=k\mid z_i)\) from the downstream logistic regression (not totalVI posterior uncertainty). Latent posterior dispersion: mean posterior variance across latent dimensions \(d\) (representation-level dispersion; **not** formal epistemic uncertainty). Diagnostics included error-detection AUROC/AUPRC, NLL, multiclass Brier score, ECE (equal-frequency and equal-width; equal-frequency treated as primary where both reported), full AURC, \(\mathrm{pAURC}_{[0.5,1.0]}\), and E-AURC. Selective prediction ordered cells by increasing uncertainty and retained the lowest-uncertainty fraction; primary risk used accuracy-based risk among retained cells. Macro-F1 and balanced accuracy versus coverage are reported as Supplementary sensitivities because the primary selective metric is accuracy-based.

### Training stochasticity and variance decomposition

Balanced computational design: **5** source subsets × **5** model seeds × **2** run replicates (Direction A focus in PHASE 12A), with paired Δ = totalVI − scVI_matched. Variance components for source, seed, source×seed interaction, and residual run-level variation were estimated by method-of-moments expected mean squares. Residual denotes computational/run-level variation under this design, not biological variance.

### Lawlor donor-level external validation

Five deterministic folds held out **two** donors each so that every donor was a target exactly once; **10** computational seeds per deep-model fold yielded 50 matched scVI_matched–totalVI fold×seed comparisons (computational consistency, e.g. latent AUROC 38/50). Donor-level present-class macro-F1 used seed-averaged within-donor scores. Simple baselines (RNA/protein/concat PCA) used source-donor-only fits applied to held-out target donors. Donor summary statistics (mean/median/SD/IQR/min/max of paired Δ; exact sign test) are reported in Supplementary Results. Fold-mean protein PCA ≈ **0.927** and concat ≈ **0.945**; donor-mean deep scores were **0.794** (scVI_matched) and **0.884** (totalVI), with mean paired Δ ≈ **+0.090** (SD ≈ 0.010; 10/10 donors).

### Software and reproducibility

Analyses used Python **3.11.7**, scvi-tools **1.3.3**, scikit-learn **1.5.2**, scanpy **1.10.3**, anndata **0.11.4**, PyTorch **2.14.0** (CUDA builds used for PHASE 10 server runs). Environment files: `environment.yml`, `requirements.txt`. Code: https://github.com/Jizhuy/robust-multiomics-integration. An immutable archival release and DOI will accompany the submitted version.

## Data Availability

The Lawlor Human Cell Atlas dataset is available under HCA project efea6426-510a-4b60-9a19-277e52bfa815 and ENA accessions PRJEB40376 / ERP124005 [13].

PBMC10k/PBMC5k analyses used the public scvi-tools CITE-seq example AnnData files `pbmc_10k_protein_v3.h5ad` and `pbmc_5k_protein_v3.h5ad` (scverse exampledata URLs; YosefLab/scVI-data fallback). No GEO/SRA accession is claimed.

## Code Availability

Analysis code and reproducible workflows are publicly available at https://github.com/Jizhuy/robust-multiomics-integration.git. An immutable archival release and DOI will accompany the submitted version of the manuscript.

## Acknowledgements

I thank Bo Li for valuable discussions, methodological suggestions, and constructive feedback on the framing and interpretation of this work.

## Figure Legends

**Figure 1.** Study logic for separating protein-access and residual totalVI-associated contrasts. RNA PCA defines the RNA-only reference, concatenated RNA-protein PCA defines the simple multimodal feature-access reference, and totalVI defines the more complex multimodal representation. The quantities \(G_{\mathrm{access}}\), \(G_{\mathrm{assoc}}\), and \(R_{\mathrm{access}}\) are descriptive contrasts.

**Figure 2.** Internal development-cohort contrasts with bootstrap intervals. RNA PCA reached 0.666, concat reached 0.745, and totalVI reached 0.779. Cell-level conditional paired-bootstrap 95% intervals with n=10000 were \(G_{\mathrm{access}}=+0.079\) [0.043, 0.132], \(G_{\mathrm{assoc}}=+0.034\) [−0.012, 0.073], and total gap +0.113 [0.084, 0.154]. The \(G_{\mathrm{assoc}}\) interval crossed zero. These intervals are conditional cell-resampling intervals, not biological confidence intervals.

**Figure 3.** Transfer comparison across source-only inductive concat, feature-access-matched transductive concat control, scVI, and totalVI. Inductive concat reached 0.735/0.734, feature-access transductive concat reached 0.636/0.767, scVI reached 0.656/0.601, and totalVI reached 0.652/0.654 for Directions A/B. totalVI minus transductive concat was approximately +0.016 in Direction A and −0.113 in Direction B. The feature-access control matches unlabeled target features but not the learning objective or adaptation algorithm of scArches.

**Figure 4.** Protein degradation experiments. Training-time degradation shows mean±SD across five mask/retraining replicates at nonzero corruption levels, with model_seed fixed at 0 and a single clean run at 0.00. Values were 0.779, 0.771±0.010, 0.758±0.012, 0.760±0.012, 0.750±0.016, 0.736±0.020, and 0.755±0.011 across increasing corruption. High-corruption non-monotonicity is interpreted as within mask/retraining variability. Post-adaptation Direction B totalVI declined from 0.654 to 0.621, with clean scVI at 0.601.

**Figure 5.** Uncertainty diagnostics, training stochasticity, and selective prediction. Classifier predictive uncertainty ranked errors with AUROC around 0.84-0.86 internally, while latent posterior variance measured representation-level dispersion. Risk-coverage metrics worsened under Direction B protein corruption, and selective prediction reduced retained-cell error while producing cell-type-dependent retention. Stochasticity summaries show substantial residual run-level variation for small contrasts.

**Figure 6.** Lawlor donor-level external validation. The Lawlor analysis used 10 biological donors and protein-gated author labels. All 10 donors showed a positive totalVI-minus-scVI_matched difference (exact two-sided sign test, *P*=0.00195; mean Δ≈+0.090). Protein PCA (~0.927) and concat (~0.945) exceeded totalVI (0.884) under protein-gated labels. The result supports multimodal measurement value under protein-gated labels and emphasizes the importance of simple protein-aligned baselines.

## Tables

### Table 1. Internal development-cohort methods

| Representation | Modalities | Macro-F1 | Role |
|---|---:|---:|---|
| RNA PCA | RNA | 0.666 | RNA-only reference |
| Protein PCA | Protein | 0.512 | Protein-only secondary baseline |
| Concatenated RNA-protein PCA | RNA + protein | 0.745 | Simple multimodal protein-access reference |
| MOFA+ | RNA + protein | 0.549 | Multimodal factor-analysis secondary baseline |
| scVI | RNA | 0.722 | RNA-only deep representation |
| totalVI | RNA + protein | 0.779 | Multimodal deep representation |

### Table 2. Core internal contrasts with bootstrap intervals

| Contrast | Definition | Estimate | 95% interval | prop>0 | Interpretation |
|---|---|---:|---:|---:|---|
| \(G_{\mathrm{access}}\) | \(F1_{\mathrm{concat}} - F1_{\mathrm{RNA\,PCA}}\) | +0.079 | [0.043, 0.132] | 1.0 | Protein-access gain |
| \(G_{\mathrm{assoc}}\) | \(F1_{\mathrm{totalVI}} - F1_{\mathrm{concat}}\) | +0.034 | [−0.012, 0.073] | 0.9045 | Residual totalVI-associated contrast; interval crosses zero |
| \(Gap_{\mathrm{total}}\) | \(F1_{\mathrm{totalVI}} - F1_{\mathrm{RNA\,PCA}}\) | +0.113 | [0.084, 0.154] | NA | Total RNA-to-totalVI gap |
| \(R_{\mathrm{access}}\) | \(G_{\mathrm{access}} / Gap_{\mathrm{total}}\) | ~0.70 | [0.405, 1.114] | NA | Descriptive recovery ratio |

Footnote: Intervals are cell-level conditional paired-bootstrap 95% intervals with n=10000 resamples over evaluated cells. They are not biological confidence intervals. The \(R_{\mathrm{access}}\) denominator was stable, with minimum gap 0.056 and never ≤0.

### Table 3. Transfer contrasts

| Direction | Inductive concat | Feature-access transductive concat | scVI | totalVI | totalVI - scVI | totalVI - transductive concat |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.735 | 0.636 | 0.656 | 0.652 | -0.004 | +0.016 |
| B | 0.734 | 0.767 | 0.601 | 0.654 | +0.052 | -0.113 |

### Table 4. Protein degradation training-time results

| Protein corruption fraction | totalVI macro-F1 |
|---:|---:|
| 0.00 | 0.779 (n=1) |
| 0.10 | 0.771±0.010 |
| 0.25 | 0.758±0.012 |
| 0.40 | 0.760±0.012 |
| 0.55 | 0.750±0.016 |
| 0.70 | 0.736±0.020 |
| 0.85 | 0.755±0.011 |

### Table 5. Uncertainty, stochasticity, and selective prediction

| Diagnostic | Value |
|---|---:|
| Predictive error AUROC in PBMC transfer | ~0.84-0.86 |
| Direction B predictive AUROC under corruption | 0.856 to 0.841 |
| Direction A latent AUROC, scVI / totalVI | 0.544±0.054 / 0.669±0.021 |
| Direction B latent AUROC, scVI / totalVI | 0.469±0.075 / 0.610±0.038 |
| Direction B full AURC under corruption | 0.060 to 0.079 |
| Direction B pAURC\(_{[0.5,1.0]}\) under corruption | 0.112 to 0.144 |
| Direction B E-AURC under corruption | 0.0350 to 0.0448 |
| PBMC 50% retention CV | ~0.88-1.02 |
| Lawlor 50% retention CV | ~0.57-0.59 |
| scVI variance fractions | source 2.31%; seed 0.00%; interaction 35.45%; residual 62.23% |
| totalVI variance fractions | source 0.00%; seed 5.62%; interaction 3.24%; residual 91.13% |
| totalVI-scVI variance fractions | source 5.09%; seed 12.04%; interaction 14.11%; residual 68.75% |

### Table 6. Lawlor donor-level validation

| Metric | Value |
|---|---:|
| Biological donors | 10 |
| totalVI - scVI donor-level difference | ~+0.090 |
| Donors favoring totalVI over scVI_matched | 10/10 |
| Exact two-sided donor-level sign test | *P*=0.00195 |
| scVI_matched donor-mean present-class macro-F1 | 0.794 |
| totalVI donor-mean present-class macro-F1 | 0.884 |
| Protein PCA macro-F1 | ~0.927 |
| Concat PCA macro-F1 | ~0.945 |
| Latent consistency favoring totalVI | 38/50 computational comparisons |

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

15. Naeini MP, Cooper GF, Hauskrecht M. Obtaining well calibrated probabilities using Bayesian binning. Proc AAAI Conf Artif Intell. 2015;29:2901-2907.

16. Brier GW. Verification of forecasts expressed in terms of probability. Mon Weather Rev. 1950;78:1-3. doi:10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2.
