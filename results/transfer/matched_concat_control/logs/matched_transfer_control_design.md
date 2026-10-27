# Matched transfer control design

- timestamp: 2026-09-22T14:34:49Z
- purpose: match ACCESS TO UNLABELED TARGET FEATURES for simple concat PCA
- NOT a fully matched architecture comparison with scArches

## Designs

1. **Source-fitted inductive concat PCA** (preserved): HVG/PCA/scaling fit on SOURCE only; TARGET transformed; classifier on SOURCE labels.
2. **Target-data-access-matched transductive concat PCA** (new): HVG/PCA/scaling fit on SOURCE + unlabeled TARGET features jointly; classifier still SOURCE labels only; TARGET labels used only at final scoring.
3. **scVI / totalVI**: frozen PHASE11B scArches unsupervised transductive adaptation (20 seeds).

## Integrity

- Target labels never enter HVG selection, scaling, PCA, concat StandardScaler, or LogisticRegression training.
- Primary endpoint: macro-F1; secondary: accuracy, balanced accuracy.

## Results snapshot

### direction_A
- RNA_PCA_10 | none_inductive_source_fit_only: macro-F1=0.6994
- concat_PCA_10_10 | none_inductive_source_fit_only: macro-F1=0.7350
- concat_PCA_10_10 | transductive_joint_unsupervised_feature_fit: macro-F1=0.6363
- scVI_matched | scArches_unsupervised_transductive: macro-F1=0.6556 ±0.0175
- totalVI | scArches_unsupervised_transductive: macro-F1=0.6524 ±0.0146

### direction_B
- RNA_PCA_10 | none_inductive_source_fit_only: macro-F1=0.6265
- concat_PCA_10_10 | none_inductive_source_fit_only: macro-F1=0.7339
- concat_PCA_10_10 | transductive_joint_unsupervised_feature_fit: macro-F1=0.7670
- scVI_matched | scArches_unsupervised_transductive: macro-F1=0.6015 ±0.0195
- totalVI | scArches_unsupervised_transductive: macro-F1=0.6540 ±0.0254

- direction_A: totalVI − transductive_concat = +0.0161
- direction_B: totalVI − transductive_concat = -0.1131
