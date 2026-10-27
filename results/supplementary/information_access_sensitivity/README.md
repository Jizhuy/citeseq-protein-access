# SECONDARY INFORMATION-ACCESS SENSITIVITY

This folder contains the **historical feature-transductive** internal analysis. Here scVI_matched and totalVI were fitted on the features (not the labels) of all development cells, including the cells that were later held out for classification. It is reported only as a secondary sensitivity analysis of information access. It is **not** the primary analysis; the primary analysis is the strict train-only evaluation in `results/primary/`.

| File | Content | Generating script |
| --- | --- | --- |
| `feature_transductive_bootstrap/core_contrast_bootstrap_summary.csv` | Paired-bootstrap contrasts under feature-transductive access | `scripts/run_feature_transductive_bootstrap.py` |
| `feature_transductive_bootstrap/core_contrast_bootstrap_audit.md` | Bootstrap audit | same |
| `class_stratified_bootstrap/*` | Global vs class-stratified bootstrap of the feature-transductive contrasts | `scripts/cache_feature_transductive_predictions.py`, `scripts/run_class_stratified_bootstrap.py` |
| `information_access_sensitivity_comparison.csv` | Side-by-side strict train-only vs feature-transductive estimates | Frozen derived artifact; no standalone generator retained (hand-formatted from the frozen tables of both analyses; see `docs/RESULT_MANIFEST.md` §9) |

The feature-transductive embeddings come from `scripts/run_feature_transductive_models.py`.
