# PRIMARY ANALYSIS — strict train-only internal evaluation

Status: **COMPLETE — PRIMARY ANALYSIS**

Every representation (RNA PCA, protein PCA, concatenated RNA–protein PCA, scVI_matched, totalVI) was fitted on 7,573 high-confidence training cells only. Held-out high-confidence cells (n = 1,921) were embedded afterwards without parameter updates. One downstream classifier was used for all representations. scVI_matched and totalVI use seed 0; the bootstrap uses B = 10,000 paired cell-level resamples (seed 20260923). The class-stratified bootstrap uses seed 20260924.

| File | Content |
| --- | --- |
| `strict_train_only/tables/strict_train_only_point_estimates.csv` | Macro-F1 per representation; G_access, G_assoc, Gap_total, R_access |
| `strict_train_only/tables/strict_train_only_bootstrap_summary.csv` | Paired-bootstrap percentile intervals and prop>0 |
| `strict_train_only/tables/strict_train_only_bootstrap_global_vs_stratified.csv` | Global vs class-stratified bootstrap |
| `strict_train_only/tables/strict_train_only_bootstrap_class_stratified.csv` | Class-stratified bootstrap summary |
| `strict_train_only/tables/strict_train_only_bootstrap_draws.csv` | Bootstrap draws |
| `strict_train_only/tables/strict_train_only_heldout_predictions.npz` | Held-out predictions used by the bootstrap |
| `strict_train_only/tables/run_meta.json` | Cell counts, software versions, device, seed |
| `STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md` | Audit of the strict train-only protocol |

Generating script: `scripts/run_strict_train_only_internal.py`.

Intervals are conditional on the fitted representations and the held-out cells. They are not biological confidence intervals.
