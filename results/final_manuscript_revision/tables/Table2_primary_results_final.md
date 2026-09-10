# Table 2. Primary reliability findings

| Setting | Model | Metric | Mean | SD | n / unit |
|---|---|---|---|---|---|
| Transfer Dir A | scVI_matched | macro-F1 | 0.656 | 0.018 | 20 seeds |
| Transfer Dir A | totalVI | macro-F1 | 0.652 | 0.015 | 20 seeds |
| Transfer Dir B | scVI_matched | macro-F1 | 0.601 | 0.020 | 20 seeds |
| Transfer Dir B | totalVI | macro-F1 | 0.654 | 0.025 | 20 seeds |
| Lawlor | scVI_matched | macro-F1 | 0.800 | 0.015 | 50 fold×seed |
| Lawlor | totalVI | macro-F1 | 0.888 | 0.012 | 50 fold×seed |
| Lawlor | scVI_matched | predictive AUROC | 0.801 | 0.020 | 50 |
| Lawlor | totalVI | predictive AUROC | 0.838 | 0.018 | 50 |
| Lawlor | scVI_matched | latent AUROC | 0.605 | 0.049 | 50 |
| Lawlor | totalVI | latent AUROC | 0.650 | 0.026 | 50 |
| Lawlor | scVI_matched | AURC | 0.087 | 0.018 | 50 |
| Lawlor | totalVI | AURC | 0.036 | 0.007 | 50 |

Footnote: 50 fold×seed pairs are descriptive consistency units, not 50 independent biological replicates.
