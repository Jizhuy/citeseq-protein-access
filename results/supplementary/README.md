# Supplementary results

| Folder | Content | Script |
| --- | --- | --- |
| `information_access_sensitivity/` | **SECONDARY INFORMATION-ACCESS SENSITIVITY** (historical feature-transductive analysis) | see folder README |
| `controls/tables/` | Simple PCA baselines (development, transfer, Lawlor), RNA-only label controls; includes the historical MOFA+ development value | `scripts/run_simple_multimodal_baselines.py` |
| `marker_channel_dropout/` | Marker-channel dropout sensitivity, **PCA only** | `scripts/run_marker_channel_dropout.py` |
| `environment/runtime_environment.json` | Software versions observed at run time | `scripts/print_environment.py` |

Not run: totalVI retraining under complete marker-channel dropout, and an additional newly trained multimodal comparator.
