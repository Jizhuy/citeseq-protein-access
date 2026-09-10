# Source-of-truth audit (PHASE 13)

Canonical artifacts for manuscript claims. **Do not cite obsolete files.**

## Global environment
| Item | Canonical path |
|---|---|
| Runtime versions | `results/logs/environment.json` |
| Pinned deps | `requirements.txt` |
| Key stack | Python 3.11.7; scvi-tools 1.3.3; scanpy 1.10.3; torch 2.14.0 (+cu126 on GPU) |

## Representation benchmark (internal PHASE 6–7)
| Use | Path |
|---|---|
| All-model metrics (PCA / MOFA+ / scVI_matched / totalVI) | `results/tables/mofa_biological_representation_metrics.csv` |
| Pairwise deltas scVI↔totalVI | `results/tables/primary_scvi_matched_vs_totalvi.csv` |
| Pairwise deltas MOFA↔totalVI | `results/tables/primary_mofa_vs_totalvi.csv` |
| Annotation provenance | PHASE 6 Seurat-v4 reference transfer logs under `results/annotation/` / phase6 figures |

**Do not use:** `results/tables/biological_representation_metrics.csv` (incomplete; missing PCA/MOFA+).

## Protein entry-wise corruption (PHASE 8)
| Use | Path |
|---|---|
| Primary delta by corruption | `results/tables/protein_sparsity_primary_delta_summary.csv` |
| Full summary JSON | `results/logs/phase8_full_summary.json` |

**Do not use:** `phase8_smoke_summary.json`.  
**Language rule:** entry-wise corruption ≠ missing modality.

## Cross-dataset transfer (PHASE 10)
| Use | Path |
|---|---|
| **CORRECTED** primary table | `results/phase10_cross_dataset/tables_corrected/cross_dataset_primary_corrected.csv` |
| **CORRECTED** delta summary | `results/phase10_cross_dataset/tables_corrected/cross_dataset_primary_delta_summary_corrected.csv` |
| Erratum | `results/phase10_cross_dataset/logs/phase10_erratum_latent_alignment.md` |

**Do not use:**
- `results/phase10_cross_dataset/tables/cross_dataset_primary*.csv` (mixed pre-/post-adaptation coordinates for totalVI seeds 1–4)
- `results/phase10_cross_dataset/logs/phase10_results.md` (uncorrected narrative)

For Direction A interpretation in the manuscript, prefer **PHASE 11B 20-seed** macro-F1 (near-null) over the corrected 5-seed wording “totalVI worse.”

## Source-size sensitivity (PHASE 10B)
| Use | Path |
|---|---|
| Summary table | `results/phase10b_source_size/tables/source_size_summary.csv` |
| Narrative | `results/phase10b_source_size/logs/phase10b_results.md` |

Not contaminated by the PHASE 10 coordinate bug.

## Uncertainty — historical vs primary
| Role | Path |
|---|---|
| Historical pilot (5-seed) | `results/phase11_uncertainty/` — **secondary / historical only** |
| **Primary uncertainty** | `results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv` |
| 20-seed summary | `.../phase11b_20seed_summary.csv` |
| Paired deltas + evidence labels | `.../phase11b_seed_delta_summary.csv` |
| Selective prediction | `.../phase11b_selective_prediction.csv` |
| Retention imbalance | `.../phase11b_retention_imbalance.csv` |
| Fragility associations | `.../phase11b_uncertainty_instability_assoc.csv` |
| Methods / versions | `.../logs/phase11b_methods.md` |

## Replicated variance (PHASE 12A)
| Use | Path |
|---|---|
| Metric variance components | `results/phase12a_variance_replication/tables/replicated_variance_components.csv` |
| Delta variance components | `.../replicated_delta_variance_components.csv` |
| Final report | `.../logs/phase12a_final_report.md` |
| Novelty audit | `.../logs/method_novelty_audit.md` |

**Do not** use PHASE 11B lumped “interaction/residual ≈ 80%” as the final interpretation.

## External validation (PHASE 12B)
| Use | Path |
|---|---|
| Primary runs | `results/phase12b_external_validation/formal_validation/tables/phase12b_primary_runs.csv` |
| Deltas | `.../phase12b_model_deltas.csv` |
| Variability | `.../phase12b_variability.csv` |
| Replication table | `.../phase12b_replication_summary.csv` |
| Report | `.../phase12b_report.md` |
| Labels | Lawlor **author** labels only |

**Do not** use Hao / Seurat-v4 / GSE164378 labels for Lawlor.

## Missing-modality API (PHASE 9)
No experimental results table. Document as Methods/Discussion limitation only.
