# PHASE 12B setup_anndata / batch covariate documentation

## Decision (frozen before performance inspection)

| Field | Role | Rationale |
|---|---|---|
| `layer="counts"` | RNA UMI counts | Integer raw counts after Ensembl→symbol aggregation by **sum** |
| `batch_key="batch"` | `obs["run_identifier"]` (10 sequencing lanes) | Technical library/lane covariate from the pooled 10x design |
| `obs["donor"]` | **Not** used as `batch_key` | Donor is the held-out biological generalization axis |
| `obs["condition"]` | Fixed to Baseline | Stimulated conditions excluded from primary analysis |
| Protein | `obsm["protein_expression"]` (12 columns) | totalVI only; canonical names |

## What we deliberately do not do

- Do not regress out donor identity in `setup_anndata` (would collapse the validation question).
- Do not treat condition as a batch (Baseline-only).
- Do not use Seurat-v4 / GSE164378 labels.

## scArches

Source model trained on 8 donors; query adaptation on 2 held-out donors.
Both source and target are re-encoded with the **adapted** query model (PHASE 10/11B coordinate rule).
