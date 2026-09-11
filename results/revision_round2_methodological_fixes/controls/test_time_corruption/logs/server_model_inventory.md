# Server model inventory — Part C

- timestamp: 2026-09-11
- SERVER_PROJECT_ROOT: `/home/jizhu/research/multiomics_robustness`

## Decision: BRANCH A

Canonical PHASE 11B adapted query models are present:

`results/phase11b_uncertainty_robustness/models/{direction_A,direction_B}/{scvi_matched,totalvi}/seed00–seed19/{source_model,adapted_query_model}/model.pt`

| metric | value |
|---|---|
| adapted+source complete units | **80 / 80** |
| missing | 0 |
| canonical status | yes (PHASE 11B) |
| loadable (smoke) | yes (CUDA) |

### Sample validated (smoke seed0)

| path | phase | model | direction | seed | loadable |
|---|---|---|---|---|---|
| `.../direction_B/totalvi/seed00/adapted_query_model` | 11B | totalVI | B | 0 | yes |
| `.../direction_B/scvi_matched/seed00/adapted_query_model` | 11B | scVI | B | 0 | yes |
| `.../direction_A/totalvi/seed00/adapted_query_model` | 11B | totalVI | A | 0 | yes |
| `.../direction_A/scvi_matched/seed00/adapted_query_model` | 11B | scVI | A | 0 | yes |

### Not used for Part C primary grid

- PHASE 10 `seed*_cuda` dirs: source-only style artifacts (no `adapted_query_model` naming)
- PHASE 12A subset `adapted_query_model`: crossed subset design — not the 20-seed full-data grid

## Do not regenerate

Historical PHASE 11B adapted models are used as-is.
