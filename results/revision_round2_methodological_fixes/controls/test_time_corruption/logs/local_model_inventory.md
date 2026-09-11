# Local model inventory — Part C

**LOCAL_PROJECT_ROOT:**  
`/Users/venus/Desktop/2026北京大学科学研究/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness`

**Search date:** 2026-09-11  
**Branch decision:** **BRANCH B** — dedicated local clean-adapted models required.

---

## Search scope

- Current project root (full recursive)
- Parent research workspace
- `results/phase10*`, `phase11*`, `phase11b*`, `phase12*`
- `results/models/`, backup-like paths

---

## PHASE 10 / 11B adapted query models

| Pattern | Found? |
|---|---|
| `adapted_query_model/` | **0** |
| `query_model/` under phase10/11b | **0** |
| `source_model/` under phase10/11b | **0** |
| `phase11b_uncertainty_robustness/models/**/model.pt` | **0** (tables/logs only; weights gitignored) |
| `phase10_cross_dataset/models/**/model.pt` | **0** |

**Canonical historical PHASE 10/11B adapted query models are NOT available locally.**

---

## Other local `model.pt` candidates (not usable as Part C adapted query)

| path | model | direction | seed | source/target | adapted_query | canonical Phase10/11B? | usable for Part C? |
|---|---|---|---|---|---|---|---|
| `results/models/scvi_matched_seed0` | scVI_matched | combined PBMC | 0 | joint train | no | no (PHASE 5) | no — not source-only + query-adapted |
| `results/models/totalvi_seed0` | totalVI | combined PBMC | 0 | joint train | no | no (PHASE 4) | no |
| `results/models/scvi_seed0` | scVI_default | combined | 0 | joint | no | no | no |
| `results/models/annotation_reference` | SCANVI ref | annotation | — | — | no | no (PHASE 6) | no |
| `results/models/stress_sparsity/corruption_*` | totalVI | training-time corruption | 0–4 | joint corrupted | no | no (PHASE 8) | no — wrong experimental design |

All stress_sparsity / combined models: `source_model_present=N/A`, `adapted_query_model_present=false`, `manifest_present` varies, `loadable` likely yes but **not** Phase 10 transfer artifacts.

---

## Decision

Proceed under **BRANCH B**:

Create new clean-adapted models under:

`results/revision_round2_methodological_fixes/controls/test_time_corruption/local_dedicated_models/`

Scientific frame:

> Dedicated repeated clean-adaptation experiment for evaluating target-only protein degradation.  
> Primary inference = **within-run Δ_from_clean**, not bitwise match to historical CUDA PHASE 10.
