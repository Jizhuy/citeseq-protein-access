# Manuscript claim audit

## Scope
Audited `manuscript/manuscript_v1.md` against canonical tables listed in `logs/source_of_truth.md` and `tables/manuscript_evidence_ledger.csv` (n=110 ledger rows).

## Numeric verification (primary)
| Claim | Manuscript | Source | Status |
|---|---|---|---|
| totalVI l2 F1 | 0.779 | primary_scvi_matched_vs_totalvi / mofa metrics 0.7792 | PASS |
| scVI_matched l2 F1 | 0.722 | 0.7219 | PASS |
| Delta + bootstrap CI | +0.057 [0.020, 0.091] | +0.0566 [0.0200, 0.0912] | PASS |
| Dir A 20-seed F1 | 0.656 / 0.652 | 0.6556 / 0.6524 | PASS |
| Dir B 20-seed F1 | 0.602 / 0.654 | 0.6015 / 0.6540 | PASS |
| Dir A latent AUROC | 0.544 / 0.669 | 0.5439 / 0.6694 | PASS |
| Dir B latent AUROC | 0.469 / 0.610 | 0.4688 / 0.6099 | PASS |
| Lawlor F1 | 0.800 / 0.888 | 0.7998 / 0.8883 | PASS |
| Lawlor latent AUROC | 0.605 / 0.650 | 0.6049 / 0.6497 | PASS |
| Lawlor paired Δ F1 | +0.088; 50/50 | model_deltas | PASS |

## Contamination checks
| Risk | Result |
|---|---|
| Obsolete PHASE 10 mixed-coordinate Δ (−0.0102 / totalVI 0.6504) | **Not present** |
| Using PHASE 11 5-seed as primary | **Not present** (20-seed used) |
| Hao/Seurat-v4 as Lawlor labels | **Not present** |
| Entry-wise corruption called missing modality | **Not present** |
| ECE treated as discrimination | **Not present** |
| “Irreducible biological noise” as claim | Present only as negation in v1 → **removed/rephrased in v2** |
| Internal PHASE numbering in body | v1 had 2 residual mentions → **removed in v2** |

## Unsupported / softened language in v2
- Removed residual internal phase references from Data availability / title note.
- Replaced “sign flip” wording even in cautionary sentence.
- Rephrased residual-variance caution without “irreducible biological noise.”
- Did **not** promote Direction A multimodal disadvantage.
- Did **not** claim causal protein→calibrated posterior mechanism.
- Did **not** claim architecture novelty / PoE novelty.

## Claims removed or never included
- Universal totalVI superiority
- Precision-weighted fusion as contribution
- Stimulation / 39-ADT results
- Fabricated cell-wise missing-protein results
- New p-values not precomputed

## Evidence ledger
- 110 mapped claims
- Support levels: mostly A; B/D/E reserved for conditional, exploratory, and limitation items

## Output
Claim-checked prose: `manuscript/manuscript_v2_claim_checked.md`
