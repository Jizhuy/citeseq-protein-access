# Controls status (PART 5–6 / P1 A–C)

## PART 5A / P1 B — RNA-only annotation labels

**Status:** DONE (frozen). PHASE 6 mapping already RNA-only; n=9494; agreement=1.0.

## PART 5B / P1 A — Simple multimodal baseline (RNA PCA + protein PCA)

**Status:** DONE (frozen). Dev HC: concat 0.745; totalVI 0.779 (Δ≈+0.034).

## PART 6 / P1 C — clean-trained → corrupted-test

**Status:** DONE on RTX 4060 (Branch A, PHASE 11B adapted models).

- 20/20 seeds both directions
- DONE marker: `controls/test_time_corruption/done/TESTTIME_CORRUPTION_DONE.json`
- Synthesis: `controls/logs/P1_CONTROL_SYNTHESIS.md`
- Distinct from historical PHASE8 training-time degradation

## Execution venue

Completed on RTX 4060 server (`ssh 4060-server`). Artifacts synced to Mac revision tree under `controls/test_time_corruption/`.
