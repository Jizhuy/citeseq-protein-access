# STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md

**Status:** PRIMARY internal analysis (final polish, 2026-09-24)

## Protocol

- Train: 7,573 HC cells only (PCA + scVI_matched + totalVI)
- Encode/evaluate: 1,921 HC held-out cells only
- Classifier: train labels only
- Bootstrap: B=10,000, seed 20260923, cell-level conditional paired; class-stratified sensitivity also run
- Artifacts: `strict_train_only/` (tables committed; models gitignored)

## Primary vs historical

| Quantity | Transductive (historical) | Strict train-only (PRIMARY) |
| --- | ---: | ---: |
| RNA PCA | 0.666 | **0.667** |
| concat | 0.745 | **0.745** |
| scVI_matched | 0.722 | **0.709** |
| totalVI | 0.779 | **0.741** |
| G_access | +0.079 [0.043, 0.132] | **+0.078 [0.042, 0.131]** |
| G_assoc | +0.034 [−0.012, 0.073] | **−0.004 [−0.060, 0.030]** (prop>0=0.355) |
| Gap_total | +0.113 | **+0.074 [0.029, 0.130]** |
| R_access | ≈0.70 | ≈1.05 [0.680, 2.304] (**secondary only**) |

## Interpretation checklist

1. G_assoc still crosses zero? **Yes** (near null).
2. R_access as recovery fraction? **No** — secondary descriptive ratio only; removed from Abstract/Fig1 emphasis.
3. Central conclusion changed? **Strengthened, not reversed** — protein-access stable; residual totalVI contrast near null under matched access.
4. Claim “concat better than totalVI”? **No** — near-null / imprecisely estimated residual contrast.
5. Lawlor unchanged? **Yes** (donor-level).
