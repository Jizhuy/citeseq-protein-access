# MANUSCRIPT_FINAL_AUDIT.md

**Date:** 2026-09-24  
**Primary protocol:** strict train-only (7,573 / 1,921)

## Primary internal numbers

| Quantity | Value |
| --- | ---: |
| RNA PCA | 0.667 |
| concat | 0.745 |
| scVI_matched | 0.709 |
| totalVI | 0.741 |
| G_access | +0.078 [0.042, 0.131] |
| G_assoc | −0.004 [−0.060, 0.030] (prop>0=0.355) |
| Gap_total | +0.074 [0.029, 0.130] |
| R_access | ≈1.05 [0.680, 2.304] (secondary descriptive only) |

## Framing checks

- [x] Abstract uses strict train-only primary values only
- [x] R_access removed from Abstract emphasis
- [x] Figure 1 emphasizes G_access / G_assoc (not R_access recovery)
- [x] Figure 2 regenerated with train-only values
- [x] Tables 1–2 primary = train-only
- [x] Old transductive values only in Suppl. information-access sensitivity
- [x] Conclusion: near-null residual contrast (not “concat superior”)
- [x] Lawlor donor-level unchanged (0.690 / 0.923 / 0.942 / 0.794 / 0.884; Δ=+0.090; P=0.00195)
- [x] Degradation 0.779/0.722 retained as historical degradation experiment (not primary internal)

## Deliverable

`manuscript/CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx`
