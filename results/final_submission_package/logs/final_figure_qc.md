# Final figure QC

## Story check (one sentence each)

| Figure | Scientific question |
|---|---|
| Figure 1 | What reliability stresses and cohorts were evaluated? |
| Figure 2 | How do representation choice and entry-wise protein degradation affect biological performance? |
| Figure 3 | How does multimodal benefit behave under reciprocal dataset shift and source-size context? |
| Figure 4 | What distinct information is carried by predictive confidence versus latent posterior variance? |
| Figure 5 | How do stochastic fragility and selective prediction trade risk against biological retention? |
| Figure 6 | Do the central findings reproduce under independent donor-held-out validation? |

**Redundancy:** No two figures answer the same primary question. Figure 4 and Figure 6 both show uncertainty metrics, but Figure 6 is the external-validation capstone (different cohort/design). Keep both.

## Technical QC

| Item | Status |
|---|---|
| Panel letters A/B/C… | Present on composites |
| PNG + PDF present | Yes (`figures/Figure1`–`Figure6`) |
| Model names consistent in legends | scVI_matched / totalVI |
| PHASE numbers in legends | Absent |
| Obsolete PHASE10 numbers | Absent |
| Entry-wise vs missing-modality language | Correct in Figure 2 legend |
| Lawlor author-label independence | Stated in Figure 6 legend |
| Uncertainty bars / n explained | Described in legends; composite panels inherit source-plot annotations |
| Single-/double-column readability | **Requires human visual inspection** (composites are large; may need reflow for journal upload) |
| Data altered? | No |

## Remaining human figure tasks
- [ ] Inspect each PDF at ~85 mm and ~170 mm width
- [ ] Confirm no clipped labels after publisher downsampling
- [ ] Optionally re-layout dense Figure 3/4 panels for typography only
