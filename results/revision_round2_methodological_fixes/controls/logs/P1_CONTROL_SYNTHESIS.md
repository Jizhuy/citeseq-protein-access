# P1 Control Synthesis

**Date:** 2026-09-11  
**Scope:** Parts A–C only. P0 frozen. No manuscript rewrite. No commit.

---

## Completion status

| Control | Status |
|---|---|
| A. Simple multimodal concat PCA | **DONE** (frozen) |
| B. RNA-only annotation control | **DONE** (frozen; provenance = already RNA-only) |
| C. Test-time protein corruption | **DONE** — Branch A on RTX 4060 using PHASE 11B adapted models |

---

## A–T answers

### A. Development HC l2 logreg macro-F1

| Representation | macro-F1 |
|---|---:|
| RNA PCA (10) | 0.666 |
| protein PCA (10) | 0.512 |
| **concat PCA 10+10** | **0.745** |
| MOFA+ | 0.549 |
| scVI_matched | 0.722 |
| totalVI | 0.779 |

### B. How much of totalVI−scVI gain does concat reproduce?

- Δ(totalVI − scVI) = **+0.057**
- Δ(concat − scVI) = **+0.023**
- Fraction of gain reproduced by concat ≈ **40%**
- Concat recovers most of the gap from RNA-only PCA toward totalVI, but not all.

### C. Does totalVI outperform concat PCA?

**Yes, modestly** on development HC: Δ(totalVI − concat) ≈ **+0.034**.

### D. RNA-only high-confidence labels

**9494 / 10849** cells (identical to historical HC).

### E. RNA-only vs historical agreement

**l1 = 1.0, l2 = 1.0** — because PHASE 6 query mapping was already RNA-only SCANVI/scArches (code-audited).

### F. Under RNA-only labels: scVI / totalVI / Δ

Identical to development (same labels): **0.722 / 0.779 / +0.057**.

### G. Under RNA-only: concat / Δ(totalVI−concat)

**0.745 / +0.034**.

### H. Does original totalVI advantage survive RNA-only annotation?

**Yes** — and it was never threatened by protein-informed *development* labels.

### I. Is development result partly explained by protein-informed label provenance?

**No for PBMC development.** Labels were RNA-derived.  
**Yes as a concern for Lawlor author labels** (protein-gated): protein PCA alone reaches macro-F1 ≈ **0.927** and concat ≈ **0.945**, exceeding P0 totalVI donor-mean present-class F1 ≈ **0.884**. That is expected when evaluation labels are protein-defined.

### J. Must the central claim change?

**Partially — STATE 2:**

- Keep: paired protein information improves cell-type discrimination.
- Soften: attributing the full gain to totalVI *architecture*.
- Development: totalVI still edges concat (+0.034) and scVI (+0.057).
- Transfer: **inductive concat exceeds both scArches-adapted VAEs** in A and B — strong caution against architecture-only storytelling.
- Lawlor: protein access dominates under protein-gated labels.

### K. Does clean-trained/fixed totalVI degrade under target-only protein corruption?

**Yes, gradually.** Mean ΔF1_from_clean at p=0.85: Dir A **−0.018**, Dir B **−0.033**. Monotonic decline; no collapse.

### L. At what corruption levels does totalVI remain above scVI?

- **Dir A:** clean totalVI already **slightly below** clean scVI (0.652 vs 0.656); remains below at all p.
- **Dir B:** totalVI **stays above** clean scVI reference (0.601) through p=0.85 (0.621).

### M. Does test-time corruption produce stronger degradation than training-time?

**Generally no** (with design caveat). Phase 8 training-time within-dataset HC ΔF1 at matched fractions is typically larger than Part C transfer test-time Δ_from_clean, except Dir B at p=0.85 (test Δ −0.033 vs train Δ −0.024). Designs are not identical.

### N. Which cell types are most sensitive?

Largest F1 drops at p=0.85 (support≥50): **CD8 Naive (Dir B, ≈−0.26)**, then CD8 TCM/TEM, CD4 Naive/TCM/TEM. Many myeloid/B subsets are nearly flat.

### O. Predictive / latent uncertainty under test-time corruption?

- Predictive median uncertainty **rises** with p (Dir B: 0.169 → 0.204).
- Predictive error AUROC **mildly declines** (Dir B: 0.856 → 0.841) — confidence becomes slightly less discriminative, not easier-to-spot errors.
- Latent error AUROC **slightly rises** (Dir B: 0.610 → 0.627); median latent uncertainty nearly flat.

### P. Does E-AURC worsen with corruption?

**Yes.** Dir A: 0.0315 → 0.0339; Dir B: 0.0350 → 0.0448. Full AURC and pAURC_[0.5,1] also worsen monotonically.

### Q. Model parameter change during corrupted inference?

**No.** `param_hash_unchanged` True for all 1514 all_runs rows; hash mismatch count 0. RNA unchanged True throughout.

### R. Target labels in preprocessing/training?

**No** for A/B PCA fits (train/source only). Class eligibility still uses target label *counts* for transfer design (same as PHASE10; not model fitting). Part C inference does not use target labels for training/adaptation.

### S. Do P1 controls materially change interpretation?

**Yes:**

1. Development advantage is **not** a protein-label artifact.
2. A large share of discriminative value is **protein access** (concat).
3. totalVI retains a **smaller residual advantage** over concat on development.
4. Transfer story must distinguish inductive multimodal features vs transductive VAE adaptation — concat can win under inductive PCA.
5. Lawlor author-label gains are partly **label–modality alignment**.
6. Post-adaptation test-time protein degradation causes **gradual**, not catastrophic, totalVI loss; reliability metrics worsen smoothly.

### T. Additional P2 required?

**No P2 experiment is required** by Part C. MultiVI / missingness / antibody dropout remain optional, not mandatory.

---

## PART C appendix (test-time corruption)

- **Venue:** RTX 4060 Laptop GPU, CUDA 12.6, torch 2.14.0+cu126, scvi-tools 1.3.3  
- **Branch:** A — frozen PHASE 11B adapted models at `results/phase11b_uncertainty_robustness/models/`  
- **Grid:** 20/20 seeds × both directions; totalVI full corruption×mask grid; scVI invariance confirmed (`max_abs=0.0`)  
- **Outputs:** `controls/test_time_corruption/{tables,logs,preview_figures,done/TESTTIME_CORRUPTION_DONE.json}`  
- **Clean baselines:** 80/80 CONSISTENT vs Phase 11B (abs_diff=0; same CUDA weights)

### Clean macro-F1 (mean ± SD, n=20)

| | scVI | totalVI |
|---|---:|---:|
| Dir A | 0.656 ± 0.018 | 0.652 ± 0.015 |
| Dir B | 0.601 ± 0.020 | 0.654 ± 0.025 |

### totalVI ΔF1 from clean (mean)

| p | Dir A | Dir B |
|---:|---:|---:|
| 0.10 | −0.0025 | −0.0023 |
| 0.25 | −0.0058 | −0.0069 |
| 0.40 | −0.0087 | −0.0108 |
| 0.55 | −0.0113 | −0.0173 |
| 0.70 | −0.0153 | −0.0253 |
| 0.85 | −0.0178 | −0.0328 |

First “meaningful” degradation (~|Δ|≥0.01): **p≈0.55 (A)** / **p≈0.40 (B)**.

---

## Transfer snapshot (inductive PCA vs transductive VAEs)

| Direction | concat | scVI | totalVI |
|---|---:|---:|---:|
| A | **0.735** | 0.656 | 0.652 |
| B | **0.734** | 0.601 | 0.654 |

Note: not an apples-to-apples architecture contest — concat has no query adaptation; VAEs use unsupervised transduction.

---

## Manuscript state classification

**STATE 2 — MULTIMODAL-INFORMATION CLAIM SURVIVES; ARCHITECTURE CLAIM WEAKENS**

Not STATE 4: test-time degradation is real but gradual, not rapid collapse.

Preferred framing:

> Paired protein measurements provide substantial discriminative value; totalVI adds selected further gains on the development HC benchmark and in reliability settings, but simple RNA+protein concatenation already captures much of the classification advantage and can dominate under inductive transfer and protein-gated external labels. After clean adaptation, target-only protein corruption degrades totalVI performance and selective-prediction metrics gradually rather than catastrophically.

---

## STOP

Per instructions: no MultiVI, no missingness, no antibody dropout, no manuscript rewrite, no final figures, no DOCX, no git commit/push.

Awaiting explicit approval for claim reframe / figures / any P2 work.
