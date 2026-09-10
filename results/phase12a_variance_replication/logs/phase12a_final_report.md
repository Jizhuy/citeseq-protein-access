# PHASE 12A Final Report (A–T)

Date: 2026-09-10  
Status: **COMPLETE**. Integrity audit **PASS**. No new model. No commit/push.

Server verification: 50/50 replicate-1 done-markers; 0 pending; 0 failures; CUDA available; PHASE 11B intact (130 done-markers). Training was already finished — no retrain. Analysis, figures, and interpretation were re-run after the integrity gate.

---

## A. Final run count and failures

**50 / 50** replicate-1 runs complete. **0 failures.**  
Combined with 50 PHASE 11B crossed observations (replicate 0) → **100 total observations**.

Integrity audit (`logs/integrity_audit.json`): PASS for A–J (models, manifests, latents, posterior, probabilities, reload validation, finite values, matching post-adaptation dims, no orphans).

## B. Exact replicate design

Direction A only: **PBMC10k → PBMC5k**

```
5 source subsets  ×  5 model seeds  ×  2 run replicates  ×  2 models
= 100 observations
```

- Replicate 0 = existing PHASE 11B path (`run_replicate_seed = None`)
- Replicate 1 = PHASE 12A (`run_replicate_seed = 90000 + 1000·rep + 10·subset + model_seed`)
- Within each subset×seed cell: identical source cells, architecture, hyperparameters, labels; only run-level stochasticity redrawn for replicate 1

## C. Replicate 1 is genuinely stochastic (not a deterministic duplicate)

**Confirmed.**

Determinism probe: same `model_seed` + identical full RNG → bit-identical (max abs latent diff = 0).  
Manifests: `model_seed` held fixed within each subset×seed pair; `run_replicate_seed` is None for rep0 and equals the designed nonzero seed for rep1 (e.g. subset0/seed0 → 91000). Architecture/hparams/subset identical. Same-seed duplicates would have produced zero residual — that did not happen.

## D–G. Metric-level variance fractions

Mean non-negative fractions over key metrics (macro-F1, error AUROC, Brier, ECE, AURC). Raw MoM estimates are preserved in the CSV (including negatives).

| Model | D source | E seed | F interaction | G residual |
|---|---:|---:|---:|---:|
| scVI_matched | **0.036** | **0.054** | **0.252** | **0.658** |
| totalVI | **0.005** | **0.071** | **0.124** | **0.800** |

Metric-specific notes:

- scVI macro-F1: interaction 0.355 (p=0.043), residual 0.622 — real interaction *and* large residual
- totalVI macro-F1: residual **0.911**, interaction 0.032 (p=0.43) — residual-dominated; non-significant interaction means no interaction of *detectable* magnitude for this design, not proof of zero interaction
- totalVI error AUROC: interaction 0.302 (p=0.037), residual 0.498 — one of the stronger interaction signals

Power caveat (preserved): Type I ≈ 0.044; ~78% power when interaction = residual; ~45% when interaction = ½ residual; ≥80% power only near interaction ≈ **2×** residual. Non-significant tests bound interaction below roughly that detectable magnitude.

## H. Comparison with PHASE 11B’s ~80% interaction/residual

**Mostly pure run-level residual stochasticity**, with a minority true subset×seed interaction.

| Model | Residual share of (interaction + residual) |
|---|---:|
| scVI_matched | **0.72** |
| totalVI | **0.87** |
| Median across model×metric | **~0.74** |

So the previous unidentified ~80% term was **not** predominantly true interaction. Metric-specific exceptions exist (scVI macro-F1 / AURC; totalVI error AUROC / AUPRC), where interaction is non-negligible and often F-significant.

## I–J. Delta (totalVI − scVI) variance decomposition

Mean over key metrics:

| Component | Fraction | Dominates? |
|---|---:|---|
| Source subset | 0.021 | no |
| Model seed | 0.090 | no |
| Interaction | 0.202 | secondary |
| **Residual** | **0.687** | **yes (J)** |

Absolute non-negative variances are in `replicated_delta_variance_components.csv` (e.g. macro-F1: residual var ≈ 4.97e-4 of total ≈ 7.23e-4). The multimodal Δ advantage is itself noisy at the run level.

## K. Per-cell median variance fractions (true-class probability)

| Model | Source | Seed | Interaction | Residual |
|---|---:|---:|---:|---:|
| scVI_matched | 0.101 | 0.011 | 0.000 | **0.776** |
| totalVI | 0.055 | 0.020 | 0.010 | **0.808** |

## L–O. Cell-type extremes (median fractions)

**L. Highest residual instability**  
- scVI: CD16 Mono (0.92), NK (0.84), CD4 TEM (0.82), Eryth (0.81), cDC2 (0.81)  
- totalVI: NK (0.92), CD8 TEM (0.88), cDC2 (0.86), CD4 Naive (0.83), CD16 Mono (0.83)

**M. Strongest interaction**  
- scVI: CD8 Naive (0.11), CD8 TEM (0.07), B memory (0.07)  
- totalVI: **CD8 TCM (0.21)**, **Eryth (0.20)**, MAIT (0.10), gdT (0.06)

**N. Most source-composition dominated**  
- scVI: B intermediate (0.36), B memory (0.29), **gdT (0.21)**, **ILC (0.16)**  
- totalVI: B intermediate (0.23), **ILC (0.22)**, B naive (0.13)

**O. Most model-seed dominated**  
- scVI: NK (0.035), CD8 TEM (0.031), **NK_CD56bright (0.030)**  
- totalVI: CD14 Mono (0.071), CD16 Mono (0.038), CD8 Naive (0.038), **ILC (0.036)**

Previously unstable classes (gdT, ILC, NK_CD56bright, Treg, CD8 TCM, Eryth): residual still large in all; CD8 TCM and Eryth show elevated *interaction* under totalVI; gdT/ILC show elevated *source* under scVI.

## P. Does 5×5×2 change the PHASE 10B / 11B interpretation?

**Yes, materially.**

1. PHASE 11B’s lumped ~80% is now identified as **mostly run-level residual**, not irreducible “subset×seed biology.”  
2. Source-subset main effects are **small** once residual is separated — PHASE 10B source-size sensitivity should not be over-read as the dominant uncertainty driver for this crossed design.  
3. Model-seed effects remain modest relative to residual.  
4. Some real interaction exists and is metric-specific; it must not be dismissed, but it is not the bulk of the old term.  
5. Practical implication: single-run point estimates (and single-seed Δ(totalVI−scVI)) are less stable than the design previously made visible.

## Q. Final novelty-audit conclusion (not re-run)

**Posterior precision weighting is not novel.** Gaussian PoE already implements it (Wu & Goodman 2018); Cobolt, scMaui, and UniVI ship it in single-cell settings. Cell-specific weighting is also prior art (including MultiVI’s tested `w_{m,c}`). Do not propose simple precision-weighted fusion as a contribution.

## R. Strongest remaining methodological / evaluative gap

**Uncertainty-guided selective prediction with explicit per-class retention reporting**, under domain shift, with calibration/reliability metrics (ECE/Brier/error-detection) reported separately from discrimination.

Secondary: source-composition robustness diagnostics remain open in the literature, but PHASE 12A shows source-subset variance is small relative to run-level residual here — so that gap is less urgent *for this empirical setting* than selective prediction / reliability reporting.

## S. Is a new model justified?

**No new fusion architecture.** Especially not precision-weighted PoE.  
A new *evaluation / reliability / selective-prediction* contribution remains justified by PHASE 8–11B failures plus the novelty audit.

## T. Scientifically justified next phase (ranked; not implemented)

| Rank | Option | Justification |
|---:|---|---|
| **1** | **A. Uncertainty-aware selective prediction with per-class retention constraints** | Directly addresses PHASE 11B’s abstention fairness failure; empty in the audited literature; residual-dominated instability makes reliability/deferral more actionable than fusion. |
| **2** | **C. Calibration-aware reliability analysis under dataset shift** | Still nearly empty for multi-omic latents; supports A; UniVI names calibration as future work. |
| **3** | **E. Stop model development and begin manuscript consolidation** | Strong option if the thesis is empirical: residual-dominated variability + precision-fusion prior art is itself a publishable negative/clarifying result. |
| **4** | **B. Source-composition robustness / fragility diagnostics** | Literature gap remains, but PHASE 12A source fractions are tiny — weaker empirical motivation *here* than A/C. |
| **5** | **D. A genuinely new architecture** | Not justified now; any future architecture must clear the six gates in `method_novelty_audit.md` and must not rediscover PoE. |

**Recommended default:** pursue **A** (with **C** as the measurement backbone), or choose **E** if consolidating the PHASE 8–12A empirical story is the priority. Wait for explicit approval before starting either.

---

## Artifacts

- Integrity: `logs/integrity_audit.json`  
- Tables: `replicated_variance_components.csv`, `replicated_delta_variance_components.csv`, `per_cell_replicated_variance.csv`, `celltype_replicated_variance.csv`, `interaction_power_check.csv`  
- Figures: `figure01_replicated_experimental_design`, `figure01b_old_vs_new_components`, `figure02`–`figure05`  
- Logs: `phase12a_interpretation.md`, `replicate_design.md`, `method_novelty_audit.md` (unchanged)  
- Replicate-0 check: max abs diff vs PHASE 11B ≈ 1e-16

## STOP

No PHASE 12B. No new architecture. No MultiVI. No RNA+ATAC. No commit. No push.
