# Revision Round 2 — Methodological Progress Report

**Tree:** `results/revision_round2_methodological_fixes/`  
**Git:** no commit/push performed.  
**Historical PHASE 6–14 artifacts:** untouched.

---

## Execution status

| Step | Status |
|---|---|
| 0 Design ledger | **DONE** `audit/experiment_design_ledger.csv` |
| 1 Leakage/transduction audit | **DONE** `audit/transduction_leakage_audit.md` |
| 2 Lawlor 50 vs 19/20 | **DONE** |
| 3 Donor-level Lawlor | **DONE** |
| 4 Variance decomposition | **DONE** (verified, not rewritten as “wrong”) |
| 6–7 Statistical units + AURC/reliability | **DONE** for Lawlor; PBMC summaries from tables |
| 5 RNA-only labels | **NOT RUN** (needs reference transfer) |
| 5B Simple concat PCA baseline | **NOT RUN** (queued) |
| 6 Test-time corruption | **NOT RUN** (needs GPU / frozen models) |
| 10 Modern baseline | **FEASIBILITY DONE — STOP before training** |
| 11 Missingness | **FEASIBILITY DONE — cell-wise unsupported** |
| 14–15 Figures/manuscript rewrite | **STOPPED** until P1 controls freeze |

---

## Final checklist answers (A–AH) — current freeze

### A. Variance decomposition mathematically wrong?
**No.** Historical fit was **per-model MoM** plus a separate **Delta** cube. Possible manuscript ambiguity only.

### B. Corrected within-model variance model
`Y_ijr^(a) = μ_a + S_i^(a) + G_j^(a) + (SG)_ij^(a) + ε_ijr^(a)`  
fit separately for scVI_matched and totalVI (5×5×2).

### C. Corrected Delta variance model
`Δ_ijr = Y_totalVI − Y_scVI = μ_Δ + S_iΔ + G_jΔ + (SG)_ijΔ + ε_ijrΔ`

### D. REML vs MoM
Primary = MoM EMS. For this balanced equal-replication design, MoM ≡ REML under normality. Numerical MixedLM with 5×5 levels is fragile; not elevated over MoM.

### E. Correct Lawlor latent sign-consistency /50
**38/50** totalVI > scVI (latent error AUROC).  
Macro-F1: **50/50**. Predictive AUROC: **45/50**.

### F. Source of historical 19/20
**Internal PBMC transfer**, not Lawlor — Direction B macro-F1 (and latent AUROC) on the **20-seed** grid. Confirmed again from `phase11b_20seed_primary.csv`.

### G. Lawlor per-donor present-class macro-F1 (seed-means)

| Donor | scVI | totalVI |
|---|---:|---:|
| Donor01 | 0.815 | 0.897 |
| Donor02 | 0.772 | 0.882 |
| Donor03 | 0.812 | 0.897 |
| Donor04 | 0.789 | 0.884 |
| Donor05 | 0.785 | 0.882 |
| Donor06 | 0.779 | 0.871 |
| Donor07 | 0.780 | 0.857 |
| Donor08 | 0.782 | 0.878 |
| Donor09 | 0.812 | 0.894 |
| Donor10 | 0.810 | 0.896 |

Primary definition = **present-class** (frozen before deltas). Fixed-ontology also saved.

### H. Lawlor per-donor ΔF1 (totalVI − scVI)
Range ≈ **+0.077 to +0.110**; mean **+0.090**; median **+0.089**.

### I. Donors favoring totalVI
**10 / 10**

### J. Between-donor variation
SD of donor ΔF1 ≈ **0.0096**

### K. Within-donor seed variation
Mean within-donor seed SD (present-class F1): scVI ≈ **0.013**; totalVI ≈ **0.0085**

### L. Target labels leaked into representation training?
**No.**

### M. Target features used transductively?
**Yes** — unlabeled target in scArches query adaptation (Lawlor + PBMC transfer).

### N. Correct terminology
“**Donor-held-out / dataset-held-out label evaluation after unsupervised query adaptation**”  
(not fully inductive unseen-target prediction).

### O. Do multimodal labels favor totalVI?
**Open — P1.** Development labels are multimodal-reference-derived; Lawlor labels are protein-gated. RNA-only control not yet run.

### P–R. RNA-only / concat PCA / totalVI vs concat
**Not yet run.**

### S. Original training-time corruption
Keep; rename to **training under progressive protein-entry degradation**.

### T. Clean-trained → corrupted-test
**Not yet run.**

### U–W. AURC (Lawlor means over 50 runs/model)
| | scVI | totalVI |
|---|---:|---:|
| Full AURC | 0.0867 | 0.0362 |
| pAURC_[0.5,1.0] | 0.1571 | 0.0670 |
| E-AURC | 0.0545 | 0.0262 |

**Important:** historical `aurc` column was already **full AURC**, not partial. Coverage grid is for tables only.

### X. NLL/Brier/ECE/AUPRC
See `calibration/full_reliability_metrics.csv` and Lawlor pairwise CSV.

### Y. Modern multimodal baseline feasibility
**MultiVI not fair** (RNA+ATAC). Do not force.

### Z. Additional modern model run?
**No.**

### AA. Honest cell-wise missing protein?
**Unsupported** under totalVI 1.3.3.

### AB. Antibody-panel dropout
**Not yet run.**

### AC–AF. Claim changes / figures / conclusion shift
Deferred until P1 controls complete. Already required language fixes: statistical units, transduction wording, AURC naming clarification, variance description clarity, 19/20 provenance.

### AG. Still REQUIRED before manuscript rewrite?
1. Simple RNA+protein concat PCA baseline (P1, low-cost)  
2. RNA-only annotation control (P1)  
3. Clean-trained → corrupted-test (P1 robustness)  
Optional P2: antibody-panel dropout; one native RNA+protein modern model if installable fairly.

### AH. Application manuscript DOCX
**Not generated yet** (Part 15 gated on frozen corrected analyses).

---

## Hard-stop evaluation

| Condition | Status |
|---|---|
| Genuine target-label leakage | **Not found** — continue |
| Mixed source/target embeddings | Not reintroduced in Lawlor/11B paths |
| Lawlor cell→donor mapping | **OK** via annotations + donor map |
| 100 Lawlor runs reconstructible | **OK** |
| Variance models unstable | MoM stable; REML optional |
| RNA-only protein-independence | **Pending** — stop before claiming O–R |
| Concat PCA target-fitted preprocess | Must enforce source-only fit when run |
| AURC reconstructible | **OK** |
| Missingness zero-fill masquerade | **Refused** |
| Modern baseline fairness | **STOP — MultiVI unfair** |

---

## Next actions (ordered)

1. On RTX 4060 (`/home/jizhu/research/multiomics_robustness`): implement & run PART 5B concat PCA + PART 6 test-time corruption with tmux/done-markers under this tree.  
2. Implement PART 5A RNA-only labels; freeze before scoring.  
3. Only then regenerate figures and application manuscript.

**Awaiting approval before git commit and before starting GPU training jobs.**
