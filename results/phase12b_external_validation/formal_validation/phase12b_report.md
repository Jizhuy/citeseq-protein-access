# PHASE 12B — Independent External CITE-seq Validation (Baseline, DESIGN A)

Lawlor et al. *Front. Immunol.* 2021 donor-held-out validation.
Primary endpoint: **macro-F1**. Class-aware selective prediction is a baseline only (not claimed novel).

**Headline:** Across 100/100 runs, totalVI improves donor-held-out biological transfer (Δmacro-F1 ≈ +0.088, 50/50 pairs), predictive uncertainty remains strongly error-informative, and totalVI latent posterior uncertainty remains *directionally* more informative than scVI — with fold variability often larger than seed variability.

---

## FINAL REPORT — A to AK

### A. Final processed Baseline cell count
**16,175** Baseline SNG cells (`lawlor_baseline_sng_shared.h5ad`).

### B. Author-labeled Baseline count
**5,207** cells with author protein-gated labels (`obs["author_cell_type"]`).
Unlabeled Baseline cells retained for unsupervised VAE training.

### C. Mapped / shared gene count
**12,776** shared genes (Ensembl GRCh37.87 → symbols; duplicate symbols summed).
External Ensembl IDs: 32,738 → mapped 32,709 (99.91%); 55 duplicate symbols aggregated.

### D. 12-protein panel confirmation
Primary panel = **12** proteins (CD15, TIGIT missing from Lawlor).
Exact: CD3, CD4, CD8a, CD14, CD16, CD56, CD19, CD25, CD45RA, CD45RO.
Aliases: PD-1 ↔ CD279_PD1; CD127 ↔ CD127_IL7Ra.
Confirmed in all 100 runs: `n_proteins=12`.

### E. Donor fold definitions
10 donors → 5 folds of consecutive pairs (frozen before performance inspection):

| Fold | Target donors | Source donors |
|---|---|---|
| 0 | Donor01–02 | 03–10 |
| 1 | Donor03–04 | others |
| 2 | Donor05–06 | others |
| 3 | Donor07–08 | others |
| 4 | Donor09–10 | others |

Each donor is target exactly once. `batch_key` = `run_identifier` (lane); **donor is not batch**.

### F. Primary eligible cell types
Eligibility: ≥20 labeled SOURCE and ≥10 labeled TARGET in **every** fold.
Eligible (7/7 author classes): **B, CD14_Mono, CD4T_Mem, CD4T_Naive, CD8T_Mem, CD8T_Naive, NK**.
Removed: none.

### G. 100-run success / failure count
**100 / 100** complete done-markers; **0** failures; integrity_ok = True for all.
Gates 1–3 validated before seeds 1–9.

### H. macro-F1 scVI by fold and overall
Overall mean **0.800** (SD 0.015).

| Fold | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| scVI | 0.803 | 0.809 | 0.790 | 0.783 | 0.815 |

### I. macro-F1 totalVI by fold and overall
Overall mean **0.888** (SD 0.012).

| Fold | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| totalVI | 0.893 | 0.895 | 0.884 | 0.871 | 0.899 |

### J. Paired macro-F1 Delta (totalVI − scVI)
Mean **+0.088**; median **+0.090**; range [0.068, 0.118].
Sign consistency: **50/50** pairs totalVI higher.

### K. Predictive error AUROC — scVI
Mean **0.801** (fold means 0.774–0.821).

### L. Predictive error AUROC — totalVI
Mean **0.838** (fold means 0.822–0.856).
Paired Δ mean **+0.037**; totalVI higher in **45/50** pairs.

### M. Predictive uncertainty replication verdict
**Replicated.** Both models show strong U_pred → error discrimination under independent donor shift (AUROC ≫ 0.5). totalVI slightly stronger on average.

### N. Latent posterior AUROC — scVI
Mean **0.605** (high seed volatility within folds; within-fold seed SD ≈ 0.047).

### O. Latent posterior AUROC — totalVI
Mean **0.650** (within-fold seed SD ≈ 0.023).

### P. totalVI − scVI latent uncertainty Delta
Mean **+0.045**; totalVI higher in **38/50** pairs (76%).
Spearman(U_latent, −log p_true): scVI 0.35; totalVI 0.54.
Spearman(U_pred, U_latent): scVI 0.37; totalVI 0.54.

### Q. Did PHASE 11B latent uncertainty result replicate?
**Directionally yes** (totalVI U_latent more error-informative than scVI), with smaller absolute AUROCs than the original PBMC transfer setting and non-trivial seed noise for scVI. Treat as **replicated in direction, attenuated in magnitude**.

### R. NLL / Brier / ECE comparison
| Metric | scVI | totalVI | Δ (t−s) |
|---|---|---|---|
| NLL | 0.580 | 0.367 | −0.213 (50/50 better) |
| Brier | 0.327 | 0.202 | −0.125 (50/50 better) |
| ECE eq-freq | 0.0344 | 0.0342 | ≈0 (mixed 29/50) |
| ECE eq-width | 0.031 | 0.034 | ≈0 |

Calibration ECE is essentially tied; proper scoring rules strongly favor totalVI.

### S. AURC comparison
scVI **0.087**; totalVI **0.036**; Δ **−0.050** (totalVI better in **50/50**).

### T. Risk–coverage results (mean across 50 runs/model)
| Coverage | scVI acc / macro-F1 / risk | totalVI acc / macro-F1 / risk |
|---|---|---|
| 100% | 0.759 / 0.800 / 0.241 | 0.864 / 0.888 / 0.136 |
| 90% | 0.795 / 0.822 / 0.205 | 0.901 / 0.916 / 0.099 |
| 80% | 0.827 / 0.840 / 0.173 | 0.928 / 0.936 / 0.072 |
| 70% | 0.859 / 0.857 / 0.141 | 0.949 / 0.950 / 0.051 |
| 60% | 0.892 / 0.874 / 0.108 | 0.963 / 0.958 / 0.037 |
| 50% | 0.925 / 0.890 / 0.075 | 0.975 / 0.961 / 0.025 |

### U. Per-class retention imbalance
At 50% coverage, mean per-run retention range ≈ **0.80–0.83**; CV ≈ **0.57–0.59**.
Abstention remains **highly class-imbalanced** on this independent cohort.

### V. Classes most frequently rejected (lowest retention_50)
- scVI: **CD8T_Naive** (0.17), CD4T_Naive (0.22), CD4T_Mem (0.35)
- totalVI: **CD4T_Naive** (0.21), CD4T_Mem (0.26), CD8T_Naive (0.33)
Best retained: CD14_Mono ≈ 1.0; B ≈ 0.94–0.95.

### W. Between-donor-fold variability (descriptive SD of fold means)
Examples (macro-F1): scVI **0.013**; totalVI **0.011**.
NLL between-fold SD notably larger (scVI 0.049; totalVI 0.035).

### X. Within-fold model-seed variability
Examples (macro-F1): scVI **0.009**; totalVI **0.006**.
Exception: scVI **latent_error_auroc** within-fold seed SD **0.047** ≫ between-fold SD 0.015.

### Y. Which variability source appears larger?
For most primary metrics (**macro-F1, predictive AUROC, NLL, Brier, AURC**): **between-fold (donor) > within-fold seed**.
For **scVI latent AUROC**: **seed > fold**.
Descriptive only (not REML); do not over-interpret fractions.

### Z. External replication table (PHASE 6–12A → 12B)
See `tables/phase12b_replication_summary.csv` / Figure 10.

| Finding | Status |
|---|---|
| 1. totalVI biological predictive gain | **Replicated** |
| 2. Predictive uncertainty detects errors | **Replicated** |
| 3. totalVI latent uncertainty more informative | **Replicated** (direction) |
| 4. Uncertainty predicts fragility under source perturbation | **Not testable** here |
| 5. Selective prediction reduces risk | **Replicated** |
| 6. Abstention rejects difficult/rare classes | **Replicated** |
| 7. Run-level stochastic variation material | **Partially replicated** (seed SD small but non-zero; fold dominates) |

### AA. Major conclusions that replicated
1. totalVI > scVI on donor-held-out macro-F1 (large, sign-consistent).
2. U_pred is a strong error detector externally.
3. Confidence abstention improves risk/AURC but with severe class retention imbalance.
4. totalVI latent posterior uncertainty is more coupled to error than scVI’s (direction).

### AB. Failed to replicate / nulls
No primary PHASE 11B claim **failed** outright. ECE advantage for totalVI did **not** appear (tied).
Source-composition fragility (PHASE 11B finding 4) was **not tested**.

### AC. Only partially testable
- Source-composition perturbation fragility.
- Formal REML variance fractions (5 folds × 10 seeds is descriptive).
- Full 14-protein / Seurat-v4 label transfer (intentionally excluded for independence).
- Stimulated conditions (excluded by DESIGN A).

### AD. Compute / runtime / RAM
- Grid: 5 folds × 2 models × 10 seeds = **100**.
- Mean runtime/run: scVI ~**134 s**; totalVI ~**229 s** (sum ≈ 5.0 GPU-hours of per-run timers).
- Wall-clock (first→last done-marker, including gates/restarts): ~**5.2 h**.
- WSL RAM ~7.6 GB; long single-process runs leaked RSS → mitigated by **chunked orchestrator** (fresh Python every 8 jobs). Available RAM typically 5+ GB after restart.
- Formal-validation tree on server ≈ **12 GB** (models dominate).

### AE. Warnings and limitations
- Author labels are coarser than PBMC `celltype.l2`; cross-study label ontology is imperfect.
- 12/14 proteins only; CD15/TIGIT absent.
- Baseline-only; stimulation shift untested.
- Latent AUROC for scVI is seed-noisy (one seed near chance).
- Donor folds are consecutive pairs (deterministic), not random; fold effects entangle biology + technical lane structure.
- Class-aware retention reporting is a **baseline**, not a novel selective algorithm.
- Do not treat GSE164378 / Seurat v4 as Lawlor labels (independence constraint).

### AF. Files created (key)
```
results/phase12b_external_validation/formal_validation/
  processed/lawlor_baseline_sng_shared.h5ad
  done/*.json                         # 100
  tables/phase12b_primary_runs.csv
  tables/phase12b_model_deltas.csv
  tables/phase12b_variability.csv
  tables/phase12b_celltype_results.csv
  tables/phase12b_replication_summary.csv
  tables/phase12b_donor_folds.csv
  tables/phase12b_primary_protein_panel.csv
  figures/figure01_design.{png,pdf} … figure10_replication_summary.{png,pdf}
  logs/{data_prep_summary.json,setup_anndata.md,donor_fold_rule.md,analysis_summary.json,...}
  models/, embeddings/, posterior/, probabilities/
```

### AG. Exact resume / rerun command
```bash
ssh 4060-server
cd /home/jizhu/research/multiomics_robustness
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate multiomics_robustness
export PATH=/usr/lib/wsl/lib:$PATH LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH TMPDIR=$HOME/tmp
# resumes via done-markers; skips completed
python scripts/12b_orchestrate.py
# analysis only:
python scripts/12b_analyze.py
```

### AH. Is full 39-ADT sensitivity scientifically justified?
**Optional later sensitivity only** — not required for primary claim. Primary question was external validation under a panel close to the original 14-protein study. 39-ADT would change the modality signature and should be pre-registered separately.

### AI. Is stimulation-shift analysis scientifically justified?
**Scientifically interesting, not yet justified as automatic next step.** Would test condition shift (LPS / CD3_CD28), a different scientific question from donor generalization. Requires explicit approval and a frozen design.

### AJ. Should source-composition perturbation be repeated externally?
**Yes, if pursuing fragility claims in the manuscript** — finding 4 remains untested externally. Lower priority than consolidating Baseline donor-held-out evidence already in hand.

### AK. Sufficient external evidence for manuscript consolidation?
**Yes for the Baseline donor-held-out claims that PHASE 12B was designed to test** (totalVI transfer gain; predictive uncertainty; selective-prediction risk/retention imbalance; directional latent-uncertainty advantage).
Remaining gaps (39-ADT, stimulation, external source-composition) are **optional extensions**, not blockers for drafting the Baseline external-validation section.

---

## STOP

Completed: gated 100-run Baseline DESIGN A + analysis + Figures 01–10 + A–AK.

**Did not** run LPS / CD3_CD28 / 39-ADT / new selective algorithm / new model / commit / push.
Awaiting explicit approval for any next step.
