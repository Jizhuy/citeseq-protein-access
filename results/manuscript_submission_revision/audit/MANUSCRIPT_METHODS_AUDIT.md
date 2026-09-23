# Manuscript Methods Audit

**Date:** 2026-09-23  
**Canonical manuscript:** “Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift”  
**Canonical source (preferred):** `../CITEseq_Protein_Access_Model_Associated_Gains_Presubmission.md` (parent of repo) and synced copies under `results/final_information_vs_model_manuscript/`  
**Older manuscript (methods reservoir only):** `Reliability_CITEseq_Manuscript_v2.pdf` / `results/final_application_manuscript_v2/manuscript/Reliability_CITEseq_Application_Research_Manuscript.md`  
**Rule:** No value accepted from older text unless verified in code, config, logs, or frozen CSVs.

---

## Status legend

| Status | Meaning |
| --- | --- |
| **verified** | Confirmed in code and/or frozen result files |
| **manuscript incomplete** | Verified in repo but missing/under-specified in canonical Methods |
| **inconsistent** | Manuscript wording conflicts with code or frozen files |
| **cannot verify** | Not recoverable from available artifacts; do not invent |

---

## A. Primary numerical results (frozen)

| Item | Manuscript | Repository source | Status |
| --- | --- | --- | --- |
| RNA PCA F1 0.666 | Yes | `simple_multimodal_development.csv` 0.666495 | **verified** |
| Protein PCA 0.512 | Yes | same 0.512031 | **verified** |
| Concat 0.745 | Yes | same 0.744897 | **verified** |
| MOFA+ 0.549 | Yes | same 0.549167 | **verified** |
| scVI 0.722 | Yes | scVI_matched 0.721970 | **verified** (matched setting, not scVI default) |
| totalVI 0.779 | Yes | 0.778906 | **verified** |
| G_access +0.079 [0.043,0.132] | Yes | bootstrap summary | **verified** |
| G_assoc +0.034 [−0.012,0.073] | Yes | same | **verified** |
| Gap_total +0.113 [0.084,0.154] | Yes | same | **verified** |
| R_access ~0.70 [0.405,1.114] | Yes | point 0.697456 | **verified** (rounding) |
| Transfer inductive 0.735/0.734 | Yes | matched_transfer_control.csv | **verified** |
| Transductive 0.636/0.767 | Yes | same | **verified** |
| scVI/totalVI transfer means | Yes | phase11b_20seed_summary | **verified** |
| Dir B Δ +0.052; 19/20 | Yes | phase11b_seed_delta_summary | **verified** |
| Training-time degradation curve | Yes | degradation_audit CSVs | **verified** |
| Dir B post-adapt 0.654→0.621 | Yes | test_time_corruption_summary.csv | **verified** |
| Lawlor Δ +0.090; 10/10; SD≈0.010 | Yes | lawlor_donor_delta SUMMARY 0.090198; SD 0.009562 | **verified** |
| Lawlor protein/concat ~0.927/~0.945 | Yes | simple_multimodal_lawlor.csv fold means 0.9267 / 0.9451 | **verified** |
| Lawlor scVI/totalVI ~0.800/~0.888 | Yes | donor means **0.7936 / 0.8838** | **inconsistent** (loose ≈; prefer 0.794/0.884) |
| CD8 Naive ΔF1 ≈ −0.26 at p=0.85 Dir B | Older MS | per_class CSV mean Δ −0.2642 | **verified** (restore to Supplement) |
| 38/50 latent AUROC consistency | Yes | Lawlor pairwise consistency tables | **verified** (computational, not biological) |
| Cells 10849 / 9494 / 6855 / 3994; genes 15792; proteins 14 | Yes | audits + development CSV | **verified** |

---

## B. Dataset provenance and QC

| Methodological item | Current MS | Older MS | Repository | Status |
| --- | --- | --- | --- | --- |
| PBMC10k/5k source | “public 10x”; GEO not verified | scvi-tools example files | `config/config.yaml` URLs: scverse exampledata + YosefLab/scVI-data fallback; filenames `pbmc_*_protein_v3.h5ad` | **verified** provenance via scvi-tools files; **no GEO/SRA** in repo |
| Combined dimensions | 10849 × 15792 × 14 | same | data/README + audits | **verified** |
| Additional cell QC (min/max genes, mito, doublets) | Incomplete | Sometimes implied | `preprocessing.additional_cell_filter: false`; comment: source already filtered per totalVI reproducibility scripts | **verified: no additional QC applied** |
| Gene filters / HVG for VAE training | Incomplete | 4000 HVG in places | Config `n_top_genes: 4000` **not applied** to scVI/totalVI runners (full shared gene universe) | **manuscript incomplete** + potential confusion with PCA HVG=2000 |
| Protein panel (14 shared) | Incomplete list | lists often present | `config/protein_gene_map.yaml` (14 symbols); raw TotalSeq names in 10k file; CD45 suffix harmonization needed for raw join | **partially verified**; exact 14 column strings for processed inner: recover from h5ad when available |
| Lawlor accessions | HCA + ENA present | same | phase12b provenance logs | **verified** |
| Lawlor labeled/singlet counts | Incomplete | 16175 singlets; 5207 labeled | Older MS + phase12b; confirm from logs before asserting | **manuscript incomplete**; values in older MS need re-check vs phase12b tables |
| Healthy-control / high-confidence | 9494 HC | same | `annotation_tier_l2=="high"`; threshold 0.85 in confidence_thresholds.json | **verified** |

---

## C. PCA / concat preprocessing

| Item | Current MS | Older MS | Code | Status |
| --- | --- | --- | --- | --- |
| 2000 HVG, seurat_v3 | Incomplete | Yes | `revision2_p1_common.py` N_HVG=2000 | **verified**; **config 4000 ≠ PCA path** |
| normalize 1e4, log1p, scale clip ±10 | Incomplete | Yes | same | **verified** |
| 10 RNA + 10 protein PCs | Incomplete | Yes | PRIMARY_*_PCS=10 | **verified** |
| Protein CLR +1 pseudocount, ddof=1 z-score | Missing | Partial | `clr_rows` + std(ddof=1) | **manuscript incomplete** |
| StandardScaler on train/source only (inductive) | Partial | Yes | verified | **verified** |
| Transductive joint fit | Partial | Case C language | revision3 script | **verified**; rename “Case C” → feature-access-matched |
| Fit cells for development | Incomplete | 8679/2170 then HC | Split: test_fraction 0.20, stratify batch, seed 0 on **all** 10849; HC eval 7573/1921 | **verified**; older “8679/2170” is all-cell split implication (0.8×10849=8679.2) — **verify from obs['split'] if needed** |

---

## D. scVI / totalVI / scArches

| Item | Current MS | Older MS | Code/logs | Status |
| --- | --- | --- | --- | --- |
| Primary RNA deep model is **scVI_matched** | Calls “scVI” | Lists matched hyperparams | n_latent=20, n_hidden=256, n_layers=2, dropout=0.2, gene_likelihood=**nb** | **verified**; distinguish from scVI **default** (n_hidden=128, n_layers=1, zinb, dropout 0.1) |
| totalVI architecture | Incomplete | Yes | n_latent=20, n_hidden=256, enc 2 / dec 1, dropout 0.2, gene_likelihood=nb, protein_dispersion=protein | **verified** |
| max_epochs 200, early_stopping, patience 45, batch 256, train_size 0.9 | Missing | Partial | code | **manuscript incomplete** |
| totalVI lr=4e-3, reduce_lr_on_plateau | Missing | Partial | code | **manuscript incomplete** |
| scVI lr | Missing | often 1e-3 | not passed → **library default** | **cannot verify numeric from project code**; state “scvi-tools TrainingPlan default” |
| Batch key | Partial | batch / run_identifier | Internal: `batch`; Lawlor: `run_identifier` | **verified** |
| scArches sequence | Partial | Yes | source train → prepare_query → adapt → encode **both** with adapted model → logreg on post-adapt source | **verified** from phase10 |
| Actual early-stop epoch (baseline) | Missing | sometimes claimed | baseline manifests **absent** in checkout | **cannot verify** for PHASE 3–5; PHASE 10 logs show 200/200 epochs |

---

## E. Classifier

| Item | MS | Code | Status |
| --- | --- | --- | --- |
| LogisticRegression max_iter=2000, solver=lbfgs, random_state=seed | Incomplete | verified | **verified** |
| C=1.0, class_weight=None, multi_class default | Missing | sklearn defaults | **manuscript incomplete** (state defaults explicitly) |

---

## F. Label provenance

| Item | MS | Code | Status |
| --- | --- | --- | --- |
| RNA-only annotation audit, agreement 1.0 on 9494 | Partial (info MS) / incomplete (presubmission) | revision2_rna_only_annotation + tables | **verified**; restore dedicated subsection |
| Lawlor protein-gated labels | Yes | author labels | **verified** |

---

## G. Bootstrap

| Item | MS | Code | Status |
| --- | --- | --- | --- |
| n=10,000 cell-level paired | Yes | BOOT_SEED=20260923; percentile 2.5/97.5 | **verified** |
| Not BCa; not stratified | Incomplete | global with-replacement over test HC cells | **manuscript incomplete** |
| Class-stratified sensitivity | Absent | not yet run | **manuscript incomplete** → new analysis |

---

## H. Degradation

| Item | MS | Code/config | Status |
| --- | --- | --- | --- |
| mask_nonzero_protein_to_zero; levels; 5 seeds; model_seed=0 | Partial | config sparsity_stress | **verified** |
| Class-specific Dir B results | Missing in canonical | test_time_corruption_per_class.csv | **manuscript incomplete** → Supplement |
| Marker-channel dropout / alternative corruption | Absent | not implemented | **blocked** for totalVI retrain; concat-level sensitivity feasible |

---

## I. Uncertainty / selective prediction / stochasticity

| Item | MS | Code | Status |
| --- | --- | --- | --- |
| U_pred, U_latent definitions | Partial | phase11/11b | **verified** conceptually |
| ECE bins, NLL, Brier, AURC formulas | Incomplete | 07_09 / phase11b scripts | **manuscript incomplete** |
| Selective risk = 1−accuracy | Incomplete | phase11b_selective_prediction | verify formula in code |
| Macro-F1 vs coverage sensitivity | Absent | not run | **new analysis** |
| 5×5×2 design + MoM variance | Compressed | phase12a | **manuscript incomplete** |
| REML sensitivity | Older MS may claim | phase12a_stats notes | verify before citing |

---

## J. Lawlor design

| Item | MS | Code | Status |
| --- | --- | --- | --- |
| 10 donors, 5 folds, 2 targets/fold, 10 seeds | Incomplete | phase12b + 03_lawlor_donor_level | **verified** pattern (confirm fold table) |
| present-class macro-F1 | Yes | code | **verified** |
| Donor-level summary stats beyond mean/SD | Incomplete | recoverable from lawlor_donor_delta | **manuscript incomplete** |
| Exact sign test | Absent | can compute from 10/10 | **new optional** |

---

## K. Comparators / second cohort

| Item | Status |
| --- | --- |
| Additional multimodal method (WNN/MultiVI/Cobolt) | **blocked**: no implementation in repo; MOFA+ already secondary |
| Second external non–protein-gated cohort | **blocked**: not present in repo; write FUTURE TODO |

---

## L. Environment

| Package | Pin | Status |
| --- | --- | --- |
| Python 3.11.7 | environment.yml | **verified** |
| scvi-tools 1.3.3 | env + phase10 log | **verified** |
| scikit-learn 1.5.2 | environment.yml | **verified** (host may differ) |
| torch 2.14.0 | env; CUDA build in phase10 log | **verified** |
| scanpy 1.10.3, anndata 0.11.4 | env | **verified** |

---

## Critical inconsistencies to fix in revision

1. Prefer donor means **0.794 / 0.884** (or state “approximately 0.79 / 0.88”) instead of 0.800 / 0.888 if three-decimal precision is used.  
2. Always say **scVI_matched** hyperparameters when citing 0.722 / transfer means; do not conflate with scVI default (zinb, 128×1).  
3. PCA HVG = **2000**; config 4000 is unused by the primary PCA/VAE pipelines as implemented.  
4. Rename “residual model-associated gain” → **“residual totalVI-associated contrast”** (`G_assoc`) with non-causal definition.  
5. Restore Methods completeness from verified items only; flag cannot-verify items explicitly.
