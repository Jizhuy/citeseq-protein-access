# FINAL SCIENTIFIC REVISION REPORT

**Date:** 2026-09-11  
**Primary manuscript:** `results/final_application_manuscript_v2/manuscript/Reliability_CITEseq_Application_Research_Manuscript.docx`  
**Editable source:** `.../Reliability_CITEseq_Application_Research_Manuscript.md`

---

## 1. Files modified

- `manuscript/Reliability_CITEseq_Application_Research_Manuscript.md` (full scientific rewrite)
- `manuscript/Reliability_CITEseq_Application_Research_Manuscript.docx` (regenerated)
- `figures/Figure1_study_design.{png,pdf}` (font readability)
- `figures/Figure2_information_vs_representation.{png,pdf}` (font readability)
- `figures/Figure3_transfer_and_uncertainty.{png,pdf}` (**renumbered**; was previous Figure 4)
- `figures/Figure4_protein_degradation.{png,pdf}` (**renumbered**; was previous Figure 3)
- `figures/Figure5_stochasticity_and_selective.{png,pdf}` (font readability)
- `figures/Figure6_lawlor_donor_validation.{png,pdf}` (font readability)
- `figures/Figure_Legends.docx`
- `tables/Main_Tables.docx`
- `supplement/Supplementary_Information.docx`
- `review/Final_Scientific_Audit.docx`
- `review/Application_Reader_Summary.docx`
- `scripts/01_generate_figures.py` (swap 3↔4 + larger fonts)
- `scripts/02_build_docx_package.py` (figure map + legend parsing)
- `logs/FINAL_SCIENTIFIC_REVISION_REPORT.md` (this file)

Historical manuscript under `results/final_manuscript_revision/` was **not** modified.

---

## 2. Exact classes of scientific wording changed

- Removed graduate-application / portfolio / admissions meta-language.
- Internal ΔF1 CI now stated as **cell-level conditional paired-bootstrap** (Results) / **conditional 95% CI** (Abstract).
- RNA-only audit reframed as **label identity / consistency**, not biological validation.
- Softened “consistently valuable” → restrained protein-information wording.
- Replaced combined-uncertainty superiority claim with **complementary diagnostics** wording.
- Direction B reframed as clearest **totalVI transfer advantage**, not “multimodal robustness.”
- Training-time non-monotonicity: “consistent with … stochasticity and optimization variability.”
- Removed “practically meaningful” threshold language (no prespecified threshold found).
- Removed defensive “No external repository link is required…” sentence.
- Softened “State 2” project jargon in Background to scientific STATE-2 framing without checklist tone.
- Related-work sentence added so MultiVI/Cobolt/scMaui/atlas benchmarking/WNN refs are cited without claiming they were run.

---

## 3. Methods sections merged/restructured

From ~34 numbered micro-sections → **9** scientific subsections:

1. Datasets and cell annotations  
2. Feature processing and representation methods  
3. Internal development-cohort evaluation  
4. Bidirectional transfer and scArches adaptation  
5. Protein degradation experiments  
6. Uncertainty and selective prediction  
7. Training stochasticity and variance decomposition  
8. Lawlor donor-level validation  
9. Software and reproducibility  

---

## 4. Reproducibility details added (and sources)

| Detail | Source |
|---|---|
| Entry-wise nonzero protein → zero corruption; preexisting zeros unchanged; fractions; 5 mask seeds | Historical Methods + SI (`results/final_manuscript_revision/...`); Part C script/logs |
| Training-time retrain totalVI / fixed scVI vs post-adaptation freeze + target-only corruption | Same + P1 Part C synthesis |
| Post-adaptation coordinate rule (same adapted model; clean source classifier) | P1 Part C design / script |
| PCA: RNA HVG2000/seurat_v3/norm/log1p/scale/PCA10; protein CLR+z-score/PCA10; concat 10+10 + StandardScaler; source-only inductive fit | `revision2_simple_multimodal_baseline.py` report block |
| MoM EMS variance decomposition; separate Δ cube; residual = design residual | `logs/variance_model_audit.md`; corrected variance CSVs |
| Target-cell paired bootstrap n=500 for HC ΔF1 CI | Historical manuscript Methods |
| scvi-tools 1.3.3; scikit-learn 1.5.2; Python 3.11.7 | `results/logs/environment.json` |
| Predictive U / latent mean posterior variance d=20; full AURC / pAURC / E-AURC | Frozen calibration tables + prior Methods |

---

## 5. Requested details NOT supported / not invented

- **No prespecified “practically meaningful” ΔF1 threshold** (e.g., −0.01) documented as an analysis rule → language removed rather than retrofitted.
- **MOFA+ full hyperparameter grid** beyond frozen historical embedding / comparable latent use → not expanded inventively.
- **Bootstrap CI type** beyond “target-cell paired bootstrap, n=500” (percentile vs BCa etc.) → not guessed.
- **Exact LogisticRegression C/penalty/multi_class** beyond documented sklearn defaults for installed version → left conservative.
- No new repository URL invented.

---

## 6. Figure renumbering performed

First-citation order corrected:

| New | Content | Previous |
|---|---|---|
| Figure 1 | Study design | Figure 1 |
| Figure 2 | Information vs representation | Figure 2 |
| Figure 3 | Transfer + uncertainty | **was Figure 4** |
| Figure 4 | Protein degradation | **was Figure 3** |
| Figure 5 | Stochasticity + selective | Figure 5 |
| Figure 6 | Lawlor donor validation | Figure 6 |

All inline callouts, legends, and DOCX embeds updated.

---

## 7. Figure readability changes

- Base font size increased (≈9→10; titles/axes/ticks enlarged).
- Regenerated from **frozen tables only** (no new experiments).
- No panel content removed; no decorative redesign.

---

## 8. References added/removed/renumbered

- **No references removed.**
- Previously uncited [5,6,7,9,14] now cited once via a Background related-work sentence as **related methods / evaluation literature, not experimental comparators**.
- Citation numbering unchanged (1–16).

---

## 9–12. Confirmations

9. **No experiment was run.**  
10. **No numerical result changed** without existing evidence (canonical frozen values retained).  
11. **NO COMMIT.**  
12. **NO PUSH.**

---

## Final quality checklist (Phase 20)

- [x] No new experiment  
- [x] No unsupported number change  
- [x] No graduate-application / portfolio meta language in scientific text  
- [x] Internal CI = conditional cell-level bootstrap  
- [x] RNA-only audit = label consistency, not biological validation  
- [x] Figure numbering = citation order  
- [x] No `Figure X.**` artifacts  
- [x] No editorial “continued…” caption notes  
- [x] Figure text enlarged for A4 readability  
- [x] Methods = 9 scientific subsections  
- [x] Protein corruption mechanism described from sources  
- [x] Variance method = MoM EMS from sources  
- [x] “Practically meaningful” removed (unsupported)  
- [x] Predictive uncertainty terminology correct  
- [x] Latent = posterior variance/dispersion (not epistemic)  
- [x] full AURC / pAURC / E-AURC distinct  
- [x] Lawlor donor vs fold×seed units separated  
- [x] 19/20 only PBMC Direction B  
- [x] Concat baseline prominent  
- [x] scArches = unsupervised transductive adaptation  
- [x] Unused refs either cited appropriately or (none left unused)  
- [x] No unsupported benchmarking claim  
- [x] Repository defensive sentence removed  
- [x] Discussion tightened  
- [x] Abstract restrained (221 words)  
- [x] DOCX/MD agree  
- [x] Historical manuscript untouched  
- [x] No commit / no push  

**STOP.**
