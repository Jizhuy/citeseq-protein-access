# Manuscript Revision Plan

**Canonical paper:** protein-access vs residual totalVI-associated contrasts  
**Do not** revert to a broad “reliability benchmark” narrative.

---

## 1. Mandatory manuscript corrections

| Priority | Change |
| --- | ---: |
| P0 | Expand Methods using **verified** parameters only (PCA, scVI_matched vs default, totalVI, classifier, splits, leakage rules, Lawlor design, bootstrap definition). |
| P0 | Rename residual contrast language: **residual totalVI-associated contrast** (`G_assoc`); state it is not a causal architecture effect. |
| P0 | Clarify bootstrap = cell-level **conditional** percentile intervals; not biological CIs. |
| P0 | Fix Lawlor deep averages to **0.794 / 0.884** (file means) or explicitly “approx.” with caveat; keep Δ=+0.090, SD≈0.010, 10/10. |
| P0 | State internal “scVI” results are **scVI_matched** (nb, 256, 2 layers, dropout 0.2). |
| P0 | Document PBMC provenance via scvi-tools/10x example files + URLs; **no invented GEO/SRA**. |
| P0 | Explicit: **no additional** cell/gene QC beyond inherited processed objects (`additional_cell_filter: false`). |
| P0 | Dedicated **Annotation provenance and label-independence audit** subsection. |
| P0 | Transfer Methods as algorithms: inductive vs feature-access-matched transductive; scArches as unsupervised transduction with post-adaptation coordinate system. |
| P1 | Restore stochasticity 5×5×2 + MoM model equation. |
| P1 | Restore Lawlor fold construction, batch=`run_identifier`, 12 matched proteins (verify count), gene mapping notes if verified. |
| P1 | R_access joint-bootstrap note (not constrained to [0,1]; not causal fraction). |
| P1 | Language audit: remove “Case C”, “architecture superiority”, causal wording. |

---

## 2. Supplementary additions

- Class-stratified paired-bootstrap sensitivity (vs primary global bootstrap).  
- Cell-type–specific post-adaptation degradation (Dir B; CD8 Naive Δ≈−0.26).  
- Selective-prediction **macro-F1 / balanced error** vs coverage.  
- Full hyperparameter tables (S2–S3).  
- Software version table (S10).  
- Lawlor donor-level full table + fixed-ontology sensitivity if already computed.  
- Full calibration metrics tables from existing CSVs.  
- Seed-level transfer table (already in phase11b).  

---

## 3. New analyses (this revision)

| Analysis | Feasibility | Location |
| --- | --- | --- |
| Class-stratified bootstrap | **Yes** (no retrain) | `scripts/run_class_stratified_bootstrap.py` |
| Class-level degradation summary + figure | **Yes** (from CSV) | `scripts/summarize_class_degradation.py` |
| Lawlor donor summary + exact sign test | **Yes** | `scripts/lawlor_donor_summary_stats.py` |
| Selective macro-F1 vs coverage | **Yes if** prediction arrays available; else from phase11b selective tables where possible | check phase11b |
| Marker-channel dropout on **concat/protein PCA** | **Yes** (development HC) | `scripts/marker_channel_dropout_pca.py` |
| Marker-channel dropout + **totalVI retrain** | **Blocked** this pass (GPU/time); TODO |
| Additional multimodal comparator (WNN/MultiVI) | **Blocked** (no code); TODO |
| Second external cohort | **Blocked** (no data in repo); FUTURE TODO |

---

## 4. Blocked / unverifiable (do not invent)

- Baseline PHASE 3–5 training manifests / actual early-stop epochs.  
- Exact scVI TrainingPlan learning rate numeric value (library default).  
- GEO/SRA for PBMC10k/5k.  
- WNN/MultiVI/Cobolt results.  
- Independent non–protein-gated external cohort.  
- Full totalVI marker-dropout retraining curve.  

---

## 5. Deliverable map

| Output | Path |
| --- | --- |
| Methods audit | `audit/MANUSCRIPT_METHODS_AUDIT.md` |
| This plan | `audit/MANUSCRIPT_REVISION_PLAN.md` |
| Revised canonical MD | `manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md` |
| Supplementary Methods/Results | `supplement/` |
| Manifest | `MANUSCRIPT_OUTPUT_MANIFEST.md` |
| Final scientific audit | `audit/FINAL_SCIENTIFIC_AUDIT.md` |
| New scripts/tables/figures | `scripts/`, `tables/`, `figures/` |

---

## 6. Execution order (remaining)

1. ~~Repository audit~~  
2. ~~Verify primary numbers~~  
3. Run sensitivity scripts (bootstrap stratified, class degradation, Lawlor stats, marker dropout PCA, env printer).  
4. Expand Methods + Supplement.  
5. Update figures captions / consistency.  
6. Final scientific audit.  
7. **STOP** (no push unless instructed).
