# Final Scientific Audit

**Manuscript:** Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift  
**Revision package:** `results/manuscript_submission_revision/`  
**Date:** 2026-09-23

## Checklist (strict reviewer mode)

| ID | Question | Verdict |
| --- | --- | --- |
| A | Can the internal benchmark be reproduced? | **Yes** from frozen CSVs + scripts; VAEs not re-trained in this pass |
| B | Feature leakage? | **Controlled** for inductive/dev paths; transductive concat **intentionally** uses unlabeled target features (documented) |
| C | Labels independent of evaluated proteins internally? | **Yes** for development (RNA-only audit agreement 1.0); Lawlor labels are protein-gated (explicit) |
| D | Inductive vs transductive clear? | **Yes** after Methods rewrite |
| E | Transductive concat fairly described? | **Yes** — feature-access matched, not algorithmic scArches replica |
| F | Bootstrap intervals correctly interpreted? | **Yes** — cell-level conditional; stratified sensitivity agrees qualitatively |
| G | Biological vs computational replication? | **Yes** — donors *n*=10 vs seeds/folds (19/20, 38/50) separated |
| H | \(G_{\mathrm{assoc}}\) non-causal? | **Yes** — renamed residual totalVI-associated contrast |
| I | Lawlor limited by protein-gated labels? | **Yes**; protein/concat > totalVI retained |
| J | Corruption interpreted cautiously? | **Yes**; rebound within variability; marker dropout Supplement |
| K | Selective prediction composition? | **Yes**; CV and macro-F1 sensitivity documented |
| L | Hyperparameters documented? | **Largely yes** in Methods + Suppl. S3; scVI LR = library default (flagged) |
| M | Software versions documented? | **Yes** (pinned env + `print_environment.py`) |
| N | Figure/table → code/raw? | **Yes** via `MANUSCRIPT_OUTPUT_MANIFEST.md` |
| O | Manual number errors? | Lawlor deep means corrected **0.794/0.884** (was ≈0.800/0.888) |
| P | Unsupported claims? | Avoided universality / causal architecture language |
| Q | Repo contradictions to narrative? | Marker dropout shows concat can be buffered by RNA; does **not** overturn protein-access story |
| R | Framing remains protein-access vs residual totalVI? | **Yes** — not reverted to broad reliability benchmark |

## Conclusions changed?

**No primary conclusion change.**  
Sensitivities reinforce caution:

- Stratified bootstrap: \(G_{\mathrm{assoc}}\) CI still crosses zero.  
- Class degradation: CD8 Naive ΔF1 ≈ −0.264 under Dir B post-adapt corruption.  
- Marker dropout: protein PCA fragile; concat more RNA-buffered.  
- Lawlor sign test *P*=0.00195 on *n*=10 donors (supportive, not overclaimed).

## Remaining limitations / blocked tasks

1. totalVI retrain under marker-channel dropout — **blocked**  
2. Additional multimodal comparator (WNN/MultiVI) — **blocked**  
3. Second external non–protein-gated cohort — **blocked**  
4. Baseline early-stop epoch manifests missing — **cannot verify**  
5. DOCX regeneration / Git push — **not done** (await explicit instruction)

## Deliverables produced

- `audit/MANUSCRIPT_METHODS_AUDIT.md`  
- `audit/MANUSCRIPT_REVISION_PLAN.md`  
- `manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md` (expanded Methods)  
- `supplement/SUPPLEMENTARY_METHODS.md`  
- `supplement/SUPPLEMENTARY_RESULTS.md`  
- `MANUSCRIPT_OUTPUT_MANIFEST.md`  
- `REPRODUCIBILITY_README.md`  
- New scripts + tables/figures under this folder  
- `scripts/print_environment.py` at repo root
