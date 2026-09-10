# Supervisor summary (1–2 pages)

**Status:** For scientific review only. Do not submit. Do not release Zenodo/public repo without approval.

---

## TITLE

Reliability of deep generative single-cell multi-omics models under distribution shift and uncertainty-guided prediction

*(Alternative shortlisted: “Reliability and uncertainty of deep generative models for single-cell multi-omics integration” — less specific about shift/selective prediction.)*

## RESEARCH QUESTION

How reliable are widely used deep generative CITE-seq models (especially architecture-matched scVI vs totalVI) under protein degradation, dataset/donor shift, training stochasticity, and uncertainty-guided prediction—without inventing a new fusion architecture?

## WHY IT MATTERS

scVI/totalVI are already adopted for multi-omics integration. Average accuracy is well studied; reliability under realistic stresses (degraded proteins, shift, seed fragility, abstention-induced biological coverage loss) is less systematically characterized. Journals and users need effect sizes, variability, and external replication—not only leaderboard wins.

## DATASETS

| Cohort | Role |
|---|---|
| PBMC10k / PBMC5k CITE-seq (10x) | Development: representation, entry-wise protein corruption, reciprocal transfer, uncertainty, variance, selective prediction |
| Lawlor Baseline CITE-seq (HCA `efea6426-…`; ENA `PRJEB40376`) | Independent external validation; **author labels only**; 12 matched proteins; 5 donor-held-out folds × 10 seeds |

Development labels used a Seurat-v4 multimodal reference pathway for PBMC annotation. Those Hao labels were **not** applied to Lawlor.

## MODELS

PCA (RNA), MOFA+ (0.7.5), scVI_matched, totalVI (scvi-tools 1.3.3). Primary contrasts are architecture-matched scVI vs totalVI. scArches query adaptation with a post-adaptation coordinate rule (no mixed coordinates).

## CORE FINDINGS

1. Multimodal representations can improve fine-grained biological prediction, but do not dominate every metric.
2. Entry-wise protein corruption: average totalVI−scVI_matched macro-F1 gain remains positive even at high corruption fractions (not a missing-modality experiment).
3. Dataset shift: multimodal benefit is **direction-dependent** (near-null Dir A; Dir B Δ ≈ +0.052, 19/20 seeds).
4. Predictive uncertainty is strongly error-informative (internal AUROC ≈0.84–0.86).
5. Latent posterior variance tracks error more consistently for totalVI than matched scVI (internal ~0.61–0.67 vs ~0.47–0.54); attenuated but same direction externally.
6. Stochastic/residual run variability dominates a replicated variance decomposition (not “interaction dominates”).
7. Confidence abstention lowers risk/AURC but creates severe cell-type retention imbalance.

## EXTERNAL VALIDATION (strongest quantitative block)

Lawlor Baseline, author labels, 100 runs:

| Metric | scVI | totalVI |
|---|---:|---:|
| macro-F1 | 0.800 | 0.888 (Δ +0.088; 50/50) |
| Predictive AUROC | 0.801 | 0.838 |
| Latent AUROC | 0.605 | 0.650 |
| AURC | 0.087 | 0.036 |

## MAIN NEGATIVE / NULL-LEANING FINDINGS

- Multimodal transfer benefit is not universal across reciprocal shift directions.
- ECE is not a decisive totalVI advantage.
- Latent uncertainty is **not** a causal claim about protein supervision “calibrating” posteriors.
- Selective prediction is a diagnostic baseline, not a novel algorithm.
- No honest cell-wise missing-protein-modality experiment under totalVI 1.3.3.

## LIMITATIONS

Both cohorts PBMC; Lawlor labels coarse/partial; 12/14 protein overlap; stimulation shift and 39-ADT sensitivity not run; source-composition not repeated externally; limited variance-design levels; predictive uncertainty is downstream of latents.

## WHY NO NEW ARCHITECTURE IS CLAIMED

Precision-weighted / Product-of-Experts multimodal VAE ideas already exist (e.g., Wu & Goodman; MultiVI and related). Contribution = reliability characterization + external replication.

## TARGET JOURNAL

Primary: **Genome Biology — Research**  
Backup: **Bioinformatics — Original Paper**

## OPTIONAL EXPERIMENTS IF REVIEWERS REQUEST THEM

Marked **NOT CURRENTLY REQUIRED** (see `logs/reviewer_contingency_experiments.md`):

- A. Baseline → LPS / CD3_CD28 stimulation shift  
- B. 39-ADT sensitivity  
- C. External source-composition perturbation  

---

## Questions for supervisor feedback (do not answer for the PI)

1. Is the reliability-focused framing strong enough for Genome Biology?
2. Should stimulation-shift validation be added proactively or reserved for reviewer response?
3. Are six main figures appropriate?
4. Should stochastic-variance decomposition remain in the main paper or move partially to Supplement?
5. Is the Lawlor external validation sufficiently prominent?
6. Is the paper better framed as a reliability study or a benchmarking study?
7. Preferred author order / corresponding author / funding text?
8. Approve public code/data release plan (private → Zenodo) before journal submission?
