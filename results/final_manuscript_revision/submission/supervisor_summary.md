# Supervisor Summary

## TITLE
Reliability of scVI and totalVI for CITE-seq under distribution shift and uncertainty-guided prediction

## RESEARCH QUESTION
How reliable are probabilistic deep generative representations (scVI_matched vs totalVI) for CITE-seq when protein measurements are degraded, datasets or donors shift, stochastic training varies, and uncertainty is used for selective prediction?

## WHY IT MATTERS
Average accuracy is insufficient for trustworthy multi-omics analysis. Reliability under shift, stochasticity, and abstention-induced biological coverage loss is what users need in practice.

## DATASETS
- PBMC10k / PBMC5k (development)
- Lawlor Baseline CITE-seq (independent external; author labels; HCA/ENA accessions)

## MODELS
PCA, MOFA+, scVI_matched, totalVI (scvi-tools 1.3.3); scArches adaptation with post-adaptation coordinate rule.

## EXPERIMENTAL PROGRAMME
Representation → entry-wise protein corruption → reciprocal transfer → predictive/latent uncertainty → replicated variance → selective prediction/retention → Lawlor donor-held-out validation.

## STRONGEST FINDINGS
Lawlor macro-F1 0.800 → 0.888 (Δ +0.088; 50/50). Predictive AUROC internal ≈0.84–0.86; Lawlor 0.801 vs 0.838. Latent AUROC internal totalVI ≈0.61–0.67 vs scVI ≈0.47–0.54; Lawlor 0.650 vs 0.605. AURC Lawlor 0.036 vs 0.087.

## EXTERNAL VALIDATION
Independent Lawlor Baseline; author labels only; 5 folds × 10 seeds = 100 runs; hierarchical fold vs seed variability reported.

## NEGATIVE FINDINGS
Direction A near-null transfer; ECE not decisive; selective prediction creates cell-type retention imbalance; residual run-level variability dominates variance decomposition.

## LIMITATIONS
PBMC-only; partial Lawlor labels; 12/14 protein overlap; no stimulation/39-ADT/external source-composition in primary package; no honest cell-wise missing-modality experiment under totalVI 1.3.3.

## TARGET JOURNAL
Genome Biology — Research (backup: Bioinformatics Original Paper)

## OPTIONAL REVIEWER-CONTINGENCY EXPERIMENTS (NOT CURRENTLY REQUIRED)
A. Baseline→LPS/CD3_CD28 shift. B. 39-ADT sensitivity. C. External source-composition perturbation.

## QUESTIONS FOR SUPERVISOR
1. Is reliability framing strong enough for Genome Biology?
2. Add stimulation shift proactively or reserve for reviewers?
3. Are six main figures appropriate?
4. Keep variance decomposition in main text or move partly to Supplement?
5. Is Lawlor prominence sufficient?
6. Reliability study vs benchmarking study framing?
