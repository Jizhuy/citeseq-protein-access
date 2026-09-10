# Main figure plan (maximum 6)

## Figure 1 — Study design and reliability overview
**Question:** What reliability stresses were evaluated?
- Panels: cohorts (PBMC10k/5k; Lawlor); models; entry-wise corruption; transfer directions; uncertainty endpoints; donor-held-out external design.
- Source: design schematics compiled from PHASE 7–12B manifests.

## Figure 2 — Representation quality and protein-corruption robustness
**Question:** Does multimodal modeling improve biology, and does that survive protein degradation?
- A: model comparison schematic
- B: l2 high-conf logreg macro-F1 (PCA / MOFA+ / scVI_matched / totalVI)
- C: Delta macro-F1 vs entry-wise corruption fraction
- D: optional class-specific robustness
- Source: `mofa_biological_representation_metrics.csv`; `protein_sparsity_primary_delta_summary.csv`

## Figure 3 — Cross-dataset transfer and source-context sensitivity
**Question:** Are multimodal transfer gains direction-dependent?
- A: Dir A/B design
- B: paired macro-F1 (prefer 20-seed; show corrected 5-seed as consistency check)
- C: source-size sensitivity
- D: class-specific transfer effects (supplement-ready detail if crowded)
- Source: `tables_corrected/*`; PHASE 10B; PHASE 11B macro-F1

## Figure 4 — Predictive and latent uncertainty
**Question:** What does “uncertainty quality” mean, and is latent variance informative?
- A: predictive error AUROC
- B: latent error AUROC
- C: predictive vs latent relationship
- D: NLL / Brier / AURC (not ECE-as-discrimination)
- Source: PHASE 11B 20-seed tables

## Figure 5 — Stochastic fragility and selective prediction
**Question:** Does abstention help, and at what biological cost?
- A: uncertainty vs instability
- B: replicated variance decomposition (source / seed / interaction / residual)
- C: risk–coverage
- D: cell-type retention imbalance
- Source: PHASE 11B + 12A

## Figure 6 — Independent Lawlor donor-held-out validation
**Question:** Do key reliability findings replicate externally?
- A: cohort / 5 donor folds
- B: macro-F1
- C: predictive AUROC
- D: latent AUROC
- E: risk–coverage
- F: class retention
- Source: PHASE 12B formal_validation tables/figures
