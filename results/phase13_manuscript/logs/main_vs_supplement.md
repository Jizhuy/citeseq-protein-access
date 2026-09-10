# Main vs supplementary vs omit

## MAIN TEXT (story spine)
1. Representation quality across PCA, MOFA+, scVI_matched, totalVI (metric plurality).
2. Entry-wise protein corruption: multimodal advantage remains directionally positive.
3. Cross-dataset transfer is direction/context dependent (Dir A ≈ null; Dir B positive).
4. Predictive uncertainty detects errors under shift.
5. Latent posterior uncertainty is model-dependent (stronger for totalVI; externally attenuated but directionally replicated).
6. Uncertainty tracks stochastic fragility; run-level residual dominates variance once separated.
7. Selective prediction lowers risk but creates cell-type retention imbalance.
8. Independent Lawlor donor-held-out validation reproduces central reliability findings.

## SUPPLEMENTARY
- Default vs matched scVI geometry
- Full MOFA diagnostics / variance explained
- l1 metrics; alternate kNN
- Per-protein recovery under corruption
- Latent drift / coordinate-alignment diagnostics (PHASE 10 erratum)
- Source-composition diagnostics (PHASE 10B detail)
- 20-seed convergence / ensemble calibration (incl. ECE worsening)
- Monte Carlo confirmation of posterior variance semantics
- Complete 5×5 variance heatmaps + power analysis
- Per-class retention curves (all coverages)
- Gene/protein harmonization for Lawlor
- Per-fold Lawlor metric tables
- ECE binning sensitivity
- Runtime / RAM diagnostics
- PHASE 9 API limitation narrative

## OMIT FROM MANUSCRIPT (retain in repository)
- Chronological “phase diary” narrative
- Obsolete mixed-coordinate PHASE 10 numbers
- PHASE 11 5-seed primary claims (superseded)
- Any proposal of precision-weighted PoE as a new model
- Unrun LPS / CD3_CD28 / 39-ADT sensitivity as results
- Fabricated missing-modality experiments
- Exhaustive metric dumps that do not change the argument
