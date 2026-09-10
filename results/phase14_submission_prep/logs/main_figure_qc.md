# Main figure QC

| Figure | Panels | Source artifacts (unchanged scientifically) | Issues / notes |
|---|---|---|---|
| Figure1 | Schematic A–style overview | Newly composed schematic (no numeric axes); lists cohorts/models/experiments | Does not imply all models in all experiments (caption note) |
| Figure2 | A UMAP collage; B l2 F1; C corruption F1; D ΔF1; E recovery | phase7 F1; phase8 figs 2/3/5; PCA/MOFA/scVI/totalVI UMAPs | A is geometry illustration, not a quantitative ranking |
| Figure3 | A design; B crossed macro-F1 heatmaps; C paired deltas; D source-size F1 | phase10 design; phase11b fig07/03; phase10b fig2 | Uses corrected-coordinate era analyses; 20-seed emphasis via 11B panels |
| Figure4 | A AUROC dist; B latent correct/incorrect; C pred vs latent; D paired deltas | phase11b figs 02/04/05/03 | ECE moved to supplement intent; D mixes metrics (legend must say so) |
| Figure5 | A unc vs subset instability; B variance by model; C risk-coverage; D retention | phase11b 11/13/14; phase12a fig02 | B is replicated variance (12A), not lumped 11B residual |
| Figure6 | A–F Lawlor design, F1, pred AUROC, lat AUROC, risk, retention | phase12b formal_validation figures 01–04,07,08 | Author labels only; 12-protein panel |

## Checklist
- [x] Panel letters present
- [x] Composite PNG+PDF written
- [x] No PHASE numbers in figure files themselves (source filenames may contain phase; manuscript legends do not)
- [x] No data alteration / point removal
- [ ] Human visual check at journal column width (remaining manual)
- [ ] Possible downsampling for upload size limits (remaining manual)
- [ ] Confirm Figure3B is clearly described as crossed macro-F1 (not AUROC)

## Obsolete contamination
No mixed-coordinate PHASE10 primary tables were plotted into Figure3; design schematic is structural only. Quantitative transfer claims in text use 20-seed / corrected sources.
