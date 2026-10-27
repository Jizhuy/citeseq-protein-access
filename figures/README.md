# Figures

| File | Manuscript item | Source / generating script |
| --- | --- | --- |
| `Figure_1.tiff` | Figure 1: representation and contrast framework | `source/Figure1_representation_framework.{pdf,png,svg}`; `scripts/generate_figure1_representation_framework.py` |
| `Figure_2.tiff` | Figure 2: protein access vs model-associated contrast (strict train-only) | `source/Figure2_protein_access_vs_model.{pdf,png}`; `scripts/generate_publication_figures.py` |
| `Figure_3.tiff` | Figure 3: bidirectional transfer and feature-access-matched control | `source/Figure3_transfer_fairness.{pdf,png}`; same |
| `Figure_4.tiff` | Figure 4: protein degradation | `source/Figure4_protein_degradation.{pdf,png}`; same |
| `Figure_5.tiff` | Figure 5: uncertainty and training stochasticity | `source/Figure5_uncertainty_stochasticity.{pdf,png}`; same |
| `Figure_6.tiff` | Figure 6: Lawlor donor-level external validation | `source/Figure6_lawlor_external.{pdf,png}`; same |
| `Supplementary_Figure_S1.{pdf,png,svg}` | Supplementary Figure S1: evaluation-design schematic | `scripts/generate_supplementary_figure_s1.py` |

All data figures are drawn from frozen tables in `results/`; no statistic is recomputed. Supplementary Figure S1 is a design schematic: it plots no result, and its icons are illustrative.

Figure 3 labels the Direction A totalVI − scVI_matched difference as −0.003; the manuscript table gives −0.004. See [`../docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md`](../docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md).
