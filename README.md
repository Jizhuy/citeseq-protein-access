# CITE-seq Multimodal Reliability

Reproducibility repository for:

*Evaluating protein-access gains and residual totalVI-associated contrasts in CITE-seq representation learning under distribution shift*

**Author:** Jizhu Yang

This repository contains code, frozen result tables, publication figures, and supplementary materials supporting an empirical evaluation of protein-access gains and residual totalVI-associated contrasts under matched information access, distribution shift, protein degradation, uncertainty analysis, training stochasticity, and donor-level external validation.

This repository (`Jizhuy/citeseq-protein-access`) hosts a synchronized copy of the finalized publication and reproducibility materials. The journal-submitted Code Availability URL cites the companion repository `https://github.com/Jizhuy/citeseq-multimodal-reliability`, which remains the manuscript-cited location and is unchanged by this synchronization. Pre-synchronization development history for this repository is preserved on branch `archive/legacy-pre-publication-sync-2026-10-27` and tag `legacy-pre-sync-2026-10-27`.

---

## 1. Overview

Paired CITE-seq measurements of RNA and surface protein are often used to motivate multimodal deep generative models. A gain from such a model can come from two different sources: having access to protein measurements at all, or the specific multimodal model. This study separates the two with descriptive contrasts evaluated by a common downstream classifier (cell-type macro-F1), and stress-tests them under distribution shift, protein degradation, uncertainty analysis, training stochasticity, and donor-level external validation.

**Model positioning.**

- This study evaluates existing representations and evaluation controls. It does **not** introduce a new neural architecture.
- totalVI is an established, published model.
- scVI_matched is a capacity-matched configuration of the published scVI model (RNA only).
- G_assoc is a descriptive model-associated contrast, **not** a causal architecture effect.

## 2. Scientific question

How much of the macro-F1 advantage of a multimodal deep representation (totalVI) over an RNA-only reference is explained by **protein access** alone, as captured by a simple concatenated RNA–protein PCA, and how much remains as a **residual totalVI-associated contrast** once information access is matched?

| Contrast | Definition |
| --- | --- |
| G_access | F1(concat RNA–protein PCA) − F1(RNA PCA) |
| G_assoc | F1(totalVI) − F1(concat RNA–protein PCA) |
| Gap_total | F1(totalVI) − F1(RNA PCA) |
| R_access | G_access / Gap_total (secondary descriptive ratio; not bounded to [0, 1]) |

## 3. Primary results

Primary analysis: **strict train-only** internal evaluation on the PBMC development cohort (7,573 high-confidence training cells / 1,921 held-out high-confidence cells). Every representation, including scVI_matched and totalVI, was fitted on training cells only.

| Representation | Macro-F1 |
| --- | ---: |
| RNA PCA | 0.667 |
| Concat RNA-protein PCA | 0.745 |
| scVI_matched | 0.709 |
| totalVI | 0.741 |

| Contrast | Estimate | 95% interval |
| --- | ---: | --- |
| G_access | +0.078 | [0.042, 0.131] |
| G_assoc | -0.004 | [-0.060, 0.030] |
| Gap_total | +0.074 | [0.029, 0.130] |

Intervals are cell-level conditional paired-bootstrap percentile intervals (B = 10,000). They are not biological confidence intervals. For G_assoc, prop>0 = 0.355. R_access ≈ 1.05 is reported only as a secondary descriptive ratio.

Under strict train-only matched information access, the protein-access contrast was substantial and resampling-stable, whereas the residual totalVI-associated contrast was near null and imprecisely estimated.

Other headline results (frozen values; see [`docs/RESULT_MANIFEST.md`](docs/RESULT_MANIFEST.md) for file-level sources):

| Analysis | Result |
| --- | --- |
| Transfer, Direction A (PBMC10k → PBMC5k) | inductive concat 0.735; feature-access-matched transductive concat 0.636; scVI_matched 0.656; totalVI 0.652 |
| Transfer, Direction B (PBMC5k → PBMC10k) | inductive concat 0.734; feature-access-matched transductive concat 0.767; scVI_matched 0.601; totalVI 0.654 |
| Lawlor donor-level external validation (10 donors; donor-mean present-class macro-F1) | RNA PCA 0.690; protein PCA 0.923; concat PCA 0.942; scVI_matched 0.794; totalVI 0.884 |
| Lawlor totalVI − scVI_matched | mean donor Δ ≈ +0.090; 10/10 donors favour totalVI; exact two-sided sign test P = 0.00195 |

Lawlor labels are protein-gated. This favours protein-aligned baselines; the result supports the value of the multimodal measurement, not universal totalVI superiority.

## 4. Study design

| Component | Role | Status |
| --- | --- | --- |
| Strict train-only internal analysis | **PRIMARY** analysis | COMPLETE — PRIMARY ANALYSIS |
| Historical feature-transductive internal analysis | Secondary information-access sensitivity (deep models fitted on features of all development cells) | Secondary only |
| Bidirectional cross-dataset transfer (scArches query adaptation, 20 seeds) | Distribution-shift stress test | Complete |
| Feature-access-matched transductive concat PCA | Transfer control matching target-feature access | Complete |
| Training-time protein degradation | Antibody-signal dependence (historical feature-transductive pipeline; 5 mask replicates per level) | Complete |
| Post-adaptation test-time protein corruption | Degradation after query adaptation (20 seeds) | Complete |
| Uncertainty and selective prediction | Predictive uncertainty, latent posterior dispersion, risk–coverage | Complete |
| Training stochasticity | Replicated crossed design: source subset × model seed × run replicate | Complete |
| Lawlor donor-level external validation | Biological unit: 10 donors | Complete |

**Transfer note.** scArches query adaptation uses unlabeled target features, so it is transductive adaptation. The transductive concat baseline matches that target-feature access but not the learning objective or the architecture.

**Not run** (stated explicitly, not implied as complete):

- totalVI retraining under complete marker-channel dropout (the marker-channel dropout sensitivity is PCA-only);
- an additional newly trained multimodal comparator;
- a second external, non-protein-gated validation cohort.

Earlier development stages are described neutrally here: the *historical scVI baseline* (initial scVI, scVI_matched and totalVI fits on the development cohort), the *GPU transfer experiments* (cross-dataset scArches transfer runs), and the *stochasticity experiment* (replicated crossed design).

## 5. Repository structure

```
citeseq-protein-access/   # synchronized publication tree (companion: citeseq-multimodal-reliability)
├── README.md, CITATION.cff, LICENSE, .gitignore
├── environment.yml, requirements.txt      # pinned software environment
├── config/                                # analysis configuration (no machine paths)
├── data/annotations/                      # reference-transfer cell labels (derived); raw data not included
├── src/                                   # minimal library code used by the scripts
│   ├── data/  models/  evaluation/  experiments/  annotation/  utils/
├── scripts/                               # analysis entry points (see section 6)
├── results/
│   ├── primary/                           # strict train-only internal analysis (PRIMARY)
│   ├── transfer/                          # bidirectional transfer + feature-access-matched control
│   ├── degradation/                       # training-time and post-adaptation protein degradation
│   ├── uncertainty/                       # reliability, selective prediction, statistical units
│   ├── stochasticity/                     # replicated crossed design and variance decomposition
│   ├── donor_validation/                  # Lawlor donor-level external validation
│   └── supplementary/                     # controls, marker-channel dropout, environment,
│       └── information_access_sensitivity/  # SECONDARY information-access sensitivity
├── figures/
│   ├── Figure_1.tiff … Figure_6.tiff      # final main figures (600 dpi)
│   ├── Supplementary_Figure_S1.{pdf,png,svg}
│   └── source/                            # vector sources of Figures 1-6
├── supplementary/
│   ├── Supplementary_Information.pdf
│   └── Supplementary_Data.zip
└── docs/                                  # reproducibility, data availability, manifests, audits
```

Each `results/<analysis>/README.md` states the role and status of its tables.

## 6. Reproducing analyses

Frozen tables under `results/` are the reference outputs; you do not need to rerun anything to inspect a reported value. Full reruns need the public datasets (section 7) and substantial compute (CUDA GPU or Apple MPS for the VAE models).

| Analysis | Script(s) | Output |
| --- | --- | --- |
| Data preparation | `scripts/prepare_data.py` | `data/processed/` |
| Reference cell-type labels | `scripts/run_reference_annotation.py` | `data/annotations/pbmc_reference_annotations.csv`, annotated `.h5ad` |
| **Strict train-only internal (PRIMARY)** | `scripts/run_strict_train_only_internal.py` | `results/primary/strict_train_only/` |
| Simple PCA baselines (development, transfer, Lawlor) | `scripts/run_simple_multimodal_baselines.py` | `results/supplementary/controls/` |
| Bidirectional transfer (20 seeds) | `scripts/run_bidirectional_transfer.py`, `scripts/analyze_bidirectional_transfer.py` | `results/transfer/replicate_runs/` |
| Feature-access-matched transfer control | `scripts/run_matched_concat_transfer.py` | `results/transfer/matched_concat_control/` |
| Historical feature-transductive models | `scripts/run_feature_transductive_models.py` | `results/models/`, `results/embeddings/` (not shipped) |
| Training-time protein degradation | `scripts/run_protein_degradation.py` | `results/tables/` → frozen copies in `results/degradation/training_time/` |
| Post-adaptation protein corruption | `scripts/run_test_time_protein_corruption.py` | `results/degradation/test_time_corruption/` |
| Uncertainty / selective prediction | `scripts/run_uncertainty_selective_prediction.py`, `scripts/summarize_selective_macro_f1.py` | `results/uncertainty/` |
| Training stochasticity | `scripts/run_training_stochasticity.py`, `scripts/analyze_crossed_design.py`, `scripts/analyze_training_stochasticity.py`, `scripts/run_variance_decomposition.py` | `results/stochasticity/` |
| Lawlor donor validation | `scripts/run_lawlor_donor_validation.py`, `scripts/analyze_lawlor_donor_validation.py`, `scripts/run_lawlor_donor_level.py`, `scripts/run_lawlor_pairwise_consistency.py`, `scripts/run_lawlor_pca_per_donor.py`, `scripts/summarize_lawlor_donor_statistics.py` | `results/donor_validation/` |
| Class-stratified bootstrap (secondary sensitivity) | `scripts/cache_feature_transductive_predictions.py`, `scripts/run_class_stratified_bootstrap.py` | `results/supplementary/information_access_sensitivity/` |
| Marker-channel dropout (PCA only) | `scripts/run_marker_channel_dropout.py` | `results/supplementary/marker_channel_dropout/` |
| Summary tables from frozen per-run tables (no model run) | `scripts/build_training_time_degradation_summary.py`, `scripts/build_pbmc_reliability_summary.py` | `results/degradation/training_time/`, `results/uncertainty/calibration/` |
| Figures 1-6, Supplementary Figure S1 | `scripts/generate_figure1_representation_framework.py`, `scripts/generate_publication_figures.py`, `scripts/generate_supplementary_figure_s1.py` | `figures/` |

The figure scripts read only frozen tables and can be run without data or GPU:

```bash
conda activate multiomics_robustness
python scripts/generate_publication_figures.py
```

The full order of execution, seeds, compute notes and the analysis taxonomy are in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## 7. Datasets

| Cohort | Role | Provenance |
| --- | --- | --- |
| PBMC10k and PBMC5k CITE-seq | Development (internal) and bidirectional transfer | Public scvi-tools / 10x Genomics CITE-seq example data (`pbmc_10k_protein_v3.h5ad`, `pbmc_5k_protein_v3.h5ad`). No GEO/SRA accession is claimed. |
| Lawlor et al. PBMC CITE-seq (Baseline condition) | Donor-level external validation | Human Cell Atlas project `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005` |

## 8. Data availability

Raw public datasets are **not redistributed** in this repository. Download instructions, expected file names and local placement are in [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md). No private biological data are included.

## 9. Software environment

Verified versions: Python 3.11.7, scvi-tools 1.3.3, scikit-learn 1.5.2, scanpy 1.10.3, anndata 0.11.4, PyTorch 2.14.0, numpy 1.26.4, pandas 2.2.2, scipy 1.13.1, mofapy2 0.7.5.

```bash
conda create -n multiomics_robustness python=3.11.7 pip -y
conda activate multiomics_robustness
pip install -r requirements.txt
```

`environment.yml` lists the core pins; `requirements.txt` is the full pinned tree. `results/supplementary/environment/runtime_environment.json` records the versions observed at run time.

## 10. Supplementary material

| Item | Path |
| --- | --- |
| Supplementary Information | [`supplementary/Supplementary_Information.pdf`](supplementary/Supplementary_Information.pdf) |
| Supplementary Data (frozen tables) | [`supplementary/Supplementary_Data.zip`](supplementary/Supplementary_Data.zip). The Supplementary Information refers to this archive by its journal-package file name (`05_GPB_Supplementary_Data.zip`); the contents are identical. |
| Supplementary Figure S1 (study-design schematic; icons are illustrative; no results are plotted) | [`figures/Supplementary_Figure_S1.pdf`](figures/Supplementary_Figure_S1.pdf) |
| Figure ↔ manuscript mapping (Figures 1–6, S1) | [`figures/README.md`](figures/README.md) |
| Result ↔ file ↔ script ↔ table/figure index | [`docs/RESULT_MANIFEST.md`](docs/RESULT_MANIFEST.md) |
| Known numerical presentation notes | [`docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md`](docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md) |
| Migration manifest and audit | [`docs/MIGRATION_MANIFEST.md`](docs/MIGRATION_MANIFEST.md), [`docs/MIGRATION_AUDIT.md`](docs/MIGRATION_AUDIT.md) |

## 11. Citation

Please cite the repository using [`CITATION.cff`](CITATION.cff):

Jizhu Yang. *Evaluating protein-access gains and residual totalVI-associated contrasts in CITE-seq representation learning under distribution shift.* Manuscript-cited code repository: https://github.com/Jizhuy/citeseq-multimodal-reliability (this repository holds a synchronized copy of the same materials).

No archival DOI has been issued yet.

## 12. License

MIT License — see [`LICENSE`](LICENSE).

The MIT License applies to the code and repository-authored materials in this repository unless otherwise stated. It does not relicense, and implies no ownership of, the third-party public datasets analysed here (PBMC10k/PBMC5k, the Seurat v4 reference, Lawlor et al.). Those datasets remain subject to their original source terms; see [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

## 13. Contact

Jizhu Yang — [yangjizhu@outlook.com](mailto:yangjizhu@outlook.com) (public repository contact)
