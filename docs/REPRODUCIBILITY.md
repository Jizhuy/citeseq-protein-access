# Reproducibility

This document says how each reported analysis was produced and how to rerun it. Frozen outputs are under `results/`. A value-by-value index is in [`RESULT_MANIFEST.md`](RESULT_MANIFEST.md).

There is no single one-command end-to-end pipeline. Inspecting frozen tables and regenerating figures needs neither data nor GPU. Rerunning the deep-model analyses needs the public datasets ([`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)) and substantial compute.

## 1. Analysis taxonomy

| Label | Analysis | Status | Frozen outputs |
| --- | --- | --- | --- |
| **PRIMARY** | Strict train-only internal evaluation. All representations, including scVI_matched and totalVI, are fitted on 7,573 high-confidence training cells; 1,921 held-out cells are embedded without parameter updates. | COMPLETE — PRIMARY ANALYSIS | `results/primary/` |
| **SECONDARY** | Historical feature-transductive internal analysis. Deep models are fitted on the features (not labels) of all development cells. Information-access sensitivity only. | Secondary | `results/supplementary/information_access_sensitivity/` |
| **TRANSFER** | Bidirectional PBMC10k ↔ PBMC5k transfer with scArches query adaptation, 20 seeds per direction. scArches uses unlabeled target features, so it is transductive adaptation. | Complete | `results/transfer/replicate_runs/` |
| **TRANSDUCTIVE CONTROL** | Feature-access-matched transductive concat PCA, fitted on source + target features. It matches target-feature access, not the learning objective or the architecture. The inductive concat PCA is fitted on source cells only. | Complete | `results/transfer/matched_concat_control/` |
| **DEGRADATION (training-time)** | Protein counts masked before totalVI training, using the historical feature-transductive pipeline, with 5 mask replicates per level. | Complete | `results/degradation/training_time/` |
| **DEGRADATION (post-adaptation)** | Target protein counts corrupted after query adaptation, 20 seeds. | Complete | `results/degradation/test_time_corruption/` |
| **UNCERTAINTY** | Predictive uncertainty (classifier probabilities) and latent posterior dispersion (VAE posterior variance), reported separately; selective prediction / risk–coverage. | Complete | `results/uncertainty/` |
| **STOCHASTICITY** | Replicated crossed design: source subset × model seed × run replicate; corrected variance decomposition. | Complete | `results/stochasticity/` |
| **EXTERNAL VALIDATION** | Lawlor PBMC CITE-seq (Baseline condition), donor-held-out folds; labels are protein-gated. | Complete | `results/donor_validation/` |
| **BIOLOGICAL UNIT** | 10 donors. Donor-level paired deltas with an exact two-sided sign test (10/10, P = 0.00195). | — | `results/donor_validation/summary/` |
| **COMPUTATIONAL CONSISTENCY** | 38/50 latent-error-signal consistency over 50 matched fold × seed pairs (5 folds × 10 seeds). A computational count, not a biological replication count. | — | `results/donor_validation/statistics/` |

**Not run:**

- totalVI retraining under complete marker-channel dropout;
- an additional newly trained multimodal comparator;
- a second external, non-protein-gated validation cohort.

The marker-channel dropout sensitivity is PCA-only.

## 2. Environment

```bash
conda create -n multiomics_robustness python=3.11.7 pip -y
conda activate multiomics_robustness
pip install -r requirements.txt
python scripts/print_environment.py   # writes results/supplementary/environment/runtime_environment.json
```

- Python 3.11.7, scvi-tools 1.3.3, PyTorch 2.14.0, scikit-learn 1.5.2, scanpy 1.10.3, anndata 0.11.4, numpy 1.26.4, pandas 2.2.2, scipy 1.13.1, statsmodels 0.15.0, mofapy2 0.7.5.
- Confirmed on macOS 14.6.1 (Apple M1 Pro, MPS). The GPU transfer experiments were run on CUDA.
- `scripts/prepare_lawlor_baseline_matrices.R` requires R with the `Matrix` package.
- **CUDA library path (optional):** on WSL or other systems where the CUDA driver libraries are not on the default path, set `CUDA_LIB_DIR` (for example `export CUDA_LIB_DIR=/usr/lib/wsl/lib`) before running `scripts/run_lawlor_donor_validation.py`. If it is unset, nothing is changed. `src/utils/device.py` probes the standard WSL `nvidia-smi` location only when it exists.

All scripts resolve paths relative to the repository root. Run them from the root.

## 3. Configuration and output locations

`config/config.yaml` holds data sources, preprocessing, model hyperparameters and seeds. Its `paths` block points the config-driven pipelines (data preparation, reference annotation, historical feature-transductive models, training-time degradation) to `results/tables/`, `results/embeddings/`, `results/models/` and `results/logs/`. These intermediates are regenerated on rerun, are not shipped, and are git-ignored. The frozen copies used in the manuscript live in the analysis folders (`results/primary/`, `results/transfer/`, ...).

One deliberate exception: `scripts/run_test_time_protein_corruption.py` reads the frozen training-time summary `results/degradation/training_time/protein_sparsity_primary_delta_summary.csv` for its training-time vs test-time comparison table.

## 4. Order of execution

### 4.1 Data and labels (all analyses)

```bash
python scripts/prepare_data.py                 # PBMC10k + PBMC5k → data/processed/pbmc_cite_combined_inner.h5ad
python scripts/run_reference_annotation.py     # Seurat v4 reference label transfer → data/annotations/, *_annotated.h5ad
```

`data/annotations/pbmc_reference_annotations.csv` is shipped, so `--skip-transfer` reuses the published labels.

### 4.2 PRIMARY: strict train-only internal analysis

```bash
python scripts/run_strict_train_only_internal.py   # → results/primary/strict_train_only/
```

### 4.3 SECONDARY: information-access sensitivity (historical feature-transductive)

```bash
python scripts/run_feature_transductive_models.py --model scvi
python scripts/run_feature_transductive_models.py --model scvi_matched
python scripts/run_feature_transductive_models.py --model totalvi
python scripts/run_feature_transductive_bootstrap.py
python scripts/cache_feature_transductive_predictions.py
python scripts/run_class_stratified_bootstrap.py
```

`information_access_sensitivity_comparison.csv` places the frozen strict train-only and feature-transductive estimates side by side. It is a frozen derived artifact; no standalone generator is retained, because its approximate entries and role text were hand-formatted (see `docs/RESULT_MANIFEST.md` §9).

### 4.4 Simple baselines and controls

```bash
python scripts/run_simple_multimodal_baselines.py   # development, transfer and Lawlor PCA baselines → results/supplementary/controls/
python scripts/run_marker_channel_dropout.py        # PCA only → results/supplementary/marker_channel_dropout/
```

The historical MOFA+ development value (0.549) is frozen in `results/supplementary/controls/tables/simple_multimodal_development.csv`. MOFA+ is not re-run.

### 4.5 TRANSFER and TRANSDUCTIVE CONTROL

```bash
python scripts/run_bidirectional_transfer.py --stage gate1
python scripts/run_bidirectional_transfer.py --stage gate2
python scripts/run_bidirectional_transfer.py --stage seeds --seeds 1-19
python scripts/run_bidirectional_transfer.py --stage crossed
python scripts/analyze_bidirectional_transfer.py --stage all
python scripts/analyze_crossed_design.py --stage all
python scripts/run_matched_concat_transfer.py --run
```

### 4.6 DEGRADATION

```bash
python scripts/run_protein_degradation.py --full                      # training-time; needs the totalVI fit from 4.3
python scripts/run_test_time_protein_corruption.py --run --gate C     # post-adaptation, seeds 0-19
python scripts/summarize_class_degradation.py
```

`results/degradation/training_time/training_time_degradation_summary.csv` summarises `protein_sparsity_primary_delta.csv` (see `DEGRADATION_REPLICATE_AUDIT.md`) and is regenerated byte-identically by:

```bash
python scripts/build_training_time_degradation_summary.py
```

### 4.7 UNCERTAINTY

```bash
python scripts/run_uncertainty_selective_prediction.py
python scripts/summarize_selective_macro_f1.py
```

`results/uncertainty/calibration/pbmc_20seed_reliability_summary.csv` aggregates the frozen per-seed transfer table and is regenerated byte-identically by:

```bash
python scripts/build_pbmc_reliability_summary.py
```

### 4.8 STOCHASTICITY

```bash
python scripts/run_training_stochasticity.py
python scripts/analyze_training_stochasticity.py --stage all
python scripts/run_variance_decomposition.py
```

### 4.9 EXTERNAL VALIDATION (Lawlor)

```bash
# in data/external/lawlor/ (HCA contributor matrices; see DATA_AVAILABILITY.md)
Rscript ../../../scripts/prepare_lawlor_baseline_matrices.R     # Baseline-condition singlets → baseline_mm/
# Matrix-Market → lawlor_baseline_sng_shared.h5ad (see note below)
python scripts/run_lawlor_donor_validation.py --gate all_seeds
python scripts/analyze_lawlor_donor_validation.py
python scripts/run_lawlor_donor_level.py
python scripts/run_lawlor_pairwise_consistency.py
python scripts/run_lawlor_pca_per_donor.py
python scripts/summarize_lawlor_donor_statistics.py
```

The Python step that converts the Matrix-Market export into `results/donor_validation/lawlor_runs/processed/lawlor_baseline_sng_shared.h5ad` (gene-symbol mapping and the shared feature set) is not preserved as a standalone script. Its inputs, filters and outputs are recorded in `results/donor_validation/lawlor_runs/logs/data_prep_provenance.md` and `data_prep_summary.json`, with the gene mapping in `results/donor_validation/lawlor_dataset_audit/gene_mapping_external.csv`.

### 4.10 Figures (no data or GPU needed)

```bash
python scripts/generate_figure1_representation_framework.py   # figures/source/Figure1_*
python scripts/generate_publication_figures.py                # figures/source/Figure2-6_* from frozen tables
python scripts/generate_supplementary_figure_s1.py            # figures/Supplementary_Figure_S1.*
```

During migration these three scripts were run, and their outputs were pixel-identical to the shipped figure sources.

## 5. Seeds

| Analysis | Seed(s) |
| --- | --- |
| Development split, PCA, classifier | 0 |
| Strict train-only scVI_matched / totalVI | 0 |
| Primary paired bootstrap (B = 10,000) | 20260923 |
| Class-stratified bootstrap (B = 10,000) | 20260924 |
| Transfer deep models | 0–19 per direction |
| Lawlor donor-held-out runs | 5 folds × 10 seeds × 2 models = 100 runs |

## 6. Approximate compute

| Step | Notes |
| --- | --- |
| Figure regeneration | Seconds, CPU |
| Prediction cache + class-stratified bootstrap | ~5 min + ~2–3 min, CPU |
| Marker-channel dropout (PCA) | ~10–30 min, CPU |
| Strict train-only internal | Two VAE fits on one development set; GPU/MPS recommended |
| Transfer, degradation, stochasticity, Lawlor | Many VAE fits; GPU recommended (hours to days in total) |

## 7. Interpretation guard-rails

- Bootstrap intervals are conditional, cell-level resampling intervals, not biological confidence intervals.
- 38/50 (Lawlor fold × seed pairs) and 19/20 (PBMC transfer seeds) are computational consistency counts, not donor *n*.
- G_assoc is a descriptive model-associated contrast, not a causal architecture effect.
- Lawlor labels are protein-gated; the external result does not establish universal totalVI superiority.

## 8. Historical stage labels in frozen logs

Some frozen logs and script docstrings keep the stage labels used during development. They mean:

| Label in frozen text | Meaning |
| --- | --- |
| PHASE 1–2 | Data loading and preparation |
| PHASE 3 | Historical scVI baseline (initial scVI / scVI_matched / totalVI fits) |
| PHASE 6 | Reference-based cell-type annotation |
| PHASE 8 | Training-time protein degradation |
| PHASE 10 / 10B | Single-seed cross-dataset transfer / source-size subsampling |
| PHASE 11 | Single-seed uncertainty analysis |
| PHASE 11B | GPU transfer experiments (20-seed bidirectional transfer and crossed design) |
| PHASE 12A | Stochasticity experiment (replicated crossed design) |
| PHASE 12B | Lawlor donor-level external validation |

These labels have no meaning beyond identifying the run that produced a file. File and directory names in this repository do not use them.
