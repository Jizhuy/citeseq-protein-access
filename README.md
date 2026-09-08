# Robust and Uncertainty-Aware Deep Generative Modeling for Single-Cell Multi-Omics Integration

Stage 5A–5C research repository. The current implementation covers **PHASE 1–8**: validated PBMC CITE-seq data, scVI default, totalVI, a matched scVI fairness control, independent Seurat v4 annotation, MOFA+, RNA PCA, and a protein sparsity / corruption stress test. True missing-modality experiments are not run yet.

## 1. Research motivation

Deep generative models such as scVI and totalVI are widely used to represent single-cell RNA and CITE-seq protein measurements. Published benchmarks, including Hu et al. (2024), compare prediction and integration accuracy under relatively standard conditions. They do not fully answer whether those representations remain biologically reliable when protein measurements become sparse, when the protein modality is missing for many cells, or when datasets shift.

This project is a controlled stress test of existing models, not a new architecture and not a reproduction of the 2024 benchmark.

## 2. Main research questions

**RQ1.** Does RNA + protein multimodal modeling yield a more biologically meaningful latent representation than RNA-only modeling?

Comparisons planned after PHASE 3:

- PCA (RNA; multimodal concatenation only if RNA cannot trivially dominate protein)
- scVI (RNA only)
- MOFA+ (RNA + protein, with MOFA+-appropriate preprocessing)
- totalVI (RNA + protein counts)

**RQ2.** How does increasing protein sparsity/noise affect latent quality and downstream biological metrics?

**RQ3.** How does increasing the fraction of cells with a fully missing protein modality affect latent structure, cell-type conservation, and protein imputation?

Later stages (not implemented):

- cross-dataset generalization
- uncertainty calibration
- RNA + ATAC with MultiVI

## 3. Dataset

Development data are the public 10x PBMC CITE-seq datasets distributed by scvi-tools:

- PBMC10k protein v3 (`pbmc_10k_protein_v3.h5ad`)
- PBMC5k protein v3 (`pbmc_5k_protein_v3.h5ad`)

Validated dimensions after PHASE 2:

| Object | Cells | Genes | Proteins |
| --- | ---: | ---: | ---: |
| PBMC10k | 6855 | 16727 | 14 |
| PBMC5k | 3994 | 16581 | 29 |
| Combined inner (default analysis set) | 10849 | 15792 | 14 shared |
| Combined outer | 10849 | 15792 | 29 union |

The PBMC10k protein panel is a **strict subset** of PBMC5k (14 shared, 0 unique to 10k, 15 unique to 5k). The inner-join object is the default modeling set. Unique proteins are retained in `pbmc5k_cite.h5ad` and `pbmc_cite_combined_outer.h5ad` with NaN, not zero, for unmeasured entries.

The official files contain **no cell-type annotation**. `obs["cell_type"]` is `unknown` until an external labeled reference is added. Do not use that placeholder as a biological label.

The loader uses the same URLs as `scvi.data.pbmcs_10x_cite_seq`, but does **not** call that helper as the primary API. The helper intersects genes silently and, for an outer protein join, fills missing proteins with zero. This repository:

1. loads each file separately
2. documents RNA and protein overlap
3. writes an inner-join analysis object (shared genes and shared proteins)
4. writes an outer-join object in which unmeasured proteins are NaN, with `obsm["protein_observed_mask"]`

Raw RNA counts stay in `X` and `layers["counts"]`. Protein counts stay in `obsm["protein_counts"]`.

See `data/README.md` for the on-disk layout.

## 4. Methods

PHASE 1–2 implement data provenance and validation only.

Planned model roles (later phases):

| Model | Modalities | Latent object | Notes |
| --- | --- | --- | --- |
| PCA | RNA, optionally concatenated protein | linear factors | RNA dimension must not be allowed to dominate protein |
| scVI | RNA counts | `z_RNA` | batch covariate if present |
| totalVI | RNA + protein counts | `z_joint` | count-based VAE |
| MOFA+ | processed RNA + protein | factors | not mathematically equivalent to totalVI |

UMAP is for visualization. Quantitative representation metrics will use the latent factors, not UMAP coordinates.

## 5. Experimental design

1. **Baseline (PHASE 3–6):** train scVI, totalVI, MOFA+, and PCA on the same inner-join PBMC object; compare ARI, NMI, silhouette, and batch ASW separately. Do not collapse batch-removal and biology scores into one number.
2. **Sparsity stress (PHASE 7):** copy the AnnData; mask 0/10/25/50/75% of nonzero protein counts to zero using a shared mask across models.
3. **Missing modality (PHASE 8):** hide entire protein profiles for 0/25/50/75/90% of cells; evaluate imputation only against held-out ground truth.
4. **RNA-correlated vs uncorrelated proteins:** after imputation exists, compare accuracy by RNA–protein association. Do not copy a published threshold without documenting it.

Statistical rules:

- never evaluate protein prediction only on training cells
- preserve ground-truth matrices before perturbation
- use the same train/test split and perturbation mask across methods
- report mean, SD, and per-seed values (5 seeds)

## 6. Reproducibility

- Config: `config/config.yaml`
- Seeds: Python, NumPy, and PyTorch via `src/utils/seed.py`
- Environment snapshot: `results/logs/environment.json` after `01_prepare_data.py`
- Train/test split: `results/tables/combined_inner_train_test_split.csv`

Pinned packages will be written to `environment.yml` and `requirements.txt` after the working environment is confirmed.

### Environment notes (this machine)

Inspected before any training:

- OS: macOS 14 (Darwin 23.6.0), Apple M1 Pro, 16 GB RAM
- NVIDIA CUDA: **not available** (`nvidia-smi` missing; no CUDA GPU)
- Apple Metal: Metal 3, 14-core GPU
- Intended accelerator for later training: PyTorch MPS, falling back to CPU
- Conda environment: `multiomics_robustness`, Python 3.11.7
- scvi-tools pinned at 1.3.3 (not the latest 1.5.x, which requires Python 3.12+)

```bash
conda activate multiomics_robustness
cd multiomics_robustness
python scripts/01_prepare_data.py
```

## 7. How to run each stage

```bash
conda activate multiomics_robustness
cd multiomics_robustness

# PHASE 1-2
python scripts/01_prepare_data.py

# PHASE 3: RNA-only scVI baseline
python scripts/02_run_baseline.py --model scvi --seed 0

# PHASE 4: RNA+protein totalVI baseline (does not train MOFA+)
python scripts/02_run_baseline.py --model totalvi --seed 0

# PHASE 5A–5B: matched scVI fairness control + annotation infrastructure
python scripts/05_run_fairness_and_annotation.py
python scripts/02_run_baseline.py --model scvi_matched --seed 0

# PHASE 6: independent Seurat v4 annotation + biological evaluation
python scripts/06_run_phase6.py

# PHASE 7: MOFA+ multimodal factor baseline + RNA PCA control
python scripts/07_run_mofa.py

# PHASE 8: protein sparsity / corruption stress test (smoke, then full)
python scripts/08_run_sparsity_stress.py --smoke
python scripts/08_run_sparsity_stress.py --full

# PHASE 9 (placeholder; true missing modality — do not run yet)
python scripts/04_run_missing_modality.py
```

## 8. Output structure

```
results/
  tables/     baseline_results.csv, three-model comparison, latent geometry, training histories
  figures/    scVI / scVI_matched / totalVI UMAPs and training curves
  embeddings/ scvi_seed0_*, scvi_matched_seed0_*, totalvi_seed0_*
  models/     scvi_seed0/, scvi_matched_seed0/, totalvi_seed0/
  logs/       manifests, fairness_control_interpretation.md, annotation_source_audit.md
  exploratory/ protein-marker summaries only; not validated labels
data/
  raw/        original scvi-tools h5ad files
  processed/  validated AnnData objects with raw counts preserved
  annotations/ external label template; no fabricated rows
```

## Repository layout

Analysis logic lives in `src/`. Scripts are thin wrappers. Notebooks are exploratory only.

```
multiomics_robustness/
├── README.md
├── environment.yml
├── requirements.txt
├── config/
├── data/
├── src/
│   ├── data/
│   ├── models/
│   ├── evaluation/
│   ├── experiments/
│   └── utils/
├── scripts/
├── results/
└── notebooks/
```
