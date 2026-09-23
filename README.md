# Robust Multi-Omics Integration

Repository: [https://github.com/Jizhuy/robust-multiomics-integration.git](https://github.com/Jizhuy/robust-multiomics-integration.git)

```bash
git clone git@github.com:Jizhuy/robust-multiomics-integration.git
```

## Overview

This project evaluates how much of the observed performance improvement in CITE-seq representation learning is associated with **access to paired protein measurements**, and how much remains as a **residual model-associated advantage** beyond a simple multimodal representation—especially under distribution shift, protein degradation, and external validation.

It is an empirical reliability study of existing representations (PCA baselines, MOFA+, scVI, totalVI), not a new architecture.

## Scientific question

Central comparison on a frozen development cohort (held-out classification macro-F1):

**RNA PCA → RNA + protein concat PCA → totalVI**

with supporting RNA-only (scVI) and protein-only (protein PCA / MOFA+) references.

Terminology used throughout:

| Term | Meaning |
| --- | --- |
| **Protein-access gain** (`G_access`) | `F1_concat − F1_RNA_PCA` |
| **Residual model-associated gain** (`G_assoc`) | `F1_totalVI − F1_concat` |
| **Feature-access-matched transductive concat** | Transfer control that matches access to unlabeled target features only |

Avoid reading these as information-theoretic quantities, pure architecture effects, or causal model attributions.

## Main findings

### Development cohort (frozen point estimates)

| Representation | Macro-F1 |
| --- | ---: |
| RNA PCA | 0.666 |
| Protein PCA | 0.512 |
| Concat PCA | 0.745 |
| MOFA+ | 0.549 |
| scVI | 0.722 |
| totalVI | 0.779 |

Primary descriptive contrasts (cell-level conditional paired bootstrap, *B* = 10,000):

| Contrast | Point | 95% interval | Notes |
| --- | ---: | --- | --- |
| `G_access` = concat − RNA PCA | **+0.079** | [0.043, 0.132] | Bootstrap proportion > 0: **1.0** |
| `G_assoc` = totalVI − concat | **+0.034** | [−0.012, 0.073] | Interval **crosses zero** |
| RNA PCA → totalVI gap | **+0.113** | [0.084, 0.154] | |
| Descriptive recovery ratio ≈ `G_access / gap` | **~0.70** | [0.405, 1.114] | Descriptive only—not variance explained, not causal |

### Transfer under distribution shift

| Setting | Direction A | Direction B |
| --- | ---: | ---: |
| Source-only inductive concat | 0.735 | 0.734 |
| Feature-access-matched transductive concat | 0.636 | 0.767 |
| scVI | 0.656 | 0.601 |
| totalVI | 0.652 | 0.654 |
| totalVI − feature-access concat | ≈ +0.016 | ≈ −0.113 |

Direction B totalVI − scVI ≈ **+0.052** (19/20 seeds favor totalVI).

Matching access to unlabeled target features does **not** yield a stable totalVI advantage across transfer directions. The concat control matches **feature access only**; it does **not** match learning objective, adaptation algorithm, or optimization procedure, and is not algorithmically identical to scArches.

### Protein degradation

**Training-time degradation** (corrupt originally nonzero protein entries; preexisting zeros unchanged; retrain totalVI; five mask/retraining computational replicates at each nonzero level; `model_seed` fixed at 0; clean level *n* = 1):

| Corruption | totalVI macro-F1 |
| ---: | --- |
| 0.00 | 0.779 |
| 0.10 | 0.771 ± 0.010 |
| 0.25 | 0.758 ± 0.012 |
| 0.40 | 0.760 ± 0.012 |
| 0.55 | 0.750 ± 0.016 |
| 0.70 | 0.736 ± 0.020 |
| 0.85 | 0.755 ± 0.011 |

The high-corruption rebound (0.70 → 0.85) occurs within mask/retraining variability and is **not** evidence of improved performance under more severe corruption. These are computational replicates, not biological replicates.

**Post-adaptation target-only degradation** (Direction B totalVI): 0.654 → 0.621 (clean scVI reference: 0.601).

### Uncertainty and stochasticity

- Classifier predictive uncertainty: `U_pred = 1 − max_k p_k`
- Latent posterior variance: mean posterior variance across latent dimensions (**not** labeled epistemic uncertainty here)
- Predictive AUROC is about **0.84–0.86** internally
- Predictive and latent uncertainty respond differently under degradation
- Selective prediction reduces retained-set error but changes cell-type composition
- Stochastic training variation is large enough that small residual model-associated differences need repeated-run context

### External Lawlor validation

- **10 biological donors**; primary biological unit = donor
- totalVI − scVI donor-level present-class macro-F1 difference ≈ **+0.090** (**10/10** donors favor totalVI)
- Deep averages ≈ scVI **0.800**, totalVI **0.888**
- Stronger simple protein-aligned baselines: protein PCA ≈ **0.927**, concat PCA ≈ **0.945**
- Lawlor labels are **protein-gated**

Interpretation: Lawlor supports cross-donor value of **protein-aligned multimodal information**. It does **not** establish universal totalVI or architecture superiority.

Consistency note: Lawlor latent-AUROC computational consistency is **38/50**. The **19/20** seed count belongs only to PBMC Direction B transfer—do not interchange these.

## Datasets

| Cohort | Role | Source |
| --- | --- | --- |
| PBMC10k / PBMC5k | Development | Public 10x Genomics CITE-seq PBMC datasets (via project loaders; no unverified GEO/SRA claimed here) |
| Lawlor Baseline CITE-seq | External | HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005` |

Processed development objects and manifests live under [`data/`](data/). See [`data/README.md`](data/README.md).

## Methods / evaluated representations

| Representation | Modalities | Role |
| --- | --- | --- |
| RNA PCA | RNA | Unimodal baseline |
| Protein PCA | Protein | Unimodal protein baseline |
| Concat PCA | RNA + protein (standardized concat) | Multimodal access baseline |
| MOFA+ | RNA + protein | Factor model reference |
| scVI | RNA | RNA-only VAE |
| totalVI | RNA + protein | Joint CITE-seq VAE |

Downstream probe: logistic regression on held-out cells (development) or transfer targets; primary metric macro-F1 on present classes as specified per analysis.

## Reliability analyses

1. **Protein-access vs residual model-associated contrasts** + cell-level paired bootstrap  
2. **Inductive vs feature-access-matched transductive concat** under PBMC transfer  
3. **Training-time and post-adaptation protein degradation**  
4. **Predictive / latent uncertainty** and selective prediction  
5. **Training stochasticity** (seed replication)  
6. **Lawlor donor-held-out external validation** (protein-gated labels)

Bootstrap intervals are **conditional** on the frozen development embeddings / fits. They are not biological confidence intervals over independent cohorts.

## Repository structure

```
scripts/                         # Phase runners (data, models, transfer, uncertainty, Lawlor)
src/                             # Shared library code
data/                            # Manifests and processed analysis objects
config/                          # Run configs
results/
  revision_round2_methodological_fixes/   # Multimodal baselines, Lawlor audits, degradation controls
  revision_round3_information_vs_model/   # Feature-access-matched transfer control
  presubmission_revision/                 # Paired bootstrap + degradation replicate summary
  phase10_cross_dataset/                  # Transfer tables
  phase11_uncertainty/                    # Uncertainty tables
  phase12a_variance_replication/          # Stochasticity / variance replication
  phase12b_external_validation/           # Lawlor provenance and formal validation
  final_information_vs_model_manuscript/  # Current manuscript package + figures
requirements.txt / environment.yml        # Pinned environment
```

Historical phase outputs under `results/` are retained for provenance; prefer revision-round and presubmission paths for the current scientific framing.

## Key reproducibility scripts

| Analysis | Script |
| --- | --- |
| Internal representation comparison (concat / PCA / MOFA+ framing) | [`results/revision_round2_methodological_fixes/scripts/revision2_simple_multimodal_baseline.py`](results/revision_round2_methodological_fixes/scripts/revision2_simple_multimodal_baseline.py) |
| Transfer (historical phase runner) | [`scripts/10_run_cross_dataset.py`](scripts/10_run_cross_dataset.py) |
| Feature-access-matched concat transfer | [`results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py`](results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py) |
| Training-time protein sparsity / corruption | [`scripts/08_run_sparsity_stress.py`](scripts/08_run_sparsity_stress.py) |
| Post-adaptation target protein corruption | [`results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py`](results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py) |
| Uncertainty | [`scripts/11_run_uncertainty.py`](scripts/11_run_uncertainty.py) |
| Stochasticity / variance replication | [`scripts/12a_run_replicates.py`](scripts/12a_run_replicates.py) |
| Lawlor donor-level analysis | [`results/revision_round2_methodological_fixes/scripts/03_lawlor_donor_level.py`](results/revision_round2_methodological_fixes/scripts/03_lawlor_donor_level.py) |
| Lawlor pairwise consistency | [`results/revision_round2_methodological_fixes/scripts/02_lawlor_pairwise_consistency.py`](results/revision_round2_methodological_fixes/scripts/02_lawlor_pairwise_consistency.py) |
| Presubmission paired bootstrap | [`results/presubmission_revision/bootstrap/presubmission_core_contrast_bootstrap.py`](results/presubmission_revision/bootstrap/presubmission_core_contrast_bootstrap.py) |
| Degradation replicate summary | [`results/presubmission_revision/degradation_audit/DEGRADATION_REPLICATE_AUDIT.md`](results/presubmission_revision/degradation_audit/DEGRADATION_REPLICATE_AUDIT.md) |

Frozen tables for the matched transfer control: [`results/revision_round3_information_vs_model/transfer_control/tables/matched_transfer_control.csv`](results/revision_round3_information_vs_model/transfer_control/tables/matched_transfer_control.csv).

Bootstrap summary: [`results/presubmission_revision/bootstrap/core_contrast_bootstrap_summary.csv`](results/presubmission_revision/bootstrap/core_contrast_bootstrap_summary.csv).

## Environment

Verified development environment (see headers in [`requirements.txt`](requirements.txt) and [`environment.yml`](environment.yml)):

- Python **3.11.7**
- scvi-tools **1.3.3**
- scikit-learn **1.5.2**
- PyTorch **2.14.0** (Apple MPS available in the documented macOS setup; CUDA unavailable on that machine)

Linux CUDA variant: [`environment-linux-cuda.yml`](environment-linux-cuda.yml).

```bash
conda create -n multiomics_robustness python=3.11.7 pip -y
conda activate multiomics_robustness
pip install -r requirements.txt
```

## Reproducing the analyses

There is **no single end-to-end one-command pipeline**. Analyses were built in stages; frozen outputs under `results/` are the scientific record.

Typical workflow:

1. **Data** — prepare / validate CITE-seq objects (`scripts/01_prepare_data.py`; see `data/`).
2. **Baselines & models** — historical phase runners under `scripts/` (`02`–`08`, `10`–`12b`) write embeddings and tables into `results/`.
3. **Current framing controls** — run revision scripts above for concat baselines, matched transfer, and test-time corruption (require prior embeddings / annotated objects).
4. **Presubmission summaries** — bootstrap script reconstructs held-out predictions from frozen development embeddings (does **not** retrain VAEs); degradation audit aggregates historical sparsity tables.

Expect nontrivial runtime and GPU/MPS use for VAE stages. Prefer reading frozen CSVs/figures when verifying manuscript numbers.

## Data availability

- **PBMC10k / PBMC5k**: public 10x Genomics CITE-seq PBMC datasets used by the project loaders.
- **Lawlor**: HCA dataset [`efea6426-510a-4b60-9a19-277e52bfa815`](https://explore.data.humancellatlas.org/projects/efea6426-510a-4b60-9a19-277e52bfa815); ENA [`PRJEB40376`](https://www.ebi.ac.uk/ena/browser/view/PRJEB40376) / `ERP124005`.
- Processed analysis objects and checksums: `data/` and `results/phase12b_external_validation/`.

## Code availability

Source and analysis scripts: this repository  
[https://github.com/Jizhuy/robust-multiomics-integration.git](https://github.com/Jizhuy/robust-multiomics-integration.git)

## Manuscript status

**Working title:** Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift

**Status:** Pre-submission research manuscript

Current package (figures, tables, draft): [`results/final_information_vs_model_manuscript/`](results/final_information_vs_model_manuscript/)

## Citation / contact

Author / maintainer: **Jizhuy** ([yangjizhu@outlook.com](mailto:yangjizhu@outlook.com))

If you use this repository, please cite the manuscript once deposited and link this GitHub URL.
