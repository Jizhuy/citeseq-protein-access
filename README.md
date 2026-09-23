# Robust Multi-Omics Integration

Repository: [https://github.com/Jizhuy/robust-multiomics-integration.git](https://github.com/Jizhuy/robust-multiomics-integration.git)

```bash
git clone git@github.com:Jizhuy/robust-multiomics-integration.git
```

## Overview

This project evaluates how much of the observed performance improvement in CITE-seq representation learning is associated with **access to paired protein measurements**, and how much remains as a **residual totalVI-associated contrast** beyond a simple RNA+protein representation—especially under distribution shift, protein degradation, and external donor validation.

It is an empirical reliability study of existing representations (PCA baselines, MOFA+, scVI_matched, totalVI), not a new architecture.

## Scientific question

Central comparison (frozen development cohort, held-out macro-F1):

**RNA PCA → RNA + protein concat PCA → totalVI**

| Contrast | Definition |
| --- | --- |
| **Protein-access contrast** `G_access` | `F1_concat − F1_RNA_PCA` |
| **Residual totalVI-associated contrast** `G_assoc` | `F1_totalVI − F1_concat` |
| **Feature-access-matched transductive concat** | Transfer control matching unlabeled target **feature access** only |

These are descriptive contrasts, not causal architecture effects.

## Main findings (frozen)

### Development cohort

| Representation | Macro-F1 |
| --- | ---: |
| RNA PCA | 0.666 |
| Protein PCA | 0.512 |
| Concat PCA | 0.745 |
| MOFA+ | 0.549 |
| scVI (matched) | 0.722 |
| totalVI | 0.779 |

| Contrast | Point | 95% cell-level conditional paired-bootstrap | Notes |
| --- | ---: | --- | --- |
| `G_access` | **+0.079** | [0.043, 0.132] | prop>0 = 1.0 |
| `G_assoc` | **+0.034** | [−0.012, 0.073] | **crosses zero** |
| RNA→totalVI gap | **+0.113** | [0.084, 0.154] | |
| Descriptive recovery ≈0.70 | **~0.70** | [0.405, 1.114] | not causal / not variance explained |

### Transfer

Inductive concat 0.735 / 0.734; feature-access-matched transductive concat 0.636 / 0.767; scVI 0.656 / 0.601; totalVI 0.652 / 0.654 (Directions A/B). totalVI − feature-access concat ≈ +0.016 / −0.113. Direction B totalVI−scVI ≈ +0.052 (**19/20** seeds).

### Protein degradation

Training-time totalVI (5 mask/retraining replicates at nonzero levels; `model_seed=0`; clean *n*=1): 0.779 → 0.771±0.010 → 0.758±0.012 → 0.760±0.012 → 0.750±0.016 → 0.736±0.020 → 0.755±0.011. High-corruption rebound is within mask/retraining variability. Post-adaptation Direction B: 0.654 → 0.621 (clean scVI 0.601).

### Lawlor external validation (biological *n* = 10 donors)

- Donor-mean present-class macro-F1: scVI_matched **0.794**, totalVI **0.884**
- Mean paired Δ ≈ **+0.090**; **10/10** donors favor totalVI
- Exact two-sided donor-level sign test: **P=0.00195**
- Protein PCA ≈ **0.927**; concat ≈ **0.945** (protein-gated labels)
- Latent-AUROC computational consistency: **38/50** (not biological *n*; distinct from PBMC Direction B **19/20**)

Interpretation: supports cross-donor value of protein-aligned multimodal measurement information—not universal architecture superiority.

### Submission revision package

Frozen manuscript analyses, Methods expansion, and sensitivities live under [`results/manuscript_submission_revision/`](results/manuscript_submission_revision/):

- Canonical manuscript: [`manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md)
- Final DOCX: [`manuscript/CITEseq_Protein_Access_Model_Gains.docx`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Gains.docx)
- Audits, supplements, tables, and scripts: see [`REPRODUCIBILITY_README.md`](results/manuscript_submission_revision/REPRODUCIBILITY_README.md)

Supplementary sensitivities (do not replace primary analyses): class-stratified bootstrap; Direction B class-level degradation (CD8 Naive ΔF1≈−0.264); selective macro-F1 vs coverage; **PCA-only** marker-channel dropout (not a totalVI panel-failure test).

## Datasets

| Cohort | Role | Provenance |
| --- | --- | --- |
| PBMC10k / PBMC5k | Development / transfer | Public scvi-tools / 10x CITE-seq example AnnData files (no unverified GEO/SRA claimed) |
| Lawlor | External | HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005` |

## Repository structure (selected)

```
scripts/                         # Phase runners
src/                             # Shared library
data/                            # Processed objects / manifests
results/
  revision_round2_methodological_fixes/
  revision_round3_information_vs_model/
  presubmission_revision/        # Primary bootstrap + degradation audit
  manuscript_submission_revision/  # Frozen submission package
requirements.txt / environment.yml
```

## Key reproducibility scripts

| Analysis | Path |
| --- | --- |
| Development multimodal baselines | `results/revision_round2_methodological_fixes/scripts/revision2_simple_multimodal_baseline.py` |
| Feature-access-matched transfer | `results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py` |
| Primary paired bootstrap | `results/presubmission_revision/bootstrap/presubmission_core_contrast_bootstrap.py` |
| Class-stratified bootstrap | `results/manuscript_submission_revision/scripts/run_class_stratified_bootstrap.py` |
| Lawlor donor summaries | `results/manuscript_submission_revision/scripts/lawlor_donor_summary_stats.py` |
| Class-level degradation summary | `results/manuscript_submission_revision/scripts/summarize_class_degradation.py` |
| PCA-only marker-channel dropout | `results/manuscript_submission_revision/scripts/marker_channel_dropout_pca.py` |
| Environment printer | `scripts/print_environment.py` |
| Final DOCX builder | `results/manuscript_submission_revision/scripts/build_final_docx.py` |

## Environment

Pinned (see `environment.yml` / `requirements.txt`): Python **3.11.7**, scvi-tools **1.3.3**, scikit-learn **1.5.2**, PyTorch **2.14.0**, scanpy **1.10.3**, anndata **0.11.4**.

```bash
conda create -n multiomics_robustness python=3.11.7 pip -y
conda activate multiomics_robustness
pip install -r requirements.txt
```

## Reproducing the analyses

There is no single end-to-end one-command pipeline. Prefer frozen tables under `results/` and the submission-revision package. Full VAE stages require substantial compute (GPU/CUDA or MPS).

## Manuscript status

**Title:** Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift

**Status:** Pre-submission research manuscript (science frozen for this package)

**Author:** Jizhu Yang

## Data and code availability

Lawlor: HCA / ENA as above. PBMC: scvi-tools example CITE-seq files.  
Code: https://github.com/Jizhuy/robust-multiomics-integration.git  
An immutable archival release and DOI will accompany the submitted version.

## Citation / contact

Jizhuy — [yangjizhu@outlook.com](mailto:yangjizhu@outlook.com)
