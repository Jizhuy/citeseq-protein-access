# CITE-seq Protein Access and Representation Learning

**Manuscript title:** Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift

Repository: [https://github.com/Jizhuy/citeseq-protein-access](https://github.com/Jizhuy/citeseq-protein-access)

```bash
git clone git@github.com:Jizhuy/citeseq-protein-access.git
cd citeseq-protein-access
```

## Overview

This project evaluates how much of the observed performance improvement in CITE-seq representation learning is associated with **access to paired protein measurements**, and how much remains as a **residual totalVI-associated contrast** beyond a simple RNA+protein representation—especially under distribution shift, protein degradation, and external donor validation.

It is an empirical reliability study of existing representations (PCA baselines, MOFA+, scVI_matched, totalVI), not a new architecture.

**Primary internal protocol:** strict train-only representation learning on 7,573 high-confidence training cells, with evaluation on 1,921 high-confidence held-out cells for RNA PCA, protein PCA, concat PCA, scVI_matched, and totalVI.

## Scientific question

Central comparison (development cohort, held-out macro-F1):

**RNA PCA → RNA + protein concat PCA → totalVI**

| Contrast | Definition |
| --- | --- |
| **Protein-access contrast** `G_access` | `F1_concat − F1_RNA_PCA` |
| **Residual totalVI-associated contrast** `G_assoc` | `F1_totalVI − F1_concat` |
| **Feature-access-matched transductive concat** | Transfer control matching unlabeled target **feature access** only |

These are descriptive contrasts, not causal architecture effects.

## Main findings (primary)

### Development cohort (strict train-only)

| Representation | Macro-F1 |
| --- | ---: |
| RNA PCA | 0.667 |
| Protein PCA | 0.512 |
| Concat PCA | 0.745 |
| scVI_matched | 0.709 |
| totalVI | 0.741 |
| MOFA+ (historical secondary) | 0.549 |

| Contrast | Point | 95% cell-level conditional paired-bootstrap | Notes |
| --- | ---: | --- | --- |
| `G_access` | **+0.078** | [0.042, 0.131] | prop>0 = 1.0; resampling-stable |
| `G_assoc` | **−0.004** | [−0.060, 0.030] | near null; **crosses zero** (prop>0 = 0.355) |
| `Gap_total` | **+0.074** | [0.029, 0.130] | RNA→totalVI gap |
| `R_access` | ≈1.05 | [0.680, 2.304] | secondary descriptive ratio only (not a recovery fraction) |

### Transfer / information-access controls

Inductive concat 0.735 / 0.734; feature-access-matched transductive concat 0.636 / 0.767; scVI_matched 0.656 / 0.601; totalVI 0.652 / 0.654 (Directions A/B). totalVI − feature-access concat ≈ +0.016 / −0.113. Direction B totalVI−scVI_matched ≈ +0.052 (**19/20** seeds).

### Protein degradation, uncertainty, stochasticity

Training-time totalVI degradation, uncertainty diagnostics, selective prediction, and training-stochasticity analyses are reported in the manuscript package. These stress-test residual contrasts under protein corruption and run-level variability.

### Lawlor external validation (biological *n* = 10 donors)

Donor-mean present-class macro-F1:

| Representation | Mean |
| --- | ---: |
| RNA PCA | 0.690 |
| Protein PCA | 0.923 |
| Concat PCA | 0.942 |
| scVI_matched | 0.794 |
| totalVI | 0.884 |

- Mean paired totalVI − scVI_matched Δ ≈ **+0.090**; **10/10** donors favor totalVI
- Exact two-sided donor-level sign test: **P=0.00195**

Interpretation: supports cross-donor value of protein-aligned multimodal measurement information under protein-gated labels—not universal architecture superiority.

### Submission revision package

Frozen manuscript analyses, Methods expansion, and sensitivities live under [`results/manuscript_submission_revision/`](results/manuscript_submission_revision/):

- Canonical manuscript (MD): [`manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md)
- Canonical DOCX: [`manuscript/CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx)
- Strict train-only artifacts: [`strict_train_only/`](results/manuscript_submission_revision/strict_train_only/)
- Audits and supplements: see [`REPRODUCIBILITY_README.md`](results/manuscript_submission_revision/REPRODUCIBILITY_README.md)

Supplementary sensitivities (do not replace primary analyses): information-access (historical feature-transductive) comparison; class-stratified bootstrap; Direction B class-level degradation; selective macro-F1 vs coverage; **PCA-only** marker-channel dropout (not a totalVI panel-failure test).

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
  presubmission_revision/        # Historical bootstrap + degradation audit
  manuscript_submission_revision/  # Frozen submission package (PRIMARY)
requirements.txt / environment.yml
```

## Key reproducibility scripts

| Analysis | Path |
| --- | --- |
| Strict train-only internal (PRIMARY) | `results/manuscript_submission_revision/scripts/run_strict_train_only_internal.py` |
| Development multimodal baselines (historical) | `results/revision_round2_methodological_fixes/scripts/revision2_simple_multimodal_baseline.py` |
| Feature-access-matched transfer | `results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py` |
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

**Title:** Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift

**Status:** Pre-submission research manuscript

**Author:** Jizhu Yang

## Data and code availability

Lawlor: HCA / ENA as above. PBMC: scvi-tools example CITE-seq files.  
Code: https://github.com/Jizhuy/citeseq-protein-access.git  
An immutable archival release and DOI will accompany the submitted version.

## Citation / contact

Jizhuy — [yangjizhu@outlook.com](mailto:yangjizhu@outlook.com)
