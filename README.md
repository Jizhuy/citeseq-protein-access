# CITE-seq Protein Access and Representation Learning

**Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift**

**Author:** Jizhu Yang

**Repository:** [https://github.com/Jizhuy/citeseq-protein-access](https://github.com/Jizhuy/citeseq-protein-access)

---

## Clone this repository

HTTPS:

```bash
git clone https://github.com/Jizhuy/citeseq-protein-access.git
cd citeseq-protein-access
```

SSH:

```bash
git clone git@github.com:Jizhuy/citeseq-protein-access.git
cd citeseq-protein-access
```

---

## Overview

Paired CITE-seq measurements of RNA and surface protein can improve single-cell representation learning, but it is often unclear how much of the gain comes from **protein access** versus a more complex multimodal model.

This repository accompanies a research manuscript that separates:

1. **Protein-access gain** — improvement from adding paired protein features to a simple multimodal baseline  
2. **Residual totalVI-associated contrast** — any remaining difference between totalVI and that simple RNA–protein baseline  

Analyses stress-test these contrasts under:

- strict train-only internal evaluation  
- bidirectional transfer with feature-access-matched controls  
- protein degradation  
- uncertainty diagnostics and selective prediction  
- training stochasticity  
- Lawlor donor-level external validation  

This is an empirical reliability study of existing representations (PCA baselines, MOFA+, scVI_matched, totalVI), **not** a new architecture paper.

---

## Primary scientific finding

Under a **strict train-only** internal protocol (7,573 high-confidence training cells / 1,921 held-out cells):

| Representation | Macro-F1 |
| --- | ---: |
| RNA PCA | 0.667 |
| Concat RNA–protein PCA | 0.745 |
| scVI_matched | 0.709 |
| totalVI | 0.741 |

| Contrast | Estimate | 95% conditional paired-bootstrap CI |
| --- | ---: | --- |
| **G_access** (concat − RNA PCA) | **+0.078** | [0.042, 0.131] |
| **G_assoc** (totalVI − concat) | **−0.004** | [−0.060, 0.030] |
| **Gap_total** (totalVI − RNA PCA) | **+0.074** | [0.029, 0.130] |

**Interpretation.** Protein-access gain is substantial and resampling-stable. After matching information access, the residual totalVI-associated contrast is **near null and imprecisely estimated** (interval crosses zero). The descriptive ratio `R_access ≈ 1.05` is reported only as a secondary quantity and is **not** a bounded recovery fraction.

Secondary / historical: protein PCA = 0.512; MOFA+ = 0.549 (not re-fit under strict train-only).

---

## Analysis components

| Component | Role |
| --- | --- |
| RNA PCA / protein PCA / concat PCA | Simple baselines (train-only) |
| scVI_matched | Capacity-matched RNA-only deep representation (train-only) |
| totalVI | Multimodal deep representation (train-only) |
| Transfer + feature-access-matched concat | Distribution-shift stress test |
| Protein degradation | Antibody-signal dependence |
| Uncertainty / selective prediction | Reliability diagnostics |
| Training stochasticity | Precision of small contrasts |
| Lawlor donor validation | Biological donor unit (*n* = 10) |

### Lawlor donor-level summary (unchanged)

| Representation | Donor-mean present-class macro-F1 |
| --- | ---: |
| RNA PCA | 0.690 |
| Protein PCA | 0.923 |
| Concat PCA | 0.942 |
| scVI_matched | 0.794 |
| totalVI | 0.884 |

Paired totalVI − scVI_matched: mean Δ ≈ +0.090; 10/10 donors; exact two-sided sign test *P* = 0.00195. Protein-gated labels favor protein-aligned baselines; this supports multimodal measurement value, not universal totalVI superiority.

---

## Canonical manuscript package

All primary figures, tables, supplements, and audits for the current manuscript live under:

[`results/manuscript_submission_revision/`](results/manuscript_submission_revision/)

| Item | Path |
| --- | --- |
| Manuscript (Markdown) | [`manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md) |
| Manuscript (DOCX) | [`manuscript/CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx`](results/manuscript_submission_revision/manuscript/CITEseq_Protein_Access_Model_Gains_Submission_Ready.docx) |
| Strict train-only results | [`strict_train_only/`](results/manuscript_submission_revision/strict_train_only/) |
| Supplementary methods / results | [`supplement/`](results/manuscript_submission_revision/supplement/) |
| Reproducibility notes | [`REPRODUCIBILITY_README.md`](results/manuscript_submission_revision/REPRODUCIBILITY_README.md) |
| Internal audit | [`STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md`](results/manuscript_submission_revision/STRICT_TRAIN_ONLY_INTERNAL_AUDIT.md) |

---

## Repository layout (selected)

```
scripts/                              # Phase runners
src/                                  # Shared library code
data/                                 # Processed objects / manifests (large files not versioned)
results/
  manuscript_submission_revision/     # PRIMARY submission package
  revision_round3_information_vs_model/
  revision_round2_methodological_fixes/
  presubmission_revision/             # Historical bootstrap / degradation
environment.yml
requirements.txt
README.md
```

---

## Environment

Pinned (see `environment.yml` / `requirements.txt`):

- Python **3.11.7**
- scvi-tools **1.3.3**
- scikit-learn **1.5.2**
- PyTorch **2.14.0**
- scanpy **1.10.3**
- anndata **0.11.4**

```bash
conda create -n multiomics_robustness python=3.11.7 pip -y
conda activate multiomics_robustness
pip install -r requirements.txt
```

---

## Key scripts

| Analysis | Script |
| --- | --- |
| Strict train-only internal (**PRIMARY**) | `results/manuscript_submission_revision/scripts/run_strict_train_only_internal.py` |
| Feature-access-matched transfer | `results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py` |
| Class-stratified bootstrap | `results/manuscript_submission_revision/scripts/run_class_stratified_bootstrap.py` |
| Lawlor donor summaries | `results/manuscript_submission_revision/scripts/lawlor_donor_summary_stats.py` |
| Final DOCX builder | `results/manuscript_submission_revision/scripts/build_final_docx.py` |
| Environment printer | `scripts/print_environment.py` |

There is no single one-command end-to-end pipeline. Prefer frozen tables under `results/`. Full VAE training requires substantial compute (GPU/CUDA or MPS).

---

## Datasets

| Cohort | Role | Provenance |
| --- | --- | --- |
| PBMC10k / PBMC5k | Development / transfer | Public scvi-tools / 10x CITE-seq example AnnData (no unverified GEO/SRA claimed) |
| Lawlor blood CITE-seq | External donor validation | HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005` |

---

## Data and code availability

- **Code:** https://github.com/Jizhuy/citeseq-protein-access  
- **Lawlor:** HCA / ENA accessions above  
- **PBMC:** scvi-tools example CITE-seq files  

An immutable archival release and DOI will accompany the submitted version of the manuscript.

---

## Citation

Jizhu Yang. *Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift.* Pre-submission research manuscript.  
Repository: https://github.com/Jizhuy/citeseq-protein-access

**Contact:** [yangjizhu@outlook.com](mailto:yangjizhu@outlook.com)
