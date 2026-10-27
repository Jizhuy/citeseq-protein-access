# Data availability

All datasets analysed in this study are publicly available. This repository does **not** redistribute raw or processed single-cell data. It contains **no private or controlled-access data**.

The repository's MIT License covers its code and repository-authored materials only. The public datasets listed below are not relicensed by this repository and remain subject to their original source terms (10x Genomics / scvi-tools example data, the Seurat v4 reference, and the Human Cell Atlas / ENA release of the Lawlor et al. data).

## 1. PBMC10k and PBMC5k CITE-seq (development cohort and bidirectional transfer)

| Item | Detail |
| --- | --- |
| Origin | Public 10x Genomics PBMC CITE-seq example data, as distributed by scvi-tools |
| Accession | None claimed. These are public example datasets; no GEO/SRA accession is asserted. |
| Expected source files | `pbmc_10k_protein_v3.h5ad`, `pbmc_5k_protein_v3.h5ad` |
| Download | URLs in `config/config.yaml` (`data.sources`): `https://exampledata.scverse.org/scvi-tools/…` with a GitHub `YosefLab/scVI-data` fallback. `scripts/prepare_data.py` downloads them if they are missing. |
| Local placement | `data/raw/` |
| Derived files (not redistributed) | `data/processed/pbmc10k_cite.h5ad`, `pbmc5k_cite.h5ad`, `pbmc_cite_combined_inner.h5ad`, `pbmc_cite_combined_outer.h5ad`, `pbmc_cite_combined_inner_annotated.h5ad` |
| Shipped derived data | `data/annotations/pbmc_reference_annotations.csv`: reference-transfer cell-type labels (no expression values) |

## 2. Annotation reference (label transfer only)

| Item | Detail |
| --- | --- |
| Origin | Seurat v4 multimodal PBMC CITE-seq reference (Hao et al., 2021), loaded via `scvi.data.pbmc_seurat_v4_cite_seq` |
| Role | Label-transfer reference for the PBMC development cohort only. It is **not** an evaluation cohort and is independent of the Lawlor validation cohort (see `results/donor_validation/lawlor_dataset_audit/independence_audit.md`). |
| Local placement | `data/raw/` (downloaded automatically by `scripts/run_reference_annotation.py`) |

## 3. Lawlor et al. PBMC CITE-seq (donor-level external validation)

| Item | Detail |
| --- | --- |
| Origin | Lawlor et al. PBMC CITE-seq study (10 healthy donors; Baseline, LPS and anti-CD3/CD28 conditions). Only the **Baseline** condition and Demuxlet singlets are used. |
| Accessions | Human Cell Atlas project `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005` (raw FASTQ, not used) |
| Expected source files (HCA contributor matrices) | `CZI.PBMC.RNA.matrix.Rds`, `CZI.PBMC.ADT.matrix.Rds`, `CZI.PBMC.HTO.matrix.Rds`, `CZI.PBMC.cell.annotations.csv` (~0.53 GB in total) |
| Local placement | `data/external/lawlor/`. Run `scripts/prepare_lawlor_baseline_matrices.R` from that directory to write `baseline_mm/`. |
| Derived file (not redistributed) | `results/donor_validation/lawlor_runs/processed/lawlor_baseline_sng_shared.h5ad` |
| Provenance records (shipped) | `results/donor_validation/lawlor_dataset_audit/data_provenance.md`, `independence_audit.md`, `results/donor_validation/lawlor_runs/tables/lawlor_download_checksums.csv`, `results/donor_validation/lawlor_runs/logs/data_prep_provenance.md` |

## 4. What is and is not in this repository

| Included | Not included |
| --- | --- |
| Frozen result tables (CSV/JSON/NPZ summaries), figures, supplementary material | Raw count matrices, FASTQ, `.h5ad` / `.h5` objects |
| Reference-transfer cell labels for PBMC development cells | Trained model checkpoints (`*.pt`, `*.ckpt`, `models/`) |
| Provenance and checksum records for external data | Any private, controlled-access or personally identifying data |

The frozen tables in `results/` and `supplementary/Supplementary_Data.zip` are enough to check every reported number without downloading any dataset.
