# PHASE 12B — Data provenance audit

**Status:** audit only (no model training; no FASTQ bulk download).  
**Date:** 2026-09-10  
**Candidate:** Lawlor et al. PBMC CITE-seq (baseline / LPS / anti-CD3/CD28)

---

## Official sources (preferred order)

| Priority | Source | Accession / ID | Role |
|---|---|---|---|
| 1 | Publication | [doi:10.3389/fimmu.2021.636720](https://doi.org/10.3389/fimmu.2021.636720) — *Frontiers in Immunology* (2021) | Primary citation |
| 2 | Human Cell Atlas | Project `efea6426-510a-4b60-9a19-277e52bfa815` ([explorer](https://explore.data.humancellatlas.org/projects/efea6426-510a-4b60-9a19-277e52bfa815)) | **Preferred for processed matrices** |
| 3 | ENA / INSDC | `PRJEB40376` / `ERP124005` | Raw CITE-seq FASTQ (RNA, ADT, HTO) |
| 4 | ENA genotypes | `PRJEB40448` | Donor genotype VCF for Demuxlet |
| 5 | Author GitHub | [nlawlor/PBMC_CITEseq](https://github.com/nlawlor/PBMC_CITEseq) (MIT) | Shiny app code; points to HCA/ENA |
| 6 | Author Shiny | https://czi-pbmc-cite-seq.jax.org/ / https://thejacksonlaboratory.shinyapps.io/czi-pbmc-cite-seq | Visualization only |

**License / accessibility (HCA):** Creative Commons Attribution 4.0 (CC BY 4.0) under HCA Data Release Policy. Publicly accessible (“Access Granted” on HCA explorer).

**BioStudies accession listed by HCA:** `S-SUBS15`.

---

## Critical non-source (do not use as Lawlor)

| Source | Accession | Why rejected for PHASE 12B external validation |
|---|---|---|
| GEO | **GSE164378** | This is **Hao et al. 2021** (*Integrated analysis of multimodal single-cell data*), **not** Lawlor. It is the Seurat v4 / Azimuth multimodal PBMC reference underlying PHASE 6 label transfer. Using it as “external” would violate independence. |

Web search occasionally mis-associates GSE164378 with Lawlor; series title and contact (Yuhan Hao, NYGC) confirm it is Hao.

---

## Processed files (HCA contributor matrices)

| File | Size | Type | Raw vs processed |
|---|---|---|---|
| `CZI.PBMC.RNA.matrix.Rds` | ~481 MB | dgCMatrix 32,738 genes × 282,528 droplets | **Processed raw UMI counts** (integer) |
| `CZI.PBMC.ADT.matrix.Rds` | ~21 MB | data.frame 51 × 282,528 | **Processed raw ADT counts** (integer; 39 Abs + controls + QC rows) |
| `CZI.PBMC.HTO.matrix.Rds` | ~4.8 MB | data.frame 7 × 282,528 | HTO count matrix |
| `CZI.PBMC.cell.annotations.csv` | ~26 MB | CSV | Donor, condition, Demuxlet, author cell types |
| `CITEseqPBMCProject_metadata_18-07-2023.xlsx` | ~43 KB | XLSX | HCA ingest metadata |

**Estimated download for primary analysis (recommended):** ≈ **0.53 GB** (RNA+ADT+HTO+annotations+metadata).

**Raw FASTQ (ENA/HCA):** ≈ **2.5 TB** across 60 `fastq.gz` files — **do not download** for PHASE 12B unless processed matrices prove inadequate.

---

## Sequencing / processing (from paper Methods)

- 10 healthy donors (5F/5M, 21–32 y, Caucasian)
- Conditions: Baseline (medium), LPS, anti-CD3/CD28 (24 h culture)
- Multiplexing: Cell hashing (condition) + Demuxlet (donor genotypes)
- 10× Genomics Single Cell 3′ v2; ~25k cells/lane × 10 lanes
- RNA: Cell Ranger **2.1.0**, genome **hg19**
- ADT/HTO: CITE-seq-Count (`-cells 25000 -cbf 1 -cbl 16 -umif 17 -umil 26 --max-error 3`)
- Paper: 195,480 droplets → 55,353 Demuxlet singlets (mean ~1,932 UMI/cell) → further protein multiplet filtering

---

## Audit downloads performed (local Mac workspace)

Under `results/phase12b_external_validation/raw_audit/`:

- HCA annotations, ADT, HTO, metadata, RNA matrix (for gene-ID audit only)
- GEO GSE164378 **metadata only** (to prove it is Hao, then rejected)
- ENA filereport (sizes only; no FASTQ bodies)

No Cell Ranger reprocessing. No model training. No commits/pushes.
