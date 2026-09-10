# PHASE 12B — External dataset feature audit (Lawlor)

## Required features checklist

| Feature | Available? | Exact columns / files |
|---|---|---|
| **A. Raw RNA UMI counts** | **YES** | `CZI.PBMC.RNA.matrix.Rds` — `dgCMatrix`, 32,738 × 282,528, integer UMIs |
| **B. Raw ADT counts** | **YES** | `CZI.PBMC.ADT.matrix.Rds` — integer counts; 39 Abs + isotype controls + `bad_struct`/`no_match`/`total_reads` |
| **C. Donor identity** | **YES** | `CZI.PBMC.cell.annotations.csv` → `Donor_of_Origin` (10 Illumina array IDs); Demuxlet `Demuxlet_Classification` ∈ {SNG, DBL, AMB} |
| **D. Experimental condition** | **YES** | `HTO_Classification` ∈ {Baseline, LPS, CD3_CD28, IgM_IgG, Multiplet, Empty}; also `HTO_Barcodes` |
| **E. Author cell-type annotation** | **PARTIAL** | `Celltype_Annotation` filled for **16,382 / 282,528** droplets (7 classes). Not all Demuxlet singlets are labeled in the deposited CSV |
| **F. RNA↔ADT barcode mapping** | **YES** | Shared cell barcodes as matrix column names / annotation `barcode` (e.g. `1-Sample_CGGAGCTAGTCGAGTG`) |
| **G. Batch / run metadata** | **YES** | `Run_Identifier` (10 lanes: `1-Sample` … `10-Sample`); all donors present in all lanes (pooled design) |

---

## Cell-type label audit (critical)

### Author labels present?

**YES**, but sparse in the deposited file.

| Label | n (all conditions) | n Baseline only |
|---|---|---|
| CD4T_Naive | 4778 | (see by-condition table) |
| CD4T_Mem | 3035 | |
| CD8T_Naive | 2648 | |
| CD8T_Mem | 1778 | |
| CD14_Mono | 1653 | |
| NK | 1409 | |
| B | 1081 | |
| **Total annotated** | **16382** | **5207** |

Full tables: `tables/external_label_distribution.csv`, `external_label_distribution_by_condition.csv`.

### Annotation levels

- **One deposited level** (7 coarse classes).
- Paper describes finer within-lineage clusters (e.g. CD4 clusters with CCR6/CXCR5/CD57) in figures — **those finer names are not columns in the HCA CSV**.
- **Do not treat Leiden clusters as validated types** without author labels.

### Label derivation

**Jointly informed but protein-primary:** paper Methods annotate using ADT markers (CD19/CD20 B; CD16/CD56 NK; CD14/CD11c Mono; CD3/CD4; CD3/CD8a; CD45RA/RO naive/memory). RNA used for clustering/trajectories after protein gating.

### Rare / missing vs PHASE 6

External lacks dedicated DC, CD16 Mono, gdT, MAIT, Treg, Platelet, HSPC, etc. Primary external validation should use a **coarse harmonized ontology** (`tables/external_label_harmonization.csv`).

### Annotation coverage gap

| Subset | n cells |
|---|---|
| Matrices (all droplets) | 282,528 |
| Paper Demuxlet singlets (reported) | 55,353 |
| Deposited SNG ∩ {Baseline,LPS,CD3_CD28} | **53,707** |
| SNG ∩ Baseline | **16,175** |
| Author `Celltype_Annotation` non-null | **16,382** |
| Annotated ∩ Baseline | **5,207** |

**Implication:** For supervised metrics, either (i) restrict to annotated cells, or (ii) reproduce author protein multiplet filters + gating on all SNG Baseline cells **without** Seurat-v4 transfer. Option (ii) needs a documented reimplementation; not done in this audit.

---

## Protein panel overlap (vs current 14)

See `tables/external_protein_overlap.csv`.

| Metric | Count |
|---|---|
| Exact name overlap | **10** |
| Alias-adjusted overlap | **12** |
| Current-only (no match) | **CD15, TIGIT** |
| External ADT biological Abs | **39** |

Aliases accepted (high confidence only): `PD-1` ↔ `CD279_PD1`; `CD127` ↔ `CD127_IL7Ra`.  
No silent ambiguous matches.

---

## Gene overlap

See `tables/external_gene_overlap.csv`.

| | Current (PBMC10k∩PBMC5k inner) | Lawlor |
|---|---|---|
| n genes | 15,792 | 32,738 |
| ID system | **gene symbols** | **Ensembl ENSG\* (hg19 / Cell Ranger 2.1.0)** |
| Exact intersection | **0** | — |
| Intersection fraction | **0** | — |

**Training blocker until resolved:** need an explicit hg19 Ensembl ↔ symbol map (and likely GRCh38 symbol harmonization to current 10x genes). Do **not** uppercase or rename blindly.

---

## Donor / condition structure

Primary table: `tables/external_donor_condition_counts.csv` (SNG ∩ 3 paper conditions).

- **10 donors**, all three conditions, cells distributed across **all 10 lanes** (fully pooled).
- **Do not mix** Baseline with LPS / CD3_CD28 in the primary external-validation analysis.
- Extra deposited class **`IgM_IgG`**: 22,448 SNG cells — **not** one of the three paper-emphasized arms; **exclude** from primary design until author clarification.

### HTO matrix naming caveat

Raw HTO features are named `Bcell`, `Mono`, `Tcell`, `control` (+ QC rows), while annotations use Baseline/LPS/CD3_CD28/IgM_IgG. Prefer **`HTO_Classification`** as demultiplexed condition labels; do not assume HTO feature names equal condition names without author documentation.

---

## Design candidates (feasibility only — not implemented)

### DESIGN A — donor-held-out transfer (recommended primary)

- Cohort: **Baseline only**, Demuxlet SNG (± annotated filter).
- Train: 7–8 donors; evaluate: 2–3 held-out donors.
- Optional second split for uncertainty calibration.
- **Memory-safe** if restricted to Baseline SNG (~16k) or annotated Baseline (~5.2k) before densifying.

### DESIGN B — leave-one-donor-out / grouped folds

- 10 folds × 2 models × seeds → **~10×** DESIGN A compute.
- Feasible scientifically; expensive on RTX 4060 + 7.6 GB RAM.
- Prefer DESIGN A first; LODO as sensitivity if A succeeds.

### Computational cost (order-of-magnitude, RTX 4060)

| Design | Cells (Baseline) | Models | Rough wall-time* |
|---|---|---|---|
| A annotated-only | ~5.2k | scVI+totalVI × few seeds | hours–1 day |
| A SNG Baseline | ~16k | same | ~1–3 days |
| B LODO annotated | 10× A | | multi-day–week |

\*Depends on epochs, genes after intersection, and whether RNA is subsetted before load.

---

## Secondary biological-shift (later only)

| Shift | Feasible? | Notes |
|---|---|---|
| Baseline → LPS | **YES** | Same donors; monocyte-focused activation |
| Baseline → CD3_CD28 | **YES** | Lymphocyte activation; monocyte depletion expected |
| Mixing conditions in one “query” | **NO for primary** | Confounds donor-shift with activation-shift |

Do **not** run stimulation analyses until Baseline donor-held-out completes.

---

## Selective-prediction novelty note (PHASE 12B constraint)

Class-aware selective prediction / imbalance-aware abstention are **literature baselines**, not claimed as novel methods. Contribution frame: application under multimodal integration + donor/dataset shift + training stochasticity.
