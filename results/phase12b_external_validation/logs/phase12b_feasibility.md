# PHASE 12B — Feasibility decision (STOP after audit)

**Date:** 2026-09-10  
**Scope:** Data / feasibility audit only. No scVI/totalVI training. No FASTQ bulk download. No commits/pushes. Historical PBMC10k/PBMC5k results untouched.

---

## A–U return summary

### A. Exact external dataset identified

**Lawlor et al.** PBMC CITE-seq (baseline + LPS + anti-CD3/CD28), HCA project **CITEseqPBMCProject**.

### B. Publication and accession

- **Publication:** Lawlor N et al. *Front. Immunol.* 2021. doi:[10.3389/fimmu.2021.636720](https://doi.org/10.3389/fimmu.2021.636720)
- **HCA project:** `efea6426-510a-4b60-9a19-277e52bfa815`
- **ENA raw:** `PRJEB40376` / `ERP124005`
- **ENA genotypes:** `PRJEB40448`
- **Not Lawlor:** GEO **GSE164378** = Hao Seurat-v4 reference (rejected)

### C. Number of donors

**10** (`Donor_of_Origin` / Demuxlet).

### D. Number of cells

| Definition | n |
|---|---|
| Droplets in deposited matrices | 282,528 |
| Paper Demuxlet singlets | 55,353 |
| SNG ∩ {Baseline, LPS, CD3_CD28} | **53,707** |
| SNG ∩ Baseline (primary) | **16,175** |
| Author-annotated (all cond.) | **16,382** |
| Author-annotated ∩ Baseline | **5,207** |

### E. Conditions

**Baseline**, **LPS**, **CD3_CD28** (paper). Deposited extra: **IgM_IgG** (exclude from primary until clarified). Also Multiplet/Empty.

### F. RNA raw counts available?

**YES** — integer UMI `dgCMatrix` in `CZI.PBMC.RNA.matrix.Rds`.

### G. ADT raw counts available?

**YES** — integer ADT matrix; 39 Abs + controls.

### H. Author cell-type labels available?

**YES (partial)** — 7 protein-gated classes in `Celltype_Annotation` for 16,382 cells; not full SNG coverage.

### I. Batch/donor metadata available?

**YES** — donor IDs, Demuxlet class, `Run_Identifier` (10 lanes), HTO condition.

### J. Exact 14-protein overlap count

- **Exact:** **10 / 14**
- **Alias-adjusted:** **12 / 14** (`PD-1`→`CD279_PD1`, `CD127`→`CD127_IL7Ra`)
- **Missing:** **CD15**, **TIGIT**

### K. Gene overlap

- Current: **15,792 symbols**
- External: **32,738 Ensembl (hg19)**
- **Exact intersection: 0** → mapping required before joint training/evaluation

### L. Label compatibility with current PBMC labels

**Coarse only.** Harmonize to ~B / Mono_CD14 / CD4_Naive / CD4_Memory / CD8_Naive / CD8_Memory / NK. Fine Seurat l2 (TCM/TEM/Treg/MAIT/DC/…) **not** force-matched. See `tables/external_label_harmonization.csv`.

### M. Independence from PHASE 6 annotation reference

**Independent** if using Lawlor HCA matrices + author labels.  
**Not independent** if using GSE164378 / Seurat-v4 / Azimuth to annotate Lawlor.

### N. Expected disk usage

| Item | Size |
|---|---|
| Recommended processed download | **~0.53 GB** |
| Converted sparse h5ad (est.) | **~1–2 GB** |
| Raw FASTQ (avoid) | **~2.5 TB** |
| Audit workspace already pulled | RNA+ADT+HTO+ann (~0.53 GB) under `raw_audit/` |

### O. Expected RAM requirements

| Approach | Fit on ~7.6 GB WSL? |
|---|---|
| Load full 282k×32k dense | **NO** |
| Load full RNA Rds then subset (R peak observed ~5 GB during audit) | Marginal / risky |
| Subset barcodes **before** materializing dense arrays; keep sparse; Baseline SNG ~16k or annotated ~5.2k | **YES (recommended)** |
| Train scVI/totalVI on ~5–16k cells × ≤~10–15k mapped genes, batch 128–256 | **YES** on RTX 4060 8 GB with care |

### P. Recommended donor split

**DESIGN A:** Baseline-only; Demuxlet SNG; hold out **2–3 donors** (stratify by sex if recoverable from metadata; else random with fixed seed). Train on remaining **7–8**. Keep stimulated conditions sealed.

Prefer **annotated Baseline (5,207)** for label-based metrics in v1; optionally expand to all SNG Baseline with documented protein re-gating later.

### Q. Recommended primary external-validation design

**Donor-held-out Baseline transfer** of matched scVI vs totalVI (same recipe as PHASE 10/11 family), evaluating:

1. Fine-vs-coarse biological representation (on harmonized labels)
2. Conditional multimodal benefit under donor shift
3. Uncertainty → error / instability / abstention (findings 3–6)
4. Optional small seed panel for run-level variability (finding 7) — light, not 100-run scale

**Do not** claim novelty for class-aware selective prediction; treat as baselines.

### R. Stimulation analysis feasible?

**YES, later.** Baseline→LPS and Baseline→CD3_CD28 are well powered at SNG scale. **Not now.**

### S. Major limitations

1. Gene ID mismatch (Ensembl hg19 vs symbols) — hard prerequisite
2. Author labels cover only ~5.2k Baseline cells
3. Coarse ontology vs PHASE 6 l2
4. Missing CD15 & TIGIT vs current 14-panel
5. Fully pooled lanes → weak technical batch; donor is the shift axis
6. Unexplained `IgM_IgG` HTO class
7. HTO feature names ≠ condition names in matrix
8. RTX 4060 / 7.6 GB RAM forbids careless full-matrix loads

### T. Exact files that should be downloaded (for training prep, after approval)

**Download / retain:**

1. `CZI.PBMC.RNA.matrix.Rds`
2. `CZI.PBMC.ADT.matrix.Rds`
3. `CZI.PBMC.HTO.matrix.Rds` (optional if using annotation conditions only)
4. `CZI.PBMC.cell.annotations.csv`
5. `CITEseqPBMCProject_metadata_18-07-2023.xlsx`

**Do not download:** ENA/HCA FASTQ (~2.5 TB), GSE164378 as external cohort.

*(Audit already pulled 1–5 into `raw_audit/` on the Mac workspace; mirror to the 4060 server only after approval.)*

### U. Whether PHASE 12B training should proceed

**CONDITIONAL GO — wait for explicit approval.**

Proceed to **data prep + training** only after approving:

1. Use Lawlor HCA processed matrices (not GSE164378, not FASTQ)
2. Primary cohort = **Baseline SNG**, metrics on **author-annotated** subset (or approve protein re-gating plan)
3. Build **Ensembl↔symbol** intersection pipeline (no blind rename)
4. Protein panel = **12 shared** (drop CD15, TIGIT) or document missing-protein handling in totalVI
5. DESIGN A donor hold-out; stimulation deferred
6. Memory-safe sparse conversion on the 4060 host

**STOP here.** No training started.

---

## Findings 1–7 (replication targets once training approved)

The external cohort is suitable to test whether PHASE 6–12A conclusions generalize under **donor shift within an independent CITE-seq study**, with the caveats above. Expect partial, not universal, replication.

---

## Output index

```
results/phase12b_external_validation/
  logs/data_provenance.md
  logs/external_dataset_audit.md
  logs/independence_audit.md
  logs/phase12b_feasibility.md
  tables/external_file_inventory.csv
  tables/external_donor_condition_counts.csv
  tables/external_protein_overlap.csv
  tables/external_gene_overlap.csv
  tables/external_label_distribution.csv
  tables/external_label_distribution_by_condition.csv
  tables/external_label_harmonization.csv
  raw_audit/   # manifests + processed matrices used for audit
```
