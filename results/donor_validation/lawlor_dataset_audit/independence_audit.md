# PHASE 12B — Independence audit

**Question:** Is Lawlor et al. independent of PHASE 6 annotation transfer and of PBMC10k/PBMC5k processing?

---

## Verdict

**Lawlor HCA project `efea6426-…` is biologically and accession-independent of PHASE 6’s Seurat v4 reference.**

**GEO GSE164378 is NOT independent** — it *is* the Hao multimodal reference. It must not be used as the PHASE 12B external cohort.

---

## PHASE 6 annotation dependency chain

| Step | What we use | Source |
|---|---|---|
| Reference object | `scvi.data.pbmc_seurat_v4_cite_seq()` | Hao et al. 2021 multimodal PBMC CITE-seq atlas |
| GEO twin | GSE164378 | Same Hao study (confirmed via SOFT `!Series_title`, pubmed 34062119, contact Yuhan Hao) |
| Labels transferred | `celltype.l1` / `celltype.l2` / `celltype.l3` | Seurat v4 WNN ontology |
| Query | PBMC10k + PBMC5k (10x Genomics public CITE-seq) | Independent commercial datasets |

Code: `src/annotation/seurat_v4.py` → `download_seurat_v4_reference()`.

---

## Lawlor vs Hao

| Property | Lawlor (candidate external) | Hao GSE164378 (PHASE 6 ref) |
|---|---|---|
| Publication | Front. Immunol. 2021; 10.3389/fimmu.2021.636720 | Nat. Biotechnol. 2021; PMID 34062119 |
| Accessions | HCA efea6426…; ENA PRJEB40376 / ERP124005 | GEO GSE164378; BioProject PRJNA690251 |
| Donors | 10 healthy adults | 8 donors (P1–P8) + vaccination timepoints |
| Conditions | Baseline, LPS, anti-CD3/CD28 (+ deposited `IgM_IgG` class) | Vaccination / time-course design |
| ADT panel | ~39 Abs | 228 Abs (3′) / 54 Abs (5′) |
| Author labels | Coarse protein-gated (7 classes) | Fine Seurat l1/l2/l3 |
| Used in PHASE 6? | **No** | **Yes** (reference) |

---

## PBMC10k / PBMC5k processing

No Lawlor files appear in:

- `data/processed/pbmc*_cite*.h5ad`
- PHASE 2–5 loaders
- PHASE 6–12A training paths

PBMC10k/PBMC5k remain 10x Genomics products with Seurat-v4–transferred labels only.

---

## Label ontology caution (independence of *labels*, not of *data*)

Author Lawlor labels are **not** Seurat l1/l2 strings, but they were derived from overlapping lineage markers (CD3/CD4/CD8/CD14/…). That is expected for PBMC and does **not** create a data-dependence on Hao.

**Do not** re-annotate Lawlor with Seurat v4 / Azimuth for the primary external validation — that would reintroduce the PHASE 6 reference into the evaluation target.

---

## Residual risks

1. **Shared technology stack** (10x 3′, TotalSeq, BioLegend) — expected; not an independence failure.
2. **NYGC coauthors** (Stoeckius, Smibert) appear on both CITE-seq ecosystem papers — does not imply shared cells.
3. **Mistaken GEO download** of GSE164378 would silently destroy independence — blocked by this audit.

---

## Conclusion

Proceed with **Lawlor HCA processed matrices** as the independent external cohort.  
**Reject GSE164378** for PHASE 12B validation use.
