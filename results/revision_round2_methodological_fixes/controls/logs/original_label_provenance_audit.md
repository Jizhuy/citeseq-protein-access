# Original label provenance audit (PART B / GATE 2)

## Source of truth

Code, not manuscript prose:

- `src/annotation/seurat_v4.py::transfer_labels`
- `src/experiments/phase6_run.py`

## Reference

- Official scvi-tools loader: `scvi.data.pbmc_seurat_v4_cite_seq`
- Cached: `data/raw/pbmc_seurat_v4.h5ad`
- Stuart et al., Cell 2021 Seurat v4 PBMC CITE-seq reference
- Hierarchy: `celltype.l1` / `celltype.l2` (/l3 available on reference)

## Query mapping modalities

**RNA ONLY.**

Evidence:

1. HVGs selected on reference RNA among shared genes (`select_reference_hvgs`).
2. SCVI/SCANVI trained on gene counts only.
3. Query surgery: `prepare_query_anndata` + `load_query_data` on RNA subset;
   docstring explicitly: **"Surgery uses query RNA only."**
4. Soft max-class probabilities → predicted l2; l1 aggregated from l2 members.
5. Protein used later for **QC marker validation only**, not for assigning labels.

## Confidence

`choose_confidence_tiers`: high threshold **0.85** when q75 ≥ 0.85.
Frozen rule in later phases: `annotation_tier_l2 == "high"` ↔
`l2_confidence >= 0.85`.

## Does protein inform development labels?

**No** for the assignment step.

Implication for PART B scientific question:

The development totalVI vs scVI comparison was **already** evaluated under
RNA-derived reference labels. A separate “RNA-only remapping” would reinvent
the same pipeline. Inventing a different RNA mapper would change the question.

## Honest control implemented here

1. Freeze PHASE 6 labels as `rna_only_labels.csv` (protein-independent by construction).
2. Compare them to themselves / historical multimodal-*named* labels (identity).
3. Re-evaluate representations including concat PCA on Set A (RNA-only HC)
   and Set B (intersection with historical HC — identical here).
4. State clearly that label-provenance circularity with protein **does not**
   explain the development result.

## Lawlor note (out of scope for this development control)

Lawlor **author** labels are protein-gated by design. That is a separate
external-validation concern already addressed with author labels; it is not
solved by re-annotating PBMC development cells.
