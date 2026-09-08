# Data directory

Validated PHASE 2 contents (this machine):

| Object | Cells | Genes | Proteins | RNA sparsity | Protein sparsity |
| --- | ---: | ---: | ---: | ---: | ---: |
| PBMC10k | 6855 | 16727 | 14 | 0.9100 | 0.0018 |
| PBMC5k | 3994 | 16581 | 29 | 0.8794 | 0.0664 |
| Combined inner | 10849 | 15792 | 14 shared | 0.8933 | 0.0136 |
| Combined outer | 10849 | 15792 | 29 union | 0.8933 | 0.3518* |

\*Outer-join protein sparsity includes unmeasured PBMC10k×PBMC5k-only entries treated as missing (102,825 NaNs; 6,855 cells × 15 proteins). Do not interpret that number as biological zero-inflation.

There are **no curated cell-type labels** in the official files. `obs["cell_type"]` is an explicit `unknown` placeholder.

The official scverse download URL returned HTTP 403. Files were retrieved from the YosefLab/scVI-data GitHub fallback and match the published scvi-tools tutorial dimensions (10,849 × 15,792 combined inner join).

## `raw/`

Official scvi-tools CITE-seq example files, downloaded unchanged:

- `pbmc_10k_protein_v3.h5ad`
- `pbmc_5k_protein_v3.h5ad`

Source URLs (same as `scvi.data.pbmcs_10x_cite_seq`):

- https://exampledata.scverse.org/scvi-tools/pbmc_10k_protein_v3.h5ad
- https://exampledata.scverse.org/scvi-tools/pbmc_5k_protein_v3.h5ad

These files were previously filtered for doublets and low-quality cells following the totalVI reproducibility scripts. Phase 2 does not apply additional QC filters.

## `processed/`

Created by `python scripts/01_prepare_data.py`.

| File | Contents |
| --- | --- |
| `pbmc10k_cite.h5ad` | PBMC10k only, all measured proteins |
| `pbmc5k_cite.h5ad` | PBMC5k only, all measured proteins |
| `pbmc_cite_combined_inner.h5ad` | Shared genes + shared proteins (dense float32; not overwritten after PHASE 2) |
| `pbmc_cite_combined_inner_sparse.h5ad` | PHASE 3 CSR training cache; values verified identical to the dense inner object |
| `pbmc_cite_combined_outer.h5ad` | Shared genes + union of proteins |

### Count integrity

- `adata.X` is raw RNA counts.
- `adata.layers["counts"]` is a copy of those raw RNA counts.
- Normalized values are **not** written into `X`.
- Protein raw counts are in `adata.obsm["protein_counts"]`.
- `adata.obsm["protein_expression"]` is a scvi-tools-compatible copy of the protein matrix.
- `adata.obsm["protein_observed_mask"]` is True where a protein was actually measured.

### Protein panel join

PBMC10k and PBMC5k do not have identical ADT panels. This repository does not silently drop unique proteins:

- Inner join: analysis object using only shared proteins. Dropped proteins are listed in `results/tables/protein_panel_overlap.csv`.
- Outer join: union of proteins. Proteins absent from a dataset remain **NaN**, not zero.

The default modeling object for later phases is the inner-join combined dataset, matching the totalVI PBMC10k+PBMC5k tutorial, with the panel difference documented rather than hidden.

### Cell-type labels

The official scvi-tools PBMC10k/PBMC5k CITE-seq files are not guaranteed to contain curated cell-type annotations. Phase 2 records whatever is actually present in `obs` and writes `unknown` only as an explicit placeholder when no annotation column exists. Do not treat `unknown` as a biological label.
