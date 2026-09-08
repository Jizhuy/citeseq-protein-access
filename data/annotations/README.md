# External cell-type annotations

This directory is the only place validated biological labels should be added.

## Acceptable sources

1. Published annotation associated with the original PBMC10k / PBMC5k CITE-seq
   datasets, preferably.
2. An independent reference-based annotation procedure (for example a
   transferred atlas), documented with its source and date.

## Unacceptable sources

- Clusters from the model currently being evaluated (scVI, scVI_matched, totalVI)
- Arbitrary manual labeling of that model's own UMAP
- Unsupervised Leiden / k-means clusters treated as ground truth

Using a model's own clusters to score that model is circular evaluation.

## Expected table

Copy `pbmc_cell_annotations_template.csv` and fill one row per cell.

Required for the loader:

- `cell_id` must match `adata.obs_names` on the combined inner object
  (example: `PBMC10k_AAACCCAAGATTGTGA-1`)
- `cell_type` (or populate `cell_type_coarse` / `cell_type_fine` and map them)

Optional:

- `cell_type_coarse`
- `cell_type_fine`
- `annotation_source`
- `annotation_confidence`

Unmatched or duplicated cell IDs are reported and **not** silently dropped.

The PHASE 2 placeholder `obs["cell_type"] == "unknown"` is not a label.
