# Data

Raw and processed single-cell data are **not redistributed**. See [`../docs/DATA_AVAILABILITY.md`](../docs/DATA_AVAILABILITY.md) for sources and accessions.

Expected local layout (git-ignored):

```
data/
├── raw/
│   ├── pbmc_10k_protein_v3.h5ad
│   └── pbmc_5k_protein_v3.h5ad
├── processed/              # written by scripts/prepare_data.py and scripts/run_reference_annotation.py
└── external/
    └── lawlor/             # Lawlor HCA contributor matrices + baseline_mm/ Matrix-Market export
```

`annotations/` is shipped. It contains the reference-transfer cell-type labels for the PBMC development cells (derived data, no raw counts) and a template, as written by `scripts/run_reference_annotation.py`.
