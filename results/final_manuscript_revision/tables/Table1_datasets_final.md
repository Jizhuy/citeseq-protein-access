# Table 1. Datasets and study designs

| Dataset | Role | Cells | Genes | Proteins | Labels | Primary experiments |
|---|---|---|---|---|---|---|
| PBMC10k CITE-seq | Development / transfer | 6,855 | 16,727 | 14 | Seurat-v4–derived (dev only) | Representation; corruption; transfer; uncertainty |
| PBMC5k CITE-seq | Development / reciprocal transfer | 3,994 | 16,581 | 29 (14 shared combined) | Seurat-v4–derived (dev only) | Transfer A/B; uncertainty |
| Lawlor Baseline CITE-seq | Independent external validation | 16,175 singlets (5,207 labeled) | 12,776 shared | 12 matched | Author labels only | 5-fold donor-held-out × 10 seeds (100 runs) |

Footnote: Combined development object 10,849 cells × 15,792 genes × 14 proteins; 9,494 high-confidence cells (confidence ≥0.85).
