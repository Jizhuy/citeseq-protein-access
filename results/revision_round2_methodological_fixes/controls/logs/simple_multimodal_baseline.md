# Simple multimodal baseline (PART A)

## Predeclared primary representation

**10 RNA PCs + 10 protein PCs = 20-D concatenated representation**
(matched to scVI/totalVI latent dimensionality).

Secondary (optional): 20 RNA PCs + rank-limited protein PCs — not used for primary claims.

## Preprocessing

- RNA: counts → HVG(2000, seurat_v3) → normalize 1e4 → log1p → scale (clip 10) → PCA
- Protein: CLR → train/source z-score → PCA
- Concat: StandardScaler on concatenated PCs
- **All fits on TRAIN (development) or SOURCE (transfer/Lawlor) only**

## Development high-confidence l2 macro-F1 (logreg)

| Representation | macro-F1 |
|---|---:|
| RNA PCA (10) | 0.6665 |
| protein PCA (10) | 0.5120 |
| **concat PCA 10+10** | **0.7449** |
| MOFA+ | 0.5492 |
| scVI_matched | 0.7220 |
| totalVI | 0.7789 |

Delta(totalVI − scVI) = 0.0569
Delta(totalVI − concat) = 0.0340
Delta(concat − scVI) = 0.0229

Fraction of totalVI−scVI gain reproduced by concat:
0.403

## Transfer (inductive PCA vs transductive VAEs)

See `tables/simple_multimodal_transfer.csv`.

## Lawlor (source-fit PCA per fold)

See `tables/simple_multimodal_lawlor.csv` (deterministic; one value per fold).
