# PHASE 10 — Cross-dataset transfer / query generalization

Platform: RTX 4060 Laptop GPU, WSL2 Ubuntu 26.04, torch 2.14.0+cu126, scvi-tools 1.3.3.
All 20 source models trained on `cuda:0`. Historical Mac/MPS artifacts untouched.

## Design

| | Direction A | Direction B |
|---|---|---|
| source | PBMC10k, 6,855 cells | PBMC5k, 3,994 cells |
| target | PBMC5k, 3,994 cells | PBMC10k, 6,855 cells |
| shared genes | 15,792 (identical order) | 15,792 (identical order) |
| shared proteins | 14 (identical order) | 14 (identical order) |
| disjoint cells | yes | yes |

Query API: `prepare_query_anndata` + `load_query_data` + `train` (scArches, `ArchesMixin`),
identical for `SCVI` and `TOTALVI`. This is **query-adapted transfer, not zero-shot**.

## Primary endpoint — high-confidence l2 macro-F1 (logreg, source-trained)

Mean ± SD over 5 model seeds.

| direction | scVI_matched | totalVI | Delta | seeds with Delta>0 |
|---|---|---|---|---|
| PBMC10k → PBMC5k | 0.6606 ± 0.0157 | 0.6504 ± 0.0126 | **−0.0102 ± 0.0164** | 2/5 |
| PBMC5k → PBMC10k | 0.6056 ± 0.0280 | 0.6483 ± 0.0335 | **+0.0427 ± 0.0431** | 4/5 |

Secondary l1:

| direction | scVI_matched | totalVI | Delta |
|---|---|---|---|
| PBMC10k → PBMC5k | 0.8205 ± 0.0016 | 0.8244 ± 0.0171 | +0.0038 ± 0.0164 |
| PBMC5k → PBMC10k | 0.7619 ± 0.0418 | 0.8595 ± 0.0428 | **+0.0977 ± 0.0321** (5/5 positive) |

The combined-dataset l2 advantage of +0.057 (PHASE 6) is **not reproduced** when PBMC10k is
the source. It is approximately reproduced (+0.043) when PBMC5k is the source.

## Cross-dataset neighbor label agreement

totalVI is higher in **all 12** direction × label-level × k combinations, with small
between-seed SD. This is the most consistent signal in PHASE 10 and it does not depend on
a classifier.

| direction | level | k=15 scVI → totalVI | k=50 scVI → totalVI |
|---|---|---|---|
| A | l1 | 0.8395 → 0.8495 | 0.8166 → 0.8247 |
| A | l2 | 0.6927 → 0.7092 | 0.6576 → 0.6705 |
| B | l1 | 0.7667 → 0.8108 | 0.7496 → 0.7856 |
| B | l2 | 0.6397 → 0.6789 | 0.6110 → 0.6437 |

## Latent coordinate system

Source and target cells are both encoded by the **post-adaptation** model. Measured drift of
the source embedding relative to the pre-adaptation source model:

- scVI_matched: 0.0000 in all 10 runs. With `use_observed_lib_size=True` there is no library
  encoder, and the only trainable weights are batch one-hot columns that receive exactly zero
  gradient for the source category because no query cell carries it.
- totalVI: mean 0.016–0.020 (0.3–0.4% of the mean latent norm), from the trainable library
  encoder and protein background priors.

Both are far below between-cell-type distances, so the cross-space neighbor metrics are valid.

## Generalization gap vs PHASE 6 (descriptive only; designs differ)

| direction | level | scVI_matched | totalVI |
|---|---|---|---|
| A | l2 | −0.061 | **−0.129** |
| A | l1 | −0.075 | −0.072 |
| B | l2 | −0.116 | −0.131 |
| B | l1 | **−0.134** | −0.037 |

totalVI's within-dataset l2 lead is the part that degrades most under dataset shift.
scVI_matched is the model that collapses when the source is small (direction B, l1).

## Dataset mixing (paired with conservation, per §27)

Dataset ASW (lower = better mixed): direction A 0.0291 (scVI) vs −0.0022 (totalVI);
direction B 0.0857 vs 0.0692. totalVI mixes source and target better **while also** having
higher neighbor agreement and higher target kNN purity, so this is not over-correction.

## Cell types

Consistently helped by protein in both directions (mean Delta F1 over 5 seeds):
CD8 Naive +0.254, CD8 TCM +0.127, CD8 TEM +0.102, CD4 Naive +0.085, B memory +0.055,
CD4 TEM +0.024, CD4 TCM +0.024, CD16 Mono +0.009, pDC +0.004.

Consistently harmed in both directions:
CD4 CTL −0.085, NK_CD56bright −0.075, Eryth −0.064, Treg −0.062, CD14 Mono −0.004.

Treg and NK_CD56bright are notable: these are populations where CD25 and CD56 surface protein
would be expected to help, and they get worse under transfer. Their support is small
(30–241 source cells), so between-seed noise is large.

Direction-dependent: MAIT (−0.056 A, +0.219 B), gdT (−0.130 A, +0.061 B),
ILC (−0.163 A, +0.014 B), B intermediate (+0.059 A, −0.108 B).

## Label harmonization

l1: 8 shared classes both directions, target coverage 1.000, nothing excluded.
l2: 29 intersecting classes; after the pre-specified ≥20 source / ≥10 target rule,
22 eligible in direction A (coverage 0.9932) and 21 in direction B (coverage 0.9774).
Excluded in both: ASDC, CD4/CD8/NK Proliferating, HSPC, Plasmablast, Platelet, cDC1
(plus dnT in direction B, 14 source cells). HSPC is source-only in A and target-only in B.

## Resources

| | source epochs | query epochs | source s | query s | peak VRAM alloc | peak VRAM reserved |
|---|---|---|---|---|---|---|
| scVI_matched | 200 | 197.5 | 90.2 | 72.5 | 406 MB | 420 MB |
| totalVI | 199 | 200 | 138.8 | 136.5 | 428 MB | 450 MB |

Total 2,290 s source training + 2,090 s query adaptation. Peak process RSS 6.78 GB;
minimum available system RAM 1.69 GB. batch_size stayed 256 throughout — no OOM, no fallback.

## No-leakage audit (all 20 manifests)

`target_labels_used_in_fitting: false`, `zero_shot: false`, `source_device`/`query_device`
= `cuda:0` in every manifest. Query early stopping monitored unlabeled `elbo_validation`.
