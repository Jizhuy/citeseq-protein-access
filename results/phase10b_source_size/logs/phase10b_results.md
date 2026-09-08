# PHASE 10B — source-size sensitivity control

Aggregation only. No model was retrained for this report.

Scope: Direction A (PBMC10k → PBMC5k) with the PBMC10k source randomly
subsampled from 6,855 to 3,994 cells, the PBMC5k source size. Direction B was
not retrained; its PHASE 10 embeddings were re-scored on the shared class set.

## Completion

10/10 manifests, `model_seed = 0` throughout, `source_n = 3994` and
`target_n = 3994` in every run, `source_device = query_device = cuda:0`,
`batch_size = 256` everywhere, `target_labels_used_in_fitting = false`,
`zero_shot = false`, `phase10_artifact_overwritten = false`. Zero failures.
Source drift under adaptation: 0.000000 for scVI_matched in all 5 runs,
mean 0.0181 (max 0.0198) for totalVI — same alignment behaviour as PHASE 10.

## Common evaluation class set

19 l2 classes, 8 l1 classes, identical across all three conditions.
Built as the intersection of the 5 subsample-eligible sets, the original
Direction A eligible set (22) and the original Direction B eligible set (21),
under the pre-specified ≥20 source / ≥10 target rule.

Dropped relative to the union: **CD4 CTL** and **pDC** (both fell below 20
source cells in all five 3,994-cell subsets) and **dnT** (already ineligible in
Direction B). Target coverage on the common l2 set is 0.9702 for both Direction
A conditions and 0.9695 for Direction B; l1 coverage is 1.000 everywhere.

## Primary endpoint — common l2 logreg macro-F1

| subsample_seed | scVI_matched | totalVI | Delta | bootstrap 95% CI | verdict |
|---|---|---|---|---|---|
| 0 | 0.6512 | 0.6291 | −0.0221 | [−0.0462, +0.0018] | crosses zero |
| 1 | 0.6255 | 0.6444 | +0.0189 | [−0.0017, +0.0382] | crosses zero |
| 2 | 0.6419 | 0.6504 | +0.0084 | [−0.0172, +0.0313] | crosses zero |
| 3 | 0.6260 | 0.6738 | +0.0478 | [+0.0258, +0.0730] | entirely positive |
| 4 | 0.6159 | 0.6148 | −0.0011 | [−0.0238, +0.0199] | crosses zero |

scVI_matched 0.6321 ± 0.0142, totalVI 0.6425 ± 0.0223,
**Delta +0.0104 ± 0.0258** (range −0.0221 to +0.0478).
The SD is across source subsets, not model seeds; bootstrap CIs are per-run and
are never pooled with it.

## Three-condition comparison (all on the same 19 common l2 classes)

| condition | source n | scVI_matched | totalVI | Delta |
|---|---|---|---|---|
| A. original PBMC10k → PBMC5k | 6,855 | 0.6959 ± 0.0181 | 0.6865 ± 0.0104 | −0.0094 ± 0.0212 |
| B. size-matched PBMC10k → PBMC5k | 3,994 | 0.6321 ± 0.0142 | 0.6425 ± 0.0223 | **+0.0104 ± 0.0258** |
| C. original PBMC5k → PBMC10k | 3,994 | 0.6111 ± 0.0321 | 0.6637 ± 0.0361 | +0.0526 ± 0.0549 |

The class harmonization did not distort the reference comparison: on the
per-direction eligible sets the PHASE 10 deltas were −0.0102 (A) and +0.0427 (B);
on the common 19-class set they are −0.0094 and +0.0526.

## Size effect

- `scVI_size_effect` = 0.6321 − 0.6959 = **−0.0638**
- `totalVI_size_effect` = 0.6425 − 0.6865 = **−0.0440**
- `Delta_size_effect` = +0.0104 − (−0.0094) = **+0.0197**

scVI_matched loses more from source downsampling than totalVI, by 0.0197
macro-F1. Delta moves toward Direction B but covers only about 37% of the
distance: +0.0104 versus B's +0.0526, still +0.0422 short.

## l1 (secondary) — moves the wrong way

| condition | scVI_matched | totalVI | Delta |
|---|---|---|---|
| A original | 0.8205 ± 0.0016 | 0.8248 ± 0.0169 | +0.0043 ± 0.0161 |
| A size-matched | 0.8071 ± 0.0069 | 0.7939 ± 0.0410 | **−0.0131 ± 0.0399** |
| B original | 0.7619 ± 0.0418 | 0.8592 ± 0.0429 | +0.0974 ± 0.0322 |

At l1 the size effect is −0.0174, i.e. size matching moves Direction A *away*
from Direction B. The large l1 asymmetry in PHASE 10 (+0.0974 in B) is not
reproduced at all by shrinking the PBMC10k source, and totalVI is the model that
degrades more at l1 (−0.0309 vs −0.0135).

## Source composition

Total variation distance from the full PBMC10k source is small and uniform:
l2 0.0156–0.0237, l1 0.0081–0.0153. Seed 0 is the most distorted on both levels
but only marginally. The largest single l1 proportion shift is CD4 T in seed 0
(0.3234 → 0.3345, 0.0111); the largest l2 shift is CD14 Mono in seed 0
(0.1958 → 0.1880, 0.0077). No subset is compositionally unusual, and no subset
was rejected.

## Rare-class stability

Highest F1 SD across source subsets: Eryth (totalVI 0.212), Treg (totalVI 0.155),
CD8 TCM (totalVI 0.111), CD8 Naive (scVI 0.098), gdT (scVI 0.090). The most
variable classes are not simply the rarest, but the low-support classes dominate
the tail: NK_CD56bright falls to a mean of 22 source cells (minimum 20, right at
the eligibility threshold), gdT to 33.6, ILC to 39.6, B intermediate to 38.8.

Rare T/ILC populations remain the main negative contributors, as in PHASE 10:
gdT −0.111, ILC −0.106, NK_CD56bright −0.071. Strongest positives: CD8 TEM
+0.123, CD8 Naive +0.104, CD4 Naive +0.099, Treg +0.086, B intermediate +0.074.

Notably, downsampling *improved* totalVI's relative position on several classes
that hurt it at full source size — ILC +0.064, NK_CD56bright +0.064, Treg +0.076,
CD8 TEM +0.065, CD8 Naive +0.052 — which is what produces the +0.0197 shift.

## Neighbor label agreement

All 6 aggregate cells (2 label levels × k ∈ {15, 30, 50}) still favour totalVI
after size matching, but the margin shrinks:

| condition | l2 k=15 | l2 k=50 | l1 k=15 | l1 k=50 | mean Delta |
|---|---|---|---|---|---|
| A original | +0.0167 | +0.0138 | +0.0099 | +0.0082 | +0.0121 |
| A size-matched | +0.0080 | +0.0024 | +0.0056 | +0.0045 | **+0.0051** |
| B original | +0.0367 | +0.0309 | +0.0440 | +0.0361 | +0.0368 |

At the level of individual runs, 20 of 30 size-matched comparisons favour
totalVI, versus 5/5 in most reference cells. Neighbor agreement remains the more
stable signal in sign, but it is clearly weakened by source downsampling and
does not approach the Direction B margin.

## Secondary metrics (common set, means)

Size-matched A, l2: logreg balanced accuracy 0.6474 (scVI) vs 0.6588 (totalVI),
weighted-F1 0.7695 vs 0.7843; kNN-15 macro-F1 0.5756 vs 0.5778; kNN-30 0.5495 vs
0.5457. Target-only quality: ASW 0.1294 vs 0.1101 (scVI higher), kNN purity 0.7214
vs 0.7304 (totalVI higher). Target-only metrics were computed only for the
size-matched runs; the reference conditions reuse PHASE 10 values.

## Compute

| | source s | query s | peak VRAM alloc | peak VRAM reserved | peak RSS | min RAM avail |
|---|---|---|---|---|---|---|
| scVI_matched | 64.9 | 53.5 | 406 MB | 420 MB | 3.94 GB | 3.52 GB |
| totalVI | 99.9 | 97.2 | 428 MB | 450 MB | 4.19 GB | 3.26 GB |

Total 823.7 s source + 753.6 s query = 26.3 min of GPU work. All 10 source runs
and all 10 query runs reached the full 200 epochs; early stopping never fired.
No OOM, no batch-size fallback, 0 failures.

## Interpretation

Source size is a **partial contributor, not a sufficient explanation**.

Reducing PBMC10k to the PBMC5k source size moves the l2 Delta from −0.0094 to
+0.0104, about 37% of the way to Direction B's +0.0526, and scVI_matched is the
model that suffers more from the reduction. That is directionally consistent with
protein information mattering more when RNA training data is scarce.

But three observations argue against source size as the primary driver:

1. Four of five size-matched bootstrap CIs cross zero; only seed 3 is
   unambiguously positive. Direction B is far more consistent.
2. At l1 the effect reverses. Size matching gives −0.0131 versus B's +0.0974,
   and totalVI degrades more than scVI at l1. The largest asymmetry in PHASE 10
   is entirely unexplained by source size.
3. Between-subset SD (0.0258) is larger than the size effect itself (0.0197), so
   which source cells are available matters about as much as how many.

No causal claim is made. Remaining plausible contributors, none separable with
two datasets: RNA sequencing depth and library characteristics, protein panel
informativeness and staining efficiency in each run, cell-type composition
differences beyond what random subsampling reproduces, technical/domain shift
between the two data-generation events, and donor differences.
