# PHASE 11 methods

Formal uncertainty characterization and calibration for cross-dataset transfer.
No representation model was trained, fine-tuned, or modified.

## Compute environment

| item | value |
|---|---|
| host | RTX 4060 Laptop GPU, Ubuntu 26.04.1 LTS under WSL2 |
| VRAM / RAM | ~8188 MiB / ~7.6 GiB |
| Python | 3.11.7 |
| PyTorch | 2.14.0+cu126, CUDA 12.6 |
| scvi-tools | 1.3.3 |
| project root | `/home/jizhu/research/multiomics_robustness` |
| env | `conda activate multiomics_robustness` |
| env vars | `TMPDIR=$HOME/tmp`, `MPLBACKEND=Agg`, `PATH=/usr/lib/wsl/lib:$PATH` |

## Model artifacts consumed

PHASE 10 cross-dataset embeddings, read-only:

```
results/phase10_cross_dataset/embeddings/<direction>/
  <model>_<direction>_seed<k>_cuda_source_latent.npy   # POST-adaptation
  <model>_<direction>_seed<k>_cuda_target_latent.npy   # POST-adaptation
  <model>_<direction>_seed<k>_cuda_{source,target}_cells.csv
```

- Directions: `pbmc10k_to_pbmc5k` (A), `pbmc5k_to_pbmc10k` (B)
- Models: `scvi_matched`, `totalvi`
- Model seeds: 0, 1, 2, 3, 4  (20 embedding pairs)
- Latent dimension 20; `latent_distribution="normal"` for both models

PHASE 10B artifacts for source-subset instability:
`results/phase10b_source_size/embeddings/<model>_pbmc10k_sub3994_subsample_seed<s>_modelseed0_cuda_*`,
subsample seeds 0-4, n = 3994, Direction A only.

Cell-order identity between each saved `*_cells.csv` and the AnnData row order
is asserted before use; a mismatch raises.

### Latent coordinate rule (§4)

Source and target are both taken from the post-adaptation arrays. The
pre-adaptation array `*_source_latent_source_model.npy` is never used in a
cross-space metric. A defect in PHASE 10 that violated this for freshly-trained
totalVI runs was found and fixed during this phase; see
`../../phase10_cross_dataset/logs/phase10_erratum_latent_alignment.md`. PHASE 11
uses the corrected reference tables in
`results/phase10_cross_dataset/tables_corrected/`.

## Posterior API

`VAEMixin.get_latent_representation(..., return_dist=True)` returns
`(qz.loc, qz.scale.square())`, i.e. closed-form posterior mean and variance.
Verified on a real saved model. Monte Carlo sampling is **not** required and was
**not used**; `MC_SAMPLES` is therefore not applicable and no S=50 vs S=100
sensitivity check was run.

Target-cell latent posterior is **unsupported**: PHASE 10 persisted only the
pre-adaptation source reference model, and an untrained query model refuses
inference. Full argument in `posterior_api_validation.md`. Source-cell posterior
variance was extracted for all 20 models and is exactly valid in the PHASE 10
coordinate system for scVI only (scArches drift 0.0); totalVI reference models
define the pre-adaptation space and are flagged
`valid_in_phase10_coordinate_system = False`.

## Label set

Primary `cell_type_l2`; secondary `cell_type_l1`. Eligible classes come from
PHASE 10's own `class_plan`, unmodified: high-confidence cells
(`annotation_tier_l2 == "high"`, PHASE 6 rule, l2 confidence >= 0.85),
intersected across source and target, with `MIN_SOURCE_SUPPORT = 20` and
`MIN_TARGET_SUPPORT = 10`. No threshold was changed.

Resulting l2 evaluation sets: Direction A 22 classes / 6040 source / 3369
target; Direction B 21 classes / 3355 source / 5964 target.

## Classifier

Identical to PHASE 10:

```python
LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed, class_weight=None)
```

Fitted on source post-adaptation latents with source labels; evaluated on target
post-adaptation latents. Target labels never enter any fit. Same settings for
both models, no per-model tuning.

## Reproduction gate (§11)

Every direction x model x seed x label-level cell (40) was required to reproduce
the corrected PHASE 10 logreg macro-F1 within `1e-6`. All 40 matched to
`0.000000`. Recorded in `phase10_reproduction_check.json`.

## Uncertainty definitions

| concept | definition |
|---|---|
| B primary | `U_pred = 1 - max_c p_c` |
| B secondary | predictive entropy `-sum p log p`; margin `p_top1 - p_top2` |
| C ensemble mean | `p_bar = mean_s p^(s)` over 5 model seeds |
| C disagreement | `H(p_bar) - mean_s H(p^(s))` — "ensemble disagreement", **not** exact Bayesian mutual information |
| C variation ratio | `1 - max_class_count / 5` |
| D subset instability | across 5 PHASE 10B source subsamples: variation ratio, entropy of mean probability, SD of true-class probability, mean per-class SD |
| A latent | `U_latent = mean_j sigma_ij^2` — source cells only |

Probability vectors are aligned by class **name** before any aggregation;
mismatched class sets raise rather than silently reindex.

## Calibration

- Uncalibrated results are always primary and always reported.
- Post-hoc temperature scaling fitted on **source labels only**: a stratified
  25% source holdout (`StratifiedShuffleSplit`, `random_state=1111`), a companion
  classifier fitted on the complementary 75%, and T minimising held-out NLL over
  160 log-spaced values in [0.25, 8]. That T is then applied to the full-source
  classifier's target logits. Reported separately in `*_temperature_scaled`
  columns. Caveat: T is estimated from a 75%-data companion classifier, so it is
  an out-of-sample but not perfectly matched estimate.
- ECE: 15 **equal-frequency** bins (primary) and 15 **equal-width** bins
  (sensitivity). Both are always reported; neither was chosen after inspection.
- NLL and multiclass Brier `mean_i sum_c (p_ic - y_ic)^2`.

## Selective prediction

Coverage grid 100/90/80/70/60/50%, cells ranked by ascending `U_pred`. Risk is
`1 - accuracy`. AURC is the mean risk over all n prefixes of the confidence
ordering. Macro-F1 vs coverage is reported alongside accuracy because of l2
class imbalance.

Error enrichment strata were pre-specified: top 10%, top 20%, bottom 20%.

## Neighborhood stability

k = 15 Euclidean nearest **source** neighbors of each target cell, computed
separately per model seed. `neighbor_stability` is the mean pairwise Jaccard
overlap of the 5 neighbor sets (10 pairs). Neighbor label agreement is the
fraction of the k neighbors sharing the target cell's label, averaged and
SD'd across seeds.

## Bootstrap

500 target-cell resamples, `numpy.random.default_rng(1111)`. The **same** index
matrix is reused for scVI and totalVI so the totalVI - scVI delta is paired.
95% CIs are the 2.5/97.5 percentiles of the delta draws.

Three variability sources are reported separately and never pooled:
between-model-seed spread, target-cell bootstrap, PHASE 10B source-subset spread.

## Reproduction command

```bash
ssh 4060-server
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate multiomics_robustness
export TMPDIR="$HOME/tmp" MPLBACKEND=Agg PATH="/usr/lib/wsl/lib:$PATH"
cd /home/jizhu/research/multiomics_robustness

python scripts/11_gate.py --direction pbmc5k_to_pbmc10k --seed 0 --level l2 --probe-api
python scripts/11_gate.py --direction pbmc10k_to_pbmc5k --seed 0 --level l2
python scripts/10c_rescore_corrected.py          # PHASE 10 erratum, writes tables_corrected/
python scripts/11_run_uncertainty.py --stage main
python scripts/11_subset_latent.py --stage subset
python scripts/11_subset_latent.py --stage latent
python scripts/11_figures.py
```

Timestamp: PHASE 11 executed 2026-09-08 (UTC+10).
