# PHASE 11B methods log

Generated 2026-09-09T12:50:13.154482+00:00

## Compute environment

- GPU: NVIDIA GeForce RTX 4060 Laptop GPU
- torch 2.14.0+cu126, CUDA 12.6
- scvi-tools 1.3.3
- Python 3.11.7 on Linux-6.18.33.2-microsoft-standard-WSL2-x86_64-with-glibc2.43
- total model training + adaptation time: 7.27 h across 130 runs

## Runs completed

- full-source 20-seed grid: 80 runs (2 directions x 2 models x 20 seeds = 80 expected)
- crossed grid: 50 runs (5 source subsets x 5 model seeds x 2 models = 50 expected)
- all artifact-integrity checks passed: True

## Model configuration

Unchanged from the validated PHASE 10 configuration. No architecture was tuned.

- scVI_matched: n_latent 20, n_hidden 256, n_layers 2, dropout 0.2, gene_likelihood nb, dispersion gene, latent_distribution normal
- totalVI: n_latent 20, n_hidden 256, encoder layers 2, decoder layers 1, dropout 0.2, gene_likelihood nb, gene_dispersion gene, protein_dispersion protein
- training: batch_size 256, max_epochs 200, train_size 0.9, early stopping on elbo_validation, patience 45, mode min, check_val_every_n_epoch 1
- scArches query adaptation: prepare_query_anndata + load_query_data, max_epochs 200, plan weight_decay 0.0; totalVI additionally keeps its model-specific lr and reduce_lr_on_plateau setting
- target labels were never used in training, adaptation, early stopping, or classifier fitting

## Latent coordinate rule

Both SOURCE and TARGET are encoded with the POST-adaptation query model. 
`assert_post_adaptation` is a regression guard against the PHASE 10 defect in 
which the pre-adaptation source array was handed to the evaluator; it checks 
object identity, shape, and content before any metric is computed.

### Source latent drift, pre- to post-adaptation

```
             drift_mean           drift_median           drift_frob          
                   mean       max         mean       max       mean       max
model                                                                        
scvi_matched   0.000000  0.000000     0.000000  0.000000   0.000000  0.000000
totalvi        0.017138  0.021148     0.017037  0.020859   0.003528  0.004284
```

- scVI runs identical pre/post: True
- totalVI runs identical pre/post: False

The non-zero totalVI drift is exactly the quantity that made the PHASE 10 defect 
numerically visible for totalVI while leaving scVI unaffected.

## Artifact persistence and integrity

Every run writes `models/<direction>/<model>/seedNN/source_model/` and 
`.../adapted_query_model/`. Immediately after saving, the in-memory model is 
deleted, reloaded from disk, and its latent representation and posterior variance 
are recomputed on a fixed validation subset.

- tolerance: 1e-05
- max observed latent difference across all runs: 0.000e+00
- max observed variance difference across all runs: 0.000e+00

## Posterior uncertainty

- API: `get_latent_representation(..., return_dist=True)`
- scvi-tools 1.3.3 `VAEMixin.get_latent_representation` returns `(qz.loc, qz.scale.square())`; neither SCVI nor TOTALVI overrides it
- the second element is therefore the posterior VARIANCE, not a scale or SD
- `latent_distribution=normal`, so q(z|x) is a diagonal Gaussian and the variance is exact
- primary scalar: U_latent_i = mean_j var_ij over the 20 latent dimensions
- secondary diagnostics saved: median and max posterior variance per run

### Monte Carlo confirmation

- seeds checked: [0, 5, 10, 15]; S = 50 posterior draws per target cell
- mean empirical/closed-form variance ratio: 0.9797226153314114
- consistent with a variance: True

## Downstream classifier

Identical to PHASE 10: `LogisticRegression(max_iter=2000, solver='lbfgs', 
random_state=seed, class_weight=None)` fitted on SOURCE post-adaptation latents 
and SOURCE labels, evaluated on TARGET. Probability columns are aligned by class 
NAME before any aggregation.

## Uncertainty definitions

- U_pred = 1 - max class probability (primary)
- predictive entropy and top1-top2 margin (secondary)
- U_latent = mean posterior variance (primary latent)
- ensemble disagreement = H(mean_s p_s) - mean_s H(p_s), described as epistemic-like, not exact Bayesian mutual information
- variation ratio = 1 - modal predicted-class frequency

No composite reliability score was constructed. No temperature scaling was fitted: 
PHASE 11 found it did not reliably improve calibration, and PHASE 11B reports raw 
classifier probabilities only.

## Calibration and selective prediction

- ECE: 15 bins, equal-frequency (primary) and equal-width (sensitivity)
- NLL and multiclass Brier on the full probability vector
- coverage grid: ['100%', '90%', '80%', '70%', '60%', '50%']
- AURC computed over the full risk-coverage curve, cells ranked by U_pred

## Bootstrap

- seed-level: 10000 paired replicates over the 20 seed pairs
- target-cell: 500 paired replicates per seed, identical indices for scVI and totalVI
- base seed: 11000
- the two bootstraps are reported separately and never pooled

## Evidence classification rule

Fixed before any PHASE 11B result was inspected:

- **robust**: seed-level paired bootstrap 95% CI excludes 0 AND at least 16/20 seeds share the sign of the mean delta
- **seed_sensitive**: the majority of per-seed target-cell bootstrap CIs exclude 0, but fewer than 16/20 seeds agree in sign OR the seed-level CI includes 0: the apparent effect depends on which model seed was drawn
- **cell_sampling_sensitive**: seed-level CI excludes 0 and signs are consistent, but fewer than half of the per-seed target-cell bootstrap CIs exclude 0: the effect is stable across seeds yet small relative to target-cell sampling noise
- **inconclusive**: neither the seed-level nor the target-cell evidence separates the models

## Variance decomposition

Balanced two-factor method-of-moments decomposition on the 5x5 (source subset x model seed) tables. With one observation per grid cell the 
interaction cannot be separated from pure error, so the third component is 
reported as combined interaction/residual. Negative raw estimates are truncated 
at zero for the reported fractions and the raw values are preserved in the same 
table. With only five levels per factor these components are exploratory.

## Resource usage

```
             runtime        gpu_alloc        gpu_res            rss        
                mean    max      mean    max    mean    max    mean     max
model                                                                      
scvi_matched   149.4  206.1     309.5  309.5   318.6  320.0  3597.0  6477.4
totalvi        253.3  376.3     379.5  379.5   398.8  400.0  3473.6  6864.3
```

- peak GPU allocated across all runs: 379 MB
- peak GPU reserved across all runs: 400 MB
- peak process RSS across all runs: 6864 MB
- neighbor k (where used): 15

## Rerun / resume

```bash
ssh 4060-server
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate multiomics_robustness
cd /home/jizhu/research/multiomics_robustness
export TMPDIR="$HOME/tmp" MPLBACKEND=Agg PATH="/usr/lib/wsl/lib:$PATH"
tmux new -s phase11b

python scripts/11b_run_grid.py --stage gate1
python scripts/11b_run_grid.py --stage gate2
python scripts/11b_run_grid.py --stage seeds --seeds 1-4
python scripts/11b_run_grid.py --stage seeds --seeds 5-19
python scripts/11b_run_grid.py --stage crossed
python scripts/11b_mc_check.py --extra-s100-seed0
python scripts/11b_analyze.py --stage all
python scripts/11b_crossed.py --stage all
python scripts/11b_figures.py
python scripts/11b_methods_log.py
```

Every training run writes a done marker under `logs/done/`. Re-running any stage 
skips completed runs, so an interrupted grid resumes without retraining.
