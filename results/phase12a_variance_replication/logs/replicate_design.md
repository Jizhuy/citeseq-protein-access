# PHASE 12A replicate design

Generated 2026-09-09T13:44:54.082491+00:00

## Why a replicate was needed

The PHASE 11B crossed design had exactly one observation per
(source subset x model seed x model) cell. In a two-factor layout with n = 1
the interaction mean square has no independent error term to be compared
against, so the subset x seed interaction and the pure run-level residual are
the same number. PHASE 11B therefore reported them as a single combined
"interaction/residual" component, which accounted for roughly 80% of the
variance. That number was not interpretable as irreducible noise, and this
phase exists to split it.

## Is a replicate even definable? (probe result)

`src.utils.seed.set_global_seed` seeds Python, NumPy, Torch and CUDA, sets
`scvi.settings.seed`, and sets `cudnn.deterministic = True` with
`cudnn.benchmark = False`. The honest prior was therefore that a same-seed
rerun is bit-identical, which would make a naive "second replicate" a
duplicate rather than a replicate.

This was measured directly rather than assumed, on the real subset0 Direction A
data (n = 3994 source cells), with three runs per model:

- **A**: `set_global_seed(S)`; build; train
- **B**: `set_global_seed(S)`; build; train -- identical call to A
- **C**: `set_global_seed(S)`; build; `set_global_seed(R)`; train

```
       model  same_seed_bit_identical  same_seed_max_abs_diff  run_seed_bit_identical  run_seed_max_abs_diff  run_seed_mean_cell_displacement  elbo_val_A  elbo_val_B  elbo_val_C
scvi_matched                     True                  0.0000                   False                 9.0743                           4.5028   4451.2397   4451.2397   4466.6206
     totalvi                     True                  0.0000                   False                 6.9119                           4.3842   4524.5879   4524.5879   4507.7910
```

- all same-seed reruns bit-identical: **True**
- all run-seed runs independent: **True**

Runs A and B agree to the last bit (max absolute latent difference 0.000e+00)
and produce identical final validation ELBO for both models. **Duplicating a
PHASE 11B run at the same seed would have produced exactly zero residual
variance and a fabricated interaction estimate.** Run C differs substantially,
so an explicitly added run-level seed is the only defensible replicate
mechanism. This is the STOP condition in the phase specification, and it fired.

## Seed factorisation

Two seeds are now documented separately.

| seed | governs | held fixed within a grid cell? |
|---|---|---|
| `model_seed` | source-model weight initialisation | yes -- it is the factor level |
| `run_replicate_seed` | train/validation split, minibatch shuffling, dropout masks, and the scArches adaptation (new batch-specific weight initialisation and its optimiser path) | no -- redrawn per replicate |

The split is implemented by `phase10_cross_dataset._apply_run_seed`, called
between module construction and `model.train`. Weight initialisation happens in
the constructor, so it is fixed by `model_seed`; scvi's `DataSplitter` reads
`settings.seed` at setup time and the loader/dropout draw from the global Torch
generator, so all of those are redrawn. For the adaptation step the seed used
is `run_replicate_seed + 1`, so the two training stages of one run do not share
an identical stream.

`run_replicate_seed` is `90000 + 1000*replicate +
10*subset_seed + model_seed`, fixed before any PHASE 12A result was inspected
and never tuned. Every grid cell therefore has its own optimisation path:

```
 source_subset  model_seed  replicate    run_replicate_seed
             0           0          0 None (PHASE 11B path)
             0           0          1                 91000
             0           1          0 None (PHASE 11B path)
             0           1          1                 91001
             1           0          0 None (PHASE 11B path)
             1           0          1                 91010
             1           1          0 None (PHASE 11B path)
             1           1          1                 91011
```

## An asymmetry that is recorded rather than hidden

Replicate 0 is the existing PHASE 11B run, in which no separate run seed
existed and one global seed governed both initialisation and path. Replicate 1
holds initialisation fixed and redraws the path. So replicate 0 is the special
case `run_replicate_seed == model_seed` and replicate 1 is not.

Seeds are arbitrary labels with no intrinsic meaning, so both replicates are
valid draws from the run-level distribution conditional on the initialisation,
and the decomposition is unbiased. The asymmetry is bookkeeping only. The
alternative -- retraining all 100 runs under a symmetric scheme -- would have
cost twice the compute and discarded 50 validated runs; the phase specification
explicitly required keeping the original replicate.

What is NOT claimed: that `run_replicate_seed` isolates a single mechanism. It
bundles the split, the shuffling, the dropout masks and the adaptation. The
residual component is therefore "all run-level stochasticity downstream of
initialisation", not "minibatch order alone".

## Design executed

- direction: Direction A only (PBMC10k -> PBMC5k), matching PHASE 11B
- 5 source subsets x 5 model seeds x 2 replicates x 2 models = 100 observations
- 50 pre-existing PHASE 11B runs reused unchanged; 50 new runs trained here
- new runs completed so far: 2/50
- held fixed within a cell: source subset membership (written to
  `logs/crossed_subset*_cells.csv`), model seed, architecture, hyperparameters,
  direction, data, labels, and the common l2 class set copied verbatim from
  PHASE 11B
- no hyperparameter, architecture or class-definition change
- target labels were never used in training, adaptation, or classifier fitting

## Artifact layout

```
results/phase12a_variance_replication/models/direction_A/<model>/subset<S>_seed<NN>/replicate01/
    source_model/
    adapted_query_model/
```

This follows the PHASE 11B convention rather than the flatter example in the
phase brief, so the two replicate trees are directly comparable. Each run also
persists post-adaptation source and target embeddings, target posterior mean
and variance, and a manifest. PHASE 11B is never written to.

## Variance model

Balanced two-factor random-effects ANOVA with n replicates per cell (Searle, Casella & McCulloch, Variance Components, ch. 4). Method-of-moments estimators from expected mean squares. Interaction and pure residual are separately identified only because n > 1; with n = 1 (the PHASE 11B design) MS_E is undefined and the two collapse into one term.

Method-of-moments component estimates are unbiased on the raw scale but fluctuate around zero for a null factor. After truncation at zero a null factor therefore receives a strictly positive fraction, and with a single 5x5x2 table that upward bias is large: on pure noise the residual typically retains only about half of the truncated total. Fractions are descriptive; the F-tests, not the fractions, decide whether a component is real.

## Software

- Python 3.11.7 on Linux-6.18.33.2-microsoft-standard-WSL2-x86_64-with-glibc2.43
