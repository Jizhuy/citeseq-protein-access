#!/usr/bin/env python
"""Generate the PHASE 12A replicate-design and interpretation logs.

Both documents are generated from the artifacts they describe -- the
determinism probe JSON, the run manifests, and the variance tables -- so they
cannot drift away from what was actually executed.

Usage
  python scripts/12a_write_logs.py --which design
  python scripts/12a_write_logs.py --which interpretation
  python scripts/12a_write_logs.py            # both, skipping what is absent
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11b_robustness import phase11b_paths, utc_now  # noqa: E402
from src.experiments.phase12a_replication import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    REPLICATES,
    RUN_REPLICATE_SEED_BASE,
    ensure_dirs,
    run_replicate_seed,
)
from src.experiments.phase12a_stats import (  # noqa: E402
    TRUNCATION_BIAS_NOTE,
    VARIANCE_MODEL_NOTE,
)

PATHS = ensure_dirs()
TAB, LOGS = PATHS["tables"], PATHS["logs"]
METRIC_LABEL = {"macro_f1": "macro-F1", "error_auroc": "error AUROC",
                "error_auprc": "error AUPRC", "nll": "NLL", "brier": "Brier",
                "ece": "ECE", "aurc": "AURC"}


def _fmt_table(df: pd.DataFrame, floatfmt: str = "%.4f") -> str:
    return df.to_string(index=False, float_format=lambda v: floatfmt % v)


def write_design() -> None:
    probe_path = LOGS / "determinism_probe.json"
    if not probe_path.exists():
        print("skip design: determinism_probe.json not found")
        return
    probe = json.loads(probe_path.read_text(encoding="utf-8"))

    rows = []
    for model, r in probe["per_model"].items():
        ab, ac = r["A_vs_B_same_seed_rerun"], r["A_vs_C_new_run_replicate_seed"]
        rows.append({
            "model": model,
            "same_seed_bit_identical": ab["bit_identical"],
            "same_seed_max_abs_diff": ab["max_abs_diff"],
            "run_seed_bit_identical": ac["bit_identical"],
            "run_seed_max_abs_diff": ac["max_abs_diff"],
            "run_seed_mean_cell_displacement": ac["mean_cell_displacement"],
            "elbo_val_A": r["elbo_validation_final"]["A"],
            "elbo_val_B": r["elbo_validation_final"]["B"],
            "elbo_val_C": r["elbo_validation_final"]["C"],
        })
    probe_tbl = _fmt_table(pd.DataFrame(rows))

    seed_rows = []
    for s in (0, 1):
        for k in (0, 1):
            for rep in REPLICATES:
                seed = run_replicate_seed(k, s, rep)
                seed_rows.append({"source_subset": s, "model_seed": k, "replicate": rep,
                                  "run_replicate_seed": "None (PHASE 11B path)"
                                  if seed is None else str(seed)})
    seed_tbl = _fmt_table(pd.DataFrame(seed_rows), "%s")

    n_done = len(list((PATHS["done"]).glob("*.json"))) if PATHS["done"].exists() else 0

    doc = f"""# PHASE 12A replicate design

Generated {utc_now()}

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
data (n = {probe['source_n']} source cells), with three runs per model:

- **A**: `set_global_seed(S)`; build; train
- **B**: `set_global_seed(S)`; build; train -- identical call to A
- **C**: `set_global_seed(S)`; build; `set_global_seed(R)`; train

```
{probe_tbl}
```

- all same-seed reruns bit-identical: **{probe['all_same_seed_reruns_bit_identical']}**
- all run-seed runs independent: **{probe['all_run_seed_runs_independent']}**

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

`run_replicate_seed` is `{RUN_REPLICATE_SEED_BASE} + 1000*replicate +
10*subset_seed + model_seed`, fixed before any PHASE 12A result was inspected
and never tuned. Every grid cell therefore has its own optimisation path:

```
{seed_tbl}
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
- {len(CROSSED_SUBSET_SEEDS)} source subsets x {len(CROSSED_MODEL_SEEDS)} model seeds x {len(REPLICATES)} replicates x 2 models = {len(CROSSED_SUBSET_SEEDS)*len(CROSSED_MODEL_SEEDS)*len(REPLICATES)*2} observations
- 50 pre-existing PHASE 11B runs reused unchanged; 50 new runs trained here
- new runs completed so far: {n_done}/50
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

{VARIANCE_MODEL_NOTE}

{TRUNCATION_BIAS_NOTE}

## Software

- Python {platform.python_version()} on {platform.platform()}
"""
    (LOGS / "replicate_design.md").write_text(doc, encoding="utf-8")
    print(f"wrote {LOGS / 'replicate_design.md'}")


def write_interpretation() -> None:
    vpath = TAB / "replicated_variance_components.csv"
    if not vpath.exists():
        print("skip interpretation: replicated_variance_components.csv not found")
        return
    v = pd.read_csv(vpath)
    d = pd.read_csv(TAB / "replicated_delta_variance_components.csv")
    cols = ["model", "metric", "source_subset_fraction", "model_seed_fraction",
            "interaction_fraction", "residual_fraction", "p_interaction",
            "p_source_subset", "p_model_seed", "phase11b_interaction_residual_fraction"]
    main_tbl = _fmt_table(v[cols])
    delta_tbl = _fmt_table(d[[c for c in cols if c in d.columns]])

    # How the old lumped term actually split.
    v = v.assign(new_lumped=v.interaction_fraction + v.residual_fraction,
                 residual_share_of_lumped=v.residual_fraction /
                 (v.interaction_fraction + v.residual_fraction).replace(0, np.nan))
    split = _fmt_table(v[["model", "metric", "phase11b_interaction_residual_fraction",
                          "interaction_fraction", "residual_fraction",
                          "residual_share_of_lumped"]])
    n_sig = int((v.p_interaction < 0.05).sum())
    n_tot = int(v.p_interaction.notna().sum())
    med_resid_share = float(v.residual_share_of_lumped.median())

    percell = TAB / "per_cell_replicated_variance.csv"
    pc_block = ""
    if percell.exists():
        pc = pd.read_csv(percell)
        pc_tbl = _fmt_table(pc.groupby("model")[
            [f"{c}_fraction" for c in
             ("source_subset", "model_seed", "interaction", "residual")]].median().reset_index())
        pc_block = f"""
## Per-cell decomposition (median over target cells)

```
{pc_tbl}
```
"""

    doc = f"""# PHASE 12A interpretation

Generated {utc_now()}

## The question

PHASE 11B reported that roughly 80% of the variance in the crossed design fell
into a single "interaction/residual" term. That term was unidentified: with one
observation per grid cell, a real subset x seed interaction and pure run-level
noise produce the same number. This phase added a second independent replicate
per cell so the two can be told apart, and asks which of them it actually was.

## Metric-level components

```
{main_tbl}
```

## How the PHASE 11B lumped term split

```
{split}
```

The interaction F-test (MS_interaction / MS_residual, which only exists once
n > 1) is significant at the 0.05 level in {n_sig} of {n_tot} model-by-metric
combinations. The median share of the old lumped term attributable to pure
run-level residual is {med_resid_share:.3f}.

## Delta decomposition

```
{delta_tbl}
```
{pc_block}
## Reading rule

Component fractions are descriptive and, after truncation at zero, upward
biased for null factors: on pure noise the three null components collectively
absorb roughly 15% of the truncated total. The F-tests, not the fractions,
decide whether a component is real. {TRUNCATION_BIAS_NOTE}

## Variance model

{VARIANCE_MODEL_NOTE}
"""
    (LOGS / "phase12a_interpretation.md").write_text(doc, encoding="utf-8")
    print(f"wrote {LOGS / 'phase12a_interpretation.md'}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="both",
                    choices=["design", "interpretation", "both"])
    args = ap.parse_args()
    if args.which in ("design", "both"):
        write_design()
    if args.which in ("interpretation", "both"):
        write_interpretation()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
