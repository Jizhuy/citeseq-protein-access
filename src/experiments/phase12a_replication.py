"""PHASE 12A: replicated crossed design for variance identification.

PHASE 11B ran a 5 (source subset) x 5 (model seed) x 2 (model) crossed design
with ONE observation per cell, so the subset x seed interaction is perfectly
confounded with pure run-level residual noise. The reported ~80%
"interaction/residual" component therefore admits two very different readings:
a real interaction between which source cells were drawn and which
initialisation was used, or simply irreducible run-to-run noise.

This module adds a SECOND independent replicate of every existing cell,
producing 5 x 5 x 2 replicates x 2 models = 100 observations, of which 50
already exist in PHASE 11B and are reused unchanged.

Seed factorisation (see logs/replicate_design.md)
-------------------------------------------------
``model_seed``          fixes the source-model weight initialisation. This is a
                        genuine factor level held constant within a grid cell.
``run_replicate_seed``  redraws every downstream stochastic choice: the
                        train/validation split, minibatch shuffling, dropout
                        masks, and the scArches adaptation (its new
                        batch-specific weight initialisation and its optimiser
                        path).

Replicate 0 is the existing PHASE 11B run, in which no separate run seed was
defined and the single global seed governed both initialisation and path.
Replicate 1 holds the initialisation fixed and redraws the path. Seeds are
arbitrary labels, so both are valid draws from the run-level distribution
conditional on the initialisation; the asymmetry is in bookkeeping only and is
recorded explicitly rather than papered over.
"""

from __future__ import annotations

from pathlib import Path

from src.experiments.phase11b_robustness import (
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    MODELS,
)
from src.utils.io import resolve_path

# Direction A only, matching the PHASE 11B crossed design.
CROSSED_DIRECTION_KEY = "direction_A"
CROSSED_DIRECTION = "PBMC10k_to_PBMC5k"
CROSSED_SOURCE_BATCH = "PBMC10k"
CROSSED_TARGET_BATCH = "PBMC5k"

# Replicate 0 == PHASE 11B (already trained). Replicate 1 == this phase.
REPLICATES = (0, 1)
NEW_REPLICATES = (1,)

# Offset chosen once, before any PHASE 12A result was inspected, and never tuned.
RUN_REPLICATE_SEED_BASE = 90_000


def phase12a_root() -> Path:
    return resolve_path("results/phase12a_variance_replication")


def phase12a_paths() -> dict[str, Path]:
    root = phase12a_root()
    return {
        "root": root,
        "models": root / "models",
        "embeddings": root / "embeddings",
        "probabilities": root / "probabilities",
        "posterior": root / "posterior",
        "tables": root / "tables",
        "figures": root / "figures",
        "logs": root / "logs",
        "done": root / "logs" / "done",
    }


def ensure_dirs() -> dict[str, Path]:
    paths = phase12a_paths()
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def run_replicate_seed(model_seed: int, subset_seed: int, replicate: int) -> int | None:
    """Seed governing run-level stochasticity, or None for the PHASE 11B path.

    Distinct for every (subset, model_seed, replicate) triple so that no two
    grid cells share an optimisation path.
    """
    if replicate == 0:
        return None
    return RUN_REPLICATE_SEED_BASE + 1_000 * replicate + 10 * subset_seed + model_seed


def crossed_jobs(replicates=NEW_REPLICATES) -> list[dict]:
    """The (subset x seed x model x replicate) jobs for the given replicates."""
    jobs: list[dict] = []
    for subset_seed in CROSSED_SUBSET_SEEDS:
        for model_seed in CROSSED_MODEL_SEEDS:
            for model in MODELS:
                for replicate in replicates:
                    jobs.append({
                        "direction_key": CROSSED_DIRECTION_KEY,
                        "direction": CROSSED_DIRECTION,
                        "model": model,
                        "model_seed": model_seed,
                        "subset_seed": subset_seed,
                        "replicate": replicate,
                        "run_replicate_seed": run_replicate_seed(
                            model_seed, subset_seed, replicate),
                    })
    return jobs
