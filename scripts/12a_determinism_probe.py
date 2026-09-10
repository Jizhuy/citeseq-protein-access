#!/usr/bin/env python
"""PHASE 12A determinism probe (§2): can a legitimate replicate even be defined?

The PHASE 11B crossed design has one observation per (source subset x model
seed x model) cell, so interaction and pure run-level residual are confounded.
Adding a second replicate is only meaningful if two runs at the SAME model seed
are genuinely different. `src.utils.seed.set_global_seed` seeds python, numpy,
torch, CUDA and `scvi.settings.seed`, and sets `cudnn.deterministic=True`, so
the honest prior is that same-seed runs are bit-identical and a naive duplicate
would be FAKE replication.

This probe measures that directly, and validates a candidate replicate
mechanism, on the real subset0 Direction A data:

  A  set_global_seed(S); build; train
  B  set_global_seed(S); build; train           -> A vs B tests determinism
  C  set_global_seed(S); build; set_global_seed(R); train
                                                -> A vs C tests the mechanism

C isolates weight initialisation (model_seed S) from the stochastic
optimisation path (run_replicate_seed R): shuffling, dropout masks and the
train/validation split. If A == B bit-exactly and A != C, then C is the only
defensible replicate definition and it must be documented as an explicitly
added second seed rather than a rerun.

Usage
  python scripts/12a_determinism_probe.py --models scvi_matched totalvi
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import (  # noqa: E402
    ensure_csr,
    load_source_target,
    train_scvi_source,
    train_totalvi_source,
)
from src.experiments.phase10b_source_size import SUBSAMPLE_N, subsample_index  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    DEFAULT_HYPERPARAMETERS,
    MATCHED_HYPERPARAMETERS,
    cleanup,
    utc_now,
)
from src.experiments.phase12a_replication import ensure_dirs  # noqa: E402
from src.utils.io import load_config, save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402
from src.utils.seed import set_global_seed  # noqa: E402

PROBE_MODEL_SEED = 0
PROBE_RUN_SEED = 9000
PROBE_SUBSET = 0


def _hparams(model_name: str, cfg: dict) -> dict:
    if model_name == "scvi_matched":
        h = dict(MATCHED_HYPERPARAMETERS)
        h["early_stopping"] = bool(cfg["models"]["early_stopping"])
    else:
        h = dict(DEFAULT_HYPERPARAMETERS)
    h["max_epochs"] = int(cfg["models"]["max_epochs"])
    h["batch_size"] = int(cfg["models"]["batch_size"])
    return h


def _train_once(model_name: str, source, hparams: dict, model_seed: int,
                run_seed: int | None, logger) -> dict:
    """One source-training run. run_seed=None reproduces PHASE 11B exactly."""
    set_global_seed(model_seed)
    trainer = train_scvi_source if model_name == "scvi_matched" else train_totalvi_source
    with tempfile.TemporaryDirectory() as td:
        out = trainer(source, hparams, Path(td), logger, run_seed=run_seed)
    latent = np.asarray(out["latent"], dtype=np.float64)
    cleanup(purge_manager_store=True)
    return {"latent": latent, "epochs": out["epochs"],
            "elbo_validation_final": out["elbo_validation_final"]}


def _compare(a: np.ndarray, b: np.ndarray) -> dict:
    d = np.abs(a - b)
    denom = float(np.linalg.norm(a)) or 1.0
    return {
        "bit_identical": bool(np.array_equal(a, b)),
        "max_abs_diff": float(d.max()),
        "mean_abs_diff": float(d.mean()),
        "relative_frobenius": float(np.linalg.norm(a - b) / denom),
        "mean_cell_displacement": float(np.linalg.norm(a - b, axis=1).mean()),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["scvi_matched", "totalvi"])
    args = ap.parse_args()

    paths = ensure_dirs()
    logger = get_logger("phase12a_probe", paths["logs"] / "determinism_probe.log")
    cfg = load_config()

    logger.info("loading Direction A source/target")
    source_full, target, _ = load_source_target("PBMC10k", "PBMC5k")
    idx = subsample_index(source_full.n_obs, SUBSAMPLE_N, PROBE_SUBSET)
    source = ensure_csr(source_full[idx].copy())
    del source_full, target
    cleanup()
    logger.info("probe source subset%s n=%s", PROBE_SUBSET, source.n_obs)

    results = {}
    for model_name in args.models:
        h = _hparams(model_name, cfg)
        logger.info("=== %s: run A (model_seed=%s) ===", model_name, PROBE_MODEL_SEED)
        a = _train_once(model_name, source, h, PROBE_MODEL_SEED, None, logger)
        logger.info("=== %s: run B (identical call) ===", model_name)
        b = _train_once(model_name, source, h, PROBE_MODEL_SEED, None, logger)
        logger.info("=== %s: run C (model_seed=%s, run_replicate_seed=%s) ===",
                    model_name, PROBE_MODEL_SEED, PROBE_RUN_SEED)
        c = _train_once(model_name, source, h, PROBE_MODEL_SEED, PROBE_RUN_SEED, logger)

        ab = _compare(a["latent"], b["latent"])
        ac = _compare(a["latent"], c["latent"])
        results[model_name] = {
            "A_vs_B_same_seed_rerun": ab,
            "A_vs_C_new_run_replicate_seed": ac,
            "epochs": {"A": a["epochs"], "B": b["epochs"], "C": c["epochs"]},
            "elbo_validation_final": {
                "A": a["elbo_validation_final"],
                "B": b["elbo_validation_final"],
                "C": c["elbo_validation_final"],
            },
            "same_seed_rerun_is_fake_replication": ab["bit_identical"],
            "run_replicate_seed_produces_independent_run": not ac["bit_identical"],
        }
        logger.info("%s A vs B bit_identical=%s max_abs=%.3e | A vs C bit_identical=%s max_abs=%.3e",
                    model_name, ab["bit_identical"], ab["max_abs_diff"],
                    ac["bit_identical"], ac["max_abs_diff"])

    verdict = {
        "probe_model_seed": PROBE_MODEL_SEED,
        "probe_run_replicate_seed": PROBE_RUN_SEED,
        "probe_subset": PROBE_SUBSET,
        "source_n": int(source.n_obs),
        "all_same_seed_reruns_bit_identical": all(
            r["same_seed_rerun_is_fake_replication"] for r in results.values()),
        "all_run_seed_runs_independent": all(
            r["run_replicate_seed_produces_independent_run"] for r in results.values()),
        "per_model": results,
        "timestamp_utc": utc_now(),
    }
    save_json(verdict, paths["logs"] / "determinism_probe.json")
    print(json.dumps(verdict, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
