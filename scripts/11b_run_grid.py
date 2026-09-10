#!/usr/bin/env python
"""PHASE 11B grid runner: prospective 20-seed replication + 5x5 crossed design.

Every run trains a source model, performs scArches query adaptation, persists
BOTH models, verifies the adapted model reloads bit-faithfully, and extracts
target posterior uncertainty. Resumable via done markers.

Stages
  gate1    Direction B, both models, seed 0
  gate2    Direction A, both models, seed 0
  seeds    both directions, both models, --seeds range
  crossed  Direction A, 5 PHASE 10B source subsets x 5 model seeds x 2 models

Usage
  python scripts/11b_run_grid.py --stage gate1
  python scripts/11b_run_grid.py --stage seeds --seeds 1-4
  python scripts/11b_run_grid.py --stage seeds --seeds 5-19
  python scripts/11b_run_grid.py --stage crossed
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import (  # noqa: E402
    class_plan,
    ensure_csr,
    load_source_target,
)
from src.experiments.phase10b_source_size import SUBSAMPLE_N, subsample_index  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    DIRECTIONS,
    MODELS,
    cleanup,
    ensure_dirs,
    experiment_tag,
    ram_snapshot,
    run_experiment,
    utc_now,
)
from src.utils.io import load_config, save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

# Abort rather than let the OOM killer land mid-write on a 7.6 GB machine.
MIN_FREE_RAM_MB = 500.0


def parse_seeds(spec: str) -> list[int]:
    out: list[int] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if "-" in chunk:
            lo, hi = chunk.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(chunk))
    return sorted(set(out))


def write_class_plan(paths, direction_key, direction, source, target, logger) -> dict:
    """Eligible-class definition, fixed from labels before any model runs."""
    plans = {}
    for level, col in (("l2", "cell_type_l2"), ("l1", "cell_type_l1")):
        plan = class_plan(source, target, col)
        plans[level] = {
            k: (v.tolist() if isinstance(v, np.ndarray) else v)
            for k, v in plan.items()
            if k not in ("source_eval_index", "target_eval_index")
        }
        np.save(paths["logs"] / f"{direction_key}_{level}_source_eval_index.npy",
                plan["source_eval_index"])
        np.save(paths["logs"] / f"{direction_key}_{level}_target_eval_index.npy",
                plan["target_eval_index"])
        logger.info("%s %s: %s eligible classes, %s target cells in eval",
                    direction_key, level, plan["n_shared_classes"], plan["n_target_in_eval"])
    save_json({"direction_key": direction_key, "direction": direction, "plans": plans},
              paths["logs"] / f"{direction_key}_class_plan.json")
    return plans


def crossed_class_plan(paths, source_full, target, subset_indices, logger) -> dict:
    """Common l2/l1 class set across all five PHASE 10B subsets (fixed for all 25 cells)."""
    per_subset, common = {}, {}
    for level, col in (("l2", "cell_type_l2"), ("l1", "cell_type_l1")):
        sets = []
        for s, idx in subset_indices.items():
            plan = class_plan(source_full[idx].copy(), target, col)
            per_subset[f"{level}_subset{s}"] = plan["eligible_shared_classes"]
            sets.append(set(plan["eligible_shared_classes"]))
        common[level] = sorted(set.intersection(*sets))
        logger.info("crossed %s: %s classes common to all 5 subsets: %s",
                    level, len(common[level]), common[level])
    save_json({"common_classes": common, "per_subset": per_subset,
               "rule": "intersection of per-subset eligible shared classes; fixed for all 25 grid cells"},
              paths["logs"] / "crossed_class_plan.json")
    return common


def run_batch(jobs, paths, cfg, logger) -> tuple[int, int, list[str]]:
    """Execute jobs sequentially; one AnnData pair resident at a time."""
    done = failed = 0
    failures: list[str] = []
    by_direction: dict[str, list] = {}
    for job in jobs:
        by_direction.setdefault(job["direction_key"], []).append(job)

    for direction_key, group in by_direction.items():
        direction, source_batch, target_batch = next(
            (d, sb, tb) for dk, d, sb, tb in DIRECTIONS if dk == direction_key)
        logger.info("loading %s (%s -> %s)", direction_key, source_batch, target_batch)
        source_full, target, audit = load_source_target(source_batch, target_batch)
        save_json(audit, paths["logs"] / f"{direction_key}_data_audit.json")
        plan_path = paths["logs"] / f"{direction_key}_class_plan.json"
        if not plan_path.exists():
            write_class_plan(paths, direction_key, direction, source_full, target, logger)

        subset_cache: dict[int, object] = {}
        for job in group:
            subset_seed = job["subset_seed"]
            if subset_seed is None:
                source = source_full
            else:
                if subset_seed not in subset_cache:
                    subset_cache.clear()
                    cleanup()
                    idx = subsample_index(source_full.n_obs, SUBSAMPLE_N, subset_seed)
                    subset_cache[subset_seed] = ensure_csr(source_full[idx].copy())
                source = subset_cache[subset_seed]
            tag = experiment_tag(direction_key, job["model"], job["model_seed"], subset_seed)
            try:
                run_experiment(
                    direction_key=direction_key, direction=direction,
                    model_name=job["model"], model_seed=job["model_seed"],
                    source=source, target=target, cfg=cfg, paths=paths,
                    logger=logger, subset_seed=subset_seed,
                )
                done += 1
            except Exception as exc:  # noqa: BLE001 - record and stop the grid
                failed += 1
                failures.append(f"{tag}: {type(exc).__name__}: {exc}")
                logger.error("FAILED %s: %s\n%s", tag, exc, traceback.format_exc())
                save_json({"tag": tag, "error": f"{type(exc).__name__}: {exc}",
                           "traceback": traceback.format_exc(), "timestamp_utc": utc_now()},
                          paths["logs"] / "failures" / f"{tag}_failure.json")
                raise
            finally:
                cleanup(purge_manager_store=True)
            ram = ram_snapshot()
            if ram["available_mb"] < MIN_FREE_RAM_MB:
                raise RuntimeError(
                    f"available RAM fell to {ram['available_mb']:.0f} MB "
                    f"(floor {MIN_FREE_RAM_MB} MB) after {tag}. Stopping rather than "
                    "risking an OOM kill mid-write; completed runs are checkpointed."
                )

        del source_full, target, subset_cache
        cleanup(purge_manager_store=True)
    return done, failed, failures


def build_jobs(stage: str, seeds: list[int]) -> list[dict]:
    jobs: list[dict] = []
    if stage == "gate1":
        for model in MODELS:
            jobs.append({"direction_key": "direction_B", "model": model,
                         "model_seed": 0, "subset_seed": None})
    elif stage == "gate2":
        for model in MODELS:
            jobs.append({"direction_key": "direction_A", "model": model,
                         "model_seed": 0, "subset_seed": None})
    elif stage == "seeds":
        for direction_key, *_ in DIRECTIONS:
            for model in MODELS:
                for seed in seeds:
                    jobs.append({"direction_key": direction_key, "model": model,
                                 "model_seed": seed, "subset_seed": None})
    elif stage == "crossed":
        for subset_seed in CROSSED_SUBSET_SEEDS:
            for model in MODELS:
                for model_seed in CROSSED_MODEL_SEEDS:
                    jobs.append({"direction_key": "direction_A", "model": model,
                                 "model_seed": model_seed, "subset_seed": subset_seed})
    else:
        raise ValueError(stage)
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["gate1", "gate2", "seeds", "crossed"])
    ap.add_argument("--seeds", default="1-4")
    args = ap.parse_args()

    paths = ensure_dirs()
    (paths["logs"] / "manifests").mkdir(parents=True, exist_ok=True)
    (paths["logs"] / "failures").mkdir(parents=True, exist_ok=True)
    logger = get_logger("phase11b", paths["logs"] / f"phase11b_{args.stage}.log")
    cfg = load_config()

    seeds = parse_seeds(args.seeds) if args.stage == "seeds" else [0]
    jobs = build_jobs(args.stage, seeds)
    pending = [j for j in jobs
               if not (paths["done"] / f"{experiment_tag(j['direction_key'], j['model'], j['model_seed'], j['subset_seed'])}.json").exists()]
    logger.info("stage=%s jobs=%s pending=%s start=%s",
                args.stage, len(jobs), len(pending), utc_now())

    if args.stage == "crossed":
        source_full, target, _ = load_source_target("PBMC10k", "PBMC5k")
        idxs = {s: subsample_index(source_full.n_obs, SUBSAMPLE_N, s)
                for s in CROSSED_SUBSET_SEEDS}
        if not (paths["logs"] / "crossed_class_plan.json").exists():
            crossed_class_plan(paths, source_full, target, idxs, logger)
        for s, idx in idxs.items():
            pd.DataFrame({"cell_id": source_full.obs_names[idx].astype(str)}).to_csv(
                paths["logs"] / f"crossed_subset{s}_cells.csv", index=False)
        del source_full, target
        cleanup()

    done, failed, failures = run_batch(jobs, paths, cfg, logger)
    summary = {"stage": args.stage, "n_jobs": len(jobs), "n_completed_now": done,
               "n_failed": failed, "failures": failures, "timestamp_utc": utc_now()}
    save_json(summary, paths["logs"] / f"stage_{args.stage}_summary.json")
    logger.info("stage %s finished: %s completed, %s failed", args.stage, done, failed)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
