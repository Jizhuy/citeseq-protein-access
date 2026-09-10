#!/usr/bin/env python
"""PHASE 12A Part 1: train the second replicate of the crossed design.

Direction A only. 5 PHASE 10B source subsets x 5 model seeds x 2 models = 50
new runs, each holding the source subset, model seed, architecture, data,
labels and hyperparameters fixed and redrawing only run-level stochasticity
via an explicit `run_replicate_seed` (see phase12a_replication and
logs/replicate_design.md).

PHASE 11B artifacts are never touched: everything lands under
results/phase12a_variance_replication/. Resumable via done markers.

Usage
  python scripts/12a_run_replicates.py
  python scripts/12a_run_replicates.py --subsets 0,1     # partial
"""

from __future__ import annotations

import argparse
import shutil
import sys
import traceback
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import ensure_csr, load_source_target  # noqa: E402
from src.experiments.phase10b_source_size import SUBSAMPLE_N, subsample_index  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    cleanup,
    experiment_tag,
    phase11b_paths,
    ram_snapshot,
    run_experiment,
    utc_now,
)
from src.experiments.phase12a_replication import (  # noqa: E402
    CROSSED_SOURCE_BATCH,
    CROSSED_TARGET_BATCH,
    crossed_jobs,
    ensure_dirs,
)
from src.utils.io import load_config, save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

MIN_FREE_RAM_MB = 500.0


def _mirror_class_plan(paths, logger) -> None:
    """Reuse the PHASE 11B common-class definition verbatim.

    The replicate must score exactly the same target cells as replicate 0, so
    the class plan is copied rather than recomputed.
    """
    src = phase11b_paths()["logs"] / "crossed_class_plan.json"
    if not src.exists():
        raise FileNotFoundError(
            f"PHASE 11B crossed class plan missing at {src}; PHASE 12A must reuse it")
    dst = paths["logs"] / "crossed_class_plan.json"
    if not dst.exists():
        shutil.copy2(src, dst)
        logger.info("copied crossed class plan from PHASE 11B: %s", src)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subsets", default="", help="comma-separated subset seeds; default all")
    args = ap.parse_args()

    paths = ensure_dirs()
    (paths["logs"] / "failures").mkdir(parents=True, exist_ok=True)
    logger = get_logger("phase12a", paths["logs"] / "phase12a_replicates.log")
    cfg = load_config()
    _mirror_class_plan(paths, logger)

    jobs = crossed_jobs()
    if args.subsets:
        wanted = {int(x) for x in args.subsets.split(",") if x.strip()}
        jobs = [j for j in jobs if j["subset_seed"] in wanted]
    pending = [j for j in jobs
               if not (paths["done"] / f"{experiment_tag(j['direction_key'], j['model'], j['model_seed'], j['subset_seed'], j['replicate'])}.json").exists()]
    logger.info("PHASE 12A: %s jobs, %s pending, start %s", len(jobs), len(pending), utc_now())
    if not pending:
        logger.info("nothing to do")
        return 0

    logger.info("loading %s -> %s", CROSSED_SOURCE_BATCH, CROSSED_TARGET_BATCH)
    source_full, target, audit = load_source_target(CROSSED_SOURCE_BATCH, CROSSED_TARGET_BATCH)
    save_json(audit, paths["logs"] / "direction_A_data_audit.json")

    # Record the exact subset membership so replicate 1 is provably the same
    # source cells as replicate 0.
    for s in sorted({j["subset_seed"] for j in jobs}):
        idx = subsample_index(source_full.n_obs, SUBSAMPLE_N, s)
        pd.DataFrame({"cell_id": source_full.obs_names[idx].astype(str)}).to_csv(
            paths["logs"] / f"crossed_subset{s}_cells.csv", index=False)

    done = failed = 0
    failures: list[str] = []
    subset_cache: dict[int, object] = {}
    # Group by subset so each subsampled AnnData is built once.
    pending.sort(key=lambda j: (j["subset_seed"], j["model_seed"], j["model"]))

    for job in pending:
        s = job["subset_seed"]
        if s not in subset_cache:
            subset_cache.clear()
            cleanup(purge_manager_store=True)
            idx = subsample_index(source_full.n_obs, SUBSAMPLE_N, s)
            subset_cache[s] = ensure_csr(source_full[idx].copy())
            logger.info("built source subset %s (n=%s)", s, subset_cache[s].n_obs)
        source = subset_cache[s]
        tag = experiment_tag(job["direction_key"], job["model"], job["model_seed"],
                             s, job["replicate"])
        try:
            run_experiment(
                direction_key=job["direction_key"], direction=job["direction"],
                model_name=job["model"], model_seed=job["model_seed"],
                source=source, target=target, cfg=cfg, paths=paths, logger=logger,
                subset_seed=s, replicate=job["replicate"],
                run_replicate_seed=job["run_replicate_seed"],
            )
            done += 1
        except Exception as exc:  # noqa: BLE001
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
                f"available RAM fell to {ram['available_mb']:.0f} MB (floor "
                f"{MIN_FREE_RAM_MB} MB) after {tag}. Stopping rather than risking "
                "an OOM kill mid-write; completed runs are checkpointed."
            )

    del source_full, target, subset_cache
    cleanup(purge_manager_store=True)
    save_json({"phase": "12A", "n_jobs": len(jobs), "n_completed_now": done,
               "n_failed": failed, "failures": failures, "timestamp_utc": utc_now()},
              paths["logs"] / "stage_replicates_summary.json")
    logger.info("PHASE 12A finished: %s completed, %s failed", done, failed)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
