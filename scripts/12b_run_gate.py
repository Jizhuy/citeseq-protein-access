#!/usr/bin/env python
"""PHASE 12B gated runner.

Usage:
  python -m scripts.12b_run_gate --gate 1
  python -m scripts.12b_run_gate --gate 2
  python -m scripts.12b_run_gate --gate 3
  python -m scripts.12b_run_gate --gate all_seeds
  python -m scripts.12b_run_gate --fold 0 --model scvi_matched --seed 0
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Prefer WSL CUDA libs and private TMPDIR before importing torch/scvi.
os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
os.environ["PATH"] = "/usr/lib/wsl/lib:" + os.environ.get("PATH", "")
if "LD_LIBRARY_PATH" in os.environ:
    os.environ["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + os.environ["LD_LIBRARY_PATH"]
else:
    os.environ["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib"

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gate", type=str, default=None, help="1|2|3|all_seeds")
    parser.add_argument("--fold", type=int, default=None)
    parser.add_argument("--model", type=str, default=None, choices=["scvi_matched", "totalvi"])
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    from src.experiments.phase12b_external import (
        FOLDS,
        MODELS,
        SEEDS,
        ensure_dirs,
        run_experiment,
        run_gate,
    )
    from src.utils.logging_utils import get_logger

    paths = ensure_dirs()
    logger = get_logger("phase12b_gate", paths["logs"] / "phase12b_gate.log")
    Path(os.environ["TMPDIR"]).mkdir(parents=True, exist_ok=True)

    if args.fold is not None:
        assert args.model is not None and args.seed is not None
        jobs = [(args.fold, args.model, args.seed)]
    elif args.gate == "1":
        jobs = [(0, "scvi_matched", 0), (0, "totalvi", 0)]
    elif args.gate == "2":
        jobs = [(1, "scvi_matched", 0), (1, "totalvi", 0)]
    elif args.gate == "3":
        jobs = [(f, m, 0) for f in FOLDS for m in MODELS]
    elif args.gate == "all_seeds":
        jobs = [(f, m, s) for f in FOLDS for m in MODELS for s in SEEDS]
    else:
        raise SystemExit("Provide --gate 1|2|3|all_seeds or --fold/--model/--seed")

    logger.info("jobs=%s", jobs)
    run_gate(jobs, logger=logger)
    logger.info("gate complete")


if __name__ == "__main__":
    main()
