#!/usr/bin/env python
"""PHASE 6: independent Seurat v4 annotation and biological evaluation.

Does not train MOFA+ or modify scVI / scVI_matched / totalVI models.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiments.phase6_run import run_phase6
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--skip-train", action="store_true")
    parser.add_argument("--skip-transfer", action="store_true")
    parser.add_argument("--config", default=str(ROOT / "config" / "config.yaml"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    seed = int(cfg["reproducibility"]["default_seed"] if args.seed is None else args.seed)
    summary = run_phase6(
        config=cfg,
        seed=seed,
        skip_train=args.skip_train,
        skip_transfer=args.skip_transfer,
    )
    print("PHASE 6 finished:")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
