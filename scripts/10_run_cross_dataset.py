#!/usr/bin/env python
"""PHASE 10: cross-dataset transfer / query generalization.

Default: seed-0 gate, both directions.
Pass --full only after the gate succeeds.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
Path(os.environ["TMPDIR"]).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLBACKEND", "Agg")

from src.experiments.phase10_cross_dataset import DIRECTIONS, run_phase10
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", action="store_true", help="Seed 0, both directions (default).")
    parser.add_argument("--full", action="store_true", help="Seeds 1-4 after a successful gate.")
    parser.add_argument(
        "--direction",
        choices=["pbmc10k_to_pbmc5k", "pbmc5k_to_pbmc10k", "both"],
        default="both",
    )
    parser.add_argument("--config", default=str(ROOT / "config" / "config.yaml"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    if args.full:
        seeds = [1, 2, 3, 4]
    else:
        seeds = [0]
    if args.direction == "both":
        directions = DIRECTIONS
    else:
        directions = tuple(d for d in DIRECTIONS if d[0] == args.direction)
    summary = run_phase10(seeds=seeds, directions=directions, config=cfg)
    print("PHASE 10 finished:")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
