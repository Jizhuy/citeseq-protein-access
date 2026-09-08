#!/usr/bin/env python
"""PHASE 10B: source-size sensitivity control for PHASE 10 Direction A.

Default is the gate: subsample_seed 0 only.
Pass --full to run subsample seeds 0-4 after the gate succeeds.
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

from src.experiments.phase10b_source_size import SUBSAMPLE_SEEDS, run_phase10b
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="Subsample seeds 0-4.")
    parser.add_argument("--config", default=str(ROOT / "config" / "config.yaml"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    seeds = SUBSAMPLE_SEEDS if args.full else (0,)
    summary = run_phase10b(subsample_seeds=seeds, config=cfg)
    print("PHASE 10B finished:")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
