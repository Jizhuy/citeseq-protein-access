#!/usr/bin/env python
"""PHASE 5A–5B: matched scVI control and annotation infrastructure.

Does not train MOFA+, stress tests, or invent cell-type labels.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiments.phase5 import run_phase5_analysis
from src.models.run_scvi_matched import run_scvi_matched
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--skip-train",
        action="store_true",
        help="Reuse an existing scVI_matched model and only refresh PHASE 5 tables.",
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "config.yaml"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    seed = int(cfg["reproducibility"]["default_seed"] if args.seed is None else args.seed)
    if not args.skip_train:
        run_scvi_matched(config=cfg, seed=seed)
    summary = run_phase5_analysis(config=cfg, seed=seed)
    print("PHASE 5A–5B analysis written:")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
