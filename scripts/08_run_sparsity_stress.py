#!/usr/bin/env python
"""PHASE 8: protein sparsity / corruption stress test.

Does not run true missing-modality experiments.
Does not retrain scVI_matched.
Default is the smoke test (0/40/85%, seed 0). Pass --full only after smoke succeeds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiments.phase8_sparsity import run_phase8
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="Run 0/40/85% seed 0 only (default if neither flag is set).")
    parser.add_argument("--full", action="store_true", help="Run the full corruption grid after a successful smoke test.")
    parser.add_argument("--skip-train", action="store_true")
    parser.add_argument("--config", default=str(ROOT / "config" / "config.yaml"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    smoke = not args.full
    if args.smoke and args.full:
        raise SystemExit("Pass only one of --smoke or --full.")
    summary = run_phase8(config=cfg, smoke=smoke, skip_train=args.skip_train)
    print("PHASE 8 finished:")
    for key, value in summary.items():
        if key == "smoke_checks":
            print("  smoke_checks:")
            for row in value:
                print(f"    {row}")
        else:
            print(f"  {key}: {value}")
    if smoke and summary.get("n_failed"):
        raise SystemExit("Smoke test failed; not launching the full grid.")


if __name__ == "__main__":
    main()
