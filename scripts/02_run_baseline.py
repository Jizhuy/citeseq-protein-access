#!/usr/bin/env python
"""Run scVI, scVI_matched, or totalVI. Does not train MOFA+."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models.run_scvi import run_scvi
from src.models.run_scvi_matched import run_scvi_matched
from src.models.run_totalvi import run_totalvi
from src.utils.io import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        default="totalvi",
        help="scvi, scvi_matched, or totalvi. Other models are blocked.",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "config.yaml"),
        help="Path to config.yaml",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model = args.model.lower()
    cfg = load_config(args.config)
    seed = int(cfg["reproducibility"]["default_seed"] if args.seed is None else args.seed)
    if model == "totalvi":
        run_totalvi(config=cfg, seed=seed)
        return
    if model == "scvi":
        run_scvi(config=cfg, seed=seed)
        return
    if model in {"scvi_matched", "scvi-matched"}:
        run_scvi_matched(config=cfg, seed=seed)
        return
    print(
        f"Supported models: scvi, scvi_matched, totalvi. Received --model {args.model}. "
        "MOFA+/PCA/sparsity/missing-modality. For MOFA+ run scripts/07_run_mofa.py."
    )
    sys.exit(2)


if __name__ == "__main__":
    main()
