#!/usr/bin/env python
"""PHASE 11B §10-11: validate what return_dist=True returns.

Source inspection shows VAEMixin.get_latent_representation returns
(qz.loc, qz.scale.square()), i.e. the posterior VARIANCE. This script confirms
that numerically: it draws S posterior samples per target cell from the module's
own sampler and checks the empirical variance against the closed-form value.

If the second return were a scale/SD, the empirical-to-closed-form ratio would
land near sqrt(var) rather than 1.

Runs on seeds 0, 5, 10, 15 only; the full 20-seed grid does not need MC sampling.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import load_source_target  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    DIRECTIONS,
    MC_CHECK_SEEDS,
    MC_SAMPLES,
    MODELS,
    cleanup,
    ensure_dirs,
    experiment_tag,
    model_dir_for,
    utc_now,
)
from src.utils.io import save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

POSTERIOR_SAMPLING_SEED = 2026
MC_CELL_CAP = 2000  # sampling subset keeps the check cheap; cells are the first N in order


def mc_variance(model, indices: np.ndarray, n_samples: int, seed: int) -> np.ndarray:
    """Welford-free two-pass-free online variance over repeated stochastic encodings."""
    torch.manual_seed(seed)
    total = None
    total_sq = None
    for _ in range(n_samples):
        z = np.asarray(model.get_latent_representation(indices=indices, give_mean=False))
        total = z if total is None else total + z
        total_sq = z**2 if total_sq is None else total_sq + z**2
    mean = total / n_samples
    return (total_sq / n_samples) - mean**2, mean


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=MC_SAMPLES)
    ap.add_argument("--extra-s100-seed0", action="store_true",
                    help="also run S=100 on seed 0 as a sample-count sensitivity check")
    args = ap.parse_args()

    from scvi.model import SCVI, TOTALVI
    paths = ensure_dirs()
    logger = get_logger("phase11b_mc", paths["logs"] / "phase11b_mc_check.log")
    rows = []

    for direction_key, direction, source_batch, target_batch in DIRECTIONS:
        _, target, _ = load_source_target(source_batch, target_batch)
        for model_name in MODELS:
            cls = SCVI if model_name == "scvi_matched" else TOTALVI
            for seed in MC_CHECK_SEEDS:
                tag = experiment_tag(direction_key, model_name, seed, None)
                if not (paths["done"] / f"{tag}.json").exists():
                    logger.warning("skip missing %s", tag)
                    continue
                mdir = model_dir_for(paths, direction_key, model_name, seed, None)
                query = target.copy()
                cls.prepare_query_anndata(query, str(mdir / "source_model"))
                model = cls.load(str(mdir / "adapted_query_model"), adata=query,
                                 accelerator="gpu", device="auto")
                n = min(MC_CELL_CAP, query.n_obs)
                idx = np.arange(n)
                _, closed_var = model.get_latent_representation(indices=idx, return_dist=True)
                sample_sets = [args.samples]
                if args.extra_s100_seed0 and seed == 0:
                    sample_sets.append(100)
                for s in sample_sets:
                    emp_var, emp_mean = mc_variance(model, idx, s, POSTERIOR_SAMPLING_SEED + seed)
                    ratio = emp_var / closed_var
                    mu, _ = model.get_latent_representation(indices=idx, return_dist=True)
                    rows.append({
                        "direction": direction_key, "model": model_name, "seed": seed,
                        "n_samples": s, "n_cells": int(n),
                        "mean_closed_form_variance": float(closed_var.mean()),
                        "mean_empirical_variance": float(emp_var.mean()),
                        "mean_ratio_empirical_over_closed": float(ratio.mean()),
                        "median_ratio": float(np.median(ratio)),
                        "ratio_p05": float(np.quantile(ratio, 0.05)),
                        "ratio_p95": float(np.quantile(ratio, 0.95)),
                        "expected_mc_relative_sd": float(np.sqrt(2.0 / (s - 1))),
                        "corr_percell_mean_variance": float(np.corrcoef(
                            emp_var.mean(axis=1), closed_var.mean(axis=1))[0, 1]),
                        "max_abs_mean_shift": float(np.abs(emp_mean - mu).max()),
                        "ratio_if_second_return_were_sd": float(
                            (emp_var / np.sqrt(closed_var)).mean()),
                    })
                    logger.info("%s S=%s ratio=%.4f (expect 1.0 for variance; "
                                "%.4f if it were an SD)", tag, s,
                                rows[-1]["mean_ratio_empirical_over_closed"],
                                rows[-1]["ratio_if_second_return_were_sd"])
                del model, query
                cleanup()
        del target
        cleanup()

    frame = pd.DataFrame(rows)
    frame.to_csv(paths["tables"] / "phase11b_posterior_mc_check.csv", index=False)
    primary = frame[frame["n_samples"] == args.samples]
    verdict = {
        "quantity_returned_by_return_dist": "(qz.loc, qz.scale**2) = posterior mean and VARIANCE",
        "source": "scvi.model.base.VAEMixin.get_latent_representation, scvi-tools 1.3.3",
        "latent_distribution": "normal, so q(z|x) is a diagonal Gaussian and the variance is exact",
        "mc_samples": args.samples,
        "posterior_sampling_seed": POSTERIOR_SAMPLING_SEED,
        "mean_ratio_empirical_over_closed_form": float(
            primary["mean_ratio_empirical_over_closed"].mean()),
        "range_of_mean_ratio": [float(primary["mean_ratio_empirical_over_closed"].min()),
                                float(primary["mean_ratio_empirical_over_closed"].max())],
        "expected_mc_relative_sd_per_element": float(np.sqrt(2.0 / (args.samples - 1))),
        "consistent_with_variance": bool(
            primary["mean_ratio_empirical_over_closed"].between(0.9, 1.1).all()),
        "n_configurations_checked": int(len(primary)),
        "timestamp_utc": utc_now(),
    }
    save_json(verdict, paths["logs"] / "posterior_distribution_validation.json")
    logger.info("MC check verdict: %s", verdict)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
