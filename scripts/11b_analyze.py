#!/usr/bin/env python
"""PHASE 11B Parts I-II and V: 20-seed uncertainty stability.

Stages
  perseed    score every run, persist target probabilities, primary +
             selective-prediction tables
  delta      paired totalVI-scVI seed deltas, seed-level and target-cell
             bootstraps, evidence classification
  ensemble   20-model ensemble, ensemble-size convergence curve
  retention  cell-type retention under uncertainty-based rejection
  repro      seeds 0-4 versus corrected PHASE 10

Trains nothing. Reads only PHASE 11B artifacts (plus corrected PHASE 10 for repro).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import (  # noqa: E402
    align_probability_matrices,
    bootstrap_indices,
    calibration_metrics,
    ci_bounds,
    classification_metrics,
    coverage_table,
    ensemble_uncertainty,
    error_detection,
    phase10_reference_metrics,
    predictive_uncertainty,
    risk_coverage,
)
from src.experiments.phase11b_analysis import (  # noqa: E402
    BOOTSTRAP_METRICS,
    EVIDENCE_RULE,
    HIGHER_IS_BETTER,
    DirectionEval,
    classify_evidence,
    load_probabilities,
    load_run,
    metric_block,
    one_hot_matrix,
    save_probabilities,
    score_run,
)
from src.experiments.phase11b_robustness import (  # noqa: E402
    BOOTSTRAP_SEED,
    COVERAGE_GRID,
    DIRECTIONS,
    MODELS,
    N_BOOT_CELLS,
    N_BOOT_SEEDS,
    SEEDS,
    cleanup,
    ensure_dirs,
    experiment_tag,
    utc_now,
)
from src.utils.io import save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

DIRECTION_KEYS = tuple(d[0] for d in DIRECTIONS)
ENSEMBLE_SIZES = (1, 2, 3, 5, 10, 15, 20)
RANDOM_SUBSET_SIZES = (5, 10, 15)
N_RANDOM_DRAWS = 5
RANDOM_SUBSET_SEED = 11_202


def available_seeds(direction_key: str, model: str, paths) -> list[int]:
    return [s for s in SEEDS
            if (paths["done"] / f"{experiment_tag(direction_key, model, s, None)}.json").exists()]


def read_manifest(paths, tag: str) -> dict:
    path = paths["logs"] / "manifests" / f"{tag}_manifest.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


# ---------------------------------------------------------------------------
def stage_perseed(paths, logger) -> None:
    primary_rows, coverage_rows, drift_rows = [], [], []
    for direction_key in DIRECTION_KEYS:
        data = DirectionEval(direction_key)
        sl = data.slice("l2")
        logger.info("%s: %s eval cells, %s classes", direction_key,
                    sl["y_target"].size, len(sl["classes"]))
        for model in MODELS:
            for seed in available_seeds(direction_key, model, paths):
                tag = experiment_tag(direction_key, model, seed, None)
                run = load_run(direction_key, model, seed)
                if not np.array_equal(run["target_cells"], data.target_cells):
                    raise RuntimeError(f"{tag}: target cell order mismatch")
                if not np.array_equal(run["source_cells"], data.source_cells):
                    raise RuntimeError(f"{tag}: source cell order mismatch")
                z_src = run["source_latent"][sl["source_index"]]
                z_tgt = run["target_latent"][sl["target_index"]]
                u_lat = run["u_latent_target"][sl["target_index"]]
                res = score_run(z_src, sl["y_source"], z_tgt, sl["y_target"], u_lat, seed)
                save_probabilities(direction_key, tag, res["probs"], res["classes"],
                                   sl["target_cells"], u_lat,
                                   {"direction": direction_key, "model": model,
                                    "seed": seed, "level": "l2"})
                man = read_manifest(paths, tag)
                primary_rows.append({
                    "direction": direction_key, "model": model, "seed": seed,
                    **res["row"],
                    "source_training_epochs": man.get("source_training_epochs"),
                    "query_epochs": man.get("query_epochs"),
                    "runtime": man.get("total_runtime_seconds"),
                    "gpu_peak_mb": man.get("gpu_peak_allocated_mb"),
                    "rss_peak_mb": (man.get("ram_after") or {}).get("process_rss_mb"),
                })
                for cov in res["coverage_rows"]:
                    coverage_rows.append({"direction": direction_key, "model": model,
                                          "seed": seed, **cov,
                                          "aurc_component_if_relevant": res["row"]["aurc"]})
                drift = man.get("source_latent_drift_pre_to_post_adaptation", {})
                drift_rows.append({"direction": direction_key, "model": model, "seed": seed,
                                   **drift,
                                   "integrity_max_latent_diff":
                                       (man.get("artifact_integrity") or {}).get("max_latent_diff"),
                                   "integrity_max_variance_diff":
                                       (man.get("artifact_integrity") or {}).get("max_variance_diff")})
                logger.info("%s macroF1=%.4f errAUROC=%.4f latAUROC=%.4f ECE=%.4f",
                            tag, res["row"]["macro_f1"], res["row"]["error_auroc_predictive"],
                            res["row"]["error_auroc_latent"], res["row"]["ece_equal_frequency"])
                del run, res
                cleanup()
        del data
        cleanup()

    pd.DataFrame(primary_rows).to_csv(paths["tables"] / "phase11b_20seed_primary.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(paths["tables"] / "phase11b_selective_prediction.csv", index=False)
    pd.DataFrame(drift_rows).to_csv(paths["tables"] / "phase11b_latent_drift.csv", index=False)

    frame = pd.DataFrame(primary_rows)
    summary = (frame.groupby(["direction", "model"])
               .agg(["mean", "std"])[[
                   "macro_f1", "accuracy", "error_auroc_predictive", "error_auprc_predictive",
                   "nll", "brier", "ece_equal_frequency", "aurc", "error_auroc_latent"]]
               .round(6))
    summary.to_csv(paths["tables"] / "phase11b_20seed_summary.csv")
    logger.info("perseed done: %s rows", len(primary_rows))


# ---------------------------------------------------------------------------
def stage_delta(paths, logger) -> None:
    primary = pd.read_csv(paths["tables"] / "phase11b_20seed_primary.csv")
    metrics = ["macro_f1", "accuracy", "error_auroc_predictive", "error_auprc_predictive",
               "nll", "brier", "ece_equal_frequency", "ece_equal_width", "aurc",
               "error_auroc_latent", "error_auprc_latent"]
    rows, summary_rows = [], []
    rng_master = np.random.default_rng(BOOTSTRAP_SEED)

    for direction_key in DIRECTION_KEYS:
        sub = primary[primary["direction"] == direction_key]
        seeds = sorted(set(sub[sub["model"] == "scvi_matched"]["seed"])
                       & set(sub[sub["model"] == "totalvi"]["seed"]))
        for metric in metrics:
            deltas = []
            for seed in seeds:
                s = float(sub[(sub["model"] == "scvi_matched") & (sub["seed"] == seed)][metric].iloc[0])
                t = float(sub[(sub["model"] == "totalvi") & (sub["seed"] == seed)][metric].iloc[0])
                deltas.append(t - s)
                rows.append({"direction": direction_key, "seed": seed, "metric": metric,
                             "scvi": s, "totalvi": t, "delta": t - s})
            d = np.asarray(deltas, dtype=float)
            finite = d[np.isfinite(d)]
            # Paired bootstrap over SEED PAIRS (§20). Kept separate from cell bootstrap.
            rng = np.random.default_rng(rng_master.integers(0, 2**31))
            boots = np.array([np.mean(finite[rng.integers(0, finite.size, finite.size)])
                              for _ in range(N_BOOT_SEEDS)]) if finite.size else np.array([])
            lo, hi = ci_bounds(boots) if boots.size else (float("nan"), float("nan"))
            summary_rows.append({
                "direction": direction_key, "metric": metric,
                "higher_is_better": HIGHER_IS_BETTER.get(metric),
                "n_seed_pairs": int(finite.size),
                "mean": float(np.mean(finite)) if finite.size else float("nan"),
                "sd": float(np.std(finite, ddof=1)) if finite.size > 1 else float("nan"),
                "median": float(np.median(finite)) if finite.size else float("nan"),
                "positive_count": int(np.sum(finite > 0)),
                "negative_count": int(np.sum(finite < 0)),
                "seed_bootstrap_ci_low": lo, "seed_bootstrap_ci_high": hi,
            })

    # ---- per-seed target-cell bootstrap, identical indices for both models (§21) ----
    cell_rows = []
    for direction_key in DIRECTION_KEYS:
        data = DirectionEval(direction_key)
        y = data.slice("l2")["y_target"]
        for seed in SEEDS:
            tags = {m: experiment_tag(direction_key, m, seed, None) for m in MODELS}
            if not all((paths["done"] / f"{t}.json").exists() for t in tags.values()):
                continue
            packs = {m: load_probabilities(direction_key, tags[m]) for m in MODELS}
            ref = list(packs["scvi_matched"]["classes"])
            if list(packs["totalvi"]["classes"]) != ref:
                raise RuntimeError(f"{direction_key} seed {seed}: class order differs")
            idxs = bootstrap_indices(y.size, N_BOOT_CELLS, BOOTSTRAP_SEED + seed)
            onehot = one_hot_matrix(y, packs["scvi_matched"]["classes"])
            draws = {m: {k: [] for k in BOOTSTRAP_METRICS} for m in MODELS}
            for b in range(N_BOOT_CELLS):
                for m in MODELS:
                    blk = metric_block(packs[m]["probs"], packs[m]["classes"], y,
                                       idxs[b], onehot)
                    for k in BOOTSTRAP_METRICS:
                        draws[m][k].append(blk[k])
            for k in BOOTSTRAP_METRICS:
                d = np.asarray(draws["totalvi"][k]) - np.asarray(draws["scvi_matched"][k])
                lo, hi = ci_bounds(d)
                cell_rows.append({"direction": direction_key, "seed": seed, "metric": k,
                                  "delta_mean_bootstrap": float(np.nanmean(d)),
                                  "cell_bootstrap_ci_low": lo, "cell_bootstrap_ci_high": hi,
                                  "excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi)
                                                        and (lo > 0 or hi < 0))})
            logger.info("%s seed %s target-cell bootstrap done", direction_key, seed)
            del packs
            cleanup()
        del data
        cleanup()

    cell_frame = pd.DataFrame(cell_rows)
    cell_frame.to_csv(paths["tables"] / "phase11b_cell_bootstrap.csv", index=False)

    summary = pd.DataFrame(summary_rows)
    evidence = []
    for _, r in summary.iterrows():
        sel = cell_frame[(cell_frame["direction"] == r["direction"])
                         & (cell_frame["metric"] == r["metric"])]
        frac = float(sel["excludes_zero"].mean()) if len(sel) else 0.0
        evidence.append({
            "frac_seed_cell_ci_excluding_zero": frac,
            "evidence": classify_evidence(
                (r["seed_bootstrap_ci_low"], r["seed_bootstrap_ci_high"]),
                int(r["positive_count"]), int(r["n_seed_pairs"]), frac),
        })
    summary = pd.concat([summary, pd.DataFrame(evidence)], axis=1)

    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_seed_delta.csv", index=False)
    summary.to_csv(paths["tables"] / "phase11b_seed_delta_summary.csv", index=False)
    save_json(EVIDENCE_RULE, paths["logs"] / "evidence_rule.json")
    logger.info("delta done")


# ---------------------------------------------------------------------------
def stage_ensemble(paths, logger) -> None:
    conv_rows, ens_rows, ens_cell = [], [], {}
    for direction_key in DIRECTION_KEYS:
        data = DirectionEval(direction_key)
        sl = data.slice("l2")
        y = sl["y_target"]
        for model in MODELS:
            seeds = available_seeds(direction_key, model, paths)
            if not seeds:
                continue
            per_seed = {s: load_probabilities(direction_key,
                                              experiment_tag(direction_key, model, s, None))
                        for s in seeds}
            stacked, classes = align_probability_matrices(
                {s: {"classes": per_seed[s]["classes"], "probs": per_seed[s]["probs"]}
                 for s in seeds})
            seed_list = sorted(per_seed)

            def evaluate(sel_idx: list[int], label: str, size: int) -> dict:
                sub = stacked[sel_idx]
                ens = ensemble_uncertainty(sub, classes)
                p_bar, pred = ens["p_bar"], ens["ensemble_pred"]
                correct = (pred == y)
                err = (~correct).astype(int)
                u = 1.0 - p_bar.max(axis=1)
                det = error_detection(u, err)
                cal = calibration_metrics(p_bar, y, classes, pred)
                rc = risk_coverage(u, correct.astype(float))
                return {"direction": direction_key, "model": model, "ensemble_size": size,
                        "seed_subset_definition": label,
                        **classification_metrics(y, pred),
                        "error_auroc": det["auroc"], "error_auprc": det["auprc"],
                        "nll": cal["nll"], "brier": cal["brier"],
                        "ece": cal["ece_equal_frequency"], "aurc": rc["aurc"],
                        "median_predictive_entropy": float(np.median(ens["predictive_entropy"])),
                        "median_expected_entropy": float(np.median(ens["expected_entropy"])),
                        "median_ensemble_disagreement": float(np.median(ens["disagreement"])),
                        "median_variation_ratio": float(np.median(ens["variation_ratio"])),
                        "_ens": ens}

            for size in ENSEMBLE_SIZES:
                if size > len(seed_list):
                    continue
                sel = list(range(size))
                out = evaluate(sel, f"nested seeds {seed_list[0]}-{seed_list[size-1]}", size)
                ens_full = out.pop("_ens")
                conv_rows.append(out)
                if size == len(seed_list):
                    ens_rows.append({**out, "n_cells": int(y.size)})
                    ens_cell[(direction_key, model)] = {
                        "disagreement": ens_full["disagreement"],
                        "variation_ratio": ens_full["variation_ratio"],
                        "predictive_entropy": ens_full["predictive_entropy"],
                        "cell_ids": sl["target_cells"],
                    }

            rng = np.random.default_rng(RANDOM_SUBSET_SEED)
            for size in RANDOM_SUBSET_SIZES:
                if size >= len(seed_list):
                    continue
                for draw in range(N_RANDOM_DRAWS):
                    sel = sorted(rng.choice(len(seed_list), size=size, replace=False).tolist())
                    out = evaluate(sel, f"random draw {draw}: seeds "
                                        f"{[seed_list[i] for i in sel]}", size)
                    out.pop("_ens")
                    conv_rows.append(out)
            logger.info("%s %s ensemble over %s seeds done", direction_key, model, len(seed_list))
            del per_seed, stacked
            cleanup()
        del data
        cleanup()

    pd.DataFrame(conv_rows).to_csv(paths["tables"] / "phase11b_ensemble_convergence.csv", index=False)
    pd.DataFrame(ens_rows).to_csv(paths["tables"] / "phase11b_ensemble.csv", index=False)
    np.savez_compressed(
        paths["tables"] / "phase11b_ensemble_percell.npz",
        **{f"{d}__{m}__{k}": v for (d, m), pack in ens_cell.items()
           for k, v in pack.items()})
    logger.info("ensemble done")


# ---------------------------------------------------------------------------
def stage_retention(paths, logger) -> None:
    rows, imbalance = [], []
    for direction_key in DIRECTION_KEYS:
        data = DirectionEval(direction_key)
        sl = data.slice("l2")
        y = sl["y_target"]
        for model in MODELS:
            seeds = available_seeds(direction_key, model, paths)
            if not seeds:
                continue
            # Average retention over seeds so the pattern is not a single-seed artifact.
            per_cov: dict[float, list[pd.Series]] = {c: [] for c in COVERAGE_GRID}
            for seed in seeds:
                pack = load_probabilities(direction_key,
                                          experiment_tag(direction_key, model, seed, None))
                u = 1.0 - pack["probs"].max(axis=1)
                order = np.argsort(u, kind="mergesort")
                for cov in COVERAGE_GRID:
                    keep = max(1, int(round(cov * y.size)))
                    per_cov[cov].append(pd.Series(y[order[:keep]]).value_counts())
                del pack
            original = pd.Series(y).value_counts()
            for cov in COVERAGE_GRID:
                mean_counts = (pd.concat(per_cov[cov], axis=1).fillna(0.0).mean(axis=1))
                rates = []
                for ct in original.index:
                    retained = float(mean_counts.get(ct, 0.0))
                    rate = retained / float(original[ct])
                    rates.append(rate)
                    rows.append({"direction": direction_key, "model": model,
                                 "coverage": cov, "cell_type": ct,
                                 "original_n": int(original[ct]),
                                 "retained_n": retained, "retention_rate": rate})
                arr = np.asarray(rates)
                imbalance.append({"direction": direction_key, "model": model, "coverage": cov,
                                  "max_retention_rate": float(arr.max()),
                                  "min_retention_rate": float(arr.min()),
                                  "retention_range": float(arr.max() - arr.min()),
                                  "retention_cv": float(arr.std(ddof=1) / arr.mean())
                                  if arr.mean() > 0 else float("nan"),
                                  "n_cell_types": int(arr.size)})
            logger.info("%s %s retention done", direction_key, model)
        del data
        cleanup()
    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_celltype_retention.csv", index=False)
    pd.DataFrame(imbalance).to_csv(paths["tables"] / "phase11b_retention_imbalance.csv", index=False)
    logger.info("retention done")


# ---------------------------------------------------------------------------
def stage_repro(paths, logger) -> None:
    """Seeds 0-4 versus corrected PHASE 10. Independently retrained: no bitwise identity."""
    primary = pd.read_csv(paths["tables"] / "phase11b_20seed_primary.csv")
    ref = phase10_reference_metrics(corrected=True)
    # The corrected table already holds one common-class row per
    # (direction, model, seed, level, classifier) and uses lowercase model names.
    dir_map = {"direction_A": "pbmc10k_to_pbmc5k", "direction_B": "pbmc5k_to_pbmc10k"}
    rows = []
    for direction_key, p10_dir in dir_map.items():
        for model in MODELS:
            for seed in range(5):
                new = primary[(primary["direction"] == direction_key)
                              & (primary["model"] == model) & (primary["seed"] == seed)]
                old = ref[(ref["direction"] == p10_dir) & (ref["model"] == model)
                          & (ref["model_seed"] == seed) & (ref["label_level"] == "l2")]
                if len(old) > 1:
                    raise RuntimeError(
                        f"{p10_dir}/{model}/seed{seed}: {len(old)} PHASE 10 reference rows; "
                        "expected exactly one")
                if new.empty or old.empty:
                    continue
                for metric, col in (("macro_f1", "macro_f1"), ("accuracy", "accuracy")):
                    rows.append({"direction": direction_key, "model": model, "seed": seed,
                                 "metric": metric,
                                 "phase10_corrected": float(old[col].iloc[0]),
                                 "phase11b": float(new[metric].iloc[0]),
                                 "difference": float(new[metric].iloc[0]) - float(old[col].iloc[0])})
    frame = pd.DataFrame(rows)
    frame.to_csv(paths["tables"] / "phase11b_phase10_replication.csv", index=False)
    if not frame.empty:
        summ = (frame.groupby(["direction", "model", "metric"])["difference"]
                .agg(["mean", "std", "min", "max"]).round(6))
        summ.to_csv(paths["tables"] / "phase11b_phase10_replication_summary.csv")
        logger.info("replication difference summary:\n%s", summ.to_string())
    logger.info("repro done")


STAGES = {"perseed": stage_perseed, "delta": stage_delta, "ensemble": stage_ensemble,
          "retention": stage_retention, "repro": stage_repro}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=[*STAGES, "all"])
    args = ap.parse_args()
    paths = ensure_dirs()
    logger = get_logger("phase11b_analysis", paths["logs"] / "phase11b_analysis.log")
    todo = list(STAGES) if args.stage == "all" else [args.stage]
    for name in todo:
        logger.info("=== stage %s (%s) ===", name, utc_now())
        STAGES[name](paths, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
