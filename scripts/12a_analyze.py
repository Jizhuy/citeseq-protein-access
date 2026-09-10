#!/usr/bin/env python
"""PHASE 12A Part 1 analysis: separate interaction from pure residual.

Scores the full replicated crossed design (5 source subsets x 5 model seeds x
2 replicates x 2 models = 100 runs), then fits a balanced two-factor
random-effects model WITH replication so that the subset x seed interaction and
the pure run-level residual are separately identified for the first time.

Replicate 0 is read from the PHASE 11B artifact tree and replicate 1 from the
PHASE 12A tree; nothing under PHASE 11B is modified. Replicate 0 is rescored
from scratch and checked against the stored PHASE 11B table, so any drift in
the scoring path would be caught rather than silently absorbed.

Stages
  score     score all 100 runs on the PHASE 11B common class set
  variance  metric-level replicated variance components + Delta decomposition
  percell   per-cell and per-cell-type replicated decomposition
  all       every stage in order
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import load_source_target  # noqa: E402
from src.experiments.phase10b_source_size import SUBSAMPLE_N, subsample_index  # noqa: E402
from src.experiments.phase11_uncertainty import (  # noqa: E402
    calibration_metrics,
    classification_metrics,
    error_detection,
    fit_transfer_classifier,
    risk_coverage,
)
from src.experiments.phase11b_robustness import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    MODELS,
    cleanup,
    experiment_tag,
    phase11b_paths,
    utc_now,
)
from src.experiments.phase12a_replication import (  # noqa: E402
    CROSSED_DIRECTION_KEY,
    REPLICATES,
    ensure_dirs,
    run_replicate_seed,
)
from src.experiments.phase12a_stats import (  # noqa: E402
    TRUNCATION_BIAS_NOTE,
    VARIANCE_MODEL_NOTE,
    replicated_components,
    replicated_components_vectorised,
)
from src.utils.io import save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

METRICS = ("macro_f1", "error_auroc", "error_auprc", "nll", "brier", "ece", "aurc")
A, B, N = len(CROSSED_SUBSET_SEEDS), len(CROSSED_MODEL_SEEDS), len(REPLICATES)


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def _tree(paths, replicate: int) -> dict:
    """Replicate 0 lives in the PHASE 11B tree, replicate 1 in the 12A tree."""
    return phase11b_paths() if replicate == 0 else paths


def _rep_arg(replicate: int) -> int | None:
    """PHASE 11B tags carry no replicate suffix; PHASE 12A tags do."""
    return None if replicate == 0 else replicate


def load_replicate_run(paths, model: str, model_seed: int, subset_seed: int,
                       replicate: int) -> dict:
    tree = _tree(paths, replicate)
    tag = experiment_tag(CROSSED_DIRECTION_KEY, model, model_seed, subset_seed,
                         _rep_arg(replicate))
    emb = tree["embeddings"] / CROSSED_DIRECTION_KEY
    post = tree["posterior"] / CROSSED_DIRECTION_KEY
    z_src = np.load(emb / f"{tag}_source_latent_post.npy")
    z_tgt = np.load(emb / f"{tag}_target_latent_post.npy")
    with np.load(post / f"{tag}_posterior.npz", allow_pickle=True) as npz:
        u_latent = npz["u_latent_target"]
    if not np.isfinite(z_src).all() or not np.isfinite(z_tgt).all():
        raise RuntimeError(f"{tag}: non-finite latent")
    return {"tag": tag, "source_latent": z_src, "target_latent": z_tgt,
            "u_latent_target": u_latent}


def crossed_setup(paths, logger) -> dict:
    """Exactly the PHASE 11B evaluation set: same classes, same target cells."""
    plan_path = phase11b_paths()["logs"] / "crossed_class_plan.json"
    common = json.loads(plan_path.read_text(encoding="utf-8"))["common_classes"]["l2"]
    source_full, target, _ = load_source_target("PBMC10k", "PBMC5k")
    y_src_all = source_full.obs["cell_type_l2"].astype(str).to_numpy()
    high_src = source_full.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    y_tgt_all = target.obs["cell_type_l2"].astype(str).to_numpy()
    high_tgt = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    tgt_mask = high_tgt & np.isin(y_tgt_all, common)
    setup = {
        "common_classes": common,
        "target_index": np.where(tgt_mask)[0],
        "y_target": y_tgt_all[tgt_mask],
        "target_cells": target.obs_names.astype(str).to_numpy()[tgt_mask],
        "source_labels_full": y_src_all,
        "source_high_full": high_src,
        "n_source_full": int(source_full.n_obs),
    }
    logger.info("setup: %s common l2 classes, %s fixed target cells",
                len(common), setup["target_index"].size)
    del source_full, target
    cleanup()
    return setup


# ---------------------------------------------------------------------------
# stage: score
# ---------------------------------------------------------------------------
def stage_score(paths, logger) -> None:
    setup = crossed_setup(paths, logger)
    common = setup["common_classes"]
    rows = []
    for subset_seed in CROSSED_SUBSET_SEEDS:
        idx = subsample_index(setup["n_source_full"], SUBSAMPLE_N, subset_seed)
        y_sub = setup["source_labels_full"][idx]
        high_sub = setup["source_high_full"][idx]
        src_keep = np.where(high_sub & np.isin(y_sub, common))[0]
        y_source = y_sub[src_keep]
        for model in MODELS:
            for model_seed in CROSSED_MODEL_SEEDS:
                for replicate in REPLICATES:
                    tree = _tree(paths, replicate)
                    tag = experiment_tag(CROSSED_DIRECTION_KEY, model, model_seed,
                                         subset_seed, _rep_arg(replicate))
                    if not (tree["done"] / f"{tag}.json").exists():
                        raise FileNotFoundError(
                            f"missing run {tag}; the replicated design must be complete "
                            "before variance components are estimated")
                    run = load_replicate_run(paths, model, model_seed, subset_seed, replicate)
                    z_src = run["source_latent"][src_keep]
                    z_tgt = run["target_latent"][setup["target_index"]]
                    u_lat = run["u_latent_target"][setup["target_index"]]
                    y = setup["y_target"]
                    # Classifier seed follows PHASE 11B (model_seed) so replicate 0
                    # reproduces exactly; the replicate differs through its latents.
                    fit = fit_transfer_classifier(z_src, y_source, z_tgt, model_seed)
                    probs, classes, pred = fit["probs"], fit["classes"], fit["pred"]
                    if list(classes) != sorted(common):
                        raise RuntimeError(f"{tag}: classifier classes != common set")
                    correct = pred == y
                    u = 1.0 - probs.max(axis=1)
                    det = error_detection(u, (~correct).astype(int))
                    cal = calibration_metrics(probs, y, classes, pred)
                    rc = risk_coverage(u, correct.astype(float))
                    np.savez_compressed(
                        paths["probabilities"] / f"rep_{tag}_probs.npz",
                        probs=probs.astype(np.float32),
                        classes=np.asarray(classes, dtype=object),
                        cell_ids=np.asarray(setup["target_cells"], dtype=object),
                        u_pred=u.astype(np.float32), u_latent=u_lat.astype(np.float32),
                        true_class_probability=cal["true_class_probability"].astype(np.float32),
                        pred=np.asarray(pred, dtype=object))
                    rows.append({
                        "model": model, "source_subset_seed": subset_seed,
                        "model_seed": model_seed, "replicate": replicate,
                        "run_replicate_seed": run_replicate_seed(model_seed, subset_seed, replicate),
                        **classification_metrics(y, pred),
                        "error_auroc": det["auroc"], "error_auprc": det["auprc"],
                        "error_prevalence": det["error_prevalence"],
                        "nll": cal["nll"], "brier": cal["brier"],
                        "ece": cal["ece_equal_frequency"], "aurc": rc["aurc"],
                        "n_source_eval": int(y_source.size), "n_target_eval": int(y.size)})
                    logger.info("%s macroF1=%.4f AUROC=%.4f", tag,
                                rows[-1]["macro_f1"], rows[-1]["error_auroc"])
                    del run, fit, probs
                    cleanup()
    frame = pd.DataFrame(rows)
    frame.to_csv(paths["tables"] / "replicated_crossed_design.csv", index=False)

    # Replicate 0 must reproduce the stored PHASE 11B table exactly.
    old = pd.read_csv(phase11b_paths()["tables"] / "phase11b_crossed_design.csv")
    new0 = frame[frame.replicate == 0]
    merged = old.merge(new0, on=["model", "source_subset_seed", "model_seed"],
                       suffixes=("_11b", "_12a"))
    if len(merged) != len(old):
        raise RuntimeError(f"replicate-0 join produced {len(merged)} rows, expected {len(old)}")
    checks = {}
    for m in METRICS:
        d = float(np.abs(merged[f"{m}_11b"] - merged[f"{m}_12a"]).max())
        checks[m] = d
    save_json({"max_abs_difference_vs_phase11b": checks,
               "n_rows_compared": int(len(merged)), "timestamp_utc": utc_now()},
              paths["logs"] / "replicate0_reproduction_check.json")
    worst = max(checks.values())
    if worst > 1e-9:
        raise RuntimeError(
            f"replicate 0 did not reproduce PHASE 11B (max abs diff {worst:.3e}); "
            "the scoring path changed and the two replicates are not comparable")
    logger.info("replicate 0 reproduces PHASE 11B exactly (max abs diff %.3e)", worst)
    logger.info("score done: %s rows", len(frame))


# ---------------------------------------------------------------------------
# stage: variance
# ---------------------------------------------------------------------------
def _cube(frame: pd.DataFrame, model: str, metric: str) -> np.ndarray:
    """(subset, seed, replicate) array for one model and metric."""
    out = np.full((A, B, N), np.nan)
    sub = frame[frame.model == model]
    for i, s in enumerate(CROSSED_SUBSET_SEEDS):
        for j, k in enumerate(CROSSED_MODEL_SEEDS):
            for r, rep in enumerate(REPLICATES):
                sel = sub[(sub.source_subset_seed == s) & (sub.model_seed == k)
                          & (sub.replicate == rep)]
                if len(sel) != 1:
                    raise RuntimeError(
                        f"{model}/{metric}: expected 1 row for subset {s} seed {k} "
                        f"replicate {rep}, found {len(sel)}")
                out[i, j, r] = float(sel[metric].iloc[0])
    return out


def _phase11b_comparison(metric: str, model: str) -> dict:
    """The unreplicated PHASE 11B numbers, for the side-by-side reinterpretation."""
    old = pd.read_csv(phase11b_paths()["tables"] / "phase11b_variance_components.csv")
    key = "delta_totalvi_minus_scvi" if model == "delta" else model
    row = old[(old.model == key) & (old.metric == metric)]
    if len(row) != 1:
        return {}
    r = row.iloc[0]
    return {
        "phase11b_source_fraction": float(r["source_fraction"]),
        "phase11b_seed_fraction": float(r["seed_fraction"]),
        "phase11b_interaction_residual_fraction": float(r["interaction_fraction"]),
    }


def stage_variance(paths, logger) -> None:
    frame = pd.read_csv(paths["tables"] / "replicated_crossed_design.csv")
    cubes, rows, drows = {}, [], []
    for model in MODELS:
        for metric in METRICS:
            cube = _cube(frame, model, metric)
            cubes[f"{model}__{metric}"] = cube
            comp = replicated_components(cube)
            rows.append({"model": model, "metric": metric, **comp,
                         "matrix_min": float(cube.min()), "matrix_max": float(cube.max()),
                         "sd_over_all_100_cells": float(cube.std(ddof=1)),
                         **_phase11b_comparison(metric, model)})
            logger.info("%s %s: source %.3f seed %.3f interaction %.3f residual %.3f "
                        "(p_interaction=%.3g)", model, metric,
                        comp["source_subset_fraction"], comp["model_seed_fraction"],
                        comp["interaction_fraction"], comp["residual_fraction"],
                        comp["p_interaction"])

    # Delta = totalVI - scVI at matched (subset, seed, replicate).
    for metric in METRICS:
        d = cubes[f"totalvi__{metric}"] - cubes[f"scvi_matched__{metric}"]
        cubes[f"delta__{metric}"] = d
        comp = replicated_components(d)
        drows.append({"model": "delta_totalvi_minus_scvi", "metric": metric, **comp,
                      "matrix_min": float(d.min()), "matrix_max": float(d.max()),
                      "n_positive_of_50": int(np.sum(d > 0)),
                      "n_negative_of_50": int(np.sum(d < 0)),
                      "sd_over_all_100_cells": float(d.std(ddof=1)),
                      **_phase11b_comparison(metric, "delta")})
        logger.info("delta %s: source %.3f seed %.3f interaction %.3f residual %.3f "
                    "(p_interaction=%.3g)", metric, comp["source_subset_fraction"],
                    comp["model_seed_fraction"], comp["interaction_fraction"],
                    comp["residual_fraction"], comp["p_interaction"])

    pd.DataFrame(rows).to_csv(paths["tables"] / "replicated_variance_components.csv", index=False)
    pd.DataFrame(drows).to_csv(
        paths["tables"] / "replicated_delta_variance_components.csv", index=False)
    np.savez_compressed(paths["tables"] / "replicated_cubes.npz", **cubes)
    save_json({"variance_model": VARIANCE_MODEL_NOTE,
               "truncation_bias": TRUNCATION_BIAS_NOTE,
               "design": {"source_subsets": A, "model_seeds": B, "replicates": N,
                          "observations_per_model": A * B * N},
               "timestamp_utc": utc_now()},
              paths["logs"] / "variance_model.json")
    logger.info("variance done")


# ---------------------------------------------------------------------------
# stage: percell
# ---------------------------------------------------------------------------
def _percell_cube(paths, model: str) -> tuple[np.ndarray, np.ndarray]:
    """(n_cells, subset, seed, replicate) true-class probability, plus cell ids."""
    cells = None
    cube = None
    for i, s in enumerate(CROSSED_SUBSET_SEEDS):
        for j, k in enumerate(CROSSED_MODEL_SEEDS):
            for r, rep in enumerate(REPLICATES):
                tag = experiment_tag(CROSSED_DIRECTION_KEY, model, k, s, _rep_arg(rep))
                with np.load(paths["probabilities"] / f"rep_{tag}_probs.npz",
                             allow_pickle=True) as npz:
                    tcp = npz["true_class_probability"].astype(float)
                    ids = np.asarray(npz["cell_ids"], dtype=object)
                if cells is None:
                    cells = ids
                    cube = np.empty((tcp.size, A, B, N))
                elif not np.array_equal(cells, ids):
                    raise RuntimeError(f"{tag}: target cell order differs across runs")
                cube[:, i, j, r] = tcp
    return cube, cells


def stage_percell(paths, logger) -> None:
    setup = crossed_setup(paths, logger)
    y_target = setup["y_target"]
    frames = []
    for model in MODELS:
        cube, cells = _percell_cube(paths, model)
        if not np.array_equal(cells, setup["target_cells"]):
            raise RuntimeError(f"{model}: probability cell ids differ from the setup order")
        comp = replicated_components_vectorised(cube)
        df = pd.DataFrame({
            "cell_id": cells, "cell_type": y_target, "model": model,
            "mean_true_class_probability": cube.mean(axis=(1, 2, 3)),
            "total_variance": comp["total_variance"],
        })
        for name in ("source_subset", "model_seed", "interaction", "residual"):
            df[f"{name}_variance"] = comp[f"{name}_variance"]
            df[f"{name}_fraction"] = comp[f"{name}_fraction"]
        frames.append(df)
        logger.info("%s per-cell medians: source %.3f seed %.3f interaction %.3f residual %.3f",
                    model, df.source_subset_fraction.median(), df.model_seed_fraction.median(),
                    df.interaction_fraction.median(), df.residual_fraction.median())
    percell = pd.concat(frames, ignore_index=True)
    percell.to_csv(paths["tables"] / "per_cell_replicated_variance.csv", index=False)

    agg = (percell.groupby(["model", "cell_type"])
           .agg(n_cells=("cell_id", "size"),
                mean_true_class_probability=("mean_true_class_probability", "mean"),
                median_total_variance=("total_variance", "median"),
                median_source_fraction=("source_subset_fraction", "median"),
                median_seed_fraction=("model_seed_fraction", "median"),
                median_interaction_fraction=("interaction_fraction", "median"),
                median_residual_fraction=("residual_fraction", "median"),
                median_source_variance=("source_subset_variance", "median"),
                median_seed_variance=("model_seed_variance", "median"),
                median_interaction_variance=("interaction_variance", "median"),
                median_residual_variance=("residual_variance", "median"))
           .reset_index())
    agg.to_csv(paths["tables"] / "celltype_replicated_variance.csv", index=False)
    logger.info("percell done: %s cells x %s models", len(percell) // len(MODELS), len(MODELS))


STAGES = {"score": stage_score, "variance": stage_variance, "percell": stage_percell}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=[*STAGES, "all"])
    args = ap.parse_args()
    paths = ensure_dirs()
    logger = get_logger("phase12a_analysis", paths["logs"] / "phase12a_analysis.log")
    todo = list(STAGES) if args.stage == "all" else [args.stage]
    for name in todo:
        logger.info("=== stage %s (%s) ===", name, utc_now())
        STAGES[name](paths, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
