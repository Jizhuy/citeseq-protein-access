#!/usr/bin/env python
"""PHASE 11B Parts III-IV: crossed source-subset x model-seed variance decomposition.

Direction A only (PBMC10k -> PBMC5k), 5 PHASE 10B source subsets x 5 model seeds
x 2 models. The eligible l2 class set is fixed to the intersection across all
five subsets so every one of the 25 grid cells scores the same target cells.

Stages
  score     score all 50 runs on the fixed class set
  variance  5x5 metric matrices and two-factor variance components
  percell   per-cell and per-cell-type variance decomposition
  assoc     full-source uncertainty versus each instability source
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase10_cross_dataset import class_plan, load_source_target  # noqa: E402
from src.experiments.phase11_uncertainty import (  # noqa: E402
    calibration_metrics,
    classification_metrics,
    error_detection,
    fit_transfer_classifier,
    risk_coverage,
)
from src.experiments.phase11b_analysis import load_probabilities, load_run  # noqa: E402
from src.experiments.phase11b_stats import (  # noqa: E402
    two_factor_components,
    variation_ratio,
    vector_components,
)
from src.experiments.phase11b_robustness import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    MODELS,
    SEEDS,
    cleanup,
    ensure_dirs,
    experiment_tag,
    utc_now,
)
from src.utils.io import save_json  # noqa: E402
from src.utils.logging_utils import get_logger  # noqa: E402

DIRECTION_KEY = "direction_A"
CROSSED_METRICS = ("macro_f1", "error_auroc", "error_auprc", "nll", "brier", "ece", "aurc")
QUANTILE = 0.20


def crossed_setup(paths, logger):
    """Fixed target evaluation set: high-confidence cells in the common class set."""
    plan_path = paths["logs"] / "crossed_class_plan.json"
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
        "source_cells_full": source_full.obs_names.astype(str).to_numpy(),
        "n_source_full": int(source_full.n_obs),
    }
    logger.info("crossed setup: %s common l2 classes, %s fixed target cells",
                len(common), setup["target_index"].size)
    del source_full, target
    cleanup()
    return setup


def stage_score(paths, logger) -> None:
    from src.experiments.phase10b_source_size import SUBSAMPLE_N, subsample_index
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
                tag = experiment_tag(DIRECTION_KEY, model, model_seed, subset_seed)
                if not (paths["done"] / f"{tag}.json").exists():
                    logger.warning("missing %s", tag)
                    continue
                run = load_run(DIRECTION_KEY, model, model_seed, subset_seed)
                z_src = run["source_latent"][src_keep]
                z_tgt = run["target_latent"][setup["target_index"]]
                u_lat = run["u_latent_target"][setup["target_index"]]
                y = setup["y_target"]
                fit = fit_transfer_classifier(z_src, y_source, z_tgt, model_seed)
                probs, classes, pred = fit["probs"], fit["classes"], fit["pred"]
                if list(classes) != sorted(common):
                    raise RuntimeError(f"{tag}: classifier classes {list(classes)} != common set")
                correct = pred == y
                err = (~correct).astype(int)
                u = 1.0 - probs.max(axis=1)
                det = error_detection(u, err)
                cal = calibration_metrics(probs, y, classes, pred)
                rc = risk_coverage(u, correct.astype(float))
                np.savez_compressed(
                    paths["probabilities"] / f"crossed_{tag}_probs.npz",
                    probs=probs.astype(np.float32),
                    classes=np.asarray(classes, dtype=object),
                    cell_ids=np.asarray(setup["target_cells"], dtype=object),
                    u_pred=u.astype(np.float32), u_latent=u_lat.astype(np.float32),
                    true_class_probability=cal["true_class_probability"].astype(np.float32),
                    pred=np.asarray(pred, dtype=object))
                rows.append({"model": model, "source_subset_seed": subset_seed,
                             "model_seed": model_seed,
                             **classification_metrics(y, pred),
                             "error_auroc": det["auroc"], "error_auprc": det["auprc"],
                             "error_prevalence": det["error_prevalence"],
                             "nll": cal["nll"], "brier": cal["brier"],
                             "ece": cal["ece_equal_frequency"], "aurc": rc["aurc"],
                             "n_source_eval": int(y_source.size), "n_target_eval": int(y.size)})
                logger.info("%s macroF1=%.4f AUROC=%.4f", tag, rows[-1]["macro_f1"],
                            rows[-1]["error_auroc"])
                del run, fit, probs
                cleanup()
    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_crossed_design.csv", index=False)
    logger.info("score done: %s rows", len(rows))


def _matrix(frame: pd.DataFrame, model: str, metric: str) -> np.ndarray:
    piv = frame[frame["model"] == model].pivot(index="source_subset_seed",
                                               columns="model_seed", values=metric)
    return piv.reindex(index=CROSSED_SUBSET_SEEDS, columns=CROSSED_MODEL_SEEDS).to_numpy()


def stage_variance(paths, logger) -> None:
    frame = pd.read_csv(paths["tables"] / "phase11b_crossed_design.csv")
    rows, matrices = [], {}
    for model in MODELS:
        for metric in CROSSED_METRICS:
            mat = _matrix(frame, model, metric)
            if np.isnan(mat).any():
                logger.warning("incomplete matrix for %s %s", model, metric)
                continue
            matrices[f"{model}__{metric}"] = mat
            comp = two_factor_components(mat)
            rows.append({"model": model, "metric": metric, **comp,
                         "matrix_min": float(mat.min()), "matrix_max": float(mat.max()),
                         "matrix_range": float(mat.max() - mat.min()),
                         "sd_over_all_25": float(mat.std(ddof=1))})
    # Delta = totalVI - scVI on the same subset x seed cell (§32)
    for metric in CROSSED_METRICS:
        a, b = f"totalvi__{metric}", f"scvi_matched__{metric}"
        if a not in matrices or b not in matrices:
            continue
        d = matrices[a] - matrices[b]
        matrices[f"delta__{metric}"] = d
        comp = two_factor_components(d)
        rows.append({"model": "delta_totalvi_minus_scvi", "metric": metric, **comp,
                     "matrix_min": float(d.min()), "matrix_max": float(d.max()),
                     "matrix_range": float(d.max() - d.min()),
                     "sd_over_all_25": float(d.std(ddof=1)),
                     "n_positive_of_25": int(np.sum(d > 0)),
                     "n_negative_of_25": int(np.sum(d < 0))})
    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_variance_components.csv", index=False)
    np.savez_compressed(paths["tables"] / "phase11b_crossed_matrices.npz", **matrices)
    logger.info("variance done: %s rows", len(rows))


def _load_crossed_cube(paths, model: str, setup) -> dict:
    """(5 subsets, 5 seeds, n_cells) true-class probability and predicted labels."""
    tcp, preds, upred = [], [], []
    for subset_seed in CROSSED_SUBSET_SEEDS:
        r_t, r_p, r_u = [], [], []
        for model_seed in CROSSED_MODEL_SEEDS:
            tag = experiment_tag(DIRECTION_KEY, model, model_seed, subset_seed)
            with np.load(paths["probabilities"] / f"crossed_{tag}_probs.npz",
                         allow_pickle=True) as npz:
                r_t.append(npz["true_class_probability"].astype(float))
                r_p.append(np.asarray(npz["pred"], dtype=object))
                r_u.append(npz["u_pred"].astype(float))
        tcp.append(r_t); preds.append(r_p); upred.append(r_u)
    return {"tcp": np.array(tcp), "pred": np.array(preds, dtype=object),
            "u_pred": np.array(upred)}


def stage_percell(paths, logger) -> None:
    setup = crossed_setup(paths, logger)
    y = setup["y_target"]
    cells = setup["target_cells"]
    rows, ct_rows = [], []
    for model in MODELS:
        cube = _load_crossed_cube(paths, model, setup)
        tcp, pred = cube["tcp"], cube["pred"]
        comp = vector_components(tcp)
        # A: across model seeds with the source subset held fixed
        seed_vr = np.mean([variation_ratio(np.stack(list(pred[i]))) for i in range(5)], axis=0)
        seed_tcp_var = tcp.var(axis=1, ddof=1).mean(axis=0)
        # B: across source subsets with the model seed held fixed
        subset_vr = np.mean([variation_ratio(np.stack([pred[i][j] for i in range(5)]))
                             for j in range(5)], axis=0)
        subset_tcp_var = tcp.var(axis=0, ddof=1).mean(axis=0)

        full = _full_source_uncertainty(paths, model, cells, y, logger)
        rows.extend(pd.DataFrame({
            "cell_id": cells, "cell_type": y, "model": model,
            "correct_full_source": full["correct"],
            "predictive_uncertainty": full["u_pred"],
            "latent_uncertainty": full["u_latent"],
            "ensemble_disagreement": full["disagreement"],
            "source_subset_variation_ratio": subset_vr,
            "model_seed_variation_ratio": seed_vr,
            "source_subset_trueclass_var": subset_tcp_var,
            "model_seed_trueclass_var": seed_tcp_var,
            "interaction_trueclass_var": comp["interaction_var"],
            "source_variance_fraction": comp["source_fraction"],
            "seed_variance_fraction": comp["seed_fraction"],
            "interaction_variance_fraction": comp["interaction_fraction"],
        }).to_dict("records"))

        frame = pd.DataFrame({"cell_type": y, "err": (~full["correct"]).astype(float),
                              "u_pred": full["u_pred"], "u_lat": full["u_latent"],
                              "sf": comp["source_fraction"], "kf": comp["seed_fraction"],
                              "if_": comp["interaction_fraction"],
                              "sv": comp["source_var"], "kv": comp["seed_var"]})
        src_counts = pd.Series(setup["source_labels_full"][setup["source_high_full"]]).value_counts()
        for ct, g in frame.groupby("cell_type"):
            ct_rows.append({"model": model, "cell_type": ct,
                            "source_support": int(src_counts.get(ct, 0)),
                            "target_support": int(len(g)),
                            "error_rate": float(g["err"].mean()),
                            "median_predictive_uncertainty": float(g["u_pred"].median()),
                            "median_latent_uncertainty": float(g["u_lat"].median()),
                            "source_variance_fraction": float(g["sf"].median()),
                            "seed_variance_fraction": float(g["kf"].median()),
                            "interaction_fraction": float(g["if_"].median()),
                            "median_source_variance": float(g["sv"].median()),
                            "median_seed_variance": float(g["kv"].median())})
        logger.info("%s per-cell decomposition done (median source frac %.3f, seed frac %.3f)",
                    model, np.nanmedian(comp["source_fraction"]),
                    np.nanmedian(comp["seed_fraction"]))
        del cube, comp
        cleanup()

    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_cell_instability.csv", index=False)
    ct = pd.DataFrame(ct_rows)
    ct.to_csv(paths["tables"] / "phase11b_celltype_variance.csv", index=False)

    # §36 source support versus variance components / error / uncertainty
    supp = []
    for model in MODELS:
        g = ct[ct["model"] == model]
        for col in ("median_seed_variance", "median_source_variance", "error_rate",
                    "median_predictive_uncertainty", "seed_variance_fraction",
                    "source_variance_fraction"):
            r = spearmanr(g["source_support"], g[col])
            supp.append({"model": model, "against": col, "spearman_rho": float(r.statistic),
                         "p_value": float(r.pvalue), "n_cell_types": int(len(g))})
    pd.DataFrame(supp).to_csv(paths["tables"] / "phase11b_source_support_association.csv", index=False)
    logger.info("percell done")


def _full_source_uncertainty(paths, model: str, crossed_cells: np.ndarray,
                             y_target: np.ndarray, logger) -> dict:
    """Direction A full-source 20-seed uncertainty, aligned to the crossed cells.

    The full-source eval set uses the Direction A eligible-class plan while the
    crossed grid uses the narrower cross-subset common set, so cells are matched
    by ID rather than assumed to be in the same order.
    """
    seeds = [s for s in SEEDS
             if (paths["done"] / f"{experiment_tag(DIRECTION_KEY, model, s, None)}.json").exists()]
    packs = [load_probabilities(DIRECTION_KEY, experiment_tag(DIRECTION_KEY, model, s, None))
             for s in seeds]
    ids = np.asarray(packs[0]["cell_ids"], dtype=str)
    pos = pd.Series(np.arange(ids.size), index=ids)
    missing = [c for c in crossed_cells if c not in pos.index]
    if missing:
        raise RuntimeError(f"{model}: {len(missing)} crossed cells absent from full-source eval")
    sel = pos.loc[crossed_cells].to_numpy()

    ref = list(packs[0]["classes"])
    stack = []
    for p in packs:
        cls = list(p["classes"])
        if cls != ref:
            if set(cls) != set(ref):
                raise RuntimeError("class set mismatch across full-source seeds")
            p["probs"] = p["probs"][:, [cls.index(c) for c in ref]]
        stack.append(p["probs"][sel])
    stack = np.stack(stack)
    p_bar = stack.mean(axis=0)
    classes = np.asarray(ref, dtype=object)

    def ent(p):
        c = np.clip(p, 1e-12, 1.0)
        return -np.sum(c * np.log(c), axis=1)

    disagreement = ent(p_bar) - np.mean([ent(stack[i]) for i in range(stack.shape[0])], axis=0)
    u_pred = 1.0 - p_bar.max(axis=1)
    u_latent = np.mean([p["u_latent"][sel] for p in packs], axis=0)
    pred = classes[np.argmax(p_bar, axis=1)]
    correct = np.asarray(pred == y_target, dtype=bool)
    logger.info("%s full-source ensemble over %s seeds aligned to %s crossed cells "
                "(ensemble accuracy %.4f)", model, len(seeds), crossed_cells.size,
                float(correct.mean()))
    return {"u_pred": u_pred, "u_latent": u_latent, "disagreement": disagreement,
            "pred": pred, "correct": correct, "n_seeds": len(seeds)}


def stage_assoc(paths, logger) -> None:
    """Part IV: does full-source uncertainty predict either instability source?"""
    cell = pd.read_csv(paths["tables"] / "phase11b_cell_instability.csv")
    rows, quant_rows = [], []
    unc_cols = ("predictive_uncertainty", "latent_uncertainty", "ensemble_disagreement")
    inst_cols = ("source_subset_variation_ratio", "model_seed_variation_ratio",
                 "source_subset_trueclass_var", "model_seed_trueclass_var")
    for model in MODELS:
        g = cell[cell["model"] == model]
        total_vr = g["source_subset_variation_ratio"] + g["model_seed_variation_ratio"]
        for u in unc_cols:
            for inst in inst_cols:
                r = spearmanr(g[u], g[inst])
                rows.append({"model": model, "uncertainty": u, "instability": inst,
                             "spearman_rho": float(r.statistic), "p_value": float(r.pvalue),
                             "n_cells": int(len(g))})
            r = spearmanr(g[u], total_vr)
            rows.append({"model": model, "uncertainty": u, "instability": "total_crossed_variation",
                         "spearman_rho": float(r.statistic), "p_value": float(r.pvalue),
                         "n_cells": int(len(g))})
            # §42 pre-specified top/bottom 20% strata
            lo_q, hi_q = g[u].quantile(QUANTILE), g[u].quantile(1 - QUANTILE)
            low, high = g[g[u] <= lo_q], g[g[u] >= hi_q]
            for inst in inst_cols:
                stat = mannwhitneyu(high[inst], low[inst], alternative="two-sided")
                quant_rows.append({"model": model, "uncertainty": u, "instability": inst,
                                   "mean_top20": float(high[inst].mean()),
                                   "mean_bottom20": float(low[inst].mean()),
                                   "ratio_top_over_bottom": float(high[inst].mean() / low[inst].mean())
                                   if low[inst].mean() > 0 else float("inf"),
                                   "mannwhitney_p": float(stat.pvalue),
                                   "n_top": int(len(high)), "n_bottom": int(len(low))})
        # §41 correct versus incorrect instability
        if g["correct_full_source"].notna().any():
            for inst in inst_cols:
                ok = g[g["correct_full_source"] == True][inst]   # noqa: E712
                bad = g[g["correct_full_source"] == False][inst]  # noqa: E712
                if len(ok) and len(bad):
                    stat = mannwhitneyu(bad, ok, alternative="two-sided")
                    quant_rows.append({"model": model, "uncertainty": "correct_vs_incorrect",
                                       "instability": inst,
                                       "mean_top20": float(bad.mean()),
                                       "mean_bottom20": float(ok.mean()),
                                       "ratio_top_over_bottom": float(bad.mean() / ok.mean())
                                       if ok.mean() > 0 else float("inf"),
                                       "mannwhitney_p": float(stat.pvalue),
                                       "n_top": int(len(bad)), "n_bottom": int(len(ok))})
    pd.DataFrame(rows).to_csv(paths["tables"] / "phase11b_uncertainty_instability_assoc.csv", index=False)
    pd.DataFrame(quant_rows).to_csv(paths["tables"] / "phase11b_instability_quantiles.csv", index=False)
    logger.info("assoc done")


STAGES = {"score": stage_score, "variance": stage_variance,
          "percell": stage_percell, "assoc": stage_assoc}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=[*STAGES, "all"])
    args = ap.parse_args()
    paths = ensure_dirs()
    logger = get_logger("phase11b_crossed", paths["logs"] / "phase11b_crossed.log")
    for name in (list(STAGES) if args.stage == "all" else [args.stage]):
        logger.info("=== crossed stage %s (%s) ===", name, utc_now())
        STAGES[name](paths, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
