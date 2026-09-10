#!/usr/bin/env python
"""PHASE 11 main runner: uncertainty characterization and calibration.

Stages
  main      per-seed predictive uncertainty, calibration, ensemble,
            neighbor stability, cell-type tables, bootstrap comparison
  latent    source-cell latent posterior variance (target side unsupported)
  subset    PHASE 10B source-subset instability and its association with
            full-source PHASE 10 uncertainty

Trains no representation model. Reuses PHASE 10 / 10B artifacts only.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import (  # noqa: E402
    DIRECTIONS,
    LATENT_UNSUPPORTED_REASON,
    MODELS,
    NEIGHBOR_K,
    N_BOOT,
    REPRO_TOLERANCE,
    SEEDS,
    DirectionData,
    align_probability_matrices,
    bootstrap_indices,
    calibration_metrics,
    classification_metrics,
    ci_bounds,
    coverage_table,
    ensemble_uncertainty,
    ensure_dirs,
    entropy_of,
    error_detection,
    error_enrichment,
    fit_transfer_classifier,
    load_phase10_pack,
    neighbor_label_agreement,
    neighbor_stability,
    phase10_reference_metrics,
    predictive_uncertainty,
    risk_coverage,
    source_neighbor_sets,
    utc_now,
    verify_cell_alignment,
)
from src.utils.io import save_json  # noqa: E402

LEVEL = "l2"
CAL_HOLDOUT = 0.25
CAL_SPLIT_SEED = 1111


def rss_mb() -> float:
    import psutil

    return float(psutil.Process().memory_info().rss) / (1024.0 * 1024.0)


# ----------------------------------------------------------------------------
# Post-hoc temperature scaling, fitted on SOURCE labels only (§17)
# ----------------------------------------------------------------------------
def fit_temperature(z_src: np.ndarray, y_src: np.ndarray, seed: int) -> dict:
    """Hold out a stratified source split, fit T on it, never touch target labels.

    The classifier used for the primary (uncalibrated) result is fitted on the
    full source, exactly as PHASE 10 did. To estimate T out-of-sample we fit a
    companion classifier on the complementary split and optimise T on the
    held-out source cells, then transfer that T. This keeps the primary
    reproduction untouched and keeps every label used for calibration on the
    source side.
    """
    counts = pd.Series(y_src).value_counts()
    if int(counts.min()) < 2:
        return {"temperature": float("nan"), "fitted": False, "reason": "class with <2 source cells"}
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=CAL_HOLDOUT, random_state=CAL_SPLIT_SEED)
    fit_idx, hold_idx = next(splitter.split(z_src, y_src))
    companion = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    companion.fit(z_src[fit_idx], y_src[fit_idx])
    logits = companion.decision_function(z_src[hold_idx])
    if logits.ndim == 1:
        logits = np.column_stack([-logits, logits])
    classes = list(companion.classes_)
    y_hold = y_src[hold_idx]
    target_idx = np.array([classes.index(c) if c in classes else -1 for c in y_hold])
    keep = target_idx >= 0
    logits, target_idx = logits[keep], target_idx[keep]

    def nll_at(temp: float) -> float:
        scaled = logits / temp
        scaled = scaled - scaled.max(axis=1, keepdims=True)
        logsum = np.log(np.exp(scaled).sum(axis=1))
        return float(-np.mean(scaled[np.arange(target_idx.size), target_idx] - logsum))

    grid = np.exp(np.linspace(np.log(0.25), np.log(8.0), 160))
    losses = np.array([nll_at(t) for t in grid])
    best = float(grid[int(np.argmin(losses))])
    return {
        "temperature": best,
        "fitted": True,
        "holdout_n": int(keep.sum()),
        "holdout_nll_at_T1": nll_at(1.0),
        "holdout_nll_at_T": float(losses.min()),
    }


def apply_temperature(clf, z_tgt: np.ndarray, temperature: float) -> np.ndarray:
    logits = clf.decision_function(z_tgt)
    if logits.ndim == 1:
        logits = np.column_stack([-logits, logits])
    scaled = logits / temperature
    scaled = scaled - scaled.max(axis=1, keepdims=True)
    exp = np.exp(scaled)
    return exp / exp.sum(axis=1, keepdims=True)


# ----------------------------------------------------------------------------
# Stage: main
# ----------------------------------------------------------------------------
def run_main() -> None:
    paths = ensure_dirs()
    reference = phase10_reference_metrics()
    primary_rows, ensemble_rows, celltype_rows, neighbor_rows = [], [], [], []
    reliability_rows, coverage_rows, roc_rows, delta_rows = [], [], [], []
    dist_rows, support_rows = [], []
    repro_log = []

    for direction, source_batch, target_batch in DIRECTIONS:
        print(f"\n########## {direction} ##########", flush=True)
        data = DirectionData(direction, source_batch, target_batch)
        sl = data.eval_slice(LEVEL)
        y_tgt, y_src = sl["y_target"], sl["y_source"]
        n_tgt = y_tgt.size
        src_counts = pd.Series(y_src).value_counts()
        tgt_counts = pd.Series(y_tgt).value_counts()
        boot_idx = bootstrap_indices(n_tgt)
        per_model_seed: dict[str, dict[int, dict]] = {m: {} for m in MODELS}
        per_model_neighbors: dict[str, list[np.ndarray]] = {m: [] for m in MODELS}
        per_model_agreement: dict[str, list[np.ndarray]] = {m: [] for m in MODELS}

        for model in MODELS:
            for seed in SEEDS:
                pack = load_phase10_pack(direction, model, seed)
                verify_cell_alignment(data, pack)
                z_src = pack["source_latent"][sl["source_index"]]
                z_tgt = pack["target_latent"][sl["target_index"]]

                fit = fit_transfer_classifier(z_src, y_src, z_tgt, seed)
                metrics = classification_metrics(y_tgt, fit["pred"])

                ref = reference[
                    (reference["direction"] == direction)
                    & (reference["model"] == model)
                    & (reference["model_seed"] == seed)
                    & (reference["label_level"] == LEVEL)
                ]
                ref_macro = float(ref["macro_f1"].iloc[0])
                d_macro = abs(metrics["macro_f1"] - ref_macro)
                if d_macro > REPRO_TOLERANCE:
                    raise RuntimeError(
                        f"PHASE 10 reproduction failed for {direction}/{model}/seed{seed}: "
                        f"phase10={ref_macro} phase11={metrics['macro_f1']}"
                    )
                repro_log.append(
                    {
                        "direction": direction,
                        "model": model,
                        "model_seed": seed,
                        "label_level": LEVEL,
                        "phase10_macro_f1": ref_macro,
                        "phase11_macro_f1": metrics["macro_f1"],
                        "abs_delta": d_macro,
                    }
                )

                error = (fit["pred"] != y_tgt).astype(int)
                correct = 1 - error
                unc = predictive_uncertainty(fit["probs"])
                det = error_detection(unc["u_pred"], error)
                cal = calibration_metrics(fit["probs"], y_tgt, fit["classes"], fit["pred"])
                rc = risk_coverage(unc["u_pred"], correct)

                # Post-hoc temperature scaling, source-fitted only.
                temp = fit_temperature(z_src, y_src, seed)
                if temp["fitted"]:
                    probs_cal = apply_temperature(fit["classifier"], z_tgt, temp["temperature"])
                    pred_cal = fit["classes"][np.argmax(probs_cal, axis=1)]
                    cal_cal = calibration_metrics(probs_cal, y_tgt, fit["classes"], pred_cal)
                    unc_cal = predictive_uncertainty(probs_cal)
                    det_cal = error_detection(unc_cal["u_pred"], (pred_cal != y_tgt).astype(int))
                else:
                    cal_cal = {"nll": np.nan, "brier": np.nan, "ece_equal_frequency": np.nan}
                    det_cal = {"auroc": np.nan}

                primary_rows.append(
                    {
                        "direction": direction,
                        "model": model,
                        "model_seed": seed,
                        "label_level": LEVEL,
                        "n_cells": n_tgt,
                        "accuracy": metrics["accuracy"],
                        "balanced_accuracy": metrics["balanced_accuracy"],
                        "macro_f1": metrics["macro_f1"],
                        "weighted_f1": metrics["weighted_f1"],
                        "error_prevalence": det["error_prevalence"],
                        "error_detection_auroc_pred": det["auroc"],
                        "error_detection_auprc_pred": det["auprc"],
                        "error_detection_auroc_latent": np.nan,
                        "error_detection_auprc_latent": np.nan,
                        "nll": cal["nll"],
                        "brier": cal["brier"],
                        "ece_equal_frequency": cal["ece_equal_frequency"],
                        "ece_equal_width": cal["ece_equal_width"],
                        "aurc": rc["aurc"],
                        "median_predictive_uncertainty": float(np.median(unc["u_pred"])),
                        "median_latent_uncertainty": np.nan,
                        "latent_uncertainty_status": "unsupported_see_posterior_api_validation_md",
                        "error_detection_auroc_entropy": error_detection(unc["entropy"], error)["auroc"],
                        "error_detection_auroc_neg_margin": error_detection(-unc["margin"], error)["auroc"],
                        "temperature": temp.get("temperature", np.nan),
                        "nll_temperature_scaled": cal_cal["nll"],
                        "brier_temperature_scaled": cal_cal["brier"],
                        "ece_eq_freq_temperature_scaled": cal_cal["ece_equal_frequency"],
                        "error_detection_auroc_temperature_scaled": det_cal["auroc"],
                    }
                )

                for scheme, key in (
                    ("equal_frequency", "reliability_equal_frequency"),
                    ("equal_width", "reliability_equal_width"),
                ):
                    for b in cal[key]:
                        reliability_rows.append(
                            {"direction": direction, "model": model, "model_seed": seed,
                             "scheme": scheme, **b}
                        )
                for row in coverage_table(unc["u_pred"], y_tgt, fit["pred"]):
                    coverage_rows.append(
                        {"direction": direction, "model": model, "model_seed": seed, **row}
                    )
                enrich = error_enrichment(unc["u_pred"], error)
                primary_rows[-1].update(enrich)

                if seed == 0:
                    step = max(1, n_tgt // 400)
                    roc_rows.append(
                        {
                            "direction": direction,
                            "model": model,
                            "coverage": rc["coverage"][::step].tolist(),
                            "risk": rc["risk"][::step].tolist(),
                        }
                    )
                    dist_rows.extend(
                        [
                            {"direction": direction, "model": model, "group": "correct",
                             "u_pred": float(v)} for v in unc["u_pred"][correct == 1][::3]
                        ]
                        + [
                            {"direction": direction, "model": model, "group": "incorrect",
                             "u_pred": float(v)} for v in unc["u_pred"][error == 1][::3]
                        ]
                    )

                nbr = source_neighbor_sets(z_src, z_tgt, NEIGHBOR_K)
                per_model_neighbors[model].append(nbr)
                per_model_agreement[model].append(neighbor_label_agreement(nbr, y_src, y_tgt))

                per_model_seed[model][seed] = {
                    "probs": fit["probs"],
                    "classes": fit["classes"],
                    "pred": fit["pred"],
                    "u_pred": unc["u_pred"],
                    "error": error,
                    "true_class_probability": cal["true_class_probability"],
                }
                np.savez_compressed(
                    paths["arrays"] / f"probs_{direction}_{model}_seed{seed}_{LEVEL}.npz",
                    probs=fit["probs"].astype(np.float32),
                    classes=np.asarray([str(c) for c in fit["classes"]]),
                    cell_ids=sl["target_cells"].astype(str),
                    y_true=y_tgt.astype(str),
                    pred=fit["pred"].astype(str),
                    direction=np.asarray([direction]),
                    model=np.asarray([model]),
                    model_seed=np.asarray([seed]),
                )
                print(
                    f"  {model} seed{seed}: macroF1={metrics['macro_f1']:.4f} "
                    f"errAUROC={det['auroc']:.4f} ECE={cal['ece_equal_frequency']:.4f} "
                    f"AURC={rc['aurc']:.4f} T={temp.get('temperature', float('nan')):.3f}",
                    flush=True,
                )
                del pack, fit
                gc.collect()

        # ---------------- ensemble across the five model seeds ----------------
        ens_store = {}
        for model in MODELS:
            stacked, classes = align_probability_matrices(per_model_seed[model])
            ens = ensemble_uncertainty(stacked, classes)
            ens_error = (ens["ensemble_pred"] != y_tgt).astype(int)
            ens_metrics = classification_metrics(y_tgt, ens["ensemble_pred"])
            ens_cal = calibration_metrics(ens["p_bar"], y_tgt, classes, ens["ensemble_pred"])
            ens_u = predictive_uncertainty(ens["p_bar"])
            ens_det = error_detection(ens_u["u_pred"], ens_error)
            ens_rc = risk_coverage(ens_u["u_pred"], 1 - ens_error)
            dis_det = error_detection(ens["disagreement"], ens_error)
            var_det = error_detection(ens["variation_ratio"], ens_error)
            ens_store[model] = {
                "ens": ens, "error": ens_error, "u_pred": ens_u["u_pred"], "classes": classes,
            }
            ensemble_rows.append(
                {
                    "direction": direction,
                    "model": model,
                    "label_level": LEVEL,
                    "n_cells": n_tgt,
                    "ensemble_accuracy": ens_metrics["accuracy"],
                    "ensemble_balanced_accuracy": ens_metrics["balanced_accuracy"],
                    "ensemble_macro_f1": ens_metrics["macro_f1"],
                    "ensemble_weighted_f1": ens_metrics["weighted_f1"],
                    "ensemble_nll": ens_cal["nll"],
                    "ensemble_brier": ens_cal["brier"],
                    "ensemble_ece": ens_cal["ece_equal_frequency"],
                    "ensemble_error_auroc": ens_det["auroc"],
                    "ensemble_error_auprc": ens_det["auprc"],
                    "ensemble_aurc": ens_rc["aurc"],
                    "median_predictive_entropy": float(np.median(ens["predictive_entropy"])),
                    "median_expected_entropy": float(np.median(ens["expected_entropy"])),
                    "median_ensemble_disagreement": float(np.median(ens["disagreement"])),
                    "median_variation_ratio": float(np.median(ens["variation_ratio"])),
                    "disagreement_error_auroc": dis_det["auroc"],
                    "disagreement_error_auprc": dis_det["auprc"],
                    "variation_ratio_error_auroc": var_det["auroc"],
                    "mean_seed_macro_f1": float(
                        np.mean([
                            f1_score(y_tgt, per_model_seed[model][s]["pred"],
                                     average="macro", zero_division=0) for s in SEEDS
                        ])
                    ),
                    "mean_seed_ece": float(
                        np.mean([
                            r["ece_equal_frequency"] for r in primary_rows
                            if r["direction"] == direction and r["model"] == model
                        ])
                    ),
                }
            )
            print(
                f"  [ensemble] {model}: macroF1={ens_metrics['macro_f1']:.4f} "
                f"errAUROC={ens_det['auroc']:.4f} disagreeAUROC={dis_det['auroc']:.4f}",
                flush=True,
            )

        # ---------------- neighbor stability across seeds ----------------
        for model in MODELS:
            stab = neighbor_stability(per_model_neighbors[model])
            agree = np.stack(per_model_agreement[model], axis=0)
            agree_mean, agree_sd = agree.mean(axis=0), agree.std(axis=0, ddof=0)
            ens = ens_store[model]
            correct0 = 1 - per_model_seed[model][0]["error"]
            for i in range(n_tgt):
                neighbor_rows.append(
                    {
                        "direction": direction,
                        "model": model,
                        "cell_id": sl["target_cells"][i],
                        "cell_type": y_tgt[i],
                        "correct": int(correct0[i]),
                        "neighbor_jaccard_mean": float(stab[i]),
                        "neighbor_label_agreement_mean": float(agree_mean[i]),
                        "neighbor_label_agreement_sd": float(agree_sd[i]),
                        "predictive_uncertainty": float(per_model_seed[model][0]["u_pred"][i]),
                        "ensemble_disagreement": float(ens["ens"]["disagreement"][i]),
                    }
                )
            err0 = per_model_seed[model][0]["error"]
            det_stab = error_detection(-stab, err0)
            support_rows.append(
                {
                    "direction": direction,
                    "model": model,
                    "metric": "neighbor_stability_vs_error",
                    "spearman_rho": float(spearmanr(stab, err0).statistic),
                    "spearman_p": float(spearmanr(stab, err0).pvalue),
                    "auroc_low_stability_detects_error": det_stab["auroc"],
                    "spearman_stability_vs_u_pred": float(
                        spearmanr(stab, per_model_seed[model][0]["u_pred"]).statistic
                    ),
                    "spearman_stability_vs_disagreement": float(
                        spearmanr(stab, ens["ens"]["disagreement"]).statistic
                    ),
                    "mean_stability_correct": float(stab[err0 == 0].mean()),
                    "mean_stability_incorrect": float(stab[err0 == 1].mean()),
                }
            )

            # ---------------- per-cell-type table ----------------
            ens_dis = ens["ens"]["disagreement"]
            u0 = per_model_seed[model][0]["u_pred"]
            for cls in sl["classes"]:
                mask = y_tgt == cls
                if mask.sum() == 0:
                    continue
                pred_cls = per_model_seed[model][0]["pred"]
                celltype_rows.append(
                    {
                        "direction": direction,
                        "model": model,
                        "cell_type": cls,
                        "source_support": int(src_counts.get(cls, 0)),
                        "target_support": int(mask.sum()),
                        "accuracy": float(accuracy_score(y_tgt[mask], pred_cls[mask])),
                        "f1": float(
                            f1_score(y_tgt, pred_cls, labels=[cls], average="macro", zero_division=0)
                        ),
                        "error_rate": float(np.mean(pred_cls[mask] != y_tgt[mask])),
                        "median_predictive_uncertainty": float(np.median(u0[mask])),
                        "median_latent_uncertainty": np.nan,
                        "median_ensemble_disagreement": float(np.median(ens_dis[mask])),
                        "neighbor_stability": float(np.mean(stab[mask])),
                        "neighbor_label_agreement": float(np.mean(agree_mean[mask])),
                    }
                )

        # ---------------- paired bootstrap: totalVI - scVI ----------------
        def metric_values(model: str, idx: np.ndarray, which: str) -> float:
            s = per_model_seed[model][0]
            if which == "error_auroc":
                return error_detection(s["u_pred"][idx], s["error"][idx])["auroc"]
            if which == "error_auprc":
                return error_detection(s["u_pred"][idx], s["error"][idx])["auprc"]
            if which == "aurc":
                return risk_coverage(s["u_pred"][idx], 1 - s["error"][idx])["aurc"]
            probs, classes = s["probs"][idx], s["classes"]
            cal = calibration_metrics(probs, y_tgt[idx], classes, s["pred"][idx])
            return cal["ece_equal_frequency"] if which == "ece" else cal["brier"]

        for metric in ("error_auroc", "error_auprc", "ece", "brier", "aurc"):
            sc = metric_values("scvi_matched", np.arange(n_tgt), metric)
            tv = metric_values("totalvi", np.arange(n_tgt), metric)
            draws = np.array(
                [
                    metric_values("totalvi", boot_idx[b], metric)
                    - metric_values("scvi_matched", boot_idx[b], metric)
                    for b in range(N_BOOT)
                ]
            )
            lo, hi = ci_bounds(draws)
            delta_rows.append(
                {
                    "direction": direction,
                    "label_level": LEVEL,
                    "metric": metric,
                    "scvi_value": sc,
                    "totalvi_value": tv,
                    "delta_totalvi_minus_scvi": tv - sc,
                    "bootstrap_ci_low": lo,
                    "bootstrap_ci_high": hi,
                    "higher_is_better": metric in ("error_auroc", "error_auprc"),
                    "n_boot": N_BOOT,
                }
            )
            print(f"  [boot] {metric}: delta={tv - sc:+.4f} CI=[{lo:+.4f},{hi:+.4f}]", flush=True)

        # ---------------- source support vs uncertainty (§30) ----------------
        ct = pd.DataFrame([r for r in celltype_rows if r["direction"] == direction])
        for model in MODELS:
            sub = ct[ct["model"] == model]
            for col in ("f1", "error_rate", "median_predictive_uncertainty",
                        "median_ensemble_disagreement"):
                res = spearmanr(sub["source_support"], sub[col])
                support_rows.append(
                    {
                        "direction": direction,
                        "model": model,
                        "metric": f"source_support_vs_{col}",
                        "spearman_rho": float(res.statistic),
                        "spearman_p": float(res.pvalue),
                        "n_classes": int(len(sub)),
                    }
                )

        del data, per_model_seed, per_model_neighbors, per_model_agreement, ens_store
        gc.collect()
        print(f"  RSS={rss_mb():.0f} MB", flush=True)

    tables = paths["tables"]
    pd.DataFrame(primary_rows).to_csv(tables / "uncertainty_primary.csv", index=False)
    pd.DataFrame(ensemble_rows).to_csv(tables / "uncertainty_ensemble.csv", index=False)
    pd.DataFrame(celltype_rows).to_csv(tables / "uncertainty_celltype_specific.csv", index=False)
    pd.DataFrame(neighbor_rows).to_csv(tables / "uncertainty_neighbor_stability.csv", index=False)
    pd.DataFrame(delta_rows).to_csv(tables / "uncertainty_scvi_vs_totalvi_delta.csv", index=False)
    pd.DataFrame(reliability_rows).to_csv(tables / "uncertainty_reliability_bins.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(tables / "uncertainty_coverage.csv", index=False)
    pd.DataFrame(support_rows).to_csv(tables / "uncertainty_associations.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(tables / "uncertainty_distributions.csv", index=False)
    save_json(
        {
            "repro_check": repro_log,
            "all_reproduced": True,
            "tolerance": REPRO_TOLERANCE,
            "latent_status": LATENT_UNSUPPORTED_REASON,
            "timestamp_utc": utc_now(),
        },
        paths["logs"] / "phase10_reproduction_check.json",
    )
    with open(paths["logs"] / "risk_coverage_curves.json", "w", encoding="utf-8") as fh:
        json.dump(roc_rows, fh)
    print("\nstage main complete", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", default="main", choices=["main", "latent", "subset"])
    args = parser.parse_args()
    if args.stage == "main":
        run_main()
    else:
        raise SystemExit(f"stage {args.stage} lives in a separate script")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
