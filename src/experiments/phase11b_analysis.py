#!/usr/bin/env python
"""PHASE 11B analysis: 20-seed stability, posterior recovery, variance decomposition.

Consumes only PHASE 11B artifacts. Metric definitions are imported from PHASE 11
so the two phases are numerically comparable; nothing is redefined here.
"""

from __future__ import annotations

import gc
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.experiments.phase10_cross_dataset import class_plan, load_source_target
from src.experiments.phase11_uncertainty import (
    ECE_BINS,
    calibration_metrics,
    classification_metrics,
    coverage_table,
    error_detection,
    error_enrichment,
    expected_calibration_error,
    fit_transfer_classifier,
    predictive_uncertainty,
    risk_coverage,
)
from src.experiments.phase11b_stats import (
    EVIDENCE_RULE,
    classify_evidence,
    mean_pairwise_distance,
    one_hot_matrix,
    two_factor_components,
    variation_ratio,
)
from src.experiments.phase11b_robustness import (
    DIRECTIONS,
    experiment_tag,
    phase11b_paths,
)

LABEL_LEVELS = (("l2", "cell_type_l2"), ("l1", "cell_type_l1"))
PRIMARY_LEVEL = "l2"

# ---------------------------------------------------------------------------
# Artifact loading
# ---------------------------------------------------------------------------
class DirectionEval:
    """Labels and the fixed eligible-class plan for one transfer direction."""

    def __init__(self, direction_key: str) -> None:
        _, direction, source_batch, target_batch = next(
            d for d in DIRECTIONS if d[0] == direction_key)
        source, target, audit = load_source_target(source_batch, target_batch)
        self.direction_key = direction_key
        self.direction = direction
        self.audit = audit
        self.source_cells = source.obs_names.astype(str).to_numpy()
        self.target_cells = target.obs_names.astype(str).to_numpy()
        self.plans = {lv: class_plan(source, target, col) for lv, col in LABEL_LEVELS}
        self.source_labels = {lv: source.obs[col].astype(str).to_numpy()
                              for lv, col in LABEL_LEVELS}
        self.target_labels = {lv: target.obs[col].astype(str).to_numpy()
                              for lv, col in LABEL_LEVELS}
        del source, target
        gc.collect()

    def slice(self, level: str = PRIMARY_LEVEL) -> dict[str, Any]:
        plan = self.plans[level]
        si, ti = plan["source_eval_index"], plan["target_eval_index"]
        return {
            "source_index": si,
            "target_index": ti,
            "y_source": self.source_labels[level][si],
            "y_target": self.target_labels[level][ti],
            "target_cells": self.target_cells[ti],
            "classes": plan["eligible_shared_classes"],
        }


def load_run(direction_key: str, model: str, model_seed: int,
             subset_seed: int | None = None) -> dict[str, Any]:
    """Post-adaptation latents and target posterior for one PHASE 11B run."""
    paths = phase11b_paths()
    tag = experiment_tag(direction_key, model, model_seed, subset_seed)
    emb = paths["embeddings"] / direction_key
    post = paths["posterior"] / direction_key
    z_src = np.load(emb / f"{tag}_source_latent_post.npy")
    z_tgt = np.load(emb / f"{tag}_target_latent_post.npy")
    with np.load(post / f"{tag}_posterior.npz", allow_pickle=True) as npz:
        u_latent = npz["u_latent_target"]
        target_var = npz["target_var"]
    if not np.isfinite(z_src).all() or not np.isfinite(z_tgt).all():
        raise RuntimeError(f"{tag}: non-finite latent")
    return {
        "tag": tag,
        "source_latent": z_src,
        "target_latent": z_tgt,
        "u_latent_target": u_latent,
        "target_var": target_var,
        "source_cells": pd.read_csv(emb / f"{tag}_source_cells.csv")["cell_id"].astype(str).to_numpy(),
        "target_cells": pd.read_csv(emb / f"{tag}_target_cells.csv")["cell_id"].astype(str).to_numpy(),
    }


def probability_path(direction_key: str, tag: str) -> Path:
    out = phase11b_paths()["probabilities"] / direction_key
    out.mkdir(parents=True, exist_ok=True)
    return out / f"{tag}_probs.npz"


def save_probabilities(direction_key: str, tag: str, probs: np.ndarray,
                       classes: np.ndarray, cell_ids: np.ndarray,
                       u_latent: np.ndarray, meta: dict[str, Any]) -> None:
    """Class order is stored explicitly; never aggregate without checking it."""
    np.savez_compressed(
        probability_path(direction_key, tag),
        probs=probs.astype(np.float32),
        classes=np.asarray(classes, dtype=object),
        cell_ids=np.asarray(cell_ids, dtype=object),
        u_latent=np.asarray(u_latent, dtype=np.float32),
        meta=np.asarray([str(meta)], dtype=object),
    )


def load_probabilities(direction_key: str, tag: str) -> dict[str, Any]:
    with np.load(probability_path(direction_key, tag), allow_pickle=True) as npz:
        return {
            "probs": npz["probs"].astype(float),
            "classes": np.asarray(npz["classes"], dtype=object),
            "cell_ids": np.asarray(npz["cell_ids"], dtype=object),
            "u_latent": npz["u_latent"].astype(float),
        }


# ---------------------------------------------------------------------------
# Scoring one run
# ---------------------------------------------------------------------------
def score_run(z_src: np.ndarray, y_src: np.ndarray, z_tgt: np.ndarray,
              y_tgt: np.ndarray, u_latent: np.ndarray, seed: int) -> dict[str, Any]:
    """Fit the PHASE 10 classifier on SOURCE only, then score TARGET."""
    fit = fit_transfer_classifier(z_src, y_src, z_tgt, seed)
    probs, classes, pred = fit["probs"], fit["classes"], fit["pred"]
    correct = (pred == y_tgt)
    error = (~correct).astype(int)
    unc = predictive_uncertainty(probs)
    cal = calibration_metrics(probs, y_tgt, classes, pred)
    rc = risk_coverage(unc["u_pred"], correct.astype(float))
    pred_err = error_detection(unc["u_pred"], error)
    lat_err = error_detection(u_latent, error)
    neg_log_true = -np.log(np.clip(cal["true_class_probability"], 1e-12, 1.0))

    from scipy.stats import spearmanr
    rho_pred_latent = spearmanr(unc["u_pred"], u_latent)
    rho_latent_nlt = spearmanr(u_latent, neg_log_true)
    lat_rc = risk_coverage(u_latent, correct.astype(float))

    row: dict[str, Any] = {
        **classification_metrics(y_tgt, pred),
        "n_cells": int(y_tgt.size),
        "error_prevalence": pred_err["error_prevalence"],
        "error_auroc_predictive": pred_err["auroc"],
        "error_auprc_predictive": pred_err["auprc"],
        "error_auroc_latent": lat_err["auroc"],
        "error_auprc_latent": lat_err["auprc"],
        "nll": cal["nll"],
        "brier": cal["brier"],
        "ece_equal_frequency": cal["ece_equal_frequency"],
        "ece_equal_width": cal["ece_equal_width"],
        "aurc": rc["aurc"],
        "aurc_latent": lat_rc["aurc"],
        "median_predictive_uncertainty": float(np.median(unc["u_pred"])),
        "median_latent_uncertainty": float(np.median(u_latent)),
        "spearman_upred_ulatent": float(rho_pred_latent.statistic),
        "spearman_upred_ulatent_p": float(rho_pred_latent.pvalue),
        "spearman_ulatent_neglog_trueclass": float(rho_latent_nlt.statistic),
        "spearman_ulatent_neglog_trueclass_p": float(rho_latent_nlt.pvalue),
        "median_latent_uncertainty_correct": float(np.median(u_latent[correct])),
        "median_latent_uncertainty_incorrect": float(np.median(u_latent[~correct])),
        "median_predictive_uncertainty_correct": float(np.median(unc["u_pred"][correct])),
        "median_predictive_uncertainty_incorrect": float(np.median(unc["u_pred"][~correct])),
        **error_enrichment(unc["u_pred"], error),
    }
    return {
        "row": row,
        "probs": probs,
        "classes": classes,
        "pred": pred,
        "correct": correct,
        "error": error,
        "u_pred": unc["u_pred"],
        "entropy": unc["entropy"],
        "true_class_probability": cal["true_class_probability"],
        "coverage_rows": coverage_table(unc["u_pred"], y_tgt, pred),
        "reliability_equal_frequency": cal["reliability_equal_frequency"],
    }


# ---------------------------------------------------------------------------
# Bootstrap metric block used for target-cell CIs
# ---------------------------------------------------------------------------
def metric_block(probs: np.ndarray, classes: np.ndarray, y_true: np.ndarray,
                 idx: np.ndarray, onehot_full: np.ndarray | None = None) -> dict[str, float]:
    """Five paired metrics on one bootstrap resample.

    ``onehot_full`` is the one-hot over the full evaluation set. Passing it in
    keeps the per-replicate cost to an array index instead of rebuilding the
    encoding for every one of the thousands of resamples.
    """
    p = probs[idx]
    y = y_true[idx]
    pred = classes[np.argmax(p, axis=1)]
    correct = (pred == y)
    error = (~correct).astype(int)
    u = 1.0 - p.max(axis=1)
    det = error_detection(u, error)
    onehot = (one_hot_matrix(y, classes) if onehot_full is None else onehot_full[idx])
    brier = float(np.mean(np.sum((p - onehot) ** 2, axis=1)))
    ece = expected_calibration_error(p.max(axis=1), correct.astype(float),
                                     ECE_BINS, "equal_frequency")["ece"]
    aurc = risk_coverage(u, correct.astype(float))["aurc"]
    return {
        "error_auroc_predictive": det["auroc"],
        "error_auprc_predictive": det["auprc"],
        "ece_equal_frequency": ece,
        "brier": brier,
        "aurc": aurc,
    }


BOOTSTRAP_METRICS = ("error_auroc_predictive", "error_auprc_predictive",
                     "ece_equal_frequency", "brier", "aurc")
HIGHER_IS_BETTER = {"error_auroc_predictive": True, "error_auprc_predictive": True,
                    "ece_equal_frequency": False, "brier": False, "aurc": False,
                    "nll": False, "macro_f1": True, "accuracy": True,
                    "error_auroc_latent": True, "error_auprc_latent": True}

__all__ = ["EVIDENCE_RULE", "classify_evidence", "two_factor_components",
           "variation_ratio", "mean_pairwise_distance", "one_hot_matrix"]
