#!/usr/bin/env python
"""PHASE 11: formal uncertainty characterization and calibration.

Reuses the PHASE 10 cross-dataset artifacts. Trains no representation model.
The PHASE 10 latent-alignment rule is preserved: source and target are both
read from the post-adaptation embeddings written by PHASE 10.

Four uncertainty concepts are kept separate throughout:
  A. latent posterior uncertainty      (q(z|x) of one fixed trained model)
  B. downstream predictive uncertainty (source-trained classifier on target)
  C. model-seed ensemble uncertainty   (across the five PHASE 10 seeds)
  D. source-subset instability         (across the five PHASE 10B subsamples)
They are never combined into a single reliability score.
"""

from __future__ import annotations

import gc
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)
from sklearn.neighbors import NearestNeighbors

from src.experiments.phase10_cross_dataset import class_plan, load_source_target
from src.utils.io import resolve_path

# ----------------------------------------------------------------------------
# Pre-specified constants. Fixed before any result was inspected.
# ----------------------------------------------------------------------------
DIRECTIONS = (
    ("pbmc10k_to_pbmc5k", "PBMC10k", "PBMC5k"),
    ("pbmc5k_to_pbmc10k", "PBMC5k", "PBMC10k"),
)
MODELS = ("scvi_matched", "totalvi")
SEEDS = (0, 1, 2, 3, 4)
LABEL_LEVELS = (("l2", "cell_type_l2"), ("l1", "cell_type_l1"))

ECE_BINS = 15
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
ENRICHMENT_TOP = (0.10, 0.20)
ENRICHMENT_BOTTOM = 0.20
NEIGHBOR_K = 15
N_BOOT = 500
BOOTSTRAP_SEED = 1111
MC_SAMPLES = 50
POSTERIOR_SAMPLING_SEED = 2026
REPRO_TOLERANCE = 1e-6

# PHASE 10 wrote only the source reference model to disk; the scArches-adapted
# query model was discarded. See posterior_api_validation.md.
LATENT_POSTERIOR_SUPPORTED = False
LATENT_UNSUPPORTED_REASON = (
    "PHASE 10 persisted only the pre-adaptation source reference model "
    "(model.pt). The scArches-adapted query model, which defines the "
    "post-adaptation coordinate system used for every PHASE 10 metric, was "
    "deleted after its latents were written. Posterior variance for target "
    "cells in that coordinate system is therefore not recoverable without "
    "re-running query adaptation, which would be retraining."
)


def _phase11_root() -> Path:
    return resolve_path("results/phase11_uncertainty")


def _phase10_root() -> Path:
    return resolve_path("results/phase10_cross_dataset")


def phase11_paths() -> dict[str, Path]:
    root = _phase11_root()
    return {
        "root": root,
        "tables": root / "tables",
        "figures": root / "figures",
        "logs": root / "logs",
        "arrays": root / "arrays",
    }


def ensure_dirs() -> dict[str, Path]:
    paths = phase11_paths()
    for key, path in paths.items():
        if key != "root":
            path.mkdir(parents=True, exist_ok=True)
    paths["root"].mkdir(parents=True, exist_ok=True)
    return paths


# ----------------------------------------------------------------------------
# PHASE 10 artifact loading
# ----------------------------------------------------------------------------
def phase10_tag(model: str, direction: str, seed: int) -> str:
    return f"{model}_{direction}_seed{seed}_cuda"


def load_phase10_pack(direction: str, model: str, seed: int) -> dict[str, Any]:
    """Load the post-adaptation source/target latents written by PHASE 10."""
    emb = _phase10_root() / "embeddings" / direction
    tag = phase10_tag(model, direction, seed)
    src_path = emb / f"{tag}_source_latent.npy"
    tgt_path = emb / f"{tag}_target_latent.npy"
    for path in (src_path, tgt_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing PHASE 10 artifact: {path}")
    z_src = np.load(src_path)
    z_tgt = np.load(tgt_path)
    src_cells = pd.read_csv(emb / f"{tag}_source_cells.csv")["cell_id"].astype(str).to_numpy()
    tgt_cells = pd.read_csv(emb / f"{tag}_target_cells.csv")["cell_id"].astype(str).to_numpy()
    if not np.isfinite(z_src).all() or not np.isfinite(z_tgt).all():
        raise RuntimeError(f"Non-finite latent values in {tag}")
    if z_src.shape[0] != src_cells.size or z_tgt.shape[0] != tgt_cells.size:
        raise RuntimeError(f"Latent/cell-id length mismatch for {tag}")
    if z_src.shape[1] != z_tgt.shape[1]:
        raise RuntimeError(f"Source/target latent dim mismatch for {tag}")
    return {
        "source_latent": z_src,
        "target_latent": z_tgt,
        "source_cells": src_cells,
        "target_cells": tgt_cells,
        "tag": tag,
        "source_latent_path": str(src_path),
        "target_latent_path": str(tgt_path),
    }


def phase10_reference_metrics(corrected: bool = True) -> pd.DataFrame:
    """PHASE 10 logreg rows that PHASE 11 must reproduce.

    Defaults to the corrected tables. The originals scored totalVI seeds 1-4
    with pre-adaptation source coordinates; see phase10_erratum_latent_alignment.md.
    """
    corrected_path = (
        _phase10_root() / "tables_corrected" / "cross_dataset_primary_corrected.csv"
    )
    table = (
        corrected_path
        if corrected and corrected_path.exists()
        else _phase10_root() / "tables" / "cross_dataset_primary.csv"
    )
    frame = pd.read_csv(table)
    return frame[frame["classifier"] == "logreg"].copy()


class DirectionData:
    """Labels and the PHASE 10 eligible-class plan for one transfer direction."""

    def __init__(self, direction: str, source_batch: str, target_batch: str) -> None:
        source, target, audit = load_source_target(source_batch, target_batch)
        self.direction = direction
        self.audit = audit
        self.source_cells = source.obs_names.astype(str).to_numpy()
        self.target_cells = target.obs_names.astype(str).to_numpy()
        self.plans = {
            level: class_plan(source, target, col) for level, col in LABEL_LEVELS
        }
        self.source_labels = {
            level: source.obs[col].astype(str).to_numpy() for level, col in LABEL_LEVELS
        }
        self.target_labels = {
            level: target.obs[col].astype(str).to_numpy() for level, col in LABEL_LEVELS
        }
        del source, target
        gc.collect()

    def eval_slice(self, level: str) -> dict[str, Any]:
        plan = self.plans[level]
        si = plan["source_eval_index"]
        ti = plan["target_eval_index"]
        return {
            "source_index": si,
            "target_index": ti,
            "y_source": self.source_labels[level][si],
            "y_target": self.target_labels[level][ti],
            "source_cells": self.source_cells[si],
            "target_cells": self.target_cells[ti],
            "classes": plan["eligible_shared_classes"],
        }


def verify_cell_alignment(data: DirectionData, pack: dict[str, Any]) -> None:
    """PHASE 10 saved latents row-for-row with the AnnData it embedded."""
    if not np.array_equal(data.source_cells, pack["source_cells"]):
        raise RuntimeError(f"Source cell order differs from PHASE 10 for {pack['tag']}")
    if not np.array_equal(data.target_cells, pack["target_cells"]):
        raise RuntimeError(f"Target cell order differs from PHASE 10 for {pack['tag']}")


# ----------------------------------------------------------------------------
# Classifier: byte-for-byte the PHASE 10 configuration
# ----------------------------------------------------------------------------
def fit_transfer_classifier(
    z_src: np.ndarray, y_src: np.ndarray, z_tgt: np.ndarray, seed: int
) -> dict[str, Any]:
    """PHASE 10 used LogisticRegression(max_iter=2000, lbfgs, random_state=seed)."""
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed, class_weight=None)
    clf.fit(z_src, y_src)
    classes = np.asarray(clf.classes_, dtype=object)
    probs = clf.predict_proba(z_tgt)
    pred = classes[np.argmax(probs, axis=1)]
    return {"classifier": clf, "classes": classes, "probs": probs, "pred": pred}


def classification_metrics(y_true: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y_true, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "macro_f1": float(f1_score(y_true, pred, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_true, pred, average="weighted", zero_division=0)),
    }


# ----------------------------------------------------------------------------
# Predictive uncertainty (concept B)
# ----------------------------------------------------------------------------
def predictive_uncertainty(probs: np.ndarray) -> dict[str, np.ndarray]:
    """Primary U_pred = 1 - max class probability. Entropy and margin secondary."""
    clipped = np.clip(probs, 1e-12, 1.0)
    ordered = np.sort(probs, axis=1)
    p_max = ordered[:, -1]
    p_second = ordered[:, -2] if probs.shape[1] > 1 else np.zeros_like(p_max)
    return {
        "u_pred": 1.0 - p_max,
        "p_max": p_max,
        "entropy": -np.sum(clipped * np.log(clipped), axis=1),
        "margin": p_max - p_second,
    }


def entropy_of(probs: np.ndarray) -> np.ndarray:
    clipped = np.clip(probs, 1e-12, 1.0)
    return -np.sum(clipped * np.log(clipped), axis=1)


# ----------------------------------------------------------------------------
# Discrimination: can uncertainty rank unreliable cells?
# ----------------------------------------------------------------------------
def error_detection(uncertainty: np.ndarray, error: np.ndarray) -> dict[str, float]:
    """AUROC/AUPRC for flagging incorrect predictions. Prevalence always reported."""
    prevalence = float(np.mean(error))
    if prevalence <= 0.0 or prevalence >= 1.0:
        return {
            "auroc": float("nan"),
            "auprc": float("nan"),
            "error_prevalence": prevalence,
            "n": int(error.size),
        }
    return {
        "auroc": float(roc_auc_score(error, uncertainty)),
        "auprc": float(average_precision_score(error, uncertainty)),
        "error_prevalence": prevalence,
        "n": int(error.size),
    }


def risk_coverage(uncertainty: np.ndarray, correct: np.ndarray) -> dict[str, Any]:
    """Retain the most confident cells first; risk is 1 - accuracy of the kept set."""
    order = np.argsort(uncertainty, kind="mergesort")
    kept_correct = correct[order].astype(float)
    running = np.cumsum(kept_correct)
    n = correct.size
    counts = np.arange(1, n + 1)
    risks = 1.0 - running / counts
    aurc = float(np.mean(risks))
    return {"aurc": aurc, "coverage": counts / n, "risk": risks, "order": order}


def coverage_table(
    uncertainty: np.ndarray, y_true: np.ndarray, pred: np.ndarray
) -> list[dict[str, float]]:
    """Selective prediction at the pre-specified coverage grid."""
    order = np.argsort(uncertainty, kind="mergesort")
    n = y_true.size
    rows = []
    for coverage in COVERAGE_GRID:
        keep = max(1, int(round(coverage * n)))
        idx = order[:keep]
        rows.append(
            {
                "coverage": float(coverage),
                "n_retained": int(keep),
                "accuracy": float(accuracy_score(y_true[idx], pred[idx])),
                "risk": float(1.0 - accuracy_score(y_true[idx], pred[idx])),
                "macro_f1": float(f1_score(y_true[idx], pred[idx], average="macro", zero_division=0)),
                "balanced_accuracy": float(balanced_accuracy_score(y_true[idx], pred[idx])),
            }
        )
    return rows


def error_enrichment(uncertainty: np.ndarray, error: np.ndarray) -> dict[str, float]:
    """Pre-specified 10%/20% top and bottom-20% uncertainty strata."""
    n = error.size
    order = np.argsort(uncertainty, kind="mergesort")
    out: dict[str, float] = {}
    for frac in ENRICHMENT_TOP:
        k = max(1, int(round(frac * n)))
        top_idx = order[-k:]
        out[f"error_rate_top_{int(frac * 100)}pct"] = float(np.mean(error[top_idx]))
    k_bottom = max(1, int(round(ENRICHMENT_BOTTOM * n)))
    bottom_idx = order[:k_bottom]
    bottom_rate = float(np.mean(error[bottom_idx]))
    out["error_rate_bottom_20pct"] = bottom_rate
    top10 = out["error_rate_top_10pct"]
    out["error_enrichment_top10_over_bottom20"] = (
        float(top10 / bottom_rate) if bottom_rate > 0 else float("inf")
    )
    return out


# ----------------------------------------------------------------------------
# Calibration: do stated probabilities match empirical correctness?
# ----------------------------------------------------------------------------
def _one_hot(y_true: np.ndarray, classes: np.ndarray) -> np.ndarray:
    lookup = {c: i for i, c in enumerate(classes)}
    out = np.zeros((y_true.size, classes.size), dtype=float)
    for row, label in enumerate(y_true):
        col = lookup.get(label)
        if col is not None:
            out[row, col] = 1.0
    return out


def expected_calibration_error(
    confidence: np.ndarray, correct: np.ndarray, bins: int, scheme: str
) -> dict[str, Any]:
    """ECE with equal-frequency (primary) or equal-width (sensitivity) bins."""
    n = confidence.size
    if scheme == "equal_frequency":
        quantiles = np.quantile(confidence, np.linspace(0.0, 1.0, bins + 1))
        edges = np.unique(quantiles)
    elif scheme == "equal_width":
        edges = np.linspace(0.0, 1.0, bins + 1)
    else:
        raise ValueError(scheme)
    assignment = np.clip(np.digitize(confidence, edges[1:-1], right=True), 0, len(edges) - 2)
    ece = 0.0
    rows = []
    for b in range(len(edges) - 1):
        mask = assignment == b
        count = int(mask.sum())
        if count == 0:
            continue
        mean_conf = float(confidence[mask].mean())
        mean_acc = float(correct[mask].mean())
        ece += (count / n) * abs(mean_conf - mean_acc)
        rows.append(
            {
                "bin": b,
                "bin_lower": float(edges[b]),
                "bin_upper": float(edges[b + 1]),
                "n_cells": count,
                "mean_confidence": mean_conf,
                "empirical_accuracy": mean_acc,
                "gap": mean_acc - mean_conf,
            }
        )
    return {"ece": float(ece), "bins": rows, "scheme": scheme, "n_bins_used": len(rows)}


def calibration_metrics(
    probs: np.ndarray, y_true: np.ndarray, classes: np.ndarray, pred: np.ndarray
) -> dict[str, Any]:
    """NLL, multiclass Brier, and both ECE variants."""
    onehot = _one_hot(y_true, classes)
    clipped = np.clip(probs, 1e-12, 1.0)
    true_prob = np.sum(onehot * clipped, axis=1)
    nll = float(-np.mean(np.log(np.clip(true_prob, 1e-12, 1.0))))
    brier = float(np.mean(np.sum((probs - onehot) ** 2, axis=1)))
    confidence = probs.max(axis=1)
    correct = (pred == y_true).astype(float)
    eq_freq = expected_calibration_error(confidence, correct, ECE_BINS, "equal_frequency")
    eq_width = expected_calibration_error(confidence, correct, ECE_BINS, "equal_width")
    return {
        "nll": nll,
        "brier": brier,
        "true_class_probability": true_prob,
        "ece_equal_frequency": eq_freq["ece"],
        "ece_equal_width": eq_width["ece"],
        "reliability_equal_frequency": eq_freq["bins"],
        "reliability_equal_width": eq_width["bins"],
    }


# ----------------------------------------------------------------------------
# Model-seed ensemble (concept C)
# ----------------------------------------------------------------------------
def align_probability_matrices(
    per_seed: dict[int, dict[str, Any]]
) -> tuple[np.ndarray, np.ndarray]:
    """Stack per-seed probabilities after aligning columns by class NAME.

    PHASE 10 fit an independent classifier per seed, so column order is only
    guaranteed to be sklearn's sorted class order. Alignment is explicit and
    verified rather than assumed.
    """
    seeds = sorted(per_seed)
    reference = list(per_seed[seeds[0]]["classes"])
    stacked = []
    for seed in seeds:
        classes = list(per_seed[seed]["classes"])
        probs = per_seed[seed]["probs"]
        if classes == reference:
            stacked.append(probs)
            continue
        if set(classes) != set(reference):
            raise RuntimeError(
                f"Seed {seed} class set differs from reference; cannot align probabilities."
            )
        order = [classes.index(c) for c in reference]
        stacked.append(probs[:, order])
    return np.stack(stacked, axis=0), np.asarray(reference, dtype=object)


def ensemble_uncertainty(stacked: np.ndarray, classes: np.ndarray) -> dict[str, Any]:
    """Ensemble mean, predictive/expected entropy, disagreement, variation ratio."""
    p_bar = stacked.mean(axis=0)
    predictive_entropy = entropy_of(p_bar)
    member_entropy = np.stack([entropy_of(stacked[s]) for s in range(stacked.shape[0])], axis=0)
    expected_entropy = member_entropy.mean(axis=0)
    disagreement = predictive_entropy - expected_entropy
    member_pred = np.argmax(stacked, axis=2)
    n_seeds, n_cells = member_pred.shape
    variation = np.empty(n_cells, dtype=float)
    for i in range(n_cells):
        counts = np.bincount(member_pred[:, i], minlength=classes.size)
        variation[i] = 1.0 - counts.max() / n_seeds
    return {
        "p_bar": p_bar,
        "ensemble_pred": classes[np.argmax(p_bar, axis=1)],
        "predictive_entropy": predictive_entropy,
        "expected_entropy": expected_entropy,
        "disagreement": disagreement,
        "variation_ratio": variation,
    }


# ----------------------------------------------------------------------------
# Neighborhood stability across model seeds
# ----------------------------------------------------------------------------
def source_neighbor_sets(
    z_src: np.ndarray, z_tgt: np.ndarray, k: int = NEIGHBOR_K
) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(z_src)
    return nn.kneighbors(z_tgt, return_distance=False)


def neighbor_stability(neighbor_sets: list[np.ndarray]) -> np.ndarray:
    """Mean pairwise Jaccard overlap of each cell's source-neighbor set across seeds."""
    n_seeds = len(neighbor_sets)
    n_cells = neighbor_sets[0].shape[0]
    stability = np.zeros(n_cells, dtype=float)
    pairs = [(a, b) for a in range(n_seeds) for b in range(a + 1, n_seeds)]
    if not pairs:
        return np.full(n_cells, np.nan)
    for i in range(n_cells):
        total = 0.0
        for a, b in pairs:
            set_a = set(neighbor_sets[a][i].tolist())
            set_b = set(neighbor_sets[b][i].tolist())
            union = len(set_a | set_b)
            total += len(set_a & set_b) / union if union else 0.0
        stability[i] = total / len(pairs)
    return stability


def neighbor_label_agreement(
    neighbor_idx: np.ndarray, y_src: np.ndarray, y_tgt: np.ndarray
) -> np.ndarray:
    return np.array(
        [float(np.mean(y_src[neigh] == y_tgt[i])) for i, neigh in enumerate(neighbor_idx)]
    )


# ----------------------------------------------------------------------------
# Bootstrap
# ----------------------------------------------------------------------------
def paired_bootstrap(
    metric_fn,
    n: int,
    n_boot: int = N_BOOT,
    seed: int = BOOTSTRAP_SEED,
) -> np.ndarray:
    """Identical resample indices are reused across models by fixing the seed."""
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        draws[b] = metric_fn(idx)
    return draws


def bootstrap_indices(n: int, n_boot: int = N_BOOT, seed: int = BOOTSTRAP_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.stack([rng.choice(n, size=n, replace=True) for _ in range(n_boot)], axis=0)


def ci_bounds(draws: np.ndarray) -> tuple[float, float]:
    finite = draws[np.isfinite(draws)]
    if finite.size == 0:
        return float("nan"), float("nan")
    return float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()
