#!/usr/bin/env python
"""PHASE 11B pure statistics: no AnnData, no torch, no scvi-tools.

Isolated from the training and I/O layers so the numerical core can be unit
tested on any machine, including one without a GPU stack. See tests/test_phase11b.py.
"""

from __future__ import annotations

import numpy as np

# Evidence rule (§22). Fixed before any PHASE 11B result was inspected.
EVIDENCE_RULE = {
    "robust": (
        "seed-level paired bootstrap 95% CI excludes 0 AND at least 16/20 seeds "
        "share the sign of the mean delta"
    ),
    "seed_sensitive": (
        "the majority of per-seed target-cell bootstrap CIs exclude 0, but fewer "
        "than 16/20 seeds agree in sign OR the seed-level CI includes 0: the "
        "apparent effect depends on which model seed was drawn"
    ),
    "cell_sampling_sensitive": (
        "seed-level CI excludes 0 and signs are consistent, but fewer than half "
        "of the per-seed target-cell bootstrap CIs exclude 0: the effect is "
        "stable across seeds yet small relative to target-cell sampling noise"
    ),
    "inconclusive": "neither the seed-level nor the target-cell evidence separates the models",
    "sign_consistency_threshold": 16,
    "n_seeds": 20,
}


def classify_evidence(seed_ci: tuple[float, float], n_positive: int, n_seeds: int,
                      frac_cell_ci_excluding_zero: float) -> str:
    lo, hi = seed_ci
    seed_ci_excludes = np.isfinite(lo) and np.isfinite(hi) and (lo > 0 or hi < 0)
    n_agree = max(n_positive, n_seeds - n_positive)
    signs_consistent = n_agree >= EVIDENCE_RULE["sign_consistency_threshold"]
    if seed_ci_excludes and signs_consistent:
        return "robust" if frac_cell_ci_excluding_zero >= 0.5 else "cell_sampling_sensitive"
    if frac_cell_ci_excluding_zero >= 0.5:
        return "seed_sensitive"
    return "inconclusive"


# ---------------------------------------------------------------------------
# Regression guard for the PHASE 10 latent-alignment bug
# ---------------------------------------------------------------------------
def assert_post_adaptation(
    source_latent_used: np.ndarray,
    source_latent_pre_adaptation: np.ndarray,
    source_latent_post_adaptation: np.ndarray,
    tag: str,
) -> None:
    """Guarantee the evaluator receives the post-adaptation source embedding.

    PHASE 10's `run_one_model` returned the pre-adaptation array while writing
    the post-adaptation one to disk, so freshly-trained runs mixed coordinate
    systems. This raises rather than letting that recur silently.

    The identity check is the load-bearing one: it pins the exact object, so a
    same-valued copy is rejected too. When adaptation leaves the reference
    embedding untouched (scArches freezes it for scVI) pre and post coincide,
    and passing that shared object is legitimate.
    """
    if source_latent_used is not source_latent_post_adaptation:
        raise RuntimeError(
            f"{tag}: downstream source embedding is not the post-adaptation object."
        )
    if source_latent_used.shape != source_latent_post_adaptation.shape:
        raise RuntimeError(f"{tag}: post-adaptation source embedding has the wrong shape.")
    if not np.array_equal(source_latent_used, source_latent_post_adaptation):
        raise RuntimeError(f"{tag}: source embedding content differs from post-adaptation array.")
    if source_latent_pre_adaptation.shape != source_latent_post_adaptation.shape:
        raise RuntimeError(f"{tag}: pre/post adaptation source embeddings differ in shape.")


def latent_drift(pre: np.ndarray, post: np.ndarray) -> dict[str, float]:
    """Diagnostic only. Zero drift is not assumed."""
    shift = np.linalg.norm(post - pre, axis=1)
    pre_norm = float(np.linalg.norm(pre))
    return {
        "mean_cell_displacement": float(shift.mean()),
        "median_cell_displacement": float(np.median(shift)),
        "max_cell_displacement": float(shift.max()),
        "relative_frobenius_change": (
            float(np.linalg.norm(post - pre) / pre_norm) if pre_norm else float("nan")
        ),
        "identical": bool(np.array_equal(pre, post)),
    }


# ---------------------------------------------------------------------------
# Balanced two-factor variance components (§31)
# ---------------------------------------------------------------------------
def two_factor_components(matrix: np.ndarray) -> dict[str, float]:
    """Method-of-moments components for an a x b table with one observation per cell.

    Rows are source subsets, columns are model seeds. With n=1 per cell the
    interaction cannot be separated from pure error, so the third component is
    reported as combined interaction/residual. Raw (possibly negative) estimates
    are preserved alongside the non-negative truncation.
    """
    y = np.asarray(matrix, dtype=float)
    a, b = y.shape
    if a < 2 or b < 2:
        raise ValueError("need at least a 2x2 table")
    grand = y.mean()
    row_means = y.mean(axis=1)
    col_means = y.mean(axis=0)
    ms_a = b * float(np.sum((row_means - grand) ** 2)) / (a - 1)
    ms_b = a * float(np.sum((col_means - grand) ** 2)) / (b - 1)
    resid = y - row_means[:, None] - col_means[None, :] + grand
    ms_e = float(np.sum(resid**2)) / ((a - 1) * (b - 1))
    var_a_raw = (ms_a - ms_e) / b
    var_b_raw = (ms_b - ms_e) / a
    var_i_raw = ms_e
    nn = {k: max(0.0, v) for k, v in
          (("a", var_a_raw), ("b", var_b_raw), ("i", var_i_raw))}
    total = sum(nn.values())
    return {
        "source_subset_variance_raw": var_a_raw,
        "model_seed_variance_raw": var_b_raw,
        "interaction_variance_raw": var_i_raw,
        "source_subset_variance_nonnegative": nn["a"],
        "model_seed_variance_nonnegative": nn["b"],
        "interaction_variance_nonnegative": nn["i"],
        "total_variance_nonnegative": total,
        "source_fraction": nn["a"] / total if total > 0 else float("nan"),
        "seed_fraction": nn["b"] / total if total > 0 else float("nan"),
        "interaction_fraction": nn["i"] / total if total > 0 else float("nan"),
        "ms_source_subset": ms_a,
        "ms_model_seed": ms_b,
        "ms_interaction_residual": ms_e,
        "grand_mean": float(grand),
    }


def vector_components(cube: np.ndarray) -> dict[str, np.ndarray]:
    """Vectorised two-factor decomposition over the last axis.

    ``cube`` is (a source subsets, b model seeds, n cells). Equivalent to calling
    :func:`two_factor_components` once per cell, but without the Python loop.
    """
    a, b, _ = cube.shape
    y = np.moveaxis(cube, 2, 0)  # (n, a, b)
    grand = y.mean(axis=(1, 2))
    row = y.mean(axis=2)
    col = y.mean(axis=1)
    ms_a = b * np.sum((row - grand[:, None]) ** 2, axis=1) / (a - 1)
    ms_b = a * np.sum((col - grand[:, None]) ** 2, axis=1) / (b - 1)
    resid = y - row[:, :, None] - col[:, None, :] + grand[:, None, None]
    ms_e = np.sum(resid**2, axis=(1, 2)) / ((a - 1) * (b - 1))
    va = np.maximum((ms_a - ms_e) / b, 0.0)
    vb = np.maximum((ms_b - ms_e) / a, 0.0)
    vi = np.maximum(ms_e, 0.0)
    total = va + vb + vi
    with np.errstate(invalid="ignore", divide="ignore"):
        return {"source_var": va, "seed_var": vb, "interaction_var": vi,
                "source_fraction": np.where(total > 0, va / total, np.nan),
                "seed_fraction": np.where(total > 0, vb / total, np.nan),
                "interaction_fraction": np.where(total > 0, vi / total, np.nan)}


# ---------------------------------------------------------------------------
# Per-cell instability measures
# ---------------------------------------------------------------------------
def variation_ratio(pred_labels: np.ndarray) -> np.ndarray:
    """1 - modal frequency, computed along axis 0 of a (n_runs, n_cells) label array."""
    n_runs, n_cells = pred_labels.shape
    out = np.empty(n_cells, dtype=float)
    for i in range(n_cells):
        _, counts = np.unique(pred_labels[:, i], return_counts=True)
        out[i] = 1.0 - counts.max() / n_runs
    return out


def mean_pairwise_distance(prob_stack: np.ndarray) -> np.ndarray:
    """Root-mean-square pairwise distance between probability vectors, per cell.

    ``prob_stack`` is (n_runs, n_cells, n_classes). Computed from the centred
    second moment rather than by enumerating the O(n_runs^2) pairs.
    """
    n_runs = prob_stack.shape[0]
    centered = prob_stack - prob_stack.mean(axis=0, keepdims=True)
    msd = np.sum(centered**2, axis=2).mean(axis=0)
    return np.sqrt(2.0 * msd * n_runs / max(n_runs - 1, 1))


def one_hot_matrix(y_true: np.ndarray, classes: np.ndarray) -> np.ndarray:
    """Dense one-hot for the full evaluation set, built once and then resampled.

    Labels outside ``classes`` produce an all-zero row, matching the PHASE 11
    convention rather than raising.
    """
    lookup = {c: i for i, c in enumerate(classes)}
    out = np.zeros((y_true.size, classes.size), dtype=float)
    cols = np.array([lookup.get(lab, -1) for lab in y_true])
    known = cols >= 0
    out[np.arange(y_true.size)[known], cols[known]] = 1.0
    return out
