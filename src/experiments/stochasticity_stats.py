"""Pure-NumPy variance decomposition for the PHASE 12A replicated design.

PHASE 11B could only fit a two-factor model with one observation per cell, in
which the subset x seed interaction and the pure residual are the same number.
With a second replicate the design becomes a balanced two-factor random-effects
layout WITH replication, and the two are separately identified.

Model (a = source subsets, b = model seeds, n = replicates per cell):

    y_ijk = mu + A_i + B_j + (AB)_ij + e_ijk

with all effects random and independent,
A_i ~ N(0, s2_A), B_j ~ N(0, s2_B), (AB)_ij ~ N(0, s2_AB), e_ijk ~ N(0, s2_e).

Mean squares:

    MS_A  = b n SUM_i (ybar_i.. - ybar...)^2 / (a-1)
    MS_B  = a n SUM_j (ybar_.j. - ybar...)^2 / (b-1)
    MS_AB = n SUM_ij (ybar_ij. - ybar_i.. - ybar_.j. + ybar...)^2 / ((a-1)(b-1))
    MS_E  = SUM_ijk (y_ijk - ybar_ij.)^2 / (a b (n-1))

Expected mean squares for the random-effects model:

    E[MS_E]  = s2_e
    E[MS_AB] = s2_e + n s2_AB
    E[MS_A]  = s2_e + n s2_AB + b n s2_A
    E[MS_B]  = s2_e + n s2_AB + a n s2_B

giving the method-of-moments estimators

    s2_e  = MS_E
    s2_AB = (MS_AB - MS_E) / n
    s2_A  = (MS_A  - MS_AB) / (b n)
    s2_B  = (MS_B  - MS_AB) / (a n)

Raw estimates are returned unchanged, including negative values, alongside a
non-negative truncated version. Truncating silently would hide exactly the
sampling noise these estimates are meant to expose.
"""

from __future__ import annotations

import numpy as np

COMPONENT_NAMES = ("source_subset", "model_seed", "interaction", "residual")

# Below this total the cube is constant up to floating-point dust.
_ZERO_VARIANCE_EPS = 1e-24

VARIANCE_MODEL_NOTE = (
    "Balanced two-factor random-effects ANOVA with n replicates per cell "
    "(Searle, Casella & McCulloch, Variance Components, ch. 4). Method-of-"
    "moments estimators from expected mean squares. Interaction and pure "
    "residual are separately identified only because n > 1; with n = 1 (the "
    "PHASE 11B design) MS_E is undefined and the two collapse into one term."
)


TRUNCATION_BIAS_NOTE = (
    "Method-of-moments component estimates are unbiased on the raw scale but "
    "fluctuate around zero for a null factor. After truncation at zero a null "
    "factor therefore receives a strictly positive fraction, and with a single "
    "5x5x2 table that upward bias is large: on pure noise the residual "
    "typically retains only about half of the truncated total. Fractions are "
    "descriptive; the F-tests, not the fractions, decide whether a component "
    "is real."
)


def _f_sf(f: float, dfn: int, dfd: int) -> float:
    """Upper-tail F p-value; NaN if scipy is unavailable (pure-NumPy fallback)."""
    if not np.isfinite(f) or f <= 0:
        return float("nan")
    try:
        from scipy.stats import f as f_dist
    except ImportError:
        return float("nan")
    return float(f_dist.sf(f, dfn, dfd))


def replicated_components(cube: np.ndarray) -> dict[str, float]:
    """Decompose one balanced (a, b, n) array into four variance components.

    Parameters
    ----------
    cube
        Shape (a, b, n): source subsets x model seeds x replicates. n must be
        at least 2, otherwise the residual is not identified.
    """
    arr = np.asarray(cube, dtype=float)
    if arr.ndim != 3:
        raise ValueError(f"expected a 3-D (a, b, n) array, got shape {arr.shape}")
    a, b, n = arr.shape
    if n < 2:
        raise ValueError(
            "n must be >= 2 to separate interaction from residual; with n = 1 "
            "use phase11b_stats.two_factor_components, which reports the two "
            "as a single combined term"
        )
    if a < 2 or b < 2:
        raise ValueError(f"need at least 2 levels per factor, got a={a}, b={b}")
    if not np.isfinite(arr).all():
        raise ValueError("cube contains NaN or Inf")

    cell = arr.mean(axis=2)                 # (a, b) mean over replicates
    grand = arr.mean()
    row = cell.mean(axis=1)                 # (a,)
    col = cell.mean(axis=0)                 # (b,)

    ms_a = b * n * ((row - grand) ** 2).sum() / (a - 1)
    ms_b = a * n * ((col - grand) ** 2).sum() / (b - 1)
    inter = cell - row[:, None] - col[None, :] + grand
    ms_ab = n * (inter ** 2).sum() / ((a - 1) * (b - 1))
    ms_e = ((arr - cell[:, :, None]) ** 2).sum() / (a * b * (n - 1))

    raw = {
        "source_subset": (ms_a - ms_ab) / (b * n),
        "model_seed": (ms_b - ms_ab) / (a * n),
        "interaction": (ms_ab - ms_e) / n,
        "residual": ms_e,
    }
    nonneg = {k: max(v, 0.0) for k, v in raw.items()}
    total = sum(nonneg.values())
    # A constant cube leaves float dust rather than an exact zero; treat any
    # total this small as no variance at all instead of reporting fractions of
    # rounding error.
    degenerate = total <= _ZERO_VARIANCE_EPS

    out: dict[str, float] = {}
    for name in COMPONENT_NAMES:
        out[f"{name}_variance_raw"] = float(raw[name])
        out[f"{name}_variance_nonnegative"] = float(nonneg[name])
        out[f"{name}_fraction"] = 0.0 if degenerate else float(nonneg[name] / total)
    # Random-effects F-tests. The interaction test is the one PHASE 11B could
    # not run: MS_AB / MS_E is only defined once n > 1.
    f_inter = float(ms_ab / ms_e) if ms_e > 0 else float("inf")
    f_a = float(ms_a / ms_ab) if ms_ab > 0 else float("inf")
    f_b = float(ms_b / ms_ab) if ms_ab > 0 else float("inf")
    df_inter, df_e = (a - 1) * (b - 1), a * b * (n - 1)

    out.update({
        "ms_source_subset": float(ms_a),
        "ms_model_seed": float(ms_b),
        "ms_interaction": float(ms_ab),
        "ms_residual": float(ms_e),
        "df_source_subset": int(a - 1),
        "df_model_seed": int(b - 1),
        "df_interaction": int(df_inter),
        "df_residual": int(df_e),
        "f_interaction_over_residual": f_inter,
        "f_source_over_interaction": f_a,
        "f_seed_over_interaction": f_b,
        "p_interaction": _f_sf(f_inter, df_inter, df_e),
        "p_source_subset": _f_sf(f_a, a - 1, df_inter),
        "p_model_seed": _f_sf(f_b, b - 1, df_inter),
        "total_variance_nonnegative": float(total),
        "grand_mean": float(grand),
        "n_levels_source_subset": int(a),
        "n_levels_model_seed": int(b),
        "n_replicates": int(n),
        "any_raw_component_negative": bool(any(v < 0 for v in raw.values())),
    })
    return out


def replicated_components_vectorised(cube: np.ndarray) -> dict[str, np.ndarray]:
    """Same decomposition applied independently over leading axes.

    Parameters
    ----------
    cube
        Shape (..., a, b, n). Returns arrays of shape (...) for each component.
        Used for the per-cell decomposition, where the leading axis indexes
        target cells and looping in Python would be far too slow.
    """
    arr = np.asarray(cube, dtype=float)
    if arr.ndim < 3:
        raise ValueError(f"expected at least 3 dimensions, got shape {arr.shape}")
    a, b, n = arr.shape[-3:]
    if n < 2:
        raise ValueError("n must be >= 2 to separate interaction from residual")

    cell = arr.mean(axis=-1)                        # (..., a, b)
    grand = cell.mean(axis=(-2, -1))                # (...)
    row = cell.mean(axis=-1)                        # (..., a)
    col = cell.mean(axis=-2)                        # (..., b)
    g = grand[..., None]

    ms_a = b * n * ((row - g) ** 2).sum(axis=-1) / (a - 1)
    ms_b = a * n * ((col - g) ** 2).sum(axis=-1) / (b - 1)
    inter = cell - row[..., :, None] - col[..., None, :] + grand[..., None, None]
    ms_ab = n * (inter ** 2).sum(axis=(-2, -1)) / ((a - 1) * (b - 1))
    ms_e = ((arr - cell[..., None]) ** 2).sum(axis=(-3, -2, -1)) / (a * b * (n - 1))

    raw = {
        "source_subset": (ms_a - ms_ab) / (b * n),
        "model_seed": (ms_b - ms_ab) / (a * n),
        "interaction": (ms_ab - ms_e) / n,
        "residual": ms_e,
    }
    nonneg = {k: np.clip(v, 0.0, None) for k, v in raw.items()}
    total = sum(nonneg.values())
    ok = total > _ZERO_VARIANCE_EPS
    safe = np.where(ok, total, 1.0)

    out: dict[str, np.ndarray] = {}
    for name in COMPONENT_NAMES:
        out[f"{name}_variance_raw"] = raw[name]
        out[f"{name}_variance"] = nonneg[name]
        out[f"{name}_fraction"] = np.where(ok, nonneg[name] / safe, 0.0)
    out["total_variance"] = total
    out["grand_mean"] = grand
    return out
