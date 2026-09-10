"""Unit tests for the PHASE 12A replicated variance decomposition.

Pure NumPy only, so these run without torch/scvi/anndata installed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase12a_stats import (  # noqa: E402
    replicated_components,
    replicated_components_vectorised,
)


def _simulate(a, b, n, s2_a, s2_b, s2_ab, s2_e, rng):
    A = rng.normal(0, np.sqrt(s2_a), size=(a, 1, 1)) if s2_a else np.zeros((a, 1, 1))
    B = rng.normal(0, np.sqrt(s2_b), size=(1, b, 1)) if s2_b else np.zeros((1, b, 1))
    AB = rng.normal(0, np.sqrt(s2_ab), size=(a, b, 1)) if s2_ab else np.zeros((a, b, 1))
    E = rng.normal(0, np.sqrt(s2_e), size=(a, b, n))
    return 10.0 + A + B + AB + E


def test_rejects_unreplicated_design():
    """n = 1 is exactly the PHASE 11B design; it must not silently return zeros."""
    with pytest.raises(ValueError, match="n must be >= 2"):
        replicated_components(np.zeros((5, 5, 1)))


def test_rejects_bad_shape_and_nonfinite():
    with pytest.raises(ValueError, match="3-D"):
        replicated_components(np.zeros((5, 5)))
    bad = np.zeros((5, 5, 2))
    bad[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        replicated_components(bad)


def test_pure_residual_recovered_on_raw_scale():
    """No factor effects: the raw residual estimate must recover the truth."""
    est, pvals = [], []
    for s in range(200):
        rng = np.random.default_rng(s)
        out = replicated_components(_simulate(5, 5, 2, 0.0, 0.0, 0.0, 1.0, rng))
        est.append(out["residual_variance_raw"])
        pvals.append(out["p_interaction"])
    assert float(np.mean(est)) == pytest.approx(1.0, rel=0.05)
    # The F-test, not the fraction, identifies a null interaction: under H0 the
    # p-values are uniform, so only about 5% should fall below 0.05.
    assert float(np.mean(np.array(pvals) < 0.05)) < 0.12


def test_truncated_fractions_are_upward_biased_for_null_factors():
    """Documents a real caveat: on pure noise the truncated residual fraction
    is far below 1, because null components get clipped at zero from below but
    keep their positive excursions. The reported fractions must therefore never
    be read as 'this factor explains X% of variance' without the F-tests."""
    fracs, praw = [], []
    for s in range(200):
        rng = np.random.default_rng(500 + s)
        out = replicated_components(_simulate(5, 5, 2, 0.0, 0.0, 0.0, 1.0, rng))
        fracs.append(out["residual_fraction"])
        praw.append(out["source_subset_variance_raw"])
    mean_frac = float(np.mean(fracs))
    assert 0.6 < mean_frac < 0.95, (
        f"null residual fraction {mean_frac:.3f}: truncation inflates the three "
        "null components by roughly 15% of the total")
    # ...while the raw estimate of a null component is unbiased around zero.
    assert float(np.mean(praw)) == pytest.approx(0.0, abs=0.15)


def test_f_test_detects_real_interaction():
    rng = np.random.default_rng(11)
    cube = _simulate(5, 5, 2, 0.0, 0.0, 4.0, 0.05, rng)
    out = replicated_components(cube)
    assert out["p_interaction"] < 1e-6
    assert out["f_interaction_over_residual"] > 10


def test_pure_interaction_is_not_called_residual():
    """The whole point of PHASE 12A: real interaction must separate from noise."""
    rng = np.random.default_rng(1)
    cube = _simulate(5, 5, 2, 0.0, 0.0, 4.0, 0.05, rng)
    out = replicated_components(cube)
    assert out["interaction_fraction"] > 0.9
    assert out["residual_fraction"] < 0.1
    assert out["interaction_variance_raw"] == pytest.approx(4.0, rel=0.6)


def test_recovers_known_components_on_average():
    """Method-of-moments estimators are unbiased; average over many draws."""
    truth = {"source_subset": 1.0, "model_seed": 0.5, "interaction": 0.25, "residual": 2.0}
    acc = {k: [] for k in truth}
    for s in range(400):
        rng = np.random.default_rng(1000 + s)
        cube = _simulate(5, 5, 2, truth["source_subset"], truth["model_seed"],
                         truth["interaction"], truth["residual"], rng)
        out = replicated_components(cube)
        for k in truth:
            acc[k].append(out[f"{k}_variance_raw"])
    for k, want in truth.items():
        got = float(np.mean(acc[k]))
        assert got == pytest.approx(want, abs=0.25), f"{k}: got {got}, want {want}"


def test_degrees_of_freedom():
    out = replicated_components(np.random.default_rng(2).normal(size=(5, 5, 2)))
    assert out["df_source_subset"] == 4
    assert out["df_model_seed"] == 4
    assert out["df_interaction"] == 16
    assert out["df_residual"] == 25
    assert out["n_replicates"] == 2


def test_negative_raw_is_preserved_and_flagged():
    """Truncation must not erase the evidence that an estimate went negative."""
    rng = np.random.default_rng(7)
    cube = _simulate(5, 5, 2, 0.0, 0.0, 0.0, 1.0, rng)
    out = replicated_components(cube)
    if out["any_raw_component_negative"]:
        negatives = [k for k in ("source_subset", "model_seed", "interaction", "residual")
                     if out[f"{k}_variance_raw"] < 0]
        for k in negatives:
            assert out[f"{k}_variance_nonnegative"] == 0.0
            assert out[f"{k}_variance_raw"] < 0.0


def test_vectorised_matches_scalar():
    rng = np.random.default_rng(3)
    cubes = np.stack([_simulate(5, 5, 2, 1.0, 0.5, 0.25, 2.0, rng) for _ in range(12)])
    vec = replicated_components_vectorised(cubes)
    for i in range(cubes.shape[0]):
        one = replicated_components(cubes[i])
        for name in ("source_subset", "model_seed", "interaction", "residual"):
            assert vec[f"{name}_variance_raw"][i] == pytest.approx(
                one[f"{name}_variance_raw"], rel=1e-10, abs=1e-12)
            assert vec[f"{name}_fraction"][i] == pytest.approx(
                one[f"{name}_fraction"], rel=1e-10, abs=1e-12)


def test_fractions_sum_to_one():
    rng = np.random.default_rng(4)
    cube = _simulate(5, 5, 2, 1.0, 0.5, 0.25, 2.0, rng)
    out = replicated_components(cube)
    total = sum(out[f"{k}_fraction"] for k in
                ("source_subset", "model_seed", "interaction", "residual"))
    assert total == pytest.approx(1.0)


def test_constant_cube_gives_zero_variance():
    out = replicated_components(np.full((5, 5, 2), 3.7))
    assert out["total_variance_nonnegative"] == pytest.approx(0.0, abs=1e-20)
    assert out["grand_mean"] == pytest.approx(3.7)
    for k in ("source_subset", "model_seed", "interaction", "residual"):
        assert out[f"{k}_fraction"] == 0.0
