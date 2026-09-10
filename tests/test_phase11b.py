#!/usr/bin/env python
"""PHASE 11B unit tests.

Covers the pure-numpy layer that the 130-run grid depends on, in particular the
regression guard for the PHASE 10 latent-alignment defect (§5). Runs on CPU and
needs no trained artifacts.

    python tests/test_phase11b.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.experiments.phase11b_stats import (  # noqa: E402
    assert_post_adaptation,
    classify_evidence,
    latent_drift,
    mean_pairwise_distance,
    one_hot_matrix,
    two_factor_components,
    variation_ratio,
    vector_components,
)

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    if condition:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name} {detail}")
        FAILURES.append(name)


# ---------------------------------------------------------------------------
def test_phase10_regression_guard() -> None:
    """The defect: the evaluator was handed the PRE-adaptation source array."""
    print("\n[regression guard for the PHASE 10 latent-alignment defect]")
    rng = np.random.default_rng(0)
    pre = rng.normal(size=(50, 20))
    post = pre + rng.normal(scale=0.02, size=(50, 20))  # totalVI-like non-zero drift

    assert_post_adaptation(post, pre, post, "ok")
    check("post-adaptation array accepted", True)

    # This is exactly what PHASE 10 did and what must now be impossible.
    try:
        assert_post_adaptation(pre, pre, post, "bug")
        check("pre-adaptation array rejected", False, "no error raised")
    except RuntimeError:
        check("pre-adaptation array rejected", True)

    # A copy has equal content but is a different object; still rejected,
    # because only the actual post-adaptation object may reach the evaluator.
    try:
        assert_post_adaptation(post.copy(), pre, post, "copy")
        check("equal-valued copy rejected", False, "no error raised")
    except RuntimeError:
        check("equal-valued copy rejected", True)

    # scVI leaves the reference embedding untouched: zero drift is legitimate.
    assert_post_adaptation(pre, pre, pre, "zero-drift")
    check("zero-drift scVI case accepted", True)

    d = latent_drift(pre, post)
    check("drift reports non-zero displacement", d["mean_cell_displacement"] > 0)
    check("drift flags identical arrays", latent_drift(pre, pre)["identical"] is True)


def test_two_factor_components() -> None:
    """Method-of-moments decomposition on a table with a known structure."""
    print("\n[two-factor variance components]")
    a = np.array([0.0, 1.0, 2.0, 3.0, 4.0])          # strong row (source) effect
    mat = a[:, None] + np.zeros((5, 5))
    comp = two_factor_components(mat)
    check("pure row effect -> seed variance 0",
          abs(comp["model_seed_variance_nonnegative"]) < 1e-12)
    check("pure row effect -> source fraction 1", abs(comp["source_fraction"] - 1.0) < 1e-9)
    check("pure row effect -> zero residual", abs(comp["ms_interaction_residual"]) < 1e-12)

    mat = np.zeros((5, 5)) + a[None, :]               # pure column (seed) effect
    comp = two_factor_components(mat)
    check("pure column effect -> seed fraction 1", abs(comp["seed_fraction"] - 1.0) < 1e-9)

    # Additive model: variance splits between the two main effects, none left over.
    mat = a[:, None] + (2 * a)[None, :]
    comp = two_factor_components(mat)
    check("additive model leaves no residual", abs(comp["ms_interaction_residual"]) < 1e-10)
    check("additive fractions sum to 1",
          abs(comp["source_fraction"] + comp["seed_fraction"]
              + comp["interaction_fraction"] - 1.0) < 1e-9)

    # Pure noise: raw estimates may go negative and must be preserved, not clipped away.
    rng = np.random.default_rng(3)
    raws = [two_factor_components(rng.normal(size=(5, 5))) for _ in range(60)]
    check("negative raw components are retained",
          any(r["source_subset_variance_raw"] < 0 or r["model_seed_variance_raw"] < 0
              for r in raws))
    check("non-negative components never negative",
          all(r["source_subset_variance_nonnegative"] >= 0
              and r["model_seed_variance_nonnegative"] >= 0 for r in raws))

    grand = two_factor_components(np.arange(25, dtype=float).reshape(5, 5))["grand_mean"]
    check("grand mean correct", abs(grand - 12.0) < 1e-12)


def test_vectorised_matches_scalar() -> None:
    """The per-cell vectorised decomposition must equal the scalar one, cell by cell."""
    print("\n[vectorised per-cell decomposition]")
    rng = np.random.default_rng(7)
    cube = rng.random((5, 5, 40))
    vec = vector_components(cube)
    ok_src = ok_seed = True
    for i in range(cube.shape[2]):
        ref = two_factor_components(cube[:, :, i])
        ok_src &= abs(vec["source_var"][i] - ref["source_subset_variance_nonnegative"]) < 1e-12
        ok_seed &= abs(vec["seed_var"][i] - ref["model_seed_variance_nonnegative"]) < 1e-12
    check("vectorised source variance matches scalar", ok_src)
    check("vectorised seed variance matches scalar", ok_seed)
    frac = vec["source_fraction"] + vec["seed_fraction"] + vec["interaction_fraction"]
    check("per-cell fractions sum to 1", np.allclose(frac[np.isfinite(frac)], 1.0))


def test_one_hot() -> None:
    """The one-hot is built once and resampled, so it must survive indexing."""
    print("\n[one-hot precomputation]")
    classes = np.array(["A", "B", "C"], dtype=object)
    y = np.array(["A", "C", "B", "C", "A", "B"], dtype=object)
    oh = one_hot_matrix(y, classes)
    check("one-hot row sums are 1", np.allclose(oh.sum(axis=1), 1.0))
    check("one-hot places labels correctly",
          oh[0, 0] == 1.0 and oh[1, 2] == 1.0 and oh[2, 1] == 1.0)

    # A label outside the class set contributes an all-zero row rather than crashing.
    oh2 = one_hot_matrix(np.array(["A", "Z"], dtype=object), classes)
    check("unknown label yields zero row", oh2[1].sum() == 0.0)

    # Indexing the precomputed matrix must equal rebuilding it on the resample,
    # which is what lets metric_block skip the per-replicate encoding.
    rng = np.random.default_rng(11)
    n = 300
    y = rng.choice(classes, size=n)
    idx = rng.integers(0, n, n)
    check("resampled one-hot equals rebuilt one-hot",
          np.array_equal(one_hot_matrix(y, classes)[idx], one_hot_matrix(y[idx], classes)))


def test_variation_ratio_and_distance() -> None:
    print("\n[instability measures]")
    preds = np.array([["A", "A", "A"], ["A", "B", "A"], ["A", "B", "C"],
                      ["A", "B", "C"], ["A", "B", "C"]], dtype=object)
    vr = variation_ratio(preds)
    check("unanimous cell has variation ratio 0", vr[0] == 0.0)
    check("4/5 modal gives 0.2", abs(vr[1] - 0.2) < 1e-12)
    check("3/5 modal gives 0.4", abs(vr[2] - 0.4) < 1e-12)

    rng = np.random.default_rng(5)
    stack = rng.random((6, 25, 4))
    got = mean_pairwise_distance(stack)
    brute = np.zeros(25)
    for i in range(25):
        ds = [np.linalg.norm(stack[a, i] - stack[b, i])
              for a in range(6) for b in range(a + 1, 6)]
        brute[i] = np.sqrt(np.mean(np.square(ds)))
    check("mean pairwise distance matches brute force", np.allclose(got, brute))


def test_evidence_rule() -> None:
    print("\n[evidence classification rule]")
    check("CI excludes 0 + consistent signs + cell support -> robust",
          classify_evidence((0.01, 0.05), 18, 20, 0.8) == "robust")
    check("CI excludes 0 + consistent signs + weak cell support -> cell_sampling_sensitive",
          classify_evidence((0.01, 0.05), 18, 20, 0.1) == "cell_sampling_sensitive")
    check("sign flips despite strong per-seed CIs -> seed_sensitive",
          classify_evidence((-0.01, 0.05), 11, 20, 0.9) == "seed_sensitive")
    check("nothing separates the models -> inconclusive",
          classify_evidence((-0.02, 0.03), 10, 20, 0.1) == "inconclusive")
    check("consistent negative signs also count",
          classify_evidence((-0.05, -0.01), 2, 20, 0.8) == "robust")
    check("non-finite CI does not crash",
          classify_evidence((float("nan"), float("nan")), 20, 20, 0.9) == "seed_sensitive")


def main() -> int:
    print("PHASE 11B unit tests")
    test_phase10_regression_guard()
    test_two_factor_components()
    test_vectorised_matches_scalar()
    test_one_hot()
    test_variation_ratio_and_distance()
    test_evidence_rule()
    print(f"\n{'ALL TESTS PASSED' if not FAILURES else f'{len(FAILURES)} FAILURES: {FAILURES}'}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
