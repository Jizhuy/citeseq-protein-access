#!/usr/bin/env python
"""How much interaction can the PHASE 12A design actually detect?

The replicated design is 5 source subsets x 5 model seeds x 2 replicates, so
the interaction F-test has 16 and 25 degrees of freedom. That is enough to
identify the interaction at all -- which PHASE 11B could not -- but it is not
a large design, and a non-significant result must not be reported as "there is
no interaction" if the design could not have detected one.

This simulation calibrates exactly that, under the null and under a range of
true interaction-to-residual ratios, and writes the resulting power curve so
the interpretation of the real F-tests is bounded by measured power rather
than by assumption. It uses no experimental data.

Usage
  python scripts/12a_power_check.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase12a_replication import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    REPLICATES,
    ensure_dirs,
)
from src.experiments.phase12a_stats import replicated_components  # noqa: E402

RATIOS = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0)
N_SIMS = 2000
ALPHA = 0.05
SIM_SEED = 12_000


def main() -> int:
    paths = ensure_dirs()
    a, b, n = len(CROSSED_SUBSET_SEEDS), len(CROSSED_MODEL_SEEDS), len(REPLICATES)
    rows = []
    for ratio in RATIOS:
        hit, rfrac, ifrac, raw = [], [], [], []
        for s in range(N_SIMS):
            rng = np.random.default_rng(SIM_SEED + int(ratio * 1000) * N_SIMS + s)
            ab = (rng.normal(0.0, np.sqrt(ratio), size=(a, b, 1))
                  if ratio > 0 else np.zeros((a, b, 1)))
            e = rng.normal(0.0, 1.0, size=(a, b, n))
            out = replicated_components(1.0 + ab + e)
            hit.append(out["p_interaction"] < ALPHA)
            rfrac.append(out["residual_fraction"])
            ifrac.append(out["interaction_fraction"])
            raw.append(out["interaction_variance_raw"])
        rows.append({
            "true_interaction_over_residual": ratio,
            "rejection_rate": float(np.mean(hit)),
            "interpretation": "type I error" if ratio == 0 else "power",
            "median_residual_fraction": float(np.median(rfrac)),
            "median_interaction_fraction": float(np.median(ifrac)),
            "mean_interaction_variance_raw": float(np.mean(raw)),
            "n_sims": N_SIMS,
        })
        print(f"  ratio {ratio:>5.3f}: rejection {rows[-1]['rejection_rate']:.3f}  "
              f"median resid frac {rows[-1]['median_residual_fraction']:.3f}", flush=True)

    frame = pd.DataFrame(rows)
    out = paths["tables"] / "interaction_power_check.csv"
    frame.to_csv(out, index=False)
    print(f"wrote {out}")

    null = frame[frame.true_interaction_over_residual == 0].rejection_rate.iloc[0]
    print(f"\ntype I error {null:.3f} (nominal {ALPHA})")
    detectable = frame[(frame.rejection_rate >= 0.80)
                       & (frame.true_interaction_over_residual > 0)]
    if len(detectable):
        r = detectable.true_interaction_over_residual.min()
        print(f"smallest ratio detected with >=80% power: {r}")
    print("A non-significant interaction F-test therefore bounds the interaction "
          "below roughly this multiple of the residual; it does not prove zero.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
