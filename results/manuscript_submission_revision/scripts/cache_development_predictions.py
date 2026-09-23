#!/usr/bin/env python3
"""Cache development held-out predictions for bootstrap sensitivities (one-time)."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BOOT_MOD = ROOT / "results/presubmission_revision/bootstrap/presubmission_core_contrast_bootstrap.py"
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    spec = importlib.util.spec_from_file_location("core_boot", BOOT_MOD)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    y, p_rna, p_cat, p_tot, audit = mod.load_development_predictions()
    np.savez_compressed(
        OUT / "development_heldout_predictions.npz",
        y=y,
        pred_rna=p_rna,
        pred_concat=p_cat,
        pred_totalvi=p_tot,
    )
    (OUT / "development_heldout_predictions_audit.json").write_text(json.dumps(audit, indent=2, default=str))
    print(json.dumps(audit, indent=2, default=str))


if __name__ == "__main__":
    main()
