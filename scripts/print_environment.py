#!/usr/bin/env python3
"""Print pinned analysis environment versions for manuscript reproducibility."""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

OUT = (
    Path(__file__).resolve().parents[1]
    / "results/manuscript_submission_revision/tables"
)
OUT.mkdir(parents=True, exist_ok=True)

PKGS = [
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "anndata",
    "scanpy",
    "torch",
    "scvi",
    "matplotlib",
    "mudata",
]


def ver(name: str) -> str:
    try:
        if name == "sklearn":
            import sklearn

            return sklearn.__version__
        if name == "scvi":
            import scvi

            return scvi.__version__
        mod = __import__(name)
        return getattr(mod, "__version__", "unknown")
    except Exception as exc:  # noqa: BLE001
        return f"UNAVAILABLE ({exc.__class__.__name__})"


def main():
    rows = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {p: ver(p) for p in PKGS},
        "pins_file": "environment.yml / requirements.txt (project pins: Python 3.11.7, scvi-tools 1.3.3, scikit-learn 1.5.2, torch 2.14.0)",
        "note": "Runtime versions on the executing machine may differ from pinned project environment.",
    }
    out = OUT / "runtime_environment.json"
    out.write_text(json.dumps(rows, indent=2))
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
