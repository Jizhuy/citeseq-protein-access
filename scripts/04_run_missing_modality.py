#!/usr/bin/env python
"""PHASE 9: true missing-protein-modality stress test.

Training is blocked until a supported cell-level missingness encoding is approved.
Run the API validation instead:

    python scripts/09_validate_missing_modality_api.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    print(
        "PHASE 9 was stopped after API validation. "
        "scvi-tools 1.3.3 totalVI does not natively support cell-level missing "
        "protein profiles. See results/logs/phase9_api_validation.md and "
        "results/logs/phase9_stop.md. The missing-modality training grid was not launched."
    )
    sys.exit(2)


if __name__ == "__main__":
    main()
