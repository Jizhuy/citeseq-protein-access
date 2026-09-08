#!/usr/bin/env python
"""PHASE 9 API validation only. Does not train the missing-modality grid."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiments.phase9_api_validation import main


if __name__ == "__main__":
    main()
