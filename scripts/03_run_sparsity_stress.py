#!/usr/bin/env python
"""PHASE 8: protein sparsity / corruption stress test.

Use scripts/08_run_sparsity_stress.py. This wrapper keeps the original path.
"""

from __future__ import annotations

import runpy
from pathlib import Path

if __name__ == "__main__":
    target = Path(__file__).with_name("08_run_sparsity_stress.py")
    runpy.run_path(str(target), run_name="__main__")
