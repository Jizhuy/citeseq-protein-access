"""Global seeding for Python, NumPy, and PyTorch."""

from __future__ import annotations

import os
import random
from typing import Any


def set_global_seed(seed: int) -> dict[str, Any]:
    """Seed Python, NumPy, and Torch if they are importable.

    Returns a record of which backends were actually seeded so later
    experiment rows can document the RNG state provenance.
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)

    seeded = {"python": True, "numpy": False, "torch": False, "scvi": False, "seed": int(seed)}

    try:
        import numpy as np

        np.random.seed(seed)
        seeded["numpy"] = True
    except ImportError:
        pass

    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        seeded["torch"] = True
    except ImportError:
        pass

    try:
        import scvi

        scvi.settings.seed = seed
        seeded["scvi"] = True
    except ImportError:
        pass

    return seeded
