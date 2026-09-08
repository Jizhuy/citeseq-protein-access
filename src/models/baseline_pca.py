"""RNA PCA control trained with the MOFA+ RNA view. See run_mofa.py."""

from __future__ import annotations

from src.models.run_mofa import run_mofa


def run_pca(*args, **kwargs):
    """PCA_RNA is written by the MOFA+ runner on the same HVGs."""
    return run_mofa(*args, **kwargs)
