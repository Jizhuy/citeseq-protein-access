"""Protein perturbation helpers.

PHASE 8 uses nonzero-entry masking (corruption / sparsification).
True missing-modality hiding is implemented here but must not be used in PHASE 8.
These functions copy arrays. They never overwrite the caller's ground-truth matrix.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


CORRUPTION_LEVELS = (0.0, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85)
PERTURBATION_SEEDS = (0, 1, 2, 3, 4)
MODEL_SEED = 0
SMOKE_CONDITIONS = ((0.0, 0), (0.40, 0), (0.85, 0))


@dataclass
class PerturbationRecord:
    perturbation_type: str
    level: float
    seed: int
    n_cells: int
    n_proteins: int
    n_entries_changed: int


def corruption_code(fraction: float) -> str:
    return f"{int(round(float(fraction) * 100)):03d}"


def mask_path(perturbations_dir: Path, fraction: float, seed: int) -> Path:
    return perturbations_dir / f"protein_corruption_{corruption_code(fraction)}_seed{seed}.npz"


def _protein_array(protein: pd.DataFrame | np.ndarray) -> tuple[np.ndarray, pd.Index, pd.Index]:
    if isinstance(protein, pd.DataFrame):
        return protein.to_numpy(dtype=np.float64, copy=True), protein.index, protein.columns
    values = np.asarray(protein, dtype=np.float64)
    index = pd.RangeIndex(values.shape[0])
    columns = pd.Index([f"p{i}" for i in range(values.shape[1])])
    return values, index, columns


def sparsity_stats(values: np.ndarray, protein_names: list[str]) -> dict[str, Any]:
    values = np.asarray(values, dtype=np.float64)
    zeros = values <= 0
    detected = values > 0
    per_protein = {
        name: {
            "n_detected": int(detected[:, j].sum()),
            "detection_fraction": float(detected[:, j].mean()),
            "sparsity": float(zeros[:, j].mean()),
            "mean_count": float(values[:, j].mean()),
            "median_count": float(np.median(values[:, j])),
        }
        for j, name in enumerate(protein_names)
    }
    return {
        "n_cells": int(values.shape[0]),
        "n_proteins": int(values.shape[1]),
        "n_entries": int(values.size),
        "n_nonzero": int(detected.sum()),
        "n_zero": int(zeros.sum()),
        "overall_sparsity": float(zeros.mean()),
        "mean_detected_proteins_per_cell": float(detected.sum(axis=1).mean()),
        "median_detected_proteins_per_cell": float(np.median(detected.sum(axis=1))),
        "per_protein": per_protein,
        "detected_proteins_per_cell": detected.sum(axis=1).astype(int).tolist(),
    }


def create_nonzero_mask(
    original: np.ndarray,
    fraction: float,
    seed: int,
) -> np.ndarray:
    """Boolean mask of originally nonzero entries selected for zeroing."""
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be in [0, 1]")
    values = np.asarray(original, dtype=np.float64)
    mask = np.zeros(values.shape, dtype=bool)
    nonzero = np.flatnonzero(values > 0)
    n_mask = int(np.floor(fraction * nonzero.size))
    if n_mask:
        rng = np.random.default_rng(int(seed))
        chosen = rng.choice(nonzero, size=n_mask, replace=False)
        mask.ravel()[chosen] = True
    if np.any(mask & ~(values > 0)):
        raise RuntimeError("Corruption mask selected an original zero entry.")
    return mask


def apply_mask(original: np.ndarray, mask: np.ndarray) -> np.ndarray:
    corrupted = np.asarray(original, dtype=np.float64).copy()
    corrupted[np.asarray(mask, dtype=bool)] = 0.0
    return corrupted


def save_corruption_npz(
    path: Path,
    *,
    mask: np.ndarray,
    original: np.ndarray,
    cell_ids: np.ndarray,
    protein_names: list[str],
    fraction: float,
    seed: int,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    original_nonzero = np.asarray(original, dtype=np.float64) > 0
    realized = float(mask.sum() / original_nonzero.sum()) if original_nonzero.any() else 0.0
    np.savez_compressed(
        path,
        mask=np.asarray(mask, dtype=bool),
        original_nonzero=original_nonzero,
        cell_ids=np.asarray(cell_ids, dtype=object),
        protein_names=np.asarray(protein_names, dtype=object),
        fraction=np.array(float(fraction)),
        seed=np.array(int(seed)),
        n_masked=np.array(int(mask.sum())),
        realized_nonzero_masked_fraction=np.array(realized),
        experiment="protein_sparsity_corruption",
        not_missing_modality=np.array(True),
    )
    return path


def load_corruption_npz(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=True) as handle:
        return {
            "mask": np.asarray(handle["mask"], dtype=bool),
            "original_nonzero": np.asarray(handle["original_nonzero"], dtype=bool),
            "cell_ids": np.asarray(handle["cell_ids"]).astype(str),
            "protein_names": [str(x) for x in np.asarray(handle["protein_names"])],
            "fraction": float(handle["fraction"]),
            "seed": int(handle["seed"]),
            "n_masked": int(handle["n_masked"]),
            "realized_nonzero_masked_fraction": float(handle["realized_nonzero_masked_fraction"]),
        }


def create_or_load_mask(
    original: np.ndarray,
    fraction: float,
    seed: int,
    path: Path,
    cell_ids,
    protein_names: list[str],
) -> dict[str, Any]:
    """Load a saved mask if present; otherwise create and persist it.

    Downstream evaluation must use this saved mask, never a regenerated one.
    """
    original = np.asarray(original, dtype=np.float64)
    if path.exists():
        payload = load_corruption_npz(path)
        if payload["mask"].shape != original.shape:
            raise ValueError(f"Saved mask shape {payload['mask'].shape} != {original.shape}")
        if abs(payload["fraction"] - float(fraction)) > 1e-12:
            raise ValueError("Saved mask fraction does not match the requested condition.")
        if payload["seed"] != int(seed):
            raise ValueError("Saved mask seed does not match the requested condition.")
        if list(payload["cell_ids"]) != list(np.asarray(cell_ids).astype(str)):
            raise ValueError("Saved mask cell order does not match the query object.")
        return payload
    mask = create_nonzero_mask(original, fraction, seed)
    save_corruption_npz(
        path,
        mask=mask,
        original=original,
        cell_ids=np.asarray(cell_ids).astype(str),
        protein_names=protein_names,
        fraction=fraction,
        seed=seed,
    )
    return load_corruption_npz(path)


def summarize_corruption(
    original: np.ndarray,
    corrupted: np.ndarray,
    mask: np.ndarray,
    protein_names: list[str],
    fraction: float,
    seed: int,
) -> dict[str, Any]:
    original = np.asarray(original, dtype=np.float64)
    corrupted = np.asarray(corrupted, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    orig_stats = sparsity_stats(original, protein_names)
    new_stats = sparsity_stats(corrupted, protein_names)
    n_nonzero = int((original > 0).sum())
    n_masked = int(mask.sum())
    realized = (n_masked / n_nonzero) if n_nonzero else 0.0
    if np.any(original[~mask] != corrupted[~mask]):
        raise RuntimeError("Unmasked protein entries were altered.")
    if np.any(corrupted[mask] != 0):
        raise RuntimeError("Masked protein entries were not set to zero.")
    return {
        "perturbation_type": "protein_nonzero_mask_to_zero",
        "experiment_name": "protein_sparsity_stress_test",
        "not_missing_modality": True,
        "target_corruption_fraction": float(fraction),
        "perturbation_seed": int(seed),
        "n_originally_nonzero": n_nonzero,
        "n_masked": n_masked,
        "realized_nonzero_masked_fraction": float(realized),
        "original_protein_sparsity": orig_stats["overall_sparsity"],
        "final_protein_sparsity": new_stats["overall_sparsity"],
        "original_mean_detected_proteins_per_cell": orig_stats["mean_detected_proteins_per_cell"],
        "final_mean_detected_proteins_per_cell": new_stats["mean_detected_proteins_per_cell"],
        "original_per_protein": orig_stats["per_protein"],
        "final_per_protein": new_stats["per_protein"],
        "original_detected_proteins_per_cell": orig_stats["detected_proteins_per_cell"],
        "final_detected_proteins_per_cell": new_stats["detected_proteins_per_cell"],
    }


def mask_nonzero_protein_counts(
    adata: ad.AnnData,
    fraction: float,
    seed: int,
    protein_key: str = "protein_counts",
) -> tuple[ad.AnnData, PerturbationRecord]:
    """Set a random subset of nonzero protein counts to zero.

    The original object is copied. Ground-truth protein values remain in
    ``obsm['protein_counts_truth']``.
    """
    out = adata.copy()
    protein = out.obsm[protein_key].copy()
    out.obsm["protein_counts_truth"] = protein.copy()
    values, index, columns = _protein_array(protein)
    mask = create_nonzero_mask(values, fraction, seed)
    corrupted = apply_mask(values, mask)
    perturbed = pd.DataFrame(corrupted, index=index, columns=columns)
    out.obsm[protein_key] = perturbed
    if "protein_expression" in out.obsm:
        out.obsm["protein_expression"] = perturbed.copy()
    record = PerturbationRecord(
        perturbation_type="protein_nonzero_mask",
        level=float(fraction),
        seed=int(seed),
        n_cells=out.n_obs,
        n_proteins=perturbed.shape[1],
        n_entries_changed=int(mask.sum()),
    )
    return out, record


def hide_protein_modality(
    adata: ad.AnnData,
    fraction: float,
    seed: int,
    protein_key: str = "protein_counts",
) -> tuple[ad.AnnData, PerturbationRecord]:
    """Hide entire protein profiles for a random subset of cells.

    PHASE 9 helper only. Sets hidden profiles to NaN and records a cell-level
    observedness mask. This is a data-bookkeeping representation.

    It is NOT a supported totalVI 1.3.3 training encoding: NaN protein values
    produce invalid encoder latents, and ordinary zeros are treated as observed
    counts unless an entire batch/panel is all zero. Do not pass the returned
    AnnData to TOTALVI.setup_anndata without an approved missingness mapping.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be in [0, 1]")
    out = adata.copy()
    protein = out.obsm[protein_key].copy()
    out.obsm["protein_counts_truth"] = protein.copy()
    rng = np.random.default_rng(seed)
    n_hide = int(np.floor(fraction * out.n_obs))
    hide_idx = rng.choice(out.n_obs, size=n_hide, replace=False) if n_hide else np.array([], dtype=int)
    working = protein.copy()
    if n_hide:
        working.iloc[hide_idx, :] = np.nan
    mask = out.obsm.get("protein_observed_mask")
    if mask is None:
        mask = pd.DataFrame(True, index=protein.index, columns=protein.columns)
    else:
        mask = mask.copy()
    if n_hide:
        mask.iloc[hide_idx, :] = False
    out.obsm[protein_key] = working
    out.obsm["protein_expression"] = working.copy()
    out.obsm["protein_observed_mask"] = mask
    out.obs["protein_modality_hidden"] = False
    if n_hide:
        out.obs.iloc[hide_idx, out.obs.columns.get_loc("protein_modality_hidden")] = True
    record = PerturbationRecord(
        perturbation_type="missing_protein_modality",
        level=float(fraction),
        seed=int(seed),
        n_cells=out.n_obs,
        n_proteins=working.shape[1],
        n_entries_changed=int(n_hide * working.shape[1]),
    )
    return out, record


def record_to_dict(record: PerturbationRecord) -> dict[str, Any]:
    return asdict(record)
