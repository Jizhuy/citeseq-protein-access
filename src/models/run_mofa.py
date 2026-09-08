"""PHASE 7: MOFA+ on processed RNA + CLR protein. Not a totalVI architecture match."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from mofapy2.run.entry_point import entry_point
from sklearn.decomposition import PCA

from src.evaluation.plots import plot_umap_by_batch
from src.evaluation.representation_metrics import batch_silhouette, latent_diagnostics
from src.evaluation.resource_metrics import Timer, peak_rss_mb
from src.models.run_scvi import EXPECTED_BATCHES, EXPECTED_N_OBS, EXPECTED_N_VARS, ScviDataError
from src.utils.io import load_config, resolve_path, save_json
from src.utils.logging_utils import get_logger
from src.utils.seed import set_global_seed

N_HVG = 2000
N_FACTORS = 20
GROUPS = ("PBMC10k", "PBMC5k")


def clr_rows(matrix: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    """Centered log-ratio per cell. Standard CITE-seq protein transform for MOFA+."""
    values = np.asarray(matrix, dtype=np.float64)
    if np.isnan(values).any():
        raise ValueError("Protein matrix contains NaN; inner-join proteins must be observed.")
    values = np.clip(values, 0.0, None) + float(pseudocount)
    logx = np.log(values)
    return logx - logx.mean(axis=1, keepdims=True)


def _decode_h5_strings(values) -> list[str]:
    out = []
    for value in np.asarray(values):
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8"))
        else:
            out.append(str(value))
    return out


def factors_from_hdf5(outfile: str | Path, original_order: pd.Index) -> np.ndarray:
    """Read group-wise Z from a mofapy2 HDF5 file and restore original cell order.

    Saved arrays are (n_factors, n_cells_in_group).
    """
    frames = []
    with h5py.File(outfile, "r") as handle:
        for group in GROUPS:
            z = np.asarray(handle[f"expectations/Z/{group}"][:], dtype=np.float64).T
            names = _decode_h5_strings(handle[f"samples/{group}"][:])
            if z.shape != (len(names), N_FACTORS):
                raise RuntimeError(
                    f"MOFA+ HDF5 Z/{group} has shape {z.shape}, expected ({len(names)}, {N_FACTORS})"
                )
            frames.append(pd.DataFrame(z, index=names))
    latent_df = pd.concat(frames).reindex(original_order)
    if latent_df.isna().any().any():
        raise RuntimeError("MOFA+ factors could not be aligned to the original cell order.")
    latent = latent_df.to_numpy(dtype=np.float64)
    if latent.shape != (len(original_order), N_FACTORS):
        raise RuntimeError(f"MOFA+ latent shape {latent.shape}")
    return latent


def variance_explained_from_hdf5(outfile: str | Path) -> pd.DataFrame:
    rows = []
    with h5py.File(outfile, "r") as handle:
        view_names = _decode_h5_strings(handle["views/views"][:])
        for group in GROUPS:
            total = np.asarray(handle[f"variance_explained/r2_total/{group}"][:], dtype=float)
            for vi, view in enumerate(view_names[: len(total)]):
                rows.append(
                    {
                        "group": group,
                        "view": view,
                        "factor": "total",
                        "variance_explained": float(total[vi]),
                    }
                )
            per = np.asarray(
                handle[f"variance_explained/r2_per_factor/{group}"][:], dtype=float
            )
            for vi, view in enumerate(view_names[: per.shape[0]]):
                for k in range(per.shape[1]):
                    rows.append(
                        {
                            "group": group,
                            "view": view,
                            "factor": k + 1,
                            "variance_explained": float(per[vi, k]),
                        }
                    )
    return pd.DataFrame(rows)


def hvgs_from_hdf5(outfile: str | Path) -> list[str]:
    with h5py.File(outfile, "r") as handle:
        return _decode_h5_strings(handle["features/rna"][:])


def prepare_mofa_matrices(
    adata: ad.AnnData,
    n_hvg: int = N_HVG,
    seed: int = 0,
    hvgs: list[str] | None = None,
) -> dict[str, Any]:
    """MOFA+-appropriate views. Does not feed raw counts into the factor model."""
    rna = ad.AnnData(
        X=adata.layers["counts"].copy() if "counts" in adata.layers else adata.X.copy()
    )
    rna.obs_names = adata.obs_names
    rna.var_names = adata.var_names
    if hvgs is None:
        sc.pp.highly_variable_genes(rna, n_top_genes=n_hvg, flavor="seurat_v3")
        hvgs = rna.var_names[rna.var["highly_variable"]].astype(str).tolist()
    missing = [g for g in hvgs if g not in set(rna.var_names.astype(str))]
    if missing:
        raise ValueError(f"{len(missing)} saved HVGs are absent from the query RNA.")
    rna = rna[:, hvgs].copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.scale(rna, max_value=10)
    rna_mat = np.asarray(rna.X, dtype=np.float64)

    protein = adata.obsm["protein_expression"]
    if not isinstance(protein, pd.DataFrame):
        raise TypeError("protein_expression must be a DataFrame of observed counts.")
    protein_names = [str(c) for c in protein.columns]
    prot = clr_rows(protein.to_numpy())
    prot = (prot - prot.mean(axis=0)) / (prot.std(axis=0, ddof=1) + 1e-8)

    if rna_mat.shape != (adata.n_obs, len(hvgs)):
        raise ValueError(f"RNA view shape {rna_mat.shape} != ({adata.n_obs}, {len(hvgs)})")
    if prot.shape != (adata.n_obs, 14):
        raise ValueError(f"Protein view shape {prot.shape} != ({adata.n_obs}, 14)")
    if not np.isfinite(rna_mat).all() or not np.isfinite(prot).all():
        raise ValueError("MOFA+ views contain non-finite values.")
    return {
        "rna": rna_mat,
        "protein": prot,
        "hvgs": hvgs,
        "protein_names": protein_names,
        "n_hvg": int(len(hvgs)),
    }


def _split_by_group(matrix: np.ndarray, batch: np.ndarray, obs_names: pd.Index):
    data, names = [], []
    for group in GROUPS:
        mask = batch == group
        data.append(np.asarray(matrix[mask], dtype=np.float64))
        names.append(obs_names[mask].astype(str).tolist())
    return data, names


def run_mofa(
    config: dict[str, Any] | None = None,
    seed: int = 0,
    reuse_existing: bool = True,
) -> dict[str, Any]:
    """Train MOFA+ and a same-HVG RNA PCA control. Does not touch VAE models."""
    cfg = config or load_config()
    logs_dir = resolve_path(cfg["paths"]["logs_dir"])
    tables_dir = resolve_path(cfg["paths"]["tables_dir"])
    figures_dir = resolve_path(cfg["paths"]["figures_dir"])
    embeddings_dir = resolve_path(cfg["paths"]["embeddings_dir"])
    models_dir = resolve_path(cfg["paths"].get("models_dir", "results/models"))
    processed_dir = resolve_path(cfg["paths"]["processed_dir"])
    logger = get_logger("run_mofa", logs_dir / f"mofa_seed{seed}.log")

    set_global_seed(seed)
    source_path = processed_dir / "pbmc_cite_combined_inner.h5ad"
    adata = ad.read_h5ad(source_path)
    if adata.n_obs != EXPECTED_N_OBS or adata.n_vars != EXPECTED_N_VARS:
        raise ScviDataError(f"Unexpected inner-join shape {adata.shape}")
    batches = set(adata.obs["batch"].astype(str))
    if batches != EXPECTED_BATCHES:
        raise ScviDataError(f"Unexpected batches {batches}")
    original_order = adata.obs_names.astype(str).copy()
    original_split = adata.obs["split"].astype(str).copy()

    model_dir = models_dir / f"mofa_seed{seed}"
    model_dir.mkdir(parents=True, exist_ok=True)
    outfile = str(model_dir / "mofa_model.hdf5")
    hdf5_ready = Path(outfile).exists() and Path(outfile).stat().st_size > 0
    reused_hdf5 = bool(reuse_existing and hdf5_ready)
    saved_hvgs = hvgs_from_hdf5(outfile) if reused_hdf5 else None

    views = prepare_mofa_matrices(adata, n_hvg=N_HVG, seed=seed, hvgs=saved_hvgs)
    batch = adata.obs["batch"].astype(str).to_numpy()
    rna_groups, sample_names = _split_by_group(views["rna"], batch, adata.obs_names)
    prot_groups, _ = _split_by_group(views["protein"], batch, adata.obs_names)

    logger.info(
        "Training MOFA+: %s cells, RNA HVGs=%s, proteins=%s, factors=%s, groups=%s, seed=%s, reuse_hdf5=%s",
        adata.n_obs,
        views["n_hvg"],
        views["protein"].shape[1],
        N_FACTORS,
        list(GROUPS),
        seed,
        reused_hdf5,
    )
    timer = Timer()
    timer.start()
    if reused_hdf5:
        logger.info("Reusing trained MOFA+ HDF5 at %s; not refitting.", outfile)
        with h5py.File(outfile, "r") as handle:
            train_time = np.asarray(handle["training_stats/time"][:], dtype=float)
            runtime = float(np.nansum(train_time))
            n_iter_saved = int(np.isfinite(train_time).sum())
    else:
        ent = entry_point()
        ent.set_data_options(scale_views=True, scale_groups=False, center_groups=True)
        ent.set_data_matrix(
            [rna_groups, prot_groups],
            likelihoods=["gaussian", "gaussian"],
            views_names=["rna", "protein"],
            groups_names=list(GROUPS),
            samples_names=sample_names,
            features_names=[views["hvgs"], views["protein_names"]],
        )
        ent.set_model_options(
            factors=N_FACTORS,
            spikeslab_weights=True,
            ard_weights=True,
            ard_factors=True,
            spikeslab_factors=False,
        )
        ent.set_train_options(
            iter=1000,
            convergence_mode="medium",
            seed=int(seed),
            gpu_mode=False,
            weight_views=True,
            quiet=False,
            outfile=outfile,
        )
        ent.build()
        ent.run()
        ent.save(outfile=outfile, save_data=False)
        runtime = timer.stop()
        n_iter_saved = 1000

    latent = factors_from_hdf5(outfile, original_order)
    if latent.shape != (EXPECTED_N_OBS, N_FACTORS):
        raise RuntimeError(f"MOFA+ latent shape {latent.shape} != ({EXPECTED_N_OBS}, {N_FACTORS})")

    pca = PCA(n_components=N_FACTORS, random_state=seed)
    pca_latent = pca.fit_transform(views["rna"])
    if pca_latent.shape != (EXPECTED_N_OBS, N_FACTORS):
        raise RuntimeError(f"PCA latent shape {pca_latent.shape}")

    embeddings_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    r2_table = variance_explained_from_hdf5(outfile)
    r2_table.to_csv(tables_dir / f"mofa_seed{seed}_variance_explained.csv", index=False)
    r2_total = r2_table[r2_table["factor"].astype(str) == "total"]

    np.save(embeddings_dir / f"mofa_seed{seed}_latent.npy", latent)
    np.save(embeddings_dir / f"pca_rna_seed{seed}_latent.npy", pca_latent)
    cells = pd.DataFrame({"cell_id": original_order, "batch": batch, "split": original_split.to_numpy()})
    cells.to_csv(embeddings_dir / f"mofa_seed{seed}_cells.csv", index=False)
    cells.to_csv(embeddings_dir / f"pca_rna_seed{seed}_cells.csv", index=False)
    pd.Series(views["hvgs"], name="gene").to_csv(model_dir / "mofa_rna_hvgs.csv", index=False)

    adata_umap = ad.AnnData(X=latent.copy())
    sc.pp.neighbors(adata_umap, use_rep="X", n_neighbors=15, random_state=seed)
    sc.tl.umap(adata_umap, random_state=seed)
    umap = np.asarray(adata_umap.obsm["X_umap"])
    np.save(embeddings_dir / f"mofa_seed{seed}_umap.npy", umap)
    plot_umap_by_batch(
        umap,
        pd.Series(batch, index=original_order),
        figures_dir / f"mofa_seed{seed}_umap_batch",
        title="MOFA+ factor UMAP colored by batch",
    )

    pca_umap_ad = ad.AnnData(X=pca_latent.copy())
    sc.pp.neighbors(pca_umap_ad, use_rep="X", n_neighbors=15, random_state=seed)
    sc.tl.umap(pca_umap_ad, random_state=seed)
    pca_umap = np.asarray(pca_umap_ad.obsm["X_umap"])
    np.save(embeddings_dir / f"pca_rna_seed{seed}_umap.npy", pca_umap)
    plot_umap_by_batch(
        pca_umap,
        pd.Series(batch, index=original_order),
        figures_dir / f"pca_rna_seed{seed}_umap_batch",
        title="RNA PCA UMAP colored by batch",
    )

    diag_mofa = latent_diagnostics(latent)
    diag_pca = latent_diagnostics(pca_latent)
    batch_asw_mofa = float(batch_silhouette(latent, batch))
    batch_asw_pca = float(batch_silhouette(pca_latent, batch))

    info = {
        "model": "MOFA+",
        "package": "mofapy2",
        "package_version": __import__("mofapy2").__version__,
        "seed": int(seed),
        "n_cells": int(adata.n_obs),
        "n_factors": N_FACTORS,
        "n_hvg": views["n_hvg"],
        "n_proteins": int(views["protein"].shape[1]),
        "rna_transform": "seurat_v3 HVG; library-size 1e4; log1p; scale max=10",
        "protein_transform": "CLR with pseudocount 1, then per-protein z-score",
        "likelihoods": ["gaussian", "gaussian"],
        "groups": list(GROUPS),
        "scale_views": True,
        "weight_views": True,
        "ard_factors": True,
        "spikeslab_weights": True,
        "convergence_mode": "medium",
        "n_iterations": int(n_iter_saved),
        "formal_convergence": False,
        "reused_existing_hdf5": bool(reused_hdf5),
        "uses_raw_counts": False,
        "initialized_from_vae_models": False,
        "matched_to_totalVI_architecture": False,
        "runtime_seconds": float(runtime),
        "peak_rss_mb": float(peak_rss_mb()),
        "batch_silhouette": batch_asw_mofa,
        "pca_rna_batch_silhouette": batch_asw_pca,
        "pca_rna_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "variance_explained_total": r2_total.to_dict(orient="records"),
        "latent_diagnostics": diag_mofa,
        "pca_latent_diagnostics": diag_pca,
        "model_path": outfile,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "note": (
            "MOFA+ is a linear Bayesian factor model on processed views. "
            "It is not a likelihood-matched control for totalVI. "
            "weight_views=True reduces domination of the 2000-gene RNA view over 14 proteins. "
            "Training reached max_iter=1000 without formal medium-mode convergence "
            "because residual ELBO drift stayed above the 0.00005% threshold."
        ),
    }
    save_json(info, logs_dir / f"mofa_seed{seed}_manifest.json")
    save_json(info, model_dir / "mofa_model_info.json")
    logger.info(
        "PHASE 7 MOFA+ complete: runtime=%.1fs batch_silhouette=%.4f PCA_batch_silhouette=%.4f",
        runtime,
        batch_asw_mofa,
        batch_asw_pca,
    )
    if list(adata.obs_names.astype(str)) != list(original_order):
        raise RuntimeError("Query cell order changed.")
    if not adata.obs["split"].astype(str).equals(original_split):
        raise RuntimeError("PHASE 2 split was altered.")
    reloaded = ad.read_h5ad(source_path)
    if "X_MOFA" in reloaded.obsm or "cell_type_l1" in reloaded.obs:
        raise RuntimeError("Original inner h5ad must remain unmodified.")
    return info
