"""Load official 10x PBMC CITE-seq files used by scvi-tools.

This module intentionally does **not** call `scvi.data.pbmcs_10x_cite_seq`
as the primary loader. That helper:

1. silently intersects RNA features to genes shared by PBMC10k and PBMC5k
2. optionally outer-joins proteins and then fills missing values with 0

Those behaviours mix true zeros with proteins that were never measured.
We download the same official files, inspect each dataset, and construct
inner/outer combined objects with an explicit protein observed-mask.
"""

from __future__ import annotations

import re
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from src.utils.io import load_yaml, project_root, resolve_path
from src.utils.logging_utils import get_logger

logger = get_logger(__name__)

CANDIDATE_CELLTYPE_COLUMNS = (
    "cell_type",
    "celltype",
    "celltype.l1",
    "celltype.l2",
    "cell_type_l1",
    "cell_type_l2",
    "predicted.celltype.l2",
    "leiden",
    "cluster",
    "clusters",
    "annotation",
)


def download_file(url: str, dest: Path, fallback_url: str | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Using cached file: %s", dest)
        return dest

    urls = [url] + ([fallback_url] if fallback_url else [])
    last_error: Exception | None = None
    for candidate in urls:
        logger.info("Downloading %s -> %s", candidate, dest)
        try:
            request = urllib.request.Request(
                candidate,
                headers={"User-Agent": "multiomics_robustness/0.1 (research data download)"},
            )
            with urllib.request.urlopen(request) as response, dest.open("wb") as handle:
                handle.write(response.read())
            if dest.exists() and dest.stat().st_size > 0:
                return dest
        except (urllib.error.URLError, OSError) as exc:
            last_error = exc
            logger.warning("Download failed for %s: %s", candidate, exc)
            if dest.exists():
                dest.unlink()
    raise RuntimeError(f"Could not download {dest.name} from {urls}") from last_error


def _as_dense_array(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray())
    return np.asarray(matrix)


def _protein_frame_from_adata(adata: ad.AnnData) -> pd.DataFrame:
    if "protein_expression" not in adata.obsm:
        raise KeyError(
            "Expected protein counts in adata.obsm['protein_expression']. "
            f"Available obsm keys: {list(adata.obsm.keys())}"
        )
    protein = adata.obsm["protein_expression"]
    if isinstance(protein, pd.DataFrame):
        frame = protein.copy()
        frame.index = adata.obs_names
        return frame

    names = None
    if "protein_names" in adata.uns:
        names = np.asarray(adata.uns["protein_names"]).astype(str)
    elif "protein_expression" in adata.uns:
        maybe = adata.uns["protein_expression"]
        if isinstance(maybe, (list, np.ndarray, pd.Index)):
            names = np.asarray(maybe).astype(str)

    array = _as_dense_array(protein)
    if names is None:
        names = np.array([f"protein_{i}" for i in range(array.shape[1])])
        logger.warning("protein_names missing; using placeholder column names.")
    if len(names) != array.shape[1]:
        raise ValueError(
            f"protein_names length {len(names)} != protein matrix width {array.shape[1]}"
        )
    return pd.DataFrame(array, index=adata.obs_names, columns=names)


def clean_protein_name(name: str) -> str:
    """Strip common 10x/TotalSeq suffixes without renaming biology."""
    cleaned = str(name)
    cleaned = re.sub(r"_(TotalSeq[A-Z]|ADT|prot|protein)$", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"[-_]?(TotalSeq[-_]?[A-Z])$", "", cleaned, flags=re.IGNORECASE)
    return cleaned


def infer_celltype_column(obs: pd.DataFrame) -> str | None:
    lower_map = {col.lower(): col for col in obs.columns}
    for candidate in CANDIDATE_CELLTYPE_COLUMNS:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def summarize_matrix(matrix: Any, name: str) -> dict[str, Any]:
    n_obs, n_vars = matrix.shape
    if sparse.issparse(matrix):
        nnz = int(matrix.nnz)
        data = matrix.data
        row_sums = np.asarray(matrix.sum(axis=1)).ravel()
        col_detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    else:
        array = np.asarray(matrix)
        nnz = int(np.count_nonzero(array))
        data = array.ravel()
        row_sums = array.sum(axis=1)
        col_detected = (array > 0).sum(axis=0)

    total = int(n_obs * n_vars)
    sparsity = 1.0 - (nnz / total) if total else np.nan
    finite = np.asarray(data[np.isfinite(data)]) if len(data) else np.array([])
    approx_integer = False
    if finite.size:
        approx_integer = bool(np.all(np.abs(finite - np.round(finite)) < 1e-6))
    return {
        "name": name,
        "n_obs": int(n_obs),
        "n_vars": int(n_vars),
        "nnz": nnz,
        "sparsity": float(sparsity),
        "min": float(np.min(finite)) if finite.size else np.nan,
        "max": float(np.max(finite)) if finite.size else np.nan,
        "mean_counts_per_cell": float(np.mean(row_sums)) if len(row_sums) else np.nan,
        "median_counts_per_cell": float(np.median(row_sums)) if len(row_sums) else np.nan,
        "features_all_zero": int(np.sum(col_detected == 0)),
        "approx_integer": approx_integer,
        "has_negative": bool(finite.size and np.any(finite < 0)),
        "layout": "csr/sparse" if sparse.issparse(matrix) else "dense",
    }


def load_raw_dataset(source_cfg: dict[str, Any], raw_dir: Path) -> ad.AnnData:
    dest = raw_dir / source_cfg["filename"]
    download_file(source_cfg["url"], dest, source_cfg.get("fallback_url"))
    adata = ad.read_h5ad(dest)
    adata.uns["source_key"] = source_cfg["key"]
    adata.uns["source_filename"] = source_cfg["filename"]
    adata.uns["source_url"] = source_cfg["url"]
    return adata


def standardize_single_dataset(
    adata: ad.AnnData,
    batch_key_value: str,
    protein_gene_map: dict[str, str] | None = None,
) -> ad.AnnData:
    """Copy the downloaded object, keep raw RNA/protein counts, add provenance."""
    out = adata.copy()
    out.obs["batch"] = batch_key_value
    out.obs["dataset"] = batch_key_value
    out.obs_names = pd.Index([f"{batch_key_value}_{name}" for name in out.obs_names], dtype=str)
    out.obs_names_make_unique()

    if sparse.issparse(out.X):
        out.layers["counts"] = out.X.copy()
    else:
        out.layers["counts"] = np.array(out.X, copy=True)

    protein = _protein_frame_from_adata(out)
    protein.index = out.obs_names
    protein.columns = [str(c) for c in protein.columns]
    out.obsm["protein_counts"] = protein.copy()
    out.obsm["protein_expression"] = protein.copy()
    out.obsm["protein_observed_mask"] = pd.DataFrame(
        True, index=protein.index, columns=protein.columns
    )
    out.uns["protein_names"] = protein.columns.to_numpy()
    out.uns["protein_names_clean"] = np.array([clean_protein_name(c) for c in protein.columns])

    celltype_col = infer_celltype_column(out.obs)
    out.uns["celltype_column_detected"] = celltype_col
    if celltype_col is not None and celltype_col != "cell_type":
        out.obs["cell_type"] = out.obs[celltype_col].astype(str)
    elif "cell_type" not in out.obs:
        out.obs["cell_type"] = pd.Categorical(["unknown"] * out.n_obs)

    if protein_gene_map:
        mapping = {}
        for raw_name, clean_name in zip(protein.columns, out.uns["protein_names_clean"]):
            mapping[raw_name] = protein_gene_map.get(clean_name)
        out.uns["protein_to_gene"] = mapping

    out.uns["raw_counts_preserved"] = True
    out.uns["x_is_raw_counts"] = True
    return out


def _shared_and_unique(left: pd.Index, right: pd.Index) -> dict[str, list[str]]:
    shared = left.intersection(right)
    return {
        "shared": shared.tolist(),
        "left_only": left.difference(right).tolist(),
        "right_only": right.difference(left).tolist(),
    }


def combine_datasets(
    adata_10k: ad.AnnData,
    adata_5k: ad.AnnData,
    protein_join: str = "inner",
) -> ad.AnnData:
    if protein_join not in {"inner", "outer"}:
        raise ValueError("protein_join must be 'inner' or 'outer'")

    genes = _shared_and_unique(adata_10k.var_names, adata_5k.var_names)
    proteins_10k = adata_10k.obsm["protein_counts"].columns
    proteins_5k = adata_5k.obsm["protein_counts"].columns
    proteins = _shared_and_unique(proteins_10k, proteins_5k)

    rna_10k = adata_10k[:, genes["shared"]].copy()
    rna_5k = adata_5k[:, genes["shared"]].copy()
    combined = ad.concat(
        {"PBMC10k": rna_10k, "PBMC5k": rna_5k},
        join="inner",
        label="concat_batch",
        index_unique=None,
    )
    combined.obs["batch"] = pd.Categorical(combined.obs["dataset"])

    if protein_join == "inner":
        protein_cols = proteins["shared"]
        if not protein_cols:
            raise ValueError("PBMC10k and PBMC5k share no protein names.")
        protein = pd.concat(
            [
                rna_10k.obsm["protein_counts"].loc[:, protein_cols],
                rna_5k.obsm["protein_counts"].loc[:, protein_cols],
            ],
            axis=0,
        )
        mask = pd.DataFrame(True, index=protein.index, columns=protein.columns)
    else:
        protein_cols = proteins_10k.union(proteins_5k).tolist()
        protein = pd.concat(
            [
                rna_10k.obsm["protein_counts"].reindex(columns=protein_cols),
                rna_5k.obsm["protein_counts"].reindex(columns=protein_cols),
            ],
            axis=0,
        )
        mask = ~protein.isna()
        # Keep NaN for proteins that were never measured. Do not fill with 0.

    protein = protein.loc[combined.obs_names, protein_cols]
    mask = mask.loc[combined.obs_names, protein_cols]
    combined.obsm["protein_counts"] = protein
    combined.obsm["protein_observed_mask"] = mask
    # scvi-tools compatibility copy. For the inner join this is identical to
    # protein_counts. For the outer join, NaNs remain; callers must not treat
    # them as biological zeros.
    combined.obsm["protein_expression"] = protein.copy()
    combined.uns["protein_names"] = np.array(protein_cols, dtype=object)
    combined.uns["protein_join"] = protein_join
    combined.uns["gene_panel"] = genes
    combined.uns["protein_panel"] = proteins
    combined.uns["n_genes_dropped_pbmc10k"] = len(genes["left_only"])
    combined.uns["n_genes_dropped_pbmc5k"] = len(genes["right_only"])
    combined.uns["n_proteins_dropped_pbmc10k"] = (
        len(proteins["left_only"]) if protein_join == "inner" else 0
    )
    combined.uns["n_proteins_dropped_pbmc5k"] = (
        len(proteins["right_only"]) if protein_join == "inner" else 0
    )
    combined.uns["raw_counts_preserved"] = True
    combined.uns["x_is_raw_counts"] = True
    if "counts" not in combined.layers:
        combined.layers["counts"] = combined.X.copy()
    return combined


def protein_panel_table(adata_10k: ad.AnnData, adata_5k: ad.AnnData) -> pd.DataFrame:
    map_10k = {
        str(k): v
        for k, v in (adata_10k.uns.get("protein_to_gene") or {}).items()
    }
    map_5k = {
        str(k): v
        for k, v in (adata_5k.uns.get("protein_to_gene") or {}).items()
    }
    names_10k = set(adata_10k.obsm["protein_counts"].columns)
    names_5k = set(adata_5k.obsm["protein_counts"].columns)
    all_names = sorted(names_10k | names_5k)
    rows = []
    for name in all_names:
        rows.append(
            {
                "protein": name,
                "protein_clean": clean_protein_name(name),
                "in_pbmc10k": name in names_10k,
                "in_pbmc5k": name in names_5k,
                "shared": name in names_10k and name in names_5k,
                "mapped_gene": map_10k.get(name) or map_5k.get(name),
            }
        )
    return pd.DataFrame(rows)


def load_citeseq_bundle(config: dict[str, Any] | None = None) -> dict[str, ad.AnnData]:
    """Download, standardize, and combine PBMC10k/PBMC5k CITE-seq datasets."""
    from src.utils.io import load_config

    cfg = config or load_config()
    raw_dir = resolve_path(cfg["paths"]["raw_dir"])
    protein_map = load_yaml(project_root() / "config" / "protein_gene_map.yaml").get(
        "protein_to_gene", {}
    )

    sources = cfg["data"]["sources"]
    raw_10k = load_raw_dataset(sources["pbmc10k"], raw_dir)
    raw_5k = load_raw_dataset(sources["pbmc5k"], raw_dir)
    adata_10k = standardize_single_dataset(raw_10k, sources["pbmc10k"]["key"], protein_map)
    adata_5k = standardize_single_dataset(raw_5k, sources["pbmc5k"]["key"], protein_map)
    combined_inner = combine_datasets(adata_10k, adata_5k, protein_join="inner")
    combined_outer = combine_datasets(adata_10k, adata_5k, protein_join="outer")
    return {
        "pbmc10k": adata_10k,
        "pbmc5k": adata_5k,
        "combined_inner": combined_inner,
        "combined_outer": combined_outer,
        "raw_pbmc10k": raw_10k,
        "raw_pbmc5k": raw_5k,
    }
