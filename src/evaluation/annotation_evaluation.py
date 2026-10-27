"""PHASE 6 biological evaluation of existing latents. Does not retrain models."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    normalized_mutual_info_score,
    silhouette_samples,
    silhouette_score,
)
from sklearn.neighbors import NearestNeighbors

from src.data.qc_plots import BLUE, GRAY, LIGHT_BLUE, _save
from src.evaluation.representation_metrics import batch_silhouette
from src.utils.io import resolve_path

LEIDEN_RESOLUTIONS = (0.2, 0.4, 0.6, 0.8, 1.0)
PRIMARY_RESOLUTION = 0.6
KNN_KS = (15, 30, 50)
N_BOOT = 500
RARE_TEST_N = 10
MODELS = (
    ("scVI_default", "X_scVI", "RNA", "scvi"),
    ("scVI_matched", "X_scVI_matched", "RNA", "scvi_matched"),
    ("totalVI", "X_totalVI", "RNA+Protein", "totalvi"),
)


def neighborhood_purity_per_cell(latent: np.ndarray, labels: np.ndarray, k: int) -> np.ndarray:
    y = np.asarray(labels).astype(str)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(latent)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    return np.array([np.mean(y[neigh] == y[i]) for i, neigh in enumerate(idx)], dtype=float)


def neighbor_overlap_per_cell(latent_a: np.ndarray, latent_b: np.ndarray, k: int = 15) -> np.ndarray:
    nn_a = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(latent_a)
    nn_b = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(latent_b)
    idx_a = nn_a.kneighbors(return_distance=False)[:, 1:]
    idx_b = nn_b.kneighbors(return_distance=False)[:, 1:]
    return np.array(
        [len(set(a) & set(b)) / float(k) for a, b in zip(idx_a, idx_b)],
        dtype=float,
    )


def leiden_clusters(latent: np.ndarray, resolution: float, seed: int = 0) -> np.ndarray:
    import scanpy as sc

    tmp = ad.AnnData(X=latent.copy())
    sc.pp.neighbors(tmp, use_rep="X", n_neighbors=15, random_state=seed)
    try:
        sc.tl.leiden(tmp, resolution=resolution, random_state=seed, flavor="igraph")
    except TypeError:
        sc.tl.leiden(tmp, resolution=resolution, random_state=seed)
    return tmp.obs["leiden"].astype(str).to_numpy()


def load_latents(
    embeddings_dir: Path,
    obs_names: pd.Index,
    seed: int = 0,
    model_specs: tuple | None = None,
) -> dict[str, np.ndarray]:
    out = {}
    for model, key, _mod, prefix in (model_specs or MODELS):
        npy = embeddings_dir / f"{prefix}_seed{seed}_latent.npy"
        cells = pd.read_csv(embeddings_dir / f"{prefix}_seed{seed}_cells.csv")
        if list(cells["cell_id"].astype(str)) != list(obs_names.astype(str)):
            raise ValueError(f"{prefix} cell order does not match the annotated object.")
        matrix = np.load(npy)
        if matrix.shape[0] != len(obs_names) or matrix.shape[1] != 20:
            raise ValueError(f"{key} has shape {matrix.shape}, expected ({len(obs_names)}, 20).")
        out[model] = matrix
    return out


def _classification_row(clf: dict[str, Any]) -> dict[str, float]:
    return {
        "knn_accuracy": clf["knn_accuracy"],
        "knn_balanced_accuracy": clf["knn_balanced_accuracy"],
        "knn_macro_f1": clf["knn_macro_f1"],
        "knn_weighted_f1": clf["knn_weighted_f1"],
        "logreg_accuracy": clf["logreg_accuracy"],
        "logreg_balanced_accuracy": clf["logreg_balanced_accuracy"],
        "logreg_macro_f1": clf["logreg_macro_f1"],
        "logreg_weighted_f1": clf["logreg_weighted_f1"],
    }


def evaluate_subset(
    latents: dict[str, np.ndarray],
    labels: pd.Series,
    split: pd.Series,
    batch: pd.Series,
    *,
    label_level: str,
    annotation_subset: str,
    annotation_source: str,
    confidence_rule: str,
    seed: int = 0,
    model_specs: tuple | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    y = labels.astype(str)
    keep = y.notna() & ~y.isin({"Unknown", "unknown", ""})
    y = y.loc[keep]
    idx = np.where(keep.to_numpy())[0]
    split = split.astype(str).loc[keep]
    batch = batch.astype(str).loc[keep]
    rows = []
    leiden_rows = []
    extras: dict[str, Any] = {"per_cell": {}, "predictions": {}}
    class_counts = y.value_counts()
    n_classes = int(class_counts.shape[0])
    specs = model_specs or MODELS
    for model, _key, modalities, _prefix in specs:
        latent = latents[model][idx]
        asw = float(silhouette_score(latent, y.to_numpy(), metric="euclidean"))
        sil_samples = silhouette_samples(latent, y.to_numpy(), metric="euclidean")
        purities = {k: neighborhood_purity_per_cell(latent, y.to_numpy(), k) for k in KNN_KS}
        clusters = {res: leiden_clusters(latent, res, seed=seed) for res in LEIDEN_RESOLUTIONS}
        ari_primary = float(adjusted_rand_score(y, clusters[PRIMARY_RESOLUTION]))
        nmi_primary = float(normalized_mutual_info_score(y, clusters[PRIMARY_RESOLUTION]))
        for res, pred in clusters.items():
            leiden_rows.append(
                {
                    "model": model,
                    "label_level": label_level,
                    "annotation_subset": annotation_subset,
                    "resolution": res,
                    "ARI": float(adjusted_rand_score(y, pred)),
                    "NMI": float(normalized_mutual_info_score(y, pred)),
                    "n_clusters": int(pd.Series(pred).nunique()),
                    "n_cells": int(len(y)),
                }
            )
        clf = heldout_with_predictions(latent, y.to_numpy(), split.to_numpy(), seed=seed)
        extras["per_cell"][model] = {
            "asw_samples": sil_samples,
            "purity": purities,
            "index": idx,
        }
        extras["predictions"][model] = clf
        row = {
            "model": model,
            "modalities": modalities,
            "label_level": label_level,
            "annotation_subset": annotation_subset,
            "n_cells": int(len(y)),
            "n_classes": n_classes,
            "ARI": ari_primary,
            "NMI": nmi_primary,
            "leiden_resolution_for_ari_nmi": PRIMARY_RESOLUTION,
            "celltype_ASW": asw,
            "batch_silhouette": float(batch_silhouette(latent, batch)),
            "knn_purity_k15": float(purities[15].mean()),
            "knn_purity_k15_median": float(np.median(purities[15])),
            "knn_purity_k30": float(purities[30].mean()),
            "knn_purity_k30_median": float(np.median(purities[30])),
            "knn_purity_k50": float(purities[50].mean()),
            "knn_purity_k50_median": float(np.median(purities[50])),
            "annotation_source": annotation_source,
            "confidence_rule": confidence_rule,
        }
        row.update(_classification_row(clf))
        rows.append(row)
    extras["labels"] = y.to_numpy()
    extras["split"] = split.to_numpy()
    extras["class_counts"] = class_counts.to_dict()
    extras["keep_index"] = idx
    return pd.DataFrame(rows), pd.DataFrame(leiden_rows), extras


def paired_bootstrap(
    extras_by_model: dict[str, dict[str, Any]],
    *,
    seed: int = 0,
    n_boot: int = N_BOOT,
    model_a: str = "scVI_matched",
    model_b: str = "totalVI",
) -> pd.DataFrame:
    """Same bootstrap indices for every model. Report model_b - model_a."""
    rng = np.random.default_rng(seed)
    matched = extras_by_model[model_a]
    total = extras_by_model[model_b]
    n = len(matched["asw_samples"])
    asw_d, p15_d, p30_d, p50_d = [], [], [], []
    for _ in range(n_boot):
        boot = rng.choice(n, size=n, replace=True)
        asw_d.append(total["asw_samples"][boot].mean() - matched["asw_samples"][boot].mean())
        p15_d.append(total["purity"][15][boot].mean() - matched["purity"][15][boot].mean())
        p30_d.append(total["purity"][30][boot].mean() - matched["purity"][30][boot].mean())
        p50_d.append(total["purity"][50][boot].mean() - matched["purity"][50][boot].mean())

    def _clf_delta(metric: str) -> list[float]:
        y = extras_by_model["scVI_matched"]["y_test"]
        pred_m = extras_by_model["scVI_matched"][metric]
        pred_t = extras_by_model["totalVI"][metric]
        n_test = len(y)
        out = []
        for _ in range(n_boot):
            boot = rng.choice(n_test, size=n_test, replace=True)
            out.append(
                f1_score(y[boot], pred_t[boot], average="macro", zero_division=0)
                - f1_score(y[boot], pred_m[boot], average="macro", zero_division=0)
            )
        return out

    # Classification predictions are stored only as metrics, not raw preds.
    # Recompute deltas from stored per-cell scores for ASW/purity; F1 needs preds.
    rows = []
    for name, draws in (
        ("celltype_ASW", asw_d),
        ("knn_purity_k15", p15_d),
        ("knn_purity_k30", p30_d),
        ("knn_purity_k50", p50_d),
    ):
        arr = np.asarray(draws)
        rows.append(
            {
                "metric": name,
                "n_boot": n_boot,
                "delta_mean": float(arr.mean()),
                "bootstrap_ci_low": float(np.quantile(arr, 0.025)),
                "bootstrap_ci_high": float(np.quantile(arr, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def classification_bootstrap_f1(
    y_test: np.ndarray,
    pred_matched: np.ndarray,
    pred_total: np.ndarray,
    *,
    seed: int = 0,
    n_boot: int = N_BOOT,
    metric_name: str = "logreg_macro_f1",
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(y_test)
    draws = []
    for _ in range(n_boot):
        boot = rng.choice(n, size=n, replace=True)
        draws.append(
            f1_score(y_test[boot], pred_total[boot], average="macro", zero_division=0)
            - f1_score(y_test[boot], pred_matched[boot], average="macro", zero_division=0)
        )
    arr = np.asarray(draws)
    return {
        "metric": metric_name,
        "n_boot": n_boot,
        "delta_mean": float(arr.mean()),
        "bootstrap_ci_low": float(np.quantile(arr, 0.025)),
        "bootstrap_ci_high": float(np.quantile(arr, 0.975)),
    }


def heldout_with_predictions(
    latent: np.ndarray,
    labels: np.ndarray,
    split: np.ndarray,
    seed: int = 0,
) -> dict[str, Any]:
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier

    y = np.asarray(labels).astype(str)
    split = np.asarray(split).astype(str)
    train = split == "train"
    test = split == "test"
    knn = KNeighborsClassifier(n_neighbors=15)
    knn.fit(latent[train], y[train])
    knn_pred = knn.predict(latent[test])
    logreg = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    logreg.fit(latent[train], y[train])
    logreg_pred = logreg.predict(latent[test])
    y_test = y[test]
    return {
        "n_train": int(train.sum()),
        "n_test": int(test.sum()),
        "y_test": y_test,
        "knn_pred": knn_pred,
        "logreg_pred": logreg_pred,
        "knn_accuracy": float(accuracy_score(y_test, knn_pred)),
        "knn_balanced_accuracy": float(balanced_accuracy_score(y_test, knn_pred)),
        "knn_macro_f1": float(f1_score(y_test, knn_pred, average="macro", zero_division=0)),
        "knn_weighted_f1": float(f1_score(y_test, knn_pred, average="weighted", zero_division=0)),
        "logreg_accuracy": float(accuracy_score(y_test, logreg_pred)),
        "logreg_balanced_accuracy": float(balanced_accuracy_score(y_test, logreg_pred)),
        "logreg_macro_f1": float(f1_score(y_test, logreg_pred, average="macro", zero_division=0)),
        "logreg_weighted_f1": float(f1_score(y_test, logreg_pred, average="weighted", zero_division=0)),
        "test_class_counts": pd.Series(y_test).value_counts().to_dict(),
        "rare_test_classes": {
            k: int(v) for k, v in pd.Series(y_test).value_counts().items() if v < RARE_TEST_N
        },
    }


def confusion_frame(y_true: np.ndarray, y_pred: np.ndarray, model: str, label_level: str, subset: str, method: str) -> pd.DataFrame:
    labels = sorted(set(y_true) | set(y_pred))
    mat = confusion_matrix(y_true, y_pred, labels=labels)
    frame = pd.DataFrame(mat, index=labels, columns=labels)
    frame.index.name = "true"
    long = frame.reset_index().melt(id_vars="true", var_name="predicted", value_name="n")
    long["model"] = model
    long["label_level"] = label_level
    long["annotation_subset"] = subset
    long["classifier"] = method
    return long


def _binary_f1(y_true: np.ndarray, y_pred: np.ndarray, cell_type: str) -> dict[str, float]:
    mask = y_true == cell_type
    tp = int(((y_pred == cell_type) & mask).sum())
    fp = int(((y_pred == cell_type) & ~mask).sum())
    fn = int(((y_pred != cell_type) & mask).sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return {
        "n_test": int(mask.sum()),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "too_few_test_examples": bool(mask.sum() < RARE_TEST_N),
    }


def celltype_specific_purity(
    latents: dict[str, np.ndarray],
    labels: pd.Series,
    keep: np.ndarray,
    predictions: dict[str, dict[str, Any]] | None = None,
    model_specs: tuple | None = None,
) -> pd.DataFrame:
    y = labels.astype(str).to_numpy()
    rows = []
    purities = {}
    specs = model_specs or MODELS
    for model, _key, _mod, _p in specs:
        purities[model] = neighborhood_purity_per_cell(latents[model][keep], y, 15)
    for cell_type in sorted(set(y)):
        mask = y == cell_type
        row = {
            "cell_type": cell_type,
            "n_cells": int(mask.sum()),
            "label_level": labels.name or "label",
        }
        for model, _k, _m, _p in specs:
            row[f"purity_{model}"] = float(purities[model][mask].mean())
        if "totalVI" in purities and "scVI_matched" in purities:
            row["delta_purity_totalVI_minus_scVI_matched"] = (
                row["purity_totalVI"] - row["purity_scVI_matched"]
            )
            row["delta_totalVI_minus_scVI_matched"] = row["delta_purity_totalVI_minus_scVI_matched"]
        if "totalVI" in purities and "MOFA+" in purities:
            row["delta_purity_totalVI_minus_MOFA"] = row["purity_totalVI"] - row["purity_MOFA+"]
        if predictions is not None:
            for model, _k, _m, _p in specs:
                pred = predictions[model]
                for method in ("knn", "logreg"):
                    stats = _binary_f1(pred["y_test"], pred[f"{method}_pred"], cell_type)
                    row[f"{method}_n_test"] = stats["n_test"]
                    row[f"{method}_f1_{model}"] = stats["f1"]
                    row[f"{method}_recall_{model}"] = stats["recall"]
                    row[f"{method}_too_few_test"] = stats["too_few_test_examples"]
            if "logreg_f1_totalVI" in row and "logreg_f1_scVI_matched" in row:
                row["delta_logreg_f1_totalVI_minus_scVI_matched"] = (
                    row["logreg_f1_totalVI"] - row["logreg_f1_scVI_matched"]
                )
            if "logreg_f1_totalVI" in row and "logreg_f1_MOFA+" in row:
                row["delta_logreg_f1_totalVI_minus_MOFA"] = (
                    row["logreg_f1_totalVI"] - row["logreg_f1_MOFA+"]
                )
        rows.append(row)
    return pd.DataFrame(rows)
