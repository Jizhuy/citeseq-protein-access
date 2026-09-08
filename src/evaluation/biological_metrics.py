"""Biological representation metrics. Do not run without validated labels."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    f1_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

from src.annotation.annotation_validator import labels_are_validated

CELLTYPE_UNAVAILABLE_REASON = "Validated ground-truth cell-type labels unavailable"
BIOLOGICAL_METRIC_COLUMNS = [
    "model",
    "representation",
    "label_level",
    "n_cells",
    "n_classes",
    "ARI",
    "NMI",
    "celltype_ASW",
    "knn_purity_k15",
    "knn_purity_k30",
    "knn_purity_k50",
    "knn_accuracy",
    "knn_macro_f1",
    "logreg_accuracy",
    "logreg_macro_f1",
    "annotation_source",
]


class LabelsUnavailableError(ValueError):
    """Raised when biological metrics are requested without validated labels."""


def write_biological_metrics_schema(path) -> None:
    """Write header-only CSV. Do not fill zeros."""
    pd.DataFrame(columns=BIOLOGICAL_METRIC_COLUMNS).to_csv(path, index=False)
    note = path.with_name(path.stem + "_status.txt")
    note.write_text("Pending validated annotations\n", encoding="utf-8")


def clustering_agreement(
    cluster_labels: np.ndarray | pd.Series,
    true_labels: np.ndarray | pd.Series,
) -> dict[str, float]:
    """ARI/NMI between a clustering and validated labels. Not self-comparison."""
    _require_validated(true_labels)
    y_true = np.asarray(true_labels).astype(str)
    y_pred = np.asarray(cluster_labels).astype(str)
    if y_true.shape[0] != y_pred.shape[0]:
        raise ValueError("Cluster assignments and labels have different lengths.")
    return {
        "ARI": float(adjusted_rand_score(y_true, y_pred)),
        "NMI": float(normalized_mutual_info_score(y_true, y_pred)),
    }


def celltype_asw(latent: np.ndarray, labels: np.ndarray | pd.Series) -> float:
    _require_validated(labels)
    y = np.asarray(labels).astype(str)
    if latent.shape[0] != len(y):
        raise ValueError("Latent rows and labels have different lengths.")
    return float(silhouette_score(latent, y, metric="euclidean"))


def neighborhood_purity(
    latent: np.ndarray,
    labels: np.ndarray | pd.Series,
    k: int,
) -> float:
    """Fraction of k nearest neighbors sharing a cell's validated label."""
    _require_validated(labels)
    y = np.asarray(labels).astype(str)
    if latent.shape[0] != len(y):
        raise ValueError("Latent rows and labels have different lengths.")
    if k < 1:
        raise ValueError("k must be >= 1")
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(latent)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    same = np.array([np.mean(y[neigh] == y[i]) for i, neigh in enumerate(idx)])
    return float(same.mean())


def heldout_latent_classification(
    latent: np.ndarray,
    labels: np.ndarray | pd.Series,
    split: np.ndarray | pd.Series,
    *,
    seed: int = 0,
    knn_k: int = 15,
) -> dict[str, Any]:
    """Same train/test cells and same classifier settings for every representation."""
    _require_validated(labels)
    y = np.asarray(labels).astype(str)
    split = np.asarray(split).astype(str)
    if latent.shape[0] != len(y) or len(y) != len(split):
        raise ValueError("latent, labels, and split must have the same length.")
    train = split == "train"
    test = split == "test"
    if train.sum() == 0 or test.sum() == 0:
        raise ValueError("obs['split'] must contain both train and test cells.")
    x_train, x_test = latent[train], latent[test]
    y_train, y_test = y[train], y[test]

    knn = KNeighborsClassifier(n_neighbors=knn_k)
    knn.fit(x_train, y_train)
    knn_pred = knn.predict(x_test)

    logreg = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)
    logreg.fit(x_train, y_train)
    logreg_pred = logreg.predict(x_test)

    return {
        "n_train": int(train.sum()),
        "n_test": int(test.sum()),
        "knn_k": knn_k,
        "knn_accuracy": float(accuracy_score(y_test, knn_pred)),
        "knn_balanced_accuracy": float(balanced_accuracy_score(y_test, knn_pred)),
        "knn_macro_f1": float(f1_score(y_test, knn_pred, average="macro", zero_division=0)),
        "knn_weighted_f1": float(f1_score(y_test, knn_pred, average="weighted", zero_division=0)),
        "logreg_accuracy": float(accuracy_score(y_test, logreg_pred)),
        "logreg_balanced_accuracy": float(balanced_accuracy_score(y_test, logreg_pred)),
        "logreg_macro_f1": float(f1_score(y_test, logreg_pred, average="macro", zero_division=0)),
        "logreg_weighted_f1": float(
            f1_score(y_test, logreg_pred, average="weighted", zero_division=0)
        ),
        "classifier_settings_shared": True,
    }


def evaluate_biological_representation(
    latent: np.ndarray,
    labels: np.ndarray | pd.Series,
    split: np.ndarray | pd.Series,
    *,
    model: str,
    representation: str,
    label_level: str,
    annotation_source: str,
    cluster_labels: np.ndarray | pd.Series | None = None,
    seed: int = 0,
) -> dict[str, Any]:
    """Compute the PHASE 5B metric row. Raises if labels are not validated."""
    _require_validated(labels)
    y = pd.Series(labels).astype(str)
    row: dict[str, Any] = {col: np.nan for col in BIOLOGICAL_METRIC_COLUMNS}
    row.update(
        {
            "model": model,
            "representation": representation,
            "label_level": label_level,
            "n_cells": int(len(y)),
            "n_classes": int(y.nunique()),
            "annotation_source": annotation_source,
            "celltype_ASW": celltype_asw(latent, y),
            "knn_purity_k15": neighborhood_purity(latent, y, 15),
            "knn_purity_k30": neighborhood_purity(latent, y, 30),
            "knn_purity_k50": neighborhood_purity(latent, y, 50),
        }
    )
    if cluster_labels is not None:
        agree = clustering_agreement(cluster_labels, y)
        row["ARI"] = agree["ARI"]
        row["NMI"] = agree["NMI"]
    clf = heldout_latent_classification(latent, y, split, seed=seed)
    row["knn_accuracy"] = clf["knn_accuracy"]
    row["knn_macro_f1"] = clf["knn_macro_f1"]
    row["logreg_accuracy"] = clf["logreg_accuracy"]
    row["logreg_macro_f1"] = clf["logreg_macro_f1"]
    return row


def _require_validated(labels: np.ndarray | pd.Series) -> None:
    if not labels_are_validated(labels):
        raise LabelsUnavailableError(CELLTYPE_UNAVAILABLE_REASON)
