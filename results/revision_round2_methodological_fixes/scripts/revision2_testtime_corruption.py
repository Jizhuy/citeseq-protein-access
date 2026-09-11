#!/usr/bin/env python3
"""P1 PART C — Post-adaptation test-time protein degradation (Branch A/B).

Uses PHASE 11B adapted_query_model artifacts when present (Branch A).
Otherwise can train dedicated models under dedicated_models_4060/ (Branch B).

Primary analysis: Delta_from_clean within the same clean-adapted model seed.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import shutil
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from revision2_p1_common import CTRL, LOGS, ROOT, TABLES, done_ok, write_done  # noqa: E402

# Project root must be importable for src.experiments.*
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUT = CTRL / "test_time_corruption"
DONE = OUT / "done" / "TESTTIME_CORRUPTION_DONE.json"
MODELS_ROOT = ROOT / "results/phase11b_uncertainty_robustness/models"
PHASE11B_PRIMARY = ROOT / "results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv"
PHASE8_SUMMARY = ROOT / "results/tables/protein_sparsity_primary_delta_summary.csv"

FRACTIONS = (0.0, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85)
MASK_SEEDS = (0, 1, 2, 3, 4)
DIRECTIONS = (
    ("direction_A", "PBMC10k", "PBMC5k"),
    ("direction_B", "PBMC5k", "PBMC10k"),
)
MODEL_NAMES = ("scvi_matched", "totalvi")
PROTEIN_KEY = "protein_expression"
ECE_BINS = 15
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
SCVI_INVAR_TOL = 1e-5
CLEAN_F1_RANGE_PAD = 0.05  # outside historical min/max ± this → flag


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def adapted_dir(direction: str, model: str, seed: int) -> Path:
    return MODELS_ROOT / direction / model / f"seed{seed:02d}" / "adapted_query_model"


def source_dir(direction: str, model: str, seed: int) -> Path:
    return MODELS_ROOT / direction / model / f"seed{seed:02d}" / "source_model"


def models_available() -> dict[str, Any]:
    found = 0
    missing: list[str] = []
    for d, _, _ in DIRECTIONS:
        for m in MODEL_NAMES:
            for s in range(20):
                p = adapted_dir(d, m, s)
                ok = (p / "model.pt").exists()
                if ok:
                    found += 1
                else:
                    missing.append(str(p))
    return {
        "branch": "A" if found == 80 else ("A_partial" if found else "B"),
        "n_found": found,
        "n_expected": 80,
        "n_missing": len(missing),
        "missing_head": missing[:8],
        "models_root": str(MODELS_ROOT),
    }


def ensure_tree() -> None:
    for sub in (
        "tables",
        "logs",
        "manifests",
        "preview_figures",
        "done",
        "per_run",
        "dedicated_models_4060",
    ):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)


def write_gpu_env() -> Path:
    import torch
    import scvi

    lines = [
        "# GPU environment — Part C",
        "",
        f"- timestamp: {_utc()}",
        f"- torch: {torch.__version__}",
        f"- scvi-tools: {scvi.__version__}",
        f"- cuda_available: {torch.cuda.is_available()}",
        f"- cuda: {torch.version.cuda}",
    ]
    if not torch.cuda.is_available():
        lines.append("- ERROR: CUDA unavailable — STOP")
    else:
        lines += [
            f"- gpu: {torch.cuda.get_device_name(0)}",
            f"- vram_bytes: {torch.cuda.get_device_properties(0).total_memory}",
        ]
    path = OUT / "logs" / "gpu_environment.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_server_inventory(avail: dict[str, Any]) -> Path:
    rows = [
        "# Server model inventory — Part C",
        "",
        f"- timestamp: {_utc()}",
        f"- LOCAL/SERVER project root: `{ROOT}`",
        f"- branch: **{avail['branch']}**",
        f"- adapted models found: {avail['n_found']} / {avail['n_expected']}",
        f"- models_root: `{avail['models_root']}`",
        "",
        "## Canonical PHASE 11B grid",
        "",
        "| path | phase | model | direction | seed | loadable | canonical |",
        "|---|---|---|---|---|---|---|",
    ]
    for d, _, _ in DIRECTIONS:
        for m in MODEL_NAMES:
            for s in (0, 19):  # sample endpoints; full count above
                p = adapted_dir(d, m, s)
                ok = (p / "model.pt").exists()
                rows.append(
                    f"| `{p}` | phase11b | {m} | {d} | {s} | {'yes' if ok else 'no'} | yes |"
                )
    rows += [
        "",
        "## Decision",
        "",
        (
            "**BRANCH A** — use historical PHASE 11B adapted_query_model (do not regenerate)."
            if avail["n_found"] == 80
            else "**BRANCH B / partial** — see missing list; dedicated models only if absent."
        ),
        "",
    ]
    path = OUT / "logs" / "server_model_inventory.md"
    path.write_text("\n".join(rows), encoding="utf-8")
    return path


def require_cuda() -> None:
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA unavailable on this host — STOP (no CPU fallback for Part C).")


def checksum_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tensor_hash(model) -> str:
    h = hashlib.sha256()
    for name, p in model.module.state_dict().items():
        h.update(name.encode())
        h.update(p.detach().cpu().numpy().tobytes())
    return h.hexdigest()


def array_checksum(arr: np.ndarray) -> str:
    a = np.ascontiguousarray(arr)
    return hashlib.sha256(a.tobytes()).hexdigest()


def corrupt_protein_entries(protein: np.ndarray, fraction: float, rng: np.random.Generator) -> np.ndarray:
    out = np.array(protein, dtype=np.float64, copy=True, order="C")
    if fraction <= 0:
        return out
    nz = np.flatnonzero(out > 0)
    n_drop = int(round(fraction * nz.size))
    if n_drop <= 0:
        return out
    choose = rng.choice(nz, size=n_drop, replace=False)
    # Index the C-contiguous buffer directly — do NOT use ravel() assignment
    # (ravel may return a copy for non-contiguous inputs, silently no-op'ing).
    out.reshape(-1)[choose] = 0.0
    return out


def apply_protein_corruption(adata, fraction: float, mask_seed: int) -> dict[str, Any]:
    """Corrupt the field totalVI consumes: obsm['protein_expression']. Sync protein_counts."""
    prot = adata.obsm[PROTEIN_KEY]
    if isinstance(prot, pd.DataFrame):
        P0 = prot.to_numpy(dtype=np.float64)
        cols = list(prot.columns)
        index = prot.index
    else:
        P0 = np.asarray(prot, dtype=np.float64)
        cols = None
        index = adata.obs_names
    pre = array_checksum(P0)
    rng = np.random.default_rng(mask_seed)
    P1 = corrupt_protein_entries(P0, fraction, rng)
    post = array_checksum(P1)
    if cols is not None:
        frame = pd.DataFrame(P1, index=index, columns=cols)
        adata.obsm[PROTEIN_KEY] = frame
        if "protein_counts" in adata.obsm:
            adata.obsm["protein_counts"] = frame.copy()
    else:
        adata.obsm[PROTEIN_KEY] = P1
        if "protein_counts" in adata.obsm:
            adata.obsm["protein_counts"] = P1.copy()
    return {
        "protein_checksum_pre": pre,
        "protein_checksum_post": post,
        "protein_changed": pre != post,
        "n_nonzero_pre": int((P0 > 0).sum()),
        "n_nonzero_post": int((P1 > 0).sum()),
        "registered_protein_field": PROTEIN_KEY,
    }


def selective_and_reliability(
    u_pred: np.ndarray,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    probs: np.ndarray,
    classes: np.ndarray,
    u_latent: np.ndarray | None = None,
) -> dict[str, float]:
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        balanced_accuracy_score,
        f1_score,
        roc_auc_score,
    )

    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    error = (y_pred != y_true).astype(int)
    correct = (error == 0).astype(float)
    order = np.argsort(u_pred, kind="mergesort")
    kept = correct[order]
    n = len(correct)
    counts = np.arange(1, n + 1)
    risks = 1.0 - np.cumsum(kept) / counts
    coverage = counts / n
    full_aurc = float(np.mean(risks))
    paurc = float(np.mean(risks[coverage >= 0.5 - 1e-12]))
    n_correct = int(correct.sum())
    perfect = np.zeros(n)
    for k in range(n_correct + 1, n + 1):
        perfect[k - 1] = (k - n_correct) / k
    aurc_star = float(np.mean(perfect)) if n else float("nan")
    e_aurc = full_aurc - aurc_star
    risk_at = {}
    for cov in COVERAGE_GRID:
        keep = max(1, int(round(cov * n)))
        risk_at[f"risk_{int(cov * 100)}"] = float(1.0 - kept[:keep].mean())

    class_to_i = {c: i for i, c in enumerate(np.asarray(classes).astype(str))}
    true_idx = np.array([class_to_i[y] for y in y_true])
    p_true = probs[np.arange(n), true_idx]
    nll = float(-np.mean(np.log(np.clip(p_true, 1e-12, 1.0))))
    Y = np.zeros_like(probs)
    Y[np.arange(n), true_idx] = 1.0
    brier = float(np.mean(np.sum((probs - Y) ** 2, axis=1)))
    conf = probs.max(axis=1)

    def ece(confidence, correct_arr, scheme: str) -> float:
        if scheme == "equal_frequency":
            edges = np.unique(np.quantile(confidence, np.linspace(0, 1, ECE_BINS + 1)))
        else:
            edges = np.linspace(0, 1, ECE_BINS + 1)
        assignment = np.clip(np.digitize(confidence, edges[1:-1], right=True), 0, len(edges) - 2)
        val = 0.0
        for b in range(len(edges) - 1):
            m = assignment == b
            c = int(m.sum())
            if c == 0:
                continue
            val += (c / n) * abs(float(confidence[m].mean()) - float(correct_arr[m].mean()))
        return float(val)

    out: dict[str, float] = {
        "n_cells": float(n),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "weighted_f1": float(f1_score(y_true, y_pred, average="weighted", zero_division=0)),
        "error_prevalence": float(error.mean()),
        "error_auroc": float(roc_auc_score(error, u_pred)) if len(np.unique(error)) > 1 else float("nan"),
        "error_auprc": float(average_precision_score(error, u_pred))
        if len(np.unique(error)) > 1
        else float("nan"),
        "nll": nll,
        "brier": brier,
        "ece_equal_frequency": ece(conf, correct, "equal_frequency"),
        "ece_equal_width": ece(conf, correct, "equal_width"),
        "full_aurc": full_aurc,
        "paurc_0.5_1.0": paurc,
        "e_aurc": e_aurc,
        "aurc_star_perfect_ranking": aurc_star,
        "median_predictive_uncertainty": float(np.median(u_pred)),
        **risk_at,
    }
    if u_latent is not None:
        out["median_latent_uncertainty"] = float(np.median(u_latent))
        if len(np.unique(error)) > 1:
            out["latent_error_auroc"] = float(roc_auc_score(error, u_latent))
            out["latent_error_auprc"] = float(average_precision_score(error, u_latent))
        else:
            out["latent_error_auroc"] = float("nan")
            out["latent_error_auprc"] = float("nan")
    return out


def per_class_f1(y_true, y_pred, classes) -> dict[str, float]:
    from sklearn.metrics import f1_score

    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    out = {}
    for c in classes:
        out[str(c)] = float(f1_score(y_true == c, y_pred == c, zero_division=0))
    return out


def run_unit_path(direction: str, model: str, seed: int) -> Path:
    return OUT / "per_run" / direction / model / f"seed{seed:02d}"


def unit_done(direction: str, model: str, seed: int) -> bool:
    return (run_unit_path(direction, model, seed) / "DONE.json").exists()


def totalvi_unit_complete(direction: str, seed: int) -> bool:
    """Full grid: clean + 6 nonzero fractions × 5 masks."""
    mdir = run_unit_path(direction, "totalvi", seed)
    if not (mdir / "DONE.json").exists():
        return False
    runs = mdir / "runs.csv"
    if not runs.exists():
        return False
    df = pd.read_csv(runs)
    # expect clean row + 30 corrupt rows (6 fracs × 5 masks)
    n_corrupt = df[(df["corruption_fraction"] > 0) & (df["model"] == "totalvi")].shape[0] if "model" in df.columns else df[df["corruption_fraction"] > 0].shape[0]
    n_clean = df[df["corruption_fraction"] <= 0].shape[0] if "corruption_fraction" in df.columns else 0
    # also count via is_clean_baseline if present
    if "is_clean_baseline" in df.columns:
        flag = df["is_clean_baseline"].astype(str).str.lower().isin(["true", "1"])
        n_clean = int(flag.sum())
        n_corrupt = int((~flag).sum())
    return n_clean >= 1 and n_corrupt >= 30


def scvi_unit_complete(direction: str, seed: int) -> bool:
    mdir = run_unit_path(direction, "scvi_matched", seed)
    if not (mdir / "DONE.json").exists():
        return False
    runs = mdir / "runs.csv"
    if not runs.exists():
        return False
    df = pd.read_csv(runs)
    if len(df) < 1:
        return False
    if "scvi_invariance_max_abs" not in df.columns:
        return False
    inv = float(df["scvi_invariance_max_abs"].iloc[0])
    return inv <= SCVI_INVAR_TOL


def prepare_query(model_cls, target, reference_dir: Path):
    query = target.copy()
    model_cls.prepare_query_anndata(query, str(reference_dir))
    return query


def load_adapted(model_cls, adapted_path: Path, query, accelerator: str = "gpu"):
    return model_cls.load(str(adapted_path), adata=query, accelerator=accelerator, device="auto")


def encode_both(model, source, query):
    import torch

    model.module.eval()
    with torch.no_grad():
        z_tgt = np.asarray(model.get_latent_representation())
        mu_t, var_t = model.get_latent_representation(return_dist=True)
        z_src = np.asarray(model.get_latent_representation(adata=source))
        mu_s, var_s = model.get_latent_representation(adata=source, return_dist=True)
    if not np.isfinite(z_tgt).all() or not np.isfinite(z_src).all():
        raise RuntimeError("Non-finite latent")
    if not np.isfinite(var_t).all() or np.any(np.asarray(var_t) <= 0):
        raise RuntimeError("Non-positive/non-finite target posterior variance")
    return {
        "z_source": z_src,
        "z_target": z_tgt,
        "var_target": np.asarray(var_t),
        "var_source": np.asarray(var_s),
        "u_latent_target": np.asarray(var_t).mean(axis=1),
    }


def evaluate_latents(
    z_src_all,
    z_tgt_all,
    u_lat_tgt_all,
    data_slice: dict[str, Any],
    seed: int,
) -> tuple[dict[str, float], dict[str, float], np.ndarray, np.ndarray, np.ndarray]:
    from sklearn.linear_model import LogisticRegression

    si = data_slice["source_index"]
    ti = data_slice["target_index"]
    y_src = data_slice["y_source"]
    y_tgt = data_slice["y_target"]
    z_src = z_src_all[si]
    z_tgt = z_tgt_all[ti]
    u_lat = u_lat_tgt_all[ti]
    clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed, class_weight=None)
    clf.fit(z_src, y_src)
    classes = np.asarray(clf.classes_, dtype=object)
    probs = clf.predict_proba(z_tgt)
    pred = classes[np.argmax(probs, axis=1)]
    u_pred = 1.0 - probs.max(axis=1)
    metrics = selective_and_reliability(u_pred, y_tgt, pred, probs, classes, u_lat)
    pc = per_class_f1(y_tgt, pred, data_slice["classes"])
    return metrics, pc, probs, pred, u_pred


def load_direction_bundle(source_batch: str, target_batch: str):
    from src.experiments.phase10_cross_dataset import class_plan, load_source_target

    source, target, audit = load_source_target(source_batch, target_batch)
    plan = class_plan(source, target, "cell_type_l2")
    high_s = source.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    high_t = target.obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    ys = source.obs["cell_type_l2"].astype(str).to_numpy()
    yt = target.obs["cell_type_l2"].astype(str).to_numpy()
    data_slice = {
        "source_index": plan["source_eval_index"],
        "target_index": plan["target_eval_index"],
        "y_source": ys[plan["source_eval_index"]],
        "y_target": yt[plan["target_eval_index"]],
        "classes": plan["eligible_shared_classes"],
        "y_target_all_high": yt[high_t],
    }
    return source, target, audit, plan, data_slice


def historical_f1_range() -> pd.DataFrame:
    if not PHASE11B_PRIMARY.exists():
        return pd.DataFrame()
    df = pd.read_csv(PHASE11B_PRIMARY)
    return df.groupby(["direction", "model"])["macro_f1"].agg(["min", "max", "mean", "std"]).reset_index()


def classify_clean_vs_historical(direction: str, model: str, f1: float, hist: pd.DataFrame) -> str:
    if hist.empty:
        return "NO_HISTORICAL"
    row = hist[(hist["direction"] == direction) & (hist["model"] == model)]
    if row.empty:
        return "NO_HISTORICAL"
    lo = float(row["min"].iloc[0]) - CLEAN_F1_RANGE_PAD
    hi = float(row["max"].iloc[0]) + CLEAN_F1_RANGE_PAD
    return "CONSISTENT" if lo <= f1 <= hi else "OUTSIDE_HISTORICAL_RANGE"


def run_one_seed(
    direction: str,
    source_batch: str,
    target_batch: str,
    seed: int,
    models: tuple[str, ...] = MODEL_NAMES,
    smoke_only: bool = False,
) -> dict[str, Any]:
    import torch
    from scvi.model import SCVI, TOTALVI

    require_cuda()
    unit_root = run_unit_path(direction, "bundle", seed)
    # per-model paths used below

    source, target, audit, plan, data_slice = load_direction_bundle(source_batch, target_batch)
    rna_cs = array_checksum(
        source.layers["counts"].toarray() if hasattr(source.layers["counts"], "toarray") else np.asarray(source.layers["counts"])
    )
    # cheaper RNA checksum: sum + shape
    def rna_sig(ad):
        X = ad.layers["counts"]
        s = float(X.sum())
        return {"sum": s, "shape": list(X.shape), "obs_hash": hashlib.sha256("".join(ad.obs_names.astype(str)).encode()).hexdigest()[:16]}

    rna_target_sig0 = rna_sig(target)
    hist = historical_f1_range()
    results_meta: dict[str, Any] = {
        "direction": direction,
        "seed": seed,
        "source_batch": source_batch,
        "target_batch": target_batch,
        "n_eligible_classes": len(data_slice["classes"]),
        "audit": {k: audit[k] for k in ("source_n", "target_n", "n_genes", "n_proteins") if k in audit},
        "models": {},
    }

    all_rows: list[dict] = []
    per_class_rows: list[dict] = []
    integrity_rows: list[dict] = []

    for model_name in models:
        mdir = run_unit_path(direction, model_name, seed)
        mdir.mkdir(parents=True, exist_ok=True)
        done_marker = mdir / "DONE.json"
        if not smoke_only and done_marker.exists():
            complete = (
                totalvi_unit_complete(direction, seed)
                if model_name == "totalvi"
                else scvi_unit_complete(direction, seed)
            )
            if complete:
                print(f"skip done {direction} {model_name} seed{seed}")
                cached = mdir / "runs.csv"
                if cached.exists():
                    all_rows.extend(pd.read_csv(cached).to_dict("records"))
                continue
            print(f"incomplete unit — rerunning {direction} {model_name} seed{seed}")
            # remove stale DONE / partial smoke outputs so full grid can rewrite
            done_marker.unlink(missing_ok=True)
            for p in mdir.glob("frac*.json"):
                p.unlink(missing_ok=True)
            if (mdir / "runs.csv").exists():
                (mdir / "runs.csv").unlink()

        ad_path = adapted_dir(direction, model_name, seed)
        src_path = source_dir(direction, model_name, seed)
        if not (ad_path / "model.pt").exists():
            raise FileNotFoundError(ad_path)

        model_cls = SCVI if model_name == "scvi_matched" else TOTALVI
        query = prepare_query(model_cls, target, src_path)
        # registry protein field check for totalVI
        registered = PROTEIN_KEY
        if model_name == "totalvi":
            # corrupt field must match setup
            registered = PROTEIN_KEY

        model = load_adapted(model_cls, ad_path, query, accelerator="gpu")
        device = str(next(model.module.parameters()).device)
        if not device.startswith("cuda"):
            raise RuntimeError(f"{model_name} not on CUDA: {device}")
        h0 = tensor_hash(model)
        file_cs = checksum_file(ad_path / "model.pt")
        n_trainable = sum(p.numel() for p in model.module.parameters() if p.requires_grad)
        n_params = sum(p.numel() for p in model.module.parameters())

        enc = encode_both(model, source, query)
        # provenance log
        provenance = {
            "source_embedding_provenance": "post-adaptation encoding via adapted_query_model",
            "target_embedding_provenance": "post-adaptation encoding via adapted_query_model",
            "adapted_model_path": str(ad_path),
            "source_model_path": str(src_path),
            "coordinate_rule": "BOTH source and target encoded with SAME adapted query model",
            "labels_used_before_eval": False,
            "target_labels_in_training": False,
        }
        (mdir / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

        # freeze classifier on CLEAN source post-adapt
        metrics_clean, pc_clean, probs_c, pred_c, u_pred_c = evaluate_latents(
            enc["z_source"], enc["z_target"], enc["u_latent_target"], data_slice, seed
        )
        clean_status = classify_clean_vs_historical(direction, model_name, metrics_clean["macro_f1"], hist)
        if clean_status == "OUTSIDE_HISTORICAL_RANGE":
            # soft flag — stop only if grossly inconsistent (>0.15 from mean)
            if not hist.empty:
                row = hist[(hist["direction"] == direction) & (hist["model"] == model_name)]
                if not row.empty and abs(metrics_clean["macro_f1"] - float(row["mean"].iloc[0])) > 0.15:
                    raise RuntimeError(
                        f"Clean F1 grossly inconsistent: {metrics_clean['macro_f1']} vs hist mean {row['mean'].iloc[0]}"
                    )

        h1 = tensor_hash(model)
        if h1 != h0:
            raise RuntimeError("Model hash changed during clean encoding")

        base_row = {
            "direction": direction,
            "model": model_name,
            "seed": seed,
            "corruption_fraction": 0.0,
            "mask_seed": -1,
            "is_clean_baseline": True,
            "clean_vs_historical": clean_status,
            "param_hash": h0,
            "model_pt_sha256": file_cs,
            "device": device,
            "n_params": n_params,
            "n_trainable_at_load": int(n_trainable),
            "param_hash_unchanged": True,
            "rna_unchanged": True,
            "scvi_invariance_max_abs": 0.0 if model_name == "scvi_matched" else float("nan"),
            **{k: metrics_clean[k] for k in metrics_clean},
        }
        all_rows.append(base_row)
        for ct, f1 in pc_clean.items():
            per_class_rows.append(
                {
                    "direction": direction,
                    "model": model_name,
                    "seed": seed,
                    "corruption_fraction": 0.0,
                    "mask_seed": -1,
                    "cell_type": ct,
                    "f1": f1,
                    "median_pred_unc": float("nan"),
                    "median_latent_unc": float("nan"),
                }
            )

        # scVI negative control + reuse
        if model_name == "scvi_matched":
            q_corrupt = query.copy()
            corr_info = apply_protein_corruption(q_corrupt, 0.25, 0)
            import torch as _t

            model.module.eval()
            with _t.no_grad():
                z_cor = np.asarray(model.get_latent_representation(adata=q_corrupt))
            max_abs = float(np.max(np.abs(z_cor - enc["z_target"])))
            base_row["scvi_invariance_max_abs"] = max_abs
            all_rows[-1]["scvi_invariance_max_abs"] = max_abs
            if max_abs > SCVI_INVAR_TOL:
                raise RuntimeError(f"scVI changed under protein corruption: max_abs={max_abs}")
            h2 = tensor_hash(model)
            if h2 != h0:
                raise RuntimeError("scVI hash changed")
            # replicate clean metrics across fractions for table completeness (no recompute)
            if not smoke_only:
                for frac in FRACTIONS[1:]:
                    row = dict(base_row)
                    row.update(
                        {
                            "corruption_fraction": frac,
                            "mask_seed": 0,
                            "is_clean_baseline": False,
                            "note": "scVI_protein_insensitive_reused_clean",
                        }
                    )
                    all_rows.append(row)
            # integrity
            integrity_rows.append(
                {
                    "direction": direction,
                    "model": model_name,
                    "seed": seed,
                    "param_hash": h0,
                    "hash_unchanged_after_corrupt_probe": h2 == h0,
                    "scvi_invariance_max_abs": max_abs,
                    "file_sha256": file_cs,
                }
            )
            pd.DataFrame([r for r in all_rows if r["model"] == model_name and r["seed"] == seed]).to_csv(
                mdir / "runs.csv", index=False
            )
            if smoke_only:
                write_done(
                    mdir / "SMOKE_DONE.json",
                    {"ok": True, "smoke": True, "direction": direction, "model": model_name, "seed": seed, "scvi_invariance": max_abs},
                )
            else:
                write_done(
                    done_marker,
                    {"ok": True, "direction": direction, "model": model_name, "seed": seed, "scvi_invariance": max_abs},
                )
            results_meta["models"][model_name] = {
                "clean_macro_f1": metrics_clean["macro_f1"],
                "clean_vs_historical": clean_status,
                "param_hash": h0,
                "scvi_invariance_max_abs": max_abs,
            }
            del model, query, q_corrupt
            gc.collect()
            torch.cuda.empty_cache()
            continue

        # totalVI corruption grid — reuse one prepared query; swap protein matrix only
        prot0 = query.obsm[PROTEIN_KEY]
        if isinstance(prot0, pd.DataFrame):
            P_clean = prot0.to_numpy(dtype=np.float64).copy()
            prot_cols = list(prot0.columns)
            prot_index = prot0.index
        else:
            P_clean = np.asarray(prot0, dtype=np.float64).copy()
            prot_cols = None
            prot_index = query.obs_names

        if smoke_only:
            frac_list = [0.25]
            mask_list = [0]
        else:
            frac_list = list(FRACTIONS)
            mask_list = list(MASK_SEEDS)

        for frac in frac_list:
            masks = [-1] if frac == 0.0 else mask_list
            for ms in masks:
                if frac == 0.0 and not smoke_only:
                    continue  # already have clean baseline
                out_json = mdir / f"frac{int(round(frac*100)):03d}_mask{ms}.json"
                if out_json.exists() and not smoke_only:
                    all_rows.append(json.loads(out_json.read_text(encoding="utf-8")))
                    continue

                rna_before = rna_sig(query)
                if frac > 0:
                    rng = np.random.default_rng(int(ms))
                    P1 = corrupt_protein_entries(P_clean, frac, rng)
                    corr_info = {
                        "protein_checksum_pre": array_checksum(P_clean),
                        "protein_checksum_post": array_checksum(P1),
                        "protein_changed": True,
                        "registered_protein_field": registered,
                    }
                else:
                    P1 = P_clean
                    cs = array_checksum(P_clean)
                    corr_info = {
                        "protein_checksum_pre": cs,
                        "protein_checksum_post": cs,
                        "protein_changed": False,
                        "registered_protein_field": registered,
                    }
                if prot_cols is not None:
                    frame = pd.DataFrame(P1, index=prot_index, columns=prot_cols)
                    query.obsm[PROTEIN_KEY] = frame
                    if "protein_counts" in query.obsm:
                        query.obsm["protein_counts"] = frame.copy()
                else:
                    query.obsm[PROTEIN_KEY] = P1
                    if "protein_counts" in query.obsm:
                        query.obsm["protein_counts"] = np.array(P1, copy=True)

                rna_after = rna_sig(query)
                if rna_before != rna_after:
                    raise RuntimeError("RNA changed during protein corruption")

                model.module.eval()
                with torch.no_grad():
                    z_t = np.asarray(model.get_latent_representation(adata=query))
                    _, var_t = model.get_latent_representation(adata=query, return_dist=True)
                u_lat = np.asarray(var_t).mean(axis=1)
                if not np.isfinite(z_t).all() or not np.isfinite(u_lat).all():
                    raise RuntimeError("Non-finite corrupted latent")

                h_after = tensor_hash(model)
                if h_after != h0:
                    raise RuntimeError("Model parameters changed during corrupted inference")

                from sklearn.linear_model import LogisticRegression

                si = data_slice["source_index"]
                ti = data_slice["target_index"]
                clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed, class_weight=None)
                clf.fit(enc["z_source"][si], data_slice["y_source"])
                classes = np.asarray(clf.classes_, dtype=object)
                probs = clf.predict_proba(z_t[ti])
                pred = classes[np.argmax(probs, axis=1)]
                u_pred = 1.0 - probs.max(axis=1)
                metrics = selective_and_reliability(
                    u_pred, data_slice["y_target"], pred, probs, classes, u_lat[ti]
                )
                pc = per_class_f1(data_slice["y_target"], pred, data_slice["classes"])

                y_tgt = np.asarray(data_slice["y_target"]).astype(str)
                for ct in data_slice["classes"]:
                    mask_ct = y_tgt == ct
                    per_class_rows.append(
                        {
                            "direction": direction,
                            "model": model_name,
                            "seed": seed,
                            "corruption_fraction": frac,
                            "mask_seed": ms,
                            "cell_type": ct,
                            "f1": pc[str(ct)],
                            "median_pred_unc": float(np.median(u_pred[mask_ct])) if mask_ct.any() else float("nan"),
                            "median_latent_unc": float(np.median(u_lat[ti][mask_ct])) if mask_ct.any() else float("nan"),
                            "support": int(mask_ct.sum()),
                        }
                    )

                row = {
                    "direction": direction,
                    "model": model_name,
                    "seed": seed,
                    "corruption_fraction": frac,
                    "mask_seed": ms,
                    "is_clean_baseline": False,
                    "clean_vs_historical": clean_status,
                    "param_hash": h0,
                    "param_hash_after": h_after,
                    "param_hash_unchanged": h_after == h0,
                    "model_pt_sha256": file_cs,
                    "device": device,
                    "rna_unchanged": True,
                    "protein_changed": bool(corr_info.get("protein_changed")),
                    "registered_protein_field": registered,
                    "protein_checksum_pre": corr_info.get("protein_checksum_pre"),
                    "protein_checksum_post": corr_info.get("protein_checksum_post"),
                    "delta_macro_f1_from_clean": metrics["macro_f1"] - metrics_clean["macro_f1"],
                    "delta_error_auroc_from_clean": metrics["error_auroc"] - metrics_clean["error_auroc"],
                    "delta_e_aurc_from_clean": metrics["e_aurc"] - metrics_clean["e_aurc"],
                    "delta_nll_from_clean": metrics["nll"] - metrics_clean["nll"],
                    "delta_brier_from_clean": metrics["brier"] - metrics_clean["brier"],
                    "clean_macro_f1": metrics_clean["macro_f1"],
                    **{k: metrics[k] for k in metrics},
                }
                out_json.write_text(json.dumps(row, indent=2), encoding="utf-8")
                all_rows.append(row)
                integrity_rows.append(
                    {
                        "direction": direction,
                        "model": model_name,
                        "seed": seed,
                        "corruption_fraction": frac,
                        "mask_seed": ms,
                        "param_hash_unchanged": True,
                        "rna_unchanged": True,
                        "protein_changed": bool(corr_info.get("protein_changed")),
                    }
                )
                del z_t, var_t, u_lat, probs, pred, u_pred
                gc.collect()

            if smoke_only:
                break

        # restore clean protein on query
        if prot_cols is not None:
            query.obsm[PROTEIN_KEY] = pd.DataFrame(P_clean, index=prot_index, columns=prot_cols)
        else:
            query.obsm[PROTEIN_KEY] = P_clean

        pd.DataFrame([r for r in all_rows if r.get("model") == model_name and r.get("seed") == seed]).to_csv(
            mdir / "runs.csv", index=False
        )
        if smoke_only:
            write_done(
                mdir / "SMOKE_DONE.json",
                {
                    "ok": True,
                    "smoke": True,
                    "direction": direction,
                    "model": model_name,
                    "seed": seed,
                    "clean_macro_f1": metrics_clean["macro_f1"],
                    "param_hash": h0,
                },
            )
        else:
            write_done(
                done_marker,
                {
                    "ok": True,
                    "direction": direction,
                    "model": model_name,
                    "seed": seed,
                    "clean_macro_f1": metrics_clean["macro_f1"],
                    "clean_vs_historical": clean_status,
                    "param_hash": h0,
                },
            )
        results_meta["models"][model_name] = {
            "clean_macro_f1": metrics_clean["macro_f1"],
            "clean_vs_historical": clean_status,
            "param_hash": h0,
        }
        del model, query
        gc.collect()
        torch.cuda.empty_cache()

    # persist unit rows
    bundle = OUT / "per_run" / direction / f"seed{seed:02d}_rows.json"
    bundle.parent.mkdir(parents=True, exist_ok=True)
    bundle.write_text(json.dumps({"rows": all_rows, "per_class": per_class_rows, "integrity": integrity_rows}, indent=2), encoding="utf-8")
    del source, target
    gc.collect()
    return {"rows": all_rows, "per_class": per_class_rows, "integrity": integrity_rows, "meta": results_meta}


def smoke_test() -> dict[str, Any]:
    """Direction B then A, seed0, totalVI 0.25/mask0; then scVI invariance."""
    out = {}
    for direction, src, tgt in (("direction_B", "PBMC5k", "PBMC10k"), ("direction_A", "PBMC10k", "PBMC5k")):
        print(f"=== SMOKE {direction} totalVI+scVI seed0 ===", flush=True)
        res = run_one_seed(direction, src, tgt, 0, models=("totalvi", "scvi_matched"), smoke_only=True)
        out[direction] = {
            "n_rows": len(res["rows"]),
            "meta": res["meta"],
            "integrity_ok": all(i.get("param_hash_unchanged", True) for i in res["integrity"]),
        }
        scvi_rows = [r for r in res["rows"] if r["model"] == "scvi_matched"]
        if scvi_rows:
            inv = scvi_rows[0].get("scvi_invariance_max_abs", 1)
            if inv > SCVI_INVAR_TOL:
                raise RuntimeError(f"Smoke scVI invariance fail: {inv}")
            out[direction]["scvi_invariance_max_abs"] = inv
        tot = [r for r in res["rows"] if r["model"] == "totalvi" and r.get("corruption_fraction") == 0.25]
        if not tot:
            raise RuntimeError("Smoke missing totalVI 0.25 row")
        if not tot[0].get("param_hash_unchanged", False):
            raise RuntimeError("Smoke totalVI hash changed")
        out[direction]["totalvi_f1_clean"] = res["meta"]["models"]["totalvi"]["clean_macro_f1"]
        out[direction]["totalvi_f1_0.25"] = tot[0]["macro_f1"]
    write_done(OUT / "done" / "SMOKE_PASS.json", {"ok": True, **out})
    (OUT / "logs" / "smoke_test.md").write_text(
        "# Smoke test PASS\n\n```json\n" + json.dumps(out, indent=2) + "\n```\n", encoding="utf-8"
    )
    return out


def aggregate_and_write(all_rows: list[dict], per_class: list[dict], integrity: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.DataFrame(all_rows)
    if df.empty:
        raise RuntimeError("No rows to aggregate")
    if "is_clean_baseline" in df.columns:
        df["is_clean_baseline"] = df["is_clean_baseline"].astype(str).str.lower().isin(["true", "1"])
    df.to_csv(OUT / "tables" / "test_time_corruption_all_runs.csv", index=False)
    df.to_csv(TABLES / "test_time_corruption_all_runs.csv", index=False)

    # deltas from clean
    clean_cols = [
        c
        for c in [
            "direction",
            "model",
            "seed",
            "macro_f1",
            "error_auroc",
            "e_aurc",
            "nll",
            "brier",
            "full_aurc",
            "paurc_0.5_1.0",
            "latent_error_auroc",
        ]
        if c in df.columns
    ]
    clean = df[df["is_clean_baseline"]][clean_cols].rename(
        columns={
            "macro_f1": "clean_macro_f1",
            "error_auroc": "clean_error_auroc",
            "e_aurc": "clean_e_aurc",
            "nll": "clean_nll",
            "brier": "clean_brier",
            "full_aurc": "clean_full_aurc",
            "paurc_0.5_1.0": "clean_paurc",
            "latent_error_auroc": "clean_latent_error_auroc",
        }
    )
    merged = df.merge(clean, on=["direction", "model", "seed"], how="left")
    for a, b, dname in (
        ("macro_f1", "clean_macro_f1", "delta_macro_f1"),
        ("error_auroc", "clean_error_auroc", "delta_error_auroc"),
        ("e_aurc", "clean_e_aurc", "delta_e_aurc"),
        ("nll", "clean_nll", "delta_nll"),
        ("brier", "clean_brier", "delta_brier"),
        ("full_aurc", "clean_full_aurc", "delta_full_aurc"),
        ("paurc_0.5_1.0", "clean_paurc", "delta_paurc"),
        ("latent_error_auroc", "clean_latent_error_auroc", "delta_latent_error_auroc"),
    ):
        if a in merged.columns and b in merged.columns:
            merged[dname] = merged[a] - merged[b]
        elif "delta_macro_f1_from_clean" in merged.columns and dname == "delta_macro_f1":
            merged[dname] = merged["delta_macro_f1_from_clean"]
        else:
            merged[dname] = np.nan
    merged.to_csv(OUT / "tables" / "test_time_corruption_delta_from_clean.csv", index=False)

    # within seed average over masks, then across seeds
    tot = merged[(merged["model"] == "totalvi") & (~merged["is_clean_baseline"])].copy()
    # include clean as fraction 0
    tot0 = merged[(merged["model"] == "totalvi") & (merged["is_clean_baseline"])].copy()
    tot0["corruption_fraction"] = 0.0
    tot0["mask_seed"] = -1
    tot0["delta_macro_f1"] = 0.0
    tot_all = pd.concat([tot0, tot], ignore_index=True)

    metric_cols = [
        c
        for c in [
            "macro_f1",
            "delta_macro_f1",
            "error_auroc",
            "delta_error_auroc",
            "latent_error_auroc",
            "nll",
            "brier",
            "ece_equal_frequency",
            "ece_equal_width",
            "full_aurc",
            "paurc_0.5_1.0",
            "e_aurc",
            "delta_e_aurc",
            "balanced_accuracy",
            "accuracy",
            "error_prevalence",
            "median_predictive_uncertainty",
            "median_latent_uncertainty",
        ]
        if c in tot_all.columns
    ]

    within = tot_all.groupby(["direction", "model", "seed", "corruption_fraction"], as_index=False)[metric_cols].mean()

    summary_rows = []
    for (direction, model, frac), g in within.groupby(["direction", "model", "corruption_fraction"]):
        row = {"direction": direction, "model": model, "corruption_fraction": frac, "n_seeds": len(g)}
        for c in metric_cols:
            row[f"{c}_mean"] = float(g[c].mean())
            row[f"{c}_sd"] = float(g[c].std(ddof=1)) if len(g) > 1 else 0.0
            row[f"{c}_median"] = float(g[c].median())
            q1, q3 = g[c].quantile([0.25, 0.75])
            row[f"{c}_iqr"] = float(q3 - q1)
        # within-seed mask variation for this fraction (nonzero only)
        sub = tot_all[
            (tot_all["direction"] == direction)
            & (tot_all["model"] == model)
            & (tot_all["corruption_fraction"] == frac)
            & (tot_all["mask_seed"] >= 0)
        ]
        if len(sub):
            # SD of masks within each seed, then mean across seeds
            mask_sds = sub.groupby("seed")["macro_f1"].std(ddof=1)
            row["within_seed_mask_sd_macro_f1_mean"] = float(mask_sds.mean())
        else:
            row["within_seed_mask_sd_macro_f1_mean"] = 0.0
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "tables" / "test_time_corruption_summary.csv", index=False)
    summary.to_csv(TABLES / "test_time_corruption_summary.csv", index=False)

    pc = pd.DataFrame(per_class)
    if len(pc):
        pc.to_csv(OUT / "tables" / "test_time_corruption_per_class.csv", index=False)

    # reliability table = all_runs reliability columns
    rel_cols = [
        c
        for c in df.columns
        if c
        in {
            "direction",
            "model",
            "seed",
            "corruption_fraction",
            "mask_seed",
            "nll",
            "brier",
            "ece_equal_frequency",
            "ece_equal_width",
            "full_aurc",
            "paurc_0.5_1.0",
            "e_aurc",
            "error_auroc",
            "error_auprc",
            "latent_error_auroc",
            "latent_error_auprc",
            "risk_100",
            "risk_90",
            "risk_80",
            "risk_70",
            "risk_60",
            "risk_50",
        }
    ]
    df[rel_cols].to_csv(OUT / "tables" / "test_time_reliability_metrics.csv", index=False)

    # clean vs phase11b
    hist = pd.read_csv(PHASE11B_PRIMARY) if PHASE11B_PRIMARY.exists() else pd.DataFrame()
    clean_df = df[df["is_clean_baseline"]].copy()
    cmp_rows = []
    for _, r in clean_df.iterrows():
        row = {
            "direction": r["direction"],
            "model": r["model"],
            "seed": r["seed"],
            "local_clean_macro_f1": r["macro_f1"],
            "clean_vs_historical": r.get("clean_vs_historical"),
        }
        if len(hist):
            h = hist[(hist["direction"] == r["direction"]) & (hist["model"] == r["model"]) & (hist["seed"] == r["seed"])]
            if len(h):
                row["phase11b_macro_f1"] = float(h["macro_f1"].iloc[0])
                row["abs_diff"] = abs(row["local_clean_macro_f1"] - row["phase11b_macro_f1"])
        cmp_rows.append(row)
    pd.DataFrame(cmp_rows).to_csv(OUT / "tables" / "clean_baseline_vs_phase10.csv", index=False)
    # filename requested clean_baseline_vs_phase10 — content is vs phase11b primary (canonical corrected)

    # training vs test-time
    for c in [
        "delta_macro_f1_mean",
        "error_auroc_mean",
        "latent_error_auroc_mean",
        "full_aurc_mean",
        "e_aurc_mean",
        "macro_f1_mean",
    ]:
        if c not in summary.columns:
            summary[c] = np.nan
    tt = summary[summary["model"] == "totalvi"][
        [
            "direction",
            "corruption_fraction",
            "macro_f1_mean",
            "delta_macro_f1_mean",
            "error_auroc_mean",
            "latent_error_auroc_mean",
            "full_aurc_mean",
            "e_aurc_mean",
        ]
    ].copy()
    if PHASE8_SUMMARY.exists():
        p8 = pd.read_csv(PHASE8_SUMMARY)
        # phase8 is development within-dataset, not transfer — still compare numerically as requested
        cmp = []
        for _, r in tt.iterrows():
            frac = float(r["corruption_fraction"])
            p8r = p8[np.isclose(p8["corruption_fraction"], frac)]
            scvi_ref = float("nan")
            scvi_clean = clean_df[(clean_df["direction"] == r["direction"]) & (clean_df["model"] == "scvi_matched")]
            if len(scvi_clean):
                scvi_ref = float(scvi_clean["macro_f1"].mean())
            cmp.append(
                {
                    "direction": r["direction"],
                    "corruption_fraction": frac,
                    "historical_training_time_F1": float(p8r["totalVI_l2_macro_f1_mean"].iloc[0]) if len(p8r) else float("nan"),
                    "new_test_time_F1": float(r["macro_f1_mean"]),
                    "new_delta_from_clean": float(r["delta_macro_f1_mean"])
                    if pd.notna(r["delta_macro_f1_mean"])
                    else float("nan"),
                    "scVI_reference_clean_transfer": scvi_ref,
                    "predictive_AUROC": float(r["error_auroc_mean"]) if pd.notna(r["error_auroc_mean"]) else float("nan"),
                    "latent_AUROC": float(r["latent_error_auroc_mean"])
                    if pd.notna(r["latent_error_auroc_mean"])
                    else float("nan"),
                    "full_AURC": float(r["full_aurc_mean"]) if pd.notna(r["full_aurc_mean"]) else float("nan"),
                    "E_AURC": float(r["e_aurc_mean"]) if pd.notna(r["e_aurc_mean"]) else float("nan"),
                    "note": "Phase8=training-time within-dataset HC; PartC=test-time transfer — designs differ",
                }
            )
        pd.DataFrame(cmp).to_csv(OUT / "tables" / "training_vs_testtime_corruption.csv", index=False)
    else:
        tt.to_csv(OUT / "tables" / "training_vs_testtime_corruption.csv", index=False)

    # integrity log
    integ = pd.DataFrame(integrity)
    integ_path = OUT / "logs" / "frozen_model_integrity.md"
    n_bad = 0
    if len(integ) and "param_hash_unchanged" in integ.columns:
        flags = integ["param_hash_unchanged"].astype(str).str.lower().isin(["true", "1"])
        n_bad = int((~flags).sum())
    integ_path.write_text(
        f"# Frozen model integrity\n\n- rows: {len(integ)}\n- hash failures: {n_bad}\n- timestamp: {_utc()}\n",
        encoding="utf-8",
    )

    # preview figures
    figdir = OUT / "preview_figures"
    figdir.mkdir(parents=True, exist_ok=True)

    def _plot(metric_mean, ylabel, fname):
        fig, ax = plt.subplots(figsize=(6, 4))
        for direction, g in summary[summary["model"] == "totalvi"].groupby("direction"):
            g = g.sort_values("corruption_fraction")
            ax.plot(g["corruption_fraction"], g[metric_mean], marker="o", label=direction)
        ax.set_xlabel("corruption fraction")
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(ylabel)
        fig.tight_layout()
        fig.savefig(figdir / fname, dpi=150)
        plt.close(fig)

    if "macro_f1_mean" in summary.columns:
        _plot("macro_f1_mean", "totalVI macro-F1", "macro_f1_vs_corruption.png")
    if "delta_macro_f1_mean" in summary.columns:
        _plot("delta_macro_f1_mean", "Delta F1 from clean", "delta_f1_vs_corruption.png")
    if "error_auroc_mean" in summary.columns:
        _plot("error_auroc_mean", "Predictive error AUROC", "pred_auroc_vs_corruption.png")
    if "latent_error_auroc_mean" in summary.columns:
        _plot("latent_error_auroc_mean", "Latent error AUROC", "latent_auroc_vs_corruption.png")
    if "e_aurc_mean" in summary.columns:
        _plot("e_aurc_mean", "E-AURC", "e_aurc_vs_corruption.png")

    # training vs test-time plot
    tvt = OUT / "tables" / "training_vs_testtime_corruption.csv"
    if tvt.exists():
        tdf = pd.read_csv(tvt)
        fig, ax = plt.subplots(figsize=(6, 4))
        for direction, g in tdf.groupby("direction"):
            g = g.sort_values("corruption_fraction")
            ax.plot(g["corruption_fraction"], g["new_test_time_F1"], marker="o", label=f"test-time {direction}")
        if "historical_training_time_F1" in tdf.columns:
            g0 = tdf.drop_duplicates("corruption_fraction").sort_values("corruption_fraction")
            ax.plot(g0["corruption_fraction"], g0["historical_training_time_F1"], marker="s", label="train-time Phase8")
        ax.legend()
        ax.set_xlabel("corruption fraction")
        ax.set_ylabel("macro-F1")
        fig.tight_layout()
        fig.savefig(figdir / "training_vs_testtime.png", dpi=150)
        plt.close(fig)

    if len(pc):
        fig, ax = plt.subplots(figsize=(8, 5))
        sub = pc[(pc["model"] == "totalvi") & (pc["direction"] == "direction_B")]
        if len(sub):
            for ct, g in sub.groupby("cell_type"):
                gg = g.groupby("corruption_fraction")["f1"].mean().reset_index()
                ax.plot(gg["corruption_fraction"], gg["f1"], marker="o", label=ct, alpha=0.8)
            ax.legend(fontsize=6, ncol=2)
            ax.set_title("Per-class F1 (Dir B totalVI)")
            fig.tight_layout()
            fig.savefig(figdir / "per_class_f1_dirB.png", dpi=150)
            plt.close(fig)

    # report
    report = [
        "# Test-time corruption report",
        "",
        f"- timestamp: {_utc()}",
        f"- n_rows: {len(df)}",
        f"- branch: A (PHASE 11B adapted models)",
        "",
        "## Clean baselines (mean over completed seeds)",
        "",
    ]
    for (d, m), g in clean_df.groupby(["direction", "model"]):
        report.append(f"- {d} {m}: macro-F1={g['macro_f1'].mean():.4f} (n={len(g)})")
    report += ["", "## totalVI mean macro-F1 by corruption", ""]
    for direction, g in summary[summary["model"] == "totalvi"].groupby("direction"):
        report.append(f"### {direction}")
        for _, r in g.sort_values("corruption_fraction").iterrows():
            report.append(f"- p={r['corruption_fraction']}: F1={r['macro_f1_mean']:.4f}")
    (OUT / "logs" / "test_time_corruption_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    # clean baseline validation md
    cmp_df = pd.DataFrame(cmp_rows)
    lines = ["# Clean baseline validation", "", f"- timestamp: {_utc()}", ""]
    if len(cmp_df):
        lines.append(cmp_df.to_string(index=False))
    else:
        lines.append("See `tables/clean_baseline_vs_phase10.csv`.")
    lines.append("")
    (OUT / "logs" / "clean_baseline_validation.md").write_text("\n".join(lines), encoding="utf-8")


def collect_all_cached() -> tuple[list[dict], list[dict], list[dict]]:
    rows, pc, integ = [], [], []
    for d, _, _ in DIRECTIONS:
        for m in MODEL_NAMES:
            for s in range(20):
                p = run_unit_path(d, m, s) / "runs.csv"
                if p.exists():
                    rows.extend(pd.read_csv(p).to_dict("records"))
                b = OUT / "per_run" / d / f"seed{s:02d}_rows.json"
                if b.exists():
                    js = json.loads(b.read_text(encoding="utf-8"))
                    pc.extend(js.get("per_class", []))
                    integ.extend(js.get("integrity", []))
    return rows, pc, integ


def append_synthesis(summary_path: Path, state: str, bullets: list[str]) -> None:
    synth = LOGS / "P1_CONTROL_SYNTHESIS.md"
    if not synth.exists():
        synth.write_text("# P1 Control Synthesis\n", encoding="utf-8")
    text = synth.read_text(encoding="utf-8")
    if "## PART C RESULTS" in text:
        # replace from PART C onward
        text = text.split("## PART C RESULTS")[0].rstrip() + "\n\n"
    block = ["## PART C RESULTS", "", f"- completed: {_utc()}", f"- final P1 STATE: **{state}**", ""]
    block.extend(bullets)
    block.append("")
    synth.write_text(text + "\n".join(block), encoding="utf-8")


def gate_seeds(gate: str) -> list[int]:
    g = gate.upper()
    if g == "A":
        return [0]
    if g == "B":
        return list(range(5))
    if g == "C":
        return list(range(20))
    raise ValueError(gate)


def main():
    ap = argparse.ArgumentParser(description="P1 Part C test-time protein corruption")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--gate", default="A", help="A=seed0, B=0-4, C=0-19")
    ap.add_argument("--aggregate-only", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--one-seed", nargs=3, metavar=("DIRECTION", "SRC_BATCH", "SEED"), help="Worker: run one direction×seed in-process")
    args = ap.parse_args()

    ensure_tree()
    avail = models_available()

    if args.status:
        print("DONE" if done_ok(DONE) else "PENDING")
        print(json.dumps(avail, indent=2))
        return

    # Lightweight worker path — used under subprocess to bound RSS
    if args.one_seed:
        require_cuda()
        direction, src_batch, seed_s = args.one_seed
        seed = int(seed_s)
        tgt_batch = "PBMC5k" if src_batch == "PBMC10k" else "PBMC10k"
        print(f"WORKER {direction} seed{seed}", flush=True)
        run_one_seed(direction, src_batch, tgt_batch, seed, smoke_only=False)
        gc.collect()
        try:
            import torch

            torch.cuda.empty_cache()
        except Exception:
            pass
        return

    write_gpu_env()
    write_server_inventory(avail)
    require_cuda()

    if avail["n_found"] < 80:
        print("WARNING: incomplete adapted model grid", avail)
        if avail["n_found"] == 0:
            print("BRANCH B required — dedicated training not auto-launched in this entrypoint.")
            sys.exit(2)

    if args.dry_run:
        print("dry-run OK; branch", avail["branch"], "found", avail["n_found"])
        return

    if args.aggregate_only:
        rows, pc, integ = collect_all_cached()
        aggregate_and_write(rows, pc, integ)
        print("aggregated", len(rows))
        return

    if args.smoke:
        smoke = smoke_test()
        print(json.dumps(smoke, indent=2))
        if not args.run:
            return

    seeds = gate_seeds(args.gate)
    print(f"Running gate {args.gate} seeds={seeds} (subprocess-per-seed for RAM)", flush=True)

    import subprocess

    script = str(Path(__file__).resolve())
    failures = []
    for direction, src, tgt in DIRECTIONS:
        for seed in seeds:
            # skip if both models complete
            if totalvi_unit_complete(direction, seed) and scvi_unit_complete(direction, seed):
                print(f"skip complete {direction} seed{seed}", flush=True)
                continue
            print(f"=== subprocess {direction} seed{seed} ===", flush=True)
            cmd = [
                sys.executable,
                "-u",
                script,
                "--one-seed",
                direction,
                src,
                str(seed),
            ]
            env = dict(**{k: v for k, v in __import__("os").environ.items()})
            env["PYTHONPATH"] = str(ROOT)
            env["MPLBACKEND"] = "Agg"
            rc = subprocess.call(cmd, env=env)
            if rc != 0:
                failures.append((direction, seed, rc))
                print(f"FAILED {direction} seed{seed} rc={rc}", flush=True)
            gc.collect()

    rows, pc, integ = collect_all_cached()
    if not rows:
        raise RuntimeError("No cached rows after gate run")
    aggregate_and_write(rows, pc, integ)

    clean = pd.DataFrame(rows)
    if "is_clean_baseline" in clean.columns:
        clean["is_clean_baseline"] = clean["is_clean_baseline"].astype(str).str.lower().isin(["true", "1"])
        clean = clean[clean["is_clean_baseline"]]
    bullets = [
        "- Design: clean-adapted PHASE 11B models → freeze → corrupt TARGET protein only → inference.",
        f"- Adapted models: PHASE 11B (`{MODELS_ROOT}`), {avail['n_found']}/80.",
        f"- Seeds requested this gate: {seeds}",
        f"- Subprocess failures: {failures}",
    ]
    for d in ("direction_A", "direction_B"):
        for m in MODEL_NAMES:
            g = clean[(clean["direction"] == d) & (clean["model"] == m)] if len(clean) else clean
            if len(g):
                bullets.append(f"- Clean {d} {m} macro-F1 mean={g['macro_f1'].mean():.4f} (n={len(g)})")
    state = "STATE 2"
    append_synthesis(DONE, state, bullets)
    write_done(
        DONE,
        {
            "ok": len(failures) == 0,
            "gate": args.gate,
            "seeds": seeds,
            "n_rows": len(rows),
            "failures": failures,
            "state": state,
            "ts": _utc(),
        },
    )
    print("Part C gate complete:", args.gate, "rows", len(rows), "failures", failures)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
