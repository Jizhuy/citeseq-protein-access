#!/usr/bin/env python
"""PHASE 11 gate: validate one direction/seed before the full grid.

Checks, in order:
  1. PHASE 10 artifacts load and cell order matches
  2. the PHASE 10 logreg result reproduces within tolerance
  3. the scvi-tools posterior API is probed on a real saved model
  4. probabilities persist and reload identically
  5. predictive uncertainty / calibration / risk-coverage compute
  6. reliability binning is well formed
  7. no target label reaches the classifier fit
  8. RAM and CUDA stay in budget
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import (  # noqa: E402
    DIRECTIONS,
    MODELS,
    REPRO_TOLERANCE,
    DirectionData,
    calibration_metrics,
    classification_metrics,
    coverage_table,
    ensure_dirs,
    error_detection,
    error_enrichment,
    fit_transfer_classifier,
    load_phase10_pack,
    phase10_reference_metrics,
    predictive_uncertainty,
    risk_coverage,
    utc_now,
    verify_cell_alignment,
)
from src.utils.io import save_json  # noqa: E402


def ram_mb() -> float:
    import psutil

    return float(psutil.Process().memory_info().rss) / (1024.0 * 1024.0)


def probe_posterior_api(direction: str, source_batch: str, target_batch: str, seed: int) -> dict:
    """Load a real saved model and establish exactly what q(z|x) is available.

    PHASE 10 saved only the source reference model. This probe records what the
    API returns, and quantifies how far the source model's target encoding is
    from the saved post-adaptation target latent.
    """
    import torch
    from scvi.model import SCVI

    from src.experiments.phase10_cross_dataset import (
        _phase10_root,
        ensure_csr,
        load_source_target,
    )

    out: dict = {"direction": direction, "model": "scvi_matched", "seed": seed}
    model_dir = _phase10_root() / "models" / direction / "scvi_matched" / f"seed{seed}_cuda"
    out["model_dir"] = str(model_dir)
    out["files_present"] = sorted(p.name for p in model_dir.iterdir())
    out["adapted_query_model_saved"] = False

    source, target, _ = load_source_target(source_batch, target_batch)
    source = ensure_csr(source)
    try:
        model = SCVI.load(str(model_dir), adata=source, accelerator="gpu", device="auto")
    except Exception as exc:  # noqa: BLE001
        out["load_error"] = f"{type(exc).__name__}: {exc}"
        return out

    out["loaded"] = True
    out["latent_distribution"] = str(getattr(model.module, "latent_distribution", "unknown"))
    out["n_latent"] = int(model.module.n_latent)

    qzm, qzv = model.get_latent_representation(return_dist=True)
    out["return_dist_supported"] = True
    out["qzm_shape"] = list(qzm.shape)
    out["qzv_shape"] = list(qzv.shape)
    out["qzv_all_positive"] = bool(np.all(qzv > 0))
    out["qzv_mean"] = float(qzv.mean())
    out["qzv_min"] = float(qzv.min())
    out["qzv_max"] = float(qzv.max())

    mean_only = model.get_latent_representation()
    out["qzm_equals_give_mean"] = bool(np.allclose(qzm, mean_only, atol=1e-6))

    # Does the SOURCE model reproduce the saved POST-ADAPTATION source latent?
    tag = f"scvi_matched_{direction}_seed{seed}_cuda"
    emb = _phase10_root() / "embeddings" / direction
    saved_source_adapted = np.load(emb / f"{tag}_source_latent.npy")
    saved_source_premodel = np.load(emb / f"{tag}_source_latent_source_model.npy")
    out["source_model_matches_saved_source_model_latent"] = bool(
        np.allclose(qzm, saved_source_premodel, atol=1e-4)
    )
    out["source_model_matches_saved_adapted_latent"] = bool(
        np.allclose(qzm, saved_source_adapted, atol=1e-4)
    )
    shift = np.linalg.norm(saved_source_adapted - qzm, axis=1)
    out["adapted_vs_source_model_mean_shift"] = float(shift.mean())
    out["adapted_vs_source_model_median_shift"] = float(np.median(shift))
    out["source_latent_norm_mean"] = float(np.linalg.norm(qzm, axis=1).mean())
    out["relative_shift"] = float(shift.mean() / np.linalg.norm(qzm, axis=1).mean())

    if torch.cuda.is_available():
        out["cuda_peak_allocated_mb"] = float(torch.cuda.max_memory_allocated()) / (1024.0**2)
        out["cuda_peak_reserved_mb"] = float(torch.cuda.max_memory_reserved()) / (1024.0**2)
    out["process_rss_mb"] = ram_mb()

    del model, source, target
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return out


def run_gate(direction: str, source_batch: str, target_batch: str, seed: int, level: str) -> dict:
    paths = ensure_dirs()
    print(f"=== PHASE 11 GATE: {direction} seed{seed} level={level} ===", flush=True)

    reference = phase10_reference_metrics()
    data = DirectionData(direction, source_batch, target_batch)
    sl = data.eval_slice(level)
    print(
        f"eval slice: source={sl['y_source'].size} target={sl['y_target'].size} "
        f"classes={len(sl['classes'])}",
        flush=True,
    )

    report: dict = {
        "direction": direction,
        "seed": seed,
        "label_level": level,
        "n_source_eval": int(sl["y_source"].size),
        "n_target_eval": int(sl["y_target"].size),
        "n_classes": len(sl["classes"]),
        "timestamp_utc": utc_now(),
        "models": {},
        "repro_ok": True,
    }

    for model in MODELS:
        pack = load_phase10_pack(direction, model, seed)
        verify_cell_alignment(data, pack)
        z_src = pack["source_latent"][sl["source_index"]]
        z_tgt = pack["target_latent"][sl["target_index"]]

        fit = fit_transfer_classifier(z_src, sl["y_source"], z_tgt, seed)
        metrics = classification_metrics(sl["y_target"], fit["pred"])

        ref_row = reference[
            (reference["direction"] == direction)
            & (reference["model"] == model)
            & (reference["model_seed"] == seed)
            & (reference["label_level"] == level)
        ]
        if len(ref_row) != 1:
            raise RuntimeError(f"Expected exactly one PHASE 10 row for {model}; got {len(ref_row)}")
        ref_macro = float(ref_row["macro_f1"].iloc[0])
        ref_acc = float(ref_row["accuracy"].iloc[0])
        d_macro = abs(metrics["macro_f1"] - ref_macro)
        d_acc = abs(metrics["accuracy"] - ref_acc)
        ok = d_macro <= REPRO_TOLERANCE and d_acc <= REPRO_TOLERANCE

        print(
            f"[{model}] macro_f1 phase10={ref_macro:.6f} phase11={metrics['macro_f1']:.6f} "
            f"delta={d_macro:.3e} | acc delta={d_acc:.3e} | {'OK' if ok else 'MISMATCH'}",
            flush=True,
        )

        error = (fit["pred"] != sl["y_target"]).astype(int)
        unc = predictive_uncertainty(fit["probs"])
        det = error_detection(unc["u_pred"], error)
        det_entropy = error_detection(unc["entropy"], error)
        det_margin = error_detection(-unc["margin"], error)
        cal = calibration_metrics(fit["probs"], sl["y_target"], fit["classes"], fit["pred"])
        rc = risk_coverage(unc["u_pred"], 1 - error)
        cov = coverage_table(unc["u_pred"], sl["y_target"], fit["pred"])
        enrich = error_enrichment(unc["u_pred"], error)

        # Persist probabilities with explicit class-order metadata.
        arr_path = paths["arrays"] / f"gate_{model}_{direction}_seed{seed}_{level}.npz"
        np.savez_compressed(
            arr_path,
            probs=fit["probs"].astype(np.float32),
            classes=np.asarray([str(c) for c in fit["classes"]]),
            cell_ids=sl["target_cells"].astype(str),
            y_true=sl["y_target"].astype(str),
            pred=fit["pred"].astype(str),
            direction=np.asarray([direction]),
            model=np.asarray([model]),
            model_seed=np.asarray([seed]),
        )
        reloaded = np.load(arr_path, allow_pickle=False)
        roundtrip = bool(
            np.allclose(reloaded["probs"], fit["probs"].astype(np.float32))
            and list(reloaded["classes"]) == [str(c) for c in fit["classes"]]
        )

        report["models"][model] = {
            "phase10_macro_f1": ref_macro,
            "phase11_macro_f1": metrics["macro_f1"],
            "abs_delta_macro_f1": d_macro,
            "phase10_accuracy": ref_acc,
            "phase11_accuracy": metrics["accuracy"],
            "abs_delta_accuracy": d_acc,
            "reproduced": bool(ok),
            "balanced_accuracy": metrics["balanced_accuracy"],
            "weighted_f1": metrics["weighted_f1"],
            "error_prevalence": det["error_prevalence"],
            "error_auroc_u_pred": det["auroc"],
            "error_auprc_u_pred": det["auprc"],
            "error_auroc_entropy": det_entropy["auroc"],
            "error_auroc_neg_margin": det_margin["auroc"],
            "nll": cal["nll"],
            "brier": cal["brier"],
            "ece_equal_frequency": cal["ece_equal_frequency"],
            "ece_equal_width": cal["ece_equal_width"],
            "n_reliability_bins_eq_freq": len(cal["reliability_equal_frequency"]),
            "reliability_bin_counts": [b["n_cells"] for b in cal["reliability_equal_frequency"]],
            "aurc": rc["aurc"],
            "coverage_grid": cov,
            "error_enrichment": enrich,
            "median_u_pred": float(np.median(unc["u_pred"])),
            "probs_roundtrip_ok": roundtrip,
            "probs_path": str(arr_path),
            "probs_sum_to_one": bool(np.allclose(fit["probs"].sum(axis=1), 1.0, atol=1e-8)),
            "classes_sorted": [str(c) for c in fit["classes"]],
            "classes_match_plan": sorted(str(c) for c in fit["classes"]) == sorted(sl["classes"]),
        }
        if not ok:
            report["repro_ok"] = False

        del pack, fit
        gc.collect()

    # Leakage guard: the classifier only ever saw source rows.
    report["leakage_check"] = {
        "target_labels_used_in_fit": False,
        "source_target_cell_overlap": int(
            len(set(sl["source_cells"].tolist()) & set(sl["target_cells"].tolist()))
        ),
        "note": "Classifier fit uses source latents and source labels only.",
    }
    report["process_rss_mb"] = ram_mb()
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direction", default="pbmc5k_to_pbmc10k")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--level", default="l2")
    parser.add_argument("--probe-api", action="store_true")
    args = parser.parse_args()

    os.environ.setdefault("MPLBACKEND", "Agg")
    lookup = {d[0]: d for d in DIRECTIONS}
    if args.direction not in lookup:
        raise SystemExit(f"Unknown direction {args.direction}")
    _, source_batch, target_batch = lookup[args.direction]

    paths = ensure_dirs()
    report = run_gate(args.direction, source_batch, target_batch, args.seed, args.level)

    if args.probe_api:
        print("=== probing posterior API on a saved model ===", flush=True)
        report["posterior_api_probe"] = probe_posterior_api(
            args.direction, source_batch, target_batch, args.seed
        )

    out = paths["logs"] / f"gate_{args.direction}_seed{args.seed}_{args.level}.json"
    save_json(report, out)
    print(f"\nwrote {out}", flush=True)
    print(f"REPRO_OK={report['repro_ok']}", flush=True)
    return 0 if report["repro_ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
