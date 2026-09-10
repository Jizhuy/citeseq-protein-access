#!/usr/bin/env python
"""PHASE 11 stages: source-subset instability (§31-33) and source latent posterior.

subset  Uses the five PHASE 10B PBMC10k source subsamples (Direction A only) to
        measure how much each PBMC5k target prediction moves when the source
        training population is resampled, then tests whether the full-source
        PHASE 10 model already knew which cells were fragile.

latent  Extracts q(z|x) variance from the saved source reference models. Valid
        in the PHASE 10 coordinate system for scVI only (drift 0.0); for totalVI
        the reference model defines the pre-adaptation space. Target-cell
        posterior is unsupported - see posterior_api_validation.md.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import accuracy_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import (  # noqa: E402
    DIRECTIONS,
    MODELS,
    SEEDS,
    align_probability_matrices,
    classification_metrics,
    ensemble_uncertainty,
    ensure_dirs,
    entropy_of,
    error_detection,
    fit_transfer_classifier,
    predictive_uncertainty,
    utc_now,
)
from src.utils.io import resolve_path, save_json  # noqa: E402

DIRECTION_A = "pbmc10k_to_pbmc5k"
MODEL_SEED = 0
SUBSAMPLE_SEEDS = (0, 1, 2, 3, 4)
LEVEL = "l2"


def _p10b() -> Path:
    return resolve_path("results/phase10b_source_size")


def load_common_classes() -> list[str]:
    audit = json.loads((_p10b() / "logs" / "phase10b_design_audit.json").read_text())
    return list(audit["class_sets"][LEVEL]["common_classes"])


def run_subset() -> None:
    from src.experiments.phase10_cross_dataset import load_source_target

    paths = ensure_dirs()
    classes_common = load_common_classes()
    print(f"PHASE 10B common l2 classes: {len(classes_common)}", flush=True)

    source_full, target, _ = load_source_target("PBMC10k", "PBMC5k")
    s_obs = source_full.obs
    t_obs = target.obs
    src_lookup = {c: i for i, c in enumerate(source_full.obs_names.astype(str))}
    high_t = t_obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    y_t_all = t_obs["cell_type_l2"].astype(str).to_numpy()
    t_keep = np.where(high_t & np.isin(y_t_all, classes_common))[0]
    y_tgt = y_t_all[t_keep]
    tgt_ids = target.obs_names.astype(str).to_numpy()[t_keep]
    print(f"PHASE 10B target eval cells: {t_keep.size}", flush=True)

    high_s = s_obs["annotation_tier_l2"].astype(str).to_numpy() == "high"
    y_s_all = s_obs["cell_type_l2"].astype(str).to_numpy()

    emb = _p10b() / "embeddings"
    rows_out = []
    per_model_subset: dict[str, dict] = {}

    for model in MODELS:
        per_subset = {}
        for ss in SUBSAMPLE_SEEDS:
            tag = f"{model}_pbmc10k_sub3994_subsample_seed{ss}_modelseed{MODEL_SEED}_cuda"
            z_s = np.load(emb / f"{tag}_source_latent.npy")
            z_t = np.load(emb / f"{tag}_target_latent.npy")
            s_cells = pd.read_csv(emb / f"{tag}_source_cells.csv")["cell_id"].astype(str).to_numpy()
            t_cells = pd.read_csv(emb / f"{tag}_target_cells.csv")["cell_id"].astype(str).to_numpy()
            if not np.array_equal(t_cells, target.obs_names.astype(str).to_numpy()):
                raise RuntimeError(f"PHASE 10B target cell order differs for {tag}")
            rows = np.array([src_lookup[c] for c in s_cells])
            keep = high_s[rows] & np.isin(y_s_all[rows], classes_common)
            zs = z_s[keep]
            ys = y_s_all[rows][keep]
            fit = fit_transfer_classifier(zs, ys, z_t[t_keep], MODEL_SEED)
            per_subset[ss] = {"probs": fit["probs"], "classes": fit["classes"], "pred": fit["pred"]}
            print(
                f"  {model} subsample{ss}: n_src={zs.shape[0]} "
                f"acc={accuracy_score(y_tgt, fit['pred']):.4f}",
                flush=True,
            )

        stacked, classes = align_probability_matrices(per_subset)
        p_bar = stacked.mean(axis=0)
        member_pred = np.argmax(stacked, axis=2)
        n_sets, n_cells = member_pred.shape
        variation = np.array(
            [1.0 - np.bincount(member_pred[:, i], minlength=classes.size).max() / n_sets
             for i in range(n_cells)]
        )
        subset_entropy = entropy_of(p_bar)
        cls_index = {c: i for i, c in enumerate(classes)}
        true_col = np.array([cls_index.get(c, -1) for c in y_tgt])
        valid = true_col >= 0
        true_probs = np.full((n_sets, n_cells), np.nan)
        for s in range(n_sets):
            true_probs[s, valid] = stacked[s, np.arange(n_cells)[valid], true_col[valid]]
        true_sd = np.nanstd(true_probs, axis=0, ddof=0)
        dispersion = np.mean(np.std(stacked, axis=0, ddof=0), axis=1)

        per_model_subset[model] = {
            "variation": variation,
            "entropy": subset_entropy,
            "true_sd": true_sd,
            "dispersion": dispersion,
            "cell_ids": tgt_ids,
            "y": y_tgt,
        }
        print(
            f"  [{model}] median variation_ratio={np.median(variation):.4f} "
            f"median true-prob SD={np.median(true_sd):.4f}",
            flush=True,
        )

    # ---- full-source PHASE 10 Direction A uncertainty for the same cells ----
    assoc_rows = []
    for model in MODELS:
        seed_packs, ids_full, y_full = {}, None, None
        for seed in SEEDS:
            f = np.load(
                paths["arrays"] / f"probs_{DIRECTION_A}_{model}_seed{seed}_{LEVEL}.npz",
                allow_pickle=False,
            )
            seed_packs[seed] = {"probs": f["probs"].astype(float), "classes": f["classes"]}
            if ids_full is None:
                ids_full = f["cell_ids"].astype(str)
                y_full = f["y_true"].astype(str)
                pred0 = f["pred"].astype(str)
        stacked_full, classes_full = align_probability_matrices(seed_packs)
        ens_full = ensemble_uncertainty(stacked_full, classes_full)
        u_full = predictive_uncertainty(seed_packs[0]["probs"])["u_pred"]
        correct_full = (pred0 == y_full).astype(int)

        sub = per_model_subset[model]
        common_ids = np.intersect1d(ids_full, sub["cell_ids"])
        fi = {c: i for i, c in enumerate(ids_full)}
        si = {c: i for i, c in enumerate(sub["cell_ids"])}
        fidx = np.array([fi[c] for c in common_ids])
        sidx = np.array([si[c] for c in common_ids])
        print(f"  [{model}] cells common to PHASE 10 A and PHASE 10B: {common_ids.size}", flush=True)

        u = u_full[fidx]
        dis = ens_full["disagreement"][fidx]
        vr = sub["variation"][sidx]
        sd = sub["true_sd"][sidx]
        ent = sub["entropy"][sidx]

        for name, x in (("full_source_u_pred", u), ("full_source_ensemble_disagreement", dis)):
            for tname, t in (("subset_variation_ratio", vr),
                             ("subset_true_class_probability_sd", sd),
                             ("subset_probability_entropy", ent)):
                ok = np.isfinite(x) & np.isfinite(t)
                r = spearmanr(x[ok], t[ok])
                assoc_rows.append(
                    {
                        "model": model, "x": name, "y": tname,
                        "spearman_rho": float(r.statistic), "spearman_p": float(r.pvalue),
                        "n": int(ok.sum()),
                    }
                )

        order = np.argsort(u, kind="mergesort")
        k = max(1, int(round(0.20 * u.size)))
        top, bottom = order[-k:], order[:k]
        assoc_rows.append(
            {
                "model": model, "x": "stratified_top20_vs_bottom20_full_source_u_pred",
                "y": "subset_instability",
                "mean_variation_ratio_top20": float(vr[top].mean()),
                "mean_variation_ratio_bottom20": float(vr[bottom].mean()),
                "mean_true_prob_sd_top20": float(sd[top].mean()),
                "mean_true_prob_sd_bottom20": float(sd[bottom].mean()),
                "variation_ratio_enrichment": float(
                    vr[top].mean() / vr[bottom].mean()) if vr[bottom].mean() > 0 else float("inf"),
                "n": int(k),
            }
        )

        for j, cid in enumerate(common_ids):
            rows_out.append(
                {
                    "model": model,
                    "cell_id": cid,
                    "true_cell_type": sub["y"][sidx[j]],
                    "full_source_correct": int(correct_full[fidx[j]]),
                    "full_source_predicted_label": pred0[fidx[j]],
                    "full_source_predictive_uncertainty": float(u[j]),
                    "full_source_latent_uncertainty": np.nan,
                    "full_source_ensemble_disagreement": float(dis[j]),
                    "subset_prediction_variation_ratio": float(vr[j]),
                    "subset_probability_entropy": float(ent[j]),
                    "subset_true_class_probability_sd": float(sd[j]),
                }
            )

    pd.DataFrame(rows_out).to_csv(
        paths["tables"] / "uncertainty_source_subset_instability.csv", index=False)
    pd.DataFrame(assoc_rows).to_csv(
        paths["tables"] / "uncertainty_subset_associations.csv", index=False)
    print("\nstage subset complete", flush=True)


def run_latent() -> None:
    """Source-cell posterior variance from the saved reference models."""
    import torch
    from scvi.model import SCVI, TOTALVI

    from src.experiments.phase10_cross_dataset import (
        _phase10_root,
        ensure_csr,
        load_source_target,
    )

    paths = ensure_dirs()
    rows = []
    for direction, source_batch, target_batch in DIRECTIONS:
        source, target, _ = load_source_target(source_batch, target_batch)
        source = ensure_csr(source)
        del target
        gc.collect()
        saved_dir = _phase10_root() / "embeddings" / direction
        for model in MODELS:
            cls = SCVI if model == "scvi_matched" else TOTALVI
            for seed in SEEDS:
                mdir = _phase10_root() / "models" / direction / model / f"seed{seed}_cuda"
                if torch.cuda.is_available():
                    torch.cuda.reset_peak_memory_stats()
                m = cls.load(str(mdir), adata=source, accelerator="gpu", device="auto")
                qzm, qzv = m.get_latent_representation(return_dist=True)
                u_latent = qzv.mean(axis=1)
                tag = f"{model}_{direction}_seed{seed}_cuda"
                post = np.load(saved_dir / f"{tag}_source_latent.npy")
                pre = np.load(saved_dir / f"{tag}_source_latent_source_model.npy")
                rows.append(
                    {
                        "direction": direction, "model": model, "model_seed": seed,
                        "n_source_cells": int(qzv.shape[0]), "n_latent": int(qzv.shape[1]),
                        "mean_posterior_variance": float(qzv.mean()),
                        "median_u_latent_source": float(np.median(u_latent)),
                        "iqr_u_latent_source": float(
                            np.percentile(u_latent, 75) - np.percentile(u_latent, 25)),
                        "min_posterior_variance": float(qzv.min()),
                        "max_posterior_variance": float(qzv.max()),
                        "gaussian_entropy_mean": float(
                            np.mean(0.5 * np.sum(np.log(2 * np.pi * np.e * qzv), axis=1))),
                        "qzm_matches_saved_post_adaptation": bool(np.allclose(qzm, post, atol=1e-4)),
                        "qzm_matches_saved_pre_adaptation": bool(np.allclose(qzm, pre, atol=1e-4)),
                        "valid_in_phase10_coordinate_system": bool(np.allclose(qzm, post, atol=1e-4)),
                        "cuda_peak_allocated_mb": (
                            float(torch.cuda.max_memory_allocated()) / 1024**2
                            if torch.cuda.is_available() else np.nan),
                    }
                )
                np.save(paths["arrays"] / f"u_latent_source_{tag}.npy", u_latent.astype(np.float32))
                print(
                    f"  {direction} {model} seed{seed}: mean qzv={qzv.mean():.5f} "
                    f"valid_post={rows[-1]['valid_in_phase10_coordinate_system']}",
                    flush=True,
                )
                del m
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        del source
        gc.collect()

    frame = pd.DataFrame(rows)
    frame.to_csv(paths["tables"] / "uncertainty_source_latent_posterior.csv", index=False)
    save_json(
        {
            "target_latent_posterior": "UNSUPPORTED - adapted query model not persisted",
            "source_latent_posterior": "extracted via get_latent_representation(return_dist=True)",
            "valid_rows": int(frame["valid_in_phase10_coordinate_system"].sum()),
            "total_rows": int(len(frame)),
            "note": (
                "scVI source rows are exactly valid in the PHASE 10 coordinate system "
                "(scArches drift 0.0). totalVI reference models define the pre-adaptation "
                "space, so their source posterior is not in the PHASE 11 coordinate system."
            ),
            "timestamp_utc": utc_now(),
        },
        paths["logs"] / "latent_posterior_summary.json",
    )
    print("\nstage latent complete", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=["subset", "latent"])
    args = parser.parse_args()
    if args.stage == "subset":
        run_subset()
    else:
        run_latent()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
