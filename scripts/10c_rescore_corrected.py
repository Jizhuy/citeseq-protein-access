#!/usr/bin/env python
"""Re-score PHASE 10 from saved embeddings using post-adaptation source latents.

Background
----------
`run_one_model` returned ``src["latent"]`` (the PRE-adaptation source embedding)
for evaluation while writing the POST-adaptation embedding to
``*_source_latent.npy``. Freshly-trained runs were therefore scored with mixed
coordinate systems, violating the PHASE 10 latent-alignment rule. Resumed runs
read the saved file and were correct.

Effect: scVI is unaffected (scArches freezes the reference embedding, drift
exactly 0.0). totalVI seeds 1-4 are affected (drift ~0.015-0.020).

This script retrains nothing. It reloads the saved post-adaptation embeddings
and recomputes the affected metrics with PHASE 10's own functions, writing to a
NEW directory. No historical PHASE 10 artifact is modified.

Columns carried over unchanged from the original tables are those computed from
the target latent alone (target_asw, target_knn_purity_15, target_leiden_ari,
target_leiden_nmi); the bug could not touch them.
"""

from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.evaluation.representation_metrics import batch_silhouette  # noqa: E402
from src.experiments.phase10_cross_dataset import (  # noqa: E402
    DIRECTIONS,
    KNN_KS,
    WITHIN_L1_HIGH,
    WITHIN_L2_HIGH,
    WITHIN_NOTE,
    _phase10_root,
    class_plan,
    classify_transfer,
    load_source_target,
    neighbor_agreement,
    paired_bootstrap_delta,
    per_class_table,
)
from src.utils.io import save_json  # noqa: E402
from sklearn.metrics import confusion_matrix  # noqa: E402

MODELS = ("scvi_matched", "totalvi")
SEEDS = (0, 1, 2, 3, 4)
CARRY_OVER = ("target_asw", "target_knn_purity_15", "target_leiden_ari", "target_leiden_nmi")


def load_saved(direction: str, model: str, seed: int) -> tuple[np.ndarray, np.ndarray]:
    emb = _phase10_root() / "embeddings" / direction
    tag = f"{model}_{direction}_seed{seed}_cuda"
    return (
        np.load(emb / f"{tag}_source_latent.npy"),
        np.load(emb / f"{tag}_target_latent.npy"),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-dataset-asw", action="store_true",
                        help="skip the O(n^2) silhouette recomputation")
    args = parser.parse_args()

    root = _phase10_root()
    out_dir = root / "tables_corrected"
    out_dir.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(root / "tables" / "cross_dataset_primary.csv")

    primary_rows, delta_rows, cell_frames, neigh_rows, gap_rows, conf_frames = [], [], [], [], [], []
    changed = []

    for direction, source_batch, target_batch in DIRECTIONS:
        print(f"\n=== {direction} ===", flush=True)
        source, target, _ = load_source_target(source_batch, target_batch)
        plans = {
            "l2": class_plan(source, target, "cell_type_l2"),
            "l1": class_plan(source, target, "cell_type_l1"),
        }
        for seed in SEEDS:
            packs = {}
            for model in MODELS:
                z_s, z_t = load_saved(direction, model, seed)
                packs[model] = (z_s, z_t)

            dataset_asw = {}
            if not args.skip_dataset_asw:
                for model, (z_s, z_t) in packs.items():
                    joint = np.vstack([z_s, z_t])
                    ids = np.array(["source"] * z_s.shape[0] + ["target"] * z_t.shape[0])
                    dataset_asw[model] = float(batch_silhouette(joint, ids))

            for level, col in (("l2", "cell_type_l2"), ("l1", "cell_type_l1")):
                plan = plans[level]
                si, ti = plan["source_eval_index"], plan["target_eval_index"]
                classes = plan["eligible_shared_classes"]
                y_src = source.obs[col].astype(str).to_numpy()[si]
                y_tgt = target.obs[col].astype(str).to_numpy()[ti]
                src_counts = pd.Series(y_src).value_counts()
                tgt_counts = pd.Series(y_tgt).value_counts()
                preds = {}

                for model in MODELS:
                    z_s, z_t = packs[model]
                    zs, zt = z_s[si], z_t[ti]
                    clf = classify_transfer(zs, y_src, zt, y_tgt, seed)
                    preds[model] = clf

                    base_i = len(primary_rows)
                    for classifier in ("logreg", "knn15", "knn30"):
                        primary_rows.append(
                            {
                                "direction": direction,
                                "model": model,
                                "model_seed": seed,
                                "source_n": int(source.n_obs),
                                "target_n": int(target.n_obs),
                                "label_level": level,
                                "n_shared_classes": plan["n_shared_classes"],
                                "target_coverage": plan["target_coverage"],
                                "classifier": classifier,
                                "accuracy": clf[f"{classifier}_accuracy"],
                                "balanced_accuracy": clf[f"{classifier}_balanced_accuracy"],
                                "macro_f1": clf[f"{classifier}_macro_f1"],
                                "weighted_f1": clf[f"{classifier}_weighted_f1"],
                                "n_eval_source": plan["n_source_in_eval"],
                                "n_eval_target": plan["n_target_in_eval"],
                                "source_embedding_used": "adapted_model_post_adaptation",
                            }
                        )
                    orig = original[
                        (original.direction == direction)
                        & (original.model == model)
                        & (original.model_seed == seed)
                        & (original.label_level == level)
                        & (original.classifier == "logreg")
                    ]
                    for c in CARRY_OVER:
                        if c in orig.columns and len(orig):
                            primary_rows[base_i][c] = orig[c].iloc[0]
                    if model in dataset_asw:
                        primary_rows[base_i]["dataset_asw"] = dataset_asw[model]
                    primary_rows[base_i]["latent_spaces_aligned_by"] = (
                        "scArches: source and target both encoded by the post-adaptation model"
                    )

                    if len(orig):
                        old = float(orig["macro_f1"].iloc[0])
                        new = float(clf["logreg_macro_f1"])
                        if abs(new - old) > 1e-9:
                            changed.append(
                                {
                                    "direction": direction, "model": model, "model_seed": seed,
                                    "label_level": level, "original_macro_f1": old,
                                    "corrected_macro_f1": new, "delta": new - old,
                                }
                            )

                    spec = per_class_table(y_tgt, clf["logreg_pred"], classes, src_counts, tgt_counts)
                    spec["direction"], spec["model"], spec["seed"] = direction, model, seed
                    spec["label_level"], spec["subset"] = level, "high_confidence_shared"
                    cell_frames.append(spec)

                    mat = confusion_matrix(y_tgt, clf["logreg_pred"], labels=classes, normalize="true")
                    long = (
                        pd.DataFrame(mat, index=classes, columns=classes)
                        .reset_index()
                        .melt(id_vars="index", var_name="predicted", value_name="normalized")
                        .rename(columns={"index": "true"})
                    )
                    long["direction"], long["model"] = direction, model
                    long["seed"], long["label_level"] = seed, level
                    conf_frames.append(long)

                    for k in KNN_KS:
                        agr = neighbor_agreement(zs, y_src, zt, y_tgt, k)
                        neigh_rows.append(
                            {
                                "direction": direction, "model": model, "seed": seed,
                                "label_level": level, "k": k,
                                "mean_agreement": agr["mean_agreement"],
                                "median_agreement": agr["median_agreement"],
                            }
                        )

                    within = (WITHIN_L2_HIGH if level == "l2" else WITHIN_L1_HIGH)[
                        "scVI_matched" if model == "scvi_matched" else "totalVI"
                    ]
                    gap_rows.append(
                        {
                            "model": model, "direction": direction, "label_level": level,
                            "within_dataset_macro_f1": within,
                            "cross_dataset_macro_f1": clf["logreg_macro_f1"],
                            "generalization_gap": float(clf["logreg_macro_f1"] - within),
                            "within_dataset_reference": WITHIN_NOTE, "model_seed": seed,
                        }
                    )

                boot = paired_bootstrap_delta(
                    y_tgt, preds["scvi_matched"]["logreg_pred"], preds["totalvi"]["logreg_pred"],
                    seed=seed,
                )
                delta_rows.append(
                    {
                        "direction": direction, "model_seed": seed, "label_level": level,
                        "scvi_matched_macro_f1": preds["scvi_matched"]["logreg_macro_f1"],
                        "totalvi_macro_f1": preds["totalvi"]["logreg_macro_f1"],
                        "delta_totalvi_minus_scvi": float(
                            preds["totalvi"]["logreg_macro_f1"]
                            - preds["scvi_matched"]["logreg_macro_f1"]
                        ),
                        "bootstrap_ci_low": boot["bootstrap_ci_low"],
                        "bootstrap_ci_high": boot["bootstrap_ci_high"],
                        "bootstrap_delta_mean": boot["delta_mean"],
                        "target_n": plan["n_target_in_eval"],
                        "shared_class_n": plan["n_shared_classes"],
                        "target_coverage": plan["target_coverage"],
                        "classifier": "logreg",
                    }
                )
            print(f"  seed{seed} done", flush=True)
            del packs
            gc.collect()
        del source, target
        gc.collect()

    primary = pd.DataFrame(primary_rows)
    delta = pd.DataFrame(delta_rows)
    primary.to_csv(out_dir / "cross_dataset_primary_corrected.csv", index=False)
    delta.to_csv(out_dir / "cross_dataset_primary_delta_corrected.csv", index=False)
    pd.concat(cell_frames, ignore_index=True).to_csv(
        out_dir / "cross_dataset_celltype_specific_corrected.csv", index=False)
    pd.DataFrame(neigh_rows).to_csv(
        out_dir / "cross_dataset_neighbor_agreement_corrected.csv", index=False)
    pd.DataFrame(gap_rows).to_csv(
        out_dir / "cross_dataset_generalization_gap_corrected.csv", index=False)
    pd.concat(conf_frames, ignore_index=True).to_csv(
        out_dir / "cross_dataset_confusion_corrected.csv", index=False)

    summary = (
        delta.groupby(["direction", "label_level"], as_index=False)
        .agg(
            n_seeds=("model_seed", "nunique"),
            scvi_matched_macro_f1_mean=("scvi_matched_macro_f1", "mean"),
            scvi_matched_macro_f1_sd=("scvi_matched_macro_f1", "std"),
            totalvi_macro_f1_mean=("totalvi_macro_f1", "mean"),
            totalvi_macro_f1_sd=("totalvi_macro_f1", "std"),
            delta_mean=("delta_totalvi_minus_scvi", "mean"),
            delta_sd=("delta_totalvi_minus_scvi", "std"),
        )
    )
    summary.to_csv(out_dir / "cross_dataset_primary_delta_summary_corrected.csv", index=False)
    pd.DataFrame(changed).to_csv(out_dir / "correction_changed_rows.csv", index=False)
    save_json(
        {
            "n_rows_changed": len(changed),
            "models_affected": sorted({c["model"] for c in changed}),
            "max_abs_delta": float(max((abs(c["delta"]) for c in changed), default=0.0)),
            "note": "scVI unaffected: scArches freezes the reference embedding (drift 0.0).",
            "historical_tables_modified": False,
            "retraining_performed": False,
        },
        out_dir / "correction_summary.json",
    )
    print(f"\nrows changed: {len(changed)}")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
