#!/usr/bin/env python3
"""PART 3: Lawlor donor-level metrics from saved cell-level predictions.

Primary macro-F1 definition (chosen BEFORE inspecting model deltas):
  PRESENT-CLASS MACRO-F1 — average F1 over classes with >0 labeled support
  in that donor. Rationale: fixed-ontology zeros for absent classes are an
  arbitrary scoring convention, not a biological observation.

Sensitivity:
  FIXED-ONTOLOGY MACRO-F1 — average over the predeclared 7 eligible classes
  with sklearn zero_division=0 for both models.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[3]
PHASE12B = ROOT / "results" / "phase12b_external_validation" / "formal_validation"
OUT = ROOT / "results" / "revision_round2_methodological_fixes"
STATS = OUT / "statistics"
FIGS = OUT / "figures"
LOGS = OUT / "logs"

ELIGIBLE = [
    "B",
    "CD14_Mono",
    "CD4T_Mem",
    "CD4T_Naive",
    "CD8T_Mem",
    "CD8T_Naive",
    "NK",
]
ECE_BINS = 10
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)

PRIMARY_MACRO_F1 = "present_class_macro_f1"


def expected_calibration_error(confidence, correct, bins, scheme):
    n = confidence.size
    if scheme == "equal_frequency":
        quantiles = np.quantile(confidence, np.linspace(0.0, 1.0, bins + 1))
        edges = np.unique(quantiles)
    else:
        edges = np.linspace(0.0, 1.0, bins + 1)
    assignment = np.clip(np.digitize(confidence, edges[1:-1], right=True), 0, len(edges) - 2)
    ece = 0.0
    for b in range(len(edges) - 1):
        mask = assignment == b
        count = int(mask.sum())
        if count == 0:
            continue
        ece += (count / n) * abs(float(confidence[mask].mean()) - float(correct[mask].mean()))
    return float(ece)


def risk_coverage(uncertainty, correct):
    order = np.argsort(uncertainty, kind="mergesort")
    kept = correct[order].astype(float)
    n = correct.size
    counts = np.arange(1, n + 1)
    risks = 1.0 - np.cumsum(kept) / counts
    coverage = counts / n
    full_aurc = float(np.mean(risks))
    # partial AURC over coverage in [0.5, 1.0]: mean risk among prefixes with cov>=0.5
    mask = coverage >= 0.5 - 1e-12
    paurc = float(np.mean(risks[mask]))
    # E-AURC (Geifman & El-Yaniv style): full AURC minus AURC of perfect ranking
    # Perfect: all correct first, then incorrect.
    n_correct = int(correct.sum())
    n_err = n - n_correct
    if n == 0:
        e_aurc = float("nan")
        aurc_star = float("nan")
    elif n_err == 0:
        aurc_star = 0.0
        e_aurc = full_aurc - aurc_star
    else:
        # risks for perfect order: 0 for first n_correct, then (k-n_correct)/k for rest
        perfect_risks = np.zeros(n, dtype=float)
        for k in range(n_correct + 1, n + 1):
            perfect_risks[k - 1] = (k - n_correct) / k
        aurc_star = float(np.mean(perfect_risks))
        e_aurc = full_aurc - aurc_star
    risk_at = {}
    for cov in COVERAGE_GRID:
        keep = max(1, int(round(cov * n)))
        risk_at[cov] = float(1.0 - accuracy_score(correct.astype(bool)[order[:keep]], np.ones(keep, dtype=bool))) if False else float(1.0 - kept[:keep].mean())
    return {
        "full_aurc": full_aurc,
        "paurc_0.5_1.0": paurc,
        "e_aurc": float(e_aurc),
        "aurc_star": float(aurc_star),
        "risk_at": risk_at,
    }


def macro_f1_present(y_true, y_pred):
    present = sorted(set(y_true))
    if not present:
        return float("nan"), present
    return float(f1_score(y_true, y_pred, labels=present, average="macro", zero_division=0)), present


def macro_f1_fixed(y_true, y_pred, labels=ELIGIBLE):
    return float(f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0))


def safe_auroc(y, score):
    y = np.asarray(y)
    if y.size < 2 or len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, score))


def safe_auprc(y, score):
    y = np.asarray(y)
    if y.size < 2 or len(np.unique(y)) < 2:
        return float("nan")
    return float(average_precision_score(y, score))


def load_donor_map():
    ann = pd.read_csv(
        ROOT / "results/phase12b_external_validation/raw_audit/baseline_sng_annotations.csv"
    )
    counts = pd.read_csv(
        ROOT / "results/phase12b_external_validation/tables/external_donor_condition_counts.csv"
    )
    chip_to_donor = (
        counts[["donor", "Donor_of_Origin"]]
        .drop_duplicates()
        .set_index("Donor_of_Origin")["donor"]
        .to_dict()
    )
    ann = ann.copy()
    ann["donor"] = ann["Donor_of_Origin"].map(chip_to_donor)
    return ann.set_index("barcode")["donor"].to_dict()


def metrics_for_subset(y_true, y_pred, probs, classes, u_pred, u_latent):
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    classes = np.asarray(classes).astype(str)
    n = y_true.size
    if n == 0:
        return None
    error = (y_pred != y_true).astype(int)
    present_f1, present = macro_f1_present(y_true, y_pred)
    fixed_f1 = macro_f1_fixed(y_true, y_pred)
    # NLL / Brier
    class_to_i = {c: i for i, c in enumerate(classes)}
    # classes array in npz may be model class order
    true_idx = np.array([class_to_i.get(y, -1) for y in y_true], dtype=int)
    if (true_idx < 0).any():
        # drop unknowns (should not happen)
        keep = true_idx >= 0
        y_true, y_pred, probs, u_pred, u_latent, error, true_idx = (
            y_true[keep],
            y_pred[keep],
            probs[keep],
            u_pred[keep],
            u_latent[keep],
            error[keep],
            true_idx[keep],
        )
        n = y_true.size
        if n == 0:
            return None
    p_true = probs[np.arange(n), true_idx]
    nll = float(-np.mean(np.log(np.clip(p_true, 1e-12, 1.0))))
    Y = np.zeros_like(probs)
    Y[np.arange(n), true_idx] = 1.0
    brier = float(np.mean(np.sum((probs - Y) ** 2, axis=1)))
    conf = probs.max(axis=1)
    correct = (error == 0).astype(float)
    ece_f = expected_calibration_error(conf, correct, ECE_BINS, "equal_frequency")
    ece_w = expected_calibration_error(conf, correct, ECE_BINS, "equal_width")
    rc = risk_coverage(u_pred, correct)
    return {
        "n_cells": int(n),
        "n_present_classes": int(len(present)),
        "present_classes": ";".join(present),
        "present_class_macro_f1": present_f1,
        "fixed_ontology_macro_f1": fixed_f1,
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "predictive_error_auroc": safe_auroc(error, u_pred),
        "predictive_error_auprc": safe_auprc(error, u_pred),
        "latent_error_auroc": safe_auroc(error, u_latent),
        "latent_error_auprc": safe_auprc(error, u_latent),
        "nll": nll,
        "brier": brier,
        "ece_equal_frequency": ece_f,
        "ece_equal_width": ece_w,
        "full_aurc": rc["full_aurc"],
        "paurc_0.5_1.0": rc["paurc_0.5_1.0"],
        "e_aurc": rc["e_aurc"],
        "risk_100": rc["risk_at"][1.0],
        "risk_90": rc["risk_at"][0.9],
        "risk_80": rc["risk_at"][0.8],
        "risk_70": rc["risk_at"][0.7],
        "risk_60": rc["risk_at"][0.6],
        "risk_50": rc["risk_at"][0.5],
        "error_prevalence": float(error.mean()),
    }


def collect_per_donor_rows(donor_map):
    rows = []
    for fold in range(5):
        emb_dir = PHASE12B / "embeddings" / f"fold{fold}"
        prob_dir = PHASE12B / "probabilities" / f"fold{fold}"
        # cell list identical across models/seeds within fold
        target_cells = pd.read_csv(emb_dir / f"fold{fold}__scvi_matched__seed00_target_cells.csv")
        donors = target_cells["cell_id"].map(donor_map)
        if donors.isna().any():
            raise RuntimeError(f"Unmapped target cells in fold {fold}")
        for model in ("scvi_matched", "totalvi"):
            for seed in range(10):
                tag = f"fold{fold}__{model}__seed{seed:02d}"
                path = prob_dir / f"{tag}_probs.npz"
                d = np.load(path, allow_pickle=True)
                ti = d["target_eval_index"].astype(int)
                y_true = d["y_true"].astype(str)
                y_pred = d["y_pred"].astype(str)
                probs = d["probs"]
                classes = d["classes"].astype(str)
                u_pred = d["u_pred"]
                u_latent = d["u_latent"]
                eval_donors = donors.iloc[ti].to_numpy()
                assert len(eval_donors) == len(y_true)
                for donor in sorted(set(eval_donors)):
                    mask = eval_donors == donor
                    m = metrics_for_subset(
                        y_true[mask],
                        y_pred[mask],
                        probs[mask],
                        classes,
                        u_pred[mask],
                        u_latent[mask],
                    )
                    if m is None:
                        continue
                    rows.append(
                        {
                            "fold": fold,
                            "donor": donor,
                            "model": model,
                            "seed": seed,
                            **m,
                        }
                    )
    return pd.DataFrame(rows)


def donor_deltas(per_donor: pd.DataFrame) -> pd.DataFrame:
    """Average over seeds within donor×model, then totalVI - scVI."""
    # Declare primary before delta inspection is already done via PRIMARY_MACRO_F1 constant.
    metric_cols = [
        PRIMARY_MACRO_F1,
        "fixed_ontology_macro_f1",
        "accuracy",
        "balanced_accuracy",
        "predictive_error_auroc",
        "predictive_error_auprc",
        "latent_error_auroc",
        "latent_error_auprc",
        "nll",
        "brier",
        "ece_equal_frequency",
        "ece_equal_width",
        "full_aurc",
        "paurc_0.5_1.0",
        "e_aurc",
        "risk_100",
        "risk_50",
        "error_prevalence",
    ]
    # within-donor seed SD
    seed_sd = (
        per_donor.groupby(["donor", "model"])[metric_cols]
        .std(ddof=1)
        .reset_index()
        .melt(id_vars=["donor", "model"], var_name="metric", value_name="within_donor_seed_sd")
    )
    means = (
        per_donor.groupby(["donor", "model"])[metric_cols]
        .mean()
        .reset_index()
    )
    scvi = means[means["model"] == "scvi_matched"].set_index("donor")
    tot = means[means["model"] == "totalvi"].set_index("donor")
    donors = sorted(scvi.index.intersection(tot.index))
    out_rows = []
    for metric in metric_cols:
        deltas = tot.loc[donors, metric] - scvi.loc[donors, metric]
        for donor in donors:
            out_rows.append(
                {
                    "donor": donor,
                    "metric": metric,
                    "scvi_mean_over_seeds": float(scvi.loc[donor, metric]),
                    "totalvi_mean_over_seeds": float(tot.loc[donor, metric]),
                    "delta": float(deltas.loc[donor]),
                }
            )
        # summary row
        out_rows.append(
            {
                "donor": "SUMMARY",
                "metric": metric,
                "scvi_mean_over_seeds": float(scvi.loc[donors, metric].mean()),
                "totalvi_mean_over_seeds": float(tot.loc[donors, metric].mean()),
                "delta": float(deltas.mean()),
                "median_delta": float(deltas.median()),
                "sd_delta_across_donors": float(deltas.std(ddof=1)),
                "n_donors_positive": int((deltas > 0).sum()),
                "n_donors_negative": int((deltas < 0).sum()),
                "n_donors_equal": int((deltas == 0).sum()),
                "n_donors": int(len(donors)),
            }
        )
    delta_df = pd.DataFrame(out_rows)
    return delta_df, seed_sd, means


def make_figures(means: pd.DataFrame, delta_df: pd.DataFrame) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    metric = PRIMARY_MACRO_F1
    donors = sorted(means["donor"].unique())
    scvi = means[means["model"] == "scvi_matched"].set_index("donor").loc[donors, metric]
    tot = means[means["model"] == "totalvi"].set_index("donor").loc[donors, metric]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    for d in donors:
        ax.plot([0, 1], [scvi.loc[d], tot.loc[d]], "-o", color="0.55", markersize=4)
    ax.scatter(np.zeros(len(donors)), scvi.values, color="#4C72B0", zorder=3, label="scVI_matched")
    ax.scatter(np.ones(len(donors)), tot.values, color="#C44E52", zorder=3, label="totalVI")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["scVI_matched", "totalVI"])
    ax.set_ylabel("Present-class macro-F1 (seed-mean)")
    ax.set_title("Lawlor donor-paired macro-F1")
    ax.legend(frameon=False)
    ax = axes[1]
    dsum = delta_df[(delta_df["metric"] == metric) & (delta_df["donor"] != "SUMMARY")].copy()
    dsum = dsum.sort_values("donor")
    colors = ["#C44E52" if v >= 0 else "#4C72B0" for v in dsum["delta"]]
    ax.bar(dsum["donor"], dsum["delta"], color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("Δ macro-F1 (totalVI − scVI)")
    ax.set_title("Donor-level ΔF1")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(FIGS / "lawlor_donor_paired_metrics.png", dpi=200)
    fig.savefig(FIGS / "lawlor_donor_paired_metrics.pdf")
    plt.close(fig)


def main():
    STATS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    # Freeze primary definition before deltas
    (LOGS / "lawlor_macro_f1_primary_definition.md").write_text(
        f"""# Lawlor donor macro-F1 primary definition

Chosen **before** computing totalVI−scVI donor deltas:

**Primary = `{PRIMARY_MACRO_F1}` (PRESENT-CLASS)**

Reason: classes absent from a donor are not observed biological events for that
donor; assigning them F1=0 under a fixed ontology is a scoring convention that
can dilute or distort donor-level biology. Present-class macro-F1 averages only
over classes with >0 labeled support in that donor.

Sensitivity = `fixed_ontology_macro_f1` over the predeclared 7 eligible classes
with `zero_division=0` for both models.
"""
    )
    donor_map = load_donor_map()
    # integrity: 10 donors
    assert len(set(donor_map.values())) == 10
    per = collect_per_donor_rows(donor_map)
    # expected: 10 donors × 10 seeds × 2 models = 200
    assert per.groupby(["donor", "model", "seed"]).ngroups == 200, per.groupby(["donor", "model", "seed"]).ngroups
    per_path = STATS / "lawlor_per_donor_metrics.csv"
    per.to_csv(per_path, index=False)
    delta_df, seed_sd, means = donor_deltas(per)
    delta_path = STATS / "lawlor_donor_delta.csv"
    delta_df.to_csv(delta_path, index=False)
    seed_sd.to_csv(STATS / "lawlor_within_donor_seed_sd.csv", index=False)
    means.to_csv(STATS / "lawlor_donor_seedmean.csv", index=False)
    make_figures(means, delta_df)
    # concise summary
    prim = delta_df[(delta_df["metric"] == PRIMARY_MACRO_F1) & (delta_df["donor"] == "SUMMARY")].iloc[0]
    print("Primary present-class macro-F1 SUMMARY:")
    print(prim.to_string())
    print("n rows per_donor", len(per))
    print("Wrote", per_path)
    print("Wrote", delta_path)


if __name__ == "__main__":
    main()
