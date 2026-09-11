#!/usr/bin/env python3
"""PART 7–9: statistical-unit ledger, selective-prediction correction, reliability table."""
from __future__ import annotations

import json
from pathlib import Path

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
OUT = ROOT / "results" / "revision_round2_methodological_fixes"
STATS = OUT / "statistics"
CAL = OUT / "calibration"
LOGS = OUT / "logs"
PHASE12B = ROOT / "results/phase12b_external_validation/formal_validation"
PHASE11B = ROOT / "results/phase11b_uncertainty_robustness/tables"
PHASE11 = ROOT / "results/phase11_uncertainty"

ECE_BINS = 10
COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)


def ece(confidence, correct, bins, scheme):
    n = confidence.size
    if scheme == "equal_frequency":
        edges = np.unique(np.quantile(confidence, np.linspace(0, 1, bins + 1)))
    else:
        edges = np.linspace(0, 1, bins + 1)
    assignment = np.clip(np.digitize(confidence, edges[1:-1], right=True), 0, len(edges) - 2)
    val = 0.0
    for b in range(len(edges) - 1):
        m = assignment == b
        c = int(m.sum())
        if c == 0:
            continue
        val += (c / n) * abs(float(confidence[m].mean()) - float(correct[m].mean()))
    return float(val)


def selective_metrics(u_pred, y_true, y_pred, probs, classes):
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
        risk_at[f"risk_{int(cov*100)}"] = float(1.0 - kept[:keep].mean())
    # proper scores
    class_to_i = {c: i for i, c in enumerate(np.asarray(classes).astype(str))}
    true_idx = np.array([class_to_i[y] for y in y_true])
    p_true = probs[np.arange(n), true_idx]
    nll = float(-np.mean(np.log(np.clip(p_true, 1e-12, 1.0))))
    Y = np.zeros_like(probs)
    Y[np.arange(n), true_idx] = 1.0
    brier = float(np.mean(np.sum((probs - Y) ** 2, axis=1)))
    conf = probs.max(axis=1)
    out = {
        "n_cells": n,
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "error_prevalence": float(error.mean()),
        "error_auroc": float(roc_auc_score(error, u_pred)) if len(np.unique(error)) > 1 else float("nan"),
        "error_auprc": float(average_precision_score(error, u_pred)) if len(np.unique(error)) > 1 else float("nan"),
        "nll": nll,
        "brier": brier,
        "ece_equal_frequency": ece(conf, correct, ECE_BINS, "equal_frequency"),
        "ece_equal_width": ece(conf, correct, ECE_BINS, "equal_width"),
        "full_aurc": full_aurc,
        "paurc_0.5_1.0": paurc,
        "e_aurc": e_aurc,
        "aurc_star_perfect_ranking": aurc_star,
        **risk_at,
    }
    return out


def write_statistical_unit_ledger():
    rows = [
        {
            "result": "PBMC HC Delta F1 +0.057 CI",
            "point_estimate": "+0.057",
            "variation_measure": "paired bootstrap CI [0.020,0.091]",
            "unit": "target/test cell",
            "n_units": "n_test cells",
            "resampling_unit": "cell",
            "number_resamples": 500,
            "biological_unit": "cell (within one annotated cohort)",
            "computational_unit": "model seed 0 (primary)",
            "interpretation": "cell-level conditional paired bootstrap CI — NOT biological generalization",
        },
        {
            "result": "PBMC Direction A/B macro-F1 mean±SD (20-seed)",
            "point_estimate": "seed mean",
            "variation_measure": "SD across model seeds",
            "unit": "model seed",
            "n_units": 20,
            "resampling_unit": "seed (optional seed bootstrap 10000 in 11B)",
            "number_resamples": 10000,
            "biological_unit": "cell / dataset batch",
            "computational_unit": "model seed",
            "interpretation": "computational seed variability under unsupervised query adaptation",
        },
        {
            "result": "Lawlor macro-F1 paired fold×seed (historical)",
            "point_estimate": "mean over 50 pairs",
            "variation_measure": "SD over fold×seed",
            "unit": "fold×seed pair",
            "n_units": 50,
            "resampling_unit": "none primary",
            "number_resamples": 0,
            "biological_unit": "donor (but fold pairs two donors)",
            "computational_unit": "seed within fold",
            "interpretation": "NOT 50 independent biological replicates; descriptive only",
        },
        {
            "result": "Lawlor donor-level Delta F1 (revision primary)",
            "point_estimate": "mean of 10 donor deltas",
            "variation_measure": "SD across donors; within-donor seed SD separate",
            "unit": "donor",
            "n_units": 10,
            "resampling_unit": "donor bootstrap (sensitivity only)",
            "number_resamples": "optional",
            "biological_unit": "donor",
            "computational_unit": "seed (averaged within donor)",
            "interpretation": "primary biological generalization summary for Lawlor",
        },
        {
            "result": "Variance components (PHASE12A)",
            "point_estimate": "MoM VC fractions",
            "variation_measure": "component variances",
            "unit": "run (subset×seed×replicate)",
            "n_units": "5×5×2=50 per model",
            "resampling_unit": "N/A (ANOVA MoM)",
            "number_resamples": 0,
            "biological_unit": "source subset (computational subsample)",
            "computational_unit": "seed + run replicate",
            "interpretation": "exploratory; 5 levels per factor — do not overstate precision",
        },
        {
            "result": "Protein corruption deltas",
            "point_estimate": "mean Delta F1 vs scVI",
            "variation_measure": "across corruption-mask seeds",
            "unit": "corruption-mask seed",
            "n_units": 5,
            "resampling_unit": "mask seed (+ cell bootstrap secondary)",
            "number_resamples": 500,
            "biological_unit": "cell",
            "computational_unit": "perturbation seed; model seed fixed 0",
            "interpretation": "training-time protein-entry degradation variability",
        },
        {
            "result": "Selective prediction AURC",
            "point_estimate": "full AURC = mean risk over all prefixes",
            "variation_measure": "across seeds/folds",
            "unit": "inherits parent experiment",
            "n_units": "inherits",
            "resampling_unit": "inherits",
            "number_resamples": "inherits",
            "biological_unit": "inherits",
            "computational_unit": "inherits",
            "interpretation": "historical column `aurc` is FULL AURC; coverage grid tables are partial summaries",
        },
    ]
    path = STATS / "statistical_unit_ledger.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def recompute_lawlor_selective_and_reliability():
    """From saved Lawlor probabilities: full AURC, pAURC, E-AURC, reliability."""
    sel_rows = []
    rel_rows = []
    for fold in range(5):
        for model in ("scvi_matched", "totalvi"):
            for seed in range(10):
                tag = f"fold{fold}__{model}__seed{seed:02d}"
                d = np.load(PHASE12B / "probabilities" / f"fold{fold}" / f"{tag}_probs.npz", allow_pickle=True)
                m = selective_metrics(d["u_pred"], d["y_true"], d["y_pred"], d["probs"], d["classes"])
                # also latent AUROC if available
                err = (d["y_pred"].astype(str) != d["y_true"].astype(str)).astype(int)
                lat_auroc = float(roc_auc_score(err, d["u_latent"])) if len(np.unique(err)) > 1 else float("nan")
                lat_auprc = float(average_precision_score(err, d["u_latent"])) if len(np.unique(err)) > 1 else float("nan")
                row = {"experiment": "Lawlor", "fold": fold, "model": model, "seed": seed, **m, "latent_error_auroc": lat_auroc, "latent_error_auprc": lat_auprc}
                sel_rows.append(row)
                rel_rows.append(row)
    sel = pd.DataFrame(sel_rows)
    # Also seed-mean summaries for Dir A/B from PHASE11B primary if available
    pbmc_rows = []
    primary = PHASE11B / "phase11b_20seed_primary.csv"
    if primary.exists():
        df = pd.read_csv(primary)
        # expected columns vary — keep flexible
        for direction in sorted(df.get("direction", pd.Series(dtype=str)).dropna().unique()) if "direction" in df.columns else []:
            for model in sorted(df["model"].unique()):
                sub = df[(df["direction"] == direction) & (df["model"] == model)]
                # Map historical aurc as full_aurc (code-confirmed)
                rec = {
                    "experiment": f"PBMC_{direction}",
                    "model": model,
                    "n_seeds": int(len(sub)),
                    "macro_f1_mean": float(sub["macro_f1"].mean()) if "macro_f1" in sub else float("nan"),
                    "macro_f1_sd": float(sub["macro_f1"].std(ddof=1)) if "macro_f1" in sub else float("nan"),
                }
                for col, new in [
                    ("predictive_error_auroc", "error_auroc"),
                    ("error_auroc", "error_auroc"),
                    ("nll", "nll"),
                    ("brier", "brier"),
                    ("ece_equal_frequency", "ece_equal_frequency"),
                    ("aurc", "full_aurc_historical_column"),
                ]:
                    if col in sub.columns:
                        rec[f"{new}_mean"] = float(sub[col].mean())
                        rec[f"{new}_sd"] = float(sub[col].std(ddof=1))
                pbmc_rows.append(rec)
    # Recompute pAURC/E-AURC from PHASE11 saved full curves (seed0 era) as documentation
    curves_path = PHASE11 / "logs" / "risk_coverage_curves.json"
    curve_rows = []
    if curves_path.exists():
        curves = json.loads(curves_path.read_text())
        for c in curves:
            cov = np.asarray(c["coverage"], dtype=float)
            risk = np.asarray(c["risk"], dtype=float)
            full = float(np.mean(risk))
            paurc = float(np.mean(risk[cov >= 0.5 - 1e-12]))
            # Cannot compute E-AURC from curve alone without n_correct; leave nan unless recoverable
            curve_rows.append(
                {
                    "experiment": f"PBMC_phase11_{c['direction']}",
                    "model": c["model"],
                    "full_aurc_from_curve": full,
                    "paurc_0.5_1.0_from_curve": paurc,
                    "note": "PHASE11 curves; E-AURC requires cell-level correctness counts",
                }
            )

    sel_path = CAL / "selective_prediction_corrected.csv"
    sel.to_csv(sel_path, index=False)
    # Lawlor summary by model
    summary = (
        sel.groupby("model")[
            [
                "macro_f1",
                "accuracy",
                "balanced_accuracy",
                "error_prevalence",
                "error_auroc",
                "error_auprc",
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
                "risk_90",
                "risk_80",
                "risk_70",
                "risk_60",
                "risk_50",
            ]
        ]
        .agg(["mean", "std"])
        .reset_index()
    )
    # flatten columns
    summary.columns = ["_".join([c for c in col if c]).strip("_") for col in summary.columns.values]
    rel = pd.DataFrame(rel_rows)
    # Attach PBMC seed summaries
    pbmc_df = pd.DataFrame(pbmc_rows)
    curve_df = pd.DataFrame(curve_rows)
    rel_path = CAL / "full_reliability_metrics.csv"
    # Combine: lawlor run-level + pbmc summaries in one long file with schema marker
    rel["table_block"] = "lawlor_run_level"
    blocks = [rel]
    if not pbmc_df.empty:
        pbmc_df["table_block"] = "pbmc_20seed_summary_from_primary_csv"
        blocks.append(pbmc_df)
    if not curve_df.empty:
        curve_df["table_block"] = "pbmc_phase11_curve_aurc"
        blocks.append(curve_df)
    combined = pd.concat(blocks, ignore_index=True, sort=False)
    combined.to_csv(rel_path, index=False)
    summary.to_csv(CAL / "lawlor_reliability_summary_by_model.csv", index=False)

    # Naming audit for AURC
    (LOGS / "aurc_implementation_audit.md").write_text(
        """# AURC implementation audit (PART 8)

## Historical implementation

`src/experiments/phase11_uncertainty.py::risk_coverage`:

- rank cells by increasing uncertainty
- risk_k = 1 - accuracy of the first k cells
- **AURC = mean(risk_k) for k=1..n**

This is the **full risk-coverage AURC over all prefixes** (coverage grid is NOT used
for the `aurc` scalar).

`COVERAGE_GRID = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)` is used only for discrete
selective tables (accuracy/risk/macro-F1 at those coverages).

## Correction required in manuscript language

Do **not** rename historical `aurc` to partial AURC — that would be incorrect.

Instead:

- Keep calling the stored scalar **full AURC**.
- Additionally report **pAURC_[0.5,1.0]** = mean risk over prefixes with coverage ≥ 0.5.
- Report **E-AURC** = AURC − AURC* where AURC* is the AURC of a perfect uncertainty
  ranking (all correct before all incorrect), following the selective-classification
  excess-AURC idea (Geifman & El-Yaniv, 2017; standard excess over optimal ranking).

## Artifacts

- `calibration/selective_prediction_corrected.csv` (Lawlor run-level)
- `calibration/full_reliability_metrics.csv`
"""
    )
    return sel_path, rel_path


def write_feasibility_docs():
    (OUT / "baselines" / "modern_baseline_feasibility.md").write_text(
        """# Modern multimodal baseline feasibility (PART 10)

## Installed stack (project pin)

- `scvi-tools==1.3.3` (see `requirements.txt` / `environment.yml`)
- Python 3.11 target env on server; local analysis used lightweight venv

## MultiVI

| Question | Verdict |
|---|---|
| 1. Native RNA + protein? | **No / not fair.** MultiVI in scvi-tools is designed for **RNA + ATAC** (chromatin). Protein (ADT) is not a first-class MultiVI modality in 1.3.3. |
| 2. Same cells/features? | Would require forcing ADT into an ATAC-like slot — **not honest**. |
| 3. Missing protein semantics? | N/A if modality misuse. |
| 4. Latent dim match? | Possible but irrelevant if modality unfair. |
| 5. Same train/test design? | Possible mechanically, not scientifically fair. |
| 6. Evaluation without label leakage? | Same classifier discipline possible. |

**Decision: DO NOT force MultiVI into the RNA+protein benchmark.**

## Cobolt

Cobolt is primarily **RNA + ATAC** multiomic integration. Not a native CITE-seq
RNA+protein generative baseline for this comparison.

## Preferable ONE modern RNA+protein candidate (if any)

Candidates with native CITE-seq / RNA+protein support in the literature:

- **scArches-compatible totalVI** (already the multimodal deep baseline)
- **MOFA+** (already included; linear multi-view factor model)
- Possible alternatives: **scMaui**, **scAI**, or newer CITE-seq VAEs — **not currently installed**

## Verdict

**STOP before training an additional deep multimodal model** until a method with
**native RNA+protein likelihood support** is installed and can share:

- same cells/features
- same train/source vs target split
- honest missing-protein semantics if claimed
- comparable latent dimensionality
- leakage-free evaluation

**Additional modern model actually run in this revision round so far: NO.**

Recommended next step (P2, only if needed): install/evaluate one native
RNA+protein method OR elevate the **simple RNA+protein concatenated PCA**
control (PART 5B) as the architecture-vs-information contrast, which does not
require a second deep generative model.
"""
    )
    (OUT / "baselines" / "missingness_feasibility.md").write_text(
        """# Missing-modality feasibility (PART 11)

## Historical mistake (do not repeat)

Zero-filling cell-wise proteins and calling that “missing modality” is invalid
for totalVI 1.3.3 because the API does not provide arbitrary cell-level protein
missingness masks with correct likelihood semantics.

## totalVI 1.3.3

Prior project audit (`scripts/09_validate_missing_modality_api.py` / PHASE9 logs):
**no honest cell-wise missing-protein likelihood mask** for arbitrary patterns.

Therefore:

- **Cell-wise missing protein evaluation: UNSUPPORTED** for totalVI 1.3.3.
- Do not fake it.

## Valid options

A. Panel-level protein missingness **only if** API genuinely masks likelihood — not confirmed as general cell-wise.
B. **Antibody-panel dropout** (remove antibodies / random 25–50% of panel) — feasible, interpretable, not cell-wise missingness.
C. A modern model with explicit missing-modality support — not yet selected/fairly configured.

## Antibody-panel dropout status

**Designed but not yet executed** in this revision round (requires GPU training
on RTX 4060 with frozen clean models). Script scaffold to be added under
`corruption_shift/` / `controls/` before any run.

## Cell-wise missing protein

**Documented as unsupported** under the current totalVI pin.
"""
    )


def write_control_scaffolds():
    """PART 5–6 design notes / stop conditions — no training yet."""
    (OUT / "controls" / "CONTROL_STATUS.md").write_text(
        """# Controls status (PART 5–6)

## PART 5A — RNA-only annotation labels

**Status:** NOT YET RUN (requires reference RNA-only label transfer + freeze).

Hard requirements before training:

- Same ontology where possible
- **No target protein** during label assignment
- Freeze labels + confidence threshold **before** model comparison
- Evaluate RNA PCA, scVI_matched, concat PCA, MOFA+, totalVI

## PART 5B — Simple multimodal baseline (RNA PCA + protein PCA)

**Status:** NOT YET RUN (low-cost; should be first GPU/CPU control).

Hard requirements:

- Fit PCA / scaling on **train/source only** for inductive evaluations
- Transform target with train-derived operators
- Same classifier as deep models
- Dimensionality comparable where reasonable

## PART 6 — clean-trained → corrupted-test

**Status:** NOT YET RUN (needs frozen clean models + test-only protein corruption).

Distinct from historical PHASE8 training-time degradation (keep; rename).

Corruption fractions: 0, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85 on **test protein only**.
scVI should be invariant (negative control).

## Execution venue

Heavy runs: RTX 4060 server via `bash scripts/run_on_gpu.sh` inside tmux,
with done-markers under this revision tree. Local Mac workspace currently lacks
PHASE11B probability/model artifacts needed to skip retraining for some PBMC
controls.
"""
    )


def main():
    STATS.mkdir(parents=True, exist_ok=True)
    CAL.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    (OUT / "baselines").mkdir(parents=True, exist_ok=True)
    (OUT / "controls").mkdir(parents=True, exist_ok=True)
    ledger = write_statistical_unit_ledger()
    sel, rel = recompute_lawlor_selective_and_reliability()
    write_feasibility_docs()
    write_control_scaffolds()
    print("Wrote", ledger)
    print("Wrote", sel)
    print("Wrote", rel)


if __name__ == "__main__":
    main()
