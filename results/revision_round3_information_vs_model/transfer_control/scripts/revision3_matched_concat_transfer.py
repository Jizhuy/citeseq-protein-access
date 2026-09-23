#!/usr/bin/env python3
"""Revision Round 3 — target-data-access-matched transductive concat PCA control.

Preserves the existing source-only inductive concat baseline and adds a
joint unsupervised (source + unlabeled target) concat PCA baseline that
matches ACCESS TO UNLABELED TARGET FEATURES only.

Target labels are used ONLY at final scoring.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[4]
P1_SCRIPTS = ROOT / "results/revision_round2_methodological_fixes/scripts"
sys.path.insert(0, str(P1_SCRIPTS))

from revision2_p1_common import (  # noqa: E402
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    concat_and_scale,
    fit_protein_pca,
    fit_rna_pca,
    logreg_metrics,
)
from revision2_simple_multimodal_baseline import (  # noqa: E402
    _fit_prot_once,
    _fit_rna_once,
    transfer_class_plan,
)

OUT = ROOT / "results/revision_round3_information_vs_model/transfer_control"
TABLES = OUT / "tables"
LOGS = OUT / "logs"
FIGS = OUT / "figures_preview"
DONE = OUT / "done" / "MATCHED_TRANSFER_CONTROL_DONE.json"
NAVY = "#1B3A5F"
BLUE = "#2F6FAE"
LBLUE = "#9DC3E6"
GRAY = "#5A5A5A"


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_adata() -> ad.AnnData:
    path = ROOT / "data/processed/pbmc_cite_combined_inner_annotated.h5ad"
    return ad.read_h5ad(path)


def _concat_batches(source: ad.AnnData, target: ad.AnnData) -> ad.AnnData:
    """Concatenate source then target; preserve order for index split."""
    # Use AnnData concat without labels influencing anything unsupervised
    joint = ad.concat([source, target], join="inner", merge="same", index_unique=None)
    # Ensure unique names if needed
    if joint.obs_names.duplicated().any():
        joint.obs_names_make_unique()
    return joint


def inductive_concat(source, target, src_idx, tgt_idx, classes):
    rna = _fit_rna_once(source, [source, target], PRIMARY_RNA_PCS)
    prot = _fit_prot_once(source, [source, target], PRIMARY_PROT_PCS)
    X_s = np.hstack([rna["Zs"][0], prot["Zs"][0]])
    X_t = np.hstack([rna["Zs"][1], prot["Zs"][1]])
    scaler = StandardScaler().fit(X_s)  # SOURCE only
    Z_s, Z_t = scaler.transform(X_s), scaler.transform(X_t)
    y_s = source.obs["cell_type_l2"].astype(str).to_numpy()[src_idx]
    y_t = target.obs["cell_type_l2"].astype(str).to_numpy()[tgt_idx]
    keep_s = np.isin(y_s, classes)
    keep_t = np.isin(y_t, classes)
    m = logreg_metrics(Z_s[src_idx][keep_s], y_s[keep_s], Z_t[tgt_idx][keep_t], y_t[keep_t])
    return m, {
        "fit_cells": "source_only",
        "scaler_fit": "source_only",
        "target_labels_used_before_eval": False,
        "n_joint_cells": int(source.n_obs),
    }


def transductive_concat(source, target, src_idx, tgt_idx, classes):
    """Fit unsupervised transforms on source + unlabeled target features."""
    n_s = source.n_obs
    joint = _concat_batches(source, target)
    assert joint.n_obs == source.n_obs + target.n_obs

    # Fit RNA/protein PCA on joint unlabeled features; transform joint then split
    rna = fit_rna_pca(joint, joint, PRIMARY_RNA_PCS)
    prot = fit_protein_pca(joint, joint, PRIMARY_PROT_PCS)
    Z_rna_s = rna["Z_transform"][:n_s]
    Z_rna_t = rna["Z_transform"][n_s:]
    Z_prot_s = prot["Z_transform"][:n_s]
    Z_prot_t = prot["Z_transform"][n_s:]

    X_s = np.hstack([Z_rna_s, Z_prot_s])
    X_t = np.hstack([Z_rna_t, Z_prot_t])
    # Final concat scaling also unsupervised on source+target features
    X_joint = np.vstack([X_s, X_t])
    scaler = StandardScaler().fit(X_joint)
    Z_s, Z_t = scaler.transform(X_s), scaler.transform(X_t)

    y_s = source.obs["cell_type_l2"].astype(str).to_numpy()[src_idx]
    y_t = target.obs["cell_type_l2"].astype(str).to_numpy()[tgt_idx]
    keep_s = np.isin(y_s, classes)
    keep_t = np.isin(y_t, classes)

    # Integrity: target labels never enter fit — only used here for scoring
    m = logreg_metrics(Z_s[src_idx][keep_s], y_s[keep_s], Z_t[tgt_idx][keep_t], y_t[keep_t])
    return m, {
        "fit_cells": "source_plus_unlabeled_target",
        "scaler_fit": "source_plus_unlabeled_target",
        "target_labels_used_before_eval": False,
        "n_joint_cells": int(joint.n_obs),
        "n_source": int(source.n_obs),
        "n_target": int(target.n_obs),
        "rna_hvgs": len(rna["hvgs"]),
        "dim": int(Z_s.shape[1]),
    }


def inductive_modality(source, target, src_idx, tgt_idx, classes, modality: str):
    if modality == "rna":
        pack = _fit_rna_once(source, [source, target], PRIMARY_RNA_PCS)
        Zs, Zt = pack["Zs"][0], pack["Zs"][1]
        dim = PRIMARY_RNA_PCS
    else:
        pack = _fit_prot_once(source, [source, target], PRIMARY_PROT_PCS)
        Zs, Zt = pack["Zs"][0], pack["Zs"][1]
        dim = PRIMARY_PROT_PCS
    y_s = source.obs["cell_type_l2"].astype(str).to_numpy()[src_idx]
    y_t = target.obs["cell_type_l2"].astype(str).to_numpy()[tgt_idx]
    keep_s = np.isin(y_s, classes)
    keep_t = np.isin(y_t, classes)
    m = logreg_metrics(Zs[src_idx][keep_s], y_s[keep_s], Zt[tgt_idx][keep_t], y_t[keep_t])
    m["dim"] = dim
    return m


def frozen_vae_rows() -> list[dict]:
    p11 = ROOT / "results/phase11b_uncertainty_robustness/tables/phase11b_20seed_primary.csv"
    rows = []
    if not p11.exists():
        return rows
    df = pd.read_csv(p11)
    for (dname, model), sub in df.groupby(["direction", "model"]):
        label = "scVI_matched" if model == "scvi_matched" else "totalVI"
        rows.append(
            {
                "direction": dname,
                "representation": label,
                "adaptation": "scArches_unsupervised_transductive",
                "target_data_access": "unlabeled_target_features_via_scArches",
                "macro_f1": float(sub.macro_f1.mean()),
                "macro_f1_sd": float(sub.macro_f1.std(ddof=1)),
                "accuracy": float("nan"),
                "balanced_accuracy": float("nan"),
                "n_train": float("nan"),
                "n_test": float("nan"),
                "n_classes": float("nan"),
                "dim": 20,
                "n_seeds": int(len(sub)),
                "note": "Frozen PHASE11B 20-seed means; not recomputed",
            }
        )
    return rows


def run_control(adata: ad.AnnData) -> pd.DataFrame:
    directions = {
        "direction_A": ("PBMC10k", "PBMC5k"),
        "direction_B": ("PBMC5k", "PBMC10k"),
    }
    rows = []
    audits = []
    for dname, (src_batch, tgt_batch) in directions.items():
        source = adata[adata.obs["batch"].astype(str) == src_batch].copy()
        target = adata[adata.obs["batch"].astype(str) == tgt_batch].copy()
        plan = transfer_class_plan(source, target, "cell_type_l2")
        classes = plan["eligible_shared_classes"]
        src_idx = np.asarray(plan["source_eval_index"])
        tgt_idx = np.asarray(plan["target_eval_index"])

        # RNA inductive
        m_rna = inductive_modality(source, target, src_idx, tgt_idx, classes, "rna")
        rows.append(
            {
                "direction": dname,
                "representation": "RNA_PCA_10",
                "adaptation": "none_inductive_source_fit_only",
                "target_data_access": "source_only",
                "macro_f1": m_rna["macro_f1"],
                "macro_f1_sd": float("nan"),
                "accuracy": m_rna["accuracy"],
                "balanced_accuracy": m_rna["balanced_accuracy"],
                "n_train": m_rna["n_train"],
                "n_test": m_rna["n_test"],
                "n_classes": m_rna["n_classes_test"],
                "dim": m_rna.get("dim", 10),
                "n_seeds": 1,
                "note": "source-fitted inductive RNA PCA",
            }
        )

        # Concat inductive (preserve)
        m_ind, meta_ind = inductive_concat(source, target, src_idx, tgt_idx, classes)
        rows.append(
            {
                "direction": dname,
                "representation": "concat_PCA_10_10",
                "adaptation": "none_inductive_source_fit_only",
                "target_data_access": "source_only",
                "macro_f1": m_ind["macro_f1"],
                "macro_f1_sd": float("nan"),
                "accuracy": m_ind["accuracy"],
                "balanced_accuracy": m_ind["balanced_accuracy"],
                "n_train": m_ind["n_train"],
                "n_test": m_ind["n_test"],
                "n_classes": m_ind["n_classes_test"],
                "dim": 20,
                "n_seeds": 1,
                "note": "source-fitted inductive concat PCA baseline (preserved)",
            }
        )

        # Concat transductive (NEW)
        m_tr, meta_tr = transductive_concat(source, target, src_idx, tgt_idx, classes)
        rows.append(
            {
                "direction": dname,
                "representation": "concat_PCA_10_10",
                "adaptation": "transductive_joint_unsupervised_feature_fit",
                "target_data_access": "source_plus_unlabeled_target_features",
                "macro_f1": m_tr["macro_f1"],
                "macro_f1_sd": float("nan"),
                "accuracy": m_tr["accuracy"],
                "balanced_accuracy": m_tr["balanced_accuracy"],
                "n_train": m_tr["n_train"],
                "n_test": m_tr["n_test"],
                "n_classes": m_tr["n_classes_test"],
                "dim": meta_tr["dim"],
                "n_seeds": 1,
                "note": "target-data-access-matched transductive concat; NOT algorithmically matched to scArches",
            }
        )
        audits.append({"direction": dname, "inductive": meta_ind, "transductive": meta_tr, "n_classes": len(classes)})

    rows.extend(frozen_vae_rows())
    df = pd.DataFrame(rows)
    (LOGS / "matched_transfer_integrity.json").write_text(json.dumps(audits, indent=2))
    return df


def write_design_log(df: pd.DataFrame) -> None:
    lines = [
        "# Matched transfer control design",
        "",
        f"- timestamp: {_utc()}",
        "- purpose: match ACCESS TO UNLABELED TARGET FEATURES for simple concat PCA",
        "- NOT a fully matched architecture comparison with scArches",
        "",
        "## Designs",
        "",
        "1. **Source-fitted inductive concat PCA** (preserved): HVG/PCA/scaling fit on SOURCE only; TARGET transformed; classifier on SOURCE labels.",
        "2. **Target-data-access-matched transductive concat PCA** (new): HVG/PCA/scaling fit on SOURCE + unlabeled TARGET features jointly; classifier still SOURCE labels only; TARGET labels used only at final scoring.",
        "3. **scVI / totalVI**: frozen PHASE11B scArches unsupervised transductive adaptation (20 seeds).",
        "",
        "## Integrity",
        "",
        "- Target labels never enter HVG selection, scaling, PCA, concat StandardScaler, or LogisticRegression training.",
        "- Primary endpoint: macro-F1; secondary: accuracy, balanced accuracy.",
        "",
        "## Results snapshot",
        "",
    ]
    for d in ["direction_A", "direction_B"]:
        lines.append(f"### {d}")
        sub = df[df.direction == d]
        for _, r in sub.iterrows():
            lines.append(
                f"- {r.representation} | {r.adaptation}: macro-F1={r.macro_f1:.4f}"
                + (f" ±{r.macro_f1_sd:.4f}" if pd.notna(r.macro_f1_sd) else "")
            )
        lines.append("")
    # Interpretation hint
    for d in ["direction_A", "direction_B"]:
        tr = df[(df.direction == d) & (df.representation == "concat_PCA_10_10") & (df.adaptation.str.contains("transductive_joint"))]
        tot = df[(df.direction == d) & (df.representation == "totalVI")]
        if len(tr) and len(tot):
            delta = float(tot.iloc[0].macro_f1) - float(tr.iloc[0].macro_f1)
            lines.append(f"- {d}: totalVI − transductive_concat = {delta:+.4f}")
    LOGS.joinpath("matched_transfer_control_design.md").write_text("\n".join(lines) + "\n")


def plot_preview(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    order = [
        ("RNA_PCA_10", "none_inductive_source_fit_only", "RNA ind."),
        ("concat_PCA_10_10", "none_inductive_source_fit_only", "concat ind."),
        ("concat_PCA_10_10", "transductive_joint_unsupervised_feature_fit", "concat trans."),
        ("scVI_matched", "scArches_unsupervised_transductive", "scVI"),
        ("totalVI", "scArches_unsupervised_transductive", "totalVI"),
    ]
    colors = [GRAY, LBLUE, BLUE, LBLUE, NAVY]
    for ax, d in zip(axes, ["direction_A", "direction_B"]):
        vals, labs, cols, yerr = [], [], [], []
        for (rep, adapt, lab), col in zip(order, colors):
            sub = df[(df.direction == d) & (df.representation == rep) & (df.adaptation == adapt)]
            if not len(sub):
                continue
            r = sub.iloc[0]
            vals.append(float(r.macro_f1))
            labs.append(lab)
            cols.append(col)
            yerr.append(float(r.macro_f1_sd) if pd.notna(r.macro_f1_sd) else 0.0)
        x = np.arange(len(vals))
        ax.bar(x, vals, yerr=yerr, color=cols, edgecolor=NAVY, capsize=3, width=0.72)
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=20, ha="right")
        ax.set_ylabel("macro-F1")
        ax.set_title(d.replace("_", " "), fontweight="bold", color=NAVY)
        ax.set_ylim(0.45, 0.85)
        for i, v in enumerate(vals):
            ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=7, color=NAVY)
    fig.suptitle("Target-data-access-matched transfer control", color=NAVY, fontweight="bold")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "matched_transfer_control.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGS / "matched_transfer_control.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    if args.status:
        print("DONE" if DONE.exists() else "NOT_DONE", DONE)
        return

    if not args.run:
        ap.print_help()
        return

    adata = load_adata()
    df = run_control(adata)
    out_csv = TABLES / "matched_transfer_control.csv"
    df.to_csv(out_csv, index=False)
    write_design_log(df)
    plot_preview(df)

    # Case interpretation
    cases = {}
    for d in ["direction_A", "direction_B"]:
        tr = float(
            df[
                (df.direction == d)
                & (df.representation == "concat_PCA_10_10")
                & (df.adaptation == "transductive_joint_unsupervised_feature_fit")
            ].iloc[0].macro_f1
        )
        tot = float(df[(df.direction == d) & (df.representation == "totalVI")].iloc[0].macro_f1)
        scvi = float(df[(df.direction == d) & (df.representation == "scVI_matched")].iloc[0].macro_f1)
        cases[d] = {
            "transductive_concat": tr,
            "totalVI": tot,
            "scVI": scvi,
            "totalVI_minus_transductive_concat": tot - tr,
            "totalVI_minus_scVI": tot - scvi,
        }
    payload = {
        "ok": True,
        "timestamp": _utc(),
        "table": str(out_csv),
        "cases": cases,
        "integrity": "target_labels_withheld_until_eval",
    }
    DONE.parent.mkdir(parents=True, exist_ok=True)
    DONE.write_text(json.dumps(payload, indent=2))
    print(json.dumps(cases, indent=2))
    print("wrote", out_csv)
    print("DONE", DONE)


if __name__ == "__main__":
    main()
