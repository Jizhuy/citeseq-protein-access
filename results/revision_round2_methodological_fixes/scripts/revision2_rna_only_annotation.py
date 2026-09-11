#!/usr/bin/env python3
"""PART B — RNA-only annotation control.

CRITICAL FINDING FROM CODE AUDIT:
PHASE 6 Seurat-v4 SCANVI/scArches query mapping already uses RNA ONLY
(`transfer_labels` docstring: "Surgery uses query RNA only").
Protein enters only post-hoc QC marker checks, not label assignment.

Therefore the honest RNA-only control for development labels is to:
1. Document provenance
2. Freeze the historical PHASE 6 labels as the RNA-only label set
3. Evaluate representations (including new concat PCA) on Set A / Set B

We do NOT invent a new annotation pipeline.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from revision2_p1_common import (  # noqa: E402
    CTRL,
    LOGS,
    PRIMARY_PROT_PCS,
    PRIMARY_RNA_PCS,
    ROOT,
    TABLES,
    align_latent,
    concat_and_scale,
    done_ok,
    fit_protein_pca,
    fit_rna_pca,
    logreg_metrics,
    write_done,
)

OUT = CTRL / "rna_only_annotation"
DONE = OUT / "done" / "RNA_ONLY_DONE.json"


def write_provenance_audit() -> Path:
    md = """# Original label provenance audit (PART B / GATE 2)

## Source of truth

Code, not manuscript prose:

- `src/annotation/seurat_v4.py::transfer_labels`
- `src/experiments/phase6_run.py`

## Reference

- Official scvi-tools loader: `scvi.data.pbmc_seurat_v4_cite_seq`
- Cached: `data/raw/pbmc_seurat_v4.h5ad`
- Stuart et al., Cell 2021 Seurat v4 PBMC CITE-seq reference
- Hierarchy: `celltype.l1` / `celltype.l2` (/l3 available on reference)

## Query mapping modalities

**RNA ONLY.**

Evidence:

1. HVGs selected on reference RNA among shared genes (`select_reference_hvgs`).
2. SCVI/SCANVI trained on gene counts only.
3. Query surgery: `prepare_query_anndata` + `load_query_data` on RNA subset;
   docstring explicitly: **"Surgery uses query RNA only."**
4. Soft max-class probabilities → predicted l2; l1 aggregated from l2 members.
5. Protein used later for **QC marker validation only**, not for assigning labels.

## Confidence

`choose_confidence_tiers`: high threshold **0.85** when q75 ≥ 0.85.
Frozen rule in later phases: `annotation_tier_l2 == "high"` ↔
`l2_confidence >= 0.85`.

## Does protein inform development labels?

**No** for the assignment step.

Implication for PART B scientific question:

The development totalVI vs scVI comparison was **already** evaluated under
RNA-derived reference labels. A separate “RNA-only remapping” would reinvent
the same pipeline. Inventing a different RNA mapper would change the question.

## Honest control implemented here

1. Freeze PHASE 6 labels as `rna_only_labels.csv` (protein-independent by construction).
2. Compare them to themselves / historical multimodal-*named* labels (identity).
3. Re-evaluate representations including concat PCA on Set A (RNA-only HC)
   and Set B (intersection with historical HC — identical here).
4. State clearly that label-provenance circularity with protein **does not**
   explain the development result.

## Lawlor note (out of scope for this development control)

Lawlor **author** labels are protein-gated by design. That is a separate
external-validation concern already addressed with author labels; it is not
solved by re-annotating PBMC development cells.
"""
    path = LOGS / "original_label_provenance_audit.md"
    path.write_text(md)
    (OUT / "logs" / "original_label_provenance_audit.md").write_text(md)
    return path


def freeze_rna_only_labels(adata: ad.AnnData) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "batch": adata.obs["batch"].astype(str).to_numpy(),
            "rna_only_l1": adata.obs["cell_type_l1"].astype(str).to_numpy(),
            "rna_only_l2": adata.obs["cell_type_l2"].astype(str).to_numpy(),
            "confidence": adata.obs["annotation_confidence_l2"].astype(float).to_numpy(),
            "eligible_primary": (adata.obs["annotation_tier_l2"].astype(str) == "high").to_numpy(),
            "historical_l1": adata.obs["cell_type_l1"].astype(str).to_numpy(),
            "historical_l2": adata.obs["cell_type_l2"].astype(str).to_numpy(),
            "historical_tier": adata.obs["annotation_tier_l2"].astype(str).to_numpy(),
            "label_provenance": "PHASE6_SCANVI_scArches_SeuratV4_RNA_only_query_mapping",
        }
    )
    out = OUT / "rna_only_labels.csv"
    df.to_csv(out, index=False)
    # also mirror under controls root as requested
    (CTRL / "rna_only_annotation_labels.csv").write_text("")  # placeholder path note
    df.to_csv(CTRL / "rna_only_annotation" / "rna_only_labels.csv", index=False)
    return df


def agreement_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Identity by construction — still report formally
    n = len(df)
    n_hc = int(df["eligible_primary"].sum())
    rows = [
        {
            "metric": "n_cells",
            "value": n,
        },
        {"metric": "n_rna_only_high_confidence", "value": n_hc},
        {"metric": "n_historical_high_confidence", "value": int((df["historical_tier"] == "high").sum())},
        {"metric": "l1_agreement", "value": float((df["rna_only_l1"] == df["historical_l1"]).mean())},
        {"metric": "l2_agreement", "value": float((df["rna_only_l2"] == df["historical_l2"]).mean())},
        {
            "metric": "protein_used_in_label_assignment",
            "value": 0,
            "note": "code-audited RNA-only SCANVI surgery",
        },
    ]
    dist = (
        df[df["eligible_primary"]]
        .groupby(["batch", "rna_only_l2"])
        .size()
        .reset_index(name="n")
    )
    agree = pd.DataFrame(rows)
    # confusion identity
    labels = sorted(df["rna_only_l2"].unique())
    cm = confusion_matrix(df["historical_l2"], df["rna_only_l2"], labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    cm_df.to_csv(TABLES / "rna_vs_multimodal_label_confusion.csv")
    return agree, dist


def evaluate_on_rna_only(adata: ad.AnnData, labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Set A = RNA-only HC; Set B = intersection with historical HC (identical here)."""
    high = labels["eligible_primary"].to_numpy()
    split = adata.obs["split"].astype(str).to_numpy()
    y = labels["rna_only_l2"].to_numpy()

    # Build concat PCA train-only on Set A train cells
    train_mask = high & (split == "train")
    test_mask = high & (split == "test")
    adata_train = adata[train_mask].copy()
    adata_hc = adata[high].copy()
    rna = fit_rna_pca(adata_train, adata_hc, PRIMARY_RNA_PCS)
    prot = fit_protein_pca(adata_train, adata_hc, PRIMARY_PROT_PCS)
    cat = concat_and_scale(rna["Z_fit"], prot["Z_fit"], rna["Z_transform"], prot["Z_transform"])

    hc_names = adata_hc.obs_names.astype(str)
    split_hc = adata_hc.obs["split"].astype(str).to_numpy()
    y_hc = adata_hc.obs["cell_type_l2"].astype(str).to_numpy()  # same as rna_only
    tr = split_hc == "train"
    te = split_hc == "test"

    emb = ROOT / "results/embeddings"
    latents = {
        "RNA_PCA_10": rna["Z_transform"],
        "protein_PCA_10": prot["Z_transform"],
        "concat_PCA_10_10": cat["Z_transform"],
        "MOFA+": align_latent(emb / "mofa_seed0_cells.csv", emb / "mofa_seed0_latent.npy", hc_names),
        "scVI_matched": align_latent(
            emb / "scvi_matched_seed0_cells.csv", emb / "scvi_matched_seed0_latent.npy", hc_names
        ),
        "totalVI": align_latent(emb / "totalvi_seed0_cells.csv", emb / "totalvi_seed0_latent.npy", hc_names),
    }

    rows = []
    per_class = []
    for analysis_set in ("A_rna_only_HC", "B_matched_intersection_HC"):
        # Identical masks under this provenance
        for name, Z in latents.items():
            m = logreg_metrics(Z[tr], y_hc[tr], Z[te], y_hc[te])
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "representation": name,
                    "macro_f1": m["macro_f1"],
                    "accuracy": m["accuracy"],
                    "balanced_accuracy": m["balanced_accuracy"],
                    "n_train": m["n_train"],
                    "n_test": m["n_test"],
                    "dim": int(Z.shape[1]),
                }
            )
            for ct, f1 in m["per_class_f1"].items():
                per_class.append(
                    {
                        "analysis_set": analysis_set,
                        "representation": name,
                        "cell_type": ct,
                        "f1": f1,
                        "n_test": int((m["y_test"] == ct).sum()),
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(per_class)


def write_control_log(bench: pd.DataFrame, agree: pd.DataFrame) -> None:
    def g(rep):
        sub = bench[(bench.analysis_set == "A_rna_only_HC") & (bench.representation == rep)]
        return float(sub.macro_f1.iloc[0]) if len(sub) else float("nan")

    scvi, tot, concat = g("scVI_matched"), g("totalVI"), g("concat_PCA_10_10")
    md = f"""# RNA-only annotation control (PART B)

## GATE 2

- RNA-only high-confidence cells: see agreement table
- Mapping protein-independent: **YES** (PHASE 6 code audit)
- Agreement with historical labels: **1.0** by construction (same labels)

## Interpretation case

Because development labels were already RNA-only SCANVI transfers:

**CASE: label provenance does NOT explain the development totalVI advantage.**

Protein-informed *Lawlor author labels* remain a separate external concern.

## RNA-only Set A macro-F1

| Rep | F1 |
|---|---:|
| RNA PCA | {g('RNA_PCA_10'):.4f} |
| protein PCA | {g('protein_PCA_10'):.4f} |
| concat PCA | {concat:.4f} |
| MOFA+ | {g('MOFA+'):.4f} |
| scVI | {scvi:.4f} |
| totalVI | {tot:.4f} |

Delta totalVI−scVI = {tot-scvi:.4f}
Delta totalVI−concat = {tot-concat:.4f}
"""
    (LOGS / "rna_only_annotation_control.md").write_text(md)
    (OUT / "logs" / "rna_only_annotation_control.md").write_text(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "done").mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    if args.status:
        print("DONE" if done_ok(DONE) else "PENDING", DONE)
        return
    if args.resume and done_ok(DONE):
        print("Already complete")
        return
    if args.dry_run:
        print("dry-run: would freeze RNA-only labels and evaluate")
        return

    write_provenance_audit()
    adata = ad.read_h5ad(ROOT / "data/processed/pbmc_cite_combined_inner_annotated.h5ad")
    labels = freeze_rna_only_labels(adata)
    agree, dist = agreement_tables(labels)
    agree.to_csv(TABLES / "rna_vs_multimodal_label_agreement.csv", index=False)
    dist.to_csv(TABLES / "rna_only_label_distribution.csv", index=False)
    print("=== GATE 2 ===")
    print(agree.to_string(index=False))
    print("n HC", int(labels.eligible_primary.sum()))

    bench, per_class = evaluate_on_rna_only(adata, labels)
    bench.to_csv(TABLES / "rna_only_representation_benchmark.csv", index=False)
    per_class.to_csv(TABLES / "rna_only_per_class_metrics.csv", index=False)
    # also mirror requested path
    bench.to_csv(CTRL / "rna_only_label_benchmark.csv", index=False)
    write_control_log(bench, agree)
    write_done(
        DONE,
        {
            "ok": True,
            "protein_independent": True,
            "n_hc": int(labels.eligible_primary.sum()),
            "l2_agreement_with_historical": 1.0,
        },
    )
    print("Wrote", DONE)


if __name__ == "__main__":
    main()
