# Matched Transfer Code Audit

**Script:** `results/revision_round3_information_vs_model/transfer_control/scripts/revision3_matched_concat_transfer.py`  
**Table:** `results/revision_round3_information_vs_model/transfer_control/tables/matched_transfer_control.csv`  
**Audit date:** 2026-09-23  
**Scope:** Code inspection only (no retrain / no new experiment).

## Frozen numbers (verified from CSV)

| Direction | Inductive concat | Transductive concat | scVI (scArches) | totalVI (scArches) | totalVI − concat_trans |
|---|---:|---:|---:|---:|---:|
| A | 0.735 | 0.636 | 0.656 | 0.652 | ≈ +0.016 |
| B | 0.734 | 0.767 | 0.601 | 0.654 | ≈ −0.113 |

Exact CSV values: inductive 0.734952 / 0.733927; transductive 0.636313 / 0.767041; totalVI 0.652383 / 0.653966.

## Direction / modality integrity

- Direction A: source=`PBMC10k`, target=`PBMC5k`.
- Direction B: source=`PBMC5k`, target=`PBMC10k`.
- Primary concat dimensionality: 10 RNA PCs + 10 protein PCs = 20.
- Shared gene/protein ordering inherited from combined AnnData inner join.

## Target-label access audit

### Unsupervised feature fitting (transductive concat)

| Step | Uses target features? | Uses target labels? |
|---|---|---|
| Joint AnnData concat | Yes (unlabeled features) | No |
| RNA HVG / normalize / scale / PCA (`fit_rna_pca(joint, joint)`) | Yes | No |
| Protein CLR / z-score / PCA (`fit_protein_pca(joint, joint)`) | Yes | No |
| Final concat `StandardScaler.fit(X_joint)` | Yes | No |
| LogisticRegression training | Source latents only | **Source labels only** |
| Final scoring | Target latents | **Target labels only at eval** |

Inductive concat: HVG/PCA/scaling fit on **source only**; target transformed; classifier on source labels; target labels only at eval.

### Evaluation-set definition (not feature supervision)

`transfer_class_plan` uses source **and** target `cell_type_l2` / `annotation_tier_l2` solely to define eligible shared high-confidence classes and evaluation indices (min counts). This is standard transfer **scoring-set construction**. It does **not** enter HVG, PCA, scaling, component selection, classifier fitting, or hyperparameter selection.

**Verdict:** No target-label leakage into unsupervised feature transforms or classifier training. Target labels used only for evaluation-set definition and final scoring.

## Algorithmic matching caveat

- Transductive concat = joint unsupervised PCA/scaling on source + unlabeled target.
- scVI/totalVI = frozen PHASE11B scArches unsupervised query adaptation (20 seeds; not recomputed).
- Therefore: **target-data-access matched**, **not** a fully matched architecture comparison.

## Integrity artifacts

- `logs/matched_transfer_integrity.json` records `target_labels_used_before_eval: False`.
- `logs/matched_transfer_control_design.md` documents the design.
- Done marker: `done/MATCHED_TRANSFER_CONTROL_DONE.json`.

## Conclusion

**PASS.** Numbers match frozen Case C. No accidental label access into representation fitting or classifier training.
