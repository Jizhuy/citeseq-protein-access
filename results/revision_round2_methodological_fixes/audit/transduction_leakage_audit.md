# Transduction / leakage audit (code-grounded)

**Rule applied:** unlabeled target features in scArches query adaptation = **unsupervised transductive adaptation**, not automatic label leakage.

## Summary verdict

| Experiment | Target used to train representation? | Unlabeled target in scArches? | Target labels before final eval? | Genuine label leakage? |
|---|---|---|---|---|
| Representation benchmark (PHASE6 VAEs) | Yes — all query cells jointly | N/A (no query surgery for VAEs) | Labels used for classifier/eval after transfer annotation | **No** for VAE; annotation is separate transfer |
| Protein corruption (PHASE8) | Yes — all cells with corrupted protein | No | PHASE6 labels for clf/eval | **No** |
| Cross-dataset A/B (PHASE10/11B) | **No** (source batch only) | **Yes** | Class eligibility uses target label counts; clf fit source-only | **No label leakage into model**; evaluation-design uses target labels |
| Source-size (PHASE10B) | **No** (subsampled source) | **Yes** | Same as PHASE10 | **No** (same caveat) |
| Variance (PHASE12A) | **No** (source subset) | **Yes** | Same | **No** |
| Lawlor (PHASE12B) | **No** (source donors) | **Yes** | Eligibility rule uses target supports every fold; clf source-only | **No** (same caveat) |
| Selective prediction | inherits | inherits | labels only at risk scoring | **No** |

## Seven-question checklist (Lawlor + transfer)

1. **Was the target dataset used to train the representation?**  
   Transfer/Lawlor: **No** (source only). PHASE6/8 joint VAEs: **Yes** (all cells).

2. **Was unlabeled target data used in scArches query adaptation?**  
   Transfer/Lawlor/11B/12A: **Yes**. This is unsupervised transduction.

3. **Were target labels used anywhere before final evaluation?**  
   - Training/adaptation/classifier fit: **No**.  
   - **Evaluation design:** eligible class sets use target label counts (`class_plan` / Lawlor `>=20 source & >=10 target in every fold`).  
   This is **target-label-informed evaluation design**, not parameter leakage.

4. **Were feature-selection steps fitted on train/source only or on all cells?**  
   VAE transfer/Lawlor: **no HVG** (fixed shared gene universe).  
   PHASE6 annotation HVGs: reference-only.  
   PHASE7 MOFA+/PCA_RNA: fitted on **all combined cells** (joint unsupervised).

5. **Were PCA / scaling / HVG / normalization parameters fitted on evaluation cells?**  
   VAE pipelines: **No** such steps.  
   PHASE7 PCA/MOFA: **Yes** (all cells) — note for baseline fairness discussions.

6. **Was classifier preprocessing fitted only on classifier-training cells?**  
   **Yes vacuously** — logistic regression on raw latents; no StandardScaler.

7. **Were class eligibility thresholds defined using target labels?**  
   **Yes** for transfer and Lawlor.

## Correct manuscript terminology

Use: **donor-held-out / dataset-held-out label evaluation after unsupervised query adaptation**.

Do **not** imply fully inductive prediction on completely unseen target feature distributions without adaptation.

Do **not** call valid unsupervised transduction “leakage.”

## Hard-stop status

- Genuine target-label leakage into representation training/adaptation/classifier fit: **NOT found**.  
- Source/target embeddings mixed (historical PHASE10 bug): corrected in prior phases; Lawlor/11B use post-adaptation coordinates consistently.  
- Continue with statistical corrections.

## Key code citations

- `src/experiments/phase10_cross_dataset.py` — `query_adapt`, `class_plan`, `target_labels_used_in_fitting: False`
- `src/experiments/phase12b_external.py` — donor folds, `evaluate_run`, unlabeled query adapt
- `src/experiments/phase12a_stats.py` — per-model MoM variance (not pooled without model term)
- `src/experiments/phase11_uncertainty.py` — `risk_coverage` full-prefix AURC; `COVERAGE_GRID` for tables
