# PHASE 10 erratum — source latent coordinate system

Discovered during the PHASE 11 §11 reproduction check. No model was retrained
to establish or correct this. No historical PHASE 10 artifact was modified;
corrected values are written to `results/phase10_cross_dataset/tables_corrected/`.

## 1. Defect

`run_one_model` in `src/experiments/phase10_cross_dataset.py` wrote the
**post-adaptation** source embedding to disk but returned the
**pre-adaptation** embedding to the evaluator:

```python
source_latent = q["source_latent_adapted"]      # post-adaptation
_save_latent(latent_src_path, source_latent, source)   # written to disk
np.save(emb_dir / f"{tag}_source_latent_source_model.npy", src["latent"])
...
return {
    "source_latent": src["latent"],             # <-- PRE-adaptation returned
    "target_latent": q["latent"],
    ...
}
```

`evaluate_direction` consumes `pack["source_latent"]`, so a freshly-trained run
scored **pre-adaptation source coordinates against post-adaptation target
coordinates**. This is the exact mixing the PHASE 10 design forbids and that the
module docstring claims to have fixed.

The resume path was correct. `_load_completed_run` reads
`*_source_latent.npy` — the post-adaptation array — so any resumed run scored
correctly.

## 2. Scope

| model | drift (mean L2 source shift) | affected |
|---|---|---|
| scVI_matched | exactly 0.000000, arrays byte-identical | **no** |
| totalVI | 0.0150 – 0.0200 | **yes** |

scVI is immune because scvi-tools scArches masks gradients to the query batch's
one-hot columns, so the reference embedding is exactly frozen. totalVI's adapted
model does move the source embedding slightly, so the two arrays differ.

By seed:

| rows | reproduce the published table using |
|---|---|
| totalVI seeds 1–4 (16 rows, l1+l2) | pre-adaptation source |
| totalVI seed 0 (4 rows, l1+l2) | post-adaptation source |
| all scVI rows (20) | either (identical) |

Seed 0 was run twice — once fresh during gating and again in a later invocation
that resumed from saved artifacts — and `write_tables` merges with
`keep="last"`, so its rows were rewritten correctly. Seeds 1–4 ran fresh only.

All 40 grid cells were verified to match the published tables **exactly
(|Δ| = 0.000000)** under this account, which confirms both the diagnosis and
that the PHASE 11 classifier is byte-identical to PHASE 10's.

## 3. Magnitude

Per-row macro-F1 changes range from 0.000029 to 0.003853 (totalVI only).
Five-seed means, l2:

| Direction | quantity | published | corrected |
|---|---|---|---|
| A: PBMC10k→PBMC5k | scVI_matched | 0.660625 | 0.660625 |
| A | totalVI | 0.650438 | 0.652796 |
| A | **Delta** | **−0.010186** | **−0.007829** |
| B: PBMC5k→PBMC10k | scVI_matched | 0.605566 | 0.605566 |
| B | totalVI | 0.648312 | 0.648792 |
| B | **Delta** | **+0.042746** | **+0.043226** |

## 4. Effect on PHASE 10 conclusions

None of the qualitative findings change:

- Direction A Delta remains small and negative; Direction B Delta remains
  clearly positive. The directional asymmetry — the central PHASE 10 result and
  the motivation for PHASE 10B and PHASE 11 — is unchanged and marginally
  strengthened.
- The correction is well inside the between-seed spread (Delta SD ≈ 0.02–0.03),
  so no significance statement flips.
- PHASE 10B is **not** affected: `phase10b_source_size.py` returns
  `source_latent` (the adapted array) and never had the defect. Its conclusions
  stand as reported.

## 5. Fix

The return statement now passes the post-adaptation array:

```python
return {
    "source_latent": source_latent,
    ...
}
```

With this fix, fresh and resumed runs agree, and the in-code comment stating the
source-model embedding "is never used in a cross-space metric" becomes true.

## 6. What was regenerated

`scripts/10c_rescore_corrected.py` reloads the saved post-adaptation embeddings
and recomputes the affected metrics with PHASE 10's own functions
(`classify_transfer`, `neighbor_agreement`, `per_class_table`,
`paired_bootstrap_delta`). Nothing is retrained. Columns derived from the target
latent alone (`target_asw`, `target_knn_purity_15`, `target_leiden_ari`,
`target_leiden_nmi`) are carried over unchanged, since the defect could not
reach them.

Outputs, all suffixed `_corrected`, in `tables_corrected/`:
`cross_dataset_primary`, `cross_dataset_primary_delta`,
`cross_dataset_primary_delta_summary`, `cross_dataset_celltype_specific`,
`cross_dataset_neighbor_agreement`, `cross_dataset_generalization_gap`,
`cross_dataset_confusion`, plus `correction_changed_rows.csv` and
`correction_summary.json`.

PHASE 11 uses the corrected, §4-compliant coordinates throughout.
