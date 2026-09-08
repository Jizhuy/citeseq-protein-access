# PHASE 10 methods note

This experiment is **cross-dataset transfer / query generalization (scArches)**.
It is **not** strict zero-shot generalization: unlabeled target RNA (and the 14
observed proteins for totalVI) participate in query adaptation.

## PHASE 9 status

PHASE 9 was not empirically executed because native totalVI in scvi-tools 1.3.3
does not support arbitrary cell-level missing-protein masking without conflating
missingness with observed zeros or introducing unmatched-panel approximations.
That is a model/API limitation discovered during validation, not a failed
scientific result. PHASE 9 is not reopened here.

## Query API (scvi-tools 1.3.3)

Both `SCVI` and `TOTALVI` inherit `ArchesMixin`. Documented workflow:

1. Train the source model on source cells only.
2. `Model.prepare_query_anndata(query, reference_dir)`
3. `Model.load_query_data(query, reference_dir, accelerator="gpu", device="auto")`
4. `model.train(..., plan_kwargs={"weight_decay": 0.0})`
5. `model.get_latent_representation()` on query cells
6. `model.get_latent_representation(source)` on source cells

## Latent coordinate system

Both source and target cells are encoded by the **post-adaptation** model.
This matters: scArches leaves `l_encoder` trainable and its gradient hook keeps
gradients on every batch one-hot column of the first encoder layer, including
the source batch column. The source embedding therefore shifts during
adaptation, and an embedding taken from the pre-adaptation source model would
not be in the same coordinate system as the query embedding. Cross-space
metrics (source→target classifier, neighbor label agreement, dataset ASW) all
use the adapted-model embeddings. The pre-adaptation source embedding is saved
as `*_source_latent_source_model.npy` and used only to report how far source
cells moved (`adapted_source_drift_vs_source_model` in each manifest).

A true zero-shot mode without target parameter updates was **not** added:
`load_query_data` expands batch embeddings for unseen categories with random
padding. Embedding query cells without adaptation would use untrained batch
parameters.

## No-leakage rules

- Source models are not trained on target cells.
- Target labels are evaluation-only.
- Query early stopping uses unlabeled ELBO, not cell-type performance.
- Shared-class support thresholds (20 source / 10 target high-confidence cells)
  were pre-specified before reviewing model-specific results.

## Historical artifacts

New CUDA models are stored under `results/phase10_cross_dataset/` with `_cuda`
in the name. Historical Mac/MPS artifacts such as `scvi_matched_seed0` and
`totalvi_seed0` are not overwritten.
