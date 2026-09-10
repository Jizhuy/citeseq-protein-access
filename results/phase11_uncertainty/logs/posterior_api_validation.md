# PHASE 11 — posterior API validation (scvi-tools 1.3.3)

Written before any uncertainty metric was computed, per PHASE 11 §6.
Environment: RTX 4060 Laptop, CUDA 12.6, torch 2.14.0+cu126, scvi-tools 1.3.3,
Python 3.11.7.

## 1. What was inspected

The installed implementation was read directly rather than assumed from a
tutorial. `SCVI.get_latent_representation` and `TOTALVI.get_latent_representation`
both resolve to the same inherited implementation:

```
scvi.model.base.VAEMixin.get_latent_representation(
    adata=None, indices=None, give_mean=True, mc_samples=5000,
    batch_size=None, return_dist=False, dataloader=None,
) -> npt.NDArray | tuple[npt.NDArray, npt.NDArray]
```

MRO confirms the shared source:

- `SCVI` → `EmbeddingMixin, RNASeqMixin, VAEMixin, ArchesMixin, ...`
- `TOTALVI` → `RNASeqMixin, VAEMixin, ArchesMixin, ...`

## 2. How q(z|x) is represented

From the installed source of `VAEMixin.get_latent_representation`:

```python
outputs = self.module.inference(**self.module._get_inference_input(tensors))
if MODULE_KEYS.QZ_KEY in outputs:
    qz   = outputs.get(MODULE_KEYS.QZ_KEY)
    qzm  = qz.loc
    qzv  = qz.scale.square()
else:
    qzm = outputs.get(MODULE_KEYS.QZM_KEY)
    qzv = outputs.get(MODULE_KEYS.QZV_KEY)
    qz  = Normal(qzm, qzv.sqrt())
if return_dist:
    qz_means.append(qzm.cpu()); qz_vars.append(qzv.cpu()); continue
...
if return_dist:
    return torch.cat(qz_means).numpy(), torch.cat(qz_vars).numpy()
```

Conclusions:

- The encoder returns a `torch.distributions.Distribution` under `QZ_KEY`.
- `return_dist=True` returns the tuple **(posterior mean, posterior variance)**,
  where variance is `qz.scale ** 2`. This is a closed-form variance, not a
  sampled estimate.
- When `return_dist=True`, `give_mean` and `mc_samples` are ignored entirely.
- Both PHASE 10 models were built with `latent_distribution="normal"`
  (`MATCHED_HYPERPARAMETERS` and `DEFAULT_HYPERPARAMETERS`, both `n_latent=20`),
  so `qz` is a diagonal Gaussian. Posterior variance is well defined and
  Gaussian differential entropy would be mathematically appropriate. The `"ln"`
  (logistic-normal) branch that would require Monte Carlo is **not** taken.

## 3. Direct probe on a real saved model

`SCVI.load` on
`results/phase10_cross_dataset/models/pbmc5k_to_pbmc10k/scvi_matched/seed0_cuda`
with the PBMC5k source AnnData:

| property | value |
|---|---|
| loaded | True |
| `latent_distribution` | `normal` |
| `n_latent` | 20 |
| `return_dist=True` supported | True |
| `qzm` shape | (3994, 20) |
| `qzv` shape | (3994, 20) |
| all `qzv > 0` | True |
| `qzv` mean / min / max | 0.07336 / 0.000601 / 0.41241 |
| `qzm` equals `give_mean=True` output | True |
| CUDA peak allocated | 161.7 MiB |

So the API itself is fully sufficient: **direct posterior variance is available,
no Monte Carlo sampling is required.**

## 4. The blocking problem: the adapted model was not persisted

PHASE 11 §4 requires that source and target both be encoded by the
**post-adaptation** model. Each PHASE 10 seed directory contains exactly one
model:

```
model.pt
query_training_history.csv
source_training_history.csv
```

That `model.pt` is the **pre-adaptation source reference model**, saved inside
`train_scvi_source` / `train_totalvi_source`. In `query_adapt`, the scArches
model is built, trained for up to 200 epochs, its latents are written, and then:

```python
del qmodel, query
```

The adapted encoder weights were never written to disk. All 20 `model.pt` files
(2 directions × 2 models × 5 seeds) are source reference models.

### 4.1 The reference embedding is genuinely frozen

A useful and verified property: the saved "post-adaptation" source latent is
**byte-identical** to the pre-adaptation source latent.

- `np.array_equal(source_latent.npy, source_latent_source_model.npy)` → `True`
- max abs difference → `0.0`
- PHASE 10's own manifest recorded `mean_euclidean_shift: 0.0` at run time.

This is correct scArches behaviour: scvi-tools masks gradients so only the
query batch's one-hot columns update, leaving the reference embedding fixed.
(The PHASE 10 source comment asserting the source embedding "shifts during
adaptation" is misleading, but the recorded diagnostic was honest and the
artifacts are internally consistent.)

Consequence: **source-cell posterior variance in the exact PHASE 10 coordinate
system is recoverable** from the saved reference model. The probe confirmed
`source_model_matches_saved_adapted_latent = True`.

### 4.2 Target-cell posterior is not recoverable

Rebuilding the query model without training and encoding the target fails:

```
SCVI.prepare_query_anndata(query, model_dir)     # 100.0% reference vars found
qmodel = SCVI.load_query_data(query, model_dir)  # builds fine
qmodel.get_latent_representation()
RuntimeError: Trying to query inferred values from an untrained model.
              Please train the model first.
```

scvi-tools refuses inference from an untrained query model, so there is no
substitute encoding. Recovering the target posterior would require re-running
query adaptation, which is (a) retraining, explicitly forbidden, and (b) not
reproducible — the RNG state entering query adaptation depends on the whole
preceding source-training run, so a rerun would produce different latents and
would break the PHASE 10 reproduction requirement of §11.

## 5. Determination under §7

The §7 stop condition is invoked, on the third branch:

- Direct posterior variance **is** available from the API — verified, not assumed.
- Reproducible posterior sampling is likewise available in principle.
- But neither can be applied to **target** cells in the validated
  post-adaptation coordinate system, because the required artifact does not
  exist and cannot be regenerated without retraining.

Therefore:

- **Latent posterior uncertainty for target cells is marked UNSUPPORTED.**
  No value is fabricated. `error_detection_auroc_latent`,
  `error_detection_auprc_latent` and `median_latent_uncertainty` are written as
  `NA` in `uncertainty_primary.csv` with this document as the recorded reason.
- Monte Carlo posterior sampling (§9, §25) is **not** performed: it would
  inherit the same missing-artifact problem, so S=50 is not used and no
  sample-count sensitivity check is run.
- **Source-cell** latent posterior uncertainty **is** extracted and reported as
  a descriptive secondary result, since it is exactly valid in the PHASE 10
  coordinate system (§4.1). It cannot be compared against target latent
  uncertainty, so the latent arm of §26 is reported as source-only.
- PHASE 11 proceeds on the three remaining concepts, exactly as §7 directs:
  classifier predictive uncertainty (B), model-seed ensemble uncertainty (C),
  and PHASE 10B source-subset instability (D).

## 6. Recommendation for future phases

Persisting the adapted query model (`qmodel.save(...)`) costs roughly the same
as the reference model (~80 MB for scVI, ~68 MB for totalVI; ~1.5 GB for a full
20-model grid) and would make latent posterior uncertainty available. This is
recorded as a limitation of PHASE 10's artifact policy, not a limitation of
scvi-tools 1.3.3.
