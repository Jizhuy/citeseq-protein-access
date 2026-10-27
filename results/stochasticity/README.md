# Training stochasticity

Replicated crossed design: source subset × model seed × run replicate, for scVI_matched and totalVI under transfer.

| Folder | Content | Script(s) |
| --- | --- | --- |
| `replicate_runs/tables/` | Crossed design, replicated variance components, per-cell and cell-type variance, power check | `scripts/run_training_stochasticity.py`, `scripts/analyze_training_stochasticity.py` |
| `variance_decomposition/statistics/` | Corrected variance fractions (source subset / model seed / subset × seed interaction / residual; run replicates estimate the residual) for each model and for the totalVI − scVI_matched delta | `scripts/run_variance_decomposition.py` |
| `variance_decomposition/logs/` | Variance-model audit | same |
