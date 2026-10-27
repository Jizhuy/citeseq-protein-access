# Protein degradation

Two distinct degradation analyses:

1. **Training-time protein degradation** (`training_time/`). Protein counts were masked before totalVI training, using the historical feature-transductive development pipeline, with 5 mask replicates per level. Script: `scripts/run_protein_degradation.py`. It needs the reference totalVI fit from `scripts/run_feature_transductive_models.py`. `training_time_degradation_summary.csv` summarises `protein_sparsity_primary_delta.csv` (see `DEGRADATION_REPLICATE_AUDIT.md`). `scripts/build_training_time_degradation_summary.py` regenerates it byte-identically.
2. **Post-adaptation test-time protein corruption** (`test_time_corruption/`). Target protein counts were corrupted after scArches query adaptation (20 seeds). Script: `scripts/run_test_time_protein_corruption.py`.

`class_level/` holds class-level degradation summaries (`scripts/summarize_class_degradation.py`).

Not run: totalVI retraining under complete marker-channel dropout. The marker-channel dropout sensitivity in `results/supplementary/marker_channel_dropout/` is PCA-only.
