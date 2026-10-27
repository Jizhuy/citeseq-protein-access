# Bidirectional cross-dataset transfer

- Direction A: PBMC10k (source) → PBMC5k (target).
- Direction B: PBMC5k (source) → PBMC10k (target).
- Deep models: scVI_matched and totalVI with scArches query adaptation, 20 seeds per direction.

scArches query adaptation uses **unlabeled target features**, so it is transductive adaptation. The **feature-access-matched transductive concat PCA** control fits PCA on source and target features together. It matches that target-feature access, but not the learning objective or the architecture. The **inductive concat PCA** baseline is fitted on source cells only.

| Folder | Content | Script(s) |
| --- | --- | --- |
| `replicate_runs/tables/` | 20-seed transfer summaries, per-seed deltas, crossed design, retention/imbalance, selective prediction, variance components | `scripts/run_bidirectional_transfer.py`, `scripts/analyze_bidirectional_transfer.py`, `scripts/analyze_crossed_design.py` |
| `matched_concat_control/tables/matched_transfer_control.csv` | Inductive vs feature-access-matched transductive concat PCA, alongside the 20-seed deep-model means | `scripts/run_matched_concat_transfer.py` |
| `matched_concat_control/logs/` | Control design and integrity checks | same |
