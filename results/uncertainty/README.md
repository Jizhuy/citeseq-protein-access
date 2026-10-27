# Uncertainty and selective prediction

Two different uncertainty signals are reported separately:

- **Predictive uncertainty**: classifier-probability-based error detection and risk–coverage.
- **Latent posterior dispersion**: VAE posterior variance used as an error signal.

| Folder | Content | Script(s) |
| --- | --- | --- |
| `calibration/` | PBMC 20-seed reliability summary (predictive and latent error AUROC), full reliability metrics, corrected selective prediction, Lawlor reliability summary | `scripts/run_uncertainty_selective_prediction.py` |
| `selective_prediction/` | Macro-F1 by coverage | `scripts/summarize_selective_macro_f1.py` |
| `statistics/statistical_unit_ledger.csv` | Statistical unit used for each inference (cells, seeds, donors) | `scripts/run_uncertainty_selective_prediction.py` |

`calibration/pbmc_20seed_reliability_summary.csv` aggregates the frozen per-seed transfer table `results/transfer/replicate_runs/tables/transfer_replicates_20seed_primary.csv`. `scripts/build_pbmc_reliability_summary.py` regenerates it byte-identically.
