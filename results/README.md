# Frozen results

All files under `results/` are frozen outputs copied byte-for-byte from the analysis runs that support the manuscript (file names were neutralised; file contents were not changed, apart from one machine-specific path normalised in `donor_validation/lawlor_runs/logs/data_prep_summary.json`). No value was recomputed during migration.

| Folder | Content | Status |
| --- | --- | --- |
| `primary/` | Strict train-only internal analysis | **PRIMARY** |
| `transfer/` | Bidirectional PBMC10k ↔ PBMC5k transfer (20 seeds) and the feature-access-matched transductive concat control | Complete |
| `degradation/` | Training-time protein degradation and post-adaptation test-time protein corruption | Complete |
| `uncertainty/` | Predictive uncertainty, latent posterior dispersion, selective prediction, statistical-unit ledger | Complete |
| `stochasticity/` | Replicated crossed design and corrected variance decomposition | Complete |
| `donor_validation/` | Lawlor donor-level external validation (10 donors) | Complete |
| `supplementary/` | Simple-baseline controls, marker-channel dropout (PCA only), runtime environment, and the **secondary** information-access sensitivity | Supplementary |

A value-by-value index (value → file → script → manuscript figure/table) is in [`../docs/RESULT_MANIFEST.md`](../docs/RESULT_MANIFEST.md).

Pipelines that use `config/config.yaml` write fresh outputs to `results/tables/`, `results/embeddings/`, `results/models/` and `results/logs/`. These intermediates are not shipped and are git-ignored; frozen copies live in the analysis folders above.
