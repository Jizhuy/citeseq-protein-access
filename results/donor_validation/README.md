# Lawlor donor-level external validation

- Cohort: Lawlor et al. PBMC CITE-seq, Baseline condition, Demuxlet singlets. HCA `efea6426-510a-4b60-9a19-277e52bfa815`; ENA `PRJEB40376` / `ERP124005`.
- **Biological unit: 10 donors.** Inference uses donor-level paired deltas (exact two-sided sign test).
- Labels are protein-gated. This favours protein-aligned baselines.

| Folder | Content | Script(s) |
| --- | --- | --- |
| `lawlor_dataset_audit/` | Provenance, independence audit, gene/protein overlap, label harmonisation | dataset audit (provenance in `data_provenance.md`) |
| `lawlor_runs/tables/`, `lawlor_runs/logs/` | Donor-held-out runs, folds, replication summary, data-preparation provenance | `scripts/run_lawlor_donor_validation.py`, `scripts/analyze_lawlor_donor_validation.py` |
| `statistics/` | Donor-level deltas, per-donor metrics, pairwise consistency | `scripts/run_lawlor_donor_level.py`, `scripts/run_lawlor_pairwise_consistency.py` |
| `summary/` | All-method donor-level summary, sign test, PCA per donor | `scripts/run_lawlor_pca_per_donor.py`, `scripts/summarize_lawlor_donor_statistics.py` |
| `logs/` | Macro-F1 definition and inconsistency audit | — |

The 38/50 latent-error-signal consistency in `statistics/lawlor_pairwise_consistency_corrected.csv` is a **computational** consistency count over 50 matched fold × seed pairs (5 donor-held-out folds × 10 seeds). It is not a biological replication count; the biological unit remains the 10 donors.

Not run: a second external, non-protein-gated validation cohort.
