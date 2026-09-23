# Submission-revision reproducibility notes

Canonical scientific framing remains **protein-access vs residual totalVI-associated contrasts**, not a broad reliability benchmark.

## Key paths

| Path | Role |
| --- | --- |
| `manuscript/CITEseq_Protein_Access_Model_Associated_Gains.md` | Revised canonical manuscript source |
| `audit/MANUSCRIPT_METHODS_AUDIT.md` | Code↔manuscript audit |
| `audit/MANUSCRIPT_REVISION_PLAN.md` | Mandatory vs blocked items |
| `supplement/SUPPLEMENTARY_METHODS.md` | Full Methods detail |
| `supplement/SUPPLEMENTARY_RESULTS.md` | Sensitivities |
| `scripts/` | Deterministic new analyses |
| `tables/`, `figures/` | New outputs |

## New analyses in this revision

```bash
conda activate multiomics_robustness
# Lawlor donor summaries + class degradation + selective macro-F1 (fast)
python results/manuscript_submission_revision/scripts/summarize_class_degradation.py
python results/manuscript_submission_revision/scripts/lawlor_donor_summary_stats.py
python results/manuscript_submission_revision/scripts/summarize_selective_macroF1.py
python scripts/print_environment.py

# Bootstrap sensitivity (requires cached predictions)
python results/manuscript_submission_revision/scripts/cache_development_predictions.py
python results/manuscript_submission_revision/scripts/run_class_stratified_bootstrap.py

# Marker-channel dropout (PCA only; slow RNA HVG once)
python results/manuscript_submission_revision/scripts/marker_channel_dropout_pca.py
```

## Approximate compute

| Step | Notes |
| --- | --- |
| Cache predictions | ~5 min CPU (HVG + PCA + logreg) |
| Stratified bootstrap 10k | ~2–3 min after cache |
| Marker dropout | ~10–30 min (one RNA PCA + several protein PCAs) |
| Full VAE retrains | Not re-run in this revision |

## Seed policy

| Analysis | Seed |
| --- | --- |
| Development split / PCA | 0 |
| Primary bootstrap | 20260923 |
| Stratified bootstrap | 20260924 |
| Transfer deep models | 0–19 |

## Do not

- Treat bootstrap CIs as biological CIs  
- Treat 38/50 or 19/20 as donor *n*  
- Claim causal architecture effects from \(G_{\mathrm{assoc}}\)  
- Invent GEO accessions for PBMC example files
