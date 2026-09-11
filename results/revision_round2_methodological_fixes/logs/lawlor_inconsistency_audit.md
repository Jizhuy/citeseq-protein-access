# Lawlor inconsistency audit: 50 vs 19/20

## Design reminder

Lawlor formal validation = 5 donor folds × 10 seeds × 2 models = **100 runs**.
Matched comparisons = **50** fold×seed pairs.

## Canonical counts (from `phase12b_primary_runs.csv`)

Computed in `statistics/lawlor_pairwise_consistency_corrected.csv`.

Headline historical report values (phase12b_report.md) already used fold×seed pairs:

| Metric | Historical report | Nature |
|---|---|---|
| macro-F1 | **50/50** totalVI higher | Lawlor |
| predictive AUROC | **45/50** | Lawlor |
| latent AUROC | **38/50** | Lawlor |
| NLL / Brier / AURC | **50/50** totalVI better | Lawlor |
| ECE | mixed (~29/50) | Lawlor |

## Where did **19/20** come from?

**Not from Lawlor.** It is an **internal PBMC transfer** sign-consistency count:

1. **Direction B macro-F1 (PBMC5k→PBMC10k), 20-seed grid:** totalVI favored in **19/20** seeds  
   — see `results/phase13_manuscript/manuscript/manuscript_v2_claim_checked.md`,  
   `results/phase14_submission_prep/manuscript/manuscript_genome_biology_v1.md`,  
   figure legends for cross-dataset transfer.

2. **Latent error AUROC in PBMC transfer:** also written as **19/20 seeds per direction**  
   — same manuscript family (PHASE 11B transfer, not Lawlor).

3. Lawlor latent sign consistency was correctly reported elsewhere as **38/50**, not 19/20.

## Root cause of the apparent inconsistency

The **19/20** figure is an **internal-transfer seed-consistency** number that coexists in the same manuscripts as Lawlor **50/50** / **38/50**.  
If prose or a table ever juxtaposed 19/20 with Lawlor without naming the experiment, that is a **cross-experiment copy/context error**, not a miscount of the Lawlor 50 pairs.

No Lawlor retraining is required for this correction.

## Action

- Keep Lawlor pairwise denominators as **/50**.
- Keep PBMC 20-seed denominators as **/20**.
- Never mix them in a single sentence without naming the experiment.
