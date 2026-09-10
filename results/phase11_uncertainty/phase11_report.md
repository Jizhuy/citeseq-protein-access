# PHASE 11 — formal uncertainty characterization and calibration

Answers A–AM. No representation model was trained, modified, or retrained.
Primary endpoint and all cutoffs were pre-specified before results were seen.

Headline: **predictive uncertainty is strongly discriminative and reasonably
well calibrated in both transfer directions, and it predicts which target cells
are fragile to source-population resampling. But totalVI does not provide
reliably better uncertainty than scVI_matched once model-seed variability is
accounted for.**

---

## A. Posterior uncertainty API / extraction method

`scvi.model.base.VAEMixin.get_latent_representation(..., return_dist=True)`,
shared by `SCVI` and `TOTALVI`. It calls `module.inference()` and returns
`(qz.loc, qz.scale.square())` — closed-form posterior mean and variance. Both
models were built with `latent_distribution="normal"` and `n_latent=20`, so
`qz` is a diagonal Gaussian and the `"ln"` Monte-Carlo branch is never taken.
Verified by direct probe on a saved model: `qzv` shape (3994, 20), all
positive, mean 0.0734. Full record in `logs/posterior_api_validation.md`.

## B. Was direct posterior variance available?

**Yes from the API, no for the cells that matter.** The API exposes it exactly.
But PHASE 10 persisted only the *pre-adaptation source reference* model; the
scArches-adapted query model — which defines the post-adaptation coordinate
system used by every PHASE 10 metric — was deleted after its latents were
written. Rebuilding it untrained and encoding the target fails with
`RuntimeError: Trying to query inferred values from an untrained model`, and
training it would be retraining and non-reproducible.

Per §7, **target-cell latent posterior uncertainty is marked UNSUPPORTED**. No
value was fabricated; `error_detection_auroc_latent`,
`error_detection_auprc_latent` and `median_latent_uncertainty` are `NA` with
`latent_uncertainty_status` recording the reason.

Source-cell posterior variance *was* extracted for all 20 models. It is exactly
valid in the PHASE 10 coordinate system for **scVI only** — scArches freezes the
reference embedding (drift 0.000000) — while totalVI reference models define the
pre-adaptation space and are flagged `valid_in_phase10_coordinate_system=False`.

| direction | model | mean posterior variance | median U_latent (source) | valid |
|---|---|---|---|---|
| A | scVI_matched | 0.1282 | 0.1284 | yes |
| A | totalVI | 0.0236 | 0.0234 | no (pre-adaptation space) |
| B | scVI_matched | 0.0807 | 0.0795 | yes |
| B | totalVI | 0.0316 | 0.0307 | no (pre-adaptation space) |

## C / D. Posterior sampling and MC sample count

**Not used.** Direct variance is available, so sampling was unnecessary; and it
would have inherited the identical missing-artifact problem for target cells.
S=50 was not drawn and no S=50 vs S=100 check was run. §25 (posterior-latent
predictive uncertainty) is consequently not reportable.

## E. PHASE 10 reproduction check

**All 40 cells reproduce exactly (|Δ| = 0.000000)** — 2 directions × 2 models ×
5 seeds × 2 label levels, tolerance 1e-6.

Reaching that required finding and fixing a real PHASE 10 defect. `run_one_model`
returned `src["latent"]` (pre-adaptation source) to the evaluator while writing
the post-adaptation array to disk, so freshly-trained runs scored **pre-adaptation
source against post-adaptation target** — the exact mixing §4 forbids. scVI was
numerically immune (drift 0.0); totalVI seeds 1–4 were affected (drift
0.015–0.020). Seed 0 was resumed from saved artifacts and so was already correct.

The bug was patched, PHASE 10 was re-scored from saved embeddings into
`results/phase10_cross_dataset/tables_corrected/` (no retraining, historical
files untouched), and 16 rows changed. Conclusions are unchanged:

| l2 Delta (totalVI − scVI), 5-seed mean | published | corrected |
|---|---|---|
| Direction A | −0.0102 | −0.0078 |
| Direction B | +0.0427 | +0.0432 |

Details: `results/phase10_cross_dataset/logs/phase10_erratum_latent_alignment.md`.

## F–I. Primary endpoint — l2 error-detection AUROC using U_pred = 1 − p_max

Mean ± SD over the five model seeds (seed 0 in brackets):

| | scVI_matched | totalVI |
|---|---|---|
| **F/G. Direction A** (PBMC10k→PBMC5k) | **0.8451 ± 0.0118** [0.8270] | **0.8472 ± 0.0134** [0.8682] |
| **H/I. Direction B** (PBMC5k→PBMC10k) | **0.8575 ± 0.0113** [0.8627] | **0.8553 ± 0.0285** [0.8733] |

All four are far above chance: uncertainty ranks errors well in every condition.
The two models are indistinguishable on this endpoint.

## J. AUPRC and error prevalence

| direction | model | error prevalence | AUPRC (mean) | AUPRC / prevalence |
|---|---|---|---|---|
| A | scVI_matched | 0.203 | 0.526 | 2.59× |
| A | totalVI | 0.186 | 0.506 | 2.72× |
| B | scVI_matched | 0.237 | 0.610 | 2.58× |
| B | totalVI | 0.213 | 0.570 | 2.68× |

AUPRC is 2.6–2.7× the prevalence baseline in every cell.

## K. Latent posterior uncertainty vs error

**Not evaluable for target cells** (see B). Error-detection AUROC/AUPRC against
`U_latent`, the Spearman correlation with true-class negative log probability,
and the correct-vs-incorrect latent distributions are all `NA` with a documented
reason rather than estimated from a substitute encoding.

## L. NLL / Brier / ECE for every model and direction

Mean over five seeds:

| direction | model | NLL | Brier | ECE (15 eq-freq) | ECE (15 eq-width) | T | ECE after T |
|---|---|---|---|---|---|---|---|
| A | scVI_matched | 0.5600 | 0.2844 | 0.0271 | 0.0298 | 0.933 | 0.0358 |
| A | totalVI | 0.5169 | 0.2628 | 0.0220 | 0.0216 | 0.881 | 0.0416 |
| B | scVI_matched | 0.7122 | 0.3279 | 0.0354 | 0.0347 | 0.878 | 0.0550 |
| B | totalVI | 0.5608 | 0.2903 | 0.0399 | 0.0403 | 0.812 | 0.0437 |

Equal-frequency and equal-width ECE agree closely, so the binning choice does
not drive the conclusion. ECE of 0.022–0.040 is genuinely small for a 21–22-class
problem under domain shift.

Temperature scaling **did not help** and mostly hurt: every fitted T < 1
(0.81–0.93), meaning the source-holdout wanted *sharper* probabilities, but
applying that to shifted target data increased ECE in three of four cells. The
uncalibrated probabilities are the better ones, and are the ones reported as
primary. This is an honest negative result about post-hoc calibration
transferring across a domain shift.

## M. Reliability-diagram findings

(`figures/phase11_figure2_reliability.*`) Both models sit slightly **below** the
diagonal through the mid-confidence range — mild overconfidence — and converge
onto it above ~0.9 confidence, where most cells live. Bin occupancy is even by
construction (15 equal-frequency bins, 397–398 cells each in Direction B). The
deviation is modest and there is no severe miscalibration in either direction.
In Direction A totalVI tracks the diagonal marginally better; in Direction B it
is visibly worse in the 0.75–0.85 band, which is what drives its higher ECE.

## N. Risk-coverage / AURC

Mean AURC (lower better): A scVI 0.0575, A totalVI 0.0511; B scVI 0.0684,
B totalVI 0.0606. totalVI is nominally better in both, but see AB — the margin
is inside the seed spread.

## O. Does rejecting high-uncertainty cells help?

**Yes, substantially and monotonically for accuracy.** Mean over five seeds:

| coverage | A scVI acc / macro-F1 | A totalVI | B scVI | B totalVI |
|---|---|---|---|---|
| 100% | 0.797 / 0.661 | 0.814 / 0.653 | 0.763 / 0.606 | 0.787 / 0.649 |
| 90% | 0.839 / 0.683 | 0.856 / 0.665 | 0.812 / 0.617 | 0.833 / 0.659 |
| 70% | 0.909 / 0.715 | 0.923 / 0.667 | 0.896 / 0.639 | 0.910 / 0.671 |
| 50% | 0.970 / 0.747 | 0.975 / 0.736 | 0.965 / 0.660 | 0.964 / 0.704 |

Abstaining on the most uncertain half lifts accuracy from ~0.78 to ~0.97.
Macro-F1 also improves but far more slowly (0.661→0.747), because rejection
removes disproportionately many rare-class cells — accuracy gains overstate the
biological benefit.

Error enrichment (pre-specified strata) is extreme: in Direction A the bottom-20%
uncertainty stratum contains essentially zero errors, giving enrichment ratios of
372× to ∞. Reported as ∞ where the denominator is exactly zero, not silently
clipped.

## P. Five-seed ensemble results

| direction | model | ens. macro-F1 | mean seed macro-F1 | ens. ECE | mean seed ECE | ens. error AUROC | ens. AURC |
|---|---|---|---|---|---|---|---|
| A | scVI_matched | 0.6873 | 0.6606 | 0.0230 | 0.0271 | 0.8508 | 0.0476 |
| A | totalVI | 0.6595 | 0.6528 | 0.0217 | 0.0220 | 0.8699 | 0.0424 |
| B | scVI_matched | 0.6121 | 0.6056 | 0.0329 | 0.0354 | 0.8801 | 0.0574 |
| B | totalVI | 0.6839 | 0.6488 | 0.0859 | 0.0399 | 0.8872 | 0.0380 |

Ensembling improves macro-F1 in all four cells (+0.007 to +0.035) and improves
error-detection AUROC and AURC over the seed mean. It is **not** uniformly
superior: Direction B totalVI ensemble ECE is 0.0859, more than double its mean
single-seed ECE of 0.0399. Averaging five disagreeing totalVI seeds produces
flatter probabilities that are then underconfident.

## Q. Ensemble disagreement vs prediction error

Ensemble disagreement is a good but not superior error detector: AUROC 0.8119
(A scVI), 0.8269 (A totalVI), 0.8409 (B scVI), **0.7014 (B totalVI)**. Variation
ratio is weaker still (0.731–0.772). Both are below the simple `U_pred` baseline
in every cell, so the extra cost of five seeds buys no better error *ranking*
than one model's own confidence.

The B/totalVI outlier is informative: totalVI's Direction B seeds have the
highest median disagreement (0.0685 vs 0.0282 for scVI) yet the *least*
error-informative disagreement. Its seeds disagree a lot, and largely for
reasons unrelated to correctness.

Median variation ratio is 0.0 in all four cells — for most target cells all five
seeds agree on the label.

## R. Individual-seed vs ensemble calibration

Ensembling improved ECE for A scVI (0.0271→0.0230), left A totalVI unchanged
(0.0220→0.0217), improved B scVI (0.0354→0.0329) and markedly degraded B totalVI
(0.0399→0.0859). So ensembling is a reliable way to improve *accuracy* but not a
reliable way to improve *calibration*.

## S. Source vs target uncertainty shift

Source uncertainty computed **out of sample** by 5-fold stratified CV on source
cells; in-sample source probabilities would be optimistically biased.

| direction | model | median U source | IQR | median U target | IQR |
|---|---|---|---|---|---|
| A | scVI_matched | 0.0294 | 0.146 | 0.0992 | 0.336 |
| A | totalVI | 0.0277 | 0.132 | 0.0915 | 0.272 |
| B | scVI_matched | 0.0994 | 0.309 | 0.1339 | 0.351 |
| B | totalVI | 0.0825 | 0.286 | 0.1490 | 0.326 |

Target uncertainty is higher than source in all four cells, and the *relative*
jump is much larger in Direction A (~3.4×) than Direction B (~1.4–1.8×). Stated
descriptively: this is not by itself evidence of useful domain-shift awareness.

## T. Rare-class uncertainty findings

Uncertainty tracks difficulty across the l2 classes. Direction A totalVI, sorted
by median `U_pred`, runs from CD8 TCM (0.500, error 0.765), Treg (0.416, 0.797),
CD4 CTL (0.406, 0.829), gdT (0.396, 0.930), ILC (0.358, 0.929) and dnT (0.342,
error 1.000) at the top, down to NK (0.004, error 0.004), CD14 Mono (0.007,
0.007) and pDC (0.026, 0.000) at the bottom.

Every historically unstable class — gdT, ILC, NK_CD56bright, Treg, CD8 TCM,
Eryth — is in the high-uncertainty tail with high error rates, so the model
correctly signals low confidence on exactly the populations PHASE 10 flagged. It
signals them but still gets them wrong: uncertainty is diagnostic, not curative.
All eligible classes are reported in
`tables/uncertainty_celltype_specific.csv`; none were cherry-picked.

## U. Source class support vs uncertainty and error

Spearman across eligible l2 classes:

| direction | model | vs F1 | vs error rate | vs median U_pred | vs disagreement |
|---|---|---|---|---|---|
| A | scVI_matched | +0.318 (p=0.15) | −0.379 (0.082) | −0.214 (0.34) | −0.220 (0.32) |
| A | totalVI | +0.447 (0.037) | −0.532 (0.011) | −0.457 (0.033) | −0.416 (0.054) |
| B | scVI_matched | +0.404 (0.069) | −0.466 (0.033) | −0.396 (0.076) | −0.546 (0.011) |
| B | totalVI | +0.435 (0.049) | −0.547 (0.010) | −0.566 (0.008) | −0.486 (0.026) |

All eight support-vs-uncertainty and support-vs-error correlations have the
expected sign: rarer source classes are more uncertain and more often wrong. The
relationship is consistently stronger for totalVI. Association only — no causal
claim; support is confounded with intrinsic class separability.

## V. PHASE 10B source-subset instability

19 common l2 classes, 3291 shared PBMC5k target cells, five independent PBMC10k
source subsamples (n=3994), model seed 0. Median variation ratio is 0.0 for both
models — most cells keep their label under source resampling — while median
true-class probability SD is 0.0570 (scVI) and 0.0732 (totalVI). Instability is
concentrated in a minority of cells rather than spread evenly.

## W. Full-model uncertainty vs subset instability

**This is the strongest result in PHASE 11.** Spearman correlation between
full-source Direction A uncertainty and PHASE 10B instability (n=3291, all
p < 1e-300):

| full-source signal | vs subset variation ratio | vs true-class prob SD | vs subset entropy |
|---|---|---|---|
| scVI U_pred | 0.563 | 0.755 | 0.868 |
| scVI ensemble disagreement | 0.582 | 0.771 | 0.890 |
| totalVI U_pred | 0.627 | 0.731 | 0.859 |
| totalVI ensemble disagreement | 0.584 | 0.733 | 0.860 |

Pre-specified stratification (§33):

| model | mean variation ratio top-20% / bottom-20% | mean true-prob SD top / bottom |
|---|---|---|
| scVI_matched | 0.2261 / 0.0006 (372×) | 0.1249 / 0.0057 (22×) |
| totalVI | 0.2666 / 0.0000 (∞) | 0.1483 / 0.0073 (20×) |

A single full-source model's own confidence already identifies which target
cells would change prediction if the source cohort were resampled. Cells it is
confident about are essentially immune; cells it is unsure about are where the
source-subset instability that motivated PHASE 11 lives.

## X. Neighbor stability result

Mean pairwise Jaccard overlap of k=15 source-neighbor sets across the five seeds
is **low in absolute terms**: 0.171 (A scVI), 0.138 (A totalVI), 0.297 (B scVI),
0.197 (B totalVI). Independently trained seeds place target cells next to largely
different source cells even when they agree on the label.

**totalVI neighborhoods are consistently less stable than scVI's in both
directions**, despite PHASE 10 finding totalVI had *higher* neighbor label
agreement.

## Y. Neighbor stability vs prediction error

**Essentially no relationship.** Spearman(stability, error) is −0.042, +0.021,
−0.115, +0.007 across the four cells; AUROC for low stability detecting errors is
0.530, 0.485, 0.581, 0.495 — at or near chance. Mean stability for correct vs
incorrect cells is nearly identical (e.g. A scVI 0.1711 vs 0.1638). Correlations
with `U_pred` (−0.062 to +0.031) and with ensemble disagreement (−0.196 to
+0.046) are similarly negligible.

Geometric neighborhood reproducibility and predictive reliability are
**different things**. Neighbor instability is not a usable error signal.

## Z. Direction A vs Direction B

- Error-detection AUROC is slightly higher in B (0.855–0.858 vs 0.845–0.847).
- Calibration is worse in B for both models (NLL 0.56–0.71 vs 0.52–0.56; ECE up
  to 0.040 vs 0.027).
- Ensemble disagreement is higher in B for totalVI (0.0685) and less informative
  there (AUROC 0.701).
- Rare-class support correlations are stronger in B.
- Source→target uncertainty inflation is relatively larger in A.

So neither direction is uniformly harder: **B is harder to calibrate but not
harder to rank**. Uncertainty does *not* mirror the macro-F1 asymmetry — totalVI
gains ~+0.043 macro-F1 in B and loses ~0.008 in A, yet its uncertainty quality is
no better than scVI's in either.

## AA. scVI_matched vs totalVI uncertainty comparison

totalVI has consistently lower NLL (−0.043 in A, −0.151 in B) and Brier (−0.022,
−0.038) across seeds, so its probability *vectors* are better. On the pre-specified
primary endpoint — error-detection AUROC — the two are indistinguishable
(+0.0021 ± 0.0229 in A, −0.0022 ± 0.0275 in B). ECE shows no consistent winner.

## AB. Paired bootstrap confidence intervals — and a necessary caveat

500 paired target-cell resamples with identical indices, seed 1111, conditional
on model seed 0. Reporting those CIs alone would be misleading, so across-seed
spread is carried in the same table:

| direction | metric | seed-0 Δ | bootstrap 95% CI | across-seed Δ (mean ± SD) | evidence |
|---|---|---|---|---|---|
| A | error AUROC | +0.0412 | [+0.0245, +0.0567] | +0.0021 ± 0.0229 | seed-dependent |
| A | error AUPRC | +0.0408 | [−0.0159, +0.0961] | −0.0192 ± 0.0414 | no reliable difference |
| A | ECE | −0.0102 | [−0.0217, +0.0042] | −0.0051 ± 0.0121 | no reliable difference |
| A | Brier | −0.0443 | [−0.0554, −0.0344] | −0.0216 ± 0.0193 | **robust** |
| A | AURC | −0.0197 | [−0.0242, −0.0156] | −0.0064 ± 0.0083 | seed-dependent |
| B | error AUROC | +0.0106 | [+0.0002, +0.0216] | −0.0022 ± 0.0275 | seed-dependent |
| B | error AUPRC | +0.0066 | [−0.0275, +0.0458] | −0.0400 ± 0.0520 | no reliable difference |
| B | ECE | +0.0176 | [+0.0060, +0.0267] | +0.0045 ± 0.0182 | seed-dependent |
| B | Brier | −0.0109 | [−0.0184, −0.0033] | −0.0376 ± 0.0389 | seed-dependent |
| B | AURC | −0.0034 | [−0.0068, −0.0004] | −0.0078 ± 0.0217 | seed-dependent |

**Only one of ten comparisons is robust** (Direction A Brier, favouring totalVI).
Six have bootstrap CIs excluding zero while the five seeds disagree in sign.
Between-model-seed variability dominates target-cell sampling variability for
these metrics, exactly the §38 hazard. Any claim resting on a single seed's
bootstrap would be an artifact.

## AC. CUDA / runtime / RAM

Peak CUDA allocated 161.7 MiB, peak reserved 162.0 MiB — trivial against the
8188 MiB budget, because PHASE 11 is inference and aggregation only. Peak process
RSS 1105–1624 MiB against ~7.6 GiB. Models were loaded strictly one at a time,
with `gc.collect()` and `torch.cuda.empty_cache()` between loads, and
`reset_peak_memory_stats()` before each. Approximate wall times: gates ~25 s
each; PHASE 10 re-scoring ~9 min (dominated by the O(n²) silhouette); main stage
~4 min; subset ~13 s; latent (20 model loads) ~30 s; figures ~27 s.

## AD. Warnings and limitations

1. **Target latent posterior is unavailable** (§7 invoked), so the four-way
   taxonomy is realised three ways. Concept A is source-only and scVI-only.
2. **A PHASE 10 defect was found and corrected here.** PHASE 11 runs on corrected
   coordinates; PHASE 10's published totalVI seed 1–4 numbers need the erratum.
3. **Between-seed variability dominates** most scVI-vs-totalVI comparisons. Five
   seeds is few for stable metric-level inference.
4. **Temperature scaling degraded calibration** in three of four cells; post-hoc
   calibration fitted on source does not transfer cleanly under this shift.
   T was also estimated from a 75%-data companion classifier, a slight mismatch.
5. **Source-subset instability uses model seed 0 only** — PHASE 10B trained one
   model seed per subsample, so concepts C and D cannot be fully crossed.
6. Ensemble disagreement is *not* exact Bayesian mutual information, and
   model-seed and source-subset spread are not posterior uncertainty.
7. §25 and the S=50/S=100 check are not reportable; §49 and §50 were not run.
8. Error-enrichment ratios are ∞ where the bottom-20% stratum has zero errors.
9. `balanced_accuracy`/`f1` emit `y_pred contains classes not in y_true`
   warnings when a rare class is never predicted; zero_division=0 throughout.

## AE. Files created

Tables (`results/phase11_uncertainty/tables/`): `uncertainty_primary.csv`,
`uncertainty_ensemble.csv`, `uncertainty_celltype_specific.csv`,
`uncertainty_source_subset_instability.csv`, `uncertainty_neighbor_stability.csv`,
`uncertainty_scvi_vs_totalvi_delta.csv`, `uncertainty_reliability_bins.csv`,
`uncertainty_coverage.csv`, `uncertainty_associations.csv`,
`uncertainty_subset_associations.csv`, `uncertainty_distributions.csv`,
`uncertainty_source_latent_posterior.csv`, `uncertainty_source_vs_target.csv`.

Figures (PNG + PDF): `phase11_figure1_taxonomy` … `phase11_figure13_source_vs_target`.

Logs: `posterior_api_validation.md`, `phase11_methods.md`,
`phase10_reproduction_check.json`, `latent_posterior_summary.json`,
`risk_coverage_curves.json`, gate JSONs, console logs.

Arrays: per-seed probability matrices with cell IDs and explicit class order;
source `U_latent` vectors.

Code: `src/experiments/phase11_uncertainty.py`, `scripts/11_gate.py`,
`scripts/11_run_uncertainty.py`, `scripts/11_subset_latent.py`,
`scripts/11_figures.py`, `scripts/11_augment_delta.py`,
`scripts/10c_rescore_corrected.py`; one-line fix in
`src/experiments/phase10_cross_dataset.py`.

PHASE 10 erratum: `results/phase10_cross_dataset/tables_corrected/*` and
`logs/phase10_erratum_latent_alignment.md`.

## AF. Exact rerun command

```bash
ssh 4060-server
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate multiomics_robustness
export TMPDIR="$HOME/tmp" MPLBACKEND=Agg PATH="/usr/lib/wsl/lib:$PATH"
cd /home/jizhu/research/multiomics_robustness

python scripts/11_gate.py --direction pbmc5k_to_pbmc10k --seed 0 --level l2 --probe-api
python scripts/11_gate.py --direction pbmc10k_to_pbmc5k --seed 0 --level l2
python scripts/10c_rescore_corrected.py
python scripts/11_run_uncertainty.py --stage main
python scripts/11_subset_latent.py --stage subset
python scripts/11_subset_latent.py --stage latent
python scripts/11_figures.py
python scripts/11_augment_delta.py
```

---

## AG. Are totalVI uncertainty estimates informative under dataset shift?

**Yes.** Error-detection AUROC 0.847 (A) and 0.855 (B), AUPRC ~2.7× prevalence,
ECE 0.022–0.040, and rejecting the most uncertain 50% raises accuracy from ~0.80
to ~0.97. They carry real information about which predictions to distrust.

## AH. Are totalVI estimates calibrated, discriminative, both, or neither?

**Both, with a direction-dependent qualification on calibration.** Discrimination
is strong and consistent (AUROC ~0.85 in both directions). Calibration is good in
Direction A (ECE 0.0220, the best of the four cells) and mediocre in Direction B
(ECE 0.0399, worse than scVI's 0.0354). Mild overconfidence in the mid-confidence
range in both. So: **discriminative in both directions; well calibrated in
Direction A, less well calibrated in Direction B.**

## AI. Does totalVI provide more useful uncertainty than scVI_matched?

**No — not reliably.** Its probability vectors are better in aggregate (lower NLL
and Brier across seeds), and Direction A Brier is the one robust win. But on the
pre-specified primary endpoint the models are indistinguishable
(Δ AUROC +0.002 ± 0.023 in A, −0.002 ± 0.028 in B), ECE has no consistent winner,
and totalVI's Direction B ensemble is *worse* calibrated and has markedly less
informative disagreement. Multimodal information improved *representation
quality* in PHASE 10 without translating into better *reliability estimates*.

## AJ. Does uncertainty help explain the PHASE 10 directional asymmetry?

**Only partially, and not in the way one might hope.** Uncertainty quality does
not track the macro-F1 asymmetry: totalVI's large Direction B accuracy advantage
comes with no uncertainty advantage and slightly worse calibration. What
uncertainty *does* explain is the rare-class mechanism behind the asymmetry —
support-vs-uncertainty correlations are stronger in Direction B (ρ = −0.57 for
totalVI), consistent with the smaller PBMC5k source giving thin rare-class support
and driving both the errors and the instability. And §36's question resolves
negatively: totalVI's better neighbor *label agreement* coexists with worse
neighbor *stability* in both directions, so purer neighborhoods are not more
reproducible ones.

## AK. Can full-model uncertainty predict sensitivity to source-subset changes?

**Yes, strongly — the clearest positive finding in PHASE 11.** Spearman ρ = 0.73–0.89
between full-source uncertainty and PHASE 10B instability, and the top-20%
uncertainty stratum shows 372× (scVI) to ∞ (totalVI) the prediction variation
ratio of the bottom-20%. A single trained model already knows which target cells
are fragile to reasonable perturbation of the source cohort, without needing the
five source subsets to be trained at all.

## AL. Is there sufficient evidence for a later uncertainty-aware multimodal method?

**Partially, and not the version originally implied.** Evidence *for* an
uncertainty-aware selective-prediction or abstention layer is strong: rejection
works, uncertainty forecasts source-subset fragility, and it flags exactly the
rare classes that fail.

Evidence *against* uncertainty-aware **modality weighting** specifically: totalVI's
uncertainty is not reliably better than scVI's, so weighting modalities by totalVI
uncertainty is not supported by these data. Evidence against composite reliability
scores: the components behave very differently — neighbor stability is
uninformative about error while `U_pred` is highly informative — so pooling them
would dilute signal. §48 is respected; no composite score was built.

## AM. What should the next phase be?

The highest-value next step is **not** a new model. Ranked:

1. **Fix the artifact policy and settle the latent-posterior question.** Persist
   adapted query models (~1.5 GB for the 20-model grid). This is the only reason
   concept A is missing, and re-running adaptation alone is far cheaper than any
   new method.
2. **Increase model seeds from 5 to ~20.** Between-seed variability dominated six
   of ten primary comparisons. Almost every scVI-vs-totalVI question in PHASE 11
   is currently seed-limited, and no new method should be judged against a
   baseline this noisy.
3. **Cross concepts C and D.** PHASE 10B has one model seed per subsample, so
   model-seed and source-subset variation cannot be separated. A small
   seeds × subsets grid would let both be estimated properly.
4. **Then, if 1–3 hold up, a selective-prediction / abstention method** with a
   pre-registered risk-coverage target — justified by AK — rather than
   uncertainty-weighted modality fusion, which these results do not support.

Applying the same diagnostics to the existing PHASE 8 corruption models (§49)
would be worthwhile and cheap, since it reuses this exact pipeline; it was not
run and awaits approval. MultiVI and RNA+ATAC remain out of scope.
