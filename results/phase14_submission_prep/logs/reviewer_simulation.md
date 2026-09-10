# Simulated peer review (no manuscript changes from these opinions)

## Reviewer 1 — single-cell multi-omics methods

**Summary.** A careful reliability study of scVI/totalVI for CITE-seq with corruption, transfer, uncertainty, and independent donor-held-out validation. Not an architecture paper.

**Strengths.** Clear claim discipline; post-adaptation coordinate hygiene; external Lawlor validation with author labels; separation of uncertainty facets; abstention coverage trade-off.

**Major concerns.**
1. PBMC-only biology may limit general interest.
2. Why not include MultiVI/other multimodal VAEs as comparators?
3. Entry-wise corruption ≠ missing modality — is the robustness claim overstated?

**Minor.** Figure composites need human polish; some panels denser than ideal.

**Challenged claims.** “Reliability” as general property; latent-uncertainty superiority outside tested settings.

**Requested experiments.** Non-PBMC cohort; missing-modality experiment; broader model panel.

**Necessary?** Non-PBMC: high-value optional. Missing-modality: not honestly supported by totalVI API used → rewrite already addresses. Broader models: optional supplement, not required for present claims.

**Likely recommendation.** Minor revision / accept after tightening scope language.

## Reviewer 2 — statistics / uncertainty

**Summary.** Distinguishes discrimination, scoring rules, calibration, and selective risk; reports seed/fold variability; replicated residual decomposition.

**Strengths.** 20-seed transfer; paired sign consistency; residual vs interaction separation; external attenuation reported honestly.

**Major concerns.**
1. Limited factor levels for variance ANOVA/power.
2. Downstream logistic regression uncertainty ≠ model posterior predictive.
3. Possible multiple-comparisons narrative across many metrics.

**Minor.** ECE binning sensitivity mostly in supplement.

**Challenged claims.** Causal reading of multimodal training → informative latent variance.

**Requested experiments.** End-to-end classifiers; more replicates; formal mixed models.

**Necessary?** Mostly B/C (rewrite + existing supplement). New end-to-end model: optional, not required for latent/predictive distinction already drawn.

**Likely recommendation.** Minor revision.

## Reviewer 3 — computational biology journal editor/reviewer

**Summary.** Fits Genome Biology if framed as reliability characterization with biological coverage consequences; desk risk if it reads as endless benchmarking.

**Strengths.** External validation capstone; concrete Δ=+0.088 Lawlor result; abstention trade-off is biologically interpretable.

**Major concerns.**
1. Contribution vs “just benchmarking.”
2. Stimulation conditions unused in Lawlor.
3. Software/resource packaging incomplete for some readers.

**Minor.** Title/abstract must lead with reliability payoff.

**Requested experiments.** LPS/CD3_CD28; 39-ADT; packaged tutorial.

**Necessary?** Stimulation: high-value optional, not required. 39-ADT: low-value optional for primary claim. Tutorial: B/C packaging, not science gap.

**Likely recommendation.** Send to review; minor/major revision depending on framing polish.

## Consensus predicted concerns (top 5)
1. PBMC-only / limited biological breadth
2. Missing-modality confusion with entry-wise corruption
3. Incomplete multimodal method panel
4. Downstream classifier vs generative uncertainty
5. Optional stimulation shift not shown
