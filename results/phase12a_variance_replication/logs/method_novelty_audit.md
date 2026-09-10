# PHASE 12A Part 2 — Method novelty audit

Purpose: establish, before any architecture is written, whether
per-cell posterior-precision-weighted multimodal latent fusion is novel.

**Verdict: it is not novel. It is Gaussian Product-of-Experts, which has been
the standard multimodal-VAE fusion rule since 2018, and it is already
implemented in at least three single-cell integration methods, two of them
peer-reviewed. The per-cell character of the weighting is not novel either: it
comes free with PoE. This idea must not be proposed as a contribution.**

Two independent audits were run, one covering the single-cell methods named in
the phase brief and one covering the general multimodal-VAE / product-of-experts
literature. They agree without exception. Source code was inspected where
available rather than relying on paper prose.

---

## 1. The mathematical core

For cell *i*, modality *m*, with a diagonal Gaussian posterior
q_m(z | x_im) = N(mu_im, diag(sigma^2_im)) and precision tau_im = 1/sigma^2_im,
the product of Gaussian experts is itself Gaussian with, per latent dimension d:

    tau_joint_id = SUM_m tau_imd
    var_joint_id = 1 / SUM_m tau_imd
    mu_joint_id  = ( SUM_m mu_imd * tau_imd ) / ( SUM_m tau_imd )

This is stated verbatim in Wu & Goodman, "Multimodal Generative Models for
Scalable Weakly-Supervised Learning", NeurIPS 2018, section 2.1 and Lemma B.1
(arXiv:1802.05335), and independently in Cao & Fleet, "Generalized Product of
Experts...", 2014 (arXiv:1410.7827), equations 2-3.

The implied weight on modality *m* for cell *i* is therefore exactly

    w_im = tau_im / SUM_m' tau_im'

**Precision weighting is not a way of doing Gaussian PoE; it is what Gaussian
PoE mechanically is.** Gaussians are an exponential family whose natural
parameters are the precision and the precision-weighted mean, so multiplying
the experts adds precisions by definition. There is no mathematical content in
"posterior-precision weighted fusion" beyond the definition of Gaussian PoE.

Two corollaries that also cannot be claimed:

- **Cell-specific weighting is free.** In any amortised model,
  (mu_im, sigma^2_im) = Enc_m(x_im), so the precisions are functions of that
  cell's own measurements and the weights already differ per cell. Cao & Fleet
  write the precision as T_i(x), explicitly input-dependent.
- **Per-latent-dimension weighting is free.** Diagonal PoE already weights each
  latent dimension separately. A scalar per-cell w_im is a *restriction* of
  what Cobolt and scMaui already compute, not a generalisation of it.

---

## 2. Method matrix

Machine-readable version: `tables/method_novelty_matrix.csv`.
1 = yes, 0 = no.

| Method | Year | Per-modality posteriors | Fusion rule | Precision-weighted | Cell-specific weights | Calibration evaluated | Domain shift | Selective prediction | Source-composition perturbation |
|---|---|---|---|---|---|---|---|---|---|
| MVAE (PoE) | 2018 | yes | Gaussian PoE | **yes** | yes (free) | no | no | no | no |
| MMVAE | 2019 | yes | MoE, alpha=1/M | no | no | no | no | no | no |
| MoPoE-VAE | 2021 | yes | mixture over subset PoEs | **yes** | yes | no | no | no | no |
| gPoE-VAE (Joshi) | 2022 | yes | PoE, per-sample learned alpha | **yes** | yes | no | no | no | no |
| gPoE-normVAE | 2023 | yes | PoE, per-dim learned alpha | **yes** | no | no | no | no | no |
| COLD Fusion | 2023 | yes | calibrated inverse-variance weights | **yes** | yes | **yes** | no | no | no |
| totalVI | 2021 | **no** (single joint encoder) | input concatenation | n/a | n/a | partial (coverage) | no | no | no |
| MultiVI | 2023 | yes | weighted average of mu and var | no (learned) | **yes, tested** | no (monotonicity only) | partial | no | no |
| Cobolt | 2021 | yes | Gaussian PoE | **yes** | yes | no | structural | no | no |
| scMM | 2021 | yes | MoE, alpha=1/M | no | no | no | no | no | no |
| scMaui | 2024 | yes | Gaussian PoE | **yes** | yes | no | no | no | no |
| UniVI | 2025 (preprint) | yes | precision-weighted PoE + gate | **yes** | yes | no (named future work) | **yes, extensive** | qualitative only | no |
| MOFA+ | 2020 | no (feature-level precisions) | shared factors | no | no | no | no | no | no |

---

## 3. Per-method records

### totalVI — Gayoso et al., Nature Methods 18:272-282 (2021), doi:10.1038/s41592-020-01050-x
RNA + surface protein. **One joint encoder over concatenated RNA and protein**;
there are no modality-specific posteriors, therefore no fusion rule over
posteriors and nothing to precision-weight. Auxiliary latents (protein
background beta_n, RNA size factor l_n) are nuisance variables, not
modality-specific cell-state posteriors. Reports a genuine held-out calibration
error (Methods eq. 31: squared deviation of posterior-predictive interval
coverage from 1) — but over counts, not for a classifier, and not ECE or Brier.
Supports missing modalities at dataset level. No cross-dataset transfer
evaluation of uncertainty, no abstention, no source-perturbation study
(only stability across 5 initialisations).

### MultiVI — Ashuach et al., Nature Methods 20:1222-1231 (2023), doi:10.1038/s41592-023-01909-9
RNA + ATAC (+ protein). Modality-specific encoders with per-modality variances.
Fusion is a **weighted arithmetic average of both means and variances** with
free learnable weights, plus a symmetric-KL or MMD coupling penalty. This is
*not* precision weighting: the weights never read the variances, and the joint
variance is a weighted sum rather than 1/SUM(tau), so it does not contract as
evidence accumulates.

**Direct threat to any cell-specific-weighting claim.** The paper implements and
tests `weighting="cell"`, documented as "Learn weights across modalities and
cells. w_{m,c}", and reports that global weighting, cell-specific weighting and
MMD "all... yielded highly similar results" on DOGMA-seq and TEA-seq. Any future
proposal involving per-cell modality weights must engage this published negative
result directly. Two legitimate distinctions exist: their weights are learned
free parameters whereas precision-derived weights need none and transfer to
unseen cells; and they evaluated only LISI/scib integration quality, never error
detection, calibration or abstention.

### Cobolt — Gong et al., Genome Biology 22:351 (2021), doi:10.1186/s13059-021-02556-z
mRNA + ATAC; framework modality-general. Ships a literal `ProductOfExperts`
class computing `T = 1./var; pd_mu = sum(mu*T)/sum(T); pd_var = 1./sum(T)` —
exactly the precision-weighted mean and reciprocal-summed-precision variance,
per cell and per latent dimension. This is the **canonical peer-reviewed
single-cell instance of the proposed fusion rule**. MultiVI's own Methods
corroborates it as "the classical product of experts calculation used by
Cobolt".

### scMaui — Jeong et al., BMC Bioinformatics 25:257 (2024), doi:10.1186/s12859-024-05880-w
Assay-agnostic, benchmarked on CITE-seq, RNA+ATAC and RNA+methylation.
Self-describes as a variational product-of-experts autoencoder. Code
(`layers.py`, `JointMean`/`JointSigma`) computes
`jointmean = SUM(mask * mu * exp(-sigma)) * reciprocal(SUM(mask * exp(-sigma)))`
where `sigma` is the log-variance, so `exp(-sigma)` is the precision.
Per-cell, per-dimension precision weighting with masked missing modalities.
No calibration, domain-shift, abstention or source-perturbation evaluation.

### UniVI — bioRxiv 2025.02.28.640429 (preprint, not peer-reviewed)
**Closest prior art, and the most damaging.** Three separate overlaps:
1. `mixture_of_experts()` computes `precisions = exp(-logvar)`,
   `mu_comb = SUM(mu*prec)/SUM(prec)`, `var_comb = 1/SUM(prec)` — precision-weighted PoE.
2. `evaluation.py::encode_moe_gates_from_tensors` exposes
   `kind="effective_precision"` computing `s_m = SUM_d exp(-logvar_m[d])`,
   `w_m = s_m / SUM_j s_j` — **the proposed scalar w_im verbatim**.
3. It uses that weight for **the same purpose we would**: a per-cell "support
   map" indicating "which regions should be interpreted cautiously when one
   modality is weak or absent".

Also the most thorough domain-shift evaluation in the set (parameter-frozen
projection of unimodal cohorts, mosaic bridging across independent AML cohorts,
inductive vs transductive framing, overlap sweeps). Explicitly names "more
explicit uncertainty calibration" as future work — a citable admission of the
gap that remains.

### scMM — Minoura et al., Cell Reports Methods 1:100071 (2021), doi:10.1016/j.crmeth.2021.100071
CITE-seq and SHARE-seq. Mixture-of-experts with **fixed uniform alpha_m = 1/M**,
independent of posterior variance. The clean negative control: modality-specific
posteriors exist but are not precision-weighted, and a mixture does not contract
variance the way a product does.

### MOFA+ — Argelaguet et al., Genome Biology 21:111 (2020), doi:10.1186/s13059-020-02015-1
Linear Bayesian multi-view factor analysis; not an autoencoder. Has view- and
feature-specific noise precisions tau^m_d, which do inverse-variance-weight each
view's residual contribution — but this is **per feature, per view, per group,
never per cell**. No per-cell modality weight exists.

### General multimodal-VAE prior art
- **MVAE**, Wu & Goodman, NeurIPS 2018 (arXiv:1802.05335) — origin of the
  Gaussian PoE VAE, including a prior expert so fused precision is 1 + SUM tau.
- **MMVAE**, Shi et al., NeurIPS 2019 (arXiv:1911.03393) — MoE alternative;
  source of the standard criticism that PoE precisions are miscalibrated.
- **MoPoE-VAE**, Sutter et al., ICLR 2021 (arXiv:2105.02470) — unifies PoE and
  MoE as abstract means; both MVAE and MMVAE are special cases.
- **gPoE**, Cao & Fleet 2014 (arXiv:1410.7827) — rescaled precisions
  SUM alpha_m tau_m, motivated explicitly by expert overconfidence.
- **gPoE-VAE**, Joshi et al., ICMI 2022 (doi:10.1145/3536221.3556596) —
  per-sample, per-latent-dimension learned alpha from auxiliary networks.
- **gPoE-normVAE**, Lawry Aguila et al., MICCAI 2023
  (doi:10.1007/978-3-031-43907-0_41) — per-dimension learned alpha, motivated
  verbatim by "an overconfident miscalibrated expert".
- **COLD Fusion**, Tellamekala et al., IEEE TPAMI 2023
  (doi:10.1109/TPAMI.2023.3325770) — learns *calibrated and ordinal* latent
  variances via a distributional matching loss, then fuses with normalised
  inverse-variance-norm weights. Supervised and discriminative, not a VAE.
- **PoE non-identifiability**: Kurle et al., AAAI 2019
  (doi:10.1609/aaai.v33i01.33014114) — "the natural parameters of the individual
  normal distributions are not uniquely identifiable by the natural parameters
  of the integrated normal distribution". Wu & Goodman say the same.

---

## 4. Answer to the critical novelty question (brief §9)

**Does any existing method already use w_im proportional to posterior precision,
or equivalent Gaussian precision weighting? YES — verified in source code in
three single-cell methods and originating in the general literature in 2018.**

- Cobolt (Genome Biology 2021) — `ProductOfExperts`, peer-reviewed.
- scMaui (BMC Bioinformatics 2024) — `JointMean`/`JointSigma`, peer-reviewed.
- UniVI (bioRxiv 2025) — precision-weighted PoE **plus** an explicit
  `effective_precision` diagnostic computing w_im exactly, used as a per-cell
  reliability map.
- Originating form: Wu & Goodman, NeurIPS 2018.

The method mathematically closest to the proposed idea is **UniVI**, which
matches on the fusion rule, on the exact scalar weight, and on the motivation.
For peer-reviewed prior art, **Cobolt** is the canonical citation.

This must therefore never be described as novel. If precision weighting is used
at all, the correct framing is "we adopt the standard Gaussian PoE fusion of
Wu & Goodman (2018), as used in single-cell by Cobolt and scMaui", cited as
prior work in the methods section.

---

## 5. A correction that matters for our own framing

**totalVI has no modality-specific encoders.** RNA and protein counts are
jointly transformed by a single encoder. Our PHASE 11B finding — that totalVI's
latent posterior variance detects downstream error (AUROC 0.61-0.67) while
scVI's does not (0.47-0.54) — is therefore a statement about a **single joint
posterior**, not about a fusible per-modality posterior.

Consequences:
- We must not write "totalVI's protein-modality posterior variance".
- There is nothing inside totalVI to precision-weight, so our own empirical
  result does not, by itself, motivate a fusion rule at all. It motivates
  treating posterior width as a *reliability signal*.
- A per-modality comparison in the CITE-seq setting would require MultiVI's
  protein extension, scMM, or UniVI — not totalVI.

---

## 6. Where the real gap is (brief §10)

Assessed against the six candidate gaps in the phase brief.

**A. Calibration-aware multimodal fusion — largely occupied.** COLD Fusion
(TPAMI 2023) learns calibrated variances and fuses with them; the gPoE-VAE
family rescales precisions to fix overconfidence. Cao & Fleet stated the
motivation in 2014. Remaining space is narrow and mechanism-level only:
a *label-free* calibration diagnostic and correction applied to posterior
precisions inside a self-supervised model. Neither auditor found that exact
combination, but both flagged the surrounding space as dense and rated the
chance that something close exists as moderate-to-high. Not a safe primary claim.

**B. Domain-shift-aware reliability estimation — mostly open.** UniVI evaluates
domain shift thoroughly but does not evaluate uncertainty quality under it, and
names calibration as future work. No method measures whether latent posterior
variance predicts cross-dataset error.

**C. Uncertainty-guided selective prediction with class-retention control —
fully open.** Zero of the thirteen methods report risk-coverage curves,
accuracy-versus-retention, per-class retention, or any formal deferral rule.
UniVI's "support map" is qualitative advice, not a quantified evaluation. This
is directly motivated by our own PHASE 11B failure finding: abstention cuts risk
6-7x but retention CV reaches 0.88-1.02 at 50% coverage, i.e. it works by
discarding rare cell types.

**D. Source-perturbation fragility prediction — fully open.** No method
perturbs the *composition of the source training set*. The existing perturbation
studies (MultiVI's unpairing sweep, UniVI's overlap sweep and per-population
dropout, scMaui's imputation sweep) all perturb *input modality availability*,
which is a different thing. PHASE 11B found full-model predictive uncertainty
predicts crossed-design instability at rho ~= 0.86.

**E. Fusion optimising representation and calibration jointly — partly
occupied** by COLD Fusion, and weakened for us by MultiVI's published negative
result that per-cell modality weights did not change integration quality.

**F. Other — the strongest option: reconciling MultiVI's negative result.**
MultiVI found w_{m,c} made no difference *on LISI/scib integration metrics*.
Integration metrics are blind to error detection, calibration and abstention.
The defensible claim is that modality weighting may be neutral for embedding
quality while being informative for reliability — a testable reframing that
turns a published null into a motivation.

**Recommendation: C and D, framed as B, with F as the argument.** The
contribution should be evaluative and diagnostic, not architectural: a
reliability and deferral framework for cross-dataset single-cell transfer, with
per-class retention as a first-class outcome. This is where all thirteen audited
methods are silent.

*Caveat: the relative weight of C versus D depends on the pending PHASE 12A
variance result. If the replicated decomposition shows the source-subset
component is negligible once pure run-level noise is separated out, then "source
fragility" is the wrong target and the finding becomes "run-level irreproducibility",
which shifts the emphasis onto C and onto reporting practice. This section will
be revisited once the 100-run design is scored.*

---

## 7. Requirements any future method must satisfy (brief §11)

1. **Must not duplicate precision weighting.** Ruled out by sections 1 and 4.
2. **Must differ mathematically from totalVI / MultiVI / UniVI / scMaui /
   Cobolt.** Note that totalVI has no per-modality posterior at all, and that
   Cobolt/scMaui already do per-dimension precision weighting, which is strictly
   more general than a per-cell scalar.
3. **Must be motivated by an empirical failure found in PHASE 8-11B.** The two
   qualifying failures are: latent posterior variance is uninformative for scVI
   but informative for totalVI (PHASE 11B), and abstention disproportionately
   removes rare classes (PHASE 11B, retention CV 0.88-1.02 at 50% coverage).
4. **Must have a falsifying baseline.** Non-negotiable baselines: plain PoE
   (MVAE/Cobolt), MoE (MMVAE/scMM), MoPoE, gPoE with learned alpha in both the
   per-sample (Joshi) and per-dimension (Lawry Aguila) variants, MultiVI with
   `weighting="cell"`, and a 10-seed ensemble. Omitting the gPoE baselines is
   the single most likely cause of desk rejection on novelty grounds.
5. **Must report calibration separately from discrimination.** These dissociate
   throughout PHASE 11B.
6. **Must report rare-class and per-class consequences.**

---

## 8. Required citations

Cobolt (Gong et al. 2021) and scMaui (Jeong et al. 2024) for peer-reviewed
single-cell PoE precision weighting; UniVI (2025) as closest prior art including
its explicit `effective_precision` weight; MultiVI (Ashuach et al. 2023) for the
`w_{m,c}` negative result; Wu & Goodman (NeurIPS 2018) for the origin of the
PoE-VAE form; Cao & Fleet (2014) for generalised precision rescaling; Shi et al.
(2019) and Kurle et al. (2019) for PoE precision miscalibration and
non-identifiability; totalVI eq. 31 for the one real calibration metric in the
single-cell set; COLD Fusion (TPAMI 2023) for calibrated variance-based fusion
weights.

Adjacent work a reviewer will know and that should be cited proactively:
uncertainty quantification for atlas-level cell-type transfer
(arXiv:2211.03793, reports ECE and calibration curves under dataset shift) and
conformal prediction for single-cell annotation with rejection
(doi:10.1093/bioinformatics/btaf521; Khatri & Bonn, PMLR v179, 2022).

## 9. Stated limitations of this audit

- UniVI is an unrefereed preprint; its claims rest on its public repository.
- Cobolt's own supplementary derivation was not read; the PoE was confirmed from
  the main-text MVAE citation, the source code, and MultiVI's independent
  description.
- MultiVI's Supplementary Fig. 3a,b was not retrieved, so the numerical size of
  the `w_{m,c}` versus equal-weight difference is known only from the authors'
  summary that results were "highly similar".
- The claim that no prior work applies a label-free calibration correction to
  posterior precisions before fusion inside a self-supervised VAE is **not
  verified as absent, only as not found**. The adjacent literature is dense.
