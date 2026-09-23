# FINAL HUMAN-STYLE MANUSCRIPT AUDIT

**Canonical manuscript:** `results/final_information_vs_model_manuscript/manuscript/CITEseq_Information_vs_Model_Reliability_Research_Manuscript.docx`  
**Pre-audit backup:** `.../CITEseq_Information_vs_Model_Reliability_Research_Manuscript_pre_audit.docx`  
**Markdown mirror:** updated in place (pre-audit backup `.md` saved)  
**Date:** 2026-09-23  
**Mode:** Writing/interpretation audit only. No new experiments, datasets, training, commit, or push.

---

## 3-minute application reader test

**Question:** What is this paper about?

**Answer a careful reader should give:**  
It separates the benefit of having protein measurements from the additional benefit of a complex multimodal model, and shows that the information benefit is strong while the residual model-specific benefit is smaller and context-dependent.

**Verdict after Title / Abstract / Figures 1–3 / Discussion ¶1 / Conclusion:** **PASS** (framing no longer reads as a generic scVI–totalVI robustness benchmark).

---

## Task 1 — Old-narrative search (selected findings)

| Location | Sentence / phrase | Class | Action |
|---|---|---|---|
| Abstract (pre) | “Additional totalVI-specific gains existed… but were not universal” | SOFTEN | Rewrote Conclusions to “Protein-information gain was strong; residual model-specific gain was smaller and context-dependent.” |
| Results §1 | “totalVI remained the highest internal representation…” | SOFTEN | Rewrote to factual ranking + emphasize \(G_{\mathrm{info}}\) vs \(G_{\mathrm{model}}\) sizes |
| Results §2 (pre) | Numbers present but pattern easy to miss | REWRITE | Explicit inductive / transductive / VAE blocks + Case C bullets |
| Results §3 heading | “reduced multimodal model advantage” | SOFTEN | “narrowed…”; gradual-not-catastrophic language |
| Results §4 (pre) | Predictive/latent terminology loose | REWRITE | Classifier predictive uncertainty; latent posterior variance; anti-“uncertainty improved” |
| Results §7 (pre) | “Uncertainty… also favored totalVI” | REWRITE | Diagnostics under protein-gated labels; required Lawlor wording |
| Discussion ¶1 (pre) | Good but not unmistakable | REWRITE | Inserted recommended “most reproducible multimodal benefit…” sentence |
| Methods \(G_{\mathrm{model}}\) | Descriptive but missing exact non-causal sentence | REWRITE | Added required Methods + Discussion causal-decomposition disclaimer |
| Captions / SI | Negations of “universal superiority” | KEEP | Correct cautionary language |
| “epistemic / robust to corruption / uncertainty improved” | Appear only in **negation** | KEEP | Explicitly rejected |

**No KEEP→DELETE of valid secondary experiments.**

---

## Tasks 2–13 — Section audits (summary)

| Check | Status |
|---|---|
| Abstract information-vs-model + 70% descriptive + Case C + protein-gated Lawlor | PASS (~218 words) |
| Results logic chain R1→R7 answers central question | PASS after reconnects |
| Transfer: inductive + transductive both reported; Case C explicit | PASS |
| Figure 3 Case C annotations (+0.016 / −0.113); y-axis from 0 | PASS (regenerated) |
| Figure 2 RNA→concat→totalVI hierarchy | PASS |
| \(G_{\mathrm{model}}\) non-causal in Methods + Discussion | PASS |
| Degradation framed as narrowing advantage | PASS |
| Uncertainty subordinated; terminology corrected | PASS |
| Stochasticity linked to small-effect interpretation | PASS |
| Lawlor protein-gated + protein/concat > totalVI prominent in Results | PASS |
| Discussion dominant interpretation | PASS |
| Limitations checklist (PBMC, panels, classification focus, protein-gated, scArches≠PCA, no protein-independent external, no cell-wise missingness, no universal architecture comparison) | PASS |

---

## Task 14–15 — Matched transfer fairness / code

See `review/MATCHED_TRANSFER_CODE_AUDIT.md`.

- Target labels: **no leakage** into HVG/PCA/scaling/classifier fit.
- Eligible shared classes use target labels only for **evaluation-set definition**.
- Case C confirmed: A ≈ +0.016; B ≈ −0.113.
- Target-data matching **does** change central interpretation (Case C).

---

## Task 17 — Figure audit

| Figure | Verdict |
|---|---|
| 1 | Central logic clear in ~10s |
| 2 | Primary trio dominates; secondaries de-emphasized |
| 3 | Case C immediately readable after annotation + full y-axis |
| 4 | Progressive Direction B narrowing visible |
| 5 | Secondary diagnostics; acceptable density |
| 6 | 10/10 donor Δ **and** protein/concat > totalVI; protein-gated annotated |

---

## Task 20 — Three reviewer simulations

### A. Computational single-cell reviewer

**CRITICAL:** None remaining after Case C + Lawlor label–modality wording.  
**MAJOR:** (1) Concat vs scArches still not algorithmically matched — now stated clearly. (2) PBMC / panel scope.  
**MINOR:** Future protein-independent-label cohort.  
**STRENGTHS:** Clear information-vs-model decomposition; honest simple baselines; matched target-feature control.

### B. Statistical-methods reviewer

**CRITICAL:** None. \(G_{\mathrm{info}}/G_{\mathrm{model}}\) and ~70% explicitly descriptive/non-causal.  
**MAJOR:** Internal bootstrap remains cell-level conditional (stated). Donor n=10 primary for Lawlor.  
**MINOR:** Notation mixed between unicode and LaTeX in MD→DOCX conversion.  
**STRENGTHS:** Case C direction dependence reported without soft-pedaling; variance residual correctly named.

### C. Graduate admissions faculty reader

**CRITICAL:** None for application use.  
**MAJOR:** Document remains unpublished research manuscript (appropriate).  
**MINOR:** Word layout is professional single-column A4 (prior package convention; two-column journal typesetting not forced).  
**STRENGTHS:** 3-minute takeaway matches intended scientific identity; intellectual honesty about Lawlor.

---

## Fixes applied (material rewrites)

1. Abstract Conclusions + Case C / degradation / Lawlor emphasis.  
2. Transfer Results: explicit A/B number pattern; non-uniform benefit of unlabeled target access.  
3. Degradation Results: narrowed-advantage framing; Dir B vs clean scVI.  
4. Uncertainty Results: classifier predictive vs latent posterior variance; anti-“improved.”  
5. Stochasticity: linked to \(G_{\mathrm{model}}\) interpretation.  
6. Lawlor Results: required protein-marker / label–modality alignment paragraph.  
7. Discussion ¶1 + non-causal model-specific language; expanded limitations.  
8. Methods: equations for \(G_{\mathrm{info}}/G_{\mathrm{model}}\); fairness audit paragraph; uncertainty notation.  
9. Figure 3 regenerated with Case C callouts and non-truncatedating y-axis.

---

## Final report answers (1–22)

1. **Old “totalVI > scVI” narrative remaining?** No as a general claim. Direction B / Lawlor contrasts remain as scoped evidence.  
2. **Exact sentences materially rewritten:** Abstract Conclusions; Result 1 secondary ranking; entire Result 2 transfer block; Result 3 heading+interpretation; Result 4 uncertainty block; Result 5 stochasticity open; Result 7 Lawlor block; Discussion ¶1–limitations; Methods G/transfer/uncertainty; Fig 3–4 legends.  
3. **Abstract distinguishes information vs model?** YES.  
4. **~70% marked descriptive?** YES (“without causal attribution”).  
5. **Model-specific gain non-causal/descriptive?** YES (Methods + Discussion required sentence).  
6. **Transfer Case C clearly stated?** YES.  
7. **Source-only and transductive concat both reported?** YES.  
8. **Figure 3 communicates Case C?** YES (after revision).  
9. **Target-label leakage?** NONE into fitting; eval-set definition only.  
10. **Target-data matching changes interpretation?** YES (Case C).  
11. **Degradation framed as narrowing model advantage?** YES.  
12. **Uncertainty subordinated?** YES.  
13. **Stochasticity linked to small effects?** YES.  
14. **Lawlor protein-gated prominent?** YES (Results, not only Limitations).  
15. **Protein PCA / concat > totalVI visible?** YES (text + Figure 6).  
16. **10 donors as biological units?** YES.  
17. **Universal architecture-superiority claims?** NONE remaining.  
18. **Remaining CRITICAL?** None.  
19. **Remaining MAJOR?** Scope limits (PBMC-only; protein-gated Lawlor; concat≠scArches algorithmically) — already stated, not writing defects.  
20. **New experiment required?** NO.  
21. **Ready as graduate-application research manuscript?** YES, pending human review of the audited DOCX.  
22. **Final DOCX path:** `results/final_information_vs_model_manuscript/manuscript/CITEseq_Information_vs_Model_Reliability_Research_Manuscript.docx`

---

## STOP

No experiments · no new external cohort · no model training · no commit · no push · no journal submission.  
Awaiting human review.
