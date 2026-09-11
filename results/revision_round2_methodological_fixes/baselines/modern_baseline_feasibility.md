# Modern multimodal baseline feasibility (PART 10)

## Installed stack (project pin)

- `scvi-tools==1.3.3` (see `requirements.txt` / `environment.yml`)
- Python 3.11 target env on server; local analysis used lightweight venv

## MultiVI

| Question | Verdict |
|---|---|
| 1. Native RNA + protein? | **No / not fair.** MultiVI in scvi-tools is designed for **RNA + ATAC** (chromatin). Protein (ADT) is not a first-class MultiVI modality in 1.3.3. |
| 2. Same cells/features? | Would require forcing ADT into an ATAC-like slot — **not honest**. |
| 3. Missing protein semantics? | N/A if modality misuse. |
| 4. Latent dim match? | Possible but irrelevant if modality unfair. |
| 5. Same train/test design? | Possible mechanically, not scientifically fair. |
| 6. Evaluation without label leakage? | Same classifier discipline possible. |

**Decision: DO NOT force MultiVI into the RNA+protein benchmark.**

## Cobolt

Cobolt is primarily **RNA + ATAC** multiomic integration. Not a native CITE-seq
RNA+protein generative baseline for this comparison.

## Preferable ONE modern RNA+protein candidate (if any)

Candidates with native CITE-seq / RNA+protein support in the literature:

- **scArches-compatible totalVI** (already the multimodal deep baseline)
- **MOFA+** (already included; linear multi-view factor model)
- Possible alternatives: **scMaui**, **scAI**, or newer CITE-seq VAEs — **not currently installed**

## Verdict

**STOP before training an additional deep multimodal model** until a method with
**native RNA+protein likelihood support** is installed and can share:

- same cells/features
- same train/source vs target split
- honest missing-protein semantics if claimed
- comparable latent dimensionality
- leakage-free evaluation

**Additional modern model actually run in this revision round so far: NO.**

Recommended next step (P2, only if needed): install/evaluate one native
RNA+protein method OR elevate the **simple RNA+protein concatenated PCA**
control (PART 5B) as the architecture-vs-information contrast, which does not
require a second deep generative model.
