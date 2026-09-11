# Test-time corruption — HARD STOP / BLOCKED

## Reason

PHASE 11B `adapted_query_model` directories are **not available in this workspace**.

Availability check:
```json
{
  "n_found": 0,
  "n_expected": 80,
  "n_missing": 80,
  "missing_head": [
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed00/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed01/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed02/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed03/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed04/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed05/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed06/adapted_query_model",
    "/Users/venus/Desktop/2026\u5317\u4eac\u5927\u5b66\u79d1\u5b66\u7814\u7a76/Deep Generative Modeling for Single-Cell Multi-Omics Integration/multiomics_robustness/results/phase11b_uncertainty_robustness/models/direction_A/scvi_matched/seed07/adapted_query_model"
  ]
}
```

Expected path pattern:
`results/phase11b_uncertainty_robustness/models/{direction}/{model}/seedXX/adapted_query_model/`

These weights live on the RTX 4060 server (`/home/jizhu/research/multiomics_robustness`)
and are gitignored locally.

## SSH status from this agent environment

`ssh 4060-server` to `192.168.1.162` failed (host unreachable / operation not permitted).

## What is already prepared

- Script: `scripts/revision2_testtime_corruption.py`
- Design: clean source → scArches on CLEAN target → freeze → corrupt TARGET protein only → re-encode → frozen source classifier
- Corruption fractions: [0.0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85]
- Mask seeds: [0, 1, 2, 3, 4]
- Model seeds: 0–19, Directions A/B

## To complete PART C

On the 4060 server:

```bash
ssh 4060-server
cd /home/jizhu/research/multiomics_robustness
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate multiomics_robustness
export TMPDIR="$HOME/tmp" MPLBACKEND=Agg
export PATH="/usr/lib/wsl/lib:$PATH"
tmux new -s revision2_p1
python results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py --smoke
python results/revision_round2_methodological_fixes/scripts/revision2_testtime_corruption.py --run
```

Do NOT fake results. Do NOT retrain. Do NOT skip freeze checks.
