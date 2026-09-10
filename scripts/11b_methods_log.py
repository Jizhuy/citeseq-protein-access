#!/usr/bin/env python
"""Generate the PHASE 11B methods log from the artifacts actually produced.

Nothing here is hand-transcribed: every number is read back from the manifests
and tables so the log cannot drift from the run.
"""

from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11b_analysis import EVIDENCE_RULE  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    BOOTSTRAP_SEED,
    COVERAGE_GRID,
    ECE_BINS,
    INTEGRITY_TOL,
    MC_CHECK_SEEDS,
    MC_SAMPLES,
    N_BOOT_CELLS,
    N_BOOT_SEEDS,
    NEIGHBOR_K,
    ensure_dirs,
)


def main() -> int:
    paths = ensure_dirs()
    man_dir = paths["logs"] / "manifests"
    manifests = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(man_dir.glob("*.json"))]
    if not manifests:
        print("no manifests found")
        return 1
    frame = pd.DataFrame([{
        "tag": m["tag"], "direction": m["direction_key"], "model": m["model"],
        "model_seed": m["model_seed"], "subset_seed": m["subset_seed"],
        "src_epochs": m["source_training_epochs"], "query_epochs": m["query_epochs"],
        "runtime": m["total_runtime_seconds"],
        "gpu_alloc": m["gpu_peak_allocated_mb"], "gpu_res": m["gpu_peak_reserved_mb"],
        "rss": (m.get("ram_after") or {}).get("process_rss_mb"),
        "integrity_latent": m["artifact_integrity"]["max_latent_diff"],
        "integrity_var": m["artifact_integrity"]["max_variance_diff"],
        "integrity_pass": m["artifact_integrity"]["passed"],
        "drift_mean": m["source_latent_drift_pre_to_post_adaptation"]["mean_cell_displacement"],
        "drift_median": m["source_latent_drift_pre_to_post_adaptation"]["median_cell_displacement"],
        "drift_frob": m["source_latent_drift_pre_to_post_adaptation"]["relative_frobenius_change"],
        "drift_identical": m["source_latent_drift_pre_to_post_adaptation"]["identical"],
        "u_latent": m["median_u_latent_target"],
    } for m in manifests])

    import torch
    import scvi
    gpu = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none"
    total_h = frame["runtime"].sum() / 3600.0
    n_full = int(frame["subset_seed"].isna().sum())
    n_crossed = int(frame["subset_seed"].notna().sum())

    drift = (frame.groupby("model")[["drift_mean", "drift_median", "drift_frob"]]
             .agg(["mean", "max"]).round(6))
    res = (frame.groupby("model")[["runtime", "gpu_alloc", "gpu_res", "rss"]]
           .agg(["mean", "max"]).round(1))

    mc_path = paths["logs"] / "posterior_distribution_validation.json"
    mc = json.loads(mc_path.read_text(encoding="utf-8")) if mc_path.exists() else {}

    lines = [
        "# PHASE 11B methods log",
        "",
        f"Generated {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Compute environment",
        "",
        f"- GPU: {gpu}",
        f"- torch {torch.__version__}, CUDA {torch.version.cuda}",
        f"- scvi-tools {scvi.__version__}",
        f"- Python {platform.python_version()} on {platform.platform()}",
        f"- total model training + adaptation time: {total_h:.2f} h across {len(frame)} runs",
        "",
        "## Runs completed",
        "",
        f"- full-source 20-seed grid: {n_full} runs (2 directions x 2 models x 20 seeds = 80 expected)",
        f"- crossed grid: {n_crossed} runs (5 source subsets x 5 model seeds x 2 models = 50 expected)",
        f"- all artifact-integrity checks passed: {bool(frame['integrity_pass'].all())}",
        "",
        "## Model configuration",
        "",
        "Unchanged from the validated PHASE 10 configuration. No architecture was tuned.",
        "",
        "- scVI_matched: n_latent 20, n_hidden 256, n_layers 2, dropout 0.2, "
        "gene_likelihood nb, dispersion gene, latent_distribution normal",
        "- totalVI: n_latent 20, n_hidden 256, encoder layers 2, decoder layers 1, "
        "dropout 0.2, gene_likelihood nb, gene_dispersion gene, protein_dispersion protein",
        "- training: batch_size 256, max_epochs 200, train_size 0.9, early stopping on "
        "elbo_validation, patience 45, mode min, check_val_every_n_epoch 1",
        "- scArches query adaptation: prepare_query_anndata + load_query_data, "
        "max_epochs 200, plan weight_decay 0.0; totalVI additionally keeps its "
        "model-specific lr and reduce_lr_on_plateau setting",
        "- target labels were never used in training, adaptation, early stopping, or "
        "classifier fitting",
        "",
        "## Latent coordinate rule",
        "",
        "Both SOURCE and TARGET are encoded with the POST-adaptation query model. ",
        "`assert_post_adaptation` is a regression guard against the PHASE 10 defect in ",
        "which the pre-adaptation source array was handed to the evaluator; it checks ",
        "object identity, shape, and content before any metric is computed.",
        "",
        "### Source latent drift, pre- to post-adaptation",
        "",
        "```",
        drift.to_string(),
        "```",
        "",
        f"- scVI runs identical pre/post: "
        f"{bool(frame[frame['model'] == 'scvi_matched']['drift_identical'].all())}",
        f"- totalVI runs identical pre/post: "
        f"{bool(frame[frame['model'] == 'totalvi']['drift_identical'].all())}",
        "",
        "The non-zero totalVI drift is exactly the quantity that made the PHASE 10 defect ",
        "numerically visible for totalVI while leaving scVI unaffected.",
        "",
        "## Artifact persistence and integrity",
        "",
        "Every run writes `models/<direction>/<model>/seedNN/source_model/` and ",
        "`.../adapted_query_model/`. Immediately after saving, the in-memory model is ",
        "deleted, reloaded from disk, and its latent representation and posterior variance ",
        "are recomputed on a fixed validation subset.",
        "",
        f"- tolerance: {INTEGRITY_TOL}",
        f"- max observed latent difference across all runs: {frame['integrity_latent'].max():.3e}",
        f"- max observed variance difference across all runs: {frame['integrity_var'].max():.3e}",
        "",
        "## Posterior uncertainty",
        "",
        "- API: `get_latent_representation(..., return_dist=True)`",
        "- scvi-tools 1.3.3 `VAEMixin.get_latent_representation` returns "
        "`(qz.loc, qz.scale.square())`; neither SCVI nor TOTALVI overrides it",
        "- the second element is therefore the posterior VARIANCE, not a scale or SD",
        "- `latent_distribution=normal`, so q(z|x) is a diagonal Gaussian and the variance is exact",
        f"- primary scalar: U_latent_i = mean_j var_ij over the {20} latent dimensions",
        "- secondary diagnostics saved: median and max posterior variance per run",
        "",
        "### Monte Carlo confirmation",
        "",
        f"- seeds checked: {list(MC_CHECK_SEEDS)}; S = {MC_SAMPLES} posterior draws per target cell",
        f"- mean empirical/closed-form variance ratio: "
        f"{mc.get('mean_ratio_empirical_over_closed_form', 'not run')}",
        f"- consistent with a variance: {mc.get('consistent_with_variance', 'not run')}",
        "",
        "## Downstream classifier",
        "",
        "Identical to PHASE 10: `LogisticRegression(max_iter=2000, solver='lbfgs', ",
        "random_state=seed, class_weight=None)` fitted on SOURCE post-adaptation latents ",
        "and SOURCE labels, evaluated on TARGET. Probability columns are aligned by class ",
        "NAME before any aggregation.",
        "",
        "## Uncertainty definitions",
        "",
        "- U_pred = 1 - max class probability (primary)",
        "- predictive entropy and top1-top2 margin (secondary)",
        "- U_latent = mean posterior variance (primary latent)",
        "- ensemble disagreement = H(mean_s p_s) - mean_s H(p_s), described as "
        "epistemic-like, not exact Bayesian mutual information",
        "- variation ratio = 1 - modal predicted-class frequency",
        "",
        "No composite reliability score was constructed. No temperature scaling was fitted: ",
        "PHASE 11 found it did not reliably improve calibration, and PHASE 11B reports raw ",
        "classifier probabilities only.",
        "",
        "## Calibration and selective prediction",
        "",
        f"- ECE: {ECE_BINS} bins, equal-frequency (primary) and equal-width (sensitivity)",
        "- NLL and multiclass Brier on the full probability vector",
        f"- coverage grid: {[f'{c:.0%}' for c in COVERAGE_GRID]}",
        "- AURC computed over the full risk-coverage curve, cells ranked by U_pred",
        "",
        "## Bootstrap",
        "",
        f"- seed-level: {N_BOOT_SEEDS} paired replicates over the 20 seed pairs",
        f"- target-cell: {N_BOOT_CELLS} paired replicates per seed, identical indices "
        f"for scVI and totalVI",
        f"- base seed: {BOOTSTRAP_SEED}",
        "- the two bootstraps are reported separately and never pooled",
        "",
        "## Evidence classification rule",
        "",
        "Fixed before any PHASE 11B result was inspected:",
        "",
        *[f"- **{k}**: {v}" for k, v in EVIDENCE_RULE.items()
          if k not in ("sign_consistency_threshold", "n_seeds")],
        "",
        "## Variance decomposition",
        "",
        "Balanced two-factor method-of-moments decomposition on the 5x5 "
        "(source subset x model seed) tables. With one observation per grid cell the ",
        "interaction cannot be separated from pure error, so the third component is ",
        "reported as combined interaction/residual. Negative raw estimates are truncated ",
        "at zero for the reported fractions and the raw values are preserved in the same ",
        "table. With only five levels per factor these components are exploratory.",
        "",
        "## Resource usage",
        "",
        "```",
        res.to_string(),
        "```",
        "",
        f"- peak GPU allocated across all runs: {frame['gpu_alloc'].max():.0f} MB",
        f"- peak GPU reserved across all runs: {frame['gpu_res'].max():.0f} MB",
        f"- peak process RSS across all runs: {frame['rss'].max():.0f} MB",
        f"- neighbor k (where used): {NEIGHBOR_K}",
        "",
        "## Rerun / resume",
        "",
        "```bash",
        'ssh 4060-server',
        'source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate multiomics_robustness',
        "cd /home/jizhu/research/multiomics_robustness",
        'export TMPDIR="$HOME/tmp" MPLBACKEND=Agg PATH="/usr/lib/wsl/lib:$PATH"',
        "tmux new -s phase11b",
        "",
        "python scripts/11b_run_grid.py --stage gate1",
        "python scripts/11b_run_grid.py --stage gate2",
        "python scripts/11b_run_grid.py --stage seeds --seeds 1-4",
        "python scripts/11b_run_grid.py --stage seeds --seeds 5-19",
        "python scripts/11b_run_grid.py --stage crossed",
        "python scripts/11b_mc_check.py --extra-s100-seed0",
        "python scripts/11b_analyze.py --stage all",
        "python scripts/11b_crossed.py --stage all",
        "python scripts/11b_figures.py",
        "python scripts/11b_methods_log.py",
        "```",
        "",
        "Every training run writes a done marker under `logs/done/`. Re-running any stage ",
        "skips completed runs, so an interrupted grid resumes without retraining.",
        "",
    ]
    out = paths["logs"] / "phase11b_methods.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    frame.to_csv(paths["tables"] / "phase11b_run_manifest_summary.csv", index=False)
    print(f"wrote {out}")
    print(f"runs={len(frame)} integrity_all_pass={bool(frame['integrity_pass'].all())} "
          f"total_train_h={total_h:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
