#!/usr/bin/env python3
"""PART 0–1: experiment design ledger + transduction/leakage audit (code-grounded)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "revision_round2_methodological_fixes"
AUDIT = OUT / "audit"
LOGS = OUT / "logs"


ROWS = [
    {
        "experiment": "representation_benchmark",
        "dataset": "PBMC_CITE_combined_inner",
        "source_dataset": "PBMC10k+PBMC5k (joint)",
        "target_dataset": "same joint object (held-out split for classifier)",
        "representation_model": "scVI_default / scVI_matched / totalVI (+ MOFA+/PCA_RNA in PHASE7)",
        "representation_training_cells": "all 10849 combined cells (unsupervised; internal ELBO train_size=0.9 unrelated to obs.split)",
        "query_adaptation_cells": "N/A for VAE reps; SCANVI annotation surgery uses all query cells unlabeled",
        "classifier_training_cells": "obs.split==train within high-confidence annotated cells",
        "evaluation_cells": "obs.split==test within high-confidence annotated cells",
        "labels_used_for_training": "none for VAE; SCANVI uses reference labels then transfers; classifier uses transferred query labels on train split",
        "labels_used_for_evaluation": "PHASE6 transferred cell_type_l1/l2 (high-confidence)",
        "label_source": "Seurat v4 PBMC CITE-seq reference via SCANVI/scArches",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "N/A (VAE trained on all cells jointly — unsupervised transduction over full matrix)",
        "biological_unit": "cell",
        "computational_unit": "model seed (primary seed 0)",
        "seed_definition": "torch/numpy/scvi model seed",
        "run_replicate_definition": "single primary run",
        "bootstrap_unit": "test cells (N_BOOT=500) for Delta F1",
        "metric": "macro-F1 (l2 high-conf), ASW/purity secondary",
        "primary_statistical_unit": "target/test cell (conditional bootstrap)",
        "transductive_or_inductive": "unsupervised_transductive_representation + inductive_classifier",
        "artifact_source": "results/tables/biological_representation_metrics.csv; results/embeddings/*_seed0_*; scripts/06_run_phase6.py; src/experiments/phase6*.py",
    },
    {
        "experiment": "protein_corruption",
        "dataset": "PBMC_CITE_combined_inner",
        "source_dataset": "same",
        "target_dataset": "same (train/test split)",
        "representation_model": "totalVI retrained under corruption; scVI_matched latents reused (seed0)",
        "representation_training_cells": "all cells with corrupted protein matrix (training-time degradation)",
        "query_adaptation_cells": "N/A",
        "classifier_training_cells": "high-conf train split",
        "evaluation_cells": "high-conf test split",
        "labels_used_for_training": "none for VAE; classifier uses PHASE6 labels on train",
        "labels_used_for_evaluation": "PHASE6 high-confidence labels",
        "label_source": "PHASE6 Seurat-v4 transfer",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "N/A; corruption applied to protein entries before totalVI retrain",
        "biological_unit": "cell",
        "computational_unit": "(corruption_fraction, perturbation_seed); model_seed fixed 0",
        "seed_definition": "perturbation/mask seed for entry zeroing",
        "run_replicate_definition": "5 perturbation seeds × corruption fractions",
        "bootstrap_unit": "test cells",
        "metric": "macro-F1 Delta vs scVI; corruption recovery",
        "primary_statistical_unit": "corruption-mask seed (descriptive) + cell bootstrap",
        "transductive_or_inductive": "training_time_protein_degradation (not pure test-time shift)",
        "artifact_source": "results/perturbations/; results/tables/protein_sparsity_*.csv; src/experiments/phase8_sparsity.py",
    },
    {
        "experiment": "cross_dataset_Direction_A",
        "dataset": "PBMC_CITE_combined_inner_annotated",
        "source_dataset": "PBMC10k",
        "target_dataset": "PBMC5k",
        "representation_model": "scVI_matched / totalVI + scArches query adaptation",
        "representation_training_cells": "all source-batch cells only",
        "query_adaptation_cells": "all target-batch cells (unlabeled features)",
        "classifier_training_cells": "high-conf source cells in eligible classes",
        "evaluation_cells": "high-conf target cells in eligible classes",
        "labels_used_for_training": "source high-conf labels for classifier only",
        "labels_used_for_evaluation": "target high-conf labels",
        "label_source": "PHASE6 transferred labels",
        "target_labels_seen_by_representation_model": "No (labels not passed to VAE/adaptation)",
        "target_features_seen_during_adaptation": "Yes — unsupervised scArches query adaptation",
        "biological_unit": "cell (dataset shift; not donor)",
        "computational_unit": "(direction, model, model_seed)",
        "seed_definition": "model seed 0–4 (PHASE10); 0–19 (PHASE11B primary)",
        "run_replicate_definition": "one run per seed",
        "bootstrap_unit": "target eval cells (PHASE10); also seed bootstrap in PHASE11B",
        "metric": "l2 macro-F1; uncertainty metrics in PHASE11/11B",
        "primary_statistical_unit": "model seed (20-seed primary for transfer claims)",
        "transductive_or_inductive": "unsupervised_transductive_query_adaptation + inductive_classifier",
        "artifact_source": "results/phase10_cross_dataset/; results/phase11b_uncertainty_robustness/; src/experiments/phase10_cross_dataset.py",
    },
    {
        "experiment": "cross_dataset_Direction_B",
        "dataset": "PBMC_CITE_combined_inner_annotated",
        "source_dataset": "PBMC5k",
        "target_dataset": "PBMC10k",
        "representation_model": "scVI_matched / totalVI + scArches query adaptation",
        "representation_training_cells": "all source-batch cells only",
        "query_adaptation_cells": "all target-batch cells (unlabeled features)",
        "classifier_training_cells": "high-conf source cells in eligible classes",
        "evaluation_cells": "high-conf target cells in eligible classes",
        "labels_used_for_training": "source high-conf labels for classifier only",
        "labels_used_for_evaluation": "target high-conf labels",
        "label_source": "PHASE6 transferred labels",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "Yes — unsupervised scArches query adaptation",
        "biological_unit": "cell",
        "computational_unit": "(direction, model, model_seed)",
        "seed_definition": "model seed 0–4 / 0–19",
        "run_replicate_definition": "one run per seed",
        "bootstrap_unit": "target eval cells / seed bootstrap",
        "metric": "l2 macro-F1 (+ uncertainty)",
        "primary_statistical_unit": "model seed",
        "transductive_or_inductive": "unsupervised_transductive_query_adaptation + inductive_classifier",
        "artifact_source": "results/phase10_cross_dataset/; results/phase11b_uncertainty_robustness/",
    },
    {
        "experiment": "source_size_analysis",
        "dataset": "PBMC_CITE_combined_inner_annotated",
        "source_dataset": "PBMC10k subsampled to 3994 (Direction A size-match)",
        "target_dataset": "PBMC5k",
        "representation_model": "scVI_matched / totalVI + scArches",
        "representation_training_cells": "subsampled source cells (label-blind draw)",
        "query_adaptation_cells": "all target cells unlabeled",
        "classifier_training_cells": "high-conf source eligible",
        "evaluation_cells": "high-conf target eligible (common class set)",
        "labels_used_for_training": "source labels for classifier only",
        "labels_used_for_evaluation": "target labels",
        "label_source": "PHASE6",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "Yes",
        "biological_unit": "cell",
        "computational_unit": "(subsample_seed, model)",
        "seed_definition": "subsample seed 0–4; MODEL_SEED=0 for new runs",
        "run_replicate_definition": "5 subsample seeds",
        "bootstrap_unit": "target cells",
        "metric": "macro-F1 under common eligible classes",
        "primary_statistical_unit": "subsample seed",
        "transductive_or_inductive": "unsupervised_transductive_query_adaptation + inductive_classifier",
        "artifact_source": "results/phase10b_source_size/; src/experiments/phase10b_source_size.py",
    },
    {
        "experiment": "20seed_uncertainty",
        "dataset": "PBMC_CITE_combined_inner_annotated",
        "source_dataset": "PBMC10k or PBMC5k by direction",
        "target_dataset": "reciprocal batch",
        "representation_model": "scVI_matched / totalVI + persisted scArches adaptation (PHASE11B)",
        "representation_training_cells": "source batch (full or crossed subsets)",
        "query_adaptation_cells": "all unlabeled target cells",
        "classifier_training_cells": "source eligible high-conf",
        "evaluation_cells": "target eligible high-conf",
        "labels_used_for_training": "source labels for classifier; class_plan uses target counts for eligibility design",
        "labels_used_for_evaluation": "target high-conf labels",
        "label_source": "PHASE6",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "Yes",
        "biological_unit": "cell",
        "computational_unit": "model seed 0–19 (full-source); crossed 5 subsets × 5 seeds",
        "seed_definition": "model seed",
        "run_replicate_definition": "one run per seed (crossed design reused by PHASE12A as rep0)",
        "bootstrap_unit": "seed-level (10000) AND target-cell (500) kept separate",
        "metric": "macro-F1, error AUROC, NLL/Brier/ECE, AURC",
        "primary_statistical_unit": "model seed",
        "transductive_or_inductive": "unsupervised_transductive_query_adaptation + inductive_classifier",
        "artifact_source": "results/phase11b_uncertainty_robustness/; NOTE: PHASE11 was 5 seeds only",
    },
    {
        "experiment": "variance_decomposition",
        "dataset": "PBMC Direction A crossed design",
        "source_dataset": "PBMC10k source subsets (5)",
        "target_dataset": "PBMC5k",
        "representation_model": "scVI_matched and totalVI SEPARATELY (+ delta cube)",
        "representation_training_cells": "source subset cells",
        "query_adaptation_cells": "all unlabeled target",
        "classifier_training_cells": "source eligible",
        "evaluation_cells": "target eligible",
        "labels_used_for_training": "source labels for classifier",
        "labels_used_for_evaluation": "target labels",
        "label_source": "PHASE6",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "Yes",
        "biological_unit": "cell (metric is scalar per run)",
        "computational_unit": "5 subsets × 5 seeds × 2 run replicates × 2 models",
        "seed_definition": "model seed",
        "run_replicate_definition": "rep0=PHASE11B reuse; rep1=PHASE12A retrain",
        "bootstrap_unit": "N/A (MoM ANOVA on run scalars)",
        "metric": "macro_f1, error_auroc, nll, brier, ece, aurc, ...",
        "primary_statistical_unit": "run (subset×seed×replicate) within model; delta separately",
        "transductive_or_inductive": "unsupervised_transductive_query_adaptation + inductive_classifier",
        "artifact_source": "results/phase12a_variance_replication/; src/experiments/phase12a_stats.py (MoM; per-model + delta)",
    },
    {
        "experiment": "Lawlor_external_validation",
        "dataset": "Lawlor Baseline SNG shared genes/proteins",
        "source_dataset": "8 source donors (per fold)",
        "target_dataset": "2 target donors (per fold)",
        "representation_model": "scVI_matched / totalVI + scArches",
        "representation_training_cells": "all source-donor cells (labeled+unlabeled)",
        "query_adaptation_cells": "all target-donor cells unlabeled features",
        "classifier_training_cells": "source labeled cells in eligible classes",
        "evaluation_cells": "target labeled cells in eligible classes present in both",
        "labels_used_for_training": "source author labels for classifier only",
        "labels_used_for_evaluation": "target author labels",
        "label_source": "author protein-gated author_cell_type",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "Yes — unsupervised query adaptation",
        "biological_unit": "donor (10 donors; folds are 2-donor pairs)",
        "computational_unit": "5 folds × 10 seeds × 2 models = 100 runs",
        "seed_definition": "model seed 0–9",
        "run_replicate_definition": "one run per fold×model×seed",
        "bootstrap_unit": "historically none as primary; revision uses donor-level deltas",
        "metric": "macro-F1 primary; predictive/latent AUROC; NLL/Brier/ECE; AURC",
        "primary_statistical_unit": "HISTORICAL: fold×seed pair (n=50) descriptive; REVISION PRIMARY: donor",
        "transductive_or_inductive": "donor-held-out label evaluation after unsupervised query adaptation",
        "artifact_source": "results/phase12b_external_validation/formal_validation/",
    },
    {
        "experiment": "selective_prediction",
        "dataset": "PBMC transfer (11/11B) and Lawlor (12B)",
        "source_dataset": "(inherits parent experiment)",
        "target_dataset": "(inherits parent experiment)",
        "representation_model": "(inherits)",
        "representation_training_cells": "(inherits)",
        "query_adaptation_cells": "(inherits)",
        "classifier_training_cells": "(inherits)",
        "evaluation_cells": "target eval cells ranked by U_pred=1-p_max",
        "labels_used_for_training": "none for ranking",
        "labels_used_for_evaluation": "target labels for risk/accuracy",
        "label_source": "(inherits)",
        "target_labels_seen_by_representation_model": "No",
        "target_features_seen_during_adaptation": "(inherits)",
        "biological_unit": "(inherits)",
        "computational_unit": "(inherits)",
        "seed_definition": "(inherits)",
        "run_replicate_definition": "(inherits)",
        "bootstrap_unit": "(inherits)",
        "metric": "FULL AURC = mean risk over all prefixes; discrete coverage table at {1.0..0.5}",
        "primary_statistical_unit": "(inherits)",
        "transductive_or_inductive": "(inherits)",
        "artifact_source": "src/experiments/phase11_uncertainty.py::risk_coverage; COVERAGE_GRID for tables only",
    },
]


def write_ledger() -> Path:
    df = pd.DataFrame(ROWS)
    path = AUDIT / "experiment_design_ledger.csv"
    df.to_csv(path, index=False)
    return path


def write_transduction_audit() -> Path:
    md = """# Transduction / leakage audit (code-grounded)

**Rule applied:** unlabeled target features in scArches query adaptation = **unsupervised transductive adaptation**, not automatic label leakage.

## Summary verdict

| Experiment | Target used to train representation? | Unlabeled target in scArches? | Target labels before final eval? | Genuine label leakage? |
|---|---|---|---|---|
| Representation benchmark (PHASE6 VAEs) | Yes — all query cells jointly | N/A (no query surgery for VAEs) | Labels used for classifier/eval after transfer annotation | **No** for VAE; annotation is separate transfer |
| Protein corruption (PHASE8) | Yes — all cells with corrupted protein | No | PHASE6 labels for clf/eval | **No** |
| Cross-dataset A/B (PHASE10/11B) | **No** (source batch only) | **Yes** | Class eligibility uses target label counts; clf fit source-only | **No label leakage into model**; evaluation-design uses target labels |
| Source-size (PHASE10B) | **No** (subsampled source) | **Yes** | Same as PHASE10 | **No** (same caveat) |
| Variance (PHASE12A) | **No** (source subset) | **Yes** | Same | **No** |
| Lawlor (PHASE12B) | **No** (source donors) | **Yes** | Eligibility rule uses target supports every fold; clf source-only | **No** (same caveat) |
| Selective prediction | inherits | inherits | labels only at risk scoring | **No** |

## Seven-question checklist (Lawlor + transfer)

1. **Was the target dataset used to train the representation?**  
   Transfer/Lawlor: **No** (source only). PHASE6/8 joint VAEs: **Yes** (all cells).

2. **Was unlabeled target data used in scArches query adaptation?**  
   Transfer/Lawlor/11B/12A: **Yes**. This is unsupervised transduction.

3. **Were target labels used anywhere before final evaluation?**  
   - Training/adaptation/classifier fit: **No**.  
   - **Evaluation design:** eligible class sets use target label counts (`class_plan` / Lawlor `>=20 source & >=10 target in every fold`).  
   This is **target-label-informed evaluation design**, not parameter leakage.

4. **Were feature-selection steps fitted on train/source only or on all cells?**  
   VAE transfer/Lawlor: **no HVG** (fixed shared gene universe).  
   PHASE6 annotation HVGs: reference-only.  
   PHASE7 MOFA+/PCA_RNA: fitted on **all combined cells** (joint unsupervised).

5. **Were PCA / scaling / HVG / normalization parameters fitted on evaluation cells?**  
   VAE pipelines: **No** such steps.  
   PHASE7 PCA/MOFA: **Yes** (all cells) — note for baseline fairness discussions.

6. **Was classifier preprocessing fitted only on classifier-training cells?**  
   **Yes vacuously** — logistic regression on raw latents; no StandardScaler.

7. **Were class eligibility thresholds defined using target labels?**  
   **Yes** for transfer and Lawlor.

## Correct manuscript terminology

Use: **donor-held-out / dataset-held-out label evaluation after unsupervised query adaptation**.

Do **not** imply fully inductive prediction on completely unseen target feature distributions without adaptation.

Do **not** call valid unsupervised transduction “leakage.”

## Hard-stop status

- Genuine target-label leakage into representation training/adaptation/classifier fit: **NOT found**.  
- Source/target embeddings mixed (historical PHASE10 bug): corrected in prior phases; Lawlor/11B use post-adaptation coordinates consistently.  
- Continue with statistical corrections.

## Key code citations

- `src/experiments/phase10_cross_dataset.py` — `query_adapt`, `class_plan`, `target_labels_used_in_fitting: False`
- `src/experiments/phase12b_external.py` — donor folds, `evaluate_run`, unlabeled query adapt
- `src/experiments/phase12a_stats.py` — per-model MoM variance (not pooled without model term)
- `src/experiments/phase11_uncertainty.py` — `risk_coverage` full-prefix AURC; `COVERAGE_GRID` for tables
"""
    path = AUDIT / "transduction_leakage_audit.md"
    path.write_text(md)
    return path


def main() -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    ledger = write_ledger()
    audit = write_transduction_audit()
    (LOGS / "part0_part1_done.txt").write_text(
        f"ledger={ledger}\naudit={audit}\n"
    )
    print("Wrote", ledger)
    print("Wrote", audit)


if __name__ == "__main__":
    main()
