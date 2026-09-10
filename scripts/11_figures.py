#!/usr/bin/env python
"""PHASE 11 figures 1-13. Saves PNG and PDF."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import precision_recall_curve, roc_curve  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11_uncertainty import (  # noqa: E402
    DIRECTIONS,
    MODELS,
    SEEDS,
    DirectionData,
    align_probability_matrices,
    ensemble_uncertainty,
    ensure_dirs,
    load_phase10_pack,
    predictive_uncertainty,
    risk_coverage,
)

LEVEL = "l2"
DIR_LABEL = {
    "pbmc10k_to_pbmc5k": "A: PBMC10k$\\rightarrow$PBMC5k",
    "pbmc5k_to_pbmc10k": "B: PBMC5k$\\rightarrow$PBMC10k",
}
MODEL_LABEL = {"scvi_matched": "scVI_matched", "totalvi": "totalVI"}
COLOR = {"scvi_matched": "#3B6FB6", "totalvi": "#D1495B"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def save(fig, figures: Path, name: str) -> None:
    fig.savefig(figures / f"{name}.png", dpi=300)
    fig.savefig(figures / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}", flush=True)


def load_seed_arrays(paths, direction, model):
    packs, meta = {}, {}
    for seed in SEEDS:
        f = np.load(paths["arrays"] / f"probs_{direction}_{model}_seed{seed}_{LEVEL}.npz",
                    allow_pickle=False)
        packs[seed] = {"probs": f["probs"].astype(float), "classes": f["classes"]}
        if not meta:
            meta = {"y": f["y_true"].astype(str), "pred0": f["pred"].astype(str),
                    "ids": f["cell_ids"].astype(str)}
    return packs, meta


def main() -> int:
    paths = ensure_dirs()
    figures, tables = paths["figures"], paths["tables"]
    primary = pd.read_csv(tables / "uncertainty_primary.csv")
    ensemble = pd.read_csv(tables / "uncertainty_ensemble.csv")
    celltype = pd.read_csv(tables / "uncertainty_celltype_specific.csv")
    reliability = pd.read_csv(tables / "uncertainty_reliability_bins.csv")
    coverage = pd.read_csv(tables / "uncertainty_coverage.csv")
    neighbor = pd.read_csv(tables / "uncertainty_neighbor_stability.csv")
    subset = pd.read_csv(tables / "uncertainty_source_subset_instability.csv")

    cache = {}
    for direction, _, _ in DIRECTIONS:
        for model in MODELS:
            packs, meta = load_seed_arrays(paths, direction, model)
            stacked, classes = align_probability_matrices(packs)
            ens = ensemble_uncertainty(stacked, classes)
            u = predictive_uncertainty(packs[0]["probs"])["u_pred"]
            err = (meta["pred0"] != meta["y"]).astype(int)
            cache[(direction, model)] = {"u": u, "err": err, "ens": ens, "meta": meta}

    # ---- Figure 1: uncertainty taxonomy -----------------------------------
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    ax.axis("off")
    boxes = [
        (0.03, "A. Latent posterior\nq(z|x), one fixed model",
         "UNSUPPORTED for target cells\nadapted query model not saved\n(source-side scVI only)", "#BBBBBB"),
        (0.27, "B. Predictive\nsource-trained classifier",
         "PRIMARY\n$U_{pred}=1-p_{max}$\nentropy, margin", "#3B6FB6"),
        (0.51, "C. Model-seed ensemble\n5 PHASE 10 seeds",
         "$H(\\bar{p})-\\overline{H(p^{(s)})}$\nvariation ratio", "#7A9E7E"),
        (0.75, "D. Source-subset\n5 PHASE 10B subsamples",
         "variation ratio\ntrue-class prob. SD", "#D1495B"),
    ]
    for x, title, body, color in boxes:
        ax.add_patch(plt.Rectangle((x, 0.42), 0.21, 0.42, facecolor=color, alpha=0.20,
                                   edgecolor=color, linewidth=2))
        ax.text(x + 0.105, 0.78, title, ha="center", va="top", fontsize=9, fontweight="bold")
        ax.text(x + 0.105, 0.60, body, ha="center", va="center", fontsize=7.6)
    ax.text(0.5, 0.30, "kept separate throughout — never combined into one reliability index (§48)",
            ha="center", fontsize=8.5, style="italic")
    ax.text(0.5, 0.94, "PHASE 11 uncertainty taxonomy", ha="center", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0.2, 1)
    save(fig, figures, "phase11_figure1_taxonomy")

    # ---- Figure 2: reliability diagrams -----------------------------------
    rel = reliability[(reliability.scheme == "equal_frequency") & (reliability.model_seed == 0)]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="perfect calibration")
        for model in MODELS:
            r = rel[(rel.direction == direction) & (rel.model == model)]
            ax.plot(r["mean_confidence"], r["empirical_accuracy"], "o-", ms=4,
                    color=COLOR[model], label=MODEL_LABEL[model])
            e = primary[(primary.direction == direction) & (primary.model == model)
                        & (primary.model_seed == 0)]["ece_equal_frequency"].iloc[0]
            ax.text(0.05, 0.95 - 0.07 * list(MODELS).index(model),
                    f"{MODEL_LABEL[model]} ECE={e:.4f}", transform=ax.transAxes,
                    fontsize=7.5, color=COLOR[model])
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
        ax.set_xlabel("mean predicted confidence")
        ax.legend(fontsize=7.5, loc="lower right")
    axes[0].set_ylabel("empirical accuracy")
    fig.suptitle("Reliability diagrams, l2, seed 0, 15 equal-frequency bins", fontsize=10.5)
    save(fig, figures, "phase11_figure2_reliability")

    # ---- Figure 3: error-detection ROC ------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        ax.plot([0, 1], [0, 1], "k--", lw=1)
        for model in MODELS:
            c = cache[(direction, model)]
            fpr, tpr, _ = roc_curve(c["err"], c["u"])
            a = primary[(primary.direction == direction) & (primary.model == model)
                        & (primary.model_seed == 0)]["error_detection_auroc_pred"].iloc[0]
            ax.plot(fpr, tpr, color=COLOR[model], lw=1.6,
                    label=f"{MODEL_LABEL[model]} AUROC={a:.3f}")
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
        ax.set_xlabel("false positive rate"); ax.legend(fontsize=7.5, loc="lower right")
    axes[0].set_ylabel("true positive rate")
    fig.suptitle("Error detection using $U_{pred}=1-p_{max}$ (seed 0, l2)", fontsize=10.5)
    save(fig, figures, "phase11_figure3_error_roc")

    # ---- Figure 4: error-detection PR -------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        for model in MODELS:
            c = cache[(direction, model)]
            pr, rc_, _ = precision_recall_curve(c["err"], c["u"])
            a = primary[(primary.direction == direction) & (primary.model == model)
                        & (primary.model_seed == 0)]["error_detection_auprc_pred"].iloc[0]
            ax.plot(rc_, pr, color=COLOR[model], lw=1.6,
                    label=f"{MODEL_LABEL[model]} AUPRC={a:.3f}")
            ax.axhline(c["err"].mean(), color=COLOR[model], ls=":", lw=1)
        ax.set_title(f"{DIR_LABEL[direction]}\n(dotted = error prevalence baseline)", fontsize=9)
        ax.set_xlabel("recall"); ax.legend(fontsize=7.5, loc="upper right")
    axes[0].set_ylabel("precision")
    fig.suptitle("Error-detection precision-recall (seed 0, l2)", fontsize=10.5)
    save(fig, figures, "phase11_figure4_error_pr")

    # ---- Figure 5: risk-coverage ------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        for model in MODELS:
            c = cache[(direction, model)]
            rc_ = risk_coverage(c["u"], 1 - c["err"])
            step = max(1, rc_["coverage"].size // 500)
            ax.plot(rc_["coverage"][::step], rc_["risk"][::step], color=COLOR[model], lw=1.6,
                    label=f"{MODEL_LABEL[model]} AURC={rc_['aurc']:.4f}")
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
        ax.set_xlabel("coverage (fraction retained)"); ax.legend(fontsize=7.5)
    axes[0].set_ylabel("risk (1 - accuracy)")
    fig.suptitle("Risk-coverage: lower is better (seed 0, l2)", fontsize=10.5)
    save(fig, figures, "phase11_figure5_risk_coverage")

    # ---- Figure 6: macro-F1 vs coverage -----------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        for model in MODELS:
            g = (coverage[(coverage.direction == direction) & (coverage.model == model)]
                 .groupby("coverage", as_index=False)
                 .agg(macro_f1=("macro_f1", "mean"), sd=("macro_f1", "std"),
                      acc=("accuracy", "mean")))
            ax.errorbar(g["coverage"], g["macro_f1"], yerr=g["sd"], marker="o", ms=4,
                        color=COLOR[model], lw=1.5, capsize=2,
                        label=f"{MODEL_LABEL[model]} macro-F1")
            ax.plot(g["coverage"], g["acc"], "--", color=COLOR[model], lw=1, alpha=0.6,
                    label=f"{MODEL_LABEL[model]} accuracy")
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
        ax.set_xlabel("coverage retained"); ax.invert_xaxis(); ax.legend(fontsize=7)
    axes[0].set_ylabel("score (mean $\\pm$ SD over 5 seeds)")
    fig.suptitle("Selective prediction: macro-F1 and accuracy vs coverage", fontsize=10.5)
    save(fig, figures, "phase11_figure6_macrof1_coverage")

    # ---- Figure 7: correct vs incorrect distributions ---------------------
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.4), sharex=True)
    for r, (direction, _, _) in enumerate(DIRECTIONS):
        for c_, model in enumerate(MODELS):
            ax = axes[r, c_]
            c = cache[(direction, model)]
            ax.hist(c["u"][c["err"] == 0], bins=40, range=(0, 1), alpha=0.65,
                    color="#4C9F70", density=True, label="correct")
            ax.hist(c["u"][c["err"] == 1], bins=40, range=(0, 1), alpha=0.65,
                    color="#C1436D", density=True, label="incorrect")
            ax.set_title(f"{DIR_LABEL[direction]} — {MODEL_LABEL[model]}", fontsize=8.5)
            ax.legend(fontsize=7)
            if r == 1:
                ax.set_xlabel("$U_{pred} = 1 - p_{max}$")
            if c_ == 0:
                ax.set_ylabel("density")
    fig.suptitle("Predictive uncertainty, correct vs incorrect target cells (seed 0, l2)",
                 fontsize=10.5)
    save(fig, figures, "phase11_figure7_correct_incorrect")

    # ---- Figure 8: cell-type uncertainty heatmap --------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.2))
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        sub = celltype[celltype.direction == direction]
        piv = sub.pivot_table(index="cell_type", columns="model",
                              values="median_predictive_uncertainty")
        piv = piv.reindex(sorted(piv.index, key=lambda c: -piv.mean(axis=1)[c]))
        piv.columns = [MODEL_LABEL[c] for c in piv.columns]
        im = ax.imshow(piv.values, aspect="auto", cmap="magma_r", vmin=0, vmax=0.8)
        ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, fontsize=8)
        ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=7)
        ax.set_title(DIR_LABEL[direction], fontsize=9.5); ax.grid(False)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                ax.text(j, i, f"{piv.values[i, j]:.2f}", ha="center", va="center",
                        fontsize=6, color="white" if piv.values[i, j] > 0.4 else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, label="median $U_{pred}$")
    fig.suptitle("Cell-type-specific predictive uncertainty (l2, seed 0)", fontsize=10.5)
    save(fig, figures, "phase11_figure8_celltype_heatmap")

    # ---- Figure 9: source support vs uncertainty/error --------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4))
    for ax, ycol, ylabel in zip(
        axes, ["median_predictive_uncertainty", "error_rate"],
        ["median $U_{pred}$", "error rate"],
    ):
        for direction, _, _ in DIRECTIONS:
            for model in MODELS:
                s = celltype[(celltype.direction == direction) & (celltype.model == model)]
                ax.scatter(s["source_support"], s[ycol], s=26, alpha=0.75,
                           color=COLOR[model],
                           marker="o" if direction == "pbmc10k_to_pbmc5k" else "^",
                           label=f"{DIR_LABEL[direction][0]} {MODEL_LABEL[model]}")
        ax.set_xscale("log"); ax.set_xlabel("source class support (log)")
        ax.set_ylabel(ylabel)
    handles, labels = axes[0].get_legend_handles_labels()
    seen, h2, l2_ = set(), [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l); h2.append(h); l2_.append(l)
    axes[0].legend(h2, l2_, fontsize=6.5, loc="upper right")
    fig.suptitle("Source class support vs uncertainty and error (per l2 class)", fontsize=10.5)
    save(fig, figures, "phase11_figure9_support")

    # ---- Figure 10: ensemble disagreement by direction --------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3))
    data, labels, colors = [], [], []
    for direction, _, _ in DIRECTIONS:
        for model in MODELS:
            data.append(cache[(direction, model)]["ens"]["disagreement"])
            labels.append(f"{DIR_LABEL[direction][0]}\n{MODEL_LABEL[model]}")
            colors.append(COLOR[model])
    bp = axes[0].boxplot(data, labels=labels, showfliers=False, patch_artist=True)
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.5)
    axes[0].set_ylabel("ensemble disagreement  $H(\\bar{p})-\\overline{H}$")
    axes[0].set_title("Disagreement distribution", fontsize=9.5)
    e = ensemble.copy()
    e["lab"] = e["direction"].str[0].str.upper() + " " + e["model"].map(MODEL_LABEL)
    x = np.arange(len(e))
    axes[1].bar(x - 0.2, e["disagreement_error_auroc"], 0.4, label="disagreement", color="#7A9E7E")
    axes[1].bar(x + 0.2, e["variation_ratio_error_auroc"], 0.4, label="variation ratio",
                color="#C6862B")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"{DIR_LABEL[d][0]}\n{MODEL_LABEL[m]}"
                             for d, m in zip(e["direction"], e["model"])], fontsize=7.5)
    axes[1].axhline(0.5, color="k", ls="--", lw=1)
    axes[1].set_ylabel("error-detection AUROC"); axes[1].set_ylim(0.4, 1.0)
    axes[1].set_title("Ensemble signals as error detectors", fontsize=9.5)
    axes[1].legend(fontsize=7.5)
    fig.suptitle("Model-seed ensemble disagreement by transfer direction", fontsize=10.5)
    save(fig, figures, "phase11_figure10_ensemble_disagreement")

    # ---- Figure 11: neighbor stability vs correctness ---------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    for ax, (direction, _, _) in zip(axes, DIRECTIONS):
        pos, ticks = [], []
        for i, model in enumerate(MODELS):
            n = neighbor[(neighbor.direction == direction) & (neighbor.model == model)]
            for j, corr in enumerate((1, 0)):
                v = n[n.correct == corr]["neighbor_jaccard_mean"].to_numpy()
                p = i * 2.5 + j
                bp = ax.boxplot([v], positions=[p], widths=0.7, showfliers=False,
                                patch_artist=True)
                bp["boxes"][0].set_facecolor("#4C9F70" if corr else "#C1436D")
                bp["boxes"][0].set_alpha(0.6)
                pos.append(p); ticks.append(f"{MODEL_LABEL[model]}\n{'correct' if corr else 'wrong'}")
        ax.set_xticks(pos); ax.set_xticklabels(ticks, fontsize=7)
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
    axes[0].set_ylabel(f"neighbor-set Jaccard across 5 seeds (k={15})")
    fig.suptitle("Neighborhood stability vs prediction correctness", fontsize=10.5)
    save(fig, figures, "phase11_figure11_neighbor_stability")

    # ---- Figure 12: full-source uncertainty vs subset instability ---------
    assoc = pd.read_csv(tables / "uncertainty_subset_associations.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), sharey=True)
    for ax, model in zip(axes, MODELS):
        s = subset[subset.model == model]
        ax.scatter(s["full_source_predictive_uncertainty"],
                   s["subset_true_class_probability_sd"], s=5, alpha=0.25, color=COLOR[model])
        rho = assoc[(assoc.model == model) & (assoc.x == "full_source_u_pred")
                    & (assoc.y == "subset_true_class_probability_sd")]["spearman_rho"].iloc[0]
        ax.set_title(f"{MODEL_LABEL[model]}  Spearman $\\rho$={rho:.3f}", fontsize=9.5)
        ax.set_xlabel("full-source $U_{pred}$ (PHASE 10 Direction A)")
    axes[0].set_ylabel("subset true-class probability SD\n(5 PHASE 10B source subsamples)")
    fig.suptitle("Does the full model know which cells are fragile to source resampling?",
                 fontsize=10.5)
    save(fig, figures, "phase11_figure12_subset_instability")

    # ---- Figure 13: source vs target uncertainty --------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.3), sharey=True)
    src_stats = []
    for ax, (direction, sb, tb) in zip(axes, DIRECTIONS):
        data_d = DirectionData(direction, sb, tb)
        sl = data_d.eval_slice(LEVEL)
        pos, ticks, cols = [], [], []
        for i, model in enumerate(MODELS):
            pack = load_phase10_pack(direction, model, 0)
            z_s = pack["source_latent"][sl["source_index"]]
            y_s = sl["y_source"]
            # out-of-sample source uncertainty via 5-fold CV (in-sample would be biased)
            u_src = np.empty(y_s.size)
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
            for tr, te in skf.split(z_s, y_s):
                clf = LogisticRegression(max_iter=2000, solver="lbfgs", random_state=0)
                clf.fit(z_s[tr], y_s[tr])
                u_src[te] = 1.0 - clf.predict_proba(z_s[te]).max(axis=1)
            u_tgt = cache[(direction, model)]["u"]
            for j, (v, lab) in enumerate(((u_src, "source\n(5-fold CV)"), (u_tgt, "target"))):
                p = i * 2.5 + j
                bp = ax.boxplot([v], positions=[p], widths=0.7, showfliers=False,
                                patch_artist=True)
                bp["boxes"][0].set_facecolor(COLOR[model])
                bp["boxes"][0].set_alpha(0.35 if j == 0 else 0.75)
                pos.append(p); ticks.append(f"{MODEL_LABEL[model]}\n{lab}")
            src_stats.append(
                {"direction": direction, "model": model,
                 "median_u_source_cv": float(np.median(u_src)),
                 "iqr_u_source_cv": float(np.percentile(u_src, 75) - np.percentile(u_src, 25)),
                 "median_u_target": float(np.median(u_tgt)),
                 "iqr_u_target": float(np.percentile(u_tgt, 75) - np.percentile(u_tgt, 25))}
            )
            del pack
        ax.set_xticks(pos); ax.set_xticklabels(ticks, fontsize=7)
        ax.set_title(DIR_LABEL[direction], fontsize=9.5)
        del data_d
    axes[0].set_ylabel("$U_{pred}$")
    fig.suptitle("Source vs target predictive uncertainty (seed 0, l2)", fontsize=10.5)
    save(fig, figures, "phase11_figure13_source_vs_target")
    pd.DataFrame(src_stats).to_csv(tables / "uncertainty_source_vs_target.csv", index=False)

    print("\nall figures written", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
