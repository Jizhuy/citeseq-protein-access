#!/usr/bin/env python
"""PHASE 11B figures 1-14. Saves PNG and PDF.

Each figure is skipped with a notice if its input table is absent, so the script
can be run repeatedly as stages complete.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11b_analysis import DirectionEval, load_probabilities  # noqa: E402
from src.experiments.phase11b_robustness import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    MODELS,
    SEEDS,
    ensure_dirs,
    experiment_tag,
)

DIR_LABEL = {"direction_A": "A: PBMC10k$\\rightarrow$PBMC5k",
             "direction_B": "B: PBMC5k$\\rightarrow$PBMC10k"}
MODEL_LABEL = {"scvi_matched": "scVI_matched", "totalvi": "totalVI"}
COLOR = {"scvi_matched": "#3B6FB6", "totalvi": "#D1495B"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
PATHS = ensure_dirs()
TAB, FIG = PATHS["tables"], PATHS["figures"]


def save(fig, name: str) -> None:
    fig.savefig(FIG / f"{name}.png", dpi=300)
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}", flush=True)


def read(name: str) -> pd.DataFrame | None:
    path = TAB / name
    if not path.exists():
        print(f"  skip: {name} not found", flush=True)
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
def fig1_design() -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off"); ax.grid(False)
    ax.set_title("Figure 1. PHASE 11B experimental design", fontsize=11, weight="bold")

    def box(x, y, w, h, text, fc):
        ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec="#333", lw=1.0, alpha=0.85))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.2)

    box(0.2, 8.1, 9.6, 1.3,
        "PART I  Prospective 20-seed replication\n"
        "2 directions x 2 models x 20 model seeds = 80 runs\n"
        "each run persists source model AND scArches-adapted query model", "#DCE9F7")
    box(0.2, 6.2, 4.6, 1.6,
        "Standardised per-run pipeline\n1. train source   2. save source\n"
        "3. query adapt   4. save adapted model\n"
        "5. encode source+target with ADAPTED model", "#EAF3E4")
    box(5.2, 6.2, 4.6, 1.6,
        "Recovered in PHASE 11B, missing in PHASE 11\n"
        "target q(z|x) posterior variance\n"
        "$U_{latent}=\\mathrm{mean}_j\\,\\sigma^2_{ij}$\n"
        "from get_latent_representation(return_dist=True)", "#F7E9DC")
    box(0.2, 4.3, 9.6, 1.5,
        "PART III  Crossed design (Direction A only)\n"
        "5 PHASE 10B source subsets x 5 model seeds x 2 models = 50 runs\n"
        "separates source-data sampling from model-initialisation stochasticity", "#F2E4F7")
    box(0.2, 2.4, 3.0, 1.5, "Model-seed\nvariability\n(fixed source subset)", "#DCE9F7")
    box(3.5, 2.4, 3.0, 1.5, "Source-subset\nvariability\n(fixed model seed)", "#F7E9DC")
    box(6.8, 2.4, 3.0, 1.5, "Interaction /\nresidual\n(n=1 per cell)", "#E8E8E8")
    box(0.2, 0.4, 9.6, 1.4,
        "PART IV-V  Does full-source uncertainty predict BOTH instability sources?\n"
        "Selective prediction stability across 20 seeds + cell-type retention imbalance", "#EAF3E4")
    for x in (1.7, 5.0, 8.3):
        ax.add_patch(FancyArrowPatch((x, 4.3), (x, 3.9), arrowstyle="-|>",
                                     mutation_scale=11, color="#333"))
    save(fig, "figure01_experimental_design")


def fig2_auroc_distribution(primary: pd.DataFrame) -> None:
    dirs = sorted(primary["direction"].unique())
    fig, axes = plt.subplots(1, len(dirs), figsize=(4.6 * len(dirs), 4.0), squeeze=False)
    for ax, d in zip(axes[0], dirs):
        sub = primary[primary["direction"] == d]
        data, labels, colors = [], [], []
        for m in MODELS:
            v = sub[sub["model"] == m]["error_auroc_predictive"].dropna().to_numpy()
            if v.size:
                data.append(v); labels.append(MODEL_LABEL[m]); colors.append(COLOR[m])
        bp = ax.boxplot(data, labels=labels, widths=0.5, patch_artist=True,
                        medianprops={"color": "black"})
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.45)
        for i, (v, c) in enumerate(zip(data, colors), start=1):
            ax.scatter(np.random.default_rng(0).normal(i, 0.055, v.size), v,
                       s=14, color=c, zorder=3, edgecolor="white", linewidth=0.4)
        ax.set_title(DIR_LABEL[d]); ax.set_ylabel("error-detection AUROC ($U_{pred}$)")
        ax.text(0.02, 0.02, f"n = {len(data[0])} seeds", transform=ax.transAxes, fontsize=7.5)
    fig.suptitle("Figure 2. 20-seed distribution of predictive error-detection AUROC",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure02_20seed_auroc_distribution")


def fig3_paired_deltas(delta: pd.DataFrame, summary: pd.DataFrame | None) -> None:
    metrics = ["error_auroc_predictive", "error_auprc_predictive", "ece_equal_frequency",
               "brier", "nll", "aurc", "macro_f1"]
    metrics = [m for m in metrics if m in set(delta["metric"])]
    dirs = sorted(delta["direction"].unique())
    fig, axes = plt.subplots(1, len(dirs), figsize=(5.4 * len(dirs), 4.4), squeeze=False)
    for ax, d in zip(axes[0], dirs):
        for i, metric in enumerate(metrics):
            v = delta[(delta["direction"] == d) & (delta["metric"] == metric)]["delta"].to_numpy()
            ax.scatter(np.random.default_rng(1).normal(i, 0.07, v.size), v, s=16,
                       color="#666", alpha=0.65, edgecolor="none")
            ax.plot([i - 0.28, i + 0.28], [v.mean()] * 2, color="#D1495B", lw=2.2)
            if summary is not None:
                s = summary[(summary["direction"] == d) & (summary["metric"] == metric)]
                if len(s):
                    r = s.iloc[0]
                    ax.plot([i, i], [r["seed_bootstrap_ci_low"], r["seed_bootstrap_ci_high"]],
                            color="#D1495B", lw=1.3)
        ax.axhline(0, color="black", lw=0.9, ls="--")
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels([m.replace("_predictive", "").replace("_equal_frequency", "")
                            for m in metrics], rotation=35, ha="right", fontsize=7.5)
        ax.set_ylabel("$\\Delta$ = totalVI $-$ scVI_matched")
        ax.set_title(DIR_LABEL[d])
    fig.suptitle("Figure 3. Paired per-seed totalVI $-$ scVI deltas (red: mean and seed-bootstrap 95% CI)",
                 fontsize=10.5, weight="bold")
    fig.tight_layout()
    save(fig, "figure03_paired_seed_deltas")


def _percell_pack(direction: str, model: str, seed: int):
    tag = experiment_tag(direction, model, seed, None)
    if not (PATHS["probabilities"] / direction / f"{tag}_probs.npz").exists():
        return None
    pack = load_probabilities(direction, tag)
    return pack


def fig4_5_latent(data_cache: dict) -> None:
    seed = 0
    rows = []
    for direction in ("direction_A", "direction_B"):
        if direction not in data_cache:
            data_cache[direction] = DirectionEval(direction)
        y = data_cache[direction].slice("l2")["y_target"]
        for model in MODELS:
            pack = _percell_pack(direction, model, seed)
            if pack is None:
                continue
            pred = np.asarray(pack["classes"], dtype=object)[np.argmax(pack["probs"], axis=1)]
            rows.append({"direction": direction, "model": model,
                         "u_latent": pack["u_latent"],
                         "u_pred": 1.0 - pack["probs"].max(axis=1),
                         "correct": pred == y})
    if not rows:
        print("  skip: no probability packs for figures 4-5", flush=True)
        return

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6))
    for ax, r in zip(axes.ravel(), rows):
        ok, bad = r["u_latent"][r["correct"]], r["u_latent"][~r["correct"]]
        bins = np.linspace(min(ok.min(), bad.min()), np.quantile(r["u_latent"], 0.995), 45)
        ax.hist(ok, bins=bins, alpha=0.6, density=True, color="#3B8C6E", label="correct")
        ax.hist(bad, bins=bins, alpha=0.6, density=True, color="#D1495B", label="incorrect")
        ax.axvline(np.median(ok), color="#3B8C6E", ls="--", lw=1.2)
        ax.axvline(np.median(bad), color="#D1495B", ls="--", lw=1.2)
        ax.set_title(f"{DIR_LABEL[r['direction']]}  {MODEL_LABEL[r['model']]}", fontsize=9)
        ax.set_xlabel("$U_{latent}$ = mean posterior variance"); ax.set_ylabel("density")
        ax.legend(fontsize=7.5)
    fig.suptitle("Figure 4. Target latent posterior uncertainty: correct vs incorrect (seed 0)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure04_latent_uncertainty_correct_incorrect")

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6))
    for ax, r in zip(axes.ravel(), rows):
        from scipy.stats import spearmanr
        rho = spearmanr(r["u_pred"], r["u_latent"]).statistic
        ax.scatter(r["u_pred"], r["u_latent"], s=4, alpha=0.22,
                   color=COLOR[r["model"]], edgecolor="none")
        ax.set_xlabel("$U_{pred}=1-p_{max}$"); ax.set_ylabel("$U_{latent}$")
        ax.set_title(f"{DIR_LABEL[r['direction']]}  {MODEL_LABEL[r['model']]}\n"
                     f"Spearman $\\rho$ = {rho:.3f}", fontsize=8.5)
    fig.suptitle("Figure 5. Predictive versus latent posterior uncertainty (seed 0)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure05_predictive_vs_latent_uncertainty")


def fig6_convergence(conv: pd.DataFrame) -> None:
    nested = conv[conv["seed_subset_definition"].str.startswith("nested")]
    random = conv[conv["seed_subset_definition"].str.startswith("random")]
    metrics = [("macro_f1", "ensemble macro-F1"), ("error_auroc", "error AUROC"),
               ("ece", "ECE"), ("brier", "Brier"), ("aurc", "AURC"), ("nll", "NLL")]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.4))
    for ax, (col, label) in zip(axes.ravel(), metrics):
        for d in sorted(nested["direction"].unique()):
            for m in MODELS:
                s = nested[(nested["direction"] == d) & (nested["model"] == m)].sort_values("ensemble_size")
                if s.empty:
                    continue
                ax.plot(s["ensemble_size"], s[col], marker="o", ms=3.6,
                        color=COLOR[m], ls="-" if d == "direction_A" else "--",
                        lw=1.4, label=f"{DIR_LABEL[d][0]} {MODEL_LABEL[m]}")
                rs = random[(random["direction"] == d) & (random["model"] == m)]
                if not rs.empty:
                    ax.scatter(rs["ensemble_size"], rs[col], s=9, color=COLOR[m],
                               alpha=0.32, edgecolor="none", zorder=1)
        ax.set_xlabel("ensemble size (number of model seeds)"); ax.set_ylabel(label)
        ax.set_xticks([1, 2, 3, 5, 10, 15, 20])
    axes[0, 0].legend(fontsize=6.3, ncol=2)
    fig.suptitle("Figure 6. Ensemble-size convergence (lines: nested seed sets; "
                 "faint points: random subsets)", fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure06_ensemble_convergence")


def _heat(ax, mat, title, cmap, fmt="{:.3f}", center=None):
    if center is not None:
        lim = np.nanmax(np.abs(mat - center))
        im = ax.imshow(mat, cmap=cmap, vmin=center - lim, vmax=center + lim)
    else:
        im = ax.imshow(mat, cmap=cmap)
    ax.set_xticks(range(len(CROSSED_MODEL_SEEDS)), CROSSED_MODEL_SEEDS)
    ax.set_yticks(range(len(CROSSED_SUBSET_SEEDS)), CROSSED_SUBSET_SEEDS)
    ax.set_xlabel("model seed"); ax.set_ylabel("source subset seed")
    ax.set_title(title, fontsize=9); ax.grid(False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, fmt.format(mat[i, j]), ha="center", va="center", fontsize=6.8,
                    color="black")
    return im


def fig7_8_heatmaps() -> None:
    path = TAB / "phase11b_crossed_matrices.npz"
    if not path.exists():
        print("  skip: crossed matrices not found", flush=True)
        return
    z = np.load(path)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    for ax, m in zip(axes, MODELS):
        key = f"{m}__macro_f1"
        if key in z:
            im = _heat(ax, z[key], f"{MODEL_LABEL[m]} macro-F1", "viridis")
            fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("Figure 7. Crossed source-subset x model-seed l2 macro-F1 (Direction A)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure07_crossed_macrof1_heatmaps")

    keys = [k for k in ("delta__macro_f1", "delta__error_auroc") if k in z]
    if keys:
        fig, axes = plt.subplots(1, len(keys), figsize=(4.9 * len(keys), 4.2), squeeze=False)
        for ax, k in zip(axes[0], keys):
            im = _heat(ax, z[k], f"$\\Delta$ {k.split('__')[1]}", "RdBu_r",
                       fmt="{:+.3f}", center=0.0)
            fig.colorbar(im, ax=ax, fraction=0.046)
        fig.suptitle("Figure 8. $\\Delta$ = totalVI $-$ scVI_matched across the crossed grid",
                     fontsize=11, weight="bold")
        fig.tight_layout()
        save(fig, "figure08_crossed_delta_heatmap")


def fig9_variance(vc: pd.DataFrame) -> None:
    models = [m for m in (*MODELS, "delta_totalvi_minus_scvi") if m in set(vc["model"])]
    metrics = ["macro_f1", "error_auroc", "ece", "aurc"]
    fig, axes = plt.subplots(1, len(models), figsize=(4.2 * len(models), 4.0), squeeze=False)
    for ax, m in zip(axes[0], models):
        sub = vc[vc["model"] == m].set_index("metric").reindex(metrics).dropna(how="all")
        if sub.empty:
            continue
        bottom = np.zeros(len(sub))
        for col, lab, c in (("source_fraction", "source subset", "#E08A3C"),
                            ("seed_fraction", "model seed", "#3B6FB6"),
                            ("interaction_fraction", "interaction / residual", "#9E9E9E")):
            vals = sub[col].to_numpy(dtype=float)
            ax.bar(range(len(sub)), vals, bottom=bottom, label=lab, color=c, width=0.62)
            bottom += np.nan_to_num(vals)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub.index, rotation=30, ha="right", fontsize=7.5)
        ax.set_ylim(0, 1); ax.set_ylabel("fraction of total variance")
        ax.set_title(MODEL_LABEL.get(m, "$\\Delta$ totalVI $-$ scVI"), fontsize=9)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Figure 9. Two-factor variance decomposition (non-negative truncation)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure09_variance_components")


def fig10_celltype(ct: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, len(MODELS), figsize=(6.6 * len(MODELS), 6.2), squeeze=False)
    for ax, m in zip(axes[0], MODELS):
        sub = ct[ct["model"] == m].sort_values("source_variance_fraction")
        if sub.empty:
            continue
        idx = np.arange(len(sub))
        ax.barh(idx, sub["source_variance_fraction"], color="#E08A3C", label="source subset")
        ax.barh(idx, sub["seed_variance_fraction"], left=sub["source_variance_fraction"],
                color="#3B6FB6", label="model seed")
        ax.barh(idx, sub["interaction_fraction"],
                left=sub["source_variance_fraction"] + sub["seed_variance_fraction"],
                color="#9E9E9E", label="interaction / residual")
        ax.set_yticks(idx); ax.set_yticklabels(sub["cell_type"], fontsize=7)
        ax.set_xlabel("median per-cell variance fraction"); ax.set_xlim(0, 1)
        ax.set_title(MODEL_LABEL[m], fontsize=9)
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Figure 10. Cell-type-specific variance contributions (true-class probability)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure10_celltype_variance")


def fig11_12_instability(cell: pd.DataFrame) -> None:
    from scipy.stats import spearmanr
    for num, inst, title in ((11, "source_subset_trueclass_var", "source-subset instability"),
                             (12, "model_seed_trueclass_var", "model-seed instability")):
        fig, axes = plt.subplots(1, len(MODELS), figsize=(5.2 * len(MODELS), 4.2), squeeze=False)
        for ax, m in zip(axes[0], MODELS):
            g = cell[cell["model"] == m]
            if g.empty:
                continue
            rho = spearmanr(g["predictive_uncertainty"], g[inst]).statistic
            rho_l = spearmanr(g["latent_uncertainty"], g[inst]).statistic
            ax.scatter(g["predictive_uncertainty"], g[inst], s=4, alpha=0.2,
                       color=COLOR[m], edgecolor="none")
            q = pd.qcut(g["predictive_uncertainty"], 20, duplicates="drop")
            binned = g.groupby(q, observed=True).agg(
                x=("predictive_uncertainty", "median"), y=(inst, "median"))
            ax.plot(binned["x"], binned["y"], color="black", lw=1.6, marker="o", ms=3)
            ax.set_xlabel("full-source 20-seed ensemble $U_{pred}$")
            ax.set_ylabel(f"{title}\n(true-class probability variance)")
            ax.set_title(f"{MODEL_LABEL[m]}   $\\rho_{{pred}}$={rho:.3f}, "
                         f"$\\rho_{{latent}}$={rho_l:.3f}", fontsize=8.5)
        fig.suptitle(f"Figure {num}. Full-model uncertainty versus {title} "
                     f"(black: 20-quantile medians)", fontsize=10.5, weight="bold")
        fig.tight_layout()
        save(fig, f"figure{num}_uncertainty_vs_"
                  f"{'source_subset' if num == 11 else 'model_seed'}_instability")


def fig13_risk_coverage(sel: pd.DataFrame) -> None:
    dirs = sorted(sel["direction"].unique())
    fig, axes = plt.subplots(2, len(dirs), figsize=(4.9 * len(dirs), 6.8), squeeze=False)
    for col, d in enumerate(dirs):
        for row, (metric, label) in enumerate((("accuracy", "accuracy"), ("macro_f1", "macro-F1"))):
            ax = axes[row][col]
            for m in MODELS:
                g = sel[(sel["direction"] == d) & (sel["model"] == m)]
                if g.empty:
                    continue
                stat = g.groupby("coverage")[metric].agg(["mean", "std"]).sort_index()
                ax.errorbar(stat.index, stat["mean"], yerr=stat["std"], marker="o", ms=4,
                            capsize=3, color=COLOR[m], label=MODEL_LABEL[m], lw=1.5)
            ax.set_xlabel("coverage (fraction retained)"); ax.set_ylabel(f"{label} of retained cells")
            ax.invert_xaxis()
            ax.set_title(f"{DIR_LABEL[d]}", fontsize=9)
            ax.legend(fontsize=7.5)
    fig.suptitle("Figure 13. Selective prediction across 20 seeds (mean $\\pm$ SD)",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure13_risk_coverage_20seed")


def fig14_retention(ret: pd.DataFrame) -> None:
    dirs = sorted(ret["direction"].unique())
    fig, axes = plt.subplots(len(dirs), len(MODELS),
                             figsize=(6.4 * len(MODELS), 4.6 * len(dirs)), squeeze=False)
    for i, d in enumerate(dirs):
        for j, m in enumerate(MODELS):
            ax = axes[i][j]
            g = ret[(ret["direction"] == d) & (ret["model"] == m)]
            if g.empty:
                continue
            piv = g.pivot_table(index="cell_type", columns="coverage", values="retention_rate")
            piv = piv.reindex(columns=sorted(piv.columns, reverse=True))
            piv = piv.sort_values(piv.columns[-1])
            im = ax.imshow(piv.to_numpy(), aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
            ax.set_xticks(range(len(piv.columns)), [f"{c:.0%}" for c in piv.columns], fontsize=7)
            ax.set_yticks(range(len(piv.index)), piv.index, fontsize=6.5)
            ax.set_xlabel("coverage"); ax.grid(False)
            ax.set_title(f"{DIR_LABEL[d]}  {MODEL_LABEL[m]}", fontsize=9)
            fig.colorbar(im, ax=ax, fraction=0.046, label="retention rate")
    fig.suptitle("Figure 14. Cell-type retention under uncertainty-based rejection",
                 fontsize=11, weight="bold")
    fig.tight_layout()
    save(fig, "figure14_celltype_retention")


def main() -> int:
    print("PHASE 11B figures", flush=True)
    fig1_design()
    cache: dict = {}
    primary = read("phase11b_20seed_primary.csv")
    if primary is not None:
        fig2_auroc_distribution(primary)
        fig4_5_latent(cache)
    delta = read("phase11b_seed_delta.csv")
    if delta is not None:
        fig3_paired_deltas(delta, read("phase11b_seed_delta_summary.csv"))
    conv = read("phase11b_ensemble_convergence.csv")
    if conv is not None:
        fig6_convergence(conv)
    fig7_8_heatmaps()
    vc = read("phase11b_variance_components.csv")
    if vc is not None:
        fig9_variance(vc)
    ct = read("phase11b_celltype_variance.csv")
    if ct is not None:
        fig10_celltype(ct)
    cell = read("phase11b_cell_instability.csv")
    if cell is not None:
        fig11_12_instability(cell)
    sel = read("phase11b_selective_prediction.csv")
    if sel is not None:
        fig13_risk_coverage(sel)
    ret = read("phase11b_celltype_retention.csv")
    if ret is not None:
        fig14_retention(ret)
    print("figures done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
