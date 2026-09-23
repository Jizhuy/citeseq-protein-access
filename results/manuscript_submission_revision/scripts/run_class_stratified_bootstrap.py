#!/usr/bin/env python3
"""Class-stratified paired bootstrap from cached held-out predictions."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "tables"
FIG = Path(__file__).resolve().parents[1] / "figures"
CACHE = OUT / "development_heldout_predictions.npz"
PRIMARY = ROOT / "results/presubmission_revision/bootstrap/core_contrast_bootstrap_summary.csv"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

N_BOOT = 10_000
STRAT_SEED = 20260924


def macro_f1(y_true, y_pred) -> float:
    return float(f1_score(y_true, y_pred, average="macro", zero_division=0))


def stratified_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    parts = []
    for lab in np.unique(y):
        idx = np.flatnonzero(y == lab)
        parts.append(rng.choice(idx, size=len(idx), replace=True))
    return np.concatenate(parts)


def summarize(df: pd.DataFrame, point: dict, method: str) -> pd.DataFrame:
    rows = []
    for key, pe in point.items():
        s = df[key]
        rows.append(
            {
                "method": method,
                "quantity": key,
                "point_estimate": pe,
                "boot_mean": float(s.mean()),
                "boot_median": float(s.median()),
                "boot_sd": float(s.std(ddof=1)),
                "ci_2.5": float(np.nanpercentile(s, 2.5)),
                "ci_97.5": float(np.nanpercentile(s, 97.5)),
                "prop_gt_0": float(np.mean(s > 0))
                if key in ("G_access", "G_assoc", "Gap_total")
                else np.nan,
                "n_boot": int(len(s)),
                "n_nan": int(s.isna().sum()),
            }
        )
    return pd.DataFrame(rows)


def main():
    if not CACHE.exists():
        raise SystemExit(f"Missing cache {CACHE}; run cache_development_predictions.py first")
    z = np.load(CACHE, allow_pickle=True)
    y = z["y"]
    p_rna = z["pred_rna"]
    p_cat = z["pred_concat"]
    p_tot = z["pred_totalvi"]
    # Prefer frozen manuscript-aligned point estimates from primary summary when available
    primary = pd.read_csv(PRIMARY)
    pe_map = {r.quantity: float(r.point_estimate) for _, r in primary.iterrows()}
    point = {
        "F1_RNA_PCA": pe_map.get("F1_RNA_PCA", macro_f1(y, p_rna)),
        "F1_concat": pe_map.get("F1_concat", macro_f1(y, p_cat)),
        "F1_totalVI": pe_map.get("F1_totalVI", macro_f1(y, p_tot)),
        "G_access": pe_map.get("G_access"),
        "G_assoc": pe_map.get("G_assoc"),
        "Gap_total": pe_map.get("Gap_total"),
        "R_access": pe_map.get("R_access"),
    }
    if point["G_access"] is None:
        point["G_access"] = point["F1_concat"] - point["F1_RNA_PCA"]
        point["G_assoc"] = point["F1_totalVI"] - point["F1_concat"]
        point["Gap_total"] = point["F1_totalVI"] - point["F1_RNA_PCA"]
        point["R_access"] = point["G_access"] / point["Gap_total"]

    rng = np.random.default_rng(STRAT_SEED)
    rows = []
    for b in range(N_BOOT):
        idx = stratified_indices(y, rng)
        yt = y[idx]
        fr = macro_f1(yt, p_rna[idx])
        fc = macro_f1(yt, p_cat[idx])
        ft = macro_f1(yt, p_tot[idx])
        g_a = fc - fr
        g_m = ft - fc
        gap = ft - fr
        rows.append(
            {
                "replicate": b,
                "F1_RNA_PCA": fr,
                "F1_concat": fc,
                "F1_totalVI": ft,
                "G_access": g_a,
                "G_assoc": g_m,
                "Gap_total": gap,
                "R_access": g_a / gap if gap != 0 else np.nan,
            }
        )
    strat_df = pd.DataFrame(rows)
    strat_sum = summarize(strat_df, point, "class_stratified_paired")
    strat_sum.to_csv(OUT / "bootstrap_class_stratified_summary.csv", index=False)

    primary_out = primary[primary.quantity.isin(point.keys())].copy()
    primary_out["method"] = "global_paired"
    keep = [
        "method",
        "quantity",
        "point_estimate",
        "boot_mean",
        "boot_median",
        "boot_sd",
        "ci_2.5",
        "ci_97.5",
        "prop_gt_0",
        "n_boot",
        "n_nan",
    ]
    if "n_nan" not in primary_out.columns:
        primary_out["n_nan"] = 0
    compare = pd.concat([primary_out[keep], strat_sum[keep]], ignore_index=True)
    compare.to_csv(OUT / "bootstrap_global_vs_class_stratified.csv", index=False)

    meta = {
        "n_boot": N_BOOT,
        "strat_seed": STRAT_SEED,
        "design": "within-class with-replacement resampling of held-out cells; paired predictions",
        "interval": "percentile 2.5/97.5",
        "cache": str(CACHE.relative_to(ROOT)),
        "note": "Sensitivity only; primary remains global paired bootstrap.",
    }
    (OUT / "bootstrap_class_stratified_meta.json").write_text(json.dumps(meta, indent=2))

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for ax, q in zip(axes, ["G_access", "G_assoc", "R_access"]):
        sub = compare[compare.quantity == q].reset_index(drop=True)
        for i, row in sub.iterrows():
            ax.errorbar(
                i,
                row.point_estimate,
                yerr=[[row.point_estimate - row["ci_2.5"]], [row["ci_97.5"] - row.point_estimate]],
                fmt="o",
                capsize=4,
            )
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub.method.tolist(), rotation=15, ha="right", fontsize=7)
        ax.set_title(q)
        ax.axhline(0, color="gray", lw=0.8, ls="--")
    fig.suptitle("Global vs class-stratified paired bootstrap (n=10,000; conditional)")
    fig.tight_layout()
    fig.savefig(FIG / "bootstrap_global_vs_class_stratified.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(strat_sum.to_string(index=False))


if __name__ == "__main__":
    main()
