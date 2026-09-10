#!/usr/bin/env python
"""PHASE 12A §3 integrity audit for replicate-1 grid + PHASE 11B intactness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.experiments.phase11b_robustness import (  # noqa: E402
    CROSSED_MODEL_SEEDS,
    CROSSED_SUBSET_SEEDS,
    MODELS,
    experiment_tag,
    phase11b_paths,
)
from src.experiments.phase12a_replication import (  # noqa: E402
    CROSSED_DIRECTION_KEY,
    REPLICATES,
    crossed_jobs,
    ensure_dirs,
    run_replicate_seed,
)


def main() -> int:
    paths = ensure_dirs()
    root = paths["root"]
    p11 = phase11b_paths()
    issues: list[str] = []
    ok: list[str] = []

    # A. 50 replicate-1 done markers
    rep1_tags = [
        experiment_tag(j["direction_key"], j["model"], j["model_seed"], j["subset_seed"], 1)
        for j in crossed_jobs(replicates=(1,))
    ]
    done_files = sorted(paths["done"].glob("*.json"))
    done_tags = [p.stem for p in done_files]
    print(f"A. done_markers={len(done_tags)} expected=50")
    if len(done_tags) != 50:
        issues.append(f"A: expected 50 done markers, found {len(done_tags)}")
    missing = sorted(set(rep1_tags) - set(done_tags))
    extra = sorted(set(done_tags) - set(rep1_tags))
    if missing:
        issues.append(f"A: missing done: {missing}")
    if extra:
        issues.append(f"A: unexpected done: {extra}")
    for p in done_files:
        marker = json.loads(p.read_text(encoding="utf-8"))
        if not marker.get("complete"):
            issues.append(f"A: marker not complete: {p.name}")
    if not missing and not extra and len(done_tags) == 50:
        ok.append("A: exactly 50 valid replicate-1 done-markers")

    # B. PHASE 11B crossed 50 intact
    p11_crossed = [
        experiment_tag(CROSSED_DIRECTION_KEY, model, k, s, None)
        for s in CROSSED_SUBSET_SEEDS
        for k in CROSSED_MODEL_SEEDS
        for model in MODELS
    ]
    p11_done = {p.stem for p in p11["done"].glob("*.json")}
    missing11 = [t for t in p11_crossed if t not in p11_done]
    print(f"B. phase11b crossed present={50 - len(missing11)}/50; total_p11b_done={len(p11_done)}")
    if missing11:
        issues.append(f"B: missing PHASE 11B crossed: {missing11}")
    else:
        ok.append("B: all 50 PHASE 11B crossed done-markers intact")

    # C. design size
    n_obs = (
        len(CROSSED_SUBSET_SEEDS)
        * len(CROSSED_MODEL_SEEDS)
        * len(REPLICATES)
        * len(MODELS)
    )
    print(f"C. design observations={n_obs} (expect 100)")
    if n_obs != 100:
        issues.append(f"C: design size {n_obs} != 100")
    else:
        ok.append("C: 5 x 5 x 2 x 2 = 100 observations")

    # D–I. per-run artifacts for replicate 1
    required_emb = [
        "_source_latent_post.npy",
        "_target_latent_post.npy",
        "_source_latent_pre.npy",
        "_source_cells.csv",
        "_target_cells.csv",
    ]
    n_ok = 0
    reload_ok = 0
    finite_ok = 0
    coord_ok = 0
    seed_ok = 0
    sample_keys: list[str] = []

    for j in crossed_jobs(replicates=(1,)):
        tag = experiment_tag(
            j["direction_key"], j["model"], j["model_seed"], j["subset_seed"], 1
        )
        model = j["model"]
        s, k = j["subset_seed"], j["model_seed"]
        mdir = (
            paths["models"]
            / CROSSED_DIRECTION_KEY
            / model
            / f"subset{s}_seed{k}"
            / "replicate01"
        )
        if not mdir.exists():
            mdir = (
                paths["models"]
                / CROSSED_DIRECTION_KEY
                / model
                / f"subset{s}_seed{k:02d}"
                / "replicate01"
            )
        man_path = paths["logs"] / "manifests" / f"{tag}_manifest.json"
        post_path = paths["posterior"] / CROSSED_DIRECTION_KEY / f"{tag}_posterior.npz"
        emb_dir = paths["embeddings"] / CROSSED_DIRECTION_KEY
        problems: list[str] = []

        manifest = None
        if man_path.exists():
            manifest = json.loads(man_path.read_text(encoding="utf-8"))
            if not sample_keys:
                sample_keys = sorted(manifest.keys())
            integ = manifest.get("artifact_integrity") or {}
            max_diff = integ.get("max_latent_diff")
            if max_diff is None:
                problems.append("missing artifact_integrity.max_latent_diff")
            elif float(max_diff) > 1e-4 or integ.get("passed") is False:
                problems.append(f"reload integrity fail max_latent_diff={max_diff}")
            else:
                reload_ok += 1

            expected_run = run_replicate_seed(k, s, 1)
            got_run = manifest.get("run_replicate_seed")
            got_model = manifest.get("model_seed")
            if got_model is not None and int(got_model) != int(k):
                problems.append(f"model_seed mismatch manifest={got_model} expected={k}")
            if got_run is not None and int(got_run) != int(expected_run):
                problems.append(
                    f"run_replicate_seed mismatch manifest={got_run} expected={expected_run}"
                )
            if got_model is not None and got_run is not None and int(got_run) != int(k):
                seed_ok += 1
            if manifest.get("replicate") not in (1, "1"):
                problems.append(f"replicate field={manifest.get('replicate')}")
            if manifest.get("source_embedding_used_for_metrics") != "post_adaptation":
                problems.append(
                    f"source_embedding_used_for_metrics="
                    f"{manifest.get('source_embedding_used_for_metrics')}"
                )

            src_dir = Path(manifest["source_model_dir"])
            q_dir = Path(manifest["adapted_query_model_dir"])
            src_model = src_dir / "model.pt"
            q_model = q_dir / "model.pt"
            if not src_model.exists():
                problems.append("missing source_model/model.pt")
            if not q_model.exists():
                problems.append("missing adapted_query_model/model.pt")
            if not manifest.get("adapted_query_model_saved"):
                problems.append("adapted_query_model_saved is not True")
        else:
            problems.append("missing manifest")
            src_model = mdir / "source_model" / "model.pt"
            q_model = mdir / "adapted_query_model" / "model.pt"
            if not src_model.exists():
                problems.append("missing source_model")
            if not q_model.exists():
                problems.append("missing adapted_query_model")

        if not post_path.exists():
            problems.append("missing posterior")
        for suf in required_emb:
            f = emb_dir / f"{tag}{suf}"
            if not f.exists() or f.stat().st_size == 0:
                problems.append(f"missing/empty {suf}")

        # finite + same latent dim (post-adaptation coordinate system)
        try:
            z_s = np.load(emb_dir / f"{tag}_source_latent_post.npy")
            z_t = np.load(emb_dir / f"{tag}_target_latent_post.npy")
            if not np.isfinite(z_s).all():
                problems.append("NaN/Inf in source latent post")
            if not np.isfinite(z_t).all():
                problems.append("NaN/Inf in target latent post")
            if z_s.ndim != 2 or z_t.ndim != 2 or z_s.shape[1] != z_t.shape[1]:
                problems.append(f"latent shape mismatch src{z_s.shape} tgt{z_t.shape}")
            else:
                coord_ok += 1
            post = np.load(post_path, allow_pickle=True)
            for key in ("target_mu", "target_var", "source_var", "u_latent_target"):
                if key not in post.files:
                    problems.append(f"posterior missing key {key}")
                    continue
                arr = post[key]
                if not np.issubdtype(arr.dtype, np.number) or not np.isfinite(arr).all():
                    problems.append(f"NaN/Inf or non-numeric posterior[{key}]")
            finite_ok += 1
            # target probabilities written at scoring time
            prob_file = None
            for cand in [
                paths["probabilities"] / f"rep_{tag}_probs.npz",
                paths["probabilities"] / CROSSED_DIRECTION_KEY / f"{tag}_target_probs.npy",
                paths["probabilities"] / f"{tag}_probs.npz",
            ]:
                if cand.exists():
                    prob_file = cand
                    break
            if prob_file is None:
                problems.append("missing target probabilities")
            else:
                if prob_file.suffix == ".npz":
                    pr_npz = np.load(prob_file, allow_pickle=True)
                    for key in pr_npz.files:
                        arr = pr_npz[key]
                        if np.issubdtype(arr.dtype, np.number) and arr.size and not np.isfinite(arr).all():
                            problems.append(f"NaN/Inf in probabilities[{key}]")
                else:
                    pr = np.load(prob_file)
                    if not np.isfinite(pr).all():
                        problems.append("NaN/Inf in target probabilities")
        except Exception as exc:  # noqa: BLE001
            problems.append(f"load error: {exc}")

        if problems:
            issues.append(f"D:{tag}: " + "; ".join(problems))
        else:
            n_ok += 1

    print(f"D. complete artifact runs={n_ok}/50")
    print(f"G. reload integrity pass={reload_ok}/50")
    print(f"H. finite latent/posterior loads={finite_ok}/50")
    print(f"I. matching post-adapt latent dims={coord_ok}/50")
    print(f"seed bookkeeping (run_seed!=model_seed where recorded)={seed_ok}/50")
    if n_ok != 50:
        issues.append(f"D: only {n_ok}/50 runs fully intact")
    else:
        ok.append("D: all 50 runs have models/manifests/latents/posterior")
    if reload_ok != 50:
        issues.append(f"G: reload validation passed for {reload_ok}/50")
    else:
        ok.append("G: adapted-model reload validation passed for all 50")
    if finite_ok != 50:
        issues.append(f"H: finite checks passed for {finite_ok}/50")
    else:
        ok.append("H: no NaN/Inf in required latent/posterior outputs")
    if coord_ok != 50:
        issues.append(f"I: coordinate-system dim match for {coord_ok}/50")
    else:
        ok.append("I: source/target post-adaptation latent dims match")

    # E failures
    fail_dir = paths["logs"] / "failures"
    fails = list(fail_dir.glob("*")) if fail_dir.exists() else []
    # ignore empty dirs
    fails = [p for p in fails if p.is_file()]
    print(f"E. failure files={len(fails)}")
    if fails:
        issues.append(f"E: failure files present: {[p.name for p in fails]}")
    else:
        ok.append("E: no failed runs remain")

    # F incomplete dirs without done markers
    orphan = []
    models_root = paths["models"] / CROSSED_DIRECTION_KEY
    for d in models_root.rglob("replicate01"):
        try:
            seed_part = d.parent.name
            model = d.parent.parent.name
            subset = int(seed_part.split("_")[0].replace("subset", ""))
            seed = int(seed_part.split("_")[1].replace("seed", ""))
            tag = experiment_tag(CROSSED_DIRECTION_KEY, model, seed, subset, 1)
            if not (paths["done"] / f"{tag}.json").exists():
                orphan.append(str(d))
        except Exception:
            orphan.append(str(d))
    print(f"F. orphan replicate01 dirs without done-marker={len(orphan)}")
    if orphan:
        issues.append(f"F: dirs without done-marker: {orphan}")
    else:
        ok.append("F: no incomplete orphan replicate directories")

    # J PHASE 11B unchanged at marker level
    if len(p11_done) < 130:
        issues.append(f"J: PHASE 11B done count dropped to {len(p11_done)} (<130)")
    else:
        ok.append(f"J: PHASE 11B done-markers intact (n={len(p11_done)})")

    # §4 replicate validity across a few manifests
    print("\n=== REPLICATE VALIDITY (§4) ===")
    print("sample_manifest_keys:", sample_keys[:30])
    validity_ok = True
    for s in (0, 1):
        for k in (0, 1):
            for model in MODELS:
                tag1 = experiment_tag(CROSSED_DIRECTION_KEY, model, k, s, 1)
                man1 = json.loads(
                    (paths["logs"] / "manifests" / f"{tag1}_manifest.json").read_text()
                )
                # PHASE 11B / replicate 0
                tag0 = experiment_tag(CROSSED_DIRECTION_KEY, model, k, s, None)
                man0_path = p11["logs"] / "manifests" / f"{tag0}_manifest.json"
                if not man0_path.exists():
                    # alternate naming
                    alt = list((p11["logs"] / "manifests").glob(f"{tag0}*"))
                    man0_path = alt[0] if alt else man0_path
                if not man0_path.exists():
                    print(f"WARN no phase11b manifest for {tag0}")
                    validity_ok = False
                    continue
                man0 = json.loads(man0_path.read_text(encoding="utf-8"))
                rs1 = man1.get("run_replicate_seed")
                rs0 = man0.get("run_replicate_seed")  # may be None
                ms1 = man1.get("model_seed")
                ms0 = man0.get("model_seed")
                print(
                    f"{model} subset{s} seed{k}: "
                    f"model_seed {ms0}->{ms1}; run_seed {rs0}->{rs1}; "
                    f"expected_run1={run_replicate_seed(k,s,1)}"
                )
                if ms1 is not None and int(ms1) != k:
                    validity_ok = False
                    issues.append(f"seed: model_seed not fixed for {tag1}")
                if rs1 is None or int(rs1) == k:
                    # run seed should exist and differ from model seed identity
                    if rs1 is None:
                        validity_ok = False
                        issues.append(f"seed: missing run_replicate_seed for {tag1}")
                if rs0 not in (None, "None") and rs1 is not None and int(rs0) == int(rs1):
                    validity_ok = False
                    issues.append(f"seed: run_replicate_seed identical across replicates for {tag1}")

    if validity_ok:
        ok.append("§4: model_seed fixed; run_replicate_seed differs for replicate 1")

    report = {
        "pass": len(issues) == 0,
        "n_done": len(done_tags),
        "n_artifact_ok": n_ok,
        "n_reload_ok": reload_ok,
        "n_finite_ok": finite_ok,
        "n_coord_ok": coord_ok,
        "n_failures": len(fails),
        "n_issues": len(issues),
        "issues": issues,
        "ok": ok,
        "sample_manifest_keys": sample_keys,
    }
    out = paths["logs"] / "integrity_audit.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("\n=== ISSUES ===")
    if issues:
        for item in issues:
            print("!", item)
    else:
        print("NONE — integrity PASS")
    print("\n=== OK ===")
    for item in ok:
        print("+", item)
    print(f"wrote {out}")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
