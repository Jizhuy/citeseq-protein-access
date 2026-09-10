#!/usr/bin/env python
"""PHASE 12B sequential gate orchestrator with chunked process restarts.

Runs gate1 → gate2 → gate3 → remaining seeds in chunks of CHUNK jobs per
fresh Python process (avoids WSL RSS creep). Resumes via done-markers.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("TMPDIR", str(Path.home() / "tmp"))
os.environ["PATH"] = "/usr/lib/wsl/lib:" + os.environ.get("PATH", "")
if "LD_LIBRARY_PATH" in os.environ:
    os.environ["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + os.environ["LD_LIBRARY_PATH"]
else:
    os.environ["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib"

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
DONE = ROOT / "results/phase12b_external_validation/formal_validation/done"
LOG = ROOT / "results/phase12b_external_validation/formal_validation/logs"
CHUNK = 8  # fresh interpreter every N remaining runs


def done_ok(tag: str) -> bool:
    p = DONE / f"{tag}.json"
    if not p.exists():
        return False
    try:
        return bool(json.loads(p.read_text()).get("complete"))
    except Exception:
        return False


def validate_gate(tags: list[str], name: str) -> None:
    missing = [t for t in tags if not done_ok(t)]
    if missing:
        raise RuntimeError(f"{name} incomplete: missing {missing}")
    for t in tags:
        prim = json.loads((DONE / f"{t}.json").read_text())["primary"]
        for key in ("macro_f1", "predictive_error_auroc", "latent_error_auroc"):
            v = float(prim.get(key))
            if not (0.0 <= v <= 1.0):
                raise RuntimeError(f"{t}: invalid {key}={v}")
        if not prim.get("integrity_ok", False):
            raise RuntimeError(f"{t}: integrity_ok is false")
        if prim.get("n_proteins") != 12:
            raise RuntimeError(f"{t}: n_proteins != 12")
        if prim.get("n_genes") != 12776:
            raise RuntimeError(f"{t}: unexpected n_genes={prim.get('n_genes')}")
    print(f"VALIDATE OK: {name}", flush=True)


def run_cmd(args: list[str]) -> None:
    print("RUN", " ".join(args), flush=True)
    subprocess.check_call(args, cwd=str(ROOT))


def pending_jobs() -> list[tuple[int, str, int]]:
    jobs = []
    for fold in range(5):
        for model in ("scvi_matched", "totalvi"):
            for seed in range(10):
                tag = f"fold{fold}__{model}__seed{seed:02d}"
                if not done_ok(tag):
                    jobs.append((fold, model, seed))
    return jobs


def run_job_chunk(jobs: list[tuple[int, str, int]]) -> None:
    """Run a small list of jobs in a fresh interpreter (memory hygiene)."""
    for fold, model, seed in jobs:
        run_cmd(
            [
                sys.executable,
                "scripts/12b_run_gate.py",
                "--fold",
                str(fold),
                "--model",
                model,
                "--seed",
                str(seed),
            ]
        )


def main() -> None:
    LOG.mkdir(parents=True, exist_ok=True)
    Path(os.environ["TMPDIR"]).mkdir(parents=True, exist_ok=True)

    # Gates (skip via done-markers inside gate runner)
    run_cmd([sys.executable, "scripts/12b_run_gate.py", "--gate", "1"])
    validate_gate(["fold0__scvi_matched__seed00", "fold0__totalvi__seed00"], "GATE1")

    run_cmd([sys.executable, "scripts/12b_run_gate.py", "--gate", "2"])
    validate_gate(["fold1__scvi_matched__seed00", "fold1__totalvi__seed00"], "GATE2")

    run_cmd([sys.executable, "scripts/12b_run_gate.py", "--gate", "3"])
    tags3 = [f"fold{f}__{m}__seed00" for f in range(5) for m in ("scvi_matched", "totalvi")]
    validate_gate(tags3, "GATE3")

    # Remaining seeds in chunks with process isolation
    while True:
        pending = pending_jobs()
        print(f"PENDING={len(pending)}/100", flush=True)
        if not pending:
            break
        chunk = pending[:CHUNK]
        print(f"CHUNK={chunk}", flush=True)
        run_job_chunk(chunk)

    tags_all = [
        f"fold{f}__{m}__seed{s:02d}"
        for f in range(5)
        for m in ("scvi_matched", "totalvi")
        for s in range(10)
    ]
    validate_gate(tags_all, "ALL_100")

    run_cmd([sys.executable, "scripts/12b_analyze.py"])
    print("PHASE12B ORCHESTRATOR COMPLETE", flush=True)


if __name__ == "__main__":
    main()
