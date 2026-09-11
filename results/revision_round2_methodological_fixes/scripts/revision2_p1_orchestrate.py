#!/usr/bin/env python3
"""Orchestrate P1 controls with --status / --resume."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
SCRIPTS = [
    "revision2_simple_multimodal_baseline.py",
    "revision2_rna_only_annotation.py",
    "revision2_testtime_corruption.py",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    for name in SCRIPTS:
        cmd = [PY, str(HERE / name)]
        if args.status:
            cmd.append("--status")
        elif args.dry_run:
            cmd.append("--dry-run")
        elif args.resume:
            cmd.append("--resume")
        print("==>", " ".join(cmd))
        subprocess.run(cmd, check=False)


if __name__ == "__main__":
    main()
