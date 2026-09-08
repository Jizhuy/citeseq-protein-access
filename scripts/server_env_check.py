#!/usr/bin/env python
"""Report the compute environment for the RTX 4060 server.

Does not train models. Safe to run any time.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.utils.device import detect_compute_environment, query_nvidia_smi  # noqa: E402
from src.evaluation.resource_metrics import cuda_memory_stats, peak_rss_mb  # noqa: E402


def main() -> None:
    env = detect_compute_environment()
    smi = query_nvidia_smi()
    gpu = cuda_memory_stats()
    report = {
        "hostname": env.get("hostname"),
        "cwd": os.getcwd(),
        "project_root": str(ROOT),
        "python_version": env.get("python"),
        "python_executable": env.get("python_executable"),
        "scvi_tools_version": env.get("scvi_version"),
        "torch_version": env.get("torch_version"),
        "torch_cuda_version": env.get("cuda_version"),
        "cuda_available": env.get("cuda_available"),
        "gpu_names": env.get("gpu_names"),
        "gpu_capability": env.get("gpu_capability"),
        "vram_total_mb": env.get("vram_total_mb") or smi.get("vram_total_mb"),
        "vram_free_mb": smi.get("vram_free_mb"),
        "driver_version": env.get("driver_version") or smi.get("driver_version"),
        "cuda_compatibility": env.get("cuda_compatibility") or smi.get("cuda_compatibility"),
        "accelerator": env.get("accelerator"),
        "accelerator_detail": env.get("accelerator_detail"),
        "peak_rss_mb": peak_rss_mb(),
        "cuda_memory": gpu,
        "import_errors": env.get("import_errors"),
    }
    print(json.dumps(report, indent=2, default=str))
    if not env.get("cuda_available"):
        raise SystemExit("CUDA is not available. Fix PyTorch/NVIDIA before training.")


if __name__ == "__main__":
    main()
