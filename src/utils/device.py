"""Detect CUDA, Apple MPS, and CPU without assuming a GPU exists."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from typing import Any


def _nvidia_smi_path() -> str | None:
    found = shutil.which("nvidia-smi")
    if found:
        return found
    wsl = "/usr/lib/wsl/lib/nvidia-smi"
    if os.path.isfile(wsl) and os.access(wsl, os.X_OK):
        return wsl
    return None


def query_nvidia_smi() -> dict[str, Any]:
    """Best-effort NVIDIA driver / GPU query. Never raises."""
    info: dict[str, Any] = {
        "nvidia_smi_path": _nvidia_smi_path(),
        "driver_version": None,
        "cuda_compatibility": None,
        "gpu_name": None,
        "vram_total_mb": None,
        "vram_free_mb": None,
        "gpu_utilization_pct": None,
    }
    smi = info["nvidia_smi_path"]
    if not smi:
        return info
    try:
        proc = subprocess.run(
            [
                smi,
                "--query-gpu=name,driver_version,memory.total,memory.free,utilization.gpu",
                "--format=csv,noheader,nounits",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=15,
        )
        line = proc.stdout.strip().splitlines()[0]
        name, driver, total, free, util = [p.strip() for p in line.split(",")]
        info["gpu_name"] = name
        info["driver_version"] = driver
        info["vram_total_mb"] = float(total)
        info["vram_free_mb"] = float(free)
        info["gpu_utilization_pct"] = float(util)
    except Exception as exc:  # noqa: BLE001
        info["nvidia_smi_error"] = f"{type(exc).__name__}: {exc}"
    try:
        proc = subprocess.run(
            [smi],
            check=True,
            capture_output=True,
            text=True,
            timeout=15,
        )
        for raw in proc.stdout.splitlines():
            if "CUDA Version" in raw:
                info["cuda_compatibility"] = raw.split("CUDA Version:")[-1].split()[0]
                break
    except Exception:  # noqa: BLE001
        pass
    return info


def detect_compute_environment() -> dict[str, Any]:
    info: dict[str, Any] = {
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "hostname": platform.node(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cuda_available": False,
        "cuda_version": None,
        "cuda_device_count": 0,
        "cuda_device_index": 0,
        "gpu_names": [],
        "gpu_capability": None,
        "vram_total_mb": None,
        "mps_available": False,
        "accelerator": "cpu",
        "accelerator_detail": "CPU",
        "torch_version": None,
        "scvi_version": None,
        "scanpy_version": None,
        "anndata_version": None,
        "numpy_version": None,
        "pandas_version": None,
        "scipy_version": None,
        "sklearn_version": None,
        "matplotlib_version": None,
        "driver_version": None,
        "cuda_compatibility": None,
        "import_errors": {},
    }
    info.update({f"nvidia_{k}": v for k, v in query_nvidia_smi().items()})
    if info.get("nvidia_driver_version"):
        info["driver_version"] = info["nvidia_driver_version"]
    if info.get("nvidia_cuda_compatibility"):
        info["cuda_compatibility"] = info["nvidia_cuda_compatibility"]

    try:
        import numpy as np

        info["numpy_version"] = np.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["numpy"] = f"{type(exc).__name__}: {exc}"

    try:
        import pandas as pd

        info["pandas_version"] = pd.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["pandas"] = f"{type(exc).__name__}: {exc}"

    try:
        import scipy

        info["scipy_version"] = scipy.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["scipy"] = f"{type(exc).__name__}: {exc}"

    try:
        import sklearn

        info["sklearn_version"] = sklearn.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["sklearn"] = f"{type(exc).__name__}: {exc}"

    try:
        import matplotlib

        info["matplotlib_version"] = matplotlib.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["matplotlib"] = f"{type(exc).__name__}: {exc}"

    try:
        import anndata

        info["anndata_version"] = anndata.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["anndata"] = f"{type(exc).__name__}: {exc}"

    try:
        import scanpy as sc

        info["scanpy_version"] = sc.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["scanpy"] = f"{type(exc).__name__}: {exc}"

    try:
        import scvi

        info["scvi_version"] = scvi.__version__
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["scvi"] = f"{type(exc).__name__}: {exc}"

    try:
        import torch

        info["torch_version"] = torch.__version__
        info["cuda_available"] = bool(torch.cuda.is_available())
        info["cuda_version"] = torch.version.cuda
        if torch.cuda.is_available():
            info["cuda_device_count"] = int(torch.cuda.device_count())
            info["gpu_names"] = [
                torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())
            ]
            cap = torch.cuda.get_device_capability(0)
            info["gpu_capability"] = f"{cap[0]}.{cap[1]}"
            props = torch.cuda.get_device_properties(0)
            info["vram_total_mb"] = float(props.total_memory) / (1024.0 * 1024.0)
        mps_mod = getattr(torch.backends, "mps", None)
        info["mps_available"] = bool(mps_mod is not None and mps_mod.is_available())
    except Exception as exc:  # noqa: BLE001
        info["import_errors"]["torch"] = f"{type(exc).__name__}: {exc}"

    accelerator, detail = get_accelerator(
        cuda_available=info["cuda_available"],
        gpu_names=info["gpu_names"],
        mps_available=info["mps_available"],
    )
    info["accelerator"] = accelerator
    info["accelerator_detail"] = detail
    return info


def get_accelerator(
    cuda_available: bool | None = None,
    gpu_names: list[str] | None = None,
    mps_available: bool | None = None,
) -> tuple[str, str]:
    """Return (accelerator_name, human_readable_detail).

    Preference order: CUDA GPU, Apple MPS, CPU.
    Lightning/scvi-tools use accelerator='gpu' for NVIDIA CUDA.
    """
    if cuda_available is None or mps_available is None:
        env = detect_compute_environment()
        cuda_available = env["cuda_available"]
        gpu_names = env["gpu_names"]
        mps_available = env["mps_available"]

    if cuda_available:
        names = ", ".join(gpu_names or ["CUDA GPU"])
        return "gpu", names
    if mps_available:
        return "mps", "Apple Silicon MPS"
    return "cpu", "CPU"


def requested_accelerator(
    env: dict[str, Any],
    prefer_cuda: bool = True,
    prefer_mps: bool = True,
    allow_cpu_fallback: bool = True,
) -> tuple[str, str]:
    """Select the Lightning accelerator string for future training.

    Returns 'gpu' (CUDA), 'mps', or 'cpu'. CUDA is preferred on the RTX 4060.
    """
    del allow_cpu_fallback  # callers still pass it; fallback is handled in train loops
    if prefer_cuda and env.get("cuda_available"):
        names = ", ".join(env.get("gpu_names") or ["CUDA GPU"])
        detail = f"NVIDIA CUDA ({names})"
        if env.get("cuda_version"):
            detail += f"; torch.version.cuda={env['cuda_version']}"
        if env.get("driver_version"):
            detail += f"; driver={env['driver_version']}"
        return "gpu", detail
    if prefer_mps and env.get("mps_available"):
        return "mps", "Apple Silicon MPS (explicit; scvi auto would choose CPU)"
    return "cpu", "CPU"


def lightning_devices(cfg: dict[str, Any] | None, requested: str) -> int | str | list[int]:
    """Device argument for model.train(). CUDA uses the configured GPU index."""
    if requested != "gpu":
        return "auto"
    index = 0
    if cfg is not None:
        index = int(cfg.get("device", {}).get("preferred_device", 0))
    return [index]


def training_device_log(env: dict[str, Any], requested: str, actual: str) -> dict[str, Any]:
    return {
        "requested_device": requested,
        "actual_device": actual,
        "gpu_name": (env.get("gpu_names") or [None])[0],
        "torch_version": env.get("torch_version"),
        "cuda_version": env.get("cuda_version"),
        "driver_version": env.get("driver_version"),
        "cuda_compatibility": env.get("cuda_compatibility"),
        "cuda_available": env.get("cuda_available"),
        "hostname": env.get("hostname"),
    }
