"""Runtime, process RSS, and CUDA VRAM tracking.

Process RSS and GPU VRAM are recorded separately. Do not treat one as the other.
"""

from __future__ import annotations

import resource
import sys
import time
from typing import Any


def peak_rss_mb() -> float:
    """Peak resident set size of this process.

    Darwin reports ru_maxrss in bytes; Linux reports kilobytes. This is CPU
    process RSS, not GPU memory.
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return float(rss) / (1024.0 * 1024.0)
    return float(rss) / 1024.0


def reset_cuda_peak_stats() -> None:
    """Reset CUDA peak memory counters before a training run."""
    try:
        import torch

        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.synchronize()
    except Exception:  # noqa: BLE001
        return


def cuda_memory_stats() -> dict[str, Any]:
    """Current and peak CUDA VRAM. Empty/None fields when CUDA is unavailable."""
    stats: dict[str, Any] = {
        "cuda_available": False,
        "cuda_allocated_mb": None,
        "cuda_reserved_mb": None,
        "peak_cuda_allocated_mb": None,
        "peak_cuda_reserved_mb": None,
    }
    try:
        import torch

        if not torch.cuda.is_available():
            return stats
        stats["cuda_available"] = True
        stats["cuda_allocated_mb"] = float(torch.cuda.memory_allocated()) / (1024.0 * 1024.0)
        stats["cuda_reserved_mb"] = float(torch.cuda.memory_reserved()) / (1024.0 * 1024.0)
        stats["peak_cuda_allocated_mb"] = float(torch.cuda.max_memory_allocated()) / (
            1024.0 * 1024.0
        )
        stats["peak_cuda_reserved_mb"] = float(torch.cuda.max_memory_reserved()) / (
            1024.0 * 1024.0
        )
    except Exception as exc:  # noqa: BLE001
        stats["cuda_memory_error"] = f"{type(exc).__name__}: {exc}"
    return stats


class Timer:
    def __init__(self) -> None:
        self.started: float | None = None
        self.elapsed: float | None = None

    def start(self) -> None:
        self.started = time.perf_counter()

    def stop(self) -> float:
        if self.started is None:
            raise RuntimeError("Timer was never started.")
        self.elapsed = time.perf_counter() - self.started
        return self.elapsed


def compute_resource_metrics(runtime_seconds: float | None = None) -> dict[str, Any]:
    metrics = {
        "runtime_seconds": runtime_seconds,
        "peak_rss_mb": peak_rss_mb(),
        "peak_memory_source": "resource.ru_maxrss (process RSS; not GPU VRAM)",
    }
    metrics.update({f"gpu_{k}": v for k, v in cuda_memory_stats().items()})
    return metrics
