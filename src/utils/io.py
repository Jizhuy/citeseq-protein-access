"""YAML/JSON helpers and path resolution relative to the project root."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml


def project_root() -> Path:
    """Return the repository root (`multiomics_robustness/`)."""
    return Path(__file__).resolve().parents[2]


def resolve_path(path_like: str | Path) -> Path:
    path = Path(path_like)
    if path.is_absolute():
        return path
    return project_root() / path


def load_yaml(path: str | Path) -> dict[str, Any]:
    with resolve_path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected a mapping in {path}, got {type(data)}")
    return data


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    config_path = path or (project_root() / "config" / "config.yaml")
    return load_yaml(config_path)


def save_json(obj: Any, path: str | Path) -> Path:
    out = resolve_path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, default=_json_default)
        handle.write("\n")
    return out


def save_yaml(obj: Any, path: str | Path) -> Path:
    out = resolve_path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(obj, handle, sort_keys=False)
    return out


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item"):
        return obj.item()
    return str(obj)
