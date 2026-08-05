from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a JSON configuration file."""
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_from_root(root: str | Path, value: str | Path) -> Path:
    """Resolve a repository-relative path from *root*."""
    path = Path(value)
    return path if path.is_absolute() else Path(root) / path
