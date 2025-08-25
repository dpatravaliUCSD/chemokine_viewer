# cache_manager.py
"""
Lightweight file-based cache helpers for plotting results (e.g., histogram bins).
Not strictly required, but useful to avoid recomputing heavy plots.
"""

from __future__ import annotations
import hashlib
import json
import os
import pathlib
from typing import Any, Dict, Optional

CACHE_DIR = os.getenv("CACHE_DIR", "/tmp/hist-cache")
pathlib.Path(CACHE_DIR).mkdir(parents=True, exist_ok=True)


def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-._" else "_" for c in s)


def histogram_cache_key(
    tissue: str,
    gene_x: str,
    gene_y: str,
    cell_type1: str = "",
    cell_type2: str = "",
    bins: int = 200,
    version: str = "v1"
) -> str:
    payload = f"{version}|{tissue}|{gene_x}|{gene_y}|{cell_type1}|{cell_type2}|{bins}"
    h = hashlib.sha256(payload.encode()).hexdigest()[:24]
    return f"hist_{_safe_name(tissue)}_{_safe_name(gene_x)}_{_safe_name(gene_y)}_{h}.json"


def cache_path_for(key: str) -> str:
    return os.path.join(CACHE_DIR, key)


def load_json(key: str) -> Optional[Dict[str, Any]]:
    p = cache_path_for(key)
    if os.path.exists(p):
        try:
            with open(p, "r") as f:
                return json.load(f)
        except Exception:
            return None
    return None


def save_json(key: str, data: Dict[str, Any]) -> None:
    p = cache_path_for(key)
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f)
    os.replace(tmp, p)
