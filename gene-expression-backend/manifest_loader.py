# manifest_loader.py
"""
Lazy loader for the datasets manifest, supporting both local files and S3.

Env:
  USE_S3=true|false
  S3_BUCKET=<your-bucket>
  DATA_MANIFEST=<key or path>
    - If USE_S3=true: either "s3://bucket/prefix/datasets_manifest.json" or
      a bucket-relative key like "datasets_manifest.json" or "config/manifest.json".
    - If USE_S3=false: a local path like "datasets_manifest.json".
"""

from __future__ import annotations
import json
import os
from functools import lru_cache
from typing import Dict, Any

import fsspec

USE_S3 = os.getenv("USE_S3", "false").lower() == "true"
S3_BUCKET = os.getenv("S3_BUCKET", "")
DATA_MANIFEST = os.getenv("DATA_MANIFEST", "datasets_manifest.json")


def _fs():
    return fsspec.filesystem("s3") if USE_S3 else fsspec.filesystem("file")


def _resolve_manifest_path() -> str:
    """
    Returns a path (local) or URI (s3://...) for the manifest.
    """
    dm = DATA_MANIFEST
    if USE_S3:
        if dm.startswith("s3://"):
            return dm
        if not S3_BUCKET:
            raise RuntimeError("S3 mode is enabled but S3_BUCKET is not set.")
        return f"s3://{S3_BUCKET}/{dm.lstrip('/')}"
    # local
    return dm


def load_manifest() -> Dict[str, Any]:
    """
    Loads the manifest JSON from the configured filesystem.
    """
    path = _resolve_manifest_path()
    fs = _fs()
    with fs.open(path, "rb") as f:
        return json.load(f)


# Lazy, in-process cache of the manifest (invalidate via reload_manifest()).
_MANIFEST: Dict[str, Any] | None = None


def get_manifest() -> Dict[str, Any]:
    global _MANIFEST
    if _MANIFEST is None:
        _MANIFEST = load_manifest()
    return _MANIFEST


def reload_manifest() -> Dict[str, Any]:
    global _MANIFEST
    _MANIFEST = None
    return get_manifest()


def dataset_key_for(tissue: str) -> str:
    """
    Returns the store key/path for a given tissue (e.g., "xenium_liver.zarr").
    """
    m = get_manifest()
    if tissue not in m:
        raise KeyError(f"No manifest entry for tissue '{tissue}'.")
    entry = m[tissue]
    key = entry.get("default")
    if not key:
        raise KeyError(f"Manifest entry for '{tissue}' has no 'default' key.")
    return key
