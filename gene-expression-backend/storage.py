# storage.py
"""
S3-aware data openers for AnnData. Prefer Zarr for cloud-native reading.
Optional .h5ad fallback downloads to a local cache and opens backed='r'.
"""

from __future__ import annotations
import hashlib
import os
import pathlib
from typing import Optional

import anndata as ad
import fsspec
from s3fs import S3FileSystem

# Env flags
USE_S3 = os.getenv("USE_S3", "false").lower() == "true"
S3_BUCKET = os.getenv("S3_BUCKET", "")
AWS_REGION = os.getenv("AWS_DEFAULT_REGION", "us-east-2")  # set to your bucket region

# Optional cache dir for .h5ad fallback
H5AD_CACHE_DIR = os.getenv("H5AD_CACHE_DIR", "/tmp/h5ad-cache")
pathlib.Path(H5AD_CACHE_DIR).mkdir(parents=True, exist_ok=True)


def _to_s3_uri(path_or_key: str) -> str:
    if path_or_key.startswith("s3://"):
        return path_or_key
    if not S3_BUCKET:
        raise RuntimeError("S3 mode requested but S3_BUCKET is empty.")
    return f"s3://{S3_BUCKET}/{path_or_key.lstrip('/')}"


def _s3_mapper(uri: str):
    """
    Create an fsspec mapper with tuned S3 options for Zarr.
    """
    storage_opts = {
        "anon": False,
        "client_kwargs": {"region_name": AWS_REGION},
        "config_kwargs": {"max_pool_connections": 64},
    }
    return fsspec.get_mapper(uri, **storage_opts)


def read_zarr(path_or_key: str) -> ad.AnnData:
    """
    Read an AnnData Zarr store from S3 (or local if USE_S3=false).
    """
    if USE_S3:
        uri = path_or_key if path_or_key.startswith("s3://") else _to_s3_uri(path_or_key)
        mapper = _s3_mapper(uri)
        return ad.read_zarr(mapper)
    # Local filesystem path
    return ad.read_zarr(path_or_key)


def read_h5ad_backed(path_or_key: str) -> ad.AnnData:
    """
    Near-term fallback for .h5ad: download to a local cache once, then open backed='r'.
    Useful until all datasets are converted to Zarr.

    Requires: s3fs; enough local disk for the file; beware ephemeral storage on Render.
    """
    if not USE_S3:
        return ad.read_h5ad(path_or_key, backed="r")

    uri = _to_s3_uri(path_or_key)
    h = hashlib.sha256(uri.encode()).hexdigest()[:20]
    local = os.path.join(H5AD_CACHE_DIR, f"{h}.h5ad")

    if not os.path.exists(local):
        fs = S3FileSystem(anon=False)
        # fs.get takes a bucket/key string, not an s3:// URI
        fs.get(uri.replace("s3://", ""), local)

    return ad.read_h5ad(local, backed="r")
