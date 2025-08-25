# data_access.py
"""
Small facade around manifest + storage. Import from your routes:
    from data_access import load_dataset_for, genes_for, cell_types_for
"""

from __future__ import annotations
from functools import lru_cache
from typing import List

from manifest_loader import dataset_key_for
from storage import read_zarr

# Optional fallback; keep if you still have any .h5ad keys in the manifest.
try:
    from storage import read_h5ad_backed
    HAVE_H5AD = True
except Exception:
    HAVE_H5AD = False


def load_dataset_for(tissue: str):
    key = dataset_key_for(tissue)
    if key.endswith(".zarr"):
        return read_zarr(key)
    if HAVE_H5AD:
        return read_h5ad_backed(key)
    raise RuntimeError(
        f"{tissue!r} resolves to {key!r} (not .zarr) and .h5ad fallback is disabled."
    )


@lru_cache(maxsize=32)
def genes_for(tissue: str) -> List[str]:
    A = load_dataset_for(tissue)
    try:
        return A.var_names.to_list()
    except Exception:
        return A.var.index.astype(str).to_list()


@lru_cache(maxsize=32)
def cell_types_for(tissue: str) -> List[str]:
    A = load_dataset_for(tissue)
    for cand in ["cell_type", "CellType", "celltype", "annotation", "cluster", "leiden", "louvain"]:
        if cand in A.obs.columns:
            s = A.obs[cand]
            if getattr(s.dtype, "name", "").startswith("category"):
                return [str(x) for x in s.cat.categories.to_list()]
            # non-categorical: sample to avoid heavy S3 scans
            n = min(10000, s.shape[0])
            return sorted(set(map(str, s.sample(n, random_state=0).tolist())))
    return []
