# app_dynamic.py - TRUE lazy loading version
from flask import Flask, request, send_file, render_template, jsonify
import os
import time
import json
from werkzeug.exceptions import BadRequest

# S3-aware dataset access
#from manifest_loader import MANIFEST, dataset_key_for, _fs
import fsspec
from manifest_loader import get_manifest, dataset_key_for
from storage import read_zarr # , S3_BUCKET, USE_S3

try:
    from storage import read_h5ad_backed   # optional fallback for .h5ad
    HAVE_H5AD_FALLBACK = True
except ImportError:
    HAVE_H5AD_FALLBACK = False

USE_S3   = os.getenv("USE_S3", "false").lower() == "true"
S3_BUCKET = os.getenv("S3_BUCKET", "")

def to_uri(key: str) -> str:
    if not USE_S3:
        return key
    return key if key.startswith("s3://") else f"s3://{S3_BUCKET}/{key.lstrip('/')}"

def _open_anndata_for(tissue: str):
    key = dataset_key_for(tissue)  # resolves via the (lazy) S3 manifest
    if key.endswith(".zarr"):
        return read_zarr(key)
    if HAVE_H5AD_FALLBACK:
        return read_h5ad_backed(key)
    raise RuntimeError(
        f"{tissue!r} maps to {key!r}, which isn’t a .zarr. Convert it to Zarr or add read_h5ad_backed in storage.py."
    )

# Initialize Flask app
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent  # <-- stay in gene-expression-backend
# Ensure manifest path resolves appropriately for local vs S3
if os.getenv("USE_S3", "false").lower() == "true":
	os.environ.setdefault("DATA_MANIFEST", "datasets_manifest.json")
else:
	os.environ.setdefault("DATA_MANIFEST", str(BASE_DIR / "datasets_manifest.json"))
app = Flask(
    __name__,
    template_folder=str(BASE_DIR / "templates"),
    static_folder=str(BASE_DIR / "static"),
)


print("🚀 TRUE Lazy Loading Flask app initialized - S3-ready via manifest_loader/storage")

@app.route("/health")
def health():
    return {"ok": True}, 200

@app.route("/healthz")
def healthz():
    return {"ok": True}, 200

@app.route("/")
def home():
    return render_template("index.html")

# @app.route("/api/tissues", methods=["GET"])
# def get_tissues():
#     """Return tissue names and availability from manifest; uses fsspec for size if possible."""
#     print("📊 API call: /api/tissues (manifest-based, no data loading)")
#     items = {}

#     fs = fsspec.filesystem("s3") if USE_S3 else fsspec.filesystem("file")
#     for tissue, entry in get_manifest().items():
#         key = entry.get("default")
#         rec = {"key": key}
#         try:
#             # Prefix bucket for relative keys in S3 mode
#             uri = key if (not USE_S3 or key.startswith("s3://")) else f"s3://{S3_BUCKET}/{key.lstrip('/')}"
#             # if USE_S3 and not key.startswith("s3://"):
#             #     uri = f"s3://{S3_BUCKET}/{key.lstrip('/')}"
#             info = fs.info(uri)
#             size = info.get("size")
#             if size is not None:
#                 rec["file_size_mb"] = round(size / (1024 * 1024), 1)
#             rec["available"] = True
#         except Exception as e:
#             rec["available"] = False
#             rec["error"] = str(e)
#         items[tissue] = rec
#     print(f"   ✅ {sum(1 for v in items.values() if v.get('available'))} available of {len(items)}")
#     return items

@app.route("/api/tissues", methods=["GET"])
def get_tissues():
    fs = fsspec.filesystem("s3") if USE_S3 else fsspec.filesystem("file")
    items = {}
    for tissue, entry in get_manifest().items():
        key = entry.get("default")
        uri = to_uri(key)
        rec = {"key": key}
        try:
            info = fs.info(uri)  # raises if missing/forbidden
            rec["available"] = True
            if "size" in info:
                rec["size"] = info["size"]
        except Exception as e:
            rec["available"] = False
            rec["error"] = str(e)
        items[tissue] = rec
    return items


_CELLTYPE_CANDIDATES = [
    "cell_type", "CellType", "celltype", "celltypes",
    "cluster", "annotation", "cell_type_major", "leiden", "louvain"
]

@app.route("/api/tissue-metadata", methods=["GET"])
def get_tissue_metadata():
    """Return lightweight metadata the UI needs to enable selectors."""
    tissue = request.args.get("tissue")
    if not tissue:
        raise BadRequest("Missing tissue")

    print(f"📊 /api/tissue-metadata for {tissue}")
    t0 = time.time()
    try:
        adata = _open_anndata_for(tissue)   # your S3-aware loader

        # genes (names only; stays lazy with Zarr)
        try:
            genes = adata.var_names.to_list()
        except Exception:
            genes = adata.var.index.astype(str).to_list()
        # keep payload sane (adjust if you like)
        if len(genes) > 50000:
            genes = genes[:50000]

        # optional cell types
        ct_col = next((c for c in _CELLTYPE_CANDIDATES if c in adata.obs.columns), None)
        cell_types = (
            sorted(adata.obs[ct_col].astype(str).unique().tolist()) if ct_col else []
        )

        # optional batches (your original fields preserved)
        if "batch" in adata.obs.columns:
            batches = sorted(adata.obs["batch"].astype(str).unique().tolist())
            counts = adata.obs["batch"].value_counts()
            batch_info = [{"name": b, "cells": int(counts[b])} for b in batches]
            samples = len(batches)
        else:
            batch_info = [{"name": tissue, "cells": int(adata.n_obs)}]
            samples = 1

        resp = {
            "tissue": tissue,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "genes": genes,              # <-- UI uses this to enable dropdowns
            "cell_types": cell_types,    # <-- optional; safe if empty
            "samples": samples,
            "batches": batch_info,
        }
        print(f"   ✅ metadata in {time.time()-t0:.2f}s "
              f"({resp['n_cells']:,} cells, {resp['n_genes']:,} genes)")
        return jsonify(resp)
    except Exception as e:
        print(f"   ❌ metadata error for {tissue}: {e}")
        return jsonify({"error": str(e)}), 500

@app.route("/api/genes", methods=["GET"])
def get_genes():
    """Return gene list from interacting partners - NO data loading"""
    print("🧬 API call: /api/genes (no data loading)")
    from plotting_backend import get_interacting_partners
    partners = get_interacting_partners()
    all_genes = set()
    for cytokine, receptors in partners.items():
        all_genes.add(cytokine)
        all_genes.update(receptors)
    print(f"   ✅ Returned {len(all_genes)} genes from interaction database")
    return {
        'genes': sorted(list(all_genes)),
        'count': len(all_genes)
    }

@app.route("/api/interacting-partners", methods=["GET"])
def get_interacting_partners_api():
    """Return interacting partners - NO data loading"""
    print("🔗 API call: /api/interacting-partners (no data loading)")
    from plotting_backend import get_interacting_partners
    partners = get_interacting_partners()
    return {
        'interacting_partners': partners,
        'cytokines': sorted(partners.keys()),
        'total_pairs': sum(len(receptors) for receptors in partners.values())
    }

@app.route("/api/cell-types", methods=["GET"])
def get_cell_types_for_gene():
    """Load cell types for a SPECIFIC gene+tissue combination when requested."""
    tissue = request.args.get("tissue")
    gene = request.args.get("gene")
    if not tissue or not gene:
        return {"error": "Missing tissue or gene parameter"}, 400

    print(f"🔬 API call: /api/cell-types for {gene} in {tissue} (loading specific data)")
    try:
        import numpy as np
        start_time = time.time()
        adata = _open_anndata_for(tissue)
        adata.obs['tissue_type'] = tissue
        # Filter by tissue if present
        tissue_adata = adata[adata.obs['tissue_type'] == tissue] if 'tissue_type' in adata.obs.columns else adata
        if gene not in tissue_adata.var.index:
            return {"error": f"Gene {gene} not found in {tissue}"}, 404
        gene_expr = tissue_adata[:, gene].X
        if hasattr(gene_expr, 'A'):
            expression_values = gene_expr.A.flatten()
        elif hasattr(gene_expr, 'toarray'):
            expression_values = gene_expr.toarray().flatten()
        else:
            expression_values = np.array(gene_expr).flatten()
        expressing_cells = expression_values > 0
        cell_types_column = "cell_types"
        if cell_types_column not in tissue_adata.obs.columns:
            return {"error": f"Cell types column not found in {tissue}"}, 404
        expressing_cell_types = tissue_adata.obs[cell_types_column][expressing_cells]
        if hasattr(expressing_cell_types, "cat"):
            unique_cell_types = [ct for ct in expressing_cell_types.cat.categories if ct in expressing_cell_types.values]
        else:
            unique_cell_types = sorted(expressing_cell_types.unique())
        unique_cell_types = [str(ct) for ct in unique_cell_types if str(ct) != 'nan']
        print(f"   ✅ Found {len(unique_cell_types)} cell types expressing {gene} ({time.time()-start_time:.2f}s)")
        return {
            "tissue": tissue,
            "gene": gene,
            "cell_types": unique_cell_types,
            "total_expressing_cells": int(np.sum(expressing_cells)),
            "total_cells": int(len(expression_values))
        }
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return {"error": str(e)}, 500

@app.route("/generate-plot", methods=["GET"])
def generate_plot():
    """Generate plot - load data via S3-aware storage and cache the PNG."""
    tissue = request.args.get("tissue")
    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")
    celltype1 = request.args.get("celltype1")
    celltype2 = request.args.get("celltype2")
    if not all([tissue, gene1, gene2]):
        return "Missing required parameters", 400

    print(f"🎨 PLOT REQUEST: {gene1} vs {gene2} in {tissue}")
    try:
        start_time = time.time()
        import matplotlib.pyplot as plt
        from plotting_backend import plot_histogram_for_pair_cached
        adata = _open_anndata_for(tissue)
        adata.obs['tissue_type'] = tissue
        print(f"   ✅ Data opened for plotting in {time.time()-start_time:.1f}s")
        is_cached, data = plot_histogram_for_pair_cached(
            adata, tissue, gene1, gene2, cell_type1=celltype1, cell_type2=celltype2
        )
        total_time = time.time() - start_time
        if is_cached:
            print(f"   ✅ Cached plot served in {total_time:.1f}s total")
            return send_file(str(data), mimetype="image/png")
        else:
            import io
            buf = io.BytesIO(data)
            buf.seek(0)
            print(f"   ✅ New plot generated and cached in {total_time:.1f}s total")
            return send_file(buf, mimetype="image/png")
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return f"Error: {str(e)}", 500

@app.route("/api/cache/stats", methods=["GET"])
def get_cache_stats():
    try:
        from cache_manager import cache_manager
        stats = cache_manager.get_cache_stats()
        return jsonify(stats)
    except Exception as e:
        return f"Error getting cache stats: {str(e)}", 500

@app.route("/api/cache/clear", methods=["POST"])
def clear_cache():
    try:
        from cache_manager import cache_manager
        cache_manager.clear_cache()
        return jsonify({"message": "Cache cleared successfully"})
    except Exception as e:
        return f"Error clearing cache: {str(e)}", 500

@app.route("/api/cache/cleanup", methods=["POST"])
def cleanup_cache():
    try:
        from cache_manager import cache_manager
        target_size = request.json.get("target_size_mb") if request.json else None
        cache_manager.cleanup_cache(target_size)
        return jsonify({"message": "Cache cleanup completed"})
    except Exception as e:
        return f"Error during cache cleanup: {str(e)}", 500

if __name__ == "__main__":
    print("🚀 Starting TRUE Lazy Loading Flask server with CACHING (S3-ready)...")
    print("   - Datasets resolved via manifest_loader")
    print("   - Backed/zarr access via storage.py (S3 or local)")
    app.run(debug=False, host="127.0.0.1", port=8000)