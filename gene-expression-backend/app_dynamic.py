# app_dynamic.py - TRUE lazy loading version
from flask import Flask, request, send_file, render_template, jsonify
import os
import time
import json

# Define paths for all tissue types
tissue_paths = {
    "Liver": "/mnt/sata4/Alex_Xenium_Data/20250717__210317__2025_07_17_perturb8_LCMV_run1/Liver/combined/final_object.h5ad",
    "SI": "/mnt/sata4/Alex_Xenium_Data/20250717__210317__2025_07_17_perturb8_LCMV_run1/SI/combined/guides_assigned.h5ad", 
    "MC38 tumor": "/mnt/sata1/Xenium_data/20250703__195730__perturb8_tumor_run1/combined/final_object.h5ad"
}

# Initialize Flask app
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
template_folder = os.path.join(base_dir, "templates")
static_folder = os.path.join(base_dir, "static")

app = Flask(
    __name__,
    template_folder=template_folder,
    static_folder=static_folder
)

print("🚀 TRUE Lazy Loading Flask app initialized - NO data loading until requested")

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/api/tissues", methods=["GET"])
def get_tissues():
    """Return only tissue names and basic info - NO data loading"""
    print("📊 API call: /api/tissues (no data loading)")
    
    # Return only what exists, without loading any data files
    available_tissues = {}
    for tissue_name, filepath in tissue_paths.items():
        if os.path.exists(filepath):
            # Get file size as a rough indicator without loading
            file_size = os.path.getsize(filepath)
            available_tissues[tissue_name] = {
                "available": True,
                "file_size_mb": round(file_size / (1024*1024), 1),
                "note": "Metadata will be loaded when selected"
            }
        else:
            available_tissues[tissue_name] = {
                "available": False,
                "error": "File not found"
            }
    
    print(f"   ✅ Found {len([t for t in available_tissues.values() if t['available']])} available tissues")
    return available_tissues

@app.route("/api/tissue-metadata", methods=["GET"])
def get_tissue_metadata():
    """Load metadata for a SPECIFIC tissue when selected"""
    tissue = request.args.get("tissue")
    if not tissue:
        return {"error": "Missing tissue parameter"}, 400
    
    print(f"📊 API call: /api/tissue-metadata for {tissue} (loading metadata only)")
    
    filepath = tissue_paths.get(tissue)
    if not filepath or not os.path.exists(filepath):
        return {"error": f"Tissue {tissue} not found"}, 404
    
    try:
        import scanpy as sc
        start_time = time.time()
        
        # Load ONLY metadata (obs) - much faster than full data
        print(f"   📁 Loading metadata for {tissue}...")
        adata = sc.read(filepath, first_column_names=True)
        
        # Extract metadata without loading expression data
        if 'batch' in adata.obs.columns:
            batches = sorted(adata.obs['batch'].unique())
            batch_counts = adata.obs['batch'].value_counts()
            tissue_info = {
                'samples': len(batches),
                'total_cells': adata.n_obs,
                'batches': [{'name': batch, 'cells': int(batch_counts[batch])} for batch in batches]
            }
        else:
            tissue_info = {
                'samples': 1,
                'total_cells': adata.n_obs,
                'batches': [{'name': tissue, 'cells': adata.n_obs}]
            }
        
        load_time = time.time() - start_time
        print(f"   ✅ {tissue} metadata loaded: {adata.n_obs:,} cells ({load_time:.2f}s)")
        
        return tissue_info
        
    except Exception as e:
        print(f"   ❌ Error loading {tissue} metadata: {e}")
        return {"error": str(e)}, 500

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
    """Load cell types for a SPECIFIC gene+tissue combination when requested"""
    tissue = request.args.get("tissue")
    gene = request.args.get("gene")
    
    if not tissue or not gene:
        return {"error": "Missing tissue or gene parameter"}, 400
    
    print(f"🔬 API call: /api/cell-types for {gene} in {tissue} (loading specific data)")
    
    try:
        import scanpy as sc
        import numpy as np
        start_time = time.time()
        
        # Load the specific tissue data
        file_path = tissue_paths.get(tissue)
        if not file_path or not os.path.exists(file_path):
            return {"error": f"Tissue data not found: {tissue}"}, 404
        
        print(f"   📁 Loading {tissue} data to check {gene} expression...")
        adata = sc.read(file_path)
        adata.obs['tissue_type'] = tissue
        
        # Filter by tissue
        tissue_adata = adata[adata.obs['tissue_type'] == tissue]
        
        if gene not in tissue_adata.var.index:
            return {"error": f"Gene {gene} not found in {tissue}"}, 404
        
        # Get expression data for this gene
        gene_expr = tissue_adata[:, gene].X
        if hasattr(gene_expr, 'A'):
            expression_values = gene_expr.A.flatten()
        elif hasattr(gene_expr, 'toarray'):
            expression_values = gene_expr.toarray().flatten()
        else:
            expression_values = np.array(gene_expr).flatten()
        
        # Find cells that express this gene (expression > 0)
        expressing_cells = expression_values > 0
        
        # Get cell types for expressing cells
        cell_types_column = "cell_types"  # Default column name
        if cell_types_column not in tissue_adata.obs.columns:
            return {"error": f"Cell types column not found in {tissue}"}, 404
        
        expressing_cell_types = tissue_adata.obs[cell_types_column][expressing_cells]
        
        # Get unique cell types that express this gene
        if hasattr(expressing_cell_types, "cat"):
            unique_cell_types = list(expressing_cell_types.cat.categories)
            # Filter to only those that actually appear in expressing cells
            unique_cell_types = [ct for ct in unique_cell_types if ct in expressing_cell_types.values]
        else:
            unique_cell_types = sorted(expressing_cell_types.unique())
        
        # Convert to string and remove any NaN values
        unique_cell_types = [str(ct) for ct in unique_cell_types if str(ct) != 'nan']
        
        load_time = time.time() - start_time
        print(f"   ✅ Found {len(unique_cell_types)} cell types expressing {gene} ({load_time:.2f}s)")
        
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
    """Generate plot - this is where we load full data for visualization"""
    tissue = request.args.get("tissue")
    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")
    celltype1 = request.args.get("celltype1")  # Optional cell type filter for gene1
    celltype2 = request.args.get("celltype2")  # Optional cell type filter for gene2

    if not all([tissue, gene1, gene2]):
        return "Missing required parameters", 400

    print(f"🎨 PLOT REQUEST: {gene1} vs {gene2} in {tissue}")
    if celltype1:
        print(f"   🔬 Cell type filter for {gene1}: {celltype1}")
    if celltype2:
        print(f"   🔬 Cell type filter for {gene2}: {celltype2}")
    print("   📁 NOW loading full data for plotting...")
    
    try:
        start_time = time.time()
        
        # Import heavy modules only when needed for plotting
        import scanpy as sc
        import matplotlib.pyplot as plt
        from plotting_backend import plot_histogram_for_pair
        import io
        
        # Load the specific tissue data for plotting
        file_path = tissue_paths.get(tissue)
        if not file_path or not os.path.exists(file_path):
            return f"Tissue data not found: {tissue}", 404
        
        print(f"   📁 Loading {tissue} for plotting...")
        adata = sc.read(file_path)
        adata.obs['tissue_type'] = tissue
        
        load_time = time.time() - start_time
        print(f"   ✅ Data loaded for plotting in {load_time:.1f}s")
        
        # Generate plot with separate cell type filtering for each gene
        fig = plot_histogram_for_pair(adata, tissue, gene1, gene2, cell_type1=celltype1, cell_type2=celltype2)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0)
        plt.close(fig)
        
        total_time = time.time() - start_time
        print(f"   ✅ Plot generated in {total_time:.1f}s total")

        return send_file(buf, mimetype="image/png")

    except Exception as e:
        print(f"   ❌ Error: {e}")
        return f"Error: {str(e)}", 500

if __name__ == "__main__":
    print("🚀 Starting TRUE Lazy Loading Flask server...")
    print("   - No data loaded at startup")
    print("   - Tissue metadata loaded only when tissue selected")
    print("   - Cell types loaded only when gene+tissue selected")
    print("   - Full data loaded only when generating plots")
    app.run(debug=False, host="127.0.0.1", port=8003) 