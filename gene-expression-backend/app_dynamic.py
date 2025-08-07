# app_hybrid.py - Semi-dynamic version with smart caching
from flask import Flask, request, send_file, render_template, jsonify
import os
import time
import json
import pickle

# Define paths for all tissue types
tissue_paths = {
    "Liver": "/mnt/sata4/Alex_Xenium_Data/20250717__210317__2025_07_17_perturb8_LCMV_run1/Liver/combined/final_object.h5ad",
    "SI": "/mnt/sata4/Alex_Xenium_Data/20250717__210317__2025_07_17_perturb8_LCMV_run1/SI/combined/guides_assigned.h5ad", 
    "MC38 tumor": "/mnt/sata1/Xenium_data/20250703__195730__perturb8_tumor_run1/combined/final_object.h5ad"
}

# Cache files
METADATA_CACHE_FILE = "tissue_metadata_cache.json"
GENE_CACHE_FILE = "gene_cache.json"

def get_file_modification_time(filepath):
    """Get file modification time for cache invalidation"""
    try:
        return os.path.getmtime(filepath)
    except:
        return 0

def load_tissue_metadata_cached():
    """Load tissue metadata with smart caching"""
    cache_valid = True
    
    # Check if cache exists and is newer than source files
    if os.path.exists(METADATA_CACHE_FILE):
        cache_time = os.path.getmtime(METADATA_CACHE_FILE)
        for filepath in tissue_paths.values():
            if os.path.exists(filepath) and os.path.getmtime(filepath) > cache_time:
                cache_valid = False
                break
    else:
        cache_valid = False
    
    # Load from cache if valid
    if cache_valid:
        try:
            with open(METADATA_CACHE_FILE, 'r') as f:
                print("📋 Loading tissue metadata from cache")
                return json.load(f)
        except:
            cache_valid = False
    
    # Generate fresh metadata
    print("🔄 Generating fresh tissue metadata...")
    import scanpy as sc
    
    tissue_info = {}
    for tissue_name, filepath in tissue_paths.items():
        if not os.path.exists(filepath):
            continue
            
        print(f"   📊 Scanning {tissue_name}...")
        start_time = time.time()
        
        # Load just the obs (metadata) - much faster than full data
        adata = sc.read(filepath, first_column_names=True)
        
        if 'batch' in adata.obs.columns:
            batches = sorted(adata.obs['batch'].unique())
            batch_counts = adata.obs['batch'].value_counts()
            tissue_info[tissue_name] = {
                'samples': len(batches),
                'total_cells': adata.n_obs,
                'batches': [{'name': batch, 'cells': int(batch_counts[batch])} for batch in batches]
            }
        else:
            tissue_info[tissue_name] = {
                'samples': 1,
                'total_cells': adata.n_obs,
                'batches': [{'name': tissue_name, 'cells': adata.n_obs}]
            }
        
        scan_time = time.time() - start_time
        print(f"   ✅ {tissue_name}: {adata.n_obs:,} cells ({scan_time:.1f}s)")
    
    # Save to cache
    try:
        with open(METADATA_CACHE_FILE, 'w') as f:
            json.dump(tissue_info, f, indent=2)
        print("💾 Metadata cached for future use")
    except Exception as e:
        print(f"⚠️ Could not save cache: {e}")
    
    return tissue_info

def load_genes_cached():
    """Load available genes with caching"""
    cache_valid = True
    
    # Check cache validity
    if os.path.exists(GENE_CACHE_FILE):
        cache_time = os.path.getmtime(GENE_CACHE_FILE)
        # Check if any data files are newer
        for filepath in tissue_paths.values():
            if os.path.exists(filepath) and os.path.getmtime(filepath) > cache_time:
                cache_valid = False
                break
    else:
        cache_valid = False
    
    # Load from cache if valid
    if cache_valid:
        try:
            with open(GENE_CACHE_FILE, 'r') as f:
                print("🧬 Loading genes from cache")
                return json.load(f)
        except:
            cache_valid = False
    
    # Generate fresh gene list
    print("🔄 Generating fresh gene list...")
    from plotting_backend import get_interacting_partners
    import scanpy as sc
    
    # Get genes from interacting partners
    partners = get_interacting_partners()
    interacting_genes = set()
    for cytokine, receptors in partners.items():
        interacting_genes.add(cytokine)
        interacting_genes.update(receptors)
    
    # Verify genes exist in at least one dataset
    print("   🔍 Verifying gene availability...")
    sample_tissue = next(iter(tissue_paths.keys()))
    sample_path = tissue_paths[sample_tissue]
    
    if os.path.exists(sample_path):
        adata = sc.read(sample_path, first_column_names=True)
        available_genes = set(adata.var.index)
        
        # Filter to only genes that exist
        verified_genes = interacting_genes.intersection(available_genes)
        print(f"   ✅ {len(verified_genes)}/{len(interacting_genes)} genes available")
    else:
        verified_genes = interacting_genes
        print(f"   ⚠️ Using all {len(interacting_genes)} genes (could not verify)")
    
    gene_data = {
        'genes': sorted(list(verified_genes)),
        'count': len(verified_genes),
        'interacting_partners': {k: v for k, v in partners.items() if k in verified_genes}
    }
    
    # Save to cache
    try:
        with open(GENE_CACHE_FILE, 'w') as f:
            json.dump(gene_data, f, indent=2)
        print("💾 Genes cached for future use")
    except Exception as e:
        print(f"⚠️ Could not save gene cache: {e}")
    
    return gene_data

# Initialize Flask app
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
template_folder = os.path.join(base_dir, "templates")
static_folder = os.path.join(base_dir, "static")

app = Flask(
    __name__,
    template_folder=template_folder,
    static_folder=static_folder
)

print("🚀 Hybrid Flask app initialized - Smart caching enabled")

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/api/tissues", methods=["GET"])
def get_tissues():
    """Return tissue info with smart caching"""
    print("📊 API call: /api/tissues")
    try:
        return load_tissue_metadata_cached()
    except Exception as e:
        print(f"❌ Error loading tissues: {e}")
        return {'error': str(e)}, 500

@app.route("/api/genes", methods=["GET"])
def get_genes():
    """Return gene list with smart caching"""
    print("🧬 API call: /api/genes")
    try:
        gene_data = load_genes_cached()
        return {
            'genes': gene_data['genes'],
            'count': gene_data['count']
        }
    except Exception as e:
        print(f"❌ Error loading genes: {e}")
        return {'error': str(e)}, 500

@app.route("/api/interacting-partners", methods=["GET"])
def get_interacting_partners_api():
    """Return interacting partners data with smart caching"""
    print("🔗 API call: /api/interacting-partners")
    try:
        gene_data = load_genes_cached()
        return {
            'interacting_partners': gene_data['interacting_partners'],
            'cytokines': sorted(gene_data['interacting_partners'].keys()),
            'total_pairs': sum(len(receptors) for receptors in gene_data['interacting_partners'].values())
        }
    except Exception as e:
        print(f"❌ Error loading interacting partners: {e}")
        return {'error': str(e)}, 500

@app.route("/api/refresh-cache", methods=["POST"])
def refresh_cache():
    """Manual cache refresh endpoint"""
    try:
        # Delete cache files to force refresh
        for cache_file in [METADATA_CACHE_FILE, GENE_CACHE_FILE]:
            if os.path.exists(cache_file):
                os.remove(cache_file)
        
        print("🔄 Cache cleared - will regenerate on next API call")
        return {'status': 'Cache refreshed successfully'}
    except Exception as e:
        return {'error': str(e)}, 500

@app.route("/generate-plot", methods=["GET"])
def generate_plot():
    """Generate plot - loads data on demand"""
    tissue = request.args.get("tissue")
    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")

    if not all([tissue, gene1, gene2]):
        return "Missing required parameters", 400

    print(f"🎨 PLOT REQUEST: {gene1} vs {gene2} in {tissue}")
    print("   📁 Loading tissue data...")
    
    try:
        start_time = time.time()
        
        # Import heavy modules only when needed
        import scanpy as sc
        import matplotlib.pyplot as plt
        from plotting_backend import plot_histogram_for_pair
        import io
        
        # Load the specific tissue data
        file_path = tissue_paths.get(tissue)
        if not file_path or not os.path.exists(file_path):
            return f"Tissue data not found: {tissue}", 404
        
        print(f"   📁 Loading {tissue} from {file_path}")
        adata = sc.read(file_path)
        adata.obs['tissue_type'] = tissue
        
        load_time = time.time() - start_time
        print(f"   ✅ Data loaded in {load_time:.1f}s")
        
        # Generate plot
        fig = plot_histogram_for_pair(adata, tissue, gene1, gene2)

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
    print("🚀 Starting HYBRID Flask server...")
    print("   - Smart caching for metadata and genes")
    print("   - Auto-refresh when data files change")
    print("   - Manual refresh available at /api/refresh-cache")
    app.run(debug=False, host="127.0.0.1", port=8003) 