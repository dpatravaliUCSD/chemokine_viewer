# app_minimal.py - Minimal version for fast startup
from flask import Flask, request, send_file, render_template, jsonify
import os
import time

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

print("🚀 Minimal Flask app initialized - NO data loading")

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/api/tissues", methods=["GET"])
def get_tissues():
    """Return tissue info without loading data"""
    print("📊 API call: /api/tissues")
    
    # Return basic info without loading actual data
    tissue_info = {
        "Liver": {
            "samples": 6,
            "total_cells": 241303,  # Approximate from previous run
            "batches": [
                {"name": "Liver_1_comp", "cells": 26276},
                {"name": "Liver_2_comp", "cells": 40547},
                {"name": "Liver_3", "cells": 75336},
                {"name": "Liver_3_comp", "cells": 27256},
                {"name": "Liver_4_comp", "cells": 20395},
                {"name": "Liver_5_comp", "cells": 51493}
            ]
        },
        "SI": {
            "samples": 2,
            "total_cells": 536102,  # Approximate from previous run
            "batches": [
                {"name": "SI_1", "cells": 268051},
                {"name": "SI_2", "cells": 268051}
            ]
        },
        "MC38 tumor": {
            "samples": 6,
            "total_cells": 577678,  # Approximate from previous run
            "batches": [
                {"name": "Tumor_1", "cells": 96280},
                {"name": "Tumor_2", "cells": 96280},
                {"name": "Tumor_3", "cells": 96280},
                {"name": "Tumor_4", "cells": 96280},
                {"name": "Tumor_5", "cells": 96280},
                {"name": "Tumor_6", "cells": 96278}
            ]
        }
    }
    
    return tissue_info

@app.route("/api/genes", methods=["GET"])
def get_genes():
    """Return hardcoded gene list for now"""
    print("🧬 API call: /api/genes")
    
    # Return all genes from the interacting partners dictionary
    from plotting_backend import get_interacting_partners
    partners = get_interacting_partners()
    
    all_genes = set()
    for cytokine, receptors in partners.items():
        all_genes.add(cytokine)
        all_genes.update(receptors)
    
    return {
        'genes': sorted(list(all_genes)),
        'count': len(all_genes)
    }

@app.route("/api/interacting-partners", methods=["GET"])
def get_interacting_partners_api():
    """Return interacting partners data"""
    print("🔗 API call: /api/interacting-partners")
    
    from plotting_backend import get_interacting_partners
    partners = get_interacting_partners()
    
    return {
        'interacting_partners': partners,
        'cytokines': sorted(partners.keys()),
        'total_pairs': sum(len(receptors) for receptors in partners.values())
    }

@app.route("/generate-plot", methods=["GET"])
def generate_plot():
    """Generate plot - this is where we actually load data"""
    tissue = request.args.get("tissue")
    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")

    if not all([tissue, gene1, gene2]):
        return "Missing required parameters", 400

    print(f"🎨 PLOT REQUEST: {gene1} vs {gene2} in {tissue}")
    print("   📁 NOW loading data (this may take time)...")
    
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
    print("🚀 Starting MINIMAL Flask server...")
    print("   - No data loaded at startup")
    print("   - Data loaded only when generating plots")
    app.run(debug=False, host="127.0.0.1", port=8002) 