import os
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from typing import Tuple, Union
from cache_manager import histogram_cache_key, load_json, save_json, cache_path_for
import io

def load_adata(path):
    """Load and return the AnnData object from a given path"""
    return sc.read(path)


def get_interacting_partners():
    return {
        "Ccl1": ["Ccr8"],
        "Ccl2": ["Ccr2", "Ccr4"],
        "Ccl3": ["Ccr1", "Ccr5", "Ccr4"],
        "Ccl4": ["Ccr5", "Ccr1", "Ccr8"],
        "Ccl5": ["Ccr1", "Ccr4", "Ccr5"],
        "Ccl7": ["Ccr1", "Ccr2"],
        "Ccl8": ["Ccr2", "Ccr5", "Ccr8"],
        "Ccl12": ["Ccr2"],
        "Ccl17": ["Ccr4", "Ccr8"],
        "Ccl19": ["Ccr7"],
        "Ccl20": ["Ccr6"],
        "Ccl22": ["Ccr4"],
        "Ccl25": ["Ccr9"],
        "Cx3cl1": ["Cx3cr1"],
        "Cxcl9": ["Cxcr3"],
        "Cxcl10": ["Cxcr3"],
        "Cxcl12": ["Cxcr4"],
        "Cxcl13": ["Cxcr5"],
        "Cxcl16": ["Cxcr6"],
        "Xcl1": ["Xcr1"],
        "Ccr1": ["Ccl3", "Ccl4", "Ccl5", "Ccl7"],
        "Ccr2": ["Ccl2", "Ccl7", "Ccl8", "Ccl12"],
        "Ccr4": ["Ccl2", "Ccl3", "Ccl5", "Ccl17", "Ccl22"],
        "Ccr5": ["Ccl3", "Ccl4", "Ccl5", "Ccl8"],
        "Ccr6": ["Ccl20"],
        "Ccr7": ["Ccl19"],
        "Ccr8": ["Ccl1", "Ccl4", "Ccl8", "Ccl17"],
        "Ccr9": ["Ccl25"],
        "Cx3cr1": ["Cx3cl1"],
        "Cxcr3": ["Cxcl9", "Cxcl10"],
        "Cxcr4": ["Cxcl12"],
        "Cxcr5": ["Cxcl13"],
        "Cxcr6": ["Cxcl16"],
        "Xcr1": ["Xcl1"],
        "Il2": ["Il2ra", "Il2rb", "Il2rg"],
        "Il4": ["Il4ra", "Il2rg"],
        "Il6": ["Il6ra", "Il6st"],
        "Il7": ["Il7r", "Il2rg"],
        "Il10": ["Il10ra", "Il10rb"],
        "Il15": ["Il15ra", "Il2rb", "Il2rg"],
        "Il17a": ["Il17ra", "Il17rc"],
        "Il17f": ["Il17ra", "Il17rc"],
        "Il18": ["Il18r1", "Il18rap"],
        "Ifng": ["Ifngr1", "Ifngr2"],
        "Tgfb1": ["Tgfbr1", "Tgfbr2"],
        "Tgfb2": ["Tgfbr1", "Tgfbr2"],
        "Tgfb3": ["Tgfbr1", "Tgfbr2"],
        "Flt3l": ["Flt3"],
        "Il2ra": ["Il2"],
        "Il2rb": ["Il2", "Il15"],
        "Il2rg": ["Il2", "Il4", "Il7", "Il15"],
        "Il4ra": ["Il4"],
        "Il6ra": ["Il6"],
        "Il6st": ["Il6"],
        "Il7r": ["Il7"],
        "Il10ra": ["Il10"],
        "Il10rb": ["Il10"],
        "Il15ra": ["Il15"],
        "Il17ra": ["Il17a", "Il17f"],
        "Il17rc": ["Il17a", "Il17f"],
        "Il18r1": ["Il18"],
        "Il18rap": ["Il18"],
        "Ifngr1": ["Ifng"],
        "Ifngr2": ["Ifng"],
        "Tgfbr1": ["Tgfb1", "Tgfb2", "Tgfb3"],
        "Tgfbr2": ["Tgfb1", "Tgfb2", "Tgfb3"],
        "Flt3": ["Flt3l"]
    }


def plot_histogram_for_pair(adata, tissue, gene1, gene2):
    os.makedirs('figures/temp', exist_ok=True)

    # First filter by tissue type, then find all batches within that tissue
    tissue_adata = adata[adata.obs['tissue_type'] == tissue]
    
    if tissue_adata.n_obs == 0:
        print(f"❌ No data found for tissue type '{tissue}'")
        available_tissues = sorted(adata.obs['tissue_type'].unique())
        print(f"Available tissue types: {available_tissues}")
        # Create empty plot with error message
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.text(0.5, 0.5, f"No data found for tissue: {tissue}\nAvailable: {', '.join(available_tissues)}", 
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f"Error: Tissue '{tissue}' not found")
        return fig

    # Find all batches within this tissue type
    if 'batch' in tissue_adata.obs.columns:
        tissue_batches = sorted(tissue_adata.obs['batch'].unique())
    else:
        # If no batch column, treat the entire tissue as one batch
        tissue_batches = [tissue]
        tissue_adata.obs['batch'] = tissue
    
    n_samples = len(tissue_batches)
    print(f"🧬 Found {n_samples} {tissue} samples: {tissue_batches}")
    
    # Create figure with subplots: 3 rows (gene1, gene2, overlay) x n_samples columns
    fig, axes = plt.subplots(3, n_samples, figsize=(4*n_samples, 12), squeeze=False)
    
    # If only one sample, axes needs to be reshaped
    if n_samples == 1:
        axes = axes.reshape(3, 1)
    
    for col, batch_name in enumerate(tissue_batches):
        # Filter data for this specific batch within the tissue
        sub = tissue_adata[tissue_adata.obs['batch'] == batch_name]
        print(f"🧬 {batch_name}: {sub.n_obs} cells")
        
        if sub.n_obs == 0:
            for row in range(3):
                axes[row, col].text(0.5, 0.5, f"No cells in\n{batch_name}", 
                                   ha='center', va='center', transform=axes[row, col].transAxes)
                axes[row, col].set_title(f"Error: {batch_name}")
            continue
        
        # Get spatial coordinates
        x = sub.obsm['X_spatial'][:, 0]
        y = sub.obsm['X_spatial'][:, 1]
        
        # Prepare expression arrays for both genes (handle sparse)
        def _expr_for(g):
            gX = sub[:, g].X
            if hasattr(gX, 'A'):
                return gX.A.flatten()
            if hasattr(gX, 'toarray'):
                return gX.toarray().flatten()
            return np.array(gX).flatten()
        
        # Plot each gene (rows 0 and 1) using histogram heatmaps
        for row, gene in enumerate([gene1, gene2]):
            ax = axes[row, col]
            
            if gene not in sub.var.index:
                print(f"❌ Gene '{gene}' not found in {batch_name}")
                ax.text(0.5, 0.5, f"Gene '{gene}'\nnot found", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"{gene} in {batch_name}")
                continue
                
            # Handle sparse matrix properly
            w = _expr_for(gene)
            print(f"🧬 {gene} in {batch_name}: min={w.min():.3f}, max={w.max():.3f}, mean={w.mean():.3f}, non-zero={np.sum(w > 0)}")
            
            # Create 2D histogram
            im = ax.hist2d(x, y, weights=w, bins=100, cmap='inferno')
            ax.set_title(f"{gene} in {batch_name}")
            ax.set_xlabel("Spatial X")
            ax.set_ylabel("Spatial Y")
            
            # Add colorbar for each subplot
            plt.colorbar(im[3], ax=ax, label="Expression")
        
        # Row 2: overlay categorical scatter (gene1-only, gene2-only, both)
        ax_overlay = axes[2, col]
        # Black background for overlay row for better contrast
        ax_overlay.set_facecolor("#000000")
        if (gene1 in sub.var.index) and (gene2 in sub.var.index):
            w1 = _expr_for(gene1)
            w2 = _expr_for(gene2)
            expr1 = w1 > 0
            expr2 = w2 > 0
            both = expr1 & expr2
            only1 = expr1 & (~expr2)
            only2 = expr2 & (~expr1)

            # Plot order: gene2-only (purple), gene1-only (yellow), both (orange) on top
            if np.any(only2):
                ax_overlay.scatter(x[only2], y[only2], s=1, c="#800080", alpha=0.8, label=f"{gene2} only")
            if np.any(only1):
                ax_overlay.scatter(x[only1], y[only1], s=1, c="#FFD700", alpha=0.8, label=f"{gene1} only")
            if np.any(both):
                ax_overlay.scatter(x[both], y[both], s=1, c="#FFA500", alpha=0.95, label="both")
            ax_overlay.set_title(f"Overlay in {batch_name}", color='white')
            ax_overlay.set_xlabel("Spatial X", color='white')
            ax_overlay.set_ylabel("Spatial Y", color='white')
            ax_overlay.tick_params(axis='both', colors='white')
            # Compact legend with dark frame
            leg = ax_overlay.legend(loc='upper right', fontsize=8, frameon=True)
            if leg is not None:
                leg.get_frame().set_facecolor('#000000')
                leg.get_frame().set_edgecolor('white')
                for t in leg.get_texts():
                    t.set_color('white')
        else:
            ax_overlay.text(0.5, 0.5, "One or both genes not found", ha='center', va='center', transform=ax_overlay.transAxes, color='white')
            ax_overlay.set_title(f"Overlay in {batch_name}", color='white')

    # Adjust layout to prevent overlap
    plt.tight_layout()
    return fig


def plot_histogram_for_pair_cached(adata, tissue, gene1, gene2, cell_type1: str = None, cell_type2: str = None, bins: int = 100) -> Tuple[bool, Union[str, bytes]]:
    """
    Cached wrapper for plot_histogram_for_pair.
    Returns (is_cached, data) where data is either a file path (if cached) or PNG bytes (if newly generated).
    """
    key = histogram_cache_key(
        tissue=tissue,
        gene_x=gene1,
        gene_y=gene2,
        cell_type1=cell_type1 or "",
        cell_type2=cell_type2 or "",
        bins=bins,
        version="v4"
    )
    meta = load_json(key)
    if meta and os.path.exists(meta.get("png_path", "")):
        return True, meta["png_path"]

    fig = plot_histogram_for_pair(adata, tissue, gene1, gene2)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
    plt.close(fig)
    buf.seek(0)

    png_path = cache_path_for(key).replace(".json", ".png")
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    with open(png_path, "wb") as f:
        f.write(buf.getvalue())

    save_json(key, {"png_path": png_path, "tissue": tissue, "gene1": gene1, "gene2": gene2, "bins": bins})
    return False, buf.getvalue()
