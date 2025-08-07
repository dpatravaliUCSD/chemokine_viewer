import os
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

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
    
    # Create figure with subplots: 2 rows (gene1, gene2) x n_samples columns
    fig, axes = plt.subplots(2, n_samples, figsize=(4*n_samples, 8), squeeze=False)
    
    # If only one sample, axes needs to be reshaped
    if n_samples == 1:
        axes = axes.reshape(2, 1)
    
    for col, batch_name in enumerate(tissue_batches):
        # Filter data for this specific batch within the tissue
        sub = tissue_adata[tissue_adata.obs['batch'] == batch_name]
        print(f"🧬 {batch_name}: {sub.n_obs} cells")
        
        if sub.n_obs == 0:
            for row in range(2):
                axes[row, col].text(0.5, 0.5, f"No cells in\n{batch_name}", 
                                   ha='center', va='center', transform=axes[row, col].transAxes)
                axes[row, col].set_title(f"Error: {batch_name}")
            continue
        
        # Get spatial coordinates
        x = sub.obsm['X_spatial'][:, 0]
        y = sub.obsm['X_spatial'][:, 1]
        
        # Plot each gene
        for row, gene in enumerate([gene1, gene2]):
            ax = axes[row, col]
            
            if gene not in sub.var.index:
                print(f"❌ Gene '{gene}' not found in {batch_name}")
                ax.text(0.5, 0.5, f"Gene '{gene}'\nnot found", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"{gene} in {batch_name}")
                continue
                
            # Handle sparse matrix properly
            gene_expr = sub[:, gene].X
            if hasattr(gene_expr, 'A'):
                w = gene_expr.A.flatten()
            elif hasattr(gene_expr, 'toarray'):
                w = gene_expr.toarray().flatten()
            else:
                w = np.array(gene_expr).flatten()
                
            print(f"🧬 {gene} in {batch_name}: min={w.min():.3f}, max={w.max():.3f}, mean={w.mean():.3f}, non-zero={np.sum(w > 0)}")
            
            # Create 2D histogram
            im = ax.hist2d(x, y, weights=w, bins=100, cmap='viridis')
            ax.set_title(f"{gene} in {batch_name}")
            ax.set_xlabel("Spatial X")
            ax.set_ylabel("Spatial Y")
            
            # Add colorbar for each subplot
            plt.colorbar(im[3], ax=ax, label="Expression")
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    return fig
