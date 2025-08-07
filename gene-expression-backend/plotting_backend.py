import os

def load_adata(path):
    """Load and return the AnnData object from a given path"""
    import scanpy as sc
    return sc.read(path)


def get_interacting_partners():
    
    interacting_partners = {
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
    return interacting_partners


def plot_histogram_for_pair(adata, tissue, gene1, gene2, cell_type1=None, cell_type2=None, cell_type_column="cell_types"):
    """
    Plots 2D histograms of spatial gene expression for two genes in a given tissue,
    with separate optional cell type filters for each gene.
    
    Parameters:
    - adata: AnnData object
    - tissue: tissue type to filter by
    - gene1, gene2: genes to plot
    - cell_type1: optional cell type filter for gene1 (if None, uses all cells)
    - cell_type2: optional cell type filter for gene2 (if None, uses all cells)  
    - cell_type_column: column name in adata.obs containing cell type information
    """
    # Import heavy modules only when this function is called
    import numpy as np
    import matplotlib.pyplot as plt
    
    os.makedirs('figures/temp', exist_ok=True)

    # 1) Filter by tissue
    tissue_adata = adata[adata.obs['tissue_type'] == tissue]
    if tissue_adata.n_obs == 0:
        print(f"❌ No data found for tissue type '{tissue}'")
        available_tissues = sorted(adata.obs['tissue_type'].unique())
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.text(0.5, 0.5, f"No data for tissue: {tissue}\nAvailable: {', '.join(map(str, available_tissues))}",
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(f"Error: Tissue '{tissue}' not found")
        return fig

    # 2) Batches within this tissue slice
    if 'batch' in tissue_adata.obs.columns:
        tissue_batches = sorted(tissue_adata.obs['batch'].unique())
    else:
        tissue_batches = [tissue]
        tissue_adata.obs['batch'] = tissue

    n_samples = len(tissue_batches)
    print(f"🧬 Found {n_samples} {tissue} samples: {tissue_batches}")

    fig, axes = plt.subplots(2, n_samples, figsize=(4*n_samples, 8), squeeze=False)
    if n_samples == 1:
        axes = axes.reshape(2, 1)

    for col_idx, batch_name in enumerate(tissue_batches):
        batch_adata = tissue_adata[tissue_adata.obs['batch'] == batch_name]
        print(f"🧬 {batch_name}: {batch_adata.n_obs} cells")

        if batch_adata.n_obs == 0:
            for row in range(2):
                axes[row, col_idx].text(0.5, 0.5, f"No cells in\n{batch_name}",
                                        ha='center', va='center', transform=axes[row, col_idx].transAxes)
                axes[row, col_idx].set_title(f"Error: {batch_name}")
            continue

        for row, (gene, cell_type_filter) in enumerate([(gene1, cell_type1), (gene2, cell_type2)]):
            ax = axes[row, col_idx]
            
            # Start with the batch data
            gene_adata = batch_adata.copy()
            
            # Apply cell type filter if specified
            if cell_type_filter is not None:
                if cell_type_column not in gene_adata.obs.columns:
                    ax.text(0.5, 0.5, f"Column '{cell_type_column}'\nnot found", 
                            ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f"Error: {gene} in {batch_name}")
                    continue
                
                # Handle categorical dtypes safely
                col = gene_adata.obs[cell_type_column]
                if hasattr(col, "cat"):
                    valid_types = list(col.cat.categories)
                else:
                    valid_types = sorted(map(str, col.unique()))
                
                # Filter by cell type
                gene_adata = gene_adata[col.astype(str) == str(cell_type_filter)]
                
                if gene_adata.n_obs == 0:
                    print(f"❌ No cells for {batch_name} with cell type '{cell_type_filter}' expressing {gene}")
                    ax.text(0.5, 0.5, f"No {cell_type_filter} cells\nexpressing {gene}", 
                            ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f"{gene} in {batch_name}")
                    continue

            if gene not in gene_adata.var.index:
                print(f"❌ Gene '{gene}' not found in {batch_name}")
                ax.text(0.5, 0.5, f"Gene '{gene}'\nnot found", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"{gene} in {batch_name}")
                continue

            # Get spatial coordinates and expression data
            x = gene_adata.obsm['X_spatial'][:, 0]
            y = gene_adata.obsm['X_spatial'][:, 1]
            
            gene_expr = gene_adata[:, gene].X
            if hasattr(gene_expr, 'A'):
                w = gene_expr.A.flatten()
            elif hasattr(gene_expr, 'toarray'):
                w = gene_expr.toarray().flatten()
            else:
                w = np.array(gene_expr).flatten()

            print(f"🧬 {gene} in {batch_name}" + (f" ({cell_type_filter})" if cell_type_filter else "") + 
                  f": min={w.min():.3f}, max={w.max():.3f}, mean={w.mean():.3f}, non-zero={np.sum(w > 0)}, cells={len(w)}")

            # Create histogram
            im = ax.hist2d(x, y, weights=w, bins=100, cmap='inferno')
            ax.set_title(f"{gene} in {batch_name}")
            ax.set_xlabel(""); ax.set_ylabel("")
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(im[3], ax=ax, label="Expression")

    plt.tight_layout()
    plt.subplots_adjust(left=0.04, bottom=0.08, right=0.96, top=0.95)
    fig.text(0.5, 0.02, 'Spatial X', ha='center', va='center', fontsize=14, fontweight='bold')
    fig.text(0.02, 0.5, 'Spatial Y', ha='center', va='center', fontsize=14, fontweight='bold', rotation=90)
    return fig

