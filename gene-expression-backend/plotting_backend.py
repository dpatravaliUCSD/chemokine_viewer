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
    with separate optional cell type filters for each gene, plus co-localization maps.
    
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
    from matplotlib.colors import ListedColormap
    
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

    # Create 3 rows now: gene1, gene2, and co-localization
    fig, axes = plt.subplots(3, n_samples, figsize=(4*n_samples, 12), squeeze=False)
    if n_samples == 1:
        axes = axes.reshape(3, 1)

    for col_idx, batch_name in enumerate(tissue_batches):
        batch_adata = tissue_adata[tissue_adata.obs['batch'] == batch_name]
        print(f"🧬 {batch_name}: {batch_adata.n_obs} cells")

        if batch_adata.n_obs == 0:
            for row in range(3):  # Now we have 3 rows
                axes[row, col_idx].text(0.5, 0.5, f"No cells in\n{batch_name}",
                                        ha='center', va='center', transform=axes[row, col_idx].transAxes)
                axes[row, col_idx].set_title(f"Error: {batch_name}")
            continue

        # Store expression data for co-localization
        gene1_data = None
        gene2_data = None
        spatial_coords = None

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

            # Store data for co-localization (use all cells in batch for co-localization)
            if row == 0:  # gene1
                # Get data for all cells in batch for gene1
                batch_gene1_adata = batch_adata.copy()
                if cell_type1 is not None and cell_type_column in batch_gene1_adata.obs.columns:
                    col = batch_gene1_adata.obs[cell_type_column]
                    batch_gene1_adata = batch_gene1_adata[col.astype(str) == str(cell_type1)]
                
                if gene1 in batch_gene1_adata.var.index and batch_gene1_adata.n_obs > 0:
                    gene1_x = batch_gene1_adata.obsm['X_spatial'][:, 0]
                    gene1_y = batch_gene1_adata.obsm['X_spatial'][:, 1]
                    gene1_expr = batch_gene1_adata[:, gene1].X
                    if hasattr(gene1_expr, 'A'):
                        gene1_w = gene1_expr.A.flatten()
                    elif hasattr(gene1_expr, 'toarray'):
                        gene1_w = gene1_expr.toarray().flatten()
                    else:
                        gene1_w = np.array(gene1_expr).flatten()
                    gene1_data = {'x': gene1_x, 'y': gene1_y, 'expr': gene1_w}
                    
            elif row == 1:  # gene2
                # Get data for all cells in batch for gene2
                batch_gene2_adata = batch_adata.copy()
                if cell_type2 is not None and cell_type_column in batch_gene2_adata.obs.columns:
                    col = batch_gene2_adata.obs[cell_type_column]
                    batch_gene2_adata = batch_gene2_adata[col.astype(str) == str(cell_type2)]
                
                if gene2 in batch_gene2_adata.var.index and batch_gene2_adata.n_obs > 0:
                    gene2_x = batch_gene2_adata.obsm['X_spatial'][:, 0]
                    gene2_y = batch_gene2_adata.obsm['X_spatial'][:, 1]
                    gene2_expr = batch_gene2_adata[:, gene2].X
                    if hasattr(gene2_expr, 'A'):
                        gene2_w = gene2_expr.A.flatten()
                    elif hasattr(gene2_expr, 'toarray'):
                        gene2_w = gene2_expr.toarray().flatten()
                    else:
                        gene2_w = np.array(gene2_expr).flatten()
                    gene2_data = {'x': gene2_x, 'y': gene2_y, 'expr': gene2_w}

            # Create histogram
            im = ax.hist2d(x, y, weights=w, bins=100, cmap='inferno')
            ax.set_title(f"{gene} in {batch_name}")
            ax.set_xlabel(""); ax.set_ylabel("")
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(im[3], ax=ax, label="Expression")

        # Create co-localization plot (row 2)
        colocalization_ax = axes[2, col_idx]
        
        if gene1_data is not None and gene2_data is not None:
            # Create co-localization map
            print(f"🎨 Creating co-localization map for {gene1} vs {gene2} in {batch_name}")
            
            # Get all spatial coordinates from the batch
            all_x = batch_adata.obsm['X_spatial'][:, 0]
            all_y = batch_adata.obsm['X_spatial'][:, 1]
            
            # Initialize co-localization array (0 = no expression, 1 = gene1 only, 2 = gene2 only, 3 = both)
            colocalization = np.zeros(len(all_x))
            
            # Create spatial coordinate to index mapping for the full batch
            coord_to_idx = {}
            for i, (x_coord, y_coord) in enumerate(zip(all_x, all_y)):
                coord_to_idx[(x_coord, y_coord)] = i
            
            # Mark cells expressing gene1
            if gene1_data is not None:
                for i, (x_coord, y_coord, expr) in enumerate(zip(gene1_data['x'], gene1_data['y'], gene1_data['expr'])):
                    if expr > 0:
                        if (x_coord, y_coord) in coord_to_idx:
                            idx = coord_to_idx[(x_coord, y_coord)]
                            colocalization[idx] = 1  # Gene1 only
            
            # Mark cells expressing gene2 (will overwrite to 2 for gene2 only, or add to 3 for both)
            if gene2_data is not None:
                for i, (x_coord, y_coord, expr) in enumerate(zip(gene2_data['x'], gene2_data['y'], gene2_data['expr'])):
                    if expr > 0:
                        if (x_coord, y_coord) in coord_to_idx:
                            idx = coord_to_idx[(x_coord, y_coord)]
                            if colocalization[idx] == 1:  # Already has gene1
                                colocalization[idx] = 3  # Both genes
                            elif colocalization[idx] == 0:  # No gene1
                                colocalization[idx] = 2  # Gene2 only
            
            # Filter out cells with no expression
            expressing_mask = colocalization > 0
            if np.sum(expressing_mask) > 0:
                plot_x = all_x[expressing_mask]
                plot_y = all_y[expressing_mask]
                plot_colors = colocalization[expressing_mask]
                
                # Create custom colormap: red=gene1 only, blue=gene2 only, purple=both
                colors = ['red', 'blue', 'purple']  # 1=red, 2=blue, 3=purple
                cmap = ListedColormap(colors)
                
                # Set black background to match histogram plots
                colocalization_ax.set_facecolor('black')
                
                # Create scatter plot
                scatter = colocalization_ax.scatter(plot_x, plot_y, c=plot_colors, cmap=cmap, 
                                                  s=1, alpha=0.7, vmin=1, vmax=3)
                
                # Count cells in each category
                gene1_only = np.sum(plot_colors == 1)
                gene2_only = np.sum(plot_colors == 2) 
                both_genes = np.sum(plot_colors == 3)
                
                title = f"Co-localization in {batch_name}"
                if cell_type1 or cell_type2:
                    title += f"\n({gene1_only} {gene1} only, {gene2_only} {gene2} only, {both_genes} both)"
                else:
                    title += f"\n({gene1_only} {gene1} only, {gene2_only} {gene2} only, {both_genes} both)"
                
                colocalization_ax.set_title(title)
                print(f"🎨 Co-localization: {gene1_only} {gene1} only, {gene2_only} {gene2} only, {both_genes} both")
                
                # Add custom legend
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor='red', label=f'{gene1} only'),
                                 Patch(facecolor='blue', label=f'{gene2} only'),
                                 Patch(facecolor='purple', label='Both genes')]
                colocalization_ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
                
            else:
                colocalization_ax.set_facecolor('black')
                colocalization_ax.text(0.5, 0.5, f"No cells expressing\n{gene1} or {gene2}", 
                                     ha='center', va='center', transform=colocalization_ax.transAxes, color='white')
                colocalization_ax.set_title(f"Co-localization in {batch_name}")
        else:
            colocalization_ax.set_facecolor('black')
            colocalization_ax.text(0.5, 0.5, f"Cannot create\nco-localization map", 
                                 ha='center', va='center', transform=colocalization_ax.transAxes, color='white')
            colocalization_ax.set_title(f"Co-localization in {batch_name}")
        
        colocalization_ax.set_xlabel(""); colocalization_ax.set_ylabel("")
        colocalization_ax.set_xticks([]); colocalization_ax.set_yticks([])

    plt.tight_layout()
    plt.subplots_adjust(left=0.04, bottom=0.06, right=0.96, top=0.95)
    fig.text(0.5, 0.02, 'Spatial X', ha='center', va='center', fontsize=14, fontweight='bold')
    fig.text(0.02, 0.5, 'Spatial Y', ha='center', va='center', fontsize=14, fontweight='bold', rotation=90)
    return fig

