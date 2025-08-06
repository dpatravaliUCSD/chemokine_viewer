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

    sub = adata[adata.obs['batch'] == tissue]
    x = sub.obsm['X_spatial'][:, 0]
    y = sub.obsm['X_spatial'][:, 1]

    fig, axs = plt.subplots(2, 1, figsize=(8, 10), squeeze=False)
    axs = axs.flatten()

    for i, gene in enumerate([gene1, gene2]):
        w = sub[:, gene].X.A.flatten() if hasattr(sub[:, gene].X, 'A') else sub[:, gene].X.flatten()
        axs[i].hist2d(x, y, weights=w, bins=100, cmap='viridis')
        axs[i].set_title(f"{gene} in {tissue}")
        axs[i].set_xlabel("X")
        axs[i].set_ylabel("Y")

    fig.tight_layout()
    return fig
