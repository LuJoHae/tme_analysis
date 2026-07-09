#!/usr/bin/env python3
"""
Load and visualize ImmunoCompass gene signatures using the datalair dataset class.
Generates overview plots showing:
1. Signature counts by Lineage and Broad Celltype Pathway
2. Distribution of signature sizes (number of genes)
3. Top 25 largest signatures
4. Heatmap of Jaccard overlap similarity between key immune signatures
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
from ici_datasets.other_datasets import ImmunoCompassSignatures
from gene_utils import read_gene_sets


def compute_jaccard_matrix(signatures):
    """
    Compute pairwise Jaccard similarity matrix for a dictionary of GeneSet objects.
    """
    keys = list(signatures.keys())
    n = len(keys)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            set1 = signatures[keys[i]].genes
            set2 = signatures[keys[j]].genes
            union = len(set1.union(set2))
            if union > 0:
                matrix[i, j] = len(set1.intersection(set2)) / union
            else:
                matrix[i, j] = 0.0
    return pd.DataFrame(matrix, index=keys, columns=keys)


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve output directory
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "immuno-compass-signatures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Derive dataset
    print("Deriving ImmunoCompassSignatures dataset...")
    ds = ImmunoCompassSignatures()
    lair.safe_derive(ds)
    
    filepaths = lair.get_dataset_filepaths(ds)
    txt_path = filepaths["concepts_defined.txt"]
    gmt_path = filepaths["gene_signatures.gmt"]
    
    # Load raw table for metadata plots
    df = pd.read_csv(txt_path, sep="\t")
    
    # Load GMT using package read_gene_sets
    signatures = read_gene_sets(gmt_path)
    print(f"Loaded {len(signatures)} ImmunoCompass gene signatures successfully.")
    
    # --- PLOT 1: Lineage & Pathway Distribution (2-panel horizontal) ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=300)
    
    # Lineage counts
    ax_lin = axes[0]
    lineage_counts = df['Lineage'].value_counts()
    sns.barplot(x=lineage_counts.values, y=lineage_counts.index, ax=ax_lin, palette='viridis')
    ax_lin.set_title("Gene Signature Counts by Cell Lineage", fontsize=12, fontweight='bold', color='#2C3E50')
    ax_lin.set_xlabel("Number of Signatures", fontsize=10)
    ax_lin.set_ylabel("Lineage", fontsize=10)
    ax_lin.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    for spine in ['top', 'right']:
        ax_lin.spines[spine].set_visible(False)
        
    # Broad Celltype Pathway counts (Top 15)
    ax_path = axes[1]
    pathway_counts = df['BroadCelltypePathway'].value_counts().head(15)
    sns.barplot(x=pathway_counts.values, y=pathway_counts.index, ax=ax_path, palette='mako')
    ax_path.set_title("Top 15 Broad Celltype Pathways", fontsize=12, fontweight='bold', color='#2C3E50')
    ax_path.set_xlabel("Number of Signatures", fontsize=10)
    ax_path.set_ylabel("Broad Celltype Pathway", fontsize=10)
    ax_path.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    for spine in ['top', 'right']:
        ax_path.spines[spine].set_visible(False)
        
    plt.tight_layout()
    plt.savefig(output_dir / "signature_distributions.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # --- PLOT 2: Gene Set Sizes ---
    plt.figure(figsize=(10, 5), dpi=300)
    sns.histplot(df['n_genes'], kde=True, color='#2980B9', bins=20)
    plt.title("Distribution of ImmunoCompass Gene Set Sizes", fontsize=12, fontweight='bold', color='#2C3E50')
    plt.xlabel("Number of Genes per Signature", fontsize=10)
    plt.ylabel("Count of Signatures", fontsize=10)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "gene_set_sizes_distribution.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # --- PLOT 3: Top 25 Largest Signatures ---
    plt.figure(figsize=(10, 8), dpi=300)
    df_largest = df.sort_values(by='n_genes', ascending=False).head(25)
    sns.barplot(x='n_genes', y='GeneSet', data=df_largest, palette='flare')
    plt.title("Top 25 Largest ImmunoCompass Signatures", fontsize=12, fontweight='bold', color='#2C3E50')
    plt.xlabel("Number of Genes", fontsize=10)
    plt.ylabel("Signature Name", fontsize=10)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "top_25_largest_signatures.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # --- PLOT 4: Jaccard Overlap Heatmap ---
    # Select a subset of representative cell-type signatures to show overlap
    selected_sig_names = [
        'Bcell_l_Danaher17', 'Bcell_sc', 'Plasma_sc',
        'CD4Tcell_sc', 'CD8Tcell_sc', 'RegTcell_sc', 'NKcell_sc',
        'cDC_sc', 'pDC_sc', 'Macrophage_sc', 'Monocyte_sc',
        'Mast_sc', 'Neutrophil_sc', 'Fibroblast_sc', 'Endothelial_sc'
    ]
    # Filter to signatures that actually exist
    selected_sig_names = [s for s in selected_sig_names if s in signatures]
    
    selected_sigs = {k: signatures[k] for k in selected_sig_names}
    df_jaccard = compute_jaccard_matrix(selected_sigs)
    
    plt.figure(figsize=(12, 10), dpi=300)
    sns.heatmap(df_jaccard, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Jaccard Similarity Index'}, square=True)
    plt.title("Jaccard Overlap Similarity between Key ImmunoCompass Signatures", fontsize=13, fontweight='bold', color='#2C3E50', pad=15)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_dir / "signature_overlap_heatmap.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # Save the GMT file copy to output folder for convenience
    df.to_csv(output_dir / "concepts_defined_metadata.csv", index=False)
    print(f"Overview plots and metadata successfully generated and saved to {output_dir}/")


if __name__ == "__main__":
    main()
