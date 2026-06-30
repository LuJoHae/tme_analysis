import pyensembl
import numpy as np
import pandas as pd
from collections import defaultdict
import datalair
import scanpy as sc
import anndata as ad
from ici_datasets.other_datasets import load_tcga
import importlib
import sys

# Append scripts folder
sys.path.append("scripts")
tcga_bg = importlib.import_module("TCGA-background")

def main():
    lair = datalair.Lair("/storage/halu/lair")
    
    # Load TCGA
    adata_tcga = load_tcga(lair)
    adata_tcga.var_names_make_unique()
    adata_tcga.X = adata_tcga.X.astype(np.float32)
    sc.pp.filter_genes(adata_tcga, min_cells=1)
    adata_tcga = tcga_bg.tpm_normalize_adata(adata_tcga)
    sc.pp.log1p(adata_tcga)
    
    # Map to Hugo
    ensembl = pyensembl.EnsemblRelease(release=111, species="human")
    ensembl.download()
    ensembl.index()
    gene_id_to_hugo = {}
    for gene_id in adata_tcga.var_names:
        clean_id = gene_id.split(".")[0]
        try:
            name = ensembl.gene_name_of_gene_id(clean_id)
            gene_id_to_hugo[gene_id] = name if name else gene_id
        except ValueError:
            gene_id_to_hugo[gene_id] = gene_id
    adata_tcga.var['hugo_symbol'] = [gene_id_to_hugo[g] for g in adata_tcga.var_names]
    
    # Load Bagaev signature genes
    from ici_datasets.bagaev_datasets import Signature
    from gene_utils import read_gene_sets
    ds_sig = Signature()
    lair.safe_derive(ds_sig)
    filepaths_sig = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths_sig["gene_signatures.gmt"])
    bagaev_genes = sorted(list(set.union(*[x.genes for x in bagaev_signature.values()])))
    
    # Load iAtlas
    adata_iatlas = tcga_bg.load_all_iatlas_datasets(lair)
    
    # Combined preprocessing subset (similar to run_pca_umap_combined)
    adata_tcga_sub = tcga_bg.aggregate_and_subset_by_hugo(adata_tcga, bagaev_genes)
    
    common_genes = sorted(list(
        set(adata_tcga_sub.var_names)
        .intersection(adata_iatlas.var_names)
    ))
    
    adata_tcga_sub = adata_tcga_sub[:, common_genes].copy()
    adata_iatlas_sub = adata_iatlas[:, common_genes].copy()
    
    adata_combined = ad.concat([adata_tcga_sub, adata_iatlas_sub], axis=0, join="inner")
    
    # Add source annotation
    adata_combined.obs["Source"] = ["TCGA" if "iAtlas" not in str(x) else "iAtlas" for x in adata_combined.obs["dataset"]]
    
    # Compute mean expression before scaling (log1p TPM unit-sum counts)
    mean_tcga = np.mean(adata_combined[adata_combined.obs["Source"] == "TCGA"].X, axis=0)
    mean_iatlas = np.mean(adata_combined[adata_combined.obs["Source"] == "iAtlas"].X, axis=0)
    
    # Scale features
    sc.pp.scale(adata_combined)
    
    # Compute PCA
    sc.tl.pca(adata_combined, n_comps=10, svd_solver="arpack")
    
    # Loadings of PC1
    loadings = adata_combined.varm['PCs'][:, 0]
    
    # Build summary dataframe
    df_loadings = pd.DataFrame({
        "Gene": common_genes,
        "PC1_Loading": loadings,
        "Mean_Expression_TCGA": mean_tcga,
        "Mean_Expression_iAtlas": mean_iatlas
    })
    
    # Sort by PC1 loading absolute value
    df_loadings["Abs_PC1_Loading"] = df_loadings["PC1_Loading"].abs()
    df_sorted = df_loadings.sort_values(by="PC1_Loading", ascending=False)
    
    print("\nTop 15 positive PC1 loading genes (higher in one group):")
    print(df_sorted.head(15).to_string(index=False))
    
    print("\nTop 15 negative PC1 loading genes (higher in the other group):")
    print(df_sorted.tail(15).to_string(index=False))
    
    # Let's check which genes are most different between the two sources
    df_loadings["Expression_Difference"] = df_loadings["Mean_Expression_TCGA"] - df_loadings["Mean_Expression_iAtlas"]
    df_diff_sorted = df_loadings.sort_values(by="Expression_Difference", ascending=False)
    
    print("\nTop 10 genes with highest TCGA vs iAtlas expression difference:")
    print(df_diff_sorted.head(10).to_string(index=False))
    
    print("\nTop 10 genes with highest iAtlas vs TCGA expression difference:")
    print(df_diff_sorted.tail(10).to_string(index=False))
    
    # Save the full results to a temporary CSV in output
    import os
    os.makedirs("output/analysis", exist_ok=True)
    df_loadings.to_csv("output/analysis/pc1_loadings_bagaev.csv", index=False)
    print("\nFull PC1 loadings saved to output/analysis/pc1_loadings_bagaev.csv")

if __name__ == "__main__":
    main()
