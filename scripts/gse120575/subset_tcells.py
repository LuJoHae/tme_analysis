import argparse
from pathlib import Path
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True, help="Path to input processed h5ad")
    parser.add_argument("--out-adata", required=True, help="Path to output subset h5ad")
    args = parser.parse_args()
    
    in_path = Path(args.adata)
    out_path = Path(args.out_adata)
    
    print(f"Loading {in_path}...")
    adata = ad.read_h5ad(in_path)
    
    print("Calculating mean CD3D expression per cluster in celltypist_leiden_1.0...")
    if 'CD3D' not in adata.var_names:
        raise ValueError("CD3D not found in var_names!")
        
    cd3d_expr = adata[:, 'CD3D'].X.toarray().flatten()
    df = pd.DataFrame({'cluster': adata.obs['celltypist_leiden_1.0'], 'cd3d': cd3d_expr})
    
    # Use observed=False to avoid future warnings
    means = df.groupby('cluster', observed=False)['cd3d'].mean()
    
    threshold = 0.5
    keep_clusters = means[means > threshold].index.tolist()
    drop_clusters = means[means <= threshold].index.tolist()
    
    print(f"\n--- Clusters KEEP (Mean CD3D > {threshold}) ---")
    for c in keep_clusters:
        print(f"  {c}: {means[c]:.3f}")
        
    print(f"\n--- Clusters DROP (Mean CD3D <= {threshold}) ---")
    for c in drop_clusters:
        print(f"  {c}: {means[c]:.3f}")
        
    print(f"\nSubsetting adata to {len(keep_clusters)} clusters...")
    adata = adata[adata.obs['celltypist_leiden_1.0'].isin(keep_clusters)].copy()
    print(f"Cells remaining: {adata.n_obs} (out of original {len(df)})")
    
    print("\nRecalculating Highly Variable Genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    
    print("Recalculating PCA...")
    sc.pp.pca(adata, n_comps=50)
    
    print("Recalculating Neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    
    print("Recalculating Leiden (resolutions 0.5, 1.0, 1.5, 2.0)...")
    for res in [0.5, 1.0, 1.5, 2.0]:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
        
    print("Recalculating UMAP...")
    sc.tl.umap(adata)
    
    print(f"Saving subset to {out_path}...")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    print("Done!")

if __name__ == "__main__":
    main()
