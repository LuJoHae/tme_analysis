import argparse
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np

def subset_cd8_tcells(adata_path: Path, out_path: Path) -> None:
    print(f"Loading {adata_path}...")
    adata = ad.read_h5ad(adata_path, backed='r')
    
    gene = 'CD8A'
    cluster_col = 'leiden_1.0'
    
    if gene not in adata.var_names:
        print(f"Error: {gene} not found in adata.var_names")
        sys.exit(1)
        
    print(f"Calculating mean {gene} expression per cluster in {cluster_col}...")
    expr = adata[:, gene].X.toarray().flatten()
    df = pd.DataFrame({'cluster': adata.obs[cluster_col], 'expr': expr})
    means = df.groupby('cluster', observed=False)['expr'].mean().sort_values(ascending=False)
    
    keep_clusters = means[means > 0.5].index.tolist()
    drop_clusters = means[means <= 0.5].index.tolist()
    
    print("\n--- Clusters KEEP (Mean CD8A > 0.5) ---")
    for c in keep_clusters:
        print(f"  Cluster {c}: {means[c]:.3f}")
        
    print("\n--- Clusters DROP (Mean CD8A <= 0.5) ---")
    for c in drop_clusters:
        print(f"  Cluster {c}: {means[c]:.3f}")
        
    print(f"\nSubsetting adata to {len(keep_clusters)} clusters...")
    
    # Needs to be loaded into memory to subset and recalculate PCA properly
    adata = adata.to_memory()
    subset = adata[adata.obs[cluster_col].isin(keep_clusters)].copy()
    
    print(f"Cells remaining: {subset.n_obs} (out of original {adata.n_obs})")
    
    print("\nRecalculating Highly Variable Genes...")
    # sc.pp.highly_variable_genes requires log1p transformed data, which we assume is in .X
    sc.pp.highly_variable_genes(subset, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    
    print("Recalculating PCA...")
    sc.pp.pca(subset, svd_solver='arpack')
    
    print("Recalculating Neighbors...")
    sc.pp.neighbors(subset, n_neighbors=15, n_pcs=40)
    
    print("Recalculating Leiden (resolutions 0.5, 1.0, 1.5, 2.0)...")
    for res in [0.5, 1.0, 1.5, 2.0]:
        sc.tl.leiden(subset, resolution=res, key_added=f"leiden_{res}")
        
    print("Recalculating UMAP...")
    sc.tl.umap(subset)
    
    print(f"Saving subset to {out_path}...")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(out_path)
    print("Done!")

def main() -> None:
    parser = argparse.ArgumentParser(description="Subset CD8+ T cells based on CD8A expression")
    parser.add_argument("--adata", required=True, help="Input T-cell h5ad path")
    parser.add_argument("--out-adata", required=True, help="Output CD8+ T-cell h5ad path")
    args = parser.parse_args()
    
    subset_cd8_tcells(Path(args.adata), Path(args.out_adata))

if __name__ == "__main__":
    main()
