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
    parser.add_argument("--marker", required=True, help="Marker gene to subset lineage by (e.g., CD3D, ITGAM)")
    parser.add_argument("--exclude", required=False, default=None, help="Marker gene to explicitly exclude (e.g., CD3D)")
    args = parser.parse_args()
    
    in_path = Path(args.adata)
    out_path = Path(args.out_adata)
    
    print(f"Loading {in_path}...")
    adata = ad.read_h5ad(in_path)
    
    marker = args.marker
    print(f"Calculating mean {marker} expression per cluster in celltypist_leiden_1.0...")
    if marker not in adata.var_names:
        raise ValueError(f"{marker} not found in var_names!")
        
    marker_expr = adata[:, marker].X.toarray().flatten()
    df = pd.DataFrame({'cluster': adata.obs['celltypist_leiden_1.0'], 'expr': marker_expr})
    
    if args.exclude:
        if args.exclude not in adata.var_names:
            raise ValueError(f"Exclude marker {args.exclude} not found in var_names!")
        print(f"Also calculating mean {args.exclude} expression per cluster to exclude...")
        exclude_expr = adata[:, args.exclude].X.toarray().flatten()
        df['exclude_expr'] = exclude_expr
        means = df.groupby('cluster', observed=False).mean()
        threshold = 0.5
        # Keep if marker > threshold AND exclude <= threshold
        keep_clusters = means[(means['expr'] > threshold) & (means['exclude_expr'] <= threshold)].index.tolist()
        drop_clusters = means[(means['expr'] <= threshold) | (means['exclude_expr'] > threshold)].index.tolist()
    else:
        # Use observed=False to avoid future warnings
        means = df.groupby('cluster', observed=False)['expr'].mean()
        threshold = 0.5
        keep_clusters = means[means > threshold].index.tolist()
        drop_clusters = means[means <= threshold].index.tolist()
    
    print(f"\n--- Clusters KEEP (Mean {marker} > {threshold}" + (f" and {args.exclude} <= {threshold})" if args.exclude else ") ---"))
    for c in keep_clusters:
        if args.exclude:
            print(f"  {c}: {marker}={means.loc[c, 'expr']:.3f}, {args.exclude}={means.loc[c, 'exclude_expr']:.3f}")
        else:
            print(f"  {c}: {means[c]:.3f}")
        
    print(f"\n--- Clusters DROP (Mean {marker} <= {threshold}" + (f" or {args.exclude} > {threshold})" if args.exclude else ") ---"))
    for c in drop_clusters:
        if args.exclude:
            print(f"  {c}: {marker}={means.loc[c, 'expr']:.3f}, {args.exclude}={means.loc[c, 'exclude_expr']:.3f}")
        else:
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
