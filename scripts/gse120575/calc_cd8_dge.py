import argparse
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
import scanpy as sc

def calc_dge(adata_path: Path, out_path: Path) -> None:
    print(f"Loading {adata_path}...")
    adata = ad.read_h5ad(adata_path)
    
    # Ensure response column is clean
    if "characteristics: response" not in adata.obs.columns:
        print("Error: 'characteristics: response' not found in adata.obs")
        sys.exit(1)
        
    adata.obs["response"] = adata.obs["characteristics: response"]
    
    # Filter out cells with missing response data if any
    adata = adata[~adata.obs["response"].isna()].copy()
    
    # Verify both classes exist
    unique_responses = adata.obs["response"].unique()
    if "Responder" not in unique_responses or "Non-responder" not in unique_responses:
        print(f"Error: Need both Responder and Non-responder. Found: {unique_responses}")
        sys.exit(1)
        
    print(f"Running sc.tl.rank_genes_groups (Responder vs Non-responder)...")
    sc.tl.rank_genes_groups(
        adata, 
        groupby="response", 
        reference="Non-responder", 
        method="wilcoxon",
        pts=True # Calculate fraction of cells expressing
    )
    
    # Extract results for "Responder"
    results = sc.get.rank_genes_groups_df(adata, group="Responder")
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_path, index=False)
    print(f"Successfully saved {len(results)} DEGs to {out_path}")

def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate DEGs between Responder and Non-responder CD8+ T cells")
    parser.add_argument("--adata", required=True, help="Input CD8+ T-cell h5ad path")
    parser.add_argument("--out-csv", required=True, help="Output CSV path for DGE results")
    args = parser.parse_args()
    
    calc_dge(Path(args.adata), Path(args.out_csv))

if __name__ == "__main__":
    main()
