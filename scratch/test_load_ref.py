import datalair
from single_cell_datasets._single_cell_reference import SingleCellReference
import anndata as ad
from pathlib import Path

def main():
    lair = datalair.Lair("/storage/halu/lair")
    ref_ds = SingleCellReference()
    
    print("Safe-deriving SingleCellReference...")
    lair.safe_derive(ref_ds)
    
    print("Getting file paths...")
    filepaths = lair.get_dataset_filepaths(ref_ds)
    print("Files:", filepaths)
    
    print("Reading AnnData...")
    adata = ad.read_h5ad(filepaths["adata.h5ad"])
    print("Adata shape:", adata.shape)
    print("Adata obs columns:", list(adata.obs.columns))
    print("Adata var columns:", list(adata.var.columns))
    
    # Print unique values of cell_type / celltype / cell_type_main if they exist
    for col in ["cell_type", "cell_type_main", "celltype"]:
        if col in adata.obs.columns:
            print(f"Unique values in {col}:", adata.obs[col].unique())
            
if __name__ == "__main__":
    main()
