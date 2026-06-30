import datalair
from single_cell_datasets import SingleCellDataProcessStep07
import anndata as ad

def main():
    lair = datalair.Lair("/storage/halu/lair")
    step7 = SingleCellDataProcessStep07()
    filepaths = lair.get_dataset_filepaths(step7)
    h5ad_path = filepaths["adata.h5ad"]
    print("Opening h5ad file in backed mode:", h5ad_path)
    
    adata = ad.read_h5ad(h5ad_path, backed='r')
    print("Success opening backed!")
    print("Adata shape:", adata.shape)
    print("Adata obs columns:", list(adata.obs.columns))
    print("Adata var columns:", list(adata.var.columns))
    
    # Check cell type column name
    cell_type_col = None
    for col in ["cell_type", "cell_type_main", "celltype", "cell_types"]:
        if col in adata.obs.columns:
            cell_type_col = col
            print(f"Found cell type column: {col}")
            print(adata.obs[col].value_counts())
            break
            
    if cell_type_col is None:
        print("No cell type column found!")
        
if __name__ == "__main__":
    main()
