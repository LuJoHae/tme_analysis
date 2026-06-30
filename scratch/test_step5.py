import datalair
from single_cell_datasets import SingleCellDataProcessStep05
import anndata as ad

def main():
    lair = datalair.Lair("/storage/halu/lair")
    step5 = SingleCellDataProcessStep05()
    filepaths = lair.get_dataset_filepaths(step5)
    
    for name in ["JerbyArnonCancerCellProgram2018Adata.h5ad", "BiermannDissectingTreatmentnaiveEcosystem2022Adata.h5ad"]:
        h5ad_path = filepaths[name]
        print(f"=== {name} ===")
        adata = ad.read_h5ad(h5ad_path, backed='r')
        print("Shape:", adata.shape)
        print("Obs columns:", list(adata.obs.columns))
        print("Var columns:", list(adata.var.columns))
        
        # Check cell type column values if any
        for col in ["cell_type", "cell_type_main", "celltype", "cell_types"]:
            if col in adata.obs.columns:
                print(f"Value counts for {col}:")
                print(adata.obs[col].value_counts())
                break

if __name__ == "__main__":
    main()
