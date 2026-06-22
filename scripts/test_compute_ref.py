import datalair
from single_cell_datasets import SingleCellDataProcessStep05
import anndata as ad
import numpy as np
import pandas as pd
import time

def main():
    lair = datalair.Lair("/storage/halu/lair")
    step5 = SingleCellDataProcessStep05()
    filepaths = lair.get_dataset_filepaths(step5)
    h5ad_path = filepaths["BiermannDissectingTreatmentnaiveEcosystem2022Adata.h5ad"]
    
    print("Opening Biermann in backed mode...")
    adata = ad.read_h5ad(h5ad_path, backed='r')
    print("Adata shape:", adata.shape)
    
    cell_type_col = "cell_type"
    cell_types = sorted(list(adata.obs[cell_type_col].dropna().unique()))
    print("Cell types:", cell_types)
    
    # Let's measure time to compute reference profiles
    start_time = time.time()
    ref_profiles = {}
    for ct in cell_types:
        if ct in ["Low-quality cells"]:
            continue
        print(f"Processing cell type: {ct}...")
        indices = np.where(adata.obs[cell_type_col] == ct)[0]
        # In backed mode, slicing adata[indices, :] loads that submatrix
        sub_X = adata[indices, :].X
        
        # Check if sub_X is sparse or dense
        if hasattr(sub_X, "sum"):
            sum_profile = np.asarray(sub_X.sum(axis=0)).flatten()
        else:
            # Fallback if it's something else
            sum_profile = np.sum(sub_X, axis=0)
            
        ref_profiles[ct] = sum_profile
        
    df_ref = pd.DataFrame(ref_profiles, index=adata.var_names)
    print("Reference profile shape:", df_ref.shape)
    print("Computed in:", time.time() - start_time, "seconds")
    print(df_ref.head())
    
if __name__ == "__main__":
    main()
