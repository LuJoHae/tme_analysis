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
    
    print("Opening Biermann in memory...")
    start_time = time.time()
    adata = ad.read_h5ad(h5ad_path)
    print("Loaded in:", time.time() - start_time, "seconds")
    print("Adata shape:", adata.shape)
    
    cell_type_col = "cell_type"
    cell_types = sorted(list(adata.obs[cell_type_col].dropna().unique()))
    print("Cell types:", cell_types)
    
    start_time = time.time()
    ref_profiles = {}
    for ct in cell_types:
        if ct in ["Low-quality cells"]:
            continue
        print(f"Processing cell type: {ct}...")
        sub_adata = adata[adata.obs[cell_type_col] == ct]
        sub_X = sub_adata.X
        sum_profile = np.asarray(sub_X.sum(axis=0)).flatten()
        ref_profiles[ct] = sum_profile
        
    df_ref = pd.DataFrame(ref_profiles, index=adata.var_names)
    print("Reference profile shape:", df_ref.shape)
    print("Computed in:", time.time() - start_time, "seconds")
    print(df_ref.head())
    
if __name__ == "__main__":
    main()
