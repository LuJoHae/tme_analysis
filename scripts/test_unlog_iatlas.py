import datalair
import importlib
import sys
import numpy as np
import pandas as pd

sys.path.append("scripts")
tcga_bg = importlib.import_module("TCGA-background")

def main():
    lair = datalair.Lair("/storage/halu/lair")
    
    import ici_datasets
    sys.path.append("scripts")
    iatlas_tmb = importlib.import_module("iAtlas-TMB")
    
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    ds_name = "McDermott-iAtlas"
    
    data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
    fpath = data_dir / "data_mrna_seq_tpm.txt"
    df_raw = pd.read_csv(fpath, sep='\t', index_col=0)
    
    # Let's un-log: 2^x - 1
    df_linear = np.power(2.0, df_raw) - 1.0
    
    col_sums = df_linear.sum(axis=0)
    print("Unlogged sample sums (first 5):", col_sums.iloc[:5].tolist())
    print("Unlogged sample sums min/max/mean:", col_sums.min(), col_sums.max(), col_sums.mean())
    
    # Let's check Gide-iAtlas which is data_mrna_seq_expression.txt
    ds_name_gide = "Gide-iAtlas"
    data_dir_gide = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name_gide)
    fpath_gide = data_dir_gide / "data_mrna_seq_expression.txt"
    df_raw_gide = pd.read_csv(fpath_gide, sep='\t', index_col=0)
    df_linear_gide = np.power(2.0, df_raw_gide) - 1.0
    col_sums_gide = df_linear_gide.sum(axis=0)
    print("Unlogged Gide-iAtlas sample sums min/max/mean:", col_sums_gide.min(), col_sums_gide.max(), col_sums_gide.mean())

if __name__ == "__main__":
    main()
