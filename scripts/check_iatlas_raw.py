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
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    for ds_name in iatlas_dataset_names:
        try:
            print(f"\nChecking {ds_name}...")
            data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
            
            # Find which file exists
            files = ["data_mrna_seq_expression.txt", "data_mrna_seq_tpm.txt", "data_mrna_seq_rpkm.txt"]
            for f in files:
                fpath = data_dir / f
                if fpath.exists():
                    df_raw = pd.read_csv(fpath, sep='\t', index_col=0)
                    print(f"File found: {f}, shape: {df_raw.shape}")
                    print(f"Dtype: {df_raw.dtypes.iloc[0]}")
                    print("Sample values (first 5x5):\n", df_raw.iloc[:5, :5])
                    
                    # Calculate sample sums
                    col_sums = df_raw.sum(axis=0)
                    print("Sample sums (first 5):", col_sums.iloc[:5].tolist())
                    print("Sample sums min/max/mean:", col_sums.min(), col_sums.max(), col_sums.mean())
                    
                    # Calculate mean of target genes
                    targets = ["B2M", "CD28", "MIF"]
                    for t in targets:
                        if t in df_raw.index:
                            val = df_raw.loc[t].mean()
                            print(f"Mean raw value of {t}: {val}")
                        else:
                            print(f"{t} not found in index.")
        except Exception as e:
            print(f"Error checking {ds_name}: {e}")

if __name__ == "__main__":
    main()
