#!/usr/bin/env python3
import sys
from pathlib import Path
import datalair
import anndata as ad

sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import single_cell_datasets

def main():
    lair_path = "/storage/halu/lair"
    print(f"Loading lair from {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    from single_cell_datasets._single_cell_reference import SingleCellReference
    ds = SingleCellReference()
    print("Loading SingleCellReference h5ad...")
    filepaths = lair.get_dataset_filepaths(ds)
    print("Files found:", filepaths)
    
    if "adata.h5ad" in filepaths:
        adata = ad.read_h5ad(filepaths["adata.h5ad"], backed="r")
        print("\nAnnData Summary:")
        print(adata)
        print("\nFirst few rows of obs:")
        print(adata.obs.head(10))
        print("\nColumns in obs:")
        print(list(adata.obs.columns))
        print("\nDimensions of X:", adata.X.shape)
    else:
        print("adata.h5ad not found.")

if __name__ == "__main__":
    main()
