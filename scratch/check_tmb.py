import datalair
import anndata as ad
from ici_datasets.other_datasets import HarmonizedTcgaIAtlas

def main():
    lair = datalair.Lair("/storage/halu/lair")
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    
    adata_iatlas = ad.read_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad")
    print(adata_iatlas.obs.columns)
    
    if 'TMB' in adata_iatlas.obs.columns or 'tmb' in adata_iatlas.obs.columns or 'TMB_Score' in adata_iatlas.obs.columns:
        print("TMB is already in obs!")
    else:
        print("TMB is NOT in obs!")

if __name__ == "__main__":
    main()
