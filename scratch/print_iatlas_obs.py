import datalair
from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
import anndata as ad

def main():
    lair = datalair.Lair("/storage/halu/lair")
    ds = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds)
    adata_iatlas = ad.read_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad", backed='r')
    print("iAtlas OBS Columns:", list(adata_iatlas.obs.columns))
    
    # Try to find response-related columns and print their unique values
    resp_cols = [c for c in adata_iatlas.obs.columns if 'resp' in c.lower() or 'clin' in c.lower() or 'benefit' in c.lower()]
    print("Potential response columns:", resp_cols)
    for c in resp_cols:
        print(f"Unique values for {c}:", adata_iatlas.obs[c].unique().tolist())

if __name__ == "__main__":
    main()
