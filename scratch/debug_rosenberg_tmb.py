import sys
import importlib
import pandas as pd
import numpy as np
import datalair
import ici_datasets
import anndata as ad

sys.path.append("scripts")
iatlas_tmb = importlib.import_module("iAtlas-TMB")

def main():
    lair = datalair.Lair("/storage/halu/lair")
    
    # Read adata
    from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    adata_iatlas = ad.read_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad")
    
    rosenberg_samples = adata_iatlas.obs_names[adata_iatlas.obs['dataset'] == 'Rosenberg-iAtlas']
    print(f"Rosenberg samples in adata (first 5): {rosenberg_samples[:5].tolist()}")
    
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    ds_name = "blca_iatlas_imvigor210_2017" # The actual Rosenberg cBioPortal ID
    
    data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
    clin_sample = data_dir / "data_clinical_sample.txt"
    
    df_s = pd.read_csv(clin_sample, sep='\t', comment='#')
    print(f"Rosenberg samples in clin_sample (first 5): {df_s['SAMPLE_ID'].head(5).tolist()}")
    print("Top 10 TMB_NONSYNONYMOUS values:", df_s['TMB_NONSYNONYMOUS'].head(10).tolist())
    
    # Check matching
    matches = 0
    for s in rosenberg_samples:
        clean = s
        if clean not in df_s['SAMPLE_ID'].values and '-' in clean:
            clean = clean.rsplit('-', 1)[0]
        if clean in df_s['SAMPLE_ID'].values:
            matches += 1
            
    print(f"Matches found: {matches} out of {len(rosenberg_samples)}")
    
if __name__ == "__main__":
    main()
