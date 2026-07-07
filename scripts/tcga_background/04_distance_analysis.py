from pathlib import Path
import anndata as ad
import datalair
from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
from ml_pipelines.tcga_background import analyze_cohort_distances


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    output_dir = Path(".") / "output" / "TCGA-background" / "distances"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load harmonized datasets
    print("Loading harmonized datasets from Datalair...")
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    
    # distance analysis is run on combat-corrected combined dataset
    adata_corr = ad.read_h5ad(ds_path / "tcga_harmonized_combat.h5ad")
    adata_iatlas_corr = ad.read_h5ad(ds_path / "iatlas_harmonized_combat.h5ad")
    
    print("Concatenating ComBat-corrected TCGA and iAtlas data...")
    adata_corr.obs['Source'] = 'TCGA'
    adata_iatlas_corr.obs['Source'] = 'iAtlas'
    adata_combined_full = ad.concat([adata_corr, adata_iatlas_corr], axis=0, join="inner")
    
    print("Running cohort distance matching...")
    analyze_cohort_distances(
        adata_combined_full, 
        k_values=[10, 20, 50, 200, 500], 
        pca_dims=[10, 30, 50], 
        output_dir=output_dir
    )
    print("Cohort distance and voting analysis successfully completed.")


if __name__ == "__main__":
    main()
