import json
from pathlib import Path
import anndata as ad
import datalair
from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
from ml_pipelines.tcga_background import (
    compute_combined_pca_umap,
    perform_complete_clustering_analysis
)


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    output_dir = Path(".") / "output" / "TCGA-background" / "clustering"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load harmonized datasets
    print("Loading harmonized datasets from Datalair...")
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    
    adata = ad.read_h5ad(ds_path / "tcga_harmonized_uncorrected.h5ad")
    adata_iatlas = ad.read_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad")
    adata_corr = ad.read_h5ad(ds_path / "tcga_harmonized_combat.h5ad")
    adata_iatlas_corr = ad.read_h5ad(ds_path / "iatlas_harmonized_combat.h5ad")
    
    # Load gene lists
    print("Loading gene lists...")
    with open(ds_path / "de_hugo_genes.json", "r") as f:
        de_hugo_genes = json.load(f)
        
    # Load Bagaev signatures
    print("Loading Bagaev signatures...")
    from ici_datasets.bagaev_datasets import Signature
    from gene_utils import read_gene_sets
    ds_sig = Signature()
    lair.safe_derive(ds_sig)
    filepaths_sig = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths_sig["gene_signatures.gmt"])
    bagaev_genes = sorted(list(set.union(*[x.genes for x in bagaev_signature.values()])))
    print(f"Loaded {len(bagaev_genes)} unique Bagaev signature genes.")

    # 1. Uncorrected - DE genes
    print("Clustering DE uncorrected...")
    df_pca_de, df_umap_de, _ = compute_combined_pca_umap(
        adata, adata_iatlas, de_hugo_genes
    )
    perform_complete_clustering_analysis(
        df_pca_de, df_umap_de, output_dir, "DE_Uncorrected", "DE Genes (Uncorrected)"
    )

    # 2. Uncorrected - Bagaev signatures
    print("Clustering Bagaev uncorrected...")
    df_pca_bagaev, df_umap_bagaev, _ = compute_combined_pca_umap(
        adata, adata_iatlas, bagaev_genes
    )
    perform_complete_clustering_analysis(
        df_pca_bagaev, df_umap_bagaev, output_dir, "Bagaev_Uncorrected", "Bagaev (Uncorrected)"
    )

    # 3. ComBat Corrected - DE genes
    print("Clustering DE ComBat corrected...")
    df_pca_de_corr, df_umap_de_corr, _ = compute_combined_pca_umap(
        adata_corr, adata_iatlas_corr, de_hugo_genes
    )
    perform_complete_clustering_analysis(
        df_pca_de_corr, df_umap_de_corr, output_dir, "DE_ComBat", "DE Genes (ComBat)"
    )

    # 4. ComBat Corrected - Bagaev signatures
    print("Clustering Bagaev ComBat corrected...")
    df_pca_bagaev_corr, df_umap_bagaev_corr, _ = compute_combined_pca_umap(
        adata_corr, adata_iatlas_corr, bagaev_genes
    )
    perform_complete_clustering_analysis(
        df_pca_bagaev_corr, df_umap_bagaev_corr, output_dir, "Bagaev_ComBat", "Bagaev (ComBat)"
    )

    print("Clustering analysis and reporting successfully completed.")


if __name__ == "__main__":
    main()
