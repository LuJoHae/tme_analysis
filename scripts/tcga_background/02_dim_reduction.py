import json
import shutil
from pathlib import Path
import anndata as ad
import datalair
from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
from ml_pipelines.tcga_background import (
    run_pca_umap_tcga_only,
    run_pca_umap_combined,
    run_lda_and_plot
)


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    output_dir = Path(".") / "output" / "TCGA-background" / "dimensionality_reduction"
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
    with open(ds_path / "source_de_genes.json", "r") as f:
        source_de_genes = json.load(f)
        
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

    # 1. TCGA-Only using DE genes
    print("Running TCGA-only analysis on DE genes...")
    target_de = [g for g in de_hugo_genes if g in adata.var_names]
    adata_tcga_de = adata[:, target_de].copy()
    run_pca_umap_tcga_only(
        adata_tcga_de, 
        output_dir, 
        "tcga_only_de_pca_umap.svg", 
        "TCGA-Only DE Genes"
    )
    
    # 2. TCGA-Only using Bagaev signature genes
    print("Running TCGA-only analysis on Bagaev signature genes...")
    target_bagaev = [g for g in bagaev_genes if g in adata.var_names]
    adata_tcga_bagaev = adata[:, target_bagaev].copy()
    run_pca_umap_tcga_only(
        adata_tcga_bagaev, 
        output_dir, 
        "tcga_only_bagaev_pca_umap.svg", 
        "TCGA-Only Bagaev Signature Genes"
    )

    # ---------------- UNCORRECTED RUN ---------------- #
    # 3. Combined TCGA & iAtlas using DE genes (Uncorrected)
    print("Running combined uncorrected analysis on DE genes...")
    run_pca_umap_combined(
        adata, adata_iatlas, de_hugo_genes, output_dir,
        "tcga_combined_de_pca_umap_uncorrected.svg", "TCGA & iAtlas DE Genes (Uncorrected)"
    )

    # 4. Combined TCGA & iAtlas using Bagaev signatures (Uncorrected)
    print("Running combined uncorrected analysis on Bagaev signatures...")
    run_pca_umap_combined(
        adata, adata_iatlas, bagaev_genes, output_dir,
        "tcga_combined_bagaev_pca_umap_uncorrected.svg", "TCGA & iAtlas Bagaev (Uncorrected)"
    )

    # 5. LDA on Source DE genes (Uncorrected)
    print("Running LDA on Source DE genes (Uncorrected)...")
    run_lda_and_plot(
        adata, adata_iatlas, source_de_genes, output_dir,
        "tcga_iatlas_source_de_lda_uncorrected.svg", "Source DE Genes (Uncorrected)"
    )

    # ---------------- CORRECTED RUN ---------------- #
    # 6. Combined TCGA & iAtlas using DE genes (ComBat)
    print("Running combined ComBat-corrected analysis on DE genes...")
    run_pca_umap_combined(
        adata_corr, adata_iatlas_corr, de_hugo_genes, output_dir,
        "tcga_combined_de_pca_umap_combat.svg", "TCGA & iAtlas DE Genes (ComBat)"
    )

    # 7. Combined TCGA & iAtlas using Bagaev signatures (ComBat)
    print("Running combined ComBat-corrected analysis on Bagaev signatures...")
    run_pca_umap_combined(
        adata_corr, adata_iatlas_corr, bagaev_genes, output_dir,
        "tcga_combined_bagaev_pca_umap_combat.svg", "TCGA & iAtlas Bagaev (ComBat)"
    )

    # 8. LDA on Source DE genes (ComBat)
    print("Running LDA on Source DE genes (ComBat)...")
    run_lda_and_plot(
        adata_corr, adata_iatlas_corr, source_de_genes, output_dir,
        "tcga_iatlas_source_de_lda_combat.svg", "Source DE Genes (ComBat)"
    )

    # Maintain original output file path copy to satisfy traditional checks
    print("Copying traditional PCA/UMAP SVG to tcga_pca_umap.svg...")
    shutil.copy(
        output_dir / "tcga_combined_de_pca_umap_uncorrected.svg",
        output_dir.parent / "tcga_pca_umap.svg"
    )
    print("Dimensional reduction and plotting successfully completed.")


if __name__ == "__main__":
    main()
