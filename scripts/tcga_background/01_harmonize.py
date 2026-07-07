import json
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import datalair
import tcga
from ici_datasets.other_datasets import load_tcga, HarmonizedTcgaIAtlas
from ml_pipelines.tcga_background import (
    get_ensembl_union_exon_lengths,
    compute_rpk_tcga,
    select_top_one_vs_rest_de_genes,
    load_all_iatlas_datasets,
    aggregate_and_subset_by_hugo,
    find_source_de_genes
)


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    # Safe-derive TCGA dataset before loading
    print("Safe-deriving TCGA dataset...")
    ds = tcga.AllProjectsAdata()
    lair.safe_derive(ds)
    
    # Load TCGA AnnData
    print("Loading TCGA dataset...")
    adata = load_tcga(lair)
    adata.var_names_make_unique()
    print(f"Loaded TCGA data: {adata.shape[0]} samples, {adata.shape[1]} genes.")
    
    # Convert gene expression matrix to float32 numeric type
    print("Converting gene expression matrix to float32...")
    adata.X = adata.X.astype(np.float32)
    
    # Preprocessing: Basic filtering of zero-expression genes
    print("Filtering out genes with zero expression across all samples...")
    sc.pp.filter_genes(adata, min_cells=1)
    
    # Initialize the local Ensembl database manager
    print("Initializing Ensembl release 111 database...")
    import pyensembl
    ensembl = pyensembl.EnsemblRelease(release=111, species="human")
    ensembl.download()
    ensembl.index()
    
    print("Querying Ensembl union exon lengths...")
    gene_lengths = get_ensembl_union_exon_lengths(ensembl)
    
    # Calculate RPK on TCGA
    print("Computing RPK on TCGA...")
    adata = compute_rpk_tcga(adata, gene_lengths)

    # Map TCGA Ensembl IDs to Hugo symbols and filter biotypes
    print("Mapping TCGA Ensembl IDs to Hugo symbols and filtering biotypes...")
    gene_id_to_hugo = {}
    valid_ids = []
    for gene_id in adata.var_names:
        clean_id = gene_id.split(".")[0]
        try:
            gene_obj = ensembl.gene_by_id(clean_id)
            if gene_obj.biotype != "protein_coding":
                continue
            name = ensembl.gene_name_of_gene_id(clean_id)
            if name and name.startswith("MT-"):
                continue
            gene_id_to_hugo[gene_id] = name if name else gene_id
            valid_ids.append(gene_id)
        except ValueError:
            pass
            
    print(f"Filtered to {len(valid_ids)} protein-coding, non-mitochondrial genes.")
    adata = adata[:, valid_ids].copy()
    adata.var['hugo_symbol'] = [gene_id_to_hugo[g] for g in adata.var_names]
    
    # Aggregate TCGA to valid Hugo symbols to collapse duplicates BEFORE intersecting
    valid_hugo = list(set([h for h in adata.var['hugo_symbol'] if not h.startswith("ENSG")]))
    adata = aggregate_and_subset_by_hugo(adata, valid_hugo)

    # Load all iAtlas datasets
    print("Loading all iAtlas datasets...")
    adata_iatlas = load_all_iatlas_datasets(lair, ensembl, gene_lengths)
    print(f"Loaded combined iAtlas data: {adata_iatlas.shape[0]} samples, {adata_iatlas.shape[1]} genes.")
    
    # Intersect transcriptomes to prevent compositional bias
    print("Intersecting TCGA and iAtlas transcriptomes to prevent compositional bias...")
    common_transcriptome = sorted(list(set(adata.var_names).intersection(adata_iatlas.var_names)))
    print(f"Shared transcriptome size: {len(common_transcriptome)} genes.")
    
    adata = adata[:, common_transcriptome].copy()
    adata_iatlas = adata_iatlas[:, common_transcriptome].copy()
    
    # Perform TPM normalization on shared transcriptome
    print("Performing TPM normalisation (sum to 1e6) on shared transcriptome...")
    for ad_obj in [adata, adata_iatlas]:
        X = ad_obj.X
        if sp.issparse(X):
            sample_sums = np.array(X.sum(axis=1)).flatten()
            sample_sums[sample_sums == 0] = 1.0
            diag_inv_sums = sp.diags(1e6 / sample_sums)
            ad_obj.X = diag_inv_sums.dot(X).astype(np.float32)
        else:
            sample_sums = X.sum(axis=1, keepdims=True)
            sample_sums[sample_sums == 0] = 1.0
            ad_obj.X = ((X / sample_sums) * 1e6).astype(np.float32)
            
    # Apply log1p transformation
    print("Applying log1p transformation...")
    sc.pp.log1p(adata)
    sc.pp.log1p(adata_iatlas)
    
    # Find Source DE genes (TCGA vs iAtlas)
    print("Finding Source DE genes between TCGA and iAtlas...")
    source_de_genes = find_source_de_genes(adata, adata_iatlas, n_top_genes=100)
    
    # Select top 10 differentially expressed genes for each cohort against the rest on TCGA
    print("Selecting top one-vs-rest differentially expressed genes on TCGA...")
    de_hugo_genes = select_top_one_vs_rest_de_genes(adata, groupby="dataset", n_top_genes=10, method="wilcoxon")
    print(f"Total unique DE Hugo genes selected: {len(de_hugo_genes)}")

    # ---------------- COMBAT CORRECTION ---------------- #
    print("Concatenating for ComBat correction...")
    adata.obs['Source'] = 'TCGA'
    adata_iatlas.obs['Source'] = 'iAtlas'
    adata_combined_full = ad.concat([adata, adata_iatlas], axis=0, join="inner")
    
    print("Applying ComBat batch correction...")
    if sp.issparse(adata_combined_full.X):
        adata_combined_full.X = adata_combined_full.X.toarray()
    sc.pp.combat(adata_combined_full, key='Source')
    
    print("Splitting ComBat-corrected data back into TCGA and iAtlas objects...")
    adata_corr = adata_combined_full[adata_combined_full.obs['Source'] == 'TCGA'].copy()
    adata_iatlas_corr = adata_combined_full[adata_combined_full.obs['Source'] == 'iAtlas'].copy()

    # Save Harmonized Datasets to Datalair
    print("Saving harmonized datasets to Datalair...")
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    ds_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Writing h5ad files to {ds_path}...")
    adata.write_h5ad(ds_path / "tcga_harmonized_uncorrected.h5ad")
    adata_iatlas.write_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad")
    adata_corr.write_h5ad(ds_path / "tcga_harmonized_combat.h5ad")
    adata_iatlas_corr.write_h5ad(ds_path / "iatlas_harmonized_combat.h5ad")

    print(f"Writing gene lists to {ds_path}...")
    with open(ds_path / "de_hugo_genes.json", "w") as f:
        json.dump(de_hugo_genes, f)
    with open(ds_path / "source_de_genes.json", "w") as f:
        json.dump(source_de_genes, f)
        
    print("Harmonization successfully completed and datasets saved.")


if __name__ == "__main__":
    main()
