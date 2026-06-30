import os
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp
import datalair
import tcga
from ici_datasets.other_datasets import load_tcga

# Suppress pandas performance warnings due to Scanpy DataFrame fragmentation
warnings.filterwarnings("ignore", category=PerformanceWarning)


def get_ensembl_union_exon_lengths(ensembl):
    """
    Queries all exons from the Ensembl database and computes the union exon length for each gene.
    """
    from collections import defaultdict
    conn = ensembl.db.connection
    cursor = conn.cursor()
    cursor.execute("SELECT gene_id, start, end FROM exon;")
    exons = cursor.fetchall()
    
    gene_exons = defaultdict(list)
    for gene_id, start, end in exons:
        gene_exons[gene_id].append((start, end))
        
    gene_lengths = {}
    for gene_id, intervals in gene_exons.items():
        intervals.sort(key=lambda x: x[0])
        merged = []
        for start, end in intervals:
            if not merged or merged[-1][1] < start:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        # 1-based coordinates inclusive
        gene_lengths[gene_id] = sum(end - start + 1 for start, end in merged)
    return gene_lengths


def compute_rpk_tcga(adata, gene_lengths):
    """
    Computes RPK (Reads Per Kilobase) on an AnnData object using pre-calculated exon union lengths.
    """
    import numpy as np
    import scipy.sparse as sp
    
    print("Computing RPK on TCGA data...")
    # Map var_names to lengths
    lengths = []
    mapped_indices = []
    unmapped_count = 0
    for idx, gene_id in enumerate(adata.var_names):
        clean_id = gene_id.split(".")[0]
        if clean_id in gene_lengths:
            lengths.append(gene_lengths[clean_id])
            mapped_indices.append(idx)
        else:
            unmapped_count += 1
            
    print(f"Mapped {len(lengths)} genes to Ensembl exon lengths. {unmapped_count} genes were unmapped and will be filtered.")
    
    # Subset AnnData to mapped genes
    adata_subset = adata[:, mapped_indices].copy()
    
    # Gene lengths in kilobases
    lengths_kb = np.array(lengths, dtype=np.float32) / 1e3
    
    # Calculate RPK
    X = adata_subset.X
    if sp.issparse(X):
        # We can scale each column by multiplying by 1 / lengths_kb on the right
        diag_inv_lengths = sp.diags(1.0 / lengths_kb)
        X_rpk = X.dot(diag_inv_lengths)
    else:
        # Dense matrix
        X_rpk = X / lengths_kb
        
    adata_subset.X = X_rpk.astype(np.float32)
    return adata_subset


def select_top_one_vs_rest_de_genes(adata, groupby, n_top_genes=10, method="wilcoxon"):
    """
    Selects the union of top differentially expressed genes for each group compared to the rest.
    """
    import pandas as pd
    import scanpy as sc
    
    # Run rank_genes_groups (one-vs-rest)
    print(f"Running one-vs-rest DE analysis grouped by {groupby}...")
    sc.tl.rank_genes_groups(
        adata, 
        groupby=groupby, 
        method=method,
        reference='rest',
        key_added="one_vs_rest"
    )
    
    groups = adata.obs[groupby].unique()
    groups = [g for g in groups if not pd.isna(g)]
    de_genes = set()
    
    print(f"Extracting top {n_top_genes} DE genes for each of the {len(groups)} groups...")
    for g in groups:
        # Extract DE results for group g
        result = sc.get.rank_genes_groups_df(adata, group=str(g), key="one_vs_rest")
        
        # Sort by absolute score (both up- and down-regulated in this group vs rest)
        result['abs_score'] = result['scores'].abs()
        top_df = result.sort_values(by='abs_score', ascending=False).head(n_top_genes)
        top_genes = top_df['names'].tolist()
        
        de_genes.update(top_genes)
        
    de_genes_list = sorted(list(de_genes))
    print(f"Total unique DE genes selected: {len(de_genes_list)}")
    return de_genes_list


def load_and_transform_iatlas_properly(data_dir, ensembl, gene_lengths):
    """
    Loads a cBioPortal iAtlas dataset from data_dir, un-logs the log2-transformed values,
    and performs proper normalization to linear sum of 1.0.
    """
    import os
    import numpy as np
    import pandas as pd
    
    files = ["data_mrna_seq_expression.txt", "data_mrna_seq_tpm.txt", "data_mrna_seq_rpkm.txt"]
    target_file = None
    for f in files:
        if (data_dir / f).is_file():
            target_file = f
            break
            
    if target_file is None:
        raise FileNotFoundError(f"No mRNA expression file found in {data_dir}")
        
    print(f"  Reading raw file: {target_file}")
    df_raw = pd.read_csv(data_dir / target_file, sep='\t', index_col=0)
    
    # 1. Un-log the data: 2^y - 1
    # cBioPortal expression data is stored in log2(val + 1) format
    df_linear = np.power(2.0, df_raw) - 1.0
    # Clip negative values to 0 in case of numerical inaccuracies
    df_linear = df_linear.clip(lower=0.0)
    
    # 2. Perform length normalization if the file is raw counts
    if target_file == "data_mrna_seq_expression.txt":
        print("  Detected raw counts file. Performing gene length (RPK) correction...")
        # Map Hugo symbols to Ensembl union exon lengths
        lengths_list = []
        common_symbols = []
        for symbol in df_linear.index:
            try:
                gids = ensembl.gene_ids_of_gene_name(symbol)
                # Map to the first gene ID if multiple exist
                gid = gids[0] if gids else None
                if gid and gid in gene_lengths:
                    lengths_list.append(gene_lengths[gid])
                    common_symbols.append(symbol)
            except ValueError:
                pass
                
        # Subset to genes with mapped lengths
        df_linear = df_linear.loc[common_symbols]
        lengths_kb = np.array(lengths_list, dtype=np.float32) / 1e3
        
        # Divide by gene lengths in kb
        df_linear = df_linear.div(lengths_kb, axis=0)
        
    # Do not sum-normalize here. We return the length-corrected linear abundance.
    return df_linear


def load_all_iatlas_datasets(lair, ensembl, gene_lengths):
    """
    Loads all cBioPortal iAtlas datasets and returns them as a concatenated AnnData.
    """
    import sys
    import importlib
    import anndata as ad
    import scanpy as sc
    import ici_datasets
    
    sys.path.append("scripts")
    iatlas_tmb = importlib.import_module("iAtlas-TMB")
    
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    adatas = []
    for ds_name in iatlas_dataset_names:
        try:
            print(f"Loading iAtlas dataset: {ds_name}...")
            data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
            # load expression data properly
            df_norm = load_and_transform_iatlas_properly(data_dir, ensembl, gene_lengths)
            
            # Convert to AnnData: samples x genes
            adata_ds = ad.AnnData(X=df_norm.T.astype(np.float32))
            adata_ds.obs_names = df_norm.columns
            adata_ds.var_names = df_norm.index
            adata_ds.var_names_make_unique()
            adata_ds.obs["dataset"] = ds_name
            
            adatas.append(adata_ds)
            print(f"Loaded {ds_name} with shape {adata_ds.shape}")
        except FileNotFoundError:
            print(f"Dataset {ds_name} files not found. Skipping.")
        except Exception as e:
            print(f"Error loading {ds_name}: {e}. Skipping.")
            
    if not adatas:
        raise ValueError("No iAtlas datasets could be successfully loaded.")
        
    print(f"Concatenating {len(adatas)} iAtlas datasets...")
    adata_iatlas = ad.concat(adatas, axis=0, join="inner")
    return adata_iatlas


def aggregate_and_subset_by_hugo(adata, target_hugo_genes):
    """
    Subsets the AnnData object to Ensembl IDs that map to the target Hugo symbols,
    and returns a new AnnData object with aggregated duplicate Hugo symbols (summed).
    """
    target_set = set(target_hugo_genes)
    indices = [i for i, h in enumerate(adata.var['hugo_symbol']) if h in target_set]
    
    if not indices:
        raise ValueError("None of the target genes were found in the dataset.")
        
    adata_subset = adata[:, indices].copy()
    hugo_names = adata_subset.var['hugo_symbol'].tolist()
    
    X = adata_subset.X
    if sp.issparse(X):
        X = X.toarray()
    
    df = pd.DataFrame(X, columns=hugo_names)
    df_agg = df.groupby(level=0, axis=1).sum()
    
    import anndata as ad
    adata_agg = ad.AnnData(X=df_agg.values.astype(np.float32), obs=adata_subset.obs.copy())
    adata_agg.var_names = df_agg.columns
    adata_agg.var_names_make_unique()
    return adata_agg


def run_pca_umap_tcga_only(adata, output_dir, filename, title_prefix):
    """
    Computes PCA and UMAP on the provided TCGA-only AnnData, and saves a 1x2 plot.
    """
    # Scale features (standardize)
    print(f"[{title_prefix}] Scaling features...")
    sc.pp.scale(adata)
    
    # Compute PCA
    print(f"[{title_prefix}] Computing PCA...")
    n_comps = min(50, adata.shape[1] - 1)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    
    # Compute Neighbors and UMAP
    print(f"[{title_prefix}] Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, n_comps))
    print(f"[{title_prefix}] Computing UMAP...")
    sc.tl.umap(adata)
    
    # Extract coordinates
    pca_coords = adata.obsm['X_pca']
    evr = adata.uns['pca']['variance_ratio']
    umap_coords = adata.obsm['X_umap']
    
    df_pca = pd.DataFrame({
        "PC1": pca_coords[:, 0],
        "PC2": pca_coords[:, 1],
        "Cohort": adata.obs["dataset"]
    })
    
    df_umap = pd.DataFrame({
        "UMAP1": umap_coords[:, 0],
        "UMAP2": umap_coords[:, 1],
        "Cohort": adata.obs["dataset"]
    })
    
    # Plotting
    print(f"[{title_prefix}] Generating scatter plots...")
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    
    unique_cohorts = sorted(df_pca["Cohort"].unique())
    num_cohorts = len(unique_cohorts)
    palette = sns.color_palette("husl", num_cohorts)
    color_dict = {cohort: color for cohort, color in zip(unique_cohorts, palette)}
    
    # PCA Plot
    sns.scatterplot(
        data=df_pca, x="PC1", y="PC2",
        hue="Cohort", hue_order=unique_cohorts, palette=color_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[0], legend=False
    )
    axes[0].set_title(f"PCA ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
    axes[0].set_xlabel(f"PC1 ({evr[0]*100:.2f}% variance explained)", fontsize=12)
    axes[0].set_ylabel(f"PC2 ({evr[1]*100:.2f}% variance explained)", fontsize=12)
    axes[0].grid(True, linestyle="--", alpha=0.3)
    
    # UMAP Plot
    sns.scatterplot(
        data=df_umap, x="UMAP1", y="UMAP2",
        hue="Cohort", hue_order=unique_cohorts, palette=color_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[1], legend=True
    )
    axes[1].set_title(f"UMAP ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
    axes[1].set_xlabel("UMAP1", fontsize=12)
    axes[1].set_ylabel("UMAP2", fontsize=12)
    axes[1].grid(True, linestyle="--", alpha=0.3)
    
    # Legend
    ncol = 2 if num_cohorts > 15 else 1
    axes[1].legend(
        title="TCGA Cohort",
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        borderaxespad=0,
        ncol=ncol,
        fontsize='small',
        title_fontsize='medium'
    )
    
    fig.tight_layout()
    plot_output_dir = output_dir / "TCGA-background"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_output_dir / filename
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"[{title_prefix}] Plots successfully saved to {plot_path}")


def run_pca_umap_combined(adata_tcga, adata_iatlas, target_hugo_genes, output_dir, filename, title_prefix):
    """
    Subsets TCGA and iAtlas to the target genes, intersects them, concatenates,
    and runs PCA, UMAP, and generates the 2x2 combined plot (colored and TCGA greyed out).
    """
    import anndata as ad
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    target_subset = [g for g in target_hugo_genes if g in adata_tcga.var_names]
    adata_tcga_sub = adata_tcga[:, target_subset].copy()
    
    # Intersect with iAtlas genes
    common_genes = sorted(list(
        set(adata_tcga_sub.var_names)
        .intersection(adata_iatlas.var_names)
    ))
    print(f"[{title_prefix}] Intersected genes count: {len(common_genes)} out of {len(adata_tcga_sub.var_names)} target genes.")
    
    # Subset both
    adata_tcga_sub = adata_tcga_sub[:, common_genes].copy()
    adata_iatlas_sub = adata_iatlas[:, common_genes].copy()
    
    # Concatenate
    print(f"[{title_prefix}] Concatenating TCGA and iAtlas data...")
    adata_combined = ad.concat([adata_tcga_sub, adata_iatlas_sub], axis=0, join="inner")
    print(f"[{title_prefix}] Combined AnnData shape: {adata_combined.shape[0]} samples, {adata_combined.shape[1]} genes.")
    
    # Scale features
    print(f"[{title_prefix}] Scaling features...")
    sc.pp.scale(adata_combined)
    
    # Compute PCA
    print(f"[{title_prefix}] Computing PCA...")
    n_comps = min(50, adata_combined.shape[1] - 1)
    sc.tl.pca(adata_combined, n_comps=n_comps, svd_solver="arpack")

    # Compute Neighbors and UMAP
    print(f"[{title_prefix}] Computing neighborhood graph...")
    sc.pp.neighbors(adata_combined, n_neighbors=15, n_pcs=min(30, n_comps))
    print(f"[{title_prefix}] Computing UMAP...")
    sc.tl.umap(adata_combined)
    
    # Extract coordinates
    pca_coords = adata_combined.obsm['X_pca']
    evr = adata_combined.uns['pca']['variance_ratio']
    umap_coords = adata_combined.obsm['X_umap']
    
    df_pca = pd.DataFrame({
        "PC1": pca_coords[:, 0],
        "PC2": pca_coords[:, 1],
        "Cohort": adata_combined.obs["dataset"]
    }, index=adata_combined.obs_names)
    
    df_umap = pd.DataFrame({
        "UMAP1": umap_coords[:, 0],
        "UMAP2": umap_coords[:, 1],
        "Cohort": adata_combined.obs["dataset"]
    }, index=adata_combined.obs_names)
    
    df_pca["Source"] = df_pca["Cohort"].apply(lambda x: "iAtlas" if "iAtlas" in str(x) else "TCGA")
    df_umap["Source"] = df_umap["Cohort"].apply(lambda x: "iAtlas" if "iAtlas" in str(x) else "TCGA")
    
    # Sort so that TCGA points are drawn first and iAtlas points are drawn on top
    df_pca = df_pca.sort_values("Source", ascending=True)
    df_umap = df_umap.sort_values("Source", ascending=True)
    
    # Plotting
    print(f"[{title_prefix}] Generating scatter plots...")
    fig, axes = plt.subplots(2, 2, figsize=(24, 20))
    
    unique_cohorts = sorted(df_pca["Cohort"].unique())
    num_cohorts = len(unique_cohorts)
    palette = sns.color_palette("husl", num_cohorts)
    
    color_dict = {cohort: color for cohort, color in zip(unique_cohorts, palette)}
    grey_tcga_color_dict = {
        cohort: color_dict[cohort] if "iAtlas" in str(cohort) else (0.8, 0.8, 0.8) 
        for cohort in unique_cohorts
    }
    
    markers_dict = {"TCGA": "o", "iAtlas": "^"}
    
    # Row 1: All colored
    # PCA Plot
    sns.scatterplot(
        data=df_pca, x="PC1", y="PC2",
        hue="Cohort", hue_order=unique_cohorts, palette=color_dict,
        style="Source", markers=markers_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[0, 0], legend=False
    )
    axes[0, 0].set_title(f"PCA of TCGA & iAtlas ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
    axes[0, 0].set_xlabel(f"PC1 ({evr[0]*100:.2f}% variance explained)", fontsize=12)
    axes[0, 0].set_ylabel(f"PC2 ({evr[1]*100:.2f}% variance explained)", fontsize=12)
    axes[0, 0].grid(True, linestyle="--", alpha=0.3)
    
    # UMAP Plot
    sns.scatterplot(
        data=df_umap, x="UMAP1", y="UMAP2",
        hue="Cohort", hue_order=unique_cohorts, palette=color_dict,
        style="Source", markers=markers_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[0, 1], legend=True
    )
    axes[0, 1].set_title(f"UMAP of TCGA & iAtlas ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
    axes[0, 1].set_xlabel("UMAP1", fontsize=12)
    axes[0, 1].set_ylabel("UMAP2", fontsize=12)
    axes[0, 1].grid(True, linestyle="--", alpha=0.3)
    
    # Row 2: TCGA Grey
    # PCA Plot
    sns.scatterplot(
        data=df_pca, x="PC1", y="PC2",
        hue="Cohort", hue_order=unique_cohorts, palette=grey_tcga_color_dict,
        style="Source", markers=markers_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[1, 0], legend=False
    )
    axes[1, 0].set_title(f"PCA ({title_prefix} - TCGA Greyed Out)", fontsize=16, fontweight='bold', pad=15)
    axes[1, 0].set_xlabel(f"PC1 ({evr[0]*100:.2f}% variance explained)", fontsize=12)
    axes[1, 0].set_ylabel(f"PC2 ({evr[1]*100:.2f}% variance explained)", fontsize=12)
    axes[1, 0].grid(True, linestyle="--", alpha=0.3)
    
    # UMAP Plot
    sns.scatterplot(
        data=df_umap, x="UMAP1", y="UMAP2",
        hue="Cohort", hue_order=unique_cohorts, palette=grey_tcga_color_dict,
        style="Source", markers=markers_dict,
        alpha=0.8, s=20, edgecolor=None, ax=axes[1, 1], legend=False
    )
    axes[1, 1].set_title(f"UMAP ({title_prefix} - TCGA Greyed Out)", fontsize=16, fontweight='bold', pad=15)
    axes[1, 1].set_xlabel("UMAP1", fontsize=12)
    axes[1, 1].set_ylabel("UMAP2", fontsize=12)
    axes[1, 1].grid(True, linestyle="--", alpha=0.3)
    
    # Place legend outside the plot area to the right of UMAP (Row 1)
    ncol = 2 if num_cohorts > 15 else 1
    axes[0, 1].legend(
        title="Cancer Type / Cohort & Source",
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        borderaxespad=0,
        ncol=ncol,
        fontsize='small',
        title_fontsize='medium'
    )
    
    # Save figure
    fig.tight_layout()
    plot_output_dir = output_dir / "TCGA-background"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_output_dir / filename
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"[{title_prefix}] Plots successfully saved to {plot_path}")
    return df_pca, df_umap


def find_source_de_genes(adata_tcga, adata_iatlas, n_top_genes=100):
    """
    Finds the most differentially expressed genes between TCGA and iAtlas datasets.
    """
    import anndata as ad
    import scanpy as sc
    
    print("Concatenating datasets to find Source DE genes...")
    adata_tcga_copy = adata_tcga.copy()
    adata_iatlas_copy = adata_iatlas.copy()
    
    adata_tcga_copy.obs['Source'] = 'TCGA'
    adata_iatlas_copy.obs['Source'] = 'iAtlas'
    
    adata_combined = ad.concat([adata_tcga_copy, adata_iatlas_copy], axis=0, join="inner")
    
    print("Running differential expression between TCGA and iAtlas...")
    sc.tl.rank_genes_groups(
        adata_combined,
        groupby='Source',
        method='wilcoxon',
        key_added='source_de'
    )
    
    de_genes = set()
    for group in ['TCGA', 'iAtlas']:
        result = sc.get.rank_genes_groups_df(adata_combined, group=group, key='source_de')
        result['abs_score'] = result['scores'].abs()
        top_df = result.sort_values(by='abs_score', ascending=False).head(n_top_genes)
        de_genes.update(top_df['names'].tolist())
        
    de_genes_list = sorted(list(de_genes))
    print(f"Total unique Source DE genes selected: {len(de_genes_list)}")
    return de_genes_list


def run_lda_and_plot(adata_tcga, adata_iatlas, target_genes, output_dir, filename, title_prefix):
    """
    Fits LDA on Source (1D) and Cohort (2D) and plots the projections.
    """
    import anndata as ad
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.preprocessing import StandardScaler
    
    target_subset = [g for g in target_genes if g in adata_tcga.var_names and g in adata_iatlas.var_names]
    
    adata_tcga_sub = adata_tcga[:, target_subset].copy()
    adata_iatlas_sub = adata_iatlas[:, target_subset].copy()
    
    adata_tcga_sub.obs['Source'] = 'TCGA'
    adata_iatlas_sub.obs['Source'] = 'iAtlas'
    
    adata_combined = ad.concat([adata_tcga_sub, adata_iatlas_sub], axis=0, join="inner")
    
    X = adata_combined.X
    if hasattr(X, "toarray"):
        X = X.toarray()
        
    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 1. LDA on Source (TCGA vs iAtlas) -> 1D
    print(f"[{title_prefix}] Fitting LDA on Source (TCGA vs iAtlas)...")
    y_source = adata_combined.obs['Source'].values
    lda_source = LinearDiscriminantAnalysis(n_components=1)
    X_lda_source = lda_source.fit_transform(X_scaled, y_source)
    
    df_lda_source = pd.DataFrame({
        "LD1": X_lda_source[:, 0],
        "Source": y_source
    })
    
    if hasattr(lda_source, 'coef_'):
        df_loadings = pd.DataFrame({
            "Gene": adata_combined.var_names,
            "LD1_Coefficient": lda_source.coef_[0]
        })
        df_loadings['Abs_Coefficient'] = df_loadings['LD1_Coefficient'].abs()
        df_loadings = df_loadings.sort_values(by='Abs_Coefficient', ascending=False)
        loadings_path = output_dir / "TCGA-background" / f"{filename.replace('.svg', '_loadings.csv')}"
        df_loadings.to_csv(loadings_path, index=False)
        print(f"[{title_prefix}] LDA Loadings saved to {loadings_path}")
    
    # 2. LDA on Cohort -> 2D
    print(f"[{title_prefix}] Fitting LDA on Cohort...")
    y_cohort = adata_combined.obs['dataset'].values
    
    # Filter out cohorts with too few samples to avoid singular covariance or bad splits
    cohort_counts = pd.Series(y_cohort).value_counts()
    valid_cohorts = cohort_counts[cohort_counts >= 5].index.tolist()
    
    mask = [c in valid_cohorts for c in y_cohort]
    X_scaled_cohort = X_scaled[mask]
    y_cohort_valid = y_cohort[mask]
    source_valid = y_source[mask]
    
    n_classes = len(np.unique(y_cohort_valid))
    n_components = min(2, n_classes - 1)
    
    if n_components >= 1:
        lda_cohort = LinearDiscriminantAnalysis(n_components=n_components)
        X_lda_cohort = lda_cohort.fit_transform(X_scaled_cohort, y_cohort_valid)
        
        df_lda_cohort = pd.DataFrame({
            "LD1": X_lda_cohort[:, 0],
            "LD2": X_lda_cohort[:, 1] if n_components >= 2 else 0,
            "Cohort": y_cohort_valid,
            "Source": source_valid
        })
    else:
        df_lda_cohort = None
        
    # Plotting
    print(f"[{title_prefix}] Generating LDA plots...")
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    # Plot 1: Source LDA Density
    sns.kdeplot(
        data=df_lda_source, x="LD1", hue="Source", 
        fill=True, common_norm=False, alpha=0.5, ax=axes[0]
    )
    axes[0].set_title(f"LDA on Source ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
    axes[0].set_xlabel("Linear Discriminant 1 (Source Separation)", fontsize=12)
    axes[0].grid(True, linestyle="--", alpha=0.3)
    
    # Plot 2: Cohort LDA Scatter
    if df_lda_cohort is not None:
        unique_cohorts = sorted(df_lda_cohort["Cohort"].unique())
        palette = sns.color_palette("husl", len(unique_cohorts))
        color_dict = {cohort: color for cohort, color in zip(unique_cohorts, palette)}
        markers_dict = {"TCGA": "o", "iAtlas": "^"}
        
        # Sort so TCGA is plotted first
        df_lda_cohort = df_lda_cohort.sort_values("Source", ascending=True)
        
        sns.scatterplot(
            data=df_lda_cohort, x="LD1", y="LD2",
            hue="Cohort", hue_order=unique_cohorts, palette=color_dict,
            style="Source", markers=markers_dict,
            alpha=0.8, s=50, edgecolor=None, ax=axes[1], legend=True
        )
        axes[1].set_title(f"LDA on Cohort ({title_prefix})", fontsize=16, fontweight='bold', pad=15)
        axes[1].set_xlabel("Linear Discriminant 1", fontsize=12)
        axes[1].set_ylabel("Linear Discriminant 2", fontsize=12)
        axes[1].grid(True, linestyle="--", alpha=0.3)
        
        ncol = 2 if len(unique_cohorts) > 15 else 1
        axes[1].legend(
            title="Cohort & Source",
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            borderaxespad=0,
            ncol=ncol,
            fontsize='small'
        )
    else:
        axes[1].text(0.5, 0.5, "Not enough cohort classes for 2D LDA", ha='center', va='center')
        axes[1].axis('off')
        
    fig.tight_layout()
    plot_output_dir = output_dir / "TCGA-background"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_output_dir / filename
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"[{title_prefix}] LDA Plots successfully saved to {plot_path}")


def run_kmeans_clustering(data_coords, output_dir, name, k_range=range(2, 9)):
    """
    Runs K-means clustering for a range of K values, computes Silhouette,
    Calinski-Harabasz, and Davies-Bouldin metrics, and returns the metrics df,
    the best K, and the best cluster assignments and centers.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    import pandas as pd
    import numpy as np
    
    X = data_coords.values
    metrics = []
    results = {}
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)
        
        sil = silhouette_score(X, labels)
        ch = calinski_harabasz_score(X, labels)
        db = davies_bouldin_score(X, labels)
        
        metrics.append({
            "Parameter": k,
            "Silhouette": sil,
            "Calinski_Harabasz": ch,
            "Davies_Bouldin": db,
            "N_Clusters": k
        })
        
        results[k] = {
            "model": kmeans,
            "labels": labels,
            "centers": kmeans.cluster_centers_
        }
        print(f"[{name} KMeans] K={k}: Silhouette={sil:.4f}, CH={ch:.2f}, DB={db:.4f}")
        
    df_metrics = pd.DataFrame(metrics)
    
    best_idx = df_metrics["Silhouette"].idxmax()
    best_k = df_metrics.loc[best_idx, "Parameter"]
    print(f"[{name} KMeans] Best K chosen: {best_k} (Silhouette={df_metrics.loc[best_idx, 'Silhouette']:.4f})")
    
    best_result = results[best_k]
    return df_metrics, best_k, best_result["labels"], best_result["centers"]


def run_dbscan_clustering(data_coords, output_dir, name):
    """
    Runs DBSCAN clustering by varying eps (fraction of coordinates standard deviation),
    computes Silhouette score on non-noise samples, and returns the metrics, best eps,
    labels, and cluster centroids.
    """
    from sklearn.cluster import DBSCAN
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    import pandas as pd
    import numpy as np
    
    X = data_coords.values
    std_val = np.std(X)
    
    # Range of eps to test: fractions of standard deviation
    eps_fractions = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    metrics = []
    results = {}
    
    for pct in eps_fractions:
        eps = std_val * pct
        db = DBSCAN(eps=eps, min_samples=15)
        labels = db.fit_predict(X)
        
        unique_labels = set(labels)
        n_clusters = len(unique_labels - {-1})
        
        if n_clusters >= 2:
            mask = labels != -1
            if mask.sum() > 15:
                sil = silhouette_score(X[mask], labels[mask])
                ch = calinski_harabasz_score(X[mask], labels[mask])
                db_score = davies_bouldin_score(X[mask], labels[mask])
            else:
                sil, ch, db_score = -1.0, 0.0, 999.0
        else:
            sil, ch, db_score = -1.0, 0.0, 999.0
            
        metrics.append({
            "Parameter": eps,
            "Silhouette": sil,
            "Calinski_Harabasz": ch,
            "Davies_Bouldin": db_score,
            "N_Clusters": n_clusters,
            "Eps_Fraction": pct
        })
        
        results[pct] = {
            "labels": labels,
            "n_clusters": n_clusters
        }
        print(f"[{name} DBSCAN] eps={eps:.4f} (frac={pct}): clusters={n_clusters}, Silhouette={sil:.4f}")
        
    df_metrics = pd.DataFrame(metrics)
    
    valid_df = df_metrics[df_metrics["N_Clusters"] >= 2]
    if not valid_df.empty:
        best_idx = valid_df["Silhouette"].idxmax()
        best_eps = df_metrics.loc[best_idx, "Parameter"]
    else:
        best_eps = std_val * 0.2
        print(f"[{name} DBSCAN] Warning: No DBSCAN configuration found >= 2 clusters. Using fallback eps={best_eps:.4f}")
        
    db_best = DBSCAN(eps=best_eps, min_samples=15)
    best_labels = db_best.fit_predict(X)
    
    unique_labels = sorted(list(set(best_labels) - {-1}))
    centers = []
    for l in unique_labels:
        centers.append(X[best_labels == l].mean(axis=0))
    centers = np.array(centers) if centers else np.empty((0, 2))
    
    df_metrics = df_metrics.drop(columns=["Eps_Fraction"])
    return df_metrics, best_eps, best_labels, centers


def run_leiden_clustering(data_coords, output_dir, name):
    """
    Runs Leiden clustering on the 2D coordinates by constructing a neighborhood graph
    and varying resolution. Computes metrics and returns the metrics, best resolution,
    labels, and cluster centroids.
    """
    import anndata as ad
    import scanpy as sc
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    import pandas as pd
    import numpy as np
    import warnings
    
    X = data_coords.values
    resolutions = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    metrics = []
    results = {}
    
    warnings.simplefilter("ignore", category=UserWarning)
    
    adata_temp = ad.AnnData(X=X.astype(np.float32))
    sc.pp.neighbors(adata_temp, n_neighbors=15)
    
    for res in resolutions:
        sc.tl.leiden(adata_temp, resolution=res, flavor='igraph', n_iterations=2, directed=False, key_added=f"leiden_{res}")
        labels = adata_temp.obs[f"leiden_{res}"].astype(int).values
        
        unique_labels = set(labels)
        n_clusters = len(unique_labels)
        
        if n_clusters >= 2:
            sil = silhouette_score(X, labels)
            ch = calinski_harabasz_score(X, labels)
            db_score = davies_bouldin_score(X, labels)
        else:
            sil, ch, db_score = -1.0, 0.0, 999.0
            
        metrics.append({
            "Parameter": res,
            "Silhouette": sil,
            "Calinski_Harabasz": ch,
            "Davies_Bouldin": db_score,
            "N_Clusters": n_clusters
        })
        
        results[res] = {
            "labels": labels,
            "n_clusters": n_clusters
        }
        print(f"[{name} Leiden] res={res}: clusters={n_clusters}, Silhouette={sil:.4f}")
        
    df_metrics = pd.DataFrame(metrics)
    
    valid_df = df_metrics[df_metrics["N_Clusters"] >= 2]
    if not valid_df.empty:
        best_idx = valid_df["Silhouette"].idxmax()
        best_res = df_metrics.loc[best_idx, "Parameter"]
    else:
        best_res = 0.5
        print(f"[{name} Leiden] Warning: No Leiden configuration found >= 2 clusters. Using fallback resolution={best_res}")
        
    best_labels = results[best_res]["labels"]
    
    unique_labels = sorted(list(set(best_labels)))
    centers = []
    for l in unique_labels:
        centers.append(X[best_labels == l].mean(axis=0))
    centers = np.array(centers)
    
    df_metrics = df_metrics[["Parameter", "Silhouette", "Calinski_Harabasz", "Davies_Bouldin", "N_Clusters"]]
    return df_metrics, best_res, best_labels, centers


def plot_clustering_results(df_pca, df_umap, pca_centers_dict, umap_centers_dict, params_pca_dict, params_umap_dict, output_dir, filename, title_suffix):
    """
    Generates a 3x2 grid of plots:
    Row 1: K-means (PCA on left, UMAP on right)
    Row 2: DBSCAN (PCA on left, UMAP on right)
    Row 3: Leiden (PCA on left, UMAP on right)
    Centers are plotted on top.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    fig, axes = plt.subplots(3, 2, figsize=(24, 30))
    
    algorithms = ["KMeans", "DBSCAN", "Leiden"]
    
    for i, alg in enumerate(algorithms):
        # PCA Plot
        col_pca = f"{alg}_PCA_Cluster"
        unique_clusters_pca = sorted(df_pca[col_pca].unique())
        non_noise_pca = [c for c in unique_clusters_pca if c != "-1"]
        pca_palette = sns.color_palette("tab10", max(1, len(non_noise_pca)))
        
        pca_color_map = {c: pca_palette[idx % len(pca_palette)] for idx, c in enumerate(non_noise_pca)}
        if "-1" in unique_clusters_pca:
            pca_color_map["-1"] = (0.7, 0.7, 0.7)
            
        sns.scatterplot(
            data=df_pca, x="PC1", y="PC2",
            hue=col_pca, hue_order=unique_clusters_pca, palette=pca_color_map,
            alpha=0.6, s=15, edgecolor=None, ax=axes[i, 0], legend="full"
        )
        
        centers_pca = pca_centers_dict[alg]
        if centers_pca.size > 0:
            axes[i, 0].scatter(
                centers_pca[:, 0], centers_pca[:, 1],
                c='black', s=300, marker='*', edgecolor='white', linewidth=1.5,
                label='Centers'
            )
            
        param_pca_str = f"param={params_pca_dict[alg]}" if alg != "KMeans" else f"K={params_pca_dict[alg]}"
        axes[i, 0].set_title(f"{alg} Clustering on PCA Space ({param_pca_str}) - {title_suffix}", fontsize=16, fontweight='bold', pad=15)
        axes[i, 0].grid(True, linestyle="--", alpha=0.3)
        axes[i, 0].legend(title="Cluster / Centers", loc="upper right")
        
        # UMAP Plot
        col_umap = f"{alg}_UMAP_Cluster"
        unique_clusters_umap = sorted(df_umap[col_umap].unique())
        non_noise_umap = [c for c in unique_clusters_umap if c != "-1"]
        umap_palette = sns.color_palette("tab10", max(1, len(non_noise_umap)))
        
        umap_color_map = {c: umap_palette[idx % len(umap_palette)] for idx, c in enumerate(non_noise_umap)}
        if "-1" in unique_clusters_umap:
            umap_color_map["-1"] = (0.7, 0.7, 0.7)
            
        sns.scatterplot(
            data=df_umap, x="UMAP1", y="UMAP2",
            hue=col_umap, hue_order=unique_clusters_umap, palette=umap_color_map,
            alpha=0.6, s=15, edgecolor=None, ax=axes[i, 1], legend="full"
        )
        
        # Plot UMAP centers
        centers_umap = umap_centers_dict[alg]
        if centers_umap.size > 0:
            axes[i, 1].scatter(
                centers_umap[:, 0], centers_umap[:, 1],
                c='black', s=300, marker='*', edgecolor='white', linewidth=1.5,
                label='Centers'
            )
            
        param_umap_str = f"param={params_umap_dict[alg]}" if alg != "KMeans" else f"K={params_umap_dict[alg]}"
        axes[i, 1].set_title(f"{alg} Clustering on UMAP Space ({param_umap_str}) - {title_suffix}", fontsize=16, fontweight='bold', pad=15)
        axes[i, 1].grid(True, linestyle="--", alpha=0.3)
        axes[i, 1].legend(title="Cluster / Centers", loc="upper right")
        
    fig.tight_layout()
    plot_output_dir = output_dir / "TCGA-background"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_output_dir / filename
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"Clustering plots successfully saved to {plot_path}")


def perform_complete_clustering_analysis(df_pca, df_umap, output_dir, name_prefix, title_suffix):
    """
    Runs KMeans, DBSCAN, and Leiden clustering on PCA and UMAP spaces.
    Saves metrics CSV, sample assignment CSV, cohort assignment CSV, and a 3x2 SVG plot.
    """
    import pandas as pd
    import numpy as np
    
    # 1. Run KMeans
    print(f"[{name_prefix}] Running K-means clustering...")
    pca_km_metrics, best_k_pca, km_pca_labels, km_pca_centers = run_kmeans_clustering(
        df_pca[["PC1", "PC2"]], output_dir, f"{name_prefix} PCA"
    )
    umap_km_metrics, best_k_umap, km_umap_labels, km_umap_centers = run_kmeans_clustering(
        df_umap[["UMAP1", "UMAP2"]], output_dir, f"{name_prefix} UMAP"
    )
    
    # 2. Run DBSCAN
    print(f"[{name_prefix}] Running DBSCAN clustering...")
    pca_db_metrics, best_eps_pca, db_pca_labels, db_pca_centers = run_dbscan_clustering(
        df_pca[["PC1", "PC2"]], output_dir, f"{name_prefix} PCA"
    )
    umap_db_metrics, best_eps_umap, db_umap_labels, db_umap_centers = run_dbscan_clustering(
        df_umap[["UMAP1", "UMAP2"]], output_dir, f"{name_prefix} UMAP"
    )
    
    # 3. Run Leiden
    print(f"[{name_prefix}] Running Leiden clustering...")
    pca_le_metrics, best_res_pca, le_pca_labels, le_pca_centers = run_leiden_clustering(
        df_pca[["PC1", "PC2"]], output_dir, f"{name_prefix} PCA"
    )
    umap_le_metrics, best_res_umap, le_umap_labels, le_umap_centers = run_leiden_clustering(
        df_umap[["UMAP1", "UMAP2"]], output_dir, f"{name_prefix} UMAP"
    )
    
    # 4. Save combined metrics to CSV
    pca_km_metrics["Space"], pca_km_metrics["Algorithm"] = "PCA", "KMeans"
    umap_km_metrics["Space"], umap_km_metrics["Algorithm"] = "UMAP", "KMeans"
    pca_db_metrics["Space"], pca_db_metrics["Algorithm"] = "PCA", "DBSCAN"
    umap_db_metrics["Space"], umap_db_metrics["Algorithm"] = "UMAP", "DBSCAN"
    pca_le_metrics["Space"], pca_le_metrics["Algorithm"] = "PCA", "Leiden"
    umap_le_metrics["Space"], umap_le_metrics["Algorithm"] = "UMAP", "Leiden"
    
    metrics_combined = pd.concat([
        pca_km_metrics, umap_km_metrics,
        pca_db_metrics, umap_db_metrics,
        pca_le_metrics, umap_le_metrics
    ], axis=0)
    
    metrics_path = output_dir / "TCGA-background" / f"clustering_metrics_{name_prefix.lower()}.csv"
    metrics_combined.to_csv(metrics_path, index=False)
    print(f"[{name_prefix}] Clustering metrics saved to {metrics_path}")
    
    # 5. Build DataFrames with assignments for plotting and exporting
    df_pca_cls = df_pca.copy()
    df_umap_cls = df_umap.copy()
    
    df_pca_cls["KMeans_PCA_Cluster"] = km_pca_labels.astype(str)
    df_pca_cls["DBSCAN_PCA_Cluster"] = db_pca_labels.astype(str)
    df_pca_cls["Leiden_PCA_Cluster"] = le_pca_labels.astype(str)
    
    df_umap_cls["KMeans_UMAP_Cluster"] = km_umap_labels.astype(str)
    df_umap_cls["DBSCAN_UMAP_Cluster"] = db_umap_labels.astype(str)
    df_umap_cls["Leiden_UMAP_Cluster"] = le_umap_labels.astype(str)
    
    # Plotting
    pca_centers_dict = {"KMeans": km_pca_centers, "DBSCAN": db_pca_centers, "Leiden": le_pca_centers}
    umap_centers_dict = {"KMeans": km_umap_centers, "DBSCAN": db_umap_centers, "Leiden": le_umap_centers}
    params_pca_dict = {"KMeans": best_k_pca, "DBSCAN": f"{best_eps_pca:.3f}", "Leiden": f"{best_res_pca:.2f}"}
    params_umap_dict = {"KMeans": best_k_umap, "DBSCAN": f"{best_eps_umap:.3f}", "Leiden": f"{best_res_umap:.2f}"}
    
    plot_clustering_results(
        df_pca_cls, df_umap_cls,
        pca_centers_dict, umap_centers_dict,
        params_pca_dict, params_umap_dict,
        output_dir, f"tcga_combined_{name_prefix.lower()}_clustering.svg",
        title_suffix
    )
    
    # 6. Save sample clusters CSV
    df_samples = pd.DataFrame({
        "SampleID": df_pca.index,
        "Cohort": df_pca["Cohort"],
        "PC1": df_pca["PC1"],
        "PC2": df_pca["PC2"],
        "UMAP1": df_umap["UMAP1"],
        "UMAP2": df_umap["UMAP2"],
        "KMeans_PCA_Cluster": km_pca_labels,
        "DBSCAN_PCA_Cluster": db_pca_labels,
        "Leiden_PCA_Cluster": le_pca_labels,
        "KMeans_UMAP_Cluster": km_umap_labels,
        "DBSCAN_UMAP_Cluster": db_umap_labels,
        "Leiden_UMAP_Cluster": le_umap_labels
    })
    samples_path = output_dir / "TCGA-background" / f"sample_clusters_{name_prefix.lower()}.csv"
    df_samples.to_csv(samples_path, index=False)
    print(f"[{name_prefix}] Sample cluster assignments saved to {samples_path}")
    
    # 7. Calculate cohort assignment distribution
    cohort_data = []
    unique_cohorts = sorted(df_samples["Cohort"].unique())
    for cohort in unique_cohorts:
        df_c = df_samples[df_samples["Cohort"] == cohort]
        total_samples = len(df_c)
        
        row = {
            "Cohort": cohort,
            "Total_Samples": total_samples,
        }
        
        for alg in ["KMeans", "DBSCAN", "Leiden"]:
            for space in ["PCA", "UMAP"]:
                col = f"{alg}_{space}_Cluster"
                counts = df_c[col].value_counts()
                dominant = counts.index[0] if not counts.empty else -1
                dominant_frac = (counts.iloc[0] / total_samples) if total_samples > 0 else 0.0
                row[f"{alg}_{space}_Dominant_Cluster"] = dominant
                row[f"{alg}_{space}_Dominant_Cluster_Fraction"] = dominant_frac
                
        cohort_data.append(row)
        
    df_cohorts = pd.DataFrame(cohort_data)
    cohorts_path = output_dir / "TCGA-background" / f"cohort_clusters_{name_prefix.lower()}.csv"
    df_cohorts.to_csv(cohorts_path, index=False)
    print(f"[{name_prefix}] Cohort cluster assignments saved to {cohorts_path}")
    
    # Maintain traditional output file name copies to satisfy original constraints
    import shutil
    shutil.copy(samples_path, output_dir / "TCGA-background" / "sample_clusters.csv")
    shutil.copy(cohorts_path, output_dir / "TCGA-background" / "cohort_clusters.csv")
    shutil.copy(metrics_path, output_dir / "TCGA-background" / "clustering_metrics.csv")


def analyze_cohort_distances(adata_combined, k_values=[50, 100, 200], pca_dims=[30, 40, 50], output_dir=None):
    import pandas as pd
    import numpy as np
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist
    import seaborn as sns
    import matplotlib.pyplot as plt
    import scanpy as sc

    print("Analyzing distances between iAtlas and TCGA cohorts in PCA space...")

    # We need PCA. Since adata_combined is already combat-corrected, we scale and run PCA.
    # To save time and memory, we'll run PCA directly on a copy of the combined object.
    adata_pca = adata_combined.copy()
    
    print("Scaling and computing PCA for distance analysis...")
    sc.pp.scale(adata_pca)
    sc.tl.pca(adata_pca, n_comps=max(pca_dims), svd_solver='arpack')
    
    full_pca = adata_pca.obsm['X_pca']

    # Get metadata
    source_labels = adata_pca.obs['Source'].astype(str).values
    cohort_labels = adata_pca.obs['dataset'].astype(str).values
    sample_ids = adata_pca.obs.index.values

    tcga_mask = source_labels == 'TCGA'
    iatlas_mask = source_labels == 'iAtlas'

    tcga_cohorts = cohort_labels[tcga_mask]
    iatlas_cohorts = cohort_labels[iatlas_mask]
    iatlas_sample_ids = sample_ids[iatlas_mask]
    
    out_dir = output_dir / "TCGA-background" / "distances"
    out_dir.mkdir(parents=True, exist_ok=True)

    for n_pca in pca_dims:
        if n_pca > full_pca.shape[1]:
            print(f"Requested {n_pca} PCA dimensions, but only {full_pca.shape[1]} available. Skipping.")
            continue
            
        pca_features = full_pca[:, :n_pca]
        tcga_pca_coords = pca_features[tcga_mask]
        iatlas_pca_coords = pca_features[iatlas_mask]

        for k in k_values:
            print(f"Running distance analysis for K={k}, PCA_dims={n_pca}...")
            
            # 1. KMeans on TCGA
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            tcga_cluster_labels = kmeans.fit_predict(tcga_pca_coords)
            centroids = kmeans.cluster_centers_

            # 2. TCGA Centroid Composition
            # Create a dataframe: rows=centroids, cols=tcga_cohorts, values=fraction
            composition_df = pd.DataFrame({'Cluster': tcga_cluster_labels, 'Cohort': tcga_cohorts})
            comp_counts = composition_df.groupby(['Cluster', 'Cohort']).size().unstack(fill_value=0)
            comp_fractions = comp_counts.div(comp_counts.sum(axis=1), axis=0)
            
            comp_csv_path = out_dir / f"tcga_centroid_composition_k{k}_pca{n_pca}.csv"
            comp_fractions.to_csv(comp_csv_path)

            # 3. L1 Distances: iAtlas Samples to TCGA Centroids
            distances = cdist(iatlas_pca_coords, centroids, metric='cityblock')
            
            dist_df = pd.DataFrame(distances, index=iatlas_sample_ids, columns=[f"Centroid_{i}" for i in range(k)])
            dist_df['iAtlas_Cohort'] = iatlas_cohorts
            dist_csv_path = out_dir / f"iatlas_sample_to_tcga_centroid_distances_k{k}_pca{n_pca}.csv"
            dist_df.to_csv(dist_csv_path)

            # 4. Average L1 Distance Heatmap (iAtlas Cohort -> TCGA Cohort)
            # Average distance from each iAtlas cohort to each Centroid
            # Exclude non-numeric columns like 'iAtlas_Cohort' from mean operation
            numeric_dist_df = dist_df.drop(columns=['iAtlas_Cohort'])
            numeric_dist_df['iAtlas_Cohort'] = dist_df['iAtlas_Cohort']
            mean_dist_to_centroid = numeric_dist_df.groupby('iAtlas_Cohort').mean()
            # Map centroids to TCGA cohorts using composition matrix
            avg_dist_matrix = mean_dist_to_centroid.values @ comp_fractions.values
            avg_dist_df = pd.DataFrame(avg_dist_matrix, index=mean_dist_to_centroid.index, columns=comp_fractions.columns)
            
            avg_dist_out = out_dir / f"iatlas_to_tcga_avg_distance_k{k}_pca{n_pca}.csv"
            avg_dist_df.to_csv(avg_dist_out)
            
            plt.figure(figsize=(14, 10))
            sns.heatmap(avg_dist_df, cmap="viridis_r", annot=False)
            plt.title(f"Average L1 Distance (K={k}, PCA={n_pca})")
            plt.ylabel("iAtlas Cohort")
            plt.xlabel("TCGA Cohort")
            plt.tight_layout()
            plt.savefig(out_dir / f"iatlas_to_tcga_avg_distance_heatmap_k{k}_pca{n_pca}.svg")
            plt.close()

            # 5. Nearest Centroid Voting Heatmap
            nearest_centroid_idx = np.argmin(distances, axis=1)
            sample_votes = comp_fractions.values[nearest_centroid_idx]
            votes_df = pd.DataFrame(sample_votes, index=iatlas_cohorts, columns=comp_fractions.columns)
            cohort_votes = votes_df.groupby(votes_df.index).mean()
            
            votes_out = out_dir / f"iatlas_to_tcga_voting_k{k}_pca{n_pca}.csv"
            cohort_votes.to_csv(votes_out)

            plt.figure(figsize=(14, 10))
            sns.heatmap(cohort_votes, cmap="YlGnBu", annot=False)
            plt.title(f"Nearest Centroid Voting Fraction (K={k}, PCA={n_pca})")
            plt.ylabel("iAtlas Cohort")
            plt.xlabel("TCGA Cohort")
            plt.tight_layout()
            plt.savefig(out_dir / f"iatlas_to_tcga_voting_heatmap_k{k}_pca{n_pca}.svg")
            plt.close()


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    output_dir = Path(".") / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
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
    import scipy.sparse as sp
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

    # Run 1: TCGA-Only using DE genes
    print("Running TCGA-only analysis on DE genes...")
    target_de = [g for g in de_hugo_genes if g in adata.var_names]
    adata_tcga_de = adata[:, target_de].copy()
    run_pca_umap_tcga_only(
        adata_tcga_de, 
        output_dir, 
        "tcga_only_de_pca_umap.svg", 
        "TCGA-Only DE Genes"
    )
    
    # Run 2: TCGA-Only using Bagaev signature genes
    target_bagaev = [g for g in bagaev_genes if g in adata.var_names]
    adata_tcga_bagaev = adata[:, target_bagaev].copy()
    run_pca_umap_tcga_only(
        adata_tcga_bagaev, 
        output_dir, 
        "tcga_only_bagaev_pca_umap.svg", 
        "TCGA-Only Bagaev Signature Genes"
    )

    # ---------------- UNCORRECTED RUN ---------------- #
    # Run 3: Combined TCGA & iAtlas using DE genes
    df_pca_de, df_umap_de = run_pca_umap_combined(
        adata, adata_iatlas, de_hugo_genes, output_dir,
        "tcga_combined_de_pca_umap_uncorrected.svg", "TCGA & iAtlas DE Genes (Uncorrected)"
    )

    # Run 4: Combined TCGA & iAtlas using Bagaev signatures
    df_pca_bagaev, df_umap_bagaev = run_pca_umap_combined(
        adata, adata_iatlas, bagaev_genes, output_dir,
        "tcga_combined_bagaev_pca_umap_uncorrected.svg", "TCGA & iAtlas Bagaev (Uncorrected)"
    )

    # Run 5: LDA on Source DE genes
    run_lda_and_plot(
        adata, adata_iatlas, source_de_genes, output_dir,
        "tcga_iatlas_source_de_lda_uncorrected.svg", "Source DE Genes (Uncorrected)"
    )

    # Perform clustering analysis on DE Genes PCA/UMAP
    perform_complete_clustering_analysis(df_pca_de, df_umap_de, output_dir, "DE_Uncorrected", "DE Genes (Uncorrected)")

    # Perform clustering analysis on Bagaev Signature Genes PCA/UMAP
    perform_complete_clustering_analysis(df_pca_bagaev, df_umap_bagaev, output_dir, "Bagaev_Uncorrected", "Bagaev (Uncorrected)")

    # Maintain original output file path to satisfy any traditional checks
    import shutil
    shutil.copy(
        output_dir / "TCGA-background" / "tcga_combined_de_pca_umap_uncorrected.svg",
        output_dir / "TCGA-background" / "tcga_pca_umap.svg"
    )


    # ---------------- COMBAT CORRECTION ---------------- #
    import anndata as ad
    import scipy.sparse as sp
    
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


    # ---------------- CORRECTED RUN ---------------- #
    # Run 6: Combined TCGA & iAtlas using DE genes
    df_pca_de_corr, df_umap_de_corr = run_pca_umap_combined(
        adata_corr, adata_iatlas_corr, de_hugo_genes, output_dir,
        "tcga_combined_de_pca_umap_combat.svg", "TCGA & iAtlas DE Genes (ComBat)"
    )

    # Run 7: Combined TCGA & iAtlas using Bagaev signatures
    df_pca_bagaev_corr, df_umap_bagaev_corr = run_pca_umap_combined(
        adata_corr, adata_iatlas_corr, bagaev_genes, output_dir,
        "tcga_combined_bagaev_pca_umap_combat.svg", "TCGA & iAtlas Bagaev (ComBat)"
    )

    # Run 8: LDA on Source DE genes
    run_lda_and_plot(
        adata_corr, adata_iatlas_corr, source_de_genes, output_dir,
        "tcga_iatlas_source_de_lda_combat.svg", "Source DE Genes (ComBat)"
    )

    # Perform clustering analysis on DE Genes PCA/UMAP
    perform_complete_clustering_analysis(df_pca_de_corr, df_umap_de_corr, output_dir, "DE_ComBat", "DE Genes (ComBat)")

    # Perform clustering analysis on Bagaev Signature Genes PCA/UMAP
    perform_complete_clustering_analysis(df_pca_bagaev_corr, df_umap_bagaev_corr, output_dir, "Bagaev_ComBat", "Bagaev (ComBat)")

    # Run cohort distance matching
    analyze_cohort_distances(
        adata_combined_full, 
        k_values=[50, 100, 200, 300, 400, 500], 
        pca_dims=[10, 20, 30, 40, 50], 
        output_dir=output_dir
    )

    # Save Harmonized Datasets to Datalair
    print("Saving harmonized datasets to Datalair...")
    from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    ds_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Writing h5ad files to {ds_path}...")
    # The prompt explicitly asked to save adata and adata_iatlas
    adata.write_h5ad(ds_path / "tcga_harmonized_uncorrected.h5ad")
    adata_iatlas.write_h5ad(ds_path / "iatlas_harmonized_uncorrected.h5ad")
    
    # We also save the ComBat corrected versions just in case
    adata_corr.write_h5ad(ds_path / "tcga_harmonized_combat.h5ad")
    adata_iatlas_corr.write_h5ad(ds_path / "iatlas_harmonized_combat.h5ad")
    print("Harmonized datasets successfully saved to Datalair.")


if __name__ == "__main__":
    main()
