import sys
import datalair
from pathlib import Path
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import norm
from single_cell_datasets._single_cell_datasets import SingleCellDataProcessStep07


class ReferenceSingleCell(datalair.Dataset):
    """Datalair Dataset class for all single cell datasets."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetSingleCell`."""
        super().__init__(namespace="DatasetSingleCell")


class SingleCellReference(ReferenceSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        ds = SingleCellDataProcessStep07()
        lair.get_dataset_filepaths(ds)
        adata = ad.read_h5ad(lair.get_dataset_filepaths(ds)["adata.h5ad"])
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))


# Lineage and subtype markers for hierarchical cell annotation
MAJOR_MARKERS = {
    'T_NK': ['CD3D', 'CD3E', 'CD8A', 'CD4', 'NKG7', 'NCAM1'],
    'B_Plasma': ['MS4A1', 'CD19', 'MZB1', 'SDC1'],
    'Myeloid': ['CD14', 'CD68', 'LYZ', 'CLEC9A', 'CD1C', 'LILRA4', 'FCGR3B', 'TPSAB1'],
    'Endothelial': ['PECAM1', 'VWF', 'PLVAP'],
    'Fibroblasts': ['COL1A1', 'COL3A1', 'DCN', 'LUM', 'ACTA2'],
    'Tumor_Epithelial': ['EPCAM', 'KRT18', 'KRT8', 'KRT19', 'PMEL', 'MLANA']
}

SUBTYPE_MARKERS = {
    'T_NK': {
        'CD8_Cytotoxic_T': ['CD8A', 'CD8B', 'GZMB', 'GZMK', 'PRF1'],
        'CD8_Exhausted_T': ['PDCD1', 'LAG3', 'HAVCR2', 'TOX'],
        'CD4_Helper_T': ['CD4', 'IL7R', 'CD40LG'],
        'Regulatory_T': ['FOXP3', 'IL2RA', 'BATF'],
        'Gamma_Delta_T': ['TRDV2', 'TRGV9'],
        'NK_cells': ['NCAM1', 'KLRD1', 'NKG7']
    },
    'B_Plasma': {
        'Naive_B': ['MS4A1', 'CD19', 'IGHD', 'TCL1A'],
        'Memory_B': ['MS4A1', 'CD27', 'CD83'],
        'Plasma_cells': ['MZB1', 'SDC1', 'IGHG1', 'IGHA1']
    },
    'Myeloid': {
        'Monocytes': ['CD14', 'FCGR3A', 'FCN1'],
        'M1_Macrophages': ['CD68', 'NOS2', 'TNF', 'IL1B'],
        'M2_Macrophages': ['CD68', 'CD163', 'MRC1', 'ARG1'],
        'cDC1': ['CLEC9A', 'XCR1'],
        'cDC2': ['CD1C', 'FCER1A'],
        'pDC': ['LILRA4', 'CLEC4C'],
        'Neutrophils': ['FCGR3B', 'CSF3R'],
        'Mast_cells': ['TPSAB1', 'TPSB2', 'MS4A2']
    },
    'Endothelial': {
        'Blood_Endothelial': ['PECAM1', 'VWF', 'PLVAP'],
        'Lymphatic_Endothelial': ['PROX1', 'PDPN', 'LYVE1']
    },
    'Fibroblasts': {
        'Normal_Fibroblasts': ['COL1A1', 'COL3A1', 'DCN', 'LUM'],
        'CAFs': ['ACTA2', 'TAGLN', 'FAP']
    },
    'Tumor_Epithelial': {
        'Epithelial_Tumor': ['EPCAM', 'KRT18', 'KRT8', 'KRT19'],
        'Melanoma_Tumor': ['PMEL', 'MLANA', 'MITF', 'TYR']
    }
}


def calculate_cnv_scores(adata_sub):
    """
    Calculate Chromosomal Expression Variance (CNV proxy) score.
    """
    chroms = [str(i) for i in range(1, 23)] + ['X']
    chrom_means = []
    
    for chrom in chroms:
        genes_on_chrom = adata_sub.var_names[adata_sub.var['contig'] == chrom]
        if len(genes_on_chrom) > 5:
            if isinstance(adata_sub.X, csr_matrix):
                mean_expr = np.asarray(adata_sub[:, genes_on_chrom].X.mean(axis=1)).flatten()
            else:
                mean_expr = np.asarray(adata_sub[:, genes_on_chrom].X.mean(axis=1)).flatten()
            chrom_means.append(mean_expr)
            
    if len(chrom_means) > 1:
        chrom_matrix = np.vstack(chrom_means)
        cnv_scores = np.var(chrom_matrix, axis=0)
    else:
        cnv_scores = np.zeros(adata_sub.n_obs)
        
    return cnv_scores


def annotate_clustering_resolution(adata_sub, resolution_col):
    """
    Perform hierarchical two-step annotation and tumor inference for a given clustering.
    """
    sc.tl.rank_genes_groups(adata_sub, groupby=resolution_col, method='wilcoxon')
    rank_res = adata_sub.uns['rank_genes_groups']
    
    major_annotations = {}
    subtype_annotations = {}
    tumor_classifications = {}
    
    tumor_markers = ['EPCAM', 'KRT18', 'KRT8', 'KRT19', 'PMEL', 'MLANA', 'MITF']
    valid_tumor_markers = [g for g in tumor_markers if g in adata_sub.var_names]
    if valid_tumor_markers:
        tumor_marker_matrix = adata_sub[:, valid_tumor_markers].X
        if isinstance(tumor_marker_matrix, csr_matrix):
            tumor_marker_scores = np.asarray(tumor_marker_matrix.mean(axis=1)).flatten()
        else:
            tumor_marker_scores = np.asarray(tumor_marker_matrix.mean(axis=1)).flatten()
    else:
        tumor_marker_scores = np.zeros(adata_sub.n_obs)
        
    cnv_threshold = np.percentile(adata_sub.obs['cnv_score'], 70)
    
    clusters = adata_sub.obs[resolution_col].unique()
    for cluster in clusters:
        cluster_mask = (adata_sub.obs[resolution_col] == cluster)
        top_markers = list(rank_res['names'][cluster][:30])
        
        # --- STEP 1: MAJOR CLASS INFERENCE (Dual-Method) ---
        major_scores = {}
        for major_type, genes in MAJOR_MARKERS.items():
            valid_g = [g for g in genes if g in adata_sub.var_names]
            if valid_g:
                expr_sub = adata_sub[cluster_mask, valid_g].X
                if isinstance(expr_sub, csr_matrix):
                    expr_score = float(expr_sub.mean())
                else:
                    expr_score = float(expr_sub.mean())
            else:
                expr_score = 0.0
                
            overlap_score = sum(1 for g in genes if g in top_markers) / len(genes)
            major_scores[major_type] = expr_score + overlap_score
            
        best_major = max(major_scores, key=major_scores.get)
        major_annotations[cluster] = best_major
        
        # --- STEP 2: SUBTYPE INFERENCE (Hierarchical) ---
        subtypes = SUBTYPE_MARKERS[best_major]
        subtype_scores = {}
        for sub_type, genes in subtypes.items():
            valid_g = [g for g in genes if g in adata_sub.var_names]
            if valid_g:
                expr_sub = adata_sub[cluster_mask, valid_g].X
                if isinstance(expr_sub, csr_matrix):
                    expr_score = float(expr_sub.mean())
                else:
                    expr_score = float(expr_sub.mean())
            else:
                expr_score = 0.0
                
            overlap_score = sum(1 for g in genes if g in top_markers) / len(genes)
            subtype_scores[sub_type] = expr_score + overlap_score
            
        best_subtype = max(subtype_scores, key=subtype_scores.get)
        subtype_annotations[cluster] = best_subtype
        
        # --- TUMOR CELL INFERENCE (Consensus) ---
        cluster_cnvs = adata_sub.obs.loc[cluster_mask, 'cnv_score']
        high_cnv_frac = np.mean(cluster_cnvs > cnv_threshold)
        
        cluster_marker_scores = tumor_marker_scores[cluster_mask]
        high_marker_frac = np.mean(cluster_marker_scores > 0.5)
        
        is_tumor_cluster = False
        if (high_cnv_frac > 0.5 and high_marker_frac > 0.5) or best_subtype in ['Epithelial_Tumor', 'Melanoma_Tumor']:
            is_tumor_cluster = True
            
        tumor_classifications[cluster] = "Tumor" if is_tumor_cluster else "Normal"
        
    return major_annotations, subtype_annotations, tumor_classifications


class SingleCellClustered(ReferenceSingleCell):
    """
    Datalair Dataset class for Clustered, Subsampled, and Automatically Annotated scRNA-seq reference.
    """
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        
        ds_ref = SingleCellReference()
        ref_paths = lair.get_dataset_filepaths(ds_ref)
        print(f"Loading SingleCellReference from: {ref_paths['adata.h5ad']}")
        adata = ad.read_h5ad(ref_paths["adata.h5ad"], backed="r")
        
        print("Performing stratified subsampling...")
        obs = adata.obs
        n_sample = 40000
        sampled_indices = []
        for cancer, df_sub in obs.groupby('cancer_code'):
            frac = len(df_sub) / len(obs)
            n_sub = int(np.round(frac * n_sample))
            if n_sub > 0:
                np.random.seed(42)
                chosen = np.random.choice(df_sub.index, size=min(n_sub, len(df_sub)), replace=False)
                sampled_indices.extend(chosen)
                
        sampled_indices = sorted(sampled_indices)
        
        try:
            adata_sub = adata[sampled_indices].to_memory()
        except Exception:
            adata_sub = adata[sampled_indices].copy()
            
        print("Mapping Ensembl IDs to unique gene symbols...")
        unique_symbols = []
        seen = set()
        for g in adata_sub.var['gene_name']:
            if pd.isna(g) or g == "":
                g = "Unknown"
            if g in seen:
                i = 1
                while f"{g}_{i}" in seen:
                    i += 1
                g = f"{g}_{i}"
            seen.add(g)
            unique_symbols.append(g)
        adata_sub.var_names = unique_symbols
        
        adata_sub.raw = adata_sub
        
        print("Running Scanpy library normalization and log1p...")
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
        
        print("Computing chromosomal variance CNV scores...")
        adata_sub.obs['cnv_score'] = calculate_cnv_scores(adata_sub)
        
        print("Selecting highly variable genes...")
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=2000, subset=True)
        sc.tl.pca(adata_sub, n_comps=40)
        sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata_sub)
        
        resolutions = [0.5, 1.0, 1.5, 2.0]
        print(f"Running Leiden clustering at resolutions: {resolutions}...")
        for res in resolutions:
            res_col = f"leiden_res_{res}"
            sc.tl.leiden(adata_sub, resolution=res, key_added=res_col)
            
            print(f"  Annotating clusters for resolution {res}...")
            major_map, subtype_map, tumor_map = annotate_clustering_resolution(adata_sub, res_col)
            
            adata_sub.obs[f"cell_type_major_res_{res}"] = adata_sub.obs[res_col].map(major_map)
            adata_sub.obs[f"cell_type_subtype_res_{res}"] = adata_sub.obs[res_col].map(subtype_map)
            adata_sub.obs[f"tumor_status_res_{res}"] = adata_sub.obs[res_col].map(tumor_map)
            
        print("Saving clustered AnnData object to Lair...")
        adata_sub.write_h5ad(output_dir.joinpath("adata.h5ad"))
        print(f"Dataset successfully derived and written to {output_dir}/adata.h5ad")


def calculate_witten_p_value(X, labels):
    """
    Calculate Daniela Witten's selective inference p-value for K-means split (k=2).
    Conditions on the partition boundary to correct for selection bias.
    """
    n1 = sum(labels == 0)
    n2 = sum(labels == 1)
    if n1 <= 1 or n2 <= 1:
        return 1.0
        
    X1 = X[labels == 0]
    X2 = X[labels == 1]
    
    c1 = X1.mean(axis=0)
    c2 = X2.mean(axis=0)
    
    v = c1 - c2
    v_norm = np.linalg.norm(v)
    if v_norm == 0:
        return 1.0
        
    v = v / v_norm
    
    # Project data onto separating axis
    z1 = X1 @ v
    z2 = X2 @ v
    
    m1 = z1.mean()
    m2 = z2.mean()
    D = abs(m1 - m2)
    
    # Estimate residual variance along the projection
    res_var = (np.sum((z1 - m1)**2) + np.sum((z2 - m2)**2)) / (n1 + n2 - 2)
    sigma_diff = np.sqrt(res_var * (1.0/n1 + 1.0/n2))
    if sigma_diff == 0:
        return 1.0
        
    # Truncation threshold: distance between the boundary-defining points
    if m1 > m2:
        t_trunc = max(0.0, z1.min() - z2.max())
    else:
        t_trunc = max(0.0, z2.min() - z1.max())
        
    # Truncated normal survival function (P(Z >= D | Z >= t_trunc))
    numerator = 1.0 - norm.cdf(D / sigma_diff)
    denominator = 1.0 - norm.cdf(t_trunc / sigma_diff)
    
    if denominator <= 1e-15:
        return 0.0
        
    pval = numerator / denominator
    return min(1.0, max(0.0, pval))


class SingleCellSubclustered(ReferenceSingleCell):
    """
    Datalair Dataset class for Subclustered single-cell reference.
    Runs cluster-specific PCA and K-means sub-clustering on all major cluster resolutions,
    evaluating clustering metrics and calculating Daniela Witten's selective inference p-values.
    """
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        
        # Load from SingleCellClustered dataset
        ds_clustered = SingleCellClustered()
        clustered_paths = lair.get_dataset_filepaths(ds_clustered)
        print(f"Loading SingleCellClustered from: {clustered_paths['adata.h5ad']}")
        adata = ad.read_h5ad(clustered_paths["adata.h5ad"])
        
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
        
        subclustering_records = []
        resolutions = [0.5, 1.0, 1.5, 2.0]
        
        for res in resolutions:
            res_col = f"leiden_res_{res}"
            subcluster_labels = pd.Series(index=adata.obs_names, dtype='string')
            best_k_labels = pd.Series(index=adata.obs_names, dtype='int32')
            
            clusters = adata.obs[res_col].unique()
            print(f"Processing sub-clustering for resolution {res} ({len(clusters)} clusters)...")
            
            for cluster in clusters:
                cluster_mask = (adata.obs[res_col] == cluster)
                n_cells = sum(cluster_mask)
                
                if n_cells < 10:
                    subcluster_labels.loc[cluster_mask] = f"{cluster}_sub_0"
                    best_k_labels.loc[cluster_mask] = 2
                    subclustering_records.append({
                        'Resolution': res,
                        'Cluster': cluster,
                        'n_cells': n_cells,
                        'Optimal_K': 2,
                        'Silhouette_Score': np.nan,
                        'Calinski_Harabasz_Index': np.nan,
                        'Davies_Bouldin_Index': np.nan,
                        'Witten_p_value': np.nan
                    })
                    continue
                
                # Perform cluster-specific PCA
                adata_c = adata[cluster_mask].copy()
                sc.pp.highly_variable_genes(adata_c, n_top_genes=min(500, adata_c.n_vars), subset=True)
                sc.tl.pca(adata_c, n_comps=min(10, adata_c.n_vars))
                X_pca = adata_c.obsm['X_pca']
                
                # Scan K from 2 to 5
                best_sil = -2.0
                best_k = 2
                best_labels = None
                best_ch = np.nan
                best_db = np.nan
                
                for k in [2, 3, 4, 5]:
                    if n_cells > k:
                        kmeans_k = KMeans(n_clusters=k, random_state=42, n_init=10)
                        labels_k = kmeans_k.fit_predict(X_pca)
                        sil = silhouette_score(X_pca, labels_k)
                        ch = calinski_harabasz_score(X_pca, labels_k)
                        db = davies_bouldin_score(X_pca, labels_k)
                        
                        if sil > best_sil:
                            best_sil = sil
                            best_k = k
                            best_labels = labels_k
                            best_ch = ch
                            best_db = db
                            
                # Calculate Witten's selective inference p-value for K=2 split
                kmeans_2 = KMeans(n_clusters=2, random_state=42, n_init=10)
                labels_2 = kmeans_2.fit_predict(X_pca)
                witten_p = calculate_witten_p_value(X_pca, labels_2)
                
                for idx, label in zip(adata.obs_names[cluster_mask], best_labels):
                    subcluster_labels.loc[idx] = f"{cluster}_sub_{label}"
                    best_k_labels.loc[idx] = best_k
                    
                subclustering_records.append({
                    'Resolution': res,
                    'Cluster': cluster,
                    'n_cells': n_cells,
                    'Optimal_K': best_k,
                    'Silhouette_Score': best_sil,
                    'Calinski_Harabasz_Index': best_ch,
                    'Davies_Bouldin_Index': best_db,
                    'Witten_p_value': witten_p
                })
                
            adata.obs[f"kmeans_subcluster_res_{res}"] = subcluster_labels.astype('category')
            adata.obs[f"kmeans_best_k_res_{res}"] = best_k_labels.astype('category')
            
        # Save metrics CSV
        df_metrics = pd.DataFrame(subclustering_records)
        metrics_csv_path = output_dir / "subclustering_metrics.csv"
        df_metrics.to_csv(metrics_csv_path, index=False)
        print(f"Saved sub-clustering evaluation metrics to {metrics_csv_path.name}")
        
        # Save clustered AnnData
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))
        print(f"Sub-clustered AnnData successfully written to {output_dir}/adata.h5ad")


class SingleCellClusterMeans(ReferenceSingleCell):
    """
    Datalair Dataset class for Single-Cell Cluster Means.
    Computes average expression profiles for each cluster across all Leiden resolutions and K-means sub-clusterings.
    Exports separate h5ad files for each clustering.
    """
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        
        # Load from SingleCellSubclustered dataset
        ds_sub = SingleCellSubclustered()
        sub_paths = lair.get_dataset_filepaths(ds_sub)
        print(f"Loading SingleCellSubclustered from: {sub_paths['adata.h5ad']}")
        adata = ad.read_h5ad(sub_paths["adata.h5ad"])
        
        clusterings = [
            'leiden_res_0.5', 'leiden_res_1.0', 'leiden_res_1.5', 'leiden_res_2.0',
            'kmeans_subcluster_res_0.5', 'kmeans_subcluster_res_1.0', 'kmeans_subcluster_res_1.5', 'kmeans_subcluster_res_2.0'
        ]
        
        # Determine whether to use raw expression or log-normalized counts
        use_raw = adata.raw is not None
        if use_raw:
            expr_matrix = adata.raw.X
            var_df = adata.raw.var
        else:
            expr_matrix = adata.X
            var_df = adata.var
            
        for clustering_col in clusterings:
            print(f"Averaging cluster expression profiles for: {clustering_col}...")
            # Drop cells with missing labels in this clustering
            obs_valid = adata.obs.dropna(subset=[clustering_col])
            clusters = sorted(obs_valid[clustering_col].unique())
            
            means = []
            cluster_obs = []
            
            for cluster in clusters:
                mask = (adata.obs[clustering_col] == cluster)
                sub_matrix = expr_matrix[mask.values]
                
                # Calculate mean expression profile
                if isinstance(sub_matrix, csr_matrix):
                    cluster_mean = np.asarray(sub_matrix.mean(axis=0)).flatten()
                else:
                    cluster_mean = np.asarray(sub_matrix.mean(axis=0)).flatten()
                means.append(cluster_mean)
                
                # Consolidate cluster observations (majority vote for annotations)
                cluster_subset = adata.obs[mask]
                mode_metadata = {}
                for col in adata.obs.columns:
                    if col != clustering_col:
                        # Skip other clustering columns to avoid confusion
                        if any(c in col for c in ['leiden_res', 'kmeans_subcluster_res']):
                            continue
                        mode_val = cluster_subset[col].mode()
                        if len(mode_val) > 0:
                            mode_metadata[col] = mode_val[0]
                        else:
                            mode_metadata[col] = np.nan
                            
                mode_metadata['cluster_label'] = cluster
                cluster_obs.append(mode_metadata)
                
            # Create AnnData for this clustering's means
            X_means = np.vstack(means)
            df_obs = pd.DataFrame(cluster_obs, index=clusters)
            
            adata_means = ad.AnnData(X=X_means, obs=df_obs, var=var_df.copy())
            
            # Save to separate file
            file_name = f"means_{clustering_col}.h5ad"
            adata_means.write_h5ad(output_dir.joinpath(file_name))
            print(f"  Exported {file_name} with shape {adata_means.shape}")


class SingleCellDeconvolution(ReferenceSingleCell):
    """
    Datalair Dataset class for Single-Cell iAtlas Deconvolution.
    Deconvolutes bulk expression from iAtlas datasets using cluster mean profiles from 8 clusterings as reference.
    Uses instaprism with 50 iterations.
    """
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Load SingleCellClusterMeans files
        ds_means = SingleCellClusterMeans()
        means_paths = lair.get_dataset_filepaths(ds_means)
        
        import instaprism
        import ici_datasets
        from pyensembl import EnsemblRelease
        from joblib import Parallel, delayed
        
        # Ensembl Release for TPM transformation if needed
        data_source = EnsemblRelease(110)
        
        def extract_kb_span(symbol: str) -> float | None:
            try:
                genes = data_source.genes_by_name(symbol)
                return (genes[0].end - genes[0].start + 1) / 1e3 if genes else None
            except Exception:
                return None

        def transform_expr_to_tpm(counts_df: pd.DataFrame) -> pd.DataFrame:
            gene_lengths = pd.Series({
                symbol: length
                for symbol in counts_df.index
                if (length := extract_kb_span(symbol)) is not None
            })
            common_genes = counts_df.index.intersection(gene_lengths.index)
            counts = counts_df.loc[common_genes]
            lengths = gene_lengths.loc[common_genes]
            rpk = counts.div(lengths, axis=0)
            return rpk.div(rpk.sum(axis=0), axis=1) * 1e6

        def load_bulk_data(cohort_name: str) -> pd.DataFrame:
            dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
            ds = dataset_class(name=cohort_name)
            lair.safe_derive(ds)
            
            filepaths = lair.get_dataset_filepaths(ds)
            unpacked_key = next(f for f in filepaths.keys() if not f.endswith('.tar.gz'))
            p_dir = filepaths[unpacked_key] / filepaths[unpacked_key].name
            
            dispatch_map = {
                "data_mrna_seq_expression.txt": (transform_expr_to_tpm, "none"),
                "data_mrna_seq_tpm.txt": (lambda df: df, "tpm"),
                "data_mrna_seq_rpkm.txt": (lambda df: df.div(df.sum(axis=0), axis=1) * 1e6, "rpkm"),
            }
            existing = [f for f in dispatch_map if (p_dir / f).is_file()]
            if not existing:
                raise FileNotFoundError(f"No bulk mRNA expression file found in {p_dir}")
            target = existing[0]
            df = pd.read_csv(p_dir / target, sep='\t', index_col=0)
            transform_fn, _ = dispatch_map[target]
            return transform_fn(df)
            
        # Target cohorts to deconvolute
        cohorts = [
            "Hugo-iAtlas",
            "Riaz-iAtlas",
            "Liu-iAtlas",
            "Gide-iAtlas",
            "Rosenberg-iAtlas",
            "Padron-iAtlas",
            "Anders-iAtlas",
            "McDermott-iAtlas",
            "Choueiri-iAtlas"
        ]
        
        tcga_projects = ["SKCM", "BLCA", "PAAD", "BRCA", "KIRC"]
        
        # Load all valid target expression datasets first to avoid reloading
        bulk_datasets = {}
        for cohort in cohorts:
            try:
                df = load_bulk_data(cohort)
                # Convert index to uppercase to normalize gene names
                df.index = df.index.str.upper()
                # Aggregate duplicate symbols by averaging expression
                df = df.groupby(level=0).mean()
                bulk_datasets[cohort] = df
                print(f"  Successfully loaded {cohort} bulk expression with shape {df.shape}")
            except Exception as e:
                print(f"  Could not load expression for {cohort}: {e}")
                
        # Load and preprocess TCGA cohorts
        print("  Loading TCGA datasets...")
        import tcga
        from pyensembl import EnsemblRelease
        ensembl = EnsemblRelease(release=111, species="human")
        try:
            ensembl.download()
            ensembl.index()
        except Exception:
            pass
            
        ds_tcga = tcga.AllProjectsAdata()
        paths_tcga = lair.get_dataset_filepaths(ds_tcga)
        
        from ml_pipelines.tcga_background import aggregate_and_subset_by_hugo
        
        for project in tcga_projects:
            try:
                file_key = f"{project}.h5ad"
                adata_tcga = ad.read_h5ad(paths_tcga[file_key]).T
                adata_tcga.var_names_make_unique()
                
                # Map Ensembl IDs to Hugo symbols using pyensembl
                gene_id_to_hugo = {}
                valid_ids = []
                for gene_id in adata_tcga.var_names:
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
                    except Exception:
                        pass
                        
                adata_tcga = adata_tcga[:, valid_ids].copy()
                adata_tcga.var['hugo_symbol'] = [gene_id_to_hugo[g] for g in adata_tcga.var_names]
                
                valid_hugo = list(set([h for h in adata_tcga.var['hugo_symbol'] if not h.startswith("ENSG")]))
                adata_agg = aggregate_and_subset_by_hugo(adata_tcga, valid_hugo)
                
                df_tcga = pd.DataFrame(
                    adata_agg.X.T,
                    index=adata_agg.var_names,
                    columns=adata_agg.obs_names
                )
                df_tcga.index = df_tcga.index.str.upper()
                df_tcga = df_tcga.groupby(level=0).mean()
                
                cohort_key = f"TCGA-{project}"
                bulk_datasets[cohort_key] = df_tcga
                print(f"  Successfully loaded {cohort_key} bulk expression with shape {df_tcga.shape}")
            except Exception as e:
                print(f"  Could not load/process TCGA-{project}: {e}")
                
        # Loop over each clustering means file
        for file_key, file_path in sorted(list(means_paths.items())):
            if file_key.startswith("means_"):
                clustering_name = file_key[6:].replace(".h5ad", "")
            else:
                clustering_name = file_key.replace(".h5ad", "")
            print(f"\nDeconvoluting with reference clustering: {clustering_name}...")
            
            # Load cluster means AnnData
            adata_mean = ad.read_h5ad(file_path)
            # Ensure reference gene symbols are uppercase
            adata_mean.var_names = adata_mean.var_names.str.upper()
            # Filter out genes that have zero expression across all reference clusters
            ref_sum = adata_mean.X.sum(axis=0)
            if hasattr(ref_sum, "A1"):
                ref_sum = ref_sum.A1
            valid_ref_genes = adata_mean.var_names[ref_sum > 0]
            adata_mean = adata_mean[:, valid_ref_genes].copy()
            
            for cohort, bulk_df in bulk_datasets.items():
                common_genes = bulk_df.index.intersection(adata_mean.var_names)
                if len(common_genes) == 0:
                    print(f"  Warning: No common genes between bulk {cohort} and reference {clustering_name}. Skipping.")
                    continue
                    
                print(f"  Deconvoluting {cohort} using {clustering_name} (common genes: {len(common_genes)})...")
                
                bulk_aligned = bulk_df.loc[common_genes].copy()
                ref_aligned = adata_mean[:, common_genes].copy()
                
                reference_matrix = ref_aligned.X
                if hasattr(reference_matrix, "toarray"):
                    reference_matrix = reference_matrix.toarray()
                # Normalize each cell type profile to sum to 1 across the common genes
                reference_matrix = reference_matrix / reference_matrix.sum(axis=1, keepdims=True)
                
                cell_types = adata_mean.obs['cluster_label'].values
                sample_ids = bulk_df.columns
                
                # Run instaprism deconvolution
                def deconvolute_sample(sample_id):
                    bulk_sample = bulk_aligned[sample_id].values
                    _, _, cell_fracs, _ = instaprism.insta_prism(
                        bulk=bulk_sample,
                        reference=reference_matrix,
                        n_iter=50
                    )
                    return cell_fracs
                    
                # Run in parallel using joblib
                results = Parallel(n_jobs=-1)(
                    delayed(deconvolute_sample)(sid) for sid in sample_ids
                )
                
                # Save resulting cell fractions
                deconv_df = pd.DataFrame(results, index=sample_ids, columns=cell_types)
                out_csv = output_dir / f"deconv_{cohort}_{clustering_name}.csv"
                deconv_df.to_csv(out_csv)
                print(f"  Saved cell fractions to {out_csv.name}")