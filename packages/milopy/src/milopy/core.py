import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata

def build_graph(adata, k=30, d=30, **kwargs):
    """
    Builds a kNN graph. Wraps scanpy.pp.neighbors.
    """
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA must be computed first. Run sc.pp.pca(adata)")
        
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=d, **kwargs)
    return adata

def make_nhoods(adata, prop=0.1, k=30, d=30, refined=True, random_state=None):
    """
    Defines neighbourhoods on the kNN graph.
    """
    if 'distances' not in adata.obsp:
        raise ValueError("kNN graph not found. Run build_graph(adata)")
        
    n_cells = adata.n_obs
    n_nhoods = int(np.round(n_cells * prop))
    
    if random_state is not None:
        np.random.seed(random_state)
        
    # Initial random sampling
    vertex_indices = np.random.choice(n_cells, size=n_nhoods, replace=False)
    
    knn_graph = adata.obsp['connectivities']
    
    if refined:
        # Refinement step: for each random vertex, find its neighbourhood,
        # compute the median profile in PCA space, and pick the closest vertex.
        pca_coords = adata.obsm['X_pca'][:, :d]
        refined_vertices = []
        for v in vertex_indices:
            # Find neighbors of v
            neighbors = knn_graph[v].nonzero()[1]
            # include self
            neighborhood = np.append(neighbors, v)
            
            # Compute median profile
            median_profile = np.median(pca_coords[neighborhood], axis=0)
            
            # Find vertex in neighborhood closest to median
            distances = np.linalg.norm(pca_coords[neighborhood] - median_profile, axis=1)
            closest_idx = neighborhood[np.argmin(distances)]
            refined_vertices.append(closest_idx)
            
        # Remove duplicates
        vertex_indices = np.unique(refined_vertices)
        
    # Store neighborhoods as a sparse matrix: cells x nhoods
    rows = []
    cols = []
    
    for i, v in enumerate(vertex_indices):
        neighbors = knn_graph[v].nonzero()[1]
        neighborhood = np.append(neighbors, v)
        rows.extend(neighborhood)
        cols.extend([i] * len(neighborhood))
        
    nhoods_mat = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cells, len(vertex_indices))).tocsc()
    
    adata.obsm['nhoods'] = nhoods_mat
    adata.uns['nhood_indices'] = vertex_indices
    return adata

def count_cells(adata, sample_col):
    """
    Counts cells in each neighbourhood across samples.
    """
    if 'nhoods' not in adata.obsm:
        raise ValueError("Neighborhoods not found. Run make_nhoods(adata)")
        
    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in adata.obs")
        
    samples = adata.obs[sample_col].astype('category')
    sample_categories = samples.cat.categories
    sample_codes = samples.cat.codes.values
    
    n_nhoods = adata.obsm['nhoods'].shape[1]
    n_samples = len(sample_categories)
    
    # Initialize count matrix
    counts = np.zeros((n_nhoods, n_samples))
    
    nhoods_mat = adata.obsm['nhoods'].tocsc()
    
    for i in range(n_nhoods):
        cells_in_nhood = nhoods_mat[:, i].nonzero()[0]
        # Count cells per sample in this neighborhood
        nhood_sample_codes = sample_codes[cells_in_nhood]
        counts[i, :] = np.bincount(nhood_sample_codes, minlength=n_samples)
        
    count_df = pd.DataFrame(counts, columns=sample_categories)
    
    # store in uns
    adata.uns['nhood_counts'] = count_df
    return adata

def calc_nhood_distance(adata, d=30):
    """
    Calculates distances between neighbourhoods based on overlap or PCA distance.
    """
    if 'nhoods' not in adata.obsm:
        raise ValueError("Neighborhoods not found. Run make_nhoods(adata)")
        
    nhoods_mat = adata.obsm['nhoods'].tocsc()
    pca_coords = adata.obsm['X_pca'][:, :d]
    
    n_nhoods = nhoods_mat.shape[1]
    
    # compute median of each neighborhood
    nhood_medians = np.zeros((n_nhoods, d))
    for i in range(n_nhoods):
        cells = nhoods_mat[:, i].nonzero()[0]
        nhood_medians[i, :] = np.median(pca_coords[cells, :], axis=0)
        
    # compute euclidean distance between medians
    from scipy.spatial.distance import pdist, squareform
    dists = pdist(nhood_medians, metric='euclidean')
    dist_mat = squareform(dists)
    
    adata.uns['nhood_distances'] = dist_mat
    return adata

def test_nhoods(adata, design, design_df, model_contrasts=None):
    """
    Tests for differential abundance using edgepython (pure Python edgeR port).
    `design` is a formula string like '~ condition'
    `design_df` is a pandas DataFrame with row names matching sample columns in `nhood_counts`
    """
    if 'nhood_counts' not in adata.uns:
        raise ValueError("Neighborhood counts not found. Run count_cells(adata, sample_col)")
        
    import edgepython as ep
    import patsy
    
    # Counts are nhoods x samples
    counts_df = adata.uns['nhood_counts']
    
    # Ensure design_df is aligned with counts_df columns
    design_df = design_df.loc[counts_df.columns]
    
    # edgepython expects genes x samples (rows=features, cols=samples)
    # Our counts_df is nhoods x samples, which is already the right orientation
    counts_matrix = counts_df.values.astype(float)
    
    # Create DGEList
    y = ep.make_dgelist(counts=counts_matrix)
    y = ep.calc_norm_factors(y, method="TMM")
    
    # Create design matrix using patsy
    design_mat = patsy.dmatrix(design, design_df, return_type='dataframe')
    design_array = np.asarray(design_mat, dtype=float)
    
    # Estimate dispersion
    y = ep.estimate_disp(y, design=design_array)
    
    # Fit QL GLM
    fit = ep.glm_ql_fit(y, design=design_array, robust=True)
    
    # Test
    if model_contrasts is not None:
        # Build contrast vector from string like "conditionB - conditionA"
        # For now, support simple column index contrasts
        res = ep.glm_ql_ftest(fit, contrast=model_contrasts)
    else:
        # Default to testing the last coefficient
        n_coefs = design_array.shape[1]
        res = ep.glm_ql_ftest(fit, coef=n_coefs - 1)
        
    # Get results table
    top = ep.top_tags(res, n=counts_matrix.shape[0], sort_by="none")
    res_df = top['table']
    
    # Ensure it's a DataFrame with the right index
    if not isinstance(res_df, pd.DataFrame):
        res_df = pd.DataFrame(res_df)
    
    res_df.index = counts_df.index
    adata.uns['nhood_test_results'] = res_df
    return adata


def group_nhoods(adata, da_res, max_fdr=0.1, overlap_threshold=0.0, resolution=0.05):
    """
    Groups significant neighborhoods into broader cell continuia (Milo Modules).
    Mimics miloR::groupNhoods.
    
    Parameters:
    adata: AnnData object containing neighborhood data.
    da_res: pd.DataFrame of differential abundance results (must have 'FDR' and 'logFC' columns).
    max_fdr: float, maximum FDR threshold to consider a neighborhood significant.
    overlap_threshold: float, minimum overlap fraction to draw an edge.
    resolution: float, resolution parameter for Leiden clustering.
    
    Returns:
    pd.DataFrame: da_res with an added 'NhoodGroup' column.
    """
    if 'nhoods' not in adata.obsm:
        raise ValueError("Neighborhoods not found. Run make_nhoods(adata)")
        
    nhoods_mat = adata.obsm['nhoods']
    
    # Ensure da_res is aligned with nhoods
    if len(da_res) != nhoods_mat.shape[1]:
        raise ValueError("Length of da_res does not match number of neighborhoods.")
        
    # Find significant neighborhoods
    sig_idx = np.where(da_res['FDR'].fillna(1.0) < max_fdr)[0]
    
    if len(sig_idx) == 0:
        print("No significant neighborhoods found.")
        da_res['NhoodGroup'] = np.nan
        return da_res
        
    # Calculate overlap matrix: nhoods x nhoods
    sig_nhoods_mat = nhoods_mat[:, sig_idx]
    overlap_mat = sig_nhoods_mat.T.dot(sig_nhoods_mat)
    
    # Calculate neighborhood sizes
    nhood_sizes = np.array(sig_nhoods_mat.sum(axis=0)).flatten()
    
    # Calculate overlap fraction (intersection / smaller size)
    overlap_mat_dense = overlap_mat.toarray()
    min_size_mat = np.minimum(nhood_sizes[:, None], nhood_sizes[None, :])
    
    # Avoid division by zero
    min_size_mat[min_size_mat == 0] = 1
    overlap_frac = overlap_mat_dense / min_size_mat
    
    # Filter edges between discordant logFC signs
    logfc = da_res['logFC'].values[sig_idx]
    signs = np.sign(logfc)
    sign_matrix = np.outer(signs, signs)
    valid_edges = sign_matrix > 0
    
    # Apply mask to the overlap fraction matrix
    overlap_frac[~valid_edges] = 0
    overlap_frac[overlap_frac < overlap_threshold] = 0
    
    # Zero out diagonal to prevent self-loops from affecting clustering modularity
    np.fill_diagonal(overlap_frac, 0)
    
    # Convert back to sparse
    sig_overlap = sp.csr_matrix(overlap_frac)
    
    # Create dummy AnnData for clustering
    import anndata as ad
    import scanpy as sc
    
    sub_adata = ad.AnnData(X=np.zeros((len(sig_idx), 1)))
    sub_adata.obsp['connectivities'] = sig_overlap
    
    # Run leiden
    sc.tl.leiden(sub_adata, resolution=resolution, adjacency=sig_overlap, flavor='igraph', n_iterations=2, directed=False)
    
    # Map back to da_res
    group_col = np.full(len(da_res), np.nan, dtype=object)
    group_col[sig_idx] = sub_adata.obs['leiden'].values
    
    # Ensure categoricals
    da_res_out = da_res.copy()
    da_res_out['NhoodGroup'] = pd.Categorical(group_col)
    
    return da_res_out

