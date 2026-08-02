import scanpy as sc
import pertpy as pt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp

def main():
    print("1. Loading real single-cell data (GSE120575)...")
    adata = sc.read_h5ad('/storage/halu/data/GSE120575/gse120575_processed.h5ad')
    
    # Map columns to what our script expects
    condition_col = 'characteristics: response'
    sample_col = 'patient_id'
    
    # Filter out any NA responses just in case
    adata = adata[~adata.obs[condition_col].isna()].copy()
    
    print("2. Building kNN graph...")
    # Compute neighbors (using X_pca if available, otherwise it will compute it or use X)
    if 'X_pca' not in adata.obsm:
        print("X_pca not found, computing PCA...")
        sc.pp.pca(adata)
        
    sc.pp.neighbors(adata, n_neighbors=30)
    
    print("3. Defining Milo neighborhoods manually (random sampling)...")
    n_cells = adata.n_obs
    n_index = int(n_cells * 0.1) # 10% of real cells
    np.random.seed(42)
    index_cell_indices = np.random.choice(n_cells, n_index, replace=False)
    
    knn_graph = adata.obsp["connectivities"]
    full_graph = knn_graph + sp.eye(n_cells)
    W = full_graph[:, index_cell_indices]
    W = (W > 0).astype(int)
    W_csc = W.tocsc()
    
    results = []
    
    print("4. Calculating responder and sample fractions per neighborhood...")
    condition_array = adata.obs[condition_col].values
    sample_array = adata.obs[sample_col].values
    
    for j in range(W.shape[1]):
        index_cell_idx = index_cell_indices[j]
        index_condition = condition_array[index_cell_idx]
        index_sample = sample_array[index_cell_idx]
        
        cells_in_nhood_idx = W_csc[:, j].indices
        
        # EXCLUDE the anchor cell itself from the neighborhood
        cells_in_nhood_idx = cells_in_nhood_idx[cells_in_nhood_idx != index_cell_idx]
        
        # Skip if neighborhood is empty without the anchor cell
        if len(cells_in_nhood_idx) == 0:
            continue
        
        conditions_in_nhood = condition_array[cells_in_nhood_idx]
        responder_fraction = np.sum(conditions_in_nhood == 'Responder') / len(cells_in_nhood_idx)
        
        samples_in_nhood = sample_array[cells_in_nhood_idx]
        same_sample_fraction = np.sum(samples_in_nhood == index_sample) / len(cells_in_nhood_idx)
        
        # New calculation: exclude cells from the same patient
        diff_patient_mask = samples_in_nhood != index_sample
        diff_patient_cells = conditions_in_nhood[diff_patient_mask]
        
        if len(diff_patient_cells) > 0:
            diff_patient_responder_fraction = np.sum(diff_patient_cells == 'Responder') / len(diff_patient_cells)
        else:
            diff_patient_responder_fraction = np.nan
        
        results.append({
            'nhood_id': j,
            'index_cell_condition': index_condition,
            'responder_fraction': responder_fraction,
            'same_sample_fraction': same_sample_fraction,
            'diff_patient_responder_fraction': diff_patient_responder_fraction
        })
        
    df_results = pd.DataFrame(results)
    
    print("\n5. Generating subplots...")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    sns.histplot(data=df_results, x='responder_fraction', hue='index_cell_condition', 
                 bins=15, multiple="layer", alpha=0.5, stat="probability", common_norm=False,
                 palette={'Responder': 'red', 'Non-responder': 'blue'}, ax=axes[0])
    axes[0].set_title('Responder Cells within Neighborhood\nStratified by Index Cell Condition')
    axes[0].set_xlabel('Fraction of Responder Cells')
    axes[0].set_ylabel('Probability Density')
    
    sns.histplot(data=df_results, x='same_sample_fraction', hue='index_cell_condition', 
                 bins=15, multiple="layer", alpha=0.5, stat="probability", common_norm=False,
                 palette={'Responder': 'red', 'Non-responder': 'blue'}, ax=axes[1])
    axes[1].set_title('Cells from the Same Patient\nas the Index Cell')
    axes[1].set_xlabel('Fraction of Cells from Same Patient')
    axes[1].set_ylabel('Probability Density')
    
    df_clean = df_results.dropna(subset=['diff_patient_responder_fraction'])
    sns.histplot(data=df_clean, x='diff_patient_responder_fraction', hue='index_cell_condition', 
                 bins=15, multiple="layer", alpha=0.5, stat="probability", common_norm=False,
                 palette={'Responder': 'red', 'Non-responder': 'blue'}, ax=axes[2])
    axes[2].set_title('Responder Cells (Excluding Same Patient)\nStratified by Index Cell Condition')
    axes[2].set_xlabel('Fraction of Responder Cells (Diff. Patient)')
    axes[2].set_ylabel('Probability Density')
    
    plt.tight_layout()
    
    plot_path = "neighborhood_fractions_plot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully to {plot_path}")

if __name__ == "__main__":
    main()
