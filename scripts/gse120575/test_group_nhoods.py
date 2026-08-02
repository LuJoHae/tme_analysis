import anndata as ad
import pandas as pd
import scipy.sparse as sp
import milopy

adata = ad.read_h5ad('/storage/halu/data/GSE120575/gse120575_processed.h5ad', backed='r')
obs = adata.obs[adata.obs['treatment_status'] == 'Combined'].copy()
nhoods = sp.load_npz('/storage/halu/data/output/milopy_nhoods_Combined.npz')

# Create mock adata with nhoods
test_adata = ad.AnnData(X=sp.csr_matrix((nhoods.shape[0], 1)))
test_adata.obsm['nhoods'] = nhoods

res = pd.read_csv('/storage/halu/data/output/milopy_results_Combined.csv')
if 'Unnamed: 0' in res.columns: res = res.rename(columns={'Unnamed: 0': 'Nhood'})
elif '' in res.columns: res = res.rename(columns={'': 'Nhood'})

print('Running group_nhoods...')
grouped_res = milopy.group_nhoods(test_adata, res, max_fdr=0.1, resolution=1.0)
print('Done!')
print(grouped_res['NhoodGroup'].value_counts(dropna=False))
