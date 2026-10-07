"""Unit tests for batched gene normalization and incremental H5AD caching."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from returns.result import Success

from tme_datasets.preprocessing.gene_normalization import (
    batch_normalize_to_sparse_h5ad,
    normalize_dataset_to_ensembl,
)


def test_batch_normalize_equivalence_to_standard(tmp_path: Path) -> None:
    """Verify batch_normalize_to_sparse_h5ad produces identical results to in-memory normalization."""
    genes = [
        "CD8A", "CD4", "PDCD1", "CTLA4", "FOXP3", "TP53", "EGFR", "MYC", "IL2", "IFNG",
        "TNF", "GAPDH", "ACTB", "B2M", "PTEN", "VEGFA", "CD274", "UNMAPPED_GENE_123"
    ]
    n_cells = 300
    X = sp.random(n_cells, len(genes), density=0.25, format="csr", dtype=np.float32)
    obs = pd.DataFrame(
        {"cluster": ["C1" if i < 150 else "C2" for i in range(n_cells)]},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    # 1. In-memory reference normalization
    res_std = normalize_dataset_to_ensembl(adata)
    assert isinstance(res_std, Success)
    adata_std = res_std.unwrap()

    # 2. Batched normalization to H5AD with small batch_size
    target_h5ad = tmp_path / "normalized_batched.h5ad"
    res_batch = batch_normalize_to_sparse_h5ad(adata, target_h5ad, batch_size=75)
    assert isinstance(res_batch, Success)
    assert target_h5ad.exists()

    adata_batch = ad.read_h5ad(target_h5ad)

    # Verify shapes and indexes match exactly
    assert adata_batch.shape == adata_std.shape
    assert list(adata_batch.obs_names) == list(adata_std.obs_names)
    assert list(adata_batch.var_names) == list(adata_std.var_names)
    assert sp.isspmatrix_csr(adata_batch.X)

    # Verify expression values are identical within float precision
    diff = np.abs(adata_std.X.toarray() - adata_batch.X.toarray()).max()
    assert diff < 1e-6
