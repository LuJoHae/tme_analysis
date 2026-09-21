"""Pytest fixtures for tme_datasets test suite."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def mock_single_cell_adata() -> ad.AnnData:
    """Generate synthetic single-cell AnnData for testing."""
    rng = np.random.default_rng(42)
    n_cells = 100
    n_genes = 20
    genes = [f"GENE_{i}" for i in range(n_genes)]
    # Include some key TME genes
    genes[0] = "CD8A"
    genes[1] = "CD4"
    genes[2] = "CD68"
    genes[3] = "PDCD1"

    counts = rng.negative_binomial(n=5, p=0.3, size=(n_cells, n_genes)).astype(np.float32)
    cell_types = rng.choice(["CD8_T", "CD4_T", "Macrophage", "B_cell"], size=n_cells)
    response = rng.choice([1.0, 0.0], size=n_cells, p=[0.4, 0.6])
    batch = rng.choice(["Cohort_A", "Cohort_B"], size=n_cells)

    obs = pd.DataFrame(
        {
            "cell_type": cell_types,
            "response_binary": response,
            "batch": batch,
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )

    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10))
    return adata
