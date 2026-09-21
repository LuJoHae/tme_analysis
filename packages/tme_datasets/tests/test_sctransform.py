"""Tests for SCTransform (Analytic Pearson Residuals and Regularized GLM)."""

import anndata as ad
import numpy as np
import pytest
from returns.maybe import Some
from returns.result import Success

from tme_datasets.models import SCTransformConfig
from tme_datasets.preprocessing.sctransform import normalize_sctransform
from tme_datasets.types import SCTransformFlavor


@pytest.fixture
def synthetic_adata() -> ad.AnnData:
    """Generate synthetic single-cell counts with varying cell library sizes."""
    rng = np.random.default_rng(42)
    N_cells, G_genes = 100, 40

    # Depth variation across cells
    depth_multipliers = rng.uniform(0.5, 3.0, size=(N_cells, 1))
    base_rates = rng.uniform(0.5, 20.0, size=(1, G_genes))

    lam = depth_multipliers * base_rates
    X = rng.poisson(lam=lam).astype(np.float32)

    a = ad.AnnData(X=X)
    a.obs_names = [f"Cell_{i:03d}" for i in range(N_cells)]
    a.var_names = [f"Gene_{i:03d}" for i in range(G_genes)]
    return a


def test_analytic_pearson_residuals(synthetic_adata):
    """Verify that Analytic Pearson Residuals (Lause et al.) executes and computes residuals."""
    config = SCTransformConfig(
        flavor=SCTransformFlavor.ANALYTIC,
        n_top_genes=Some(20),
        clip_residuals=True,
    )

    res = normalize_sctransform(synthetic_adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # Raw counts preserved
    assert "raw_counts" in out_adata.layers
    # Pearson residuals stored in layers and in X
    assert "pearson_residuals" in out_adata.layers
    assert out_adata.X.shape == synthetic_adata.shape

    # Residuals should have mean close to 0
    residuals = out_adata.layers["pearson_residuals"]
    assert np.abs(np.mean(residuals)) < 0.2

    # Residuals are bounded by sqrt(N_cells)
    max_bound = np.sqrt(synthetic_adata.n_obs) + 0.1
    assert np.max(np.abs(residuals)) <= max_bound

    # Highly variable genes selected
    assert "highly_variable" in out_adata.var
    assert np.sum(out_adata.var["highly_variable"]) == 20


def test_regularized_glm_residuals(synthetic_adata):
    """Verify that Regularized GLM (Hafemeister & Satija) smooths parameters and computes residuals."""
    config = SCTransformConfig(
        flavor=SCTransformFlavor.REGULARIZED_GLM,
        clip_residuals=True,
    )

    res = normalize_sctransform(synthetic_adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    assert "raw_counts" in out_adata.layers
    assert "pearson_residuals" in out_adata.layers
    assert "sct_beta0" in out_adata.var
    assert "sct_beta1" in out_adata.var

    residuals = out_adata.layers["pearson_residuals"]
    assert residuals.shape == synthetic_adata.shape
    # Bounded residuals
    assert np.max(np.abs(residuals)) <= np.sqrt(synthetic_adata.n_obs) + 0.1


def test_immutability(synthetic_adata):
    """Verify that input AnnData is not mutated in place."""
    orig_x = synthetic_adata.X.copy()
    config = SCTransformConfig(flavor=SCTransformFlavor.ANALYTIC)

    _ = normalize_sctransform(synthetic_adata, config)
    assert np.array_equal(synthetic_adata.X, orig_x)
    assert "pearson_residuals" not in synthetic_adata.layers
