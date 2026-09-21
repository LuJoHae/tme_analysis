"""Unit tests for subsampling, supersampling, and Negative Binomial perturbations."""

import anndata as ad
import numpy as np
from returns.maybe import Some
from returns.result import Success
from tme_datasets.models import NegativeBinomialConfig, SubsampleSpec
from tme_datasets.transforms.knockout import in_silico_knockout
from tme_datasets.transforms.perturbations import randomize_negative_binomial, simulate_dropout
from tme_datasets.transforms.sampling import subsample_cells, supersample_cells


def test_subsample_fraction(mock_single_cell_adata: ad.AnnData) -> None:
    spec = SubsampleSpec(n_or_fraction=0.5, seed=Some(42))
    res = subsample_cells(mock_single_cell_adata, spec)
    assert isinstance(res, Success)
    sub = res.unwrap()
    assert sub.n_obs == 50


def test_subsample_stratified(mock_single_cell_adata: ad.AnnData) -> None:
    spec = SubsampleSpec(n_or_fraction=40, stratify_by=Some("cell_type"), balanced=True, seed=Some(42))
    res = subsample_cells(mock_single_cell_adata, spec)
    assert isinstance(res, Success)
    sub = res.unwrap()
    assert sub.n_obs == 40
    # Balanced should have equal cells across the 4 types (10 each)
    counts = sub.obs["cell_type"].value_counts()
    for count in counts:
        assert count == 10


def test_supersample_cells(mock_single_cell_adata: ad.AnnData) -> None:
    res = supersample_cells(mock_single_cell_adata, n_target=150, seed=Some(42))
    assert isinstance(res, Success)
    super_adata = res.unwrap()
    assert super_adata.n_obs == 150


def test_randomize_negative_binomial(mock_single_cell_adata: ad.AnnData) -> None:
    config = NegativeBinomialConfig(dispersion=0.2, seed=Some(42))
    res = randomize_negative_binomial(mock_single_cell_adata, config)
    assert isinstance(res, Success)
    perturbed = res.unwrap()
    assert "raw_counts" in perturbed.layers
    assert "randomized_nb" in perturbed.layers
    assert perturbed.X.shape == mock_single_cell_adata.X.shape
    # Non-negative integers
    assert np.all(perturbed.X >= 0)


def test_simulate_dropout(mock_single_cell_adata: ad.AnnData) -> None:
    res = simulate_dropout(mock_single_cell_adata, rate=0.20, seed=42)
    assert isinstance(res, Success)
    dropped = res.unwrap()
    # Number of zeros in dropped should be greater than original
    orig_zeros = np.sum(mock_single_cell_adata.X == 0)
    new_zeros = np.sum(dropped.X == 0)
    assert new_zeros >= orig_zeros


def test_in_silico_knockout(mock_single_cell_adata: ad.AnnData) -> None:
    # Knockout CD8A with 100% efficiency
    res = in_silico_knockout(mock_single_cell_adata, genes=("CD8A",), efficiency=1.0)
    assert isinstance(res, Success)
    ko_adata = res.unwrap()
    cd8a_idx = list(ko_adata.var_names).index("CD8A")
    assert np.all(ko_adata.X[:, cd8a_idx] == 0.0)
