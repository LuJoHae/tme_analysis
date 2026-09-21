"""Tests for Negative Binomial parameter estimation (MoM, MLE, Empirical Bayes) and count resampling."""

import anndata as ad
import numpy as np
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success

from tme_datasets.models import NegativeBinomialConfig
from tme_datasets.transforms.nb_inference import (
    compute_size_factors,
    fit_nb_empirical_bayes,
    fit_nb_mle,
    fit_nb_moments,
    infer_dataset_nb_parameters,
)
from tme_datasets.transforms.perturbations import randomize_negative_binomial
from tme_datasets.types import NBEstimationMethod


@pytest.fixture
def synthetic_count_matrix() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate synthetic single-cell counts with known gene means and dispersions."""
    rng = np.random.default_rng(42)
    N_cells, G_genes = 200, 30

    true_means = np.linspace(1.0, 50.0, G_genes)
    true_alpha = 0.2  # true dispersion

    # Cell library size factors
    cell_depths = rng.uniform(0.5, 1.5, size=N_cells)
    size_factors = cell_depths / np.median(cell_depths)

    # Gamma-Poisson sampling:
    # shape = 1 / alpha, scale = alpha * (s_i * mu_g)
    shape = 1.0 / true_alpha
    X = np.zeros((N_cells, G_genes), dtype=np.float64)

    for i in range(N_cells):
        scales = true_alpha * size_factors[i] * true_means
        lam = rng.gamma(shape=shape, scale=scales)
        X[i] = rng.poisson(lam)

    return X, size_factors, true_means


def test_fit_nb_moments(synthetic_count_matrix):
    """Verify that Method of Moments infers means and dispersions close to ground truth."""
    X, size_factors, true_means = synthetic_count_matrix
    inferred_means, inferred_dispersions = fit_nb_moments(X, size_factors)

    # Inferred means should closely match true means
    rel_mean_error = np.abs(inferred_means - true_means) / true_means
    assert np.mean(rel_mean_error) < 0.15

    # Dispersions should be positive and in reasonable range (~0.2)
    assert np.all(inferred_dispersions > 0)
    assert np.median(inferred_dispersions) == pytest.approx(0.2, abs=0.1)


def test_fit_nb_mle(synthetic_count_matrix):
    """Verify that Newton-Raphson MLE converges and estimates dispersion accurately."""
    X, size_factors, true_means = synthetic_count_matrix
    inferred_means, inferred_dispersions = fit_nb_mle(X, size_factors, max_iter=20)

    # Mean matches ground truth
    rel_mean_error = np.abs(inferred_means - true_means) / true_means
    assert np.mean(rel_mean_error) < 0.15

    # Inferred MLE dispersion should be around true alpha = 0.2
    assert np.all(inferred_dispersions > 0)
    assert np.median(inferred_dispersions) == pytest.approx(0.2, abs=0.1)


def test_fit_nb_empirical_bayes(synthetic_count_matrix):
    """Verify Empirical Bayes fits trend and stabilizes noisy gene estimates."""
    X, size_factors, true_means = synthetic_count_matrix
    inferred_means, inferred_dispersions = fit_nb_empirical_bayes(X, size_factors)

    assert len(inferred_means) == len(true_means)
    assert len(inferred_dispersions) == len(true_means)
    assert np.all(inferred_dispersions > 0)
    assert np.all(inferred_dispersions <= 10.0)


def test_simple_mode_backward_compatibility():
    """Verify that default simple mode without estimation_method preserves raw count randomization."""
    X = np.array([[10.0, 0.0], [0.0, 20.0]], dtype=np.float32)
    adata = ad.AnnData(X=X)
    config = NegativeBinomialConfig(dispersion=0.1, seed=Some(42))

    res = randomize_negative_binomial(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # In simple mode with baseline_rate=0.0, zeros stay 0
    assert out_adata.X[0, 1] == 0.0
    assert out_adata.X[1, 0] == 0.0
    # Non-zeros are randomized
    assert out_adata.X[0, 0] > 0
    assert out_adata.X[1, 1] > 0


def test_resampling_dropout_recovery_with_inferred_parameters():
    """Verify that inferred parameters allow technical dropout zeros to recover counts > 0."""
    # Gene 0 is expressed with mean ~10 across cells, but cell 0 happens to have count 0 (dropout)
    rng = np.random.default_rng(123)
    X = rng.poisson(lam=10.0, size=(100, 5)).astype(np.float32)
    X[0, 0] = 0.0  # Force technical dropout zero in cell 0

    adata = ad.AnnData(X=X)
    config = NegativeBinomialConfig(
        estimation_method=Some(NBEstimationMethod.MOMENTS),
        seed=Some(42),
    )

    res = randomize_negative_binomial(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # Inferred parameters stored in varm
    assert "nb_means" in out_adata.varm
    assert "nb_dispersions" in out_adata.varm

    # After resampling from the inferred gene distribution, cell 0 for gene 0 has non-zero probability
    # Across 100 cells, cell 0 should now sample a plausible expression count > 0
    assert out_adata.X[0, 0] > 0


def test_resampling_biological_zero_invariant():
    """Verify that a gene with mean 0 (biological zero across all cells) remains strictly 0."""
    X = np.zeros((50, 5), dtype=np.float32)
    X[:, 1:] = 15.0  # Genes 1-4 expressed, Gene 0 completely unexpressed (all 0)

    adata = ad.AnnData(X=X)
    config = NegativeBinomialConfig(
        estimation_method=Some(NBEstimationMethod.MOMENTS),
        seed=Some(42),
    )

    res = randomize_negative_binomial(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # Gene 0 was 0 across all cells, mean = 0 => resampled counts MUST be all 0
    assert np.all(out_adata.X[:, 0] == 0.0)
    # Gene 1 was expressed => resampled counts should be > 0
    assert np.mean(out_adata.X[:, 1]) > 5.0


def test_cluster_stratified_inference():
    """Verify cluster-stratified parameter inference using infer_dataset_nb_parameters."""
    rng = np.random.default_rng(99)
    X = np.vstack([
        rng.poisson(lam=2.0, size=(50, 4)),   # Cluster A: low expression
        rng.poisson(lam=50.0, size=(50, 4)),  # Cluster B: high expression
    ]).astype(np.float32)

    obs = {"cell_type": ["CellTypeA"] * 50 + ["CellTypeB"] * 50}
    adata = ad.AnnData(X=X, obs=obs)

    config = NegativeBinomialConfig(
        estimation_method=Some(NBEstimationMethod.MOMENTS),
        cluster_key=Some("cell_type"),
        seed=Some(42),
    )

    res = randomize_negative_binomial(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # Cluster A counts should remain modest, Cluster B counts should remain high
    cluster_a_mean = np.mean(out_adata.X[:50])
    cluster_b_mean = np.mean(out_adata.X[50:])
    assert cluster_a_mean < 10.0
    assert cluster_b_mean > 30.0
