"""Unit tests for the native CIBERSORT (linear nu-SVR) adapter in bayesprism."""

import numpy as np
import polars as pl
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success, Failure

from bayesprism.adapters.cibersort import (
    CibersortConfig,
    CibersortDeconvResult,
    deconvolve_cibersort,
)


def test_cibersort_config_frozen() -> None:
    """Verify CibersortConfig is immutable and accepts Maybe parameters."""
    cfg = CibersortConfig(nu_values=(0.25, 0.5), c_param=2.0, n_cpus=Some(2))
    assert cfg.nu_values == (0.25, 0.5)
    assert cfg.c_param == 2.0
    assert cfg.n_cpus == Some(2)

    with pytest.raises(Exception):
        cfg.c_param = 3.0  # type: ignore


def test_cibersort_dimension_mismatch() -> None:
    """Verify error handling on dimension mismatch."""
    genes = tuple(f"Gene_{i}" for i in range(20))
    states = ("StateA", "StateB")
    samples = ("Sample_1",)
    phi = np.ones((3, 20))  # 3 states in matrix, but only 2 state names provided
    mix = np.ones((1, 20))

    match deconvolve_cibersort(mix, phi, states, samples, genes):
        case Failure(err):
            assert "does not match" in err.lower() or "shape" in err.lower()
        case Success(_):
            raise AssertionError("Should have failed on dimension mismatch")


def test_cibersort_synthetic_exact_recovery() -> None:
    """Test end-to-end deconvolution with synthetic mixtures."""
    rng = np.random.default_rng(42)
    G = 150
    genes = tuple(f"Gene_{i:03d}" for i in range(G))
    states = ("Lineage_A", "Lineage_B", "Lineage_C")
    sample_names = ("Sample_01", "Sample_02")

    # Distinct cell-type profiles
    s1 = rng.gamma(5.0, 1.0, G)
    s2 = rng.gamma(2.0, 2.0, G)
    s3 = rng.gamma(1.0, 4.0, G)
    phi = np.vstack([s1, s2, s3])
    phi = phi / np.sum(phi, axis=1, keepdims=True) * 1e5

    # True proportions
    true_theta = np.array([
        [0.5, 0.3, 0.2],
        [0.1, 0.7, 0.2],
    ])

    bulks = true_theta @ phi

    cfg = CibersortConfig(nu_values=(0.25, 0.5, 0.75), n_perm=0)
    match deconvolve_cibersort(bulks, phi, states, sample_names, genes, config=Some(cfg)):
        case Failure(err):
            raise AssertionError(f"deconvolve_cibersort failed: {err}")
        case Success(res):
            assert isinstance(res, CibersortDeconvResult)
            assert res.sample_names == sample_names
            assert res.cell_types == states
            assert "sample_id" in res.proportions.columns

            theta_est = res.proportions.select(list(states)).to_numpy()
            assert theta_est.shape == (2, 3)

            # Check that estimates are close to ground truth (MSE < 0.005, Corr > 0.98)
            for i in range(2):
                mse = float(np.mean((theta_est[i] - true_theta[i]) ** 2))
                corr = float(np.corrcoef(theta_est[i], true_theta[i])[0, 1])
                assert mse < 0.005, f"Sample {i} MSE too high: {mse}"
                assert corr > 0.98, f"Sample {i} Corr too low: {corr}"

            assert len(res.rmse) == 2
            assert len(res.correlation) == 2
            assert len(res.best_nu) == 2


def test_cibersort_multi_threaded() -> None:
    """Test multi-threaded batch execution."""
    rng = np.random.default_rng(123)
    G = 50
    genes = tuple(f"Gene_{i:02d}" for i in range(G))
    states = ("Type_1", "Type_2")
    sample_names = tuple(f"Bulk_{i:02d}" for i in range(6))

    phi = rng.gamma(2.0, 1.0, (2, G))
    true_theta = rng.dirichlet(np.array([1.0, 1.0]), size=6)
    bulks = true_theta @ phi

    cfg = CibersortConfig(n_cpus=Some(2))
    match deconvolve_cibersort(bulks, phi, states, sample_names, genes, config=Some(cfg)):
        case Failure(err):
            raise AssertionError(f"Parallel deconvolve_cibersort failed: {err}")
        case Success(res):
            assert len(res.proportions) == 6
            assert len(res.correlation) == 6


def test_cibersort_permutation_p_value() -> None:
    """Test permutation testing producing empirical p-values."""
    rng = np.random.default_rng(999)
    G = 40
    genes = tuple(f"Gene_{i:02d}" for i in range(G))
    states = ("StateA", "StateB")
    sample_names = ("Sample_P1",)

    phi = rng.gamma(2.0, 1.0, (2, G))
    bulks = np.array([[0.6, 0.4]]) @ phi

    cfg = CibersortConfig(n_perm=10)
    match deconvolve_cibersort(bulks, phi, states, sample_names, genes, config=Some(cfg)):
        case Failure(err):
            raise AssertionError(f"Permutation test failed: {err}")
        case Success(res):
            assert res.p_values != Nothing
            p_val = res.p_values.unwrap()[0]
            assert 0.0 <= p_val <= 1.0
