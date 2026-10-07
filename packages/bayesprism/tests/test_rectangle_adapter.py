"""Unit tests for the Rectangle DWLS-QP adapter in bayesprism."""

import numpy as np
import polars as pl
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success, Failure

from bayesprism.adapters.rectangle import (
    RectangleConfig,
    RectangleDeconvResult,
    create_signature_from_matrix,
    deconvolve_rectangle,
    is_rectangle_available,
)


def test_is_rectangle_available_bool() -> None:
    """Verify that is_rectangle_available returns True natively."""
    assert is_rectangle_available() is True


def test_rectangle_config_frozen() -> None:
    """Verify RectangleConfig is immutable and accepts Maybe parameters."""
    cfg = RectangleConfig(correct_mrna_bias=False, n_cpus=Some(2))
    assert cfg.correct_mrna_bias is False
    assert cfg.n_cpus == Some(2)


def test_create_signature_from_matrix() -> None:
    """Test creating a RectangleSignatureResult from numpy matrix."""
    rng = np.random.default_rng(42)
    genes = tuple(f"Gene_{i}" for i in range(20))
    states = ("StateA", "StateB", "StateC")
    phi = rng.gamma(2.0, 1.0, size=(3, 20))

    match create_signature_from_matrix(phi, states, genes):
        case Success(sig_obj):
            assert sig_obj is not None
            assert hasattr(sig_obj, "signature_genes")
            assert hasattr(sig_obj, "bias_factors")
            assert len(sig_obj.bias_factors) == 3
        case Failure(err):
            raise AssertionError(f"create_signature_from_matrix failed: {err}")


def test_dimension_mismatch_error() -> None:
    """Verify error handling on dimension mismatch."""
    genes = tuple(f"Gene_{i}" for i in range(20))
    states = ("StateA", "StateB")
    phi = np.ones((3, 20))  # 3 rows, but only 2 state names

    match create_signature_from_matrix(phi, states, genes):
        case Failure(err):
            assert "does not match" in err
        case Success(_):
            raise AssertionError("Should have failed on dimension mismatch")


def test_deconvolve_rectangle_synthetic() -> None:
    """Test end-to-end DWLS-QP deconvolution with synthetic mixtures."""
    rng = np.random.default_rng(123)
    G = 30
    genes = tuple(f"Gene_{i:02d}" for i in range(G))
    states = ("Lineage_A", "Lineage_B", "Lineage_C")
    sample_names = ("Sample_01", "Sample_02", "Sample_03")

    # Distinct cell-type signatures
    base_a = rng.gamma(5.0, 1.0, G)
    base_b = rng.gamma(2.0, 2.0, G)
    base_c = rng.gamma(1.0, 4.0, G)
    phi = np.vstack([base_a, base_b, base_c])
    phi = phi / np.sum(phi, axis=1, keepdims=True) * 1e5

    # True proportions
    true_theta = np.array([
        [0.6, 0.3, 0.1],
        [0.1, 0.7, 0.2],
        [0.3, 0.3, 0.4],
    ])

    # Synthetic mixtures
    bulks = true_theta @ phi

    # Build signature object
    sig_res = create_signature_from_matrix(phi, states, genes).unwrap()

    # Run deconvolution
    cfg = RectangleConfig(correct_mrna_bias=False)
    match deconvolve_rectangle(sig_res, bulks, sample_names, genes, config=Some(cfg)):
        case Failure(err):
            raise AssertionError(f"deconvolve_rectangle failed: {err}")
        case Success(res):
            assert isinstance(res, RectangleDeconvResult)
            assert res.sample_names == sample_names
            assert set(res.cell_types) == set(states)
            assert "sample_id" in res.proportions.columns
            assert "Unknown" in res.proportions.columns

            # Check that estimates are non-negative
            for state in states:
                vals = res.proportions[state].to_numpy()
                assert np.all(vals >= -1e-6)

            # Check unknown column is non-negative
            unkn_vals = res.proportions["Unknown"].to_numpy()
            assert np.all(unkn_vals >= -1e-6)

            # Check partition of unity (sum to ~1.0)
            all_cols = list(states) + ["Unknown"]
            total_sum = res.proportions.select(all_cols).to_numpy().sum(axis=1)
            assert np.allclose(total_sum, 1.0, atol=1e-3)

            # Verify proportions are NOT constant uniform fallback
            est_mat = res.proportions.select(list(states)).to_numpy()
            for n in range(3):
                state_std = float(np.std(est_mat[n]))
                assert state_std > 0.03, f"Estimates are flat for sample {n}: std={state_std}"

            # Verify strong correlation with ground truth (r > 0.95)
            corr = float(np.corrcoef(est_mat.ravel(), true_theta.ravel())[0, 1])
            assert corr > 0.95, f"Correlation with ground truth is too low: {corr:.4f}"
