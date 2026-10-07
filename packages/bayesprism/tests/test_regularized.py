"""
Unit tests for Collinearity-Aware Regularized Deconvolution Engine.
"""

import numpy as np
import polars as pl
from returns.result import Success
from returns.maybe import Some, Nothing

from bayesprism.regularized import (
    RegularizedDeconvConfig,
    RegularizedDeconvResult,
    project_onto_simplex,
    build_transcriptomic_graph,
    calibrate_spectral_lambda,
    deconvolve_collinearity_regularized,
)


def test_project_onto_simplex() -> None:
    """Test exact Euclidean projection onto probability simplex."""
    # 1. Random vector with negative entries
    v = np.array([-0.5, 1.2, 0.3, -1.0])
    w = project_onto_simplex(v)
    assert np.all(w >= 0.0)
    assert np.isclose(np.sum(w), 1.0, atol=1e-7)

    # 2. Vector already on simplex
    v_on = np.array([0.2, 0.5, 0.3])
    w_on = project_onto_simplex(v_on)
    assert np.allclose(w_on, v_on, atol=1e-7)

    # 3. Vector of equal positive entries
    v_eq = np.array([5.0, 5.0, 5.0])
    w_eq = project_onto_simplex(v_eq)
    assert np.allclose(w_eq, np.array([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]), atol=1e-7)


def test_build_transcriptomic_graph() -> None:
    """Test Graph Laplacian properties: symmetry, positive semi-definiteness, zero row sums."""
    rng = np.random.default_rng(42)
    G = 100
    S = 4
    phi = rng.gamma(2.0, 1.0, size=(S, G))
    phi /= np.sum(phi, axis=1, keepdims=True)

    weights, laplacian = build_transcriptomic_graph(phi, power=2.0)

    # Symmetry
    assert np.allclose(weights, weights.T, atol=1e-9)
    assert np.allclose(laplacian, laplacian.T, atol=1e-9)

    # Row sums are zero
    assert np.allclose(np.sum(laplacian, axis=1), 0.0, atol=1e-9)

    # Positive semi-definite
    evals = np.linalg.eigvalsh(laplacian)
    assert np.all(evals >= -1e-9)
    assert np.isclose(evals[0], 0.0, atol=1e-7)  # Constant vector in nullspace


def test_condition_number_reduction() -> None:
    """Test that Graph Laplacian regularizer strictly bounds the condition number under collinearity."""
    rng = np.random.default_rng(123)
    G = 300
    base_t1 = rng.gamma(2.0, 1.0, G)
    base_t2 = rng.gamma(2.0, 1.0, G)

    # 2 nearly identical sibling states (correlation > 0.99)
    s1 = base_t1 * (1.0 + rng.normal(0, 0.005, G))
    s2 = base_t1 * (1.0 + rng.normal(0, 0.005, G))
    s3 = base_t2

    phi = np.vstack([s1 / np.sum(s1), s2 / np.sum(s2), s3 / np.sum(s3)])

    # Nominal Poisson Hessian
    mu0 = np.mean(phi, axis=0)
    mu0_safe = np.maximum(mu0, 1e-8)
    h_nom_unreg = (phi / mu0_safe) @ phi.T
    kappa_unreg = np.linalg.cond(h_nom_unreg)
    assert kappa_unreg > 200.0, f"Expected high collinearity, got kappa={kappa_unreg}"

    weights, laplacian = build_transcriptomic_graph(phi)
    l_lap, _ = calibrate_spectral_lambda(phi, laplacian, kappa_target=30.0)

    h_nom_reg = h_nom_unreg + l_lap * laplacian
    kappa_reg = np.linalg.cond(h_nom_reg)

    assert kappa_reg < 35.0, f"Expected bounded kappa <= 35.0, got kappa={kappa_reg}"
    assert kappa_reg < kappa_unreg * 0.2, "Condition number should drop by >80%"


def test_deconvolve_collinearity_regularized_end_to_end() -> None:
    """Test end-to-end deconvolution on synthetic mixture."""
    rng = np.random.default_rng(999)
    G = 200
    S = 4
    N = 5

    base1 = rng.gamma(2.0, 1.0, G)
    base2 = rng.gamma(2.0, 1.0, G)

    s1 = base1 * (1.0 + rng.normal(0, 0.01, G))
    s2 = base1 * (1.0 + rng.normal(0, 0.01, G))
    s3 = base2 * (1.0 + rng.normal(0, 0.01, G))
    s4 = base2 * (1.0 + rng.normal(0, 0.01, G))

    phi = np.vstack([s1 / np.sum(s1), s2 / np.sum(s2), s3 / np.sum(s3), s4 / np.sum(s4)])

    # Generate synthetic mixtures
    true_theta = np.array([
        [0.4, 0.1, 0.3, 0.2],
        [0.1, 0.4, 0.2, 0.3],
        [0.25, 0.25, 0.25, 0.25],
        [0.05, 0.05, 0.45, 0.45],
        [0.5, 0.0, 0.1, 0.4],
    ])

    mixture = np.zeros((N, G))
    for n in range(N):
        p_n = true_theta[n] @ phi
        mixture[n] = rng.multinomial(50_000, p_n)

    state_names = ("State_1A", "State_1B", "State_2A", "State_2B")
    sample_names = (f"Sample_{n}" for n in range(N))

    cfg = RegularizedDeconvConfig(kappa_target=30.0, max_iter=300)
    res = deconvolve_collinearity_regularized(
        mixture=mixture,
        reference=phi,
        state_names=state_names,
        sample_names=tuple(sample_names),
        config=Some(cfg),
    )

    assert isinstance(res, Success)
    result = res.unwrap()

    # Verify simplex properties on output
    assert np.all(result.theta >= -1e-6)
    row_sums = np.sum(result.theta, axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-5)

    # Verify Polars conversion
    df = result.to_polars()
    assert isinstance(df, pl.DataFrame)
    assert df.height == N
    assert set(df.columns) == {"sample_id", "State_1A", "State_1B", "State_2A", "State_2B"}
