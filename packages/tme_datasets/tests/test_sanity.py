"""Tests for Sanity Bayesian Log-Normal Poisson normalization and denoising."""

import anndata as ad
import numpy as np
import pytest
from returns.maybe import Some
from returns.result import Success

from tme_datasets.models import NegativeBinomialConfig, SanityConfig
from tme_datasets.preprocessing.sanity import (
    _solve_qg,
    run_sanity_normalization,
    run_sanity_single_gene,
)
from tme_datasets.transforms.perturbations import randomize_negative_binomial
from tme_datasets.types import NBEstimationMethod


def test_solve_qg_consistency():
    """Verify that q_g solution satisfies sum_c f_{gc} = 1."""
    rng = np.random.default_rng(42)
    C = 100
    N_c = rng.uniform(1000, 5000, size=C)
    n_gc = rng.poisson(lam=5.0, size=C)
    n_g = float(np.sum(n_gc))
    v_g = 0.5
    y_gc = np.log(N_c) + v_g * n_gc

    q_g = _solve_qg(y_gc, v_g, n_g)

    # Check sum of f_gc equals 1
    from scipy.special import lambertw

    scale = v_g * (n_g + 1.0)
    log_z = np.clip(-q_g + y_gc + np.log(scale), -50.0, 50.0)
    w_val = np.real(lambertw(np.exp(log_z)))
    f_gc = w_val / scale
    assert np.sum(f_gc) == pytest.approx(1.0, abs=1e-3)


def test_run_sanity_single_gene_variance_recovery():
    """Verify that Sanity correctly infers high variance for variable genes and low variance for invariant genes."""
    rng = np.random.default_rng(123)
    C = 150
    N_c = np.full(C, 2000.0)
    v_grid = np.geomspace(0.01, 10.0, num=30)

    # Gene 1: High true biological variance across two bimodal subpopulations
    n_gc_variable = np.concatenate([
        rng.poisson(lam=2.0, size=75),
        rng.poisson(lam=40.0, size=75),
    ]).astype(np.float64)

    # Gene 2: Pure Poisson noise (constant baseline across all cells, true biological variance ~ 0)
    n_gc_invariant = rng.poisson(lam=10.0, size=C).astype(np.float64)

    _, _, vg_variable = run_sanity_single_gene(n_gc_variable, N_c, v_grid)
    _, _, vg_invariant = run_sanity_single_gene(n_gc_invariant, N_c, v_grid)

    assert vg_variable > vg_invariant
    # Invariant gene should have small inferred variance
    assert vg_invariant < 0.2


def test_run_sanity_normalization_anndata():
    """Verify Sanity normalization produces correct layers and metadata on AnnData."""
    rng = np.random.default_rng(42)
    N_cells, G_genes = 60, 10

    X = rng.poisson(lam=5.0, size=(N_cells, G_genes)).astype(np.float32)
    # Gene 0 is unexpressed
    X[:, 0] = 0.0

    adata = ad.AnnData(X=X)
    config = SanityConfig(n_bins=20)

    res = run_sanity_normalization(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    # Check output layers
    assert "sanity_ltq" in out_adata.layers
    assert "sanity_error" in out_adata.layers
    assert "sanity_variance" in out_adata.var
    assert "sanity_baseline_alpha" in out_adata.var

    # Dimensions match
    assert out_adata.layers["sanity_ltq"].shape == (N_cells, G_genes)
    assert out_adata.layers["sanity_error"].shape == (N_cells, G_genes)

    # Unexpressed gene has 0 LTQ and 0 error
    assert np.all(out_adata.layers["sanity_ltq"][:, 0] == 0.0)


def test_sanity_resampling_integration():
    """Verify that randomize_negative_binomial with NBEstimationMethod.SANITY executes successfully."""
    rng = np.random.default_rng(99)
    X = rng.poisson(lam=8.0, size=(50, 8)).astype(np.float32)
    adata = ad.AnnData(X=X)

    config = NegativeBinomialConfig(
        estimation_method=Some(NBEstimationMethod.SANITY),
        seed=Some(42),
    )

    res = randomize_negative_binomial(adata, config)
    assert isinstance(res, Success)
    out_adata = res.unwrap()

    assert out_adata.shape == (50, 8)
    assert "randomized_nb" in out_adata.layers
    assert "raw_counts" in out_adata.layers
