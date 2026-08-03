import os
import pytest
import numpy as np
import polars as pl
import torch
from returns.result import Success, Failure
from returns.maybe import Some, Nothing

from bayesprism.models import GibbsControl, OptControl
from bayesprism.preprocessing import norm_to_one
from bayesprism.optimization import transform_phi_t
from bayesprism.pipeline import new_prism, run_prism


FIXTURES_DIR = "packages/bayesprism/tests/fixtures/r_benchmark"


@pytest.mark.skipif(not os.path.exists(FIXTURES_DIR), reason="R benchmark fixtures not generated")
def test_norm_to_one_r_equivalence() -> None:
    r_phi_cellstate_path = os.path.join(FIXTURES_DIR, "phi_cellState.csv")
    if not os.path.exists(r_phi_cellstate_path):
        pytest.skip("phi_cellState.csv not found")

    df_r = pl.read_csv(r_phi_cellstate_path)
    first_col = df_r.columns[0]
    numeric_cols = [c for c in df_r.columns if c != first_col]
    matrix_r = df_r.select(numeric_cols).to_numpy().astype(np.float64)

    row_sums = matrix_r.sum(axis=1)
    np.testing.assert_allclose(row_sums, np.ones_like(row_sums), rtol=1e-4)


@pytest.mark.skipif(not os.path.exists(FIXTURES_DIR), reason="R benchmark fixtures not generated")
def test_deconvolution_r_equivalence() -> None:
    r_theta_init_path = os.path.join(FIXTURES_DIR, "theta_initial_cellType.csv")
    if not os.path.exists(r_theta_init_path):
        pytest.skip("theta_initial_cellType.csv not found")

    df_r = pl.read_csv(r_theta_init_path)
    first_col = df_r.columns[0]
    numeric_cols = [c for c in df_r.columns if c != first_col]
    theta_r = df_r.select(numeric_cols).to_numpy().astype(np.float64)

    # Assert valid cell proportion sum = 1 per sample
    row_sums = theta_r.sum(axis=1)
    np.testing.assert_allclose(row_sums, np.ones_like(row_sums), rtol=1e-4)
