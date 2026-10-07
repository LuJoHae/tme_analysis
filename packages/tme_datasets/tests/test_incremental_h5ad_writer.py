"""Unit tests for H5ADSparseIncrementalWriter."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from tme_datasets.storage.incremental_writer import H5ADSparseIncrementalWriter


def test_incremental_writer_multi_batch(tmp_path: Path) -> None:
    """Verify writing multiple sparse batches results in a valid AnnData H5AD."""
    out_h5ad = tmp_path / "test_incremental.h5ad"
    n_vars = 40
    genes = [f"ENSG_{i:04d}" for i in range(n_vars)]
    var_df = pd.DataFrame({"symbol": [f"GENE_{i}" for i in range(n_vars)]}, index=genes)

    batches = []
    n_batches = 4
    batch_cells = 50
    for b in range(n_batches):
        obs = pd.DataFrame(
            {"donor": f"donor_{b}", "condition": "treatment" if b % 2 == 0 else "control"},
            index=[f"b{b}_cell_{c}" for c in range(batch_cells)],
        )
        # Sparse matrix with random values
        X = sp.random(batch_cells, n_vars, density=0.3, format="csr", dtype=np.float32)
        batches.append((obs, X))

    # Stream batches
    with H5ADSparseIncrementalWriter(out_h5ad, var=var_df) as writer:
        for obs, X in batches:
            writer.append_batch(obs, X)

    assert out_h5ad.exists()
    assert not out_h5ad.with_suffix(".h5ad.tmp").exists()

    # Read back with AnnData
    adata = ad.read_h5ad(out_h5ad)
    assert adata.shape == (n_batches * batch_cells, n_vars)
    assert sp.isspmatrix_csr(adata.X)
    assert list(adata.var_names) == genes
    assert adata.obs.shape[0] == n_batches * batch_cells
    assert list(adata.obs["donor"].unique()) == [f"donor_{b}" for b in range(n_batches)]

    # Verify concatenated matrix values match expected
    expected_X = sp.vstack([b[1] for b in batches], format="csr")
    diff = np.abs(adata.X.toarray() - expected_X.toarray()).max()
    assert diff < 1e-6


def test_incremental_writer_empty_column_sanitization(tmp_path: Path) -> None:
    """Verify unnamed/empty string columns in obs and var are sanitized without h5py errors."""
    out_h5ad = tmp_path / "test_sanitization.h5ad"
    var_df = pd.DataFrame({"": ["unnamed1", "unnamed2"], "symbol": ["G1", "G2"]}, index=["E1", "E2"])
    obs_df = pd.DataFrame({"": ["idx1", "idx2"], "label": ["A", "B"]}, index=["c1", "c2"])
    X = sp.csr_matrix([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32)

    with H5ADSparseIncrementalWriter(out_h5ad, var=var_df) as writer:
        writer.append_batch(obs_df, X)

    adata = ad.read_h5ad(out_h5ad)
    assert "" not in adata.var.columns
    assert "" not in adata.obs.columns
    assert "symbol" in adata.var.columns
    assert "label" in adata.obs.columns


def test_incremental_writer_abort_cleans_tmp(tmp_path: Path) -> None:
    """Verify that an exception during streaming cleans up the staging file without creating target."""
    out_h5ad = tmp_path / "test_abort.h5ad"
    var_df = pd.DataFrame({"symbol": ["G1", "G2"]}, index=["E1", "E2"])
    obs_df = pd.DataFrame({"label": ["A"]}, index=["c1"])
    # Incompatible shape (3 columns instead of 2)
    bad_X = sp.csr_matrix([[1.0, 0.0, 3.0]], dtype=np.float32)

    with pytest.raises(ValueError, match="does not match"):
        with H5ADSparseIncrementalWriter(out_h5ad, var=var_df) as writer:
            writer.append_batch(obs_df, bad_X)

    assert not out_h5ad.exists()
    assert not out_h5ad.with_suffix(".h5ad.tmp").exists()
