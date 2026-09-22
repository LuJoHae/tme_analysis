"""Tests for Jupyter-compatible logging in tme_datasets."""

import io
import logging
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from returns.result import Failure, Success

from tme_datasets import (
    configure_logging,
    get_logger,
    load_dataset,
    normalize_sctransform,
    set_log_level,
    simulate_pseudobulk,
    PseudobulkConfig,
    SCTransformConfig,
    SCTransformFlavor,
)
from tme_datasets.logging import _TmeFlushStreamHandler


@pytest.fixture(autouse=True)
def reset_logging_fixture():
    """Reset logger configuration before each test."""
    yield
    # Restore default INFO level
    configure_logging(level=logging.INFO)


def test_default_logging_initialization():
    """Verify default logger is created and has a flushing stream handler."""
    logger = get_logger("test_module")
    pkg_logger = logging.getLogger("tme_datasets")

    assert pkg_logger.level == logging.INFO
    handlers = [h for h in pkg_logger.handlers if isinstance(h, _TmeFlushStreamHandler)]
    assert len(handlers) == 1
    assert logger.name == "tme_datasets.test_module"


def test_logging_idempotency_for_jupyter():
    """Verify multiple configure_logging calls do not accumulate duplicate handlers."""
    # Simulate a user re-running a notebook cell 5 times
    for _ in range(5):
        configure_logging(level=logging.INFO)

    pkg_logger = logging.getLogger("tme_datasets")
    flush_handlers = [h for h in pkg_logger.handlers if isinstance(h, _TmeFlushStreamHandler)]
    assert len(flush_handlers) == 1


def test_set_log_level_controls_output():
    """Verify set_log_level controls log output visibility."""
    stream = io.StringIO()
    configure_logging(level=logging.INFO, stream=stream)
    logger = get_logger("filter_test")

    # INFO message should appear
    logger.info("Informational test message")
    assert "Informational test message" in stream.getvalue()

    # Clear stream and silence via WARNING
    stream.seek(0)
    stream.truncate(0)
    set_log_level("WARNING")

    logger.info("This should be suppressed")
    assert "This should be suppressed" not in stream.getvalue()

    # WARNING message should appear
    logger.warning("Warning test message")
    assert "Warning test message" in stream.getvalue()


def test_load_dataset_logs_output():
    """Verify load_dataset emits start and completion/failure logs."""
    stream = io.StringIO()
    configure_logging(level=logging.INFO, stream=stream)

    res = load_dataset("NonExistentCohort")
    assert isinstance(res, Failure)

    output = stream.getvalue()
    assert "Loading dataset 'NonExistentCohort'" in output
    assert "is not recognized in the registry" in output


def test_sctransform_logs_progress():
    """Verify normalize_sctransform emits informational logs."""
    stream = io.StringIO()
    configure_logging(level=logging.INFO, stream=stream)

    rng = np.random.default_rng(42)
    X = rng.poisson(lam=5.0, size=(20, 30)).astype(np.float32)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(20)]),
        var=pd.DataFrame(index=[f"gene_{j}" for j in range(30)]),
    )

    cfg = SCTransformConfig(flavor=SCTransformFlavor.ANALYTIC)
    res = normalize_sctransform(adata, cfg)
    assert isinstance(res, Success)

    output = stream.getvalue()
    assert "Starting SCTransform normalization" in output
    assert "Computing Analytic Pearson Residuals" in output
    assert "Successfully completed Analytic Pearson Residuals" in output


def test_pseudobulk_simulation_logs_progress():
    """Verify simulate_pseudobulk emits progress logs."""
    stream = io.StringIO()
    configure_logging(level=logging.INFO, stream=stream)

    rng = np.random.default_rng(42)
    X = rng.poisson(lam=3.0, size=(30, 20)).astype(np.float32)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"cell_type": ["T_cell"] * 15 + ["B_cell"] * 15},
            index=[f"c_{i}" for i in range(30)],
        ),
        var=pd.DataFrame(index=[f"g_{j}" for j in range(20)]),
    )

    cfg = PseudobulkConfig(n_samples=5, cells_per_sample=50)
    res = simulate_pseudobulk(adata, cfg)
    assert isinstance(res, Success)

    output = stream.getvalue()
    assert "Simulating 5 in-silico bulk mixtures" in output
    assert "Successfully simulated 5 pseudobulk mixtures" in output
