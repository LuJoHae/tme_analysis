import pytest
import numpy as np
import scipy.sparse as sp
import torch
from returns.result import Success, Failure
from bayesprism.validation import (
    check_matrix_finite_and_positive,
    check_linear_count_space,
    validate_input,
)


def test_finite_and_positive_valid() -> None:
    mat = np.array([[10, 20], [30, 40]], dtype=np.float32)
    match check_matrix_finite_and_positive(mat):
        case Success(_):
            assert True
        case Failure(err):
            pytest.fail(f"Expected Success, got Failure: {err}")


def test_finite_and_positive_negative_rejected() -> None:
    mat = np.array([[-1, 20], [30, 40]], dtype=np.float32)
    match check_matrix_finite_and_positive(mat):
        case Failure(err):
            assert "negative values" in err
        case Success(_):
            pytest.fail("Expected Failure for negative matrix")


def test_linear_count_space_log_transformed_rejected() -> None:
    mat = np.array([[1.5, 2.3], [0.5, 3.1]], dtype=np.float32)
    match check_linear_count_space(mat):
        case Failure(err):
            assert "log-transformed" in err or "normalized" in err
        case Success(_):
            pytest.fail("Expected Failure for low max count matrix")


def test_validate_input_valid() -> None:
    mat = np.array([[100, 200], [300, 400]], dtype=np.float32)
    gene_names = ("gene1", "gene2")
    match validate_input(mat, gene_names):
        case Success(_):
            assert True
        case Failure(err):
            pytest.fail(f"Expected Success, got Failure: {err}")
