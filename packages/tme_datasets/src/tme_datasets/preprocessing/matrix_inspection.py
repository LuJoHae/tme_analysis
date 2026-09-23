"""Pure functions for verifying and classifying single-cell expression matrix representations."""

from __future__ import annotations

from enum import Enum
import anndata as ad
import numpy as np
import scipy.sparse as sp
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

from ..logging import get_logger

logger = get_logger("preprocessing.matrix_inspection")


class ExpressionType(str, Enum):
    """Categorical classification of single-cell expression matrix scale and format."""
    RAW_COUNTS = "raw_counts"
    LINEAR_NORMALIZED = "linear_normalized"
    LOG_NORMALIZED = "log_normalized"
    UNKNOWN = "unknown"


class ExpressionInspectionResult(BaseModel):
    """Immutable report characterizing expression values, scale, and integrality."""
    model_config = ConfigDict(frozen=True)

    expression_type: ExpressionType
    is_raw_counts: bool
    is_integer: bool
    min_value: float
    max_value: float
    mean_library_depth: float
    cv_library_depth: float  # coefficient of variation = std / mean
    n_zeros: int
    sparsity: float


def inspect_expression_type(adata: ad.AnnData) -> ExpressionInspectionResult:
    """Inspects an AnnData matrix to determine if it contains raw counts or normalized expression.

    Tests:
    1. Integrality: Are all non-zero entries exact non-negative integers?
    2. Library size dispersion: In normalized data (CPM / 10k), sum per cell is invariant (CV < 1e-3).
    3. Value range: Max value <= 30 and non-integer indicates log-scale expression.
    """
    X = adata.X
    n_obs, n_vars = adata.n_obs, adata.n_vars
    total_elements = n_obs * n_vars

    if sp.issparse(X):
        data = X.data
        n_zeros = int(total_elements - len(data))
    else:
        arr = np.asarray(X)
        data = arr[arr != 0]
        n_zeros = int(np.sum(arr == 0))

    if len(data) == 0:
        return ExpressionInspectionResult(
            expression_type=ExpressionType.UNKNOWN,
            is_raw_counts=False,
            is_integer=True,
            min_value=0.0,
            max_value=0.0,
            mean_library_depth=0.0,
            cv_library_depth=0.0,
            n_zeros=n_zeros,
            sparsity=1.0,
        )

    min_val = float(np.min(data))
    max_val = float(np.max(data))
    sparsity = float(n_zeros / total_elements) if total_elements > 0 else 0.0

    # Test 1: Integrality (allow tiny floating point rounding tolerances)
    is_non_negative = min_val >= -1e-6
    is_integer = is_non_negative and bool(np.allclose(data, np.round(data), atol=1e-4))

    # Test 2: Library size dispersion across cells
    if sp.issparse(X):
        cell_sums = np.asarray(X.sum(axis=1)).flatten()
    else:
        cell_sums = np.asarray(np.sum(X, axis=1)).flatten()

    mean_depth = float(np.mean(cell_sums))
    std_depth = float(np.std(cell_sums))
    cv_depth = float(std_depth / mean_depth) if mean_depth > 0 else 0.0

    # Classification logic
    if is_integer and cv_depth > 0.05:
        # Raw integer UMI counts with naturally variable cell depths
        expr_type = ExpressionType.RAW_COUNTS
        is_raw = True
    elif not is_integer and max_val <= 30.0:
        # Non-integer with max <= 30 is characteristic of log1p / log2
        expr_type = ExpressionType.LOG_NORMALIZED
        is_raw = False
    elif not is_integer and (cv_depth < 0.01 or max_val > 30.0):
        # Fractional linear CPM/TPM or library-size scaled values
        expr_type = ExpressionType.LINEAR_NORMALIZED
        is_raw = False
    else:
        expr_type = ExpressionType.RAW_COUNTS if is_integer else ExpressionType.UNKNOWN
        is_raw = is_integer

    logger.debug(
        "Expression inspection: %s (is_raw=%s, is_integer=%s, min=%.2f, max=%.2f, cv_depth=%.4f)",
        expr_type.value,
        is_raw,
        is_integer,
        min_val,
        max_val,
        cv_depth,
    )

    return ExpressionInspectionResult(
        expression_type=expr_type,
        is_raw_counts=is_raw,
        is_integer=is_integer,
        min_value=min_val,
        max_value=max_val,
        mean_library_depth=mean_depth,
        cv_library_depth=cv_depth,
        n_zeros=n_zeros,
        sparsity=sparsity,
    )


def tag_expression_metadata(adata: ad.AnnData) -> ad.AnnData:
    """Inspects adata and records expression metadata in adata.uns and adata.layers purely."""
    inspection = inspect_expression_type(adata)
    new_adata = adata.copy()

    new_adata.uns["is_raw_counts"] = inspection.is_raw_counts
    new_adata.uns["expression_type"] = inspection.expression_type.value

    if inspection.is_raw_counts:
        if "counts" not in new_adata.layers:
            new_adata.layers["counts"] = new_adata.X.copy()
    elif inspection.expression_type == ExpressionType.LOG_NORMALIZED:
        if "log1p" not in new_adata.layers:
            new_adata.layers["log1p"] = new_adata.X.copy()
        # Create linear layer expm1(X) if missing for linear deconvolution compatibility
        if "linear" not in new_adata.layers:
            if sp.issparse(new_adata.X):
                linear_data = np.expm1(new_adata.X.data)
                linear_mat = sp.csr_matrix(
                    (linear_data, new_adata.X.indices, new_adata.X.indptr),
                    shape=new_adata.X.shape,
                )
            else:
                linear_mat = np.expm1(np.asarray(new_adata.X, dtype=np.float32))
            new_adata.layers["linear"] = linear_mat
    elif inspection.expression_type == ExpressionType.LINEAR_NORMALIZED:
        if "cpm" not in new_adata.layers:
            new_adata.layers["cpm"] = new_adata.X.copy()

    return new_adata
