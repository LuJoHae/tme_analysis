from typing import Union
import numpy as np
import scipy.sparse as sp
import torch
from returns.result import Result, Success, Failure
from bayesprism.models import MatrixType


def check_matrix_finite_and_positive(mat: MatrixType) -> Result[None, str]:
    """Check that matrix contains no NaN/Inf and no negative numbers."""
    if isinstance(mat, sp.spmatrix):
        data = mat.data
    elif isinstance(mat, torch.Tensor):
        if mat.is_sparse:
            data = mat.values().detach().cpu().numpy()
        else:
            data = mat.detach().cpu().numpy()
    else:
        data = np.asarray(mat)

    if not np.all(np.isfinite(data)):
        return Failure("Error: input contains NaN or NA values.")
    if np.any(data < 0):
        return Failure("Error: input contains negative values. Please make sure your input is untransformed raw count.")
    return Success(None)


def check_linear_count_space(mat: MatrixType) -> Result[None, str]:
    """Check if matrix is in linear count space rather than log-transformed or normalized."""
    if isinstance(mat, sp.spmatrix):
        max_val = float(mat.data.max()) if mat.nnz > 0 else 0.0
    elif isinstance(mat, torch.Tensor):
        max_val = float(mat.max())
    else:
        max_val = float(np.max(mat))

    if max_val <= 1.0:
        return Failure("Error: input seems to be normalized. BayesPrism requires unnormalized raw counts.")
    if max_val < 20.0:
        return Failure("Error: input max value < 20. Data appears log-transformed. Log transformation should be avoided.")

    return Success(None)


def validate_input(
    mat: MatrixType,
    gene_names: tuple[str, ...],
) -> Result[None, str]:
    """Validate gene expression matrix and metadata."""
    if len(gene_names) == 0:
        return Failure("Error: please specify the gene names of mixture / reference!")

    match check_matrix_finite_and_positive(mat):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    match check_linear_count_space(mat):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    return Success(None)
