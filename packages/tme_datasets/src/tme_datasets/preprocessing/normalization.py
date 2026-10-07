"""Pure functional library size normalization (TPM/CPM) and log1p transformation."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success

from ..logging import get_logger

logger = get_logger("preprocessing.normalization")


def normalize_total_counts(
    adata: ad.AnnData,
    target_sum: float = 1e6,
) -> Result[ad.AnnData, str]:
    """Normalize cell/sample library sizes to a fixed target sum (e.g. CPM/TPM = 10^6)."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        is_sparse = sp.issparse(X)

        counts_per_cell = np.asarray(X.sum(axis=1)).flatten()
        counts_per_cell[counts_per_cell == 0] = 1.0  # Avoid zero division
        scale_factors = (target_sum / counts_per_cell).astype(np.float32)

        if is_sparse:
            scaling_matrix = sp.diags(scale_factors, dtype=np.float32)
            new_adata.X = (scaling_matrix @ X).astype(np.float32)
        else:
            new_adata.X = (X * scale_factors[:, np.newaxis]).astype(np.float32)

        return Success(new_adata)
    except Exception as exc:
        msg = f"Failed to normalize library sizes: {exc}"
        logger.error(msg)
        return Failure(msg)


def normalize_to_tpm(
    adata: ad.AnnData,
    target_sum: float = 1e6,
) -> Result[ad.AnnData, str]:
    """Alias for normalize_total_counts to explicitly denote Transcripts Per Million (TPM)."""
    return normalize_total_counts(adata, target_sum=target_sum)


def log1p_transform(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    """Apply natural log(1 + x) transformation to AnnData expression values in X."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        if sp.issparse(X):
            new_adata.X = X.log1p()
        else:
            new_adata.X = np.log1p(X).astype(np.float32)
        return Success(new_adata)
    except Exception as exc:
        msg = f"Failed to apply log1p transform: {exc}"
        logger.error(msg)
        return Failure(msg)


def expm1_transform(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    """Invert log1p transformation returning linear expression values: exp(x) - 1."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        if sp.issparse(X):
            new_adata.X = X.expm1()
        else:
            new_adata.X = np.expm1(X).astype(np.float32)
        return Success(new_adata)
    except Exception as exc:
        msg = f"Failed to apply expm1 transform: {exc}"
        logger.error(msg)
        return Failure(msg)


def compute_tpm_matrix(
    X: sp.spmatrix | np.ndarray,
    target_sum: float = 1e6,
) -> sp.csr_matrix:
    """Compute sparse linear TPM matrix (count / total * target_sum) as float32 CSR.

    Uses O(nnz) vectorized row scaling instead of sparse matrix multiplication
    to eliminate large working buffer allocations.
    """
    if not sp.isspmatrix_csr(X):
        csr_X = sp.csr_matrix(X, dtype=np.float32)
    else:
        csr_X = X.astype(np.float32) if X.dtype != np.float32 else X

    counts_per_cell = np.asarray(csr_X.sum(axis=1)).flatten()
    counts_per_cell[counts_per_cell == 0] = 1.0
    scale_factors = (target_sum / counts_per_cell).astype(np.float32)
    row_multipliers = np.repeat(scale_factors, np.diff(csr_X.indptr))

    return sp.csr_matrix(
        (csr_X.data * row_multipliers, csr_X.indices.copy(), csr_X.indptr.copy()),
        shape=csr_X.shape,
    )


def compute_log1p_norm_matrix(
    X: sp.spmatrix | np.ndarray,
    target_sum: float = 1e4,
) -> sp.csr_matrix:
    """Compute sparse log1p(count / total * target_sum) matrix as float32 CSR."""
    tpm_mat = compute_tpm_matrix(X, target_sum=target_sum)
    return tpm_mat.log1p().tocsr()


def standardize_processed_layers(
    adata: ad.AnnData,
    target_sum: float = 1e6,
) -> Result[ad.AnnData, str]:
    """Establish canonical multi-layer layout for single-cell processing:

    Layers established:
        - adata.layers['counts']: Raw integer count matrix (sparse CSR float32).
        - adata.layers['tpm']: Linear TPM/CPM normalized matrix (target_sum per cell).
        - adata.layers['log1p']: Natural log-transformed normalized matrix log(1 + TPM).
        - adata.X: Natural log-transformed normalized matrix (aligned with Scanpy/Seurat).
        - adata.raw: Frozen AnnData view containing raw counts.

    Args:
        adata: AnnData with raw counts in .X or .layers['counts'].
        target_sum: Target library sum per cell (default 10^6 for TPM/CPM).

    Returns:
        Success(standardized_adata) or Failure(err).
    """
    try:
        new_adata = adata.copy()
        new_adata.var_names_make_unique()

        # Handle datasets originally supplied as TPM (e.g. Smart-seq2)
        if new_adata.uns.get("expression_type") == "tpm":
            if "tpm" in new_adata.layers:
                tpm_matrix = new_adata.layers["tpm"]
            else:
                tpm_matrix = new_adata.X

            if not sp.isspmatrix_csr(tpm_matrix):
                tpm_matrix = sp.csr_matrix(tpm_matrix, dtype=np.float32)
            elif tpm_matrix.dtype != np.float32:
                tpm_matrix = tpm_matrix.astype(np.float32)

            log1p_matrix = tpm_matrix.log1p().tocsr()
            new_adata.layers["tpm"] = tpm_matrix
            new_adata.layers["normalized"] = tpm_matrix
            new_adata.layers["log1p"] = log1p_matrix
            new_adata.layers["log1p_norm"] = log1p_matrix
            new_adata.X = tpm_matrix
            new_adata.raw = new_adata.copy()
            new_adata.uns["normalization"] = {
                "target_sum": target_sum,
                "method": "original_tpm_log1p",
                "layers_created": ["tpm", "normalized", "log1p", "log1p_norm"],
            }
            logger.info("Preserved original linear TPM in .X and .layers['tpm'], generated log1p layers, omitted counts")
            return Success(new_adata)

        # Identify raw count matrix
        if "counts" in new_adata.layers:
            raw_matrix = new_adata.layers["counts"]
        else:
            raw_matrix = new_adata.X

        if not sp.isspmatrix_csr(raw_matrix):
            raw_matrix = sp.csr_matrix(raw_matrix, dtype=np.float32)
        elif raw_matrix.dtype != np.float32:
            raw_matrix = raw_matrix.astype(np.float32)

        # Sanitize nullable string types and enable writing nullable strings
        ad.settings.allow_write_nullable_strings = True
        for col in new_adata.var.columns:
            if hasattr(new_adata.var[col], "dtype") and "string" in str(new_adata.var[col].dtype).lower():
                new_adata.var[col] = new_adata.var[col].astype(object)
        for col in new_adata.obs.columns:
            if hasattr(new_adata.obs[col], "dtype") and "string" in str(new_adata.obs[col].dtype).lower():
                new_adata.obs[col] = new_adata.obs[col].astype(object)

        # Set .layers['counts']
        new_adata.layers["counts"] = raw_matrix

        # Store raw state before normalization
        new_adata.raw = new_adata

        # Compute linear TPM
        tpm_matrix = compute_tpm_matrix(raw_matrix, target_sum=target_sum)
        new_adata.layers["tpm"] = tpm_matrix
        new_adata.layers["normalized"] = tpm_matrix  # Semantic alias

        # Compute log1p(TPM)
        log1p_matrix = tpm_matrix.log1p().tocsr()
        new_adata.layers["log1p"] = log1p_matrix
        new_adata.layers["log1p_norm"] = log1p_matrix
        new_adata.X = log1p_matrix

        new_adata.uns["normalization"] = {
            "target_sum": target_sum,
            "method": "TPM_log1p",
            "layers_created": ["counts", "tpm", "normalized", "log1p"],
        }
        logger.info(
            "Successfully standardized layers: 'counts', 'tpm' (scale=%.0e), and 'log1p' in .X",
            target_sum,
        )
        return Success(new_adata)
    except Exception as exc:
        msg = f"Failed to standardize processed layers: {exc}"
        logger.error(msg)
        return Failure(msg)


def standardize_dual_layers(
    adata: ad.AnnData,
    target_sum: float = 1e4,
) -> Result[ad.AnnData, str]:
    """Standardize AnnData to contain raw counts in .X and .layers['counts'], and log1p_norm in .layers['log1p_norm']."""
    try:
        new_adata = adata.copy()
        new_adata.var_names_make_unique()

        if "counts" in new_adata.layers:
            raw_matrix = new_adata.layers["counts"]
        else:
            raw_matrix = new_adata.X

        if not sp.isspmatrix_csr(raw_matrix):
            raw_matrix = sp.csr_matrix(raw_matrix, dtype=np.float32)
        elif raw_matrix.dtype != np.float32:
            raw_matrix = raw_matrix.astype(np.float32)

        new_adata.X = raw_matrix
        new_adata.layers["counts"] = raw_matrix

        log1p_matrix = compute_log1p_norm_matrix(raw_matrix, target_sum=target_sum)
        new_adata.layers["log1p_norm"] = log1p_matrix
        new_adata.layers["normalized"] = log1p_matrix

        return Success(new_adata)
    except Exception as exc:
        msg = f"Failed to standardize dual layers: {exc}"
        logger.error(msg)
        return Failure(msg)


__all__ = [
    "compute_log1p_norm_matrix",
    "compute_tpm_matrix",
    "expm1_transform",
    "log1p_transform",
    "normalize_to_tpm",
    "normalize_total_counts",
    "standardize_dual_layers",
    "standardize_processed_layers",
]
