"""Out-of-core AnnData loading and backed matrix access."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
from returns.result import Failure, Result, Success

from ..types import StorageBackend


def load_backed(
    file_path: Path,
    backend: StorageBackend = StorageBackend.BACKED_H5AD,
    mode: str = "r",
) -> Result[ad.AnnData, str]:
    """Open an AnnData matrix in backed mode without loading the full expression array into RAM."""
    if not file_path.exists():
        return Failure(f"File not found: {file_path}")

    try:
        match backend:
            case StorageBackend.BACKED_H5AD:
                adata = ad.read_h5ad(file_path, backed=mode)
                return Success(adata)
            case StorageBackend.ZARR:
                adata = ad.read_zarr(file_path)
                return Success(adata)
            case StorageBackend.MEMORY:
                adata = ad.read_h5ad(file_path)
                return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load dataset in backed mode from {file_path}: {exc}")


def slice_backed_dataset(
    adata: ad.AnnData,
    indices: list[int] | np.ndarray,
) -> Result[ad.AnnData, str]:
    """Safely extract a cell subset from a backed AnnData into memory."""
    try:
        sub = adata[indices].to_memory()
        return Success(sub)
    except Exception as exc:
        return Failure(f"Failed to slice backed dataset: {exc}")
