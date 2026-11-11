"""Out-of-core AnnData loading and backed matrix access."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..paths import find_dataset_h5ad, find_repo_root
from ..registry import resolve_dataset_id
from ..types import StorageBackend


def _resolve_backed_target(
    dataset_or_path: str | Path,
    base_dir: Maybe[Path] | Path,
) -> Result[Path, str]:
    """Resolve either a registered dataset ID or an explicit H5AD/Zarr path."""
    if isinstance(dataset_or_path, Path):
        return (
            Success(dataset_or_path)
            if dataset_or_path.exists()
            else Failure(f"File not found: {dataset_or_path}")
        )

    candidate_path = Path(dataset_or_path)
    if candidate_path.suffix in {".h5ad", ".zarr"} and candidate_path.exists():
        return Success(candidate_path)

    root = (
        base_dir
        if isinstance(base_dir, Path)
        else base_dir.value_or(find_repo_root())
    )
    canonical_id = resolve_dataset_id(dataset_or_path)
    match find_dataset_h5ad(canonical_id, repo_root=root):
        case Some(found_path):
            return Success(found_path)
        case _:
            return Failure(
                f"Could not locate preprocessed H5AD for dataset '{dataset_or_path}' "
                f"(canonical '{canonical_id}') via tme_datasets API."
            )


def load_backed(
    dataset_or_path: str | Path,
    backend: StorageBackend = StorageBackend.BACKED_H5AD,
    mode: str = "r",
    *,
    base_dir: Maybe[Path] | Path = Nothing,
) -> Result[ad.AnnData, str]:
    """Open a registered dataset ID or AnnData path in backed mode without loading expression into RAM."""
    match _resolve_backed_target(dataset_or_path, base_dir):
        case Failure(err):
            return Failure(err)
        case Success(file_path):
            pass

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
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to load dataset in backed mode from {file_path}: {exc}")


def slice_backed_dataset(
    adata: ad.AnnData,
    indices: Sequence[int] | npt.NDArray[np.intp],
) -> Result[ad.AnnData, str]:
    """Safely extract a cell subset from a backed AnnData into memory."""
    try:
        sub = adata[indices].to_memory()
        return Success(sub)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to slice backed dataset: {exc}")
