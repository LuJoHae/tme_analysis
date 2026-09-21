"""Conversion between AnnData H5AD and chunked Zarr formats."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
from returns.result import Failure, Result, Success


def convert_to_zarr(
    adata: ad.AnnData,
    output_path: Path,
    chunk_size: tuple[int, int] = (1000, 2000),
) -> Result[Path, str]:
    """Export an AnnData object to a high-performance chunked Zarr directory."""
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_zarr(output_path, chunks=chunk_size)
        return Success(output_path)
    except Exception as exc:
        return Failure(f"Failed to export AnnData to Zarr at {output_path}: {exc}")
