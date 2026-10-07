"""Out-of-core storage, backed matrix reading, and Zarr conversion."""

from .incremental_writer import H5ADSparseIncrementalWriter
from .reader import load_backed, slice_backed_dataset
from .zarr_converter import convert_to_zarr

__all__ = [
    "load_backed",
    "slice_backed_dataset",
    "convert_to_zarr",
    "H5ADSparseIncrementalWriter",
]
