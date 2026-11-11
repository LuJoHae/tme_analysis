"""Out-of-core storage, backed matrix reading, RDS streaming, and Zarr conversion."""

from .incremental_writer import (
    H5ADSparseIncrementalWriter,
    append_csr_batch_to_group,
    init_compressed_csr_group,
    write_anndata_dataframe_elem,
)
from .rds_stream import (
    SeuratStagingManifest,
    StagedAssaySpec,
    StagedReductionSpec,
    read_r_dataframe,
    stream_r_object_to_staging_h5,
    stream_seurat_rds_to_staging_h5,
)
from .reader import load_backed, slice_backed_dataset
from .zarr_converter import convert_to_zarr

__all__ = [
    "load_backed",
    "slice_backed_dataset",
    "convert_to_zarr",
    "H5ADSparseIncrementalWriter",
    "append_csr_batch_to_group",
    "init_compressed_csr_group",
    "write_anndata_dataframe_elem",
    "SeuratStagingManifest",
    "StagedAssaySpec",
    "StagedReductionSpec",
    "read_r_dataframe",
    "stream_r_object_to_staging_h5",
    "stream_seurat_rds_to_staging_h5",
]
