"""Incremental H5AD writer for memory-efficient streaming of sparse CSR AnnData objects."""

from __future__ import annotations

import logging
from pathlib import Path
from types import TracebackType
from typing import Any, Mapping, Sequence

import anndata.io
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger("tme_datasets.storage.incremental_writer")


class H5ADSparseIncrementalWriter:
    """Context manager for incrementally writing large sparse CSR AnnData to an H5AD file.

    Features:
    - Zero full-dataset memory footprint: writes batches directly to extendable HDF5 datasets.
    - Strict CSR matrix storage with AnnData 0.1.0/0.2.0 compatibility.
    - Sanitizes observation and feature DataFrames to eliminate empty-string column crashes on Python 3.14 / h5py.
    - Atomic writes: writes to a staging `.tmp` file and replaces the destination only upon success.
    """

    def __init__(
        self,
        output_path: Path,
        var: pd.DataFrame,
        n_vars: int | None = None,
        uns: Mapping[str, Any] | None = None,
        chunk_nnz: int = 65536,
        chunk_indptr: int = 16384,
    ) -> None:
        self.output_path = Path(output_path).resolve()
        self.var = var.copy()
        self.n_vars = n_vars if n_vars is not None else len(self.var)
        self.uns = dict(uns) if uns else {}
        self.chunk_nnz = chunk_nnz
        self.chunk_indptr = chunk_indptr

        self.staging_path = self.output_path.with_suffix(self.output_path.suffix + ".tmp")
        self._file: h5py.File | None = None
        self._grp_x: h5py.Group | None = None
        self._ds_data: h5py.Dataset | None = None
        self._ds_indices: h5py.Dataset | None = None
        self._ds_indptr: h5py.Dataset | None = None

        self._total_obs: int = 0
        self._obs_chunks: list[pd.DataFrame] = []
        self._is_closed: bool = False

    def __enter__(self) -> H5ADSparseIncrementalWriter:
        self.open()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        if exc_type is not None:
            self.abort()
        else:
            self.close()

    def open(self) -> None:
        """Initialize the HDF5 staging file and create the AnnData group structure."""
        self.staging_path.parent.mkdir(parents=True, exist_ok=True)
        if self.staging_path.exists():
            self.staging_path.unlink()

        f = h5py.File(self.staging_path, "w")
        self._file = f

        # Set root AnnData specification attributes
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"

        # Sanitize var DataFrame: drop empty column names and ensure valid index name
        clean_var = self._sanitize_dataframe(self.var, default_index_name="gene_id")
        anndata.io.write_elem(f, "var", clean_var)

        # Create extendable CSR matrix group /X
        grp_x = f.create_group("X")
        grp_x.attrs["encoding-type"] = "csr_matrix"
        grp_x.attrs["encoding-version"] = "0.1.0"
        grp_x.attrs["shape"] = np.array([0, self.n_vars], dtype=np.int64)

        self._ds_data = grp_x.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=np.float32,
            chunks=(self.chunk_nnz,),
        )
        self._ds_indices = grp_x.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(None,),
            dtype=np.int32,
            chunks=(self.chunk_nnz,),
        )
        self._ds_indptr = grp_x.create_dataset(
            "indptr",
            shape=(1,),
            maxshape=(None,),
            dtype=np.int64,
            chunks=(self.chunk_indptr,),
            data=np.array([0], dtype=np.int64),
        )

        self._grp_x = grp_x
        self._total_obs = 0
        self._obs_chunks = []
        self._is_closed = False

    def append_batch(
        self,
        obs_chunk: pd.DataFrame,
        X_chunk: sp.spmatrix | np.ndarray,
    ) -> None:
        """Append a batch of observation rows and sparse matrix data.

        Args:
            obs_chunk: Metadata DataFrame for the cells in this batch.
            X_chunk: Expression matrix for this batch (cells x genes).
        """
        if self._file is None or self._is_closed:
            raise RuntimeError("Cannot append to a closed or uninitialized H5ADSparseIncrementalWriter")

        # Ensure CSR float32 format
        if not sp.isspmatrix_csr(X_chunk):
            X_csr = sp.csr_matrix(X_chunk, dtype=np.float32)
        else:
            X_csr = X_chunk.astype(np.float32) if X_chunk.dtype != np.float32 else X_chunk

        n_rows, n_cols = X_csr.shape
        if n_cols != self.n_vars:
            raise ValueError(
                f"Batch column dimension ({n_cols}) does not match initialized n_vars ({self.n_vars})"
            )
        if n_rows != len(obs_chunk):
            raise ValueError(
                f"Batch matrix rows ({n_rows}) does not match obs rows ({len(obs_chunk)})"
            )

        assert self._ds_data is not None
        assert self._ds_indices is not None
        assert self._ds_indptr is not None

        cur_nnz = self._ds_data.shape[0]
        cur_indptr_len = self._ds_indptr.shape[0]
        batch_nnz = X_csr.nnz

        # 1. Resize and append data and indices
        self._ds_data.resize((cur_nnz + batch_nnz,))
        if batch_nnz > 0:
            self._ds_data[cur_nnz:] = X_csr.data

        self._ds_indices.resize((cur_nnz + batch_nnz,))
        if batch_nnz > 0:
            self._ds_indices[cur_nnz:] = X_csr.indices

        # 2. Resize and append indptr (offsetting row pointers by cur_nnz)
        self._ds_indptr.resize((cur_indptr_len + n_rows,))
        self._ds_indptr[cur_indptr_len:] = X_csr.indptr[1:].astype(np.int64) + cur_nnz

        # 3. Buffer observation metadata
        self._obs_chunks.append(obs_chunk.copy())
        self._total_obs += n_rows

    def close(self) -> None:
        """Finalize observation metadata, root attributes, and atomically commit the file."""
        if self._is_closed or self._file is None:
            return

        f = self._file
        assert self._grp_x is not None

        # 1. Update matrix shape
        self._grp_x.attrs["shape"] = np.array([self._total_obs, self.n_vars], dtype=np.int64)

        # 2. Write combined obs metadata
        if self._obs_chunks:
            combined_obs = pd.concat(self._obs_chunks, axis=0)
        else:
            combined_obs = pd.DataFrame(index=pd.Index([], name="cell_id"))

        clean_obs = self._sanitize_dataframe(combined_obs, default_index_name="cell_id")
        anndata.io.write_elem(f, "obs", clean_obs)

        # 3. Write optional uns metadata
        if self.uns:
            anndata.io.write_elem(f, "uns", self.uns)

        f.close()
        self._file = None
        self._is_closed = True

        # 4. Atomic rename to destination path
        if self.output_path.exists():
            self.output_path.unlink()
        self.staging_path.rename(self.output_path)
        logger.debug(
            "Successfully finalized incremental H5AD: %s (%d obs x %d vars)",
            self.output_path.name,
            self._total_obs,
            self.n_vars,
        )

    def abort(self) -> None:
        """Close HDF5 file and remove incomplete staging file upon failure."""
        if self._file is not None:
            try:
                self._file.close()
            except Exception:
                pass
            self._file = None
        if self.staging_path.exists():
            try:
                self.staging_path.unlink()
            except Exception:
                pass
        self._is_closed = True

    @staticmethod
    def _sanitize_dataframe(df: pd.DataFrame, default_index_name: str) -> pd.DataFrame:
        """Sanitize a DataFrame for error-free HDF5 serialization in Python 3.14.

        - Drops unnamed/empty-string column headers.
        - Replaces missing values in object columns with empty strings.
        - Ensures the index has a non-empty name.
        """
        clean_df = df.copy()

        # Drop columns with empty names
        drop_cols = [c for c in clean_df.columns if not c or str(c).strip() == ""]
        if drop_cols:
            clean_df = clean_df.drop(columns=drop_cols)

        # Ensure index name is set
        if not clean_df.index.name or not isinstance(clean_df.index.name, str):
            clean_df.index.name = default_index_name

        # Sanitize object columns
        for col in clean_df.columns:
            if clean_df[col].dtype == object:
                clean_df[col] = clean_df[col].fillna("").astype(str)

        return clean_df
