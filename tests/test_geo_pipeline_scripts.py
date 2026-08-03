"""Unit tests for the GEO pipeline scripts (01_download_geo_data.py, 02_process_raw_to_sparse_h5.py, 03_preprocess_datasets_to_h5.py)."""

from pathlib import Path
import pytest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
from returns.result import Success, Failure

import importlib
geo_download = importlib.import_module("scripts.01_download_geo_data")
raw_to_sparse = importlib.import_module("scripts.02_process_raw_to_sparse_h5")
preprocess_h5 = importlib.import_module("scripts.03_preprocess_datasets_to_h5")

DownloadConfig = geo_download.DownloadConfig
GeoDatasetSpec = geo_download.GeoDatasetSpec
download_single_file = geo_download.download_single_file

RawProcessConfig = raw_to_sparse.RawProcessConfig
convert_df_to_sparse_adata = raw_to_sparse.convert_df_to_sparse_adata
process_tsv_gz_to_sparse_h5 = raw_to_sparse.process_tsv_gz_to_sparse_h5

PreprocessConfig = preprocess_h5.PreprocessConfig
DatasetQCSpec = preprocess_h5.DatasetQCSpec
preprocess_single_dataset = preprocess_h5.preprocess_single_dataset


def test_download_config() -> None:
    """Tests GeoDatasetSpec and DownloadConfig initialization."""
    spec = GeoDatasetSpec(accession="GSE12345", urls={"test": "https://example.com/test.txt.gz"})
    assert spec.accession == "GSE12345"
    assert "test" in spec.urls

    config = DownloadConfig(out_dir=Path("/tmp/download_test"))
    assert len(config.datasets) >= 3


def test_convert_df_to_sparse_adata() -> None:
    """Tests conversion of a pandas DataFrame expression matrix into a sparse CSR AnnData."""
    df_raw = pd.DataFrame(
        np.random.randint(10, 500, size=(100, 20)),
        index=[f"Cell_{i}" for i in range(100)],
        columns=[f"Gene_{j}" for j in range(20)],
    )

    res = convert_df_to_sparse_adata(df_raw, "TestDataset", min_cell_counts=50)
    assert isinstance(res, Success)
    adata = res.unwrap()

    assert isinstance(adata.X, csr_matrix)
    assert adata.n_obs > 0
    assert adata.n_vars == 20
    assert "original.barcode" in adata.obs.columns


def test_process_tsv_gz_to_sparse_h5(tmp_path: Path) -> None:
    """Tests parsing a TSV matrix and saving as a sparse .h5ad file."""
    df_raw = pd.DataFrame(
        np.random.randint(10, 500, size=(50, 15)),
        index=[f"Cell_{i}" for i in range(50)],
        columns=[f"Gene_{j}" for j in range(15)],
    )
    
    tsv_path = tmp_path / "raw_matrix.tsv.gz"
    df_raw.to_csv(tsv_path, sep="\t", compression="gzip")

    out_h5 = tmp_path / "output_sparse.h5ad"
    res = process_tsv_gz_to_sparse_h5(tsv_path, "TestDS", out_h5, min_cell_counts=50)
    assert isinstance(res, Success)
    assert out_h5.exists()

    loaded = ad.read_h5ad(out_h5)
    assert loaded.n_obs == 50
    assert loaded.n_vars == 15


def test_preprocess_single_dataset(tmp_path: Path) -> None:
    """Tests dataset QC filtering, normalization, and preprocessed .h5ad export."""
    X = csr_matrix(np.random.randint(50, 1000, size=(60, 20), dtype=np.int64))
    obs = pd.DataFrame({"original.barcode": [f"Cell_{i}" for i in range(60)]})
    var = pd.DataFrame({"contig": ["1"] * 18 + ["MT"] * 2}, index=[f"ENSG00000{i:06d}" for i in range(20)])
    
    sparse_h5 = tmp_path / "sparse_input.h5ad"
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(sparse_h5)

    out_h5 = tmp_path / "processed_output.h5ad"
    qc_spec = DatasetQCSpec(min_genes=5, max_genes=100, min_counts=50, max_counts=50000, max_mt_content=50.0)

    res = preprocess_single_dataset(sparse_h5, qc_spec, out_h5)
    assert isinstance(res, Success)
    assert out_h5.exists()

def test_check_and_normalize_transform_space() -> None:
    """Tests detection of prior log-transformation and continuous float values."""
    # Log-transformed AnnData
    X_log = csr_matrix(np.random.uniform(0.0, 5.0, size=(20, 10)))
    adata_log = ad.AnnData(X=X_log)
    
    preprocess_h5 = importlib.import_module("scripts.03_preprocess_datasets_to_h5")
    adata_res, status = preprocess_h5.check_and_normalize_transform_space(adata_log)

    assert status["is_log1p"] is True
    assert "log1p_original" in adata_res.layers
    # Linearized values should have larger max than log values
    assert adata_res.X.data.max() > X_log.data.max()
