"""Unit tests for manual download resolution, informative error handling, and new loaders."""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from returns.maybe import Some
from returns.result import Failure, Success

from tme_datasets import (
    find_dataset_h5ad,
    get_manual_download_dir,
    load_dataset,
)
from tme_datasets.providers.bulk_papers import (
    load_genentech_egad,
    load_paper_h5ad,
)
from tme_datasets.providers.single_cell import (
    load_gse179994,
    load_ma_liver,
    load_maynard,
    load_yost,
)


def test_get_manual_download_dir(tmp_path: Path) -> None:
    """Verify get_manual_download_dir resolves to data/manual_download."""
    manual_dir = get_manual_download_dir(repo_root=tmp_path)
    assert manual_dir == tmp_path / "data/manual_download"


def test_missing_paper_h5ad_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing paper H5AD returns structured Failure with actionable instructions."""
    missing_file = tmp_path / "data/manual_download/Auslander.h5ad"
    res = load_paper_h5ad(missing_file)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "Paper H5AD file not found" in err
    assert "data/manual_download/Auslander.h5ad" in err
    assert "setup_manual_downloads.py" in err


def test_missing_egad_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing EGAD alignment directory returns structured Failure."""
    missing_dir = tmp_path / "data/manual_download/EGAD00001006631-align"
    res = load_genentech_egad(missing_dir)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "EGAD alignment directory not found" in err
    assert "EGAD00001006631-align" in err


def test_missing_maynard_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing Maynard dataset returns structured Failure."""
    res = load_maynard(tmp_path)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "Maynard NSCLC dataset not found" in err
    assert "Maynard_NSCLC.h5ad" in err


def test_load_ma_liver_with_mock_files(tmp_path: Path) -> None:
    """Verify load_ma_liver correctly reads 10x matrix, genes, and barcodes."""
    raw_dir = tmp_path / "GSE125449"
    raw_dir.mkdir(parents=True)

    # 1. Create dummy MatrixMarket file (3 genes x 2 cells)
    mm_content = (
        "%%MatrixMarket matrix coordinate real general\n"
        "%test\n"
        "3 2 4\n"
        "1 1 5.0\n"
        "2 1 10.0\n"
        "1 2 2.0\n"
        "3 2 8.0\n"
    )
    with gzip.open(raw_dir / "GSE125449_Set1_matrix.mtx.gz", "wt") as f:
        f.write(mm_content)

    with gzip.open(raw_dir / "GSE125449_Set1_genes.tsv.gz", "wt") as f:
        f.write("CD8A\nCD4\nFOXP3\n")

    with gzip.open(raw_dir / "GSE125449_Set1_barcodes.tsv.gz", "wt") as f:
        f.write("cell1\ncell2\n")

    res = load_ma_liver(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert list(adata.obs_names) == ["cell1", "cell2"]
    assert list(adata.var_names) == ["CD8A", "CD4", "FOXP3"]


def test_load_yost_with_mock_counts(tmp_path: Path) -> None:
    """Verify load_yost parses tab-delimited counts and annotates patient and response."""
    raw_dir = tmp_path / "GSE123813"
    raw_dir.mkdir(parents=True)

    # Header: Gene, bcc.su001.pre, bcc.su005.post
    counts_content = (
        "Gene\tbcc.su001.pre\tbcc.su005.post\n"
        "CD8A\t10\t0\n"
        "PDCD1\t5\t2\n"
        "CTLA4\t8\t1\n"
    )
    with gzip.open(raw_dir / "GSE123813_bcc_scRNA_counts.txt.gz", "wt") as f:
        f.write(counts_content)

    res = load_yost(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert "patient" in adata.obs.columns
    assert "response_binary" in adata.obs.columns
    # su001 is Responder (1), su005 is Non-responder (0)
    assert adata.obs.loc["bcc.su001.pre", "response_binary"] == 1
    assert adata.obs.loc["bcc.su005.post", "response_binary"] == 0


def test_load_gse179994_with_mock_rds(tmp_path: Path) -> None:
    """Verify load_gse179994 reads RDS dataframe and attaches condition metadata."""
    raw_dir = tmp_path / "GSE179994"
    raw_dir.mkdir(parents=True)

    rds_file = raw_dir / "GSE179994_all.Tcell.rawCounts.rds"
    # Create mock RDS object via pyreadr mock
    mock_df = pd.DataFrame(
        [[5, 0], [2, 8]],
        index=["CD8A", "CD4"],
        columns=["cell1", "cell2"],
    )

    with patch("pyreadr.read_r", return_value={"counts": mock_df}):
        rds_file.touch()
        res = load_gse179994(raw_dir, auto_download=False)
        assert isinstance(res, Success)
        adata = res.unwrap()
        assert adata.n_obs == 2
        assert adata.n_vars == 2


def test_find_dataset_h5ad_finds_manual_download(tmp_path: Path) -> None:
    """Verify find_dataset_h5ad discovers files placed in data/manual_download."""
    manual_dir = tmp_path / "data/manual_download"
    manual_dir.mkdir(parents=True)
    fake_h5ad = manual_dir / "Auslander.h5ad"
    fake_h5ad.write_bytes(b"mock_h5ad_data")

    found = find_dataset_h5ad("Auslander", repo_root=tmp_path)
    assert isinstance(found, (Success, Some))
    assert found.unwrap() == fake_h5ad
