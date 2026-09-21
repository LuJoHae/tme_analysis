"""Unit tests for single_cell_immuno_datasets package."""

from pathlib import Path
import polars as pl
from returns.result import Success
from single_cell_immuno_datasets.config import TIER_1_DATASETS, DataDirectories
from single_cell_immuno_datasets.report import generate_summary_report, ALL_EVALUATED_DATASETS


def test_tier1_datasets_config() -> None:
    assert len(TIER_1_DATASETS) == 13
    accessions = [d.accession for d in TIER_1_DATASETS]
    assert "GSE120575" in accessions
    assert "Gondal2025" in accessions


def test_data_directories_default() -> None:
    dirs = DataDirectories()
    assert str(dirs.base_dir) == "/storage/halu/data"
    assert str(dirs.raw_dir) == "/storage/halu/data/raw"
    assert str(dirs.preprocessed_dir) == "/storage/halu/data/preprocessed"


def test_generate_summary_report(tmp_path: Path) -> None:
    custom_dirs = DataDirectories.with_base(tmp_path)
    res = generate_summary_report(custom_dirs)
    assert isinstance(res, Success)
    df, md_path, svg_paths = res.unwrap()
    
    assert isinstance(df, pl.DataFrame)
    assert df.height == 27
    assert md_path.exists()
    assert svg_paths[0].exists()
    assert svg_paths[1].exists()
