"""Unit tests for scripts/unify_step07_gse120575.py."""

from pathlib import Path
import pytest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
from returns.result import Success, Failure
from returns.maybe import Some, Nothing

from scripts.unify_step07_gse120575 import (
    UnifyConfig,
    IntegrationMetrics,
    align_and_harmonize_metadata,
    combine_datasets,
    evaluate_integration_quality,
    save_unified_adata,
    subsample_dataset,
)


@pytest.fixture
def mock_step07_adata() -> ad.AnnData:
    """Fixture creating a synthetic Step07 AnnData dataset."""
    X = csr_matrix(np.random.randint(0, 50, size=(100, 20), dtype=np.int64))
    obs = pd.DataFrame({
        "dataset": ["Pelka2021"] * 50 + ["Azizi2018"] * 50,
        "cell_type": ["tumor"] * 40 + ["tcell"] * 60,
        "patient": [f"P{i%5}" for i in range(100)],
        "original.barcode": [f"BC7_{i}" for i in range(100)],
    })
    var = pd.DataFrame(index=[f"ENSG00000{i:06d}" for i in range(20)])
    var["contig"] = "1"
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def mock_gse120575_adata() -> ad.AnnData:
    """Fixture creating a synthetic GSE120575 AnnData dataset."""
    X = csr_matrix(np.random.rand(80, 20).astype(np.float32) * 100.0)
    obs = pd.DataFrame({
        "dataset": ["GSE120575_SadeFeldman"] * 80,
        "cell_type": ["tumor"] * 30 + ["tcell"] * 50,
        "Patient_ID": [f"GSE_P{i%4}" for i in range(80)],
        "original.barcode": [f"BC_GSE_{i}" for i in range(80)],
        "sequencing_tech": ["Smart-seq2"] * 80,
    })
    var = pd.DataFrame(index=[f"ENSG00000{i:06d}" for i in range(20)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_align_and_harmonize_metadata(
    mock_step07_adata: ad.AnnData,
    mock_gse120575_adata: ad.AnnData,
) -> None:
    """Tests gene feature alignment and metadata column standardization."""
    result = align_and_harmonize_metadata(mock_step07_adata, mock_gse120575_adata)
    assert isinstance(result, Success)
    a7, ag = result.unwrap()

    assert a7.n_vars == 20
    assert ag.n_vars == 20
    assert "sequencing_tech" in a7.obs.columns
    assert a7.obs["sequencing_tech"].iloc[0] == "10x_UMI"
    assert ag.obs["sequencing_tech"].iloc[0] == "Smart-seq2"


def test_subsample_dataset(mock_step07_adata: ad.AnnData) -> None:
    """Tests dataset downsampling by cell type."""
    subsampled = subsample_dataset(mock_step07_adata, n_per_type=10, seed=42)
    assert subsampled.n_obs == 20  # 10 for tumor + 10 for tcell


def test_combine_datasets(
    mock_step07_adata: ad.AnnData,
    mock_gse120575_adata: ad.AnnData,
) -> None:
    """Tests dataset concatenation with subsampling configuration."""
    aligned_res = align_and_harmonize_metadata(mock_step07_adata, mock_gse120575_adata)
    a7, ag = aligned_res.unwrap()

    config = UnifyConfig(
        lair_dir=Path("/tmp/lair"),
        gse120575_dir=Path("/tmp/gse"),
        output_path=Path("/tmp/out.h5ad"),
        subsample_n_per_type=Some(15),
    )

    combined_res = combine_datasets(a7, ag, config)
    assert isinstance(combined_res, Success)
    combined = combined_res.unwrap()

    assert combined.n_vars == 20
    assert combined.n_obs == 60  # 30 from a7 (15 per type) + 30 from ag (15 per type)


def test_save_unified_adata(
    mock_step07_adata: ad.AnnData,
    tmp_path: Path,
) -> None:
    """Tests saving unified AnnData object to disk in .h5ad format."""
    out_file = tmp_path / "unified.h5ad"
    res = save_unified_adata(mock_step07_adata, out_file)
    assert isinstance(res, Success)
    assert out_file.exists()

    reloaded = ad.read_h5ad(out_file)
    assert reloaded.n_obs == mock_step07_adata.n_obs
    assert reloaded.n_vars == mock_step07_adata.n_vars


def test_evaluate_integration_quality(
    mock_step07_adata: ad.AnnData,
    mock_gse120575_adata: ad.AnnData,
) -> None:
    """Tests integration quality metric evaluation (gene count, correlation, silhouette scores)."""
    aligned_res = align_and_harmonize_metadata(mock_step07_adata, mock_gse120575_adata)
    a7, ag = aligned_res.unwrap()

    config = UnifyConfig(
        lair_dir=Path("/tmp/lair"),
        gse120575_dir=Path("/tmp/gse"),
        output_path=Path("/tmp/out.h5ad"),
    )
    combined = combine_datasets(a7, ag, config).unwrap()

    eval_res = evaluate_integration_quality(combined)
    assert isinstance(eval_res, Success)
    metrics: IntegrationMetrics = eval_res.unwrap()

    assert metrics.num_step07_cells == 100
    assert metrics.num_gse120575_cells == 80
    assert metrics.num_shared_genes == 20
    assert metrics.gene_overlap_ratio == 1.0
    assert isinstance(metrics.marker_correlation_mean, float)
