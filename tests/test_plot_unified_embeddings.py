"""Unit tests for scripts/plot_unified_embeddings.py."""

from pathlib import Path
import pytest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
from returns.result import Success, Failure

from scripts.plot_unified_embeddings import (
    PlotConfig,
    apply_gene_length_scaling,
    compute_embeddings_mode,
    extract_embedding_dataframe,
    export_mode_svgs,
    run_plotting_pipeline,
)


@pytest.fixture
def mock_unified_adata() -> ad.AnnData:
    """Fixture creating a synthetic unified AnnData dataset with obs metadata."""
    X = csr_matrix(np.random.randint(0, 50, size=(120, 30), dtype=np.int64))
    obs = pd.DataFrame({
        "dataset": ["Pelka2021"] * 40 + ["Azizi2018"] * 40 + ["GSE120575_SadeFeldman"] * 40,
        "sequencing_tech": ["10x_UMI"] * 80 + ["Smart-seq2"] * 40,
        "cell_type": ["tumor"] * 60 + ["tcell"] * 60,
        "cancer_code": ["CRC"] * 40 + ["BRCA"] * 40 + ["SKCM"] * 40,
    })
    var = pd.DataFrame({
        "contig_length": [1000 + i * 100 for i in range(30)],
    }, index=[f"ENSG00000{i:06d}" for i in range(30)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_apply_gene_length_scaling(mock_unified_adata: ad.AnnData) -> None:
    """Tests gene length scaling for Smart-seq2 cells."""
    scaled = apply_gene_length_scaling(mock_unified_adata)
    assert scaled.shape == mock_unified_adata.shape


def test_compute_embeddings_modes(mock_unified_adata: ad.AnnData) -> None:
    """Tests Uncorrected, ComBat, and Harmony embedding computations."""
    for mode in ["uncorrected", "combat", "harmony"]:
        res = compute_embeddings_mode(mock_unified_adata, mode=mode, n_hvg=20, n_pcs=10)
        assert isinstance(res, Success), f"Failed for mode {mode}: {res}"
        adata_emb = res.unwrap()

        assert "X_pca" in adata_emb.obsm or "X_pca_harmony" in adata_emb.obsm
        assert "X_umap" in adata_emb.obsm
        assert adata_emb.obsm["X_umap"].shape[0] == 120


def test_extract_embedding_dataframe(mock_unified_adata: ad.AnnData) -> None:
    """Tests extraction of PCA/UMAP coordinates and metadata into DataFrame."""
    adata_emb = compute_embeddings_mode(mock_unified_adata, mode="uncorrected", n_hvg=20, n_pcs=10).unwrap()
    res = extract_embedding_dataframe(adata_emb, "Uncorrected")
    assert isinstance(res, Success)
    df = res.unwrap()

    assert "PC1" in df.columns
    assert "UMAP1" in df.columns
    assert "GSE120575_vs_Rest" in df.columns
    assert "Correction_Mode" in df.columns

    # Verify GSE120575 vs Rest classification
    gse_count = (df["GSE120575_vs_Rest"] == "GSE120575 (Melanoma)").sum()
    assert gse_count == 40


def test_export_mode_svgs(
    mock_unified_adata: ad.AnnData,
    tmp_path: Path,
) -> None:
    """Tests generation and export of SVG plot files for a batch mode."""
    adata_emb = compute_embeddings_mode(mock_unified_adata, mode="uncorrected", n_hvg=20, n_pcs=10).unwrap()
    df = extract_embedding_dataframe(adata_emb, "Uncorrected").unwrap()

    res = export_mode_svgs(df, "uncorrected", tmp_path)
    assert isinstance(res, Success)
    paths = res.unwrap()

    assert len(paths) == 10
    for p in paths:
        assert p.exists()
        assert p.suffix == ".svg"


def test_run_plotting_pipeline(
    mock_unified_adata: ad.AnnData,
    tmp_path: Path,
) -> None:
    """Tests complete plotting pipeline with h5ad file input across all modes."""
    h5ad_file = tmp_path / "test_unified.h5ad"
    out_dir = tmp_path / "plots"
    mock_unified_adata.write_h5ad(h5ad_file)

    config = PlotConfig(
        input_h5ad_path=h5ad_file,
        output_dir=out_dir,
        n_hvg=20,
        n_pcs=10,
    )

    res = run_plotting_pipeline(config)
    assert isinstance(res, Success), f"Pipeline failed: {res}"
    paths = res.unwrap()
    assert len(paths) == 30  # 10 for uncorrected + 10 for combat + 10 for harmony
