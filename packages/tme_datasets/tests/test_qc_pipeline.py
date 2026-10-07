"""Unit tests for the Luecken & Theis (2019) single-cell processing pipeline."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
from returns.maybe import Some
from returns.result import Success
import scipy.sparse as sp

from tme_datasets import (
    QualityControlSpec,
    SingleCellProcessingResult,
    SingleCellProcessingSpec,
    apply_quality_control,
    build_qc_dashboard,
    calculate_adaptive_thresholds,
    compute_qc_covariates,
    export_qc_plots,
    extract_qc_metrics_dataframe,
    normalize_to_tpm,
    process_single_cell_dataset,
    standardize_processed_layers,
)


@pytest.fixture
def sample_raw_adata() -> ad.AnnData:
    """Create a deterministic synthetic raw counts AnnData with mitochondrial genes."""
    rng = np.random.default_rng(42)
    # 50 cells x 20 genes
    counts = rng.poisson(lam=5.0, size=(50, 20)).astype(np.float32)

    # Force some low-quality cells (cell 0: empty/low-count, cell 1: dying high-mito)
    counts[0, :] = 0.0
    counts[0, 1] = 5.0  # total count = 5, genes = 1

    # Designate genes 18 and 19 as mitochondrial
    genes = [f"GENE_{i}" for i in range(18)] + ["MT-CO1", "MT-ND1"]
    counts[1, :] = 1.0
    counts[1, 18] = 50.0  # high mito
    counts[1, 19] = 50.0  # high mito

    # Gene 0 expressed in only 1 cell
    counts[:, 0] = 0.0
    counts[2, 0] = 10.0

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(50)])
    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes
    return ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)


def test_compute_qc_covariates(sample_raw_adata: ad.AnnData) -> None:
    """Verify exact computation of count depth, gene detection, and mitochondrial %."""
    adata = compute_qc_covariates(sample_raw_adata)

    assert "total_counts" in adata.obs.columns
    assert "n_genes_by_counts" in adata.obs.columns
    assert "pct_counts_mt" in adata.obs.columns
    assert "is_mitochondrial" in adata.var.columns

    # Cell 0 has only 5 counts in 1 gene
    assert adata.obs.loc["cell_0", "total_counts"] == 5.0
    assert adata.obs.loc["cell_0", "n_genes_by_counts"] == 1
    assert adata.obs.loc["cell_0", "pct_counts_mt"] == 0.0

    # Cell 1 has ~100 mito counts out of ~118 total counts (>80% mito)
    assert adata.obs.loc["cell_1", "pct_counts_mt"] > 80.0

    # Mitochondrial gene mask
    assert bool(adata.var.loc["MT-CO1", "is_mitochondrial"]) is True
    assert bool(adata.var.loc["GENE_1", "is_mitochondrial"]) is False


def test_calculate_adaptive_thresholds(sample_raw_adata: ad.AnnData) -> None:
    """Verify adaptive MAD calculations return mathematically valid thresholds."""
    adata = compute_qc_covariates(sample_raw_adata)
    spec = calculate_adaptive_thresholds(adata, n_mads=3.0)

    assert isinstance(spec, QualityControlSpec)
    assert spec.min_counts_per_cell >= 200
    assert spec.min_genes_per_cell >= 100
    assert spec.max_genes_per_cell <= 15000
    assert 5.0 <= spec.max_pct_mitochondrial <= 20.0


def test_apply_quality_control_filtering(sample_raw_adata: ad.AnnData) -> None:
    """Verify low-quality cells and unexpressed genes are filtered."""
    spec = QualityControlSpec(
        min_counts_per_cell=20,
        min_genes_per_cell=3,
        max_pct_mitochondrial=25.0,
        min_cells_per_gene=3,
    )
    res = apply_quality_control(sample_raw_adata, qc_spec=spec)
    assert isinstance(res, Success)
    filtered = res.unwrap()

    # Cell 0 (5 counts) and Cell 1 (>80% mito) must be removed
    assert "cell_0" not in filtered.obs_names
    assert "cell_1" not in filtered.obs_names
    assert filtered.n_obs < sample_raw_adata.n_obs

    # GENE_0 (expressed in only 1 cell) must be filtered out
    assert "GENE_0" not in filtered.var_names
    assert filtered.n_vars < sample_raw_adata.n_vars

    # Metadata stored in .uns
    assert "qc_spec" in filtered.uns
    assert filtered.uns["qc_spec"]["n_cells_post_qc"] == filtered.n_obs


def test_standardize_processed_layers(sample_raw_adata: ad.AnnData) -> None:
    """Verify multi-layer standard: raw counts, linear 10^6 TPM, and log1p."""
    res = standardize_processed_layers(sample_raw_adata, target_sum=1e6)
    assert isinstance(res, Success)
    adata = res.unwrap()

    assert "counts" in adata.layers
    assert "tpm" in adata.layers
    assert "log1p" in adata.layers
    assert adata.raw is not None

    # Check TPM row sum for non-zero cell
    tpm_sums = np.asarray(adata.layers["tpm"].sum(axis=1)).ravel()
    for s in tpm_sums[1:]:  # skip cell 0 which had low count
        if s > 0:
            assert np.isclose(s, 1e6, rtol=1e-3)

    # Check that X is log1p of TPM
    tpm_entry = adata.layers["tpm"][2, 0]
    log1p_entry = adata.X[2, 0]
    assert np.isclose(log1p_entry, np.log1p(tpm_entry), rtol=1e-4)


def test_qc_dashboard_altair_and_svg_export(tmp_path: Path, sample_raw_adata: ad.AnnData) -> None:
    """Verify Altair dashboard builds and exports to valid SVG and PNG files."""
    adata = compute_qc_covariates(sample_raw_adata)
    metrics_df = extract_qc_metrics_dataframe(adata)
    spec = QualityControlSpec()

    dashboard = build_qc_dashboard(metrics_df, spec, dataset_name="TestCohort")
    assert dashboard is not None

    export_res = export_qc_plots(dashboard, tmp_path, dataset_name="TestCohort")
    assert isinstance(export_res, Success)
    svg_path, png_path = export_res.unwrap()

    assert svg_path.exists()
    assert png_path.exists()
    assert svg_path.stat().st_size > 500
    assert png_path.stat().st_size > 500

    # Verify SVG structure
    svg_content = svg_path.read_text(encoding="utf-8")
    assert "<svg" in svg_content
    assert "</svg>" in svg_content


def test_process_single_cell_dataset_end_to_end(tmp_path: Path, sample_raw_adata: ad.AnnData) -> None:
    """Verify full pipeline execution returns SingleCellProcessingResult with all artifacts."""
    spec = SingleCellProcessingSpec(
        use_adaptive_qc=False,
        qc=QualityControlSpec(
            min_counts_per_cell=20,
            min_genes_per_cell=3,
            max_pct_mitochondrial=30.0,
            min_cells_per_gene=2,
        ),
        target_sum=1e6,
        generate_plots=True,
        output_plot_dir=Some(tmp_path),
    )

    res = process_single_cell_dataset(
        sample_raw_adata,
        spec=spec,
        dataset_name="SyntheticTest",
        output_dir=tmp_path,
        display_plots=False,
    )
    assert isinstance(res, Success)
    result = res.unwrap()
    assert isinstance(result, SingleCellProcessingResult)

    assert result.n_cells_post_qc < result.n_cells_pre_qc
    assert result.pct_cells_retained > 0.0
    assert len(result.plot_paths) == 2  # SVG and PNG
    assert "counts" in result.adata.layers
    assert "tpm" in result.adata.layers
    assert "log1p" in result.adata.layers
    assert result.adata.uns["sc_processing_summary"]["dataset_name"] == "SyntheticTest"
