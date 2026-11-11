"""
Unit tests for the 3x3 Cross-Cohort Deconvolution Scatter Grid visualization.
Tests pure functions for data loading, chart assembly, and output file integrity.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import altair as alt
import polars as pl
import pytest
from returns.result import Success

import importlib

_stage4a5 = importlib.import_module(
    "scripts.manuscript_pipeline.4a5_plot_cross_cohort_deconv_scatter_grid"
)
COHORT_ORDER_AND_METADATA = _stage4a5.COHORT_ORDER_AND_METADATA
CohortPanelData = _stage4a5.CohortPanelData
PlotConfig = _stage4a5.PlotConfig
assemble_3x3_scatter_grid = _stage4a5.assemble_3x3_scatter_grid
build_cohort_scatter_panel = _stage4a5.build_cohort_scatter_panel
load_cohort_panel_data = _stage4a5.load_cohort_panel_data


@pytest.fixture
def benchmark_dir() -> Path:
    return Path("output/manuscript_pipeline/3a_cross_cohort_pseudobulk_deconv_benchmark")


@pytest.fixture
def plot_config(benchmark_dir: Path) -> PlotConfig:
    return PlotConfig(benchmark_dir=benchmark_dir)


def test_cohort_metadata_definition() -> None:
    """Verify that all 9 cohorts are registered with single grid letter 'c' on the first panel."""
    assert len(COHORT_ORDER_AND_METADATA) == 9
    expected_letters = ["c", "", "", "", "", "", "", "", ""]
    actual_letters = [m["letter"] for m in COHORT_ORDER_AND_METADATA]
    assert actual_letters == expected_letters


def test_load_all_cohorts_panel_data(benchmark_dir: Path) -> None:
    """Verify that all 9 cohorts can load 5 sweeps of baseline concordance data."""
    if not benchmark_dir.exists():
        pytest.skip(f"Benchmark dir {benchmark_dir} not available.")

    for meta in COHORT_ORDER_AND_METADATA:
        cohort_id = meta["cohort"]
        res = load_cohort_panel_data(benchmark_dir, meta)
        assert isinstance(res, Success), f"Failed to load panel data for {cohort_id}: {res}"

        panel_data = res.unwrap()
        assert isinstance(panel_data, CohortPanelData)
        assert panel_data.cohort_id == cohort_id
        assert panel_data.panel_letter == meta["letter"]
        assert panel_data.n_sweeps == 5, f"{cohort_id} should have 5 sweeps, got {panel_data.n_sweeps}"
        assert panel_data.n_states >= 4
        assert panel_data.points_df.height >= 20  # At least 5 sweeps * 4 states

        # Verify essential columns exist in points_df
        required_cols = {"cell_state", "milo_mean_logfc", "deconv_beta", "sweep_idx"}
        assert required_cols.issubset(set(panel_data.points_df.columns))

        # Verify correlation and sign concordance bounds
        assert -1.0 <= panel_data.rho_mean <= 1.0
        assert 0.0 <= panel_data.sign_mean <= 100.0


def test_build_cohort_scatter_panel(benchmark_dir: Path, plot_config: PlotConfig) -> None:
    """Verify that a single cohort panel constructs a valid Altair LayerChart."""
    meta = COHORT_ORDER_AND_METADATA[0]  # GSE120575
    res = load_cohort_panel_data(benchmark_dir, meta)
    if not isinstance(res, Success):
        pytest.skip("Benchmark data not present.")

    panel_data = res.unwrap()
    chart = build_cohort_scatter_panel(panel_data, plot_config)
    assert chart is not None
    chart_dict = chart.to_dict()
    assert "layer" in chart_dict
    assert len(chart_dict["layer"]) == 4  # h_rule, v_rule, trend_line, c_base


def test_assemble_3x3_scatter_grid(benchmark_dir: Path, plot_config: PlotConfig) -> None:
    """Verify that 9 panels assemble correctly into a 3x3 VConcatChart and fail on invalid counts."""
    meta = COHORT_ORDER_AND_METADATA[0]
    res = load_cohort_panel_data(benchmark_dir, meta)
    if not isinstance(res, Success):
        pytest.skip("Benchmark data not present.")

    panel_data = res.unwrap()
    single_panel = build_cohort_scatter_panel(panel_data, plot_config)

    # Valid assembly
    nine_panels = [single_panel] * 9
    grid = assemble_3x3_scatter_grid(nine_panels)
    assert isinstance(grid, alt.VConcatChart)
    grid_dict = grid.to_dict()
    assert "vconcat" in grid_dict
    assert len(grid_dict["vconcat"]) == 3  # 3 rows

    # Error handling when panel count is not 9
    with pytest.raises(ValueError, match="Expected exactly 9 panels"):
        assemble_3x3_scatter_grid([single_panel] * 8)


def test_output_figures_exist_and_nonempty() -> None:
    """Verify that the generated SVG and PNG figures exist, have 130mm width, and >=10pt font sizes."""
    import re

    svg_path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
    )
    png_path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.png"
    )
    article_path = Path("article/figures/milo-deconv/nine_cohorts_deconv_scatter_grid.png")

    assert svg_path.exists(), f"Missing SVG at {svg_path}"
    assert png_path.exists(), f"Missing PNG at {png_path}"
    assert article_path.exists(), f"Missing manuscript figure at {article_path}"

    assert svg_path.stat().st_size > 10_000, "SVG file suspiciously small"
    assert png_path.stat().st_size > 50_000, "PNG file suspiciously small"
    assert article_path.stat().st_size > 50_000, "Manuscript PNG file suspiciously small"

    svg_text = svg_path.read_text(encoding="utf-8")
    assert 'width="130.0mm"' in svg_text, "SVG must have physical width of 130.0mm"

    font_sizes = [float(x) for x in re.findall(r'font-size="([0-9.]+)px"', svg_text)]
    assert len(font_sizes) > 0, "Expected font-size attributes in SVG"
    assert min(font_sizes) >= 13.3, f"Minimum font-size in 130mm SVG must be >=10pt (13.33px), got {min(font_sizes)}px"

