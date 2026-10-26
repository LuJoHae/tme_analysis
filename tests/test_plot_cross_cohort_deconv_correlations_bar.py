"""
Unit tests for the Cross-Cohort Deconvolution Correlation Bar Plot and 180mm Composite Figure.
Tests pure functions for data loading, chart construction, physical mm dimensions, and font sizing.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import altair as alt
import polars as pl
import pytest
from returns.result import Success

import importlib

_stage4a6 = importlib.import_module(
    "scripts.manuscript_pipeline.4a6_plot_cross_cohort_deconv_correlations_bar"
)
BAR_WIDTH_MM = _stage4a6.BAR_WIDTH_MM
COHORT_ORDER_AND_METADATA = _stage4a6.COHORT_ORDER_AND_METADATA
COMPOSITE_WIDTH_MM = _stage4a6.COMPOSITE_WIDTH_MM
FIGURE_HEIGHT_MM = _stage4a6.FIGURE_HEIGHT_MM
MIN_FONT_PT = _stage4a6.MIN_FONT_PT
BarPlotConfig = _stage4a6.BarPlotConfig
build_correlation_bar_chart = _stage4a6.build_correlation_bar_chart
enforce_svg_physical_dimensions = _stage4a6.enforce_svg_physical_dimensions
load_cohort_correlations = _stage4a6.load_cohort_correlations
stitch_composite_grid_and_bars = _stage4a6.stitch_composite_grid_and_bars


@pytest.fixture
def benchmark_dir() -> Path:
    return Path("output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark")


@pytest.fixture
def bar_config(benchmark_dir: Path) -> BarPlotConfig:
    return BarPlotConfig(benchmark_dir=benchmark_dir)


def test_load_cohort_correlations(benchmark_dir: Path) -> None:
    """Verify loading correlation metrics across all 9 cohorts."""
    if not benchmark_dir.exists():
        pytest.skip(f"Benchmark dir {benchmark_dir} not available.")

    res = load_cohort_correlations(benchmark_dir)
    assert isinstance(res, Success), f"Failed to load correlations: {res}"

    df = res.unwrap()
    assert isinstance(df, pl.DataFrame)
    assert df.height == 9

    expected_letters = ["c", "", "", "", "", "", "", "", ""]
    assert df["letter"].to_list() == expected_letters
    assert "Sade-Feldman" in df["label"].to_list()
    assert "Raufi" in df["label"].to_list()

    for row in df.iter_rows(named=True):
        assert 0.0 <= row["rho_mean"] <= 1.0, f"Invalid rho for {row['cohort']}: {row['rho_mean']}"
        assert 0.0 <= row["rho_std"] <= 0.5
        assert row["rho_min"] <= row["rho_max"]


def test_build_correlation_bar_chart(benchmark_dir: Path, bar_config: BarPlotConfig) -> None:
    """Verify construction of Altair layered chart specification."""
    res = load_cohort_correlations(benchmark_dir)
    if not isinstance(res, Success):
        pytest.skip("Benchmark data not present.")

    df = res.unwrap()
    chart = build_correlation_bar_chart(df, bar_config)
    assert chart is not None
    chart_dict = chart.to_dict()
    assert "layer" in chart_dict
    assert len(chart_dict["layer"]) == 3  # bars + error whiskers + cross-cohort mean rule


def test_standalone_bar_svg_dimensions_and_fonts() -> None:
    """Verify standalone bar SVG has exact 50.0mm width, 179.75mm height, and font sizes >= 10pt."""
    svg_path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.svg"
    )
    assert svg_path.exists(), f"Missing SVG at {svg_path}"

    svg_text = svg_path.read_text(encoding="utf-8")
    assert f'width="{BAR_WIDTH_MM:.1f}mm"' in svg_text
    assert f'height="{FIGURE_HEIGHT_MM:.2f}mm"' in svg_text
    assert 'viewBox="0 0 188 461"' in svg_text

    # Verify minimum font size is strictly >= 10 pt (13.33 px at 96 CSS DPI)
    font_sizes = [float(x) for x in re.findall(r'font-size="([0-9.]+)px"', svg_text)]
    assert len(font_sizes) > 0, "No font-size attributes found in SVG"
    assert min(font_sizes) >= 13.3, f"Minimum font-size in 50mm SVG must be >=10pt (13.33px), got {min(font_sizes)}px"


def test_composite_svg_dimensions_and_fonts() -> None:
    """Verify combined 180mm composite SVG matches required 2-column layout and font sizes."""
    composite_svg_path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.svg"
    )
    assert composite_svg_path.exists(), f"Missing composite SVG at {composite_svg_path}"

    svg_text = composite_svg_path.read_text(encoding="utf-8")
    assert f'width="{COMPOSITE_WIDTH_MM:.1f}mm"' in svg_text
    assert f'height="{FIGURE_HEIGHT_MM:.2f}mm"' in svg_text
    assert 'viewBox="0 0 674 461"' in svg_text

    font_sizes = [float(x) for x in re.findall(r'font-size="([0-9.]+)px"', svg_text)]
    assert len(font_sizes) > 0, "No font-size attributes found in composite SVG"
    assert min(font_sizes) >= 13.3, f"Minimum font-size in composite SVG must be >=10pt (13.33px), got {min(font_sizes)}px"


def test_output_rasters_exist_and_nonempty() -> None:
    """Verify 300 DPI raster PNGs exist in results and article directories."""
    bar_png = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.png"
    )
    bar_art = Path("article/figures/milo-deconv/nine_cohorts_deconv_correlations_bar.png")
    comp_png = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.png"
    )
    comp_art = Path("article/figures/milo-deconv/nine_cohorts_deconv_composite_with_bars.png")

    for p in [bar_png, bar_art, comp_png, comp_art]:
        assert p.exists(), f"File does not exist: {p}"
        assert p.stat().st_size > 15_000, f"File {p} is suspiciously small: {p.stat().st_size} bytes"
