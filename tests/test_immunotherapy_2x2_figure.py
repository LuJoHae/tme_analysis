"""
Unit tests for the Unified 2x2 Master Immunotherapy Figure (180mm full-width).
Validates physical dimensions, typography (>=10pt), panel lettering (a, b, c, d),
and cohort symmetry across all four panels.
"""

from __future__ import annotations

from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

import polars as pl
import pytest

import importlib

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

_stage4a3 = importlib.import_module(
    "scripts.manuscript_pipeline.4a3_plot_cross_cohort_umaps_grid"
)
COHORT_ORDER_AND_METADATA = _stage4a3.COHORT_ORDER_AND_METADATA
MIN_FONT_PX = _stage4a3.MIN_FONT_PX
UMAP_WIDTH_MM = _stage4a3.TARGET_WIDTH_MM
UmapGridConfig = _stage4a3.UmapGridConfig

_stage4a4 = importlib.import_module(
    "scripts.manuscript_pipeline.4a4_plot_cross_cohort_da_bars"
)
DaBarConfig = _stage4a4.DaBarConfig
BAR_WIDTH_MM = _stage4a4.TARGET_WIDTH_MM
load_all_cohort_da_stats = _stage4a4.load_all_cohort_da_stats

_stage4a7 = importlib.import_module(
    "scripts.manuscript_pipeline.4a7_assemble_immunotherapy_2x2_figure"
)
AssemblerConfig = _stage4a7.AssemblerConfig
TOTAL_WIDTH_MM = _stage4a7.TOTAL_WIDTH_MM
VIEWBOX_TOTAL_WIDTH = _stage4a7.VIEWBOX_TOTAL_WIDTH
VIEWBOX_TOTAL_HEIGHT = _stage4a7.VIEWBOX_TOTAL_HEIGHT


@pytest.fixture
def master_svg_path() -> Path:
    return Path("results/milopy_analysis/immunotherapy_response_2x2_master.svg")


@pytest.fixture
def master_png_path() -> Path:
    return Path("results/milopy_analysis/immunotherapy_response_2x2_master.png")


@pytest.fixture
def row1_composite_svg_path() -> Path:
    return Path("results/milopy_analysis/nine_cohorts_single_cell_composite.svg")


@pytest.fixture
def panel_a_svg_path() -> Path:
    return Path("results/milopy_analysis/nine_cohorts_umaps_logfc_grid.svg")


@pytest.fixture
def panel_b_svg_path() -> Path:
    return Path("results/milopy_analysis/nine_cohorts_da_bars.svg")


@pytest.fixture
def panel_c_svg_path() -> Path:
    return Path("results/pseudobulk_deconv_benchmark/nine_cohorts_deconv_scatter_grid.svg")


@pytest.fixture
def panel_d_svg_path() -> Path:
    return Path("results/pseudobulk_deconv_benchmark/nine_cohorts_deconv_correlations_bar.svg")


def test_cohort_order_definition() -> None:
    """Verify that all 9 cohorts are correctly registered in canonical order."""
    assert len(COHORT_ORDER_AND_METADATA) == 9
    expected_cohorts = [
        "GSE120575",
        "CELLxGENE_7b20c613",
        "CELLxGENE_05a8c945",
        "CELLxGENE_6f9de485",
        "GSE207422",
        "GSE243013",
        "GSE233203",
        "GSE200996",
        "GSE316195",
    ]
    actual_cohorts = [m["cohort"] for m in COHORT_ORDER_AND_METADATA]
    assert actual_cohorts == expected_cohorts


def test_panel_a_svg_properties(panel_a_svg_path: Path) -> None:
    """Verify Panel a physical dimensions and strict font size compliance."""
    if not panel_a_svg_path.exists():
        pytest.skip(f"Panel a SVG {panel_a_svg_path} not found.")

    tree = ET.parse(panel_a_svg_path)
    root = tree.getroot()
    assert root.attrib["width"] == f"{UMAP_WIDTH_MM:.1f}mm"
    assert "123.31mm" in root.attrib["height"]

    content = panel_a_svg_path.read_text(encoding="utf-8")
    sizes = [float(m) for m in re.findall(r'font-size="([\d\.]+)px"', content)]
    assert len(sizes) > 0
    assert min(sizes) >= 13.33, f"Font size violation in Panel a: {min(sizes)}px < 13.33px"
    assert "(a)  Melanoma" in content


def test_panel_b_svg_properties(panel_b_svg_path: Path) -> None:
    """Verify Panel b physical dimensions and strict font size compliance."""
    if not panel_b_svg_path.exists():
        pytest.skip(f"Panel b SVG {panel_b_svg_path} not found.")

    tree = ET.parse(panel_b_svg_path)
    root = tree.getroot()
    assert root.attrib["width"] == f"{BAR_WIDTH_MM:.1f}mm"
    assert "123.31mm" in root.attrib["height"]

    content = panel_b_svg_path.read_text(encoding="utf-8")
    sizes = [float(m) for m in re.findall(r'font-size="([\d\.]+)px"', content)]
    assert len(sizes) > 0
    assert min(sizes) >= 13.33, f"Font size violation in Panel b: {min(sizes)}px < 13.33px"
    assert "(b)  Differential" in content


def test_master_2x2_figure_dimensions_and_fonts(master_svg_path: Path) -> None:
    """Verify the combined 2x2 master figure has 180mm width, valid viewBox, and >=10pt font."""
    if not master_svg_path.exists():
        pytest.skip(f"Master SVG {master_svg_path} not found.")

    tree = ET.parse(master_svg_path)
    root = tree.getroot()
    assert root.attrib["width"] == f"{TOTAL_WIDTH_MM:.1f}mm"
    assert root.attrib["viewBox"] == f"0 0 {VIEWBOX_TOTAL_WIDTH} {VIEWBOX_TOTAL_HEIGHT}"

    content = master_svg_path.read_text(encoding="utf-8")
    sizes = [float(m) for m in re.findall(r'font-size="([\d\.]+)px"', content)]
    assert len(sizes) > 0
    assert min(sizes) >= 13.33, f"Font size violation in Master SVG: min size {min(sizes)}px < 13.33px (10pt)"

    # Verify panel tags (a), (b), (c), (d) are all present
    assert "(a)  Melanoma" in content, "Missing Panel (a) tag"
    assert "(b)  Differential" in content, "Missing Panel (b) tag"
    assert "(c)  Melanoma" in content, "Missing Panel (c) tag"
    assert "(d)  Concordance" in content, "Missing Panel (d) tag"

    # Verify all 9 ground-truth author labels are present and no misattributed authors remain
    for author in ("Sade-Feldman", "Gondal", "Marteau", "Yan", "Hu", "Liu", "Hong", "Luoma", "Raufi"):
        assert author in content, f"Missing verified author '{author}' in Master SVG"
    for old_wrong in ("Jerby-Arnon", "Bassez", "Bockorny"):
        assert old_wrong not in content, f"Found deprecated misattributed author '{old_wrong}' in Master SVG"


def test_master_png_rasterization(master_png_path: Path) -> None:
    """Verify master PNG raster is non-empty and has high resolution."""
    if not master_png_path.exists():
        pytest.skip(f"Master PNG {master_png_path} not found.")

    stat = master_png_path.stat()
    assert stat.st_size > 100_000, f"Master PNG unexpectedly small ({stat.st_size} bytes)"
