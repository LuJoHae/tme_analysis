"""Unit tests for plotting_utils package and Altair SVG datapoint rasterization.
"""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET
import altair as alt
import polars as pl
from returns.result import Success
from plotting_utils import (
    COLOR_BG_WHITE,
    COHORT_METADATA,
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_PALETTE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
    chart_to_svg,
    export_altair_figure,
)


def test_plotting_utils_palettes_and_cohorts() -> None:
    assert len(COHORT_PANELS) == 9
    assert len(COHORT_METADATA) == 9
    assert len(OKABE_PALETTE) >= 8
    assert OKABE_BLUE == "#0072B2"
    assert OKABE_VERMILION == "#D55E00"


def test_apply_nature_methods_theme() -> None:
    df = pl.DataFrame({"x": [1, 2, 3], "y": [10, 20, 30]})
    chart = alt.Chart(df).mark_circle().encode(x="x:Q", y="y:Q")
    themed = apply_nature_methods_theme(chart)
    assert themed is not None
    chart_dict = themed.to_dict()
    assert "config" in chart_dict
    assert chart_dict["config"]["view"]["fill"] == COLOR_BG_WHITE


def test_export_altair_figure_small_stays_vector(tmp_path: Path) -> None:
    df = pl.DataFrame({"x": list(range(120)), "y": [i * 0.5 for i in range(120)]})
    chart = apply_nature_methods_theme(
        alt.Chart(df).mark_circle().encode(x="x:Q", y="y:Q")
    )
    svg_path, png_path = export_altair_figure(chart, tmp_path / "small_plot")
    assert svg_path.exists()
    assert png_path.exists()

    root = ET.parse(svg_path).getroot()
    images = [e for e in root.iter() if e.tag.endswith("image")]
    symbol_groups = [
        e for e in root.iter()
        if "mark-symbol" in e.attrib.get("class", "") and "role-mark" in e.attrib.get("class", "")
    ]
    assert len(images) == 0
    assert len(symbol_groups) == 1
    assert len(list(symbol_groups[0])) == 120


def test_export_altair_figure_large_rasterizes_datapoints_preserves_vector_axes(
    tmp_path: Path,
) -> None:
    n_pts = 1500
    df = pl.DataFrame({
        "x": [(i % 50) * 0.2 for i in range(n_pts)],
        "y": [(i // 50) * 0.3 for i in range(n_pts)],
    })
    scatter = alt.Chart(df).mark_circle(size=18, color=OKABE_BLUE).encode(
        x=alt.X("x:Q", title="UMAP 1"),
        y=alt.Y("y:Q", title="UMAP 2"),
    )
    ref_line = alt.Chart(pl.DataFrame({"y": [4.5]})).mark_rule(
        color=OKABE_VERMILION, strokeDash=[4, 4]
    ).encode(y="y:Q")
    layered = apply_nature_methods_theme(
        (scatter + ref_line).properties(title="Single-Cell Embedding")
    )

    svg_path, png_path = export_altair_figure(layered, tmp_path / "large_plot")
    assert svg_path.exists()
    assert png_path.exists()

    unrasterized_res = chart_to_svg(layered, max_vector_points=0)
    assert isinstance(unrasterized_res, Success)
    raw_svg_len = len(unrasterized_res.unwrap())
    rasterized_svg_text = svg_path.read_text(encoding="utf-8")

    # File size should be significantly smaller than 1,500 full vector <path> nodes
    assert len(rasterized_svg_text) < raw_svg_len * 0.5

    root = ET.fromstring(rasterized_svg_text)
    images = [e for e in root.iter() if e.tag.endswith("image")]
    texts = [e.text for e in root.iter() if e.tag.endswith("text") and e.text]
    rule_groups = [
        e for e in root.iter()
        if "mark-rule" in e.attrib.get("class", "") and "role-mark" in e.attrib.get("class", "")
    ]

    assert len(images) == 1
    assert "UMAP 1" in texts
    assert "UMAP 2" in texts
    assert "Single-Cell Embedding" in texts
    assert len(rule_groups) == 1
    assert len(list(rule_groups[0])) == 1


def test_chart_to_svg_multi_panel_cumulative_threshold() -> None:
    df1 = pl.DataFrame({"x": list(range(600)), "y": list(range(600))})
    df2 = pl.DataFrame({"x": list(range(600)), "y": [600 - i for i in range(600)]})
    p1 = alt.Chart(df1).mark_circle().encode(x="x:Q", y="y:Q")
    p2 = alt.Chart(df2).mark_circle().encode(x="x:Q", y="y:Q")
    grid = apply_nature_methods_theme(alt.hconcat(p1, p2))

    res = chart_to_svg(grid, max_vector_points=1000)
    assert isinstance(res, Success)
    root = ET.fromstring(res.unwrap())
    images = [e for e in root.iter() if e.tag.endswith("image")]
    assert len(images) == 2
