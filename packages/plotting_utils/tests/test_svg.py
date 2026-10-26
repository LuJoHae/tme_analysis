"""Unit and property-based tests for plotting_utils.svg module.
"""

from __future__ import annotations

import base64
from pathlib import Path
import xml.etree.ElementTree as ET
from hypothesis import HealthCheck, given, settings, strategies as st
from returns.result import Success
from plotting_utils.svg import (
    INKSCAPE_NS,
    XLINK_NS,
    count_svg_datapoints,
    rasterize_svg_datapoints,
    svg_badge,
    svg_chip,
    svg_circle,
    svg_line,
    svg_rect,
    svg_text,
    verify_and_rasterize_svg,
    wrap_inkscape_svg,
)


def test_svg_primitives() -> None:
    rect = svg_rect(10, 20, 100, 50, fill="#FFFFFF", stroke="#CBD5E1")
    assert '<rect x="10" y="20" width="100" height="50" fill="#FFFFFF" stroke="#CBD5E1"' in rect

    text = svg_text("Test Label", 15, 25, font_size=12.0)
    assert "Test Label</text>" in text
    assert 'font-size="12.0"' in text

    line = svg_line(0, 0, 100, 100, stroke="#CBD5E1")
    assert '<line x1="0" y1="0" x2="100" y2="100"' in line

    circle = svg_circle(50, 50, 10, fill="#0F172A")
    assert '<circle cx="50" cy="50" r="10" fill="#0F172A"' in circle

    badge = svg_badge(20, 20, "a")
    assert ">a</text>" in badge

    chip = svg_chip(10, 10, 200, 40, "PARAM TITLE", "Sub value")
    assert "PARAM TITLE" in chip
    assert "Sub value" in chip


def test_wrap_inkscape_svg_and_verification(tmp_path: Path) -> None:
    layers = [
        ("layer-01-header", "01_Header", svg_text("Header", 20, 30)),
        ("layer-02-body", "02_Body", svg_rect(20, 50, 200, 100)),
    ]
    svg_str = wrap_inkscape_svg(layers, width=800, height=600, title="Test Figure")
    assert 'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"' in svg_str
    assert 'id="layer-01-header"' in svg_str
    assert 'id="layer-02-body"' in svg_str

    out_svg = tmp_path / "test_fig.svg"
    out_png = tmp_path / "test_fig.png"
    res = verify_and_rasterize_svg(svg_str, out_svg, out_png, width=800, height=600, min_layers=2)
    assert isinstance(res, Success)
    assert out_svg.exists()

    tree = ET.parse(out_svg)
    root = tree.getroot()
    layer_elems = [
        e for e in root.iter()
        if e.attrib.get(f"{{{INKSCAPE_NS}}}groupmode") == "layer"
    ]
    assert len(layer_elems) == 2


def test_verify_and_rasterize_svg_auto_rasterizes_over_1000_points(tmp_path: Path) -> None:
    points_xml = "\n".join(
        svg_circle(20.0 + (i % 50) * 4.0, 40.0 + (i // 50) * 4.0, 2.0, fill="#0072B2")
        for i in range(1200)
    )
    layers = [
        ("layer-01-header", "01_Header", svg_text("High Density Scatter", 20, 25)),
        ("layer-02-points", "02_Points", f'<g id="scatter-group">{points_xml}</g>'),
    ]
    svg_str = wrap_inkscape_svg(layers, width=400, height=300, title="Large Scatter")

    out_svg = tmp_path / "large_scatter.svg"
    out_png = tmp_path / "large_scatter.png"
    res = verify_and_rasterize_svg(svg_str, out_svg, out_png, width=400, height=300, min_layers=2)
    assert isinstance(res, Success)

    saved_svg = out_svg.read_text(encoding="utf-8")
    root = ET.fromstring(saved_svg)
    circles = [e for e in root.iter() if e.tag.endswith("circle")]
    images = [e for e in root.iter() if e.tag.endswith("image")]
    texts = [e for e in root.iter() if e.tag.endswith("text")]

    assert len(circles) == 0
    assert len(images) == 1
    assert len(texts) == 1
    assert texts[0].text == "High Density Scatter"

    href = images[0].attrib.get(f"{{{XLINK_NS}}}href", "")
    assert href.startswith("data:image/png;base64,")
    png_bytes = base64.b64decode(href.split(",", 1)[1])
    assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")


@settings(
    max_examples=12,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    n_points=st.integers(min_value=10, max_value=75),
    threshold=st.integers(min_value=20, max_value=50),
)
def test_hypothesis_datapoint_rasterization_threshold_invariant(
    n_points: int,
    threshold: int,
) -> None:
    circles_xml = "\n".join(
        svg_circle(10.0 + (i % 10) * 8.0, 15.0 + (i // 10) * 8.0, 2.5, fill="#D55E00")
        for i in range(n_points)
    )
    raw_svg = wrap_inkscape_svg(
        [
            ("layer-01-title", "01_Title", svg_text("Axis Title", 15, 12)),
            ("layer-02-marks", "02_Marks", f'<g class="mark-symbol role-mark">{circles_xml}</g>'),
        ],
        width=200,
        height=200,
    )

    cnt_res = count_svg_datapoints(raw_svg, min_group_points=5)
    assert isinstance(cnt_res, Success)
    assert cnt_res.unwrap() == n_points

    rast_res = rasterize_svg_datapoints(
        raw_svg,
        max_points=threshold,
        min_group_points=5,
        scale=2.0,
    )
    assert isinstance(rast_res, Success)
    out_svg = rast_res.unwrap()
    root = ET.fromstring(out_svg)

    circles = [e for e in root.iter() if e.tag.endswith("circle")]
    images = [e for e in root.iter() if e.tag.endswith("image")]
    texts = [e for e in root.iter() if e.tag.endswith("text")]

    assert len(texts) == 1
    assert texts[0].text == "Axis Title"

    if n_points <= threshold:
        assert len(circles) == n_points
        assert len(images) == 0
    else:
        assert len(circles) == 0
        assert len(images) == 1
        assert int(images[0].attrib["width"]) >= 1
        assert int(images[0].attrib["height"]) >= 1
