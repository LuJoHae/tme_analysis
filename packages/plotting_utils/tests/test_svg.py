"""Unit tests for plotting_utils.svg module.
"""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET
from returns.result import Success
from plotting_utils.svg import (
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
        if e.attrib.get("{http://www.inkscape.org/namespaces/inkscape}groupmode") == "layer"
    ]
    assert len(layer_elems) == 2
