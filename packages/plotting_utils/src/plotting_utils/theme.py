"""Nature Methods minimalist wireframe styling for Altair charts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, TypeVar, cast
import xml.etree.ElementTree as ET
import altair as alt
from returns.result import Failure, Result, Success
import vl_convert as vlc

from .palettes import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
)
from .svg import INKSCAPE_NS, SODIPODI_NS, rasterize_svg_datapoints, verify_and_rasterize_svg

ChartT = TypeVar(
    "ChartT",
    alt.TopLevelMixin,
    alt.Chart,
    alt.HConcatChart,
    alt.VConcatChart,
    alt.LayerChart,
)


def apply_nature_methods_theme(
    chart: ChartT,
) -> ChartT:
    """Applies Nature Methods minimalist wireframe styling to an Altair chart."""
    return cast(
        ChartT,
        chart.configure(background=COLOR_BG_WHITE, font=FONT_SANS)
        .configure_view(
            strokeWidth=0.75,
            stroke=COLOR_HAIRLINE,
            fill=COLOR_BG_WHITE,
        )
        .configure_axis(
            labelFont=FONT_SANS,
            titleFont=FONT_SANS,
            labelFontSize=10,
            titleFontSize=11,
            titleFontWeight="bold",
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
            gridColor=COLOR_LIGHT_GREY,
            gridWidth=0.5,
            domainColor=COLOR_HAIRLINE,
            domainWidth=0.8,
            tickColor=COLOR_HAIRLINE,
            tickWidth=0.8,
        )
        .configure_header(
            labelFont=FONT_SANS,
            titleFont=FONT_SANS,
            labelFontSize=11,
            labelFontWeight="bold",
            labelColor=COLOR_TEXT_PRIMARY,
            titleColor=COLOR_TEXT_PRIMARY,
        )
        .configure_legend(
            labelFont=FONT_SANS,
            titleFont=FONT_SANS,
            labelFontSize=10,
            titleFontSize=10.5,
            titleFontWeight="bold",
            labelLimit=260,
            titleLimit=260,
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
        )
        .configure_title(
            font=FONT_SANS,
            color=COLOR_TEXT_PRIMARY,
            fontSize=12.5,
            fontWeight="bold",
            anchor="start",
        ),
    )


def create_nature_axis(
    title: str,
    values: list[float] | None = None,
    title_font_size: int = 11,
    label_font_size: int = 10,
    format_str: str | None = None,
) -> alt.Axis:
    """Create an Altair axis adhering strictly to Nature Methods styling."""
    kwargs: dict[str, Any] = {
        "title": title,
        "titleFontSize": title_font_size,
        "labelFontSize": label_font_size,
        "titleFont": "sans-serif",
        "labelFont": "sans-serif",
        "titleColor": COLOR_TEXT_PRIMARY,
        "labelColor": COLOR_TEXT_SECONDARY,
        "tickColor": COLOR_HAIRLINE,
        "domainColor": COLOR_HAIRLINE,
        "gridColor": COLOR_LIGHT_GREY,
    }
    if values is not None:
        kwargs["values"] = values
    if format_str is not None:
        kwargs["format"] = format_str
    return alt.Axis(**kwargs)


def chart_to_svg(
    chart: alt.Chart | alt.HConcatChart | alt.VConcatChart | alt.LayerChart,
    max_vector_points: int = 1000,
    scale: float = 3.0,
) -> Result[str, str]:
    """Purely compile an Altair chart to SVG, rasterizing datapoints when exceeding `max_vector_points`."""
    try:
        chart_dict = chart.to_dict()
        raw_svg = str(vlc.vegalite_to_svg(chart_dict))
    except Exception as exc:
        return Failure(f"Failed to compile Altair chart to SVG: {exc}")

    return rasterize_svg_datapoints(raw_svg, max_points=max_vector_points, scale=scale)


def export_altair_figure(
    chart: alt.Chart | alt.HConcatChart | alt.VConcatChart | alt.LayerChart,
    base_path: Path,
    scale: float = 2.5,
    max_vector_points: int = 1000,
) -> tuple[Path, Path]:
    """
    Purely export an Altair chart to publication vector SVG and 300 DPI raster PNG.

    Whenever the chart contains more than `max_vector_points` datapoints (default 1,000),
    dense datapoint mark groups in the SVG are automatically rasterized into embedded
    high-DPI transparent PNG images while axes, labels, titles, and legends remain vector.

    Parameters:
    - chart: Compiled Altair chart
    - base_path: Path without extension (or with .svg/.png, which will be stripped)
    - scale: Scaling factor for PNG rasterization (default 2.5 for ~300 DPI)
    - max_vector_points: Maximum total datapoints to keep as vector before rasterizing (default 1,000)

    Returns:
    - Tuple of (svg_path, png_path)
    """
    clean_stem = Path(base_path).parent / Path(base_path).stem
    svg_path = clean_stem.with_suffix(".svg")
    png_path = clean_stem.with_suffix(".png")
    svg_path.parent.mkdir(parents=True, exist_ok=True)

    chart_dict = chart.to_dict()
    raw_svg_str = str(vlc.vegalite_to_svg(chart_dict))
    match rasterize_svg_datapoints(
        raw_svg_str,
        max_points=max_vector_points,
        scale=max(scale, 3.0),
    ):
        case Success(processed_svg):
            svg_str = processed_svg
        case Failure(_):
            svg_str = raw_svg_str

    svg_path.write_text(svg_str, encoding="utf-8")

    png_bytes = vlc.vegalite_to_png(chart_dict, scale=scale)
    png_path.write_bytes(png_bytes)

    return svg_path, png_path


def _wrap_altair_svg_in_inkscape_layers(raw_svg: str) -> str:
    """Structure compiled Vega-Lite SVG elements into native Inkscape layers (`01_Background`, `02_Panels_and_Marks`)."""
    svg_ns = "http://www.w3.org/2000/svg"
    ET.register_namespace("", svg_ns)
    ET.register_namespace("inkscape", INKSCAPE_NS)
    ET.register_namespace("sodipodi", SODIPODI_NS)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

    root = ET.fromstring(raw_svg)
    root.set(f"{{{INKSCAPE_NS}}}version", "1.3")

    existing_children = list(root)
    bg_layer = ET.Element(
        f"{{{svg_ns}}}g",
        {
            f"{{{INKSCAPE_NS}}}groupmode": "layer",
            "id": "layer-01",
            f"{{{INKSCAPE_NS}}}label": "01_Background",
        },
    )
    panel_layer = ET.Element(
        f"{{{svg_ns}}}g",
        {
            f"{{{INKSCAPE_NS}}}groupmode": "layer",
            "id": "layer-02",
            f"{{{INKSCAPE_NS}}}label": "02_Panels_and_Marks",
        },
    )

    for child in existing_children:
        root.remove(child)
        tag_local = child.tag.split("}")[-1]
        if tag_local in ("defs", "namedview", "title", "desc"):
            root.append(child)
        elif tag_local == "rect" and child.attrib.get("class", "").startswith("background"):
            bg_layer.append(child)
        elif tag_local == "rect" and child.attrib.get("x", "0") == "0" and child.attrib.get("y", "0") == "0" and len(bg_layer) == 0:
            bg_layer.append(child)
        else:
            panel_layer.append(child)

    if len(bg_layer) == 0:
        w = root.attrib.get("width", "100%")
        h = root.attrib.get("height", "100%")
        bg_layer.append(
            ET.Element(
                f"{{{svg_ns}}}rect",
                {"x": "0", "y": "0", "width": str(w), "height": str(h), "fill": COLOR_BG_WHITE},
            )
        )

    root.append(bg_layer)
    root.append(panel_layer)
    return ET.tostring(root, encoding="unicode")


def export_altair_inkscape_svg(
    chart: alt.TopLevelMixin,
    base_path: Path,
    scale: float = 2.5,
    max_vector_points: int = 1500,
) -> Result[tuple[Path, Path], str]:
    """Compile an Altair chart into an Inkscape-layered SVG (`>=2` layers) and 300 DPI PNG via `verify_and_rasterize_svg`."""
    clean_stem = Path(base_path).parent / Path(base_path).stem
    svg_path = clean_stem.with_suffix(".svg")
    png_path = clean_stem.with_suffix(".png")

    raw_svg = str(vlc.vegalite_to_svg(chart.to_dict()))
    layered_svg = _wrap_altair_svg_in_inkscape_layers(raw_svg)
    verify_res = verify_and_rasterize_svg(
        layered_svg,
        output_svg_path=svg_path,
        output_png_path=png_path,
        scale=scale,
        min_layers=2,
        max_vector_points=max_vector_points,
        point_raster_scale=max(scale, 3.0),
    )
    if isinstance(verify_res, Failure):
        return Failure(verify_res.failure())
    return Success((svg_path, png_path))
