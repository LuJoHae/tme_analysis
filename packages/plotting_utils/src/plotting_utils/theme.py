"""Nature Methods minimalist wireframe styling for Altair charts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, TypeVar, cast
import altair as alt

ChartT = TypeVar(
    "ChartT",
    alt.Chart,
    alt.HConcatChart,
    alt.VConcatChart,
    alt.LayerChart,
)
from returns.result import Failure, Result, Success
import vl_convert as vlc

from .palettes import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
)
from .svg import rasterize_svg_datapoints


def apply_nature_methods_theme(
    chart: ChartT,
) -> ChartT:
    """Applies Nature Methods minimalist wireframe styling to an Altair chart."""
    return cast(
        ChartT,
        chart.configure_view(
            strokeWidth=0.75,
            stroke=COLOR_HAIRLINE,
            fill=COLOR_BG_WHITE,
        )
        .configure_axis(
            labelFont="sans-serif",
            titleFont="sans-serif",
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
            gridColor=COLOR_LIGHT_GREY,
            domainColor=COLOR_HAIRLINE,
            tickColor=COLOR_HAIRLINE,
        )
        .configure_legend(
            labelFont="sans-serif",
            titleFont="sans-serif",
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
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
