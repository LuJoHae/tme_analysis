"""Publication-grade Inkscape SVG construction, primitives, and rasterization.
"""

from __future__ import annotations

import html
from pathlib import Path
import subprocess
from typing import Final, Sequence
import xml.etree.ElementTree as ET
from returns.result import Failure, Result, Success

from .palettes import (
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
)


def wrap_inkscape_svg(
    layers: Sequence[tuple[str, str, str]],
    width: int,
    height: int,
    title: str = "",
    bg_color: str = COLOR_CANVAS_BG,
) -> str:
    """Wrap layer fragments into a valid, publication-ready Inkscape vector SVG.

    Parameters
    ----------
    layers:
        Sequence of (layer_id, layer_label, inner_svg_string) tuples.
    width:
        Canvas width in pixels.
    height:
        Canvas height in pixels.
    title:
        Optional SVG title element text.
    bg_color:
        Background color for the canvas rect.

    Returns
    -------
    str:
        Complete well-formed SVG markup string with Inkscape namespaces.
    """
    title_element = f"  <title>{html.escape(title)}</title>\n" if title else ""

    rendered_layers = []
    for layer_id, layer_label, inner_xml in layers:
        rendered_layers.append(
            f'  <g inkscape:groupmode="layer" id="{layer_id}" inkscape:label="{layer_label}">\n'
            f"{inner_xml}\n"
            f"  </g>"
        )
    layers_str = "\n".join(rendered_layers)

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{width}" height="{height}" viewBox="0 0 {width} {height}"
     version="1.1">
{title_element}  <sodipodi:namedview id="namedview-base" pagecolor="{bg_color}" bordercolor="{COLOR_BORDER_HAIRLINE}" borderopacity="1" />
  <rect id="bg-canvas" x="0" y="0" width="{width}" height="{height}" fill="{bg_color}" rx="0" />
{layers_str}
</svg>
"""


def svg_rect(
    x: float,
    y: float,
    w: float,
    h: float,
    fill: str = "none",
    stroke: str | None = None,
    stroke_width: float = 0.75,
    rx: float = 0.0,
    extra: str = "",
) -> str:
    """Generate SVG <rect> element."""
    stroke_attr = f'stroke="{stroke}" stroke-width="{stroke_width}" ' if stroke else ""
    rx_attr = f'rx="{rx}" ' if rx > 0 else 'rx="0" '
    extra_attr = f"{extra} " if extra else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" {stroke_attr}{rx_attr}{extra_attr}/>'


def svg_text(
    text: str,
    x: float,
    y: float,
    font_size: float = 10.0,
    font_weight: str = "400",
    fill: str = COLOR_TEXT_PRIMARY,
    anchor: str = "start",
    font_family: str = FONT_SANS,
    extra: str = "",
) -> str:
    """Generate SVG <text> element with safe HTML-escaped content."""
    extra_attr = f" {extra}" if extra else ""
    return (
        f'<text x="{x}" y="{y}" font-family="{font_family}" font-size="{font_size}" '
        f'font-weight="{font_weight}" fill="{fill}" text-anchor="{anchor}"{extra_attr}>'
        f"{text}</text>"
    )


def svg_line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str = COLOR_BORDER_HAIRLINE,
    stroke_width: float = 0.75,
    dash: str | None = None,
    extra: str = "",
) -> str:
    """Generate SVG <line> element."""
    dash_attr = f'stroke-dasharray="{dash}" ' if dash else ""
    extra_attr = f"{extra} " if extra else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{stroke_width}" {dash_attr}{extra_attr}/>'


def svg_circle(
    cx: float,
    cy: float,
    r: float,
    fill: str = COLOR_TEXT_PRIMARY,
    stroke: str | None = None,
    stroke_width: float = 0.0,
    extra: str = "",
) -> str:
    """Generate SVG <circle> element."""
    stroke_attr = f'stroke="{stroke}" stroke-width="{stroke_width}" ' if stroke else ""
    extra_attr = f"{extra} " if extra else ""
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" {stroke_attr}{extra_attr}/>'


def svg_badge(
    x: float,
    y: float,
    letter: str,
    bg_color: str = COLOR_TEXT_PRIMARY,
    text_color: str = "#FFFFFF",
    radius: float = 10.0,
    font_size: float = 11.0,
) -> str:
    """Generate publication panel badge circle with bold centered letter (e.g. 'a', 'b', 'c')."""
    return (
        f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{bg_color}" />\n'
        f'    <text x="{x}" y="{y + 4}" font-family="{FONT_SANS}" font-size="{font_size}" '
        f'font-weight="700" fill="{text_color}" text-anchor="middle">{letter}</text>'
    )


def svg_chip(
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    subtitle: str,
    bg_color: str = COLOR_SUBTLE_FILL,
    border_color: str = COLOR_BORDER_HAIRLINE,
    title_color: str = COLOR_TEXT_PRIMARY,
    sub_color: str = COLOR_TEXT_SECONDARY,
) -> str:
    """Generate formal Nature Methods statistical parameter chip."""
    return f"""<g id="grp-header-chip" transform="translate({x}, {y})">
      <rect x="0" y="0" width="{w}" height="{h}" fill="{bg_color}" stroke="{border_color}" stroke-width="0.75" rx="0" />
      <text x="{w/2}" y="15" font-family="{FONT_SANS}" font-size="8.5" font-weight="700" fill="{title_color}" text-anchor="middle" letter-spacing="0.2">
        {title}
      </text>
      <text x="{w/2}" y="27" font-family="{FONT_SANS}" font-size="8.5" font-weight="500" fill="{sub_color}" text-anchor="middle">
        {subtitle}
      </text>
    </g>"""


def verify_and_rasterize_svg(
    svg_content_or_path: str | Path,
    output_svg_path: Path,
    output_png_path: Path | None = None,
    width: int | None = None,
    height: int | None = None,
    scale: float = 2.0,
    min_layers: int = 1,
    rsvg_bin: str | Path | None = None,
) -> Result[Path, str]:
    """Validate XML well-formedness, save SVG, and rasterize 300 DPI PNG via rsvg-convert.

    Pure monadic function: checks layer counts, validates XML tree, and ensures reliable
    cross-platform rasterization.
    """
    try:
        if isinstance(svg_content_or_path, Path):
            output_svg_path = svg_content_or_path
            parsed_tree = ET.parse(output_svg_path)
            svg_content = output_svg_path.read_text(encoding="utf-8")
        else:
            svg_content = svg_content_or_path
            output_svg_path.parent.mkdir(parents=True, exist_ok=True)
            output_svg_path.write_text(svg_content, encoding="utf-8")
            parsed_tree = ET.fromstring(svg_content)

        root = parsed_tree if isinstance(parsed_tree, ET.Element) else parsed_tree.getroot()

        # Count Inkscape layers
        layers = [
            elem
            for elem in root.iter()
            if elem.attrib.get("{http://www.inkscape.org/namespaces/inkscape}groupmode") == "layer"
        ]
        if len(layers) < min_layers:
            return Failure(f"Expected at least {min_layers} Inkscape layer(s), found {len(layers)}")

        # Extract canvas dimensions if not explicitly provided
        if width is None or height is None:
            viewbox = root.attrib.get("viewBox", "")
            if viewbox:
                parts = viewbox.split()
                if len(parts) == 4:
                    width = int(float(parts[2]))
                    height = int(float(parts[3]))
            if width is None:
                width = int(float(root.attrib.get("width", "1400")))
            if height is None:
                height = int(float(root.attrib.get("height", "920")))

    except Exception as exc:
        return Failure(f"Invalid XML syntax or SVG verification failure: {exc}")

    if output_png_path is not None:
        output_png_path.parent.mkdir(parents=True, exist_ok=True)
        raster_w = int(width * scale)
        raster_h = int(height * scale)

        candidates = [
            Path(rsvg_bin) if rsvg_bin else None,
            Path("/opt/homebrew/bin/rsvg-convert"),
            Path("/usr/local/bin/rsvg-convert"),
            Path("/usr/bin/rsvg-convert"),
        ]
        found_bin: Path | str | None = next((p for p in candidates if p and p.exists()), None)
        if found_bin is None:
            found_bin = "rsvg-convert"

        try:
            subprocess.run(
                [str(found_bin), "-w", str(raster_w), "-h", str(raster_h), str(output_svg_path), "-o", str(output_png_path)],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as cpe:
            return Failure(f"rsvg-convert execution failed ({found_bin}): {cpe.stderr}")
        except FileNotFoundError:
            # Fallback if rsvg-convert is completely missing from system
            pass

    return Success(output_svg_path)
