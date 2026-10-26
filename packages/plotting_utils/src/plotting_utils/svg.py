"""Publication-grade Inkscape SVG construction, primitives, and rasterization.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
import html
import math
from pathlib import Path
import re
import struct
import subprocess
from typing import Final, Sequence
import xml.etree.ElementTree as ET
import zlib
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import vl_convert as vlc

from .palettes import (
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
)

SVG_NS: Final[str] = "http://www.w3.org/2000/svg"
XLINK_NS: Final[str] = "http://www.w3.org/1999/xlink"
INKSCAPE_NS: Final[str] = "http://www.inkscape.org/namespaces/inkscape"
SODIPODI_NS: Final[str] = "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"

_TRANSLATE_RE: Final[re.Pattern[str]] = re.compile(
    r"translate\(\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?:[\s,]+([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?))?\s*\)"
)
_FLOAT_RE: Final[re.Pattern[str]] = re.compile(
    r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
)
_PATH_TOKEN_RE: Final[re.Pattern[str]] = re.compile(
    r"[MmLlHhVvCcSsQqTtAaZz]|[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
)
_CLIP_URL_RE: Final[re.Pattern[str]] = re.compile(r"url\(#([^)]+)\)")

_DATAPOINT_CHILD_TAGS: Final[frozenset[str]] = frozenset(
    {"path", "circle", "rect", "ellipse", "use"}
)
_INHERITED_STYLE_ATTRS: Final[tuple[str, ...]] = (
    "fill",
    "stroke",
    "stroke-width",
    "stroke-miterlimit",
    "stroke-linecap",
    "stroke-linejoin",
)


@dataclass(frozen=True)
class SvgRasterizationConfig:
    """Immutable configuration for selective SVG datapoint mark-group rasterization."""

    max_points: int = 1000
    min_group_points: int = 50
    scale: float = 3.0
    pad_px: float = 6.0


@dataclass(frozen=True)
class BoundingBox:
    """Immutable 2D integer bounding box in local group coordinates."""

    x: int
    y: int
    width: int
    height: int


def _register_svg_namespaces() -> None:
    """Register standard SVG and Inkscape XML namespaces for clean serialization."""
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    ET.register_namespace("inkscape", INKSCAPE_NS)
    ET.register_namespace("sodipodi", SODIPODI_NS)


def _strip_ns(tag: str) -> str:
    """Strip XML namespace URI prefix from an ElementTree tag string."""
    return tag.split("}", 1)[1] if "}" in tag else tag


def _parse_float(value: str | None, default: float = 0.0) -> float:
    """Purely parse a float string, returning default on missing or malformed values."""
    if not value:
        return default
    m = _FLOAT_RE.search(value)
    return float(m.group(0)) if m else default


def _parse_svg_path_bounds(d: str) -> Maybe[tuple[float, float, float, float]]:
    """Compute (min_x, min_y, max_x, max_y) visited by an SVG path `d` command string."""
    tokens = _PATH_TOKEN_RE.findall(d)
    if not tokens:
        return Nothing

    xs: list[float] = []
    ys: list[float] = []
    cur_x = 0.0
    cur_y = 0.0
    start_x = 0.0
    start_y = 0.0
    cmd = "M"
    idx = 0
    n_tokens = len(tokens)

    while idx < n_tokens:
        tok = tokens[idx]
        if tok.isalpha():
            cmd = tok
            idx += 1
            if cmd in ("Z", "z"):
                cur_x, cur_y = start_x, start_y
            continue

        match cmd:
            case "M" | "L" | "T":
                if idx + 1 >= n_tokens:
                    break
                cur_x = float(tokens[idx])
                cur_y = float(tokens[idx + 1])
                xs.append(cur_x)
                ys.append(cur_y)
                if cmd == "M":
                    start_x, start_y = cur_x, cur_y
                    cmd = "L"
                idx += 2
            case "m" | "l" | "t":
                if idx + 1 >= n_tokens:
                    break
                cur_x += float(tokens[idx])
                cur_y += float(tokens[idx + 1])
                xs.append(cur_x)
                ys.append(cur_y)
                if cmd == "m":
                    start_x, start_y = cur_x, cur_y
                    cmd = "l"
                idx += 2
            case "H":
                cur_x = float(tokens[idx])
                xs.append(cur_x)
                ys.append(cur_y)
                idx += 1
            case "h":
                cur_x += float(tokens[idx])
                xs.append(cur_x)
                ys.append(cur_y)
                idx += 1
            case "V":
                cur_y = float(tokens[idx])
                xs.append(cur_x)
                ys.append(cur_y)
                idx += 1
            case "v":
                cur_y += float(tokens[idx])
                xs.append(cur_x)
                ys.append(cur_y)
                idx += 1
            case "C":
                if idx + 5 >= n_tokens:
                    break
                xs.extend([float(tokens[idx]), float(tokens[idx + 2]), float(tokens[idx + 4])])
                ys.extend([float(tokens[idx + 1]), float(tokens[idx + 3]), float(tokens[idx + 5])])
                cur_x = float(tokens[idx + 4])
                cur_y = float(tokens[idx + 5])
                idx += 6
            case "c":
                if idx + 5 >= n_tokens:
                    break
                xs.extend([
                    cur_x + float(tokens[idx]),
                    cur_x + float(tokens[idx + 2]),
                    cur_x + float(tokens[idx + 4]),
                ])
                ys.extend([
                    cur_y + float(tokens[idx + 1]),
                    cur_y + float(tokens[idx + 3]),
                    cur_y + float(tokens[idx + 5]),
                ])
                cur_x += float(tokens[idx + 4])
                cur_y += float(tokens[idx + 5])
                idx += 6
            case "S" | "Q":
                if idx + 3 >= n_tokens:
                    break
                xs.extend([float(tokens[idx]), float(tokens[idx + 2])])
                ys.extend([float(tokens[idx + 1]), float(tokens[idx + 3])])
                cur_x = float(tokens[idx + 2])
                cur_y = float(tokens[idx + 3])
                idx += 4
            case "s" | "q":
                if idx + 3 >= n_tokens:
                    break
                xs.extend([cur_x + float(tokens[idx]), cur_x + float(tokens[idx + 2])])
                ys.extend([cur_y + float(tokens[idx + 1]), cur_y + float(tokens[idx + 3])])
                cur_x += float(tokens[idx + 2])
                cur_y += float(tokens[idx + 3])
                idx += 4
            case "A":
                if idx + 6 >= n_tokens:
                    break
                rx = abs(float(tokens[idx]))
                ry = abs(float(tokens[idx + 1]))
                end_x = float(tokens[idx + 5])
                end_y = float(tokens[idx + 6])
                xs.extend([cur_x - rx, cur_x + rx, end_x - rx, end_x + rx])
                ys.extend([cur_y - ry, cur_y + ry, end_y - ry, end_y + ry])
                cur_x, cur_y = end_x, end_y
                idx += 7
            case "a":
                if idx + 6 >= n_tokens:
                    break
                rx = abs(float(tokens[idx]))
                ry = abs(float(tokens[idx + 1]))
                end_x = cur_x + float(tokens[idx + 5])
                end_y = cur_y + float(tokens[idx + 6])
                xs.extend([cur_x - rx, cur_x + rx, end_x - rx, end_x + rx])
                ys.extend([cur_y - ry, cur_y + ry, end_y - ry, end_y + ry])
                cur_x, cur_y = end_x, end_y
                idx += 7
            case _:
                idx += 1

    return Some((min(xs), min(ys), max(xs), max(ys))) if (xs and ys) else Nothing


def _element_local_bounds(
    elem: ET.Element[str],
) -> Maybe[tuple[float, float, float, float]]:
    """Compute local coordinate bounds (min_x, min_y, max_x, max_y) of a single child mark."""
    tag = _strip_ns(elem.tag)
    sw = _parse_float(elem.attrib.get("stroke-width"), default=1.0)
    tr = elem.attrib.get("transform", "")
    m_tr = _TRANSLATE_RE.search(tr)
    tx = float(m_tr.group(1)) if m_tr else 0.0
    ty = float(m_tr.group(2)) if (m_tr and m_tr.group(2) is not None) else 0.0

    match tag:
        case "circle":
            cx = _parse_float(elem.attrib.get("cx")) + tx
            cy = _parse_float(elem.attrib.get("cy")) + ty
            r = _parse_float(elem.attrib.get("r"), default=2.0) + sw
            return Some((cx - r, cy - r, cx + r, cy + r))
        case "ellipse":
            cx = _parse_float(elem.attrib.get("cx")) + tx
            cy = _parse_float(elem.attrib.get("cy")) + ty
            rx = _parse_float(elem.attrib.get("rx"), default=2.0) + sw
            ry = _parse_float(elem.attrib.get("ry"), default=2.0) + sw
            return Some((cx - rx, cy - ry, cx + rx, cy + ry))
        case "rect" | "use":
            x = _parse_float(elem.attrib.get("x")) + tx
            y = _parse_float(elem.attrib.get("y")) + ty
            w = _parse_float(elem.attrib.get("width"), default=2.0)
            h = _parse_float(elem.attrib.get("height"), default=2.0)
            return Some((
                min(x, x + w) - sw,
                min(y, y + h) - sw,
                max(x, x + w) + sw,
                max(y, y + h) + sw,
            ))
        case "path":
            d = elem.attrib.get("d", "")
            if m_tr:
                nums = [abs(float(v)) for v in _FLOAT_RE.findall(d)]
                r = (max(nums) if nums else 2.5) + sw
                return Some((tx - r, ty - r, tx + r, ty + r))
            match _parse_svg_path_bounds(d):
                case Some((x0, y0, x1, y1)):
                    return Some((x0 + tx - sw, y0 + ty - sw, x1 + tx + sw, y1 + ty + sw))
                case _:
                    return Nothing
        case _:
            return Nothing


def _extract_clip_rects(
    root: ET.Element[str],
) -> dict[str, tuple[float, float, float, float]]:
    """Extract (min_x, min_y, max_x, max_y) for each <clipPath id='...'> in <defs>."""
    clip_map: dict[str, tuple[float, float, float, float]] = {}
    for elem in root.iter():
        if _strip_ns(elem.tag) == "clipPath":
            clip_id = elem.attrib.get("id", "")
            if not clip_id:
                continue
            for child in elem:
                match _element_local_bounds(child):
                    case Some(bounds):
                        clip_map[clip_id] = bounds
                        break
                    case _:
                        continue
    return clip_map


def _is_datapoint_mark_group(elem: ET.Element[str], eff_min_group: int) -> bool:
    """Return True if `elem` is a `<g>` container holding dense leaf datapoint marks."""
    if _strip_ns(elem.tag) != "g":
        return False
    children = list(elem)
    if len(children) < eff_min_group:
        return False

    cls = elem.attrib.get("class", "")
    if "role-axis" in cls or "role-legend" in cls or "role-title" in cls:
        return False

    return all(_strip_ns(ch.tag) in _DATAPOINT_CHILD_TAGS and len(list(ch)) == 0 for ch in children)


def _compute_group_bounding_box(
    elem: ET.Element[str],
    clip_rects: dict[str, tuple[float, float, float, float]],
    pad_px: float,
) -> Maybe[BoundingBox]:
    """Compute integer-snapped local BoundingBox for a datapoint mark group `<g>`."""
    min_xs: list[float] = []
    min_ys: list[float] = []
    max_xs: list[float] = []
    max_ys: list[float] = []

    for child in elem:
        match _element_local_bounds(child):
            case Some((x0, y0, x1, y1)):
                min_xs.append(x0)
                min_ys.append(y0)
                max_xs.append(x1)
                max_ys.append(y1)
            case _:
                continue

    if not min_xs or not min_ys:
        return Nothing

    bx0 = min(min_xs) - pad_px
    by0 = min(min_ys) - pad_px
    bx1 = max(max_xs) + pad_px
    by1 = max(max_ys) + pad_px

    clip_attr = elem.attrib.get("clip-path", "")
    m_clip = _CLIP_URL_RE.search(clip_attr)
    if m_clip and m_clip.group(1) in clip_rects:
        cx0, cy0, cx1, cy1 = clip_rects[m_clip.group(1)]
        bx0 = max(bx0, cx0 - pad_px)
        by0 = max(by0, cy0 - pad_px)
        bx1 = min(bx1, cx1 + pad_px)
        by1 = min(by1, cy1 + pad_px)

    ix0 = math.floor(bx0)
    iy0 = math.floor(by0)
    ix1 = math.ceil(max(bx1, bx0 + 1.0))
    iy1 = math.ceil(max(by1, by0 + 1.0))

    return Some(
        BoundingBox(
            x=ix0,
            y=iy0,
            width=max(1, ix1 - ix0),
            height=max(1, iy1 - iy0),
        )
    )


def _compress_png_idat(data: bytes) -> bytes:
    """Losslessly recompress PNG IDAT chunks using zlib level 9 for compact SVG embedding."""
    png_sig = b"\x89PNG\r\n\x1a\n"
    if not data.startswith(png_sig):
        return data
    pos = 8
    chunks: list[tuple[bytes, bytes]] = []
    idat_payloads: list[bytes] = []
    data_len = len(data)
    while pos + 12 <= data_len:
        length = struct.unpack(">I", data[pos : pos + 4])[0]
        ctype = data[pos + 4 : pos + 8]
        payload = data[pos + 8 : pos + 8 + length]
        pos += 12 + length
        if ctype == b"IDAT":
            idat_payloads.append(payload)
        else:
            chunks.append((ctype, payload))
    if not idat_payloads:
        return data
    try:
        raw_scanlines = zlib.decompress(b"".join(idat_payloads))
        compressed = zlib.compress(raw_scanlines, level=9)
    except Exception:
        return data

    out: list[bytes] = [png_sig]
    for ctype, payload in chunks:
        if ctype == b"IEND":
            idat_crc = zlib.crc32(b"IDAT" + compressed) & 0xFFFFFFFF
            out.append(
                struct.pack(">I", len(compressed))
                + b"IDAT"
                + compressed
                + struct.pack(">I", idat_crc)
            )
        chunk_crc = zlib.crc32(ctype + payload) & 0xFFFFFFFF
        out.append(
            struct.pack(">I", len(payload))
            + ctype
            + payload
            + struct.pack(">I", chunk_crc)
        )
    return b"".join(out)


def count_svg_datapoints(
    svg_content: str,
    min_group_points: int = 50,
) -> Result[int, str]:
    """Purely count the total number of datapoints across dense mark groups in an SVG string."""
    try:
        root = ET.fromstring(svg_content)
    except Exception as exc:
        return Failure(f"Failed to parse SVG XML in count_svg_datapoints: {exc}")

    eff_min = max(1, min_group_points)
    total = sum(
        len(list(elem))
        for elem in root.iter()
        if _is_datapoint_mark_group(elem, eff_min)
    )
    return Success(total)


def rasterize_svg_datapoints(
    svg_content: str,
    max_points: int = 1000,
    min_group_points: int = 50,
    scale: float = 3.0,
    pad_px: float = 6.0,
) -> Result[str, str]:
    """Rasterize datapoint mark groups in an SVG when total datapoints exceed `max_points`.

    Preserves all vector axes, gridlines, titles, text labels, legends, and low-count
    overlay layers while replacing dense datapoint `<g>` children with high-DPI
    transparent PNG `<image>` elements in local coordinates.
    """
    if max_points <= 0:
        return Success(svg_content)

    _register_svg_namespaces()
    try:
        root = ET.fromstring(svg_content)
    except Exception as exc:
        return Failure(f"Invalid XML syntax in rasterize_svg_datapoints: {exc}")

    eff_min_group = max(1, min(min_group_points, max_points))
    candidate_groups = [
        elem
        for elem in root.iter()
        if _is_datapoint_mark_group(elem, eff_min_group)
    ]
    total_points = sum(len(list(elem)) for elem in candidate_groups)
    if total_points <= max_points:
        return Success(svg_content)

    clip_rects = _extract_clip_rects(root)
    defs_elements = [elem for elem in root if _strip_ns(elem.tag) == "defs"]
    defs_xml = "".join(ET.tostring(d, encoding="unicode") for d in defs_elements)

    for group_elem in candidate_groups:
        match _compute_group_bounding_box(group_elem, clip_rects, pad_px):
            case Some(bbox):
                style_attrs = " ".join(
                    f'{k}="{html.escape(group_elem.attrib[k], quote=True)}"'
                    for k in _INHERITED_STYLE_ATTRS
                    if k in group_elem.attrib
                )
                default_fill = "" if "fill" in group_elem.attrib else 'fill="none"'
                default_miter = (
                    "" if "stroke-miterlimit" in group_elem.attrib else 'stroke-miterlimit="10"'
                )
                group_open_attrs = " ".join(
                    part for part in (default_fill, default_miter, style_attrs) if part
                )
                children_xml = "".join(
                    ET.tostring(ch, encoding="unicode") for ch in list(group_elem)
                )
                isolated_svg = (
                    f'<svg xmlns="{SVG_NS}" xmlns:xlink="{XLINK_NS}" '
                    f'width="{bbox.width}" height="{bbox.height}" '
                    f'viewBox="{bbox.x} {bbox.y} {bbox.width} {bbox.height}" version="1.1">'
                    f"{defs_xml}"
                    f"<g {group_open_attrs}>{children_xml}</g>"
                    f"</svg>"
                )
                try:
                    png_bytes = _compress_png_idat(vlc.svg_to_png(isolated_svg, scale=scale))
                except Exception as exc:
                    return Failure(f"Failed to rasterize SVG datapoint group via vl_convert: {exc}")

                b64_png = base64.b64encode(png_bytes).decode("ascii")
                preserved_attrib = dict(group_elem.attrib)
                group_elem.clear()
                group_elem.attrib.update(preserved_attrib)

                image_elem = ET.SubElement(group_elem, f"{{{SVG_NS}}}image")
                image_elem.set("x", str(bbox.x))
                image_elem.set("y", str(bbox.y))
                image_elem.set("width", str(bbox.width))
                image_elem.set("height", str(bbox.height))
                image_elem.set(f"{{{XLINK_NS}}}href", f"data:image/png;base64,{b64_png}")
                image_elem.set("preserveAspectRatio", "none")
                image_elem.set("style", "image-rendering: optimizeQuality;")
            case _:
                continue

    serialized_svg = ET.tostring(root, encoding="unicode")
    has_xml_decl = svg_content.lstrip().startswith("<?xml")
    xml_header = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' if has_xml_decl else ""
    return Success(f"{xml_header}{serialized_svg}")


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
    max_vector_points: int = 1000,
    point_raster_scale: float = 3.0,
) -> Result[Path, str]:
    """Validate XML well-formedness, rasterize dense datapoint layers, save SVG, and render 300 DPI PNG.

    Pure monadic function: checks layer counts, validates XML tree, automatically rasterizes
    datapoint mark groups when total datapoints exceed `max_vector_points`, and ensures
    reliable cross-platform rasterization.
    """
    try:
        if isinstance(svg_content_or_path, Path):
            output_svg_path = svg_content_or_path
            raw_svg_content = output_svg_path.read_text(encoding="utf-8")
        else:
            raw_svg_content = svg_content_or_path

        match rasterize_svg_datapoints(
            raw_svg_content,
            max_points=max_vector_points,
            scale=point_raster_scale,
        ):
            case Failure(err):
                return Failure(err)
            case Success(processed_svg):
                svg_content = processed_svg

        output_svg_path.parent.mkdir(parents=True, exist_ok=True)
        output_svg_path.write_text(svg_content, encoding="utf-8")
        root: ET.Element[str] = ET.fromstring(svg_content)

        # Count Inkscape layers
        layers = [
            elem
            for elem in root.iter()
            if elem.attrib.get(f"{{{INKSCAPE_NS}}}groupmode") == "layer"
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
                width = int(_parse_float(root.attrib.get("width"), default=1400.0))
            if height is None:
                height = int(_parse_float(root.attrib.get("height"), default=920.0))

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
                [
                    str(found_bin),
                    "-w",
                    str(raster_w),
                    "-h",
                    str(raster_h),
                    str(output_svg_path),
                    "-o",
                    str(output_png_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as cpe:
            return Failure(f"rsvg-convert execution failed ({found_bin}): {cpe.stderr}")
        except FileNotFoundError:
            # Fallback to in-memory vl_convert if rsvg-convert is missing from system
            png_bytes = vlc.svg_to_png(svg_content, scale=scale)
            output_png_path.write_bytes(png_bytes)

    return Success(output_svg_path)
