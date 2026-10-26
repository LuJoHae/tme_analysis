#!/usr/bin/env python3
"""
Generate a conference-room-readable PowerPoint slide table of the 9 Bulk RNA-seq
Immune Checkpoint Blockade (ICB) Clinical Trial Cohorts (n = 1,097 patients).

Specifications:
  - Inner Table / Chart Area: 400 px width x 450 px height (10 rows x 45 px)
  - Main Chart Title: 24 pt
  - Axis / Column Titles: 20 pt
  - Legend Title: 18 pt
  - Axis / Row and Legend Labels (tick marks & cells): 16 pt
  - Nature Methods wireframe aesthetic with Okabe-Ito colorblind-safe palette.
"""

import html
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Final

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# ==============================================================================
# DIMENSIONS & TYPOGRAPHY CONSTANTS
# ==============================================================================
INNER_W: Final[int] = 400
INNER_H: Final[int] = 450
ROW_H: Final[int] = 45

GRID_X: Final[int] = 186
GRID_Y: Final[int] = 86

CANVAS_W: Final[int] = 600
CANVAS_H: Final[int] = 620

FONT_MAIN_TITLE: Final[int] = 24
FONT_AXIS_TITLE: Final[int] = 20
FONT_LEGEND_TITLE: Final[int] = 18
FONT_LABELS: Final[int] = 16

COLOR_CANVAS_BG: Final[str] = "#FFFFFF"
COLOR_PANEL_BG: Final[str] = "#FFFFFF"
COLOR_SUBTLE_FILL: Final[str] = "#F8FAFC"
COLOR_SUMMARY_FILL: Final[str] = "#F1F5F9"
COLOR_BORDER_HAIRLINE: Final[str] = "#CBD5E1"
COLOR_AXIS_LINE: Final[str] = "#1E293B"
COLOR_DIVIDER_RULE: Final[str] = "#E2E8F0"

COLOR_TEXT_PRIMARY: Final[str] = "#0F172A"
COLOR_TEXT_SECONDARY: Final[str] = "#334155"

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_SKY_BLUE: Final[str] = "#56B4E9"
OKABE_BLUISH_GREEN: Final[str] = "#009E73"
OKABE_ORANGE: Final[str] = "#E69F00"
OKABE_VERMILION: Final[str] = "#D55E00"
OKABE_REDDISH_PURPLE: Final[str] = "#CC79A7"

FONT_SANS: Final[str] = "system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"


# ==============================================================================
# IMMUTABLE DATA MODEL (FROZEN PYDANTIC, NO NONE)
# ==============================================================================
class SlideBulkCohort(BaseModel):
    """Immutable model for a bulk RNA-seq ICB clinical trial cohort on a PowerPoint slide."""

    model_config = ConfigDict(frozen=True)

    study: str
    cancer_type: str
    cancer_color: str
    n_total: int
    n_resp: int
    n_non_resp: int
    orr_pct: float


SLIDE_BULK_COHORTS: Final[tuple[SlideBulkCohort, ...]] = (
    SlideBulkCohort(
        study="Hugo (2016)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        n_total=27,
        n_resp=14,
        n_non_resp=13,
        orr_pct=51.9,
    ),
    SlideBulkCohort(
        study="Riaz (2017)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        n_total=107,
        n_resp=34,
        n_non_resp=73,
        orr_pct=31.8,
    ),
    SlideBulkCohort(
        study="Liu (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        n_total=122,
        n_resp=48,
        n_non_resp=74,
        orr_pct=39.3,
    ),
    SlideBulkCohort(
        study="Gide (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        n_total=91,
        n_resp=63,
        n_non_resp=28,
        orr_pct=69.2,
    ),
    SlideBulkCohort(
        study="Rosenberg (2016)",
        cancer_type="Bladder",
        cancer_color=OKABE_SKY_BLUE,
        n_total=347,
        n_resp=68,
        n_non_resp=279,
        orr_pct=19.6,
    ),
    SlideBulkCohort(
        study="McDermott (2018)",
        cancer_type="ccRCC",
        cancer_color=OKABE_BLUISH_GREEN,
        n_total=263,
        n_resp=97,
        n_non_resp=166,
        orr_pct=36.9,
    ),
    SlideBulkCohort(
        study="Choueiri (2016)",
        cancer_type="ccRCC",
        cancer_color=OKABE_BLUISH_GREEN,
        n_total=16,
        n_resp=6,
        n_non_resp=10,
        orr_pct=37.5,
    ),
    SlideBulkCohort(
        study="Padron (2022)",
        cancer_type="PDAC",
        cancer_color=OKABE_BLUE,
        n_total=93,
        n_resp=16,
        n_non_resp=77,
        orr_pct=17.2,
    ),
    SlideBulkCohort(
        study="Anders (2021)",
        cancer_type="TNBC",
        cancer_color=OKABE_REDDISH_PURPLE,
        n_total=31,
        n_resp=11,
        n_non_resp=20,
        orr_pct=35.5,
    ),
)


# ==============================================================================
# PURE FUNCTIONAL SVG RENDERERS
# ==============================================================================
def render_slide_row(idx: int, item: SlideBulkCohort) -> str:
    """Render a single cohort row (height = 45px, inner width = 400px, font-size = 16pt)."""
    y_offset = idx * ROW_H
    bar_w = 144
    raw_w_r = int(round((item.n_resp / item.n_total) * bar_w))
    w_r = max(48, min(bar_w - 58, raw_w_r))
    w_nr = bar_w - w_r
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL
    r_center = w_r // 2
    nr_center = w_r + (w_nr // 2)

    return f"""
      <!-- Cohort Row {idx + 1}: {html.escape(item.study)} -->
      <g id="slide-row-{idx + 1:02d}" transform="translate(0, {y_offset})">
        <!-- Inner 400px row background & divider -->
        <rect x="0" y="0" width="{INNER_W}" height="{ROW_H}" fill="{row_bg}" />
        <line x1="0" y1="{ROW_H}" x2="{INNER_W}" y2="{ROW_H}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />

        <!-- Outward Y-axis tick mark & Cohort Label (16pt) -->
        <line x1="-5" y1="22.5" x2="0" y2="22.5" stroke="{COLOR_AXIS_LINE}" stroke-width="1.0" />
        <text x="-10" y="28.0" font-size="{FONT_LABELS}" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">{html.escape(item.study)}</text>

        <!-- Cancer Type with Okabe-Ito indicator dot (16pt) -->
        <circle cx="13" cy="22.5" r="5.0" fill="{item.cancer_color}" />
        <text x="23" y="28.0" font-size="{FONT_LABELS}" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.cancer_type)}</text>

        <!-- Sample Size N (16pt) -->
        <text x="162" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{item.n_total}</text>

        <!-- Stacked RECIST Response Bar (R / NR, 16pt) -->
        <g transform="translate(172, 8.5)">
          <rect x="0" y="0" width="{w_r}" height="28" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{w_r}" y="0" width="{w_nr}" height="28" fill="{OKABE_VERMILION}" />
          <text x="{r_center}" y="19.5" font-size="{FONT_LABELS}" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.3">{item.n_resp}R</text>
          <text x="{nr_center}" y="19.5" font-size="{FONT_LABELS}" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.3">{item.n_non_resp}NR</text>
        </g>

        <!-- Objective Response Rate ORR % (16pt) -->
        <text x="392" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{item.orr_pct:.1f}%</text>
      </g>"""


def render_svg() -> str:
    """Declaratively render the complete PowerPoint slide table SVG."""
    total_n = sum(c.n_total for c in SLIDE_BULK_COHORTS)
    total_r = sum(c.n_resp for c in SLIDE_BULK_COHORTS)
    total_nr = sum(c.n_non_resp for c in SLIDE_BULK_COHORTS)
    total_orr = (total_r / total_n) * 100.0

    bar_w = 144
    total_w_r = max(54, int(round((total_r / total_n) * bar_w)))
    total_w_nr = bar_w - total_w_r

    rows_svg = "".join(
        render_slide_row(idx, item)
        for idx, item in enumerate(SLIDE_BULK_COHORTS)
    )

    summary_y = 9 * ROW_H

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{CANVAS_W}" height="{CANVAS_H}" viewBox="0 0 {CANVAS_W} {CANVAS_H}"
     version="1.1" id="figure-bulk-cohorts-ppt-slide"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <!-- ====================================================================== -->
  <!-- LAYER 00: CANVAS BACKGROUND & OUTER WIREFRAME                          -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{CANVAS_W}" height="{CANVAS_H}" fill="{COLOR_CANVAS_BG}" />
    <rect id="card-outer" x="4" y="4" width="{CANVAS_W - 8}" height="{CANVAS_H - 8}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 01: MAIN CHART TITLE (24 PT)                                     -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-01-main-title" inkscape:label="01_Main_Title">
    <text id="txt-main-title" x="{CANVAS_W // 2}" y="36"
          font-size="{FONT_MAIN_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}"
          text-anchor="middle" letter-spacing="-0.3">Bulk RNA-seq ICB Clinical Cohorts (n = 1,097)</text>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 02: AXIS & COLUMN TITLES (20 PT)                                 -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-02-axis-titles" inkscape:label="02_Axis_And_Column_Headers">
    <g id="grp-column-headers" transform="translate({GRID_X}, {GRID_Y})">
      <!-- Top header rule -->
      <rect x="0" y="-36" width="{INNER_W}" height="36" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="-10" y="-10" font-size="{FONT_AXIS_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">Trial Cohort</text>
      <text x="10" y="-10" font-size="{FONT_AXIS_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Cancer</text>
      <text x="162" y="-10" font-size="{FONT_AXIS_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">N</text>
      <text x="244" y="-10" font-size="{FONT_AXIS_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" letter-spacing="-0.2">RECIST</text>
      <text x="392" y="-10" font-size="{FONT_AXIS_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">ORR</text>
    </g>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 03: INNER TABLE / CHART AREA (400 x 450 PX, 16 PT LABELS)        -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-03-inner-table" inkscape:label="03_Inner_Table_400x450">
    <g id="grp-inner-table-400x450" transform="translate({GRID_X}, {GRID_Y})">
      {rows_svg}

      <!-- Summary Row 10: Overall (y = 405, height = 45px) -->
      <g id="slide-row-summary" transform="translate(0, {summary_y})">
        <rect x="0" y="0" width="{INNER_W}" height="{ROW_H}" fill="{COLOR_SUMMARY_FILL}" />
        <line x1="0" y1="0" x2="{INNER_W}" y2="0" stroke="{COLOR_AXIS_LINE}" stroke-width="1.0" />
        <line x1="-5" y1="22.5" x2="0" y2="22.5" stroke="{COLOR_AXIS_LINE}" stroke-width="1.2" />
        <text x="-10" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">Overall (9)</text>
        <text x="10" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">5 Types</text>
        <text x="162" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{total_n:,}</text>
        <g transform="translate(172, 8.5)">
          <rect x="0" y="0" width="{total_w_r}" height="28" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{total_w_r}" y="0" width="{total_w_nr}" height="28" fill="{OKABE_VERMILION}" />
          <text x="{total_w_r // 2}" y="19.5" font-size="{FONT_LABELS}" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.3">{total_r}R</text>
          <text x="{total_w_r + (total_w_nr // 2)}" y="19.5" font-size="{FONT_LABELS}" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.3">{total_nr}NR</text>
        </g>
        <text x="392" y="28.0" font-size="{FONT_LABELS}" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{total_orr:.1f}%</text>
      </g>

      <!-- Exact 400 x 450 px Inner Plot / Table Border -->
      <rect id="rect-inner-plot-400x450" x="0" y="0" width="{INNER_W}" height="{INNER_H}"
            fill="none" stroke="{COLOR_AXIS_LINE}" stroke-width="1.2" />
    </g>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 04: LEGEND BLOCK (18 PT TITLE, 16 PT LABELS)                     -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-04-legend" inkscape:label="04_Legend">
    <g id="grp-legend" transform="translate(16, 548)">
      <rect x="0" y="0" width="{CANVAS_W - 32}" height="60" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Legend Title (18pt) -->
      <text x="12" y="23" font-size="{FONT_LEGEND_TITLE}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">RECIST v1.1 Response:</text>

      <!-- Legend Labels (16pt) -->
      <rect x="12" y="35" width="15" height="15" fill="{OKABE_BLUISH_GREEN}" />
      <text x="33" y="48" font-size="{FONT_LABELS}" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">R: Responder (357, 32.5%)</text>

      <rect x="278" y="35" width="15" height="15" fill="{OKABE_VERMILION}" />
      <text x="299" y="48" font-size="{FONT_LABELS}" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">NR: Non-Responder (740, 67.5%)</text>
    </g>
  </g>
</svg>
"""


def validate_and_export(svg_content: str, out_paths: tuple[Path, ...]) -> Result[tuple[Path, ...], str]:
    """Validate XML structure, 400x450 inner dimensions, and 16-24pt font sizes, then export SVG & PNG."""
    try:
        root = ET.fromstring(svg_content)
    except ET.ParseError as exc:
        return Failure(f"Invalid SVG XML generated: {exc}")

    ns_inkscape = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
    layers = tuple(elem for elem in root.iter() if elem.attrib.get(ns_inkscape) == "layer")
    if len(layers) < 5:
        return Failure(f"Expected 5 Inkscape layers, found {len(layers)}")

    # Verify the inner table/chart area is exactly 400 x 450 px
    inner_rect = next(
        (elem for elem in root.iter() if elem.attrib.get("id") == "rect-inner-plot-400x450"),
        root,
    )
    if inner_rect is root or inner_rect.attrib.get("width") != "400" or inner_rect.attrib.get("height") != "450":
        return Failure("Inner table/chart area rect-inner-plot-400x450 does not have width=400 and height=450")

    # Verify every <text> element uses allowed conference-room font sizes {16, 18, 20, 24}
    allowed_sizes = {16.0, 18.0, 20.0, 24.0}
    for elem in root.iter("{http://www.w3.org/2000/svg}text"):
        fs_str = elem.attrib.get("font-size", "0")
        fs_val = float(fs_str.rstrip("ptx"))
        if fs_val not in allowed_sizes:
            return Failure(f"Found text element with font-size {fs_str} not in {sorted(allowed_sizes)}")

    written: list[Path] = []
    for svg_path in out_paths:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        svg_path.write_text(svg_content, encoding="utf-8")
        written.append(svg_path)

        png_path = svg_path.with_suffix(".png")
        try:
            subprocess.run(
                [
                    "rsvg-convert",
                    "-d",
                    "300",
                    "-p",
                    "300",
                    "-w",
                    "1800",
                    "-b",
                    "white",
                    str(svg_path),
                    "-o",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            written.append(png_path)
        except FileNotFoundError:
            pass
        except subprocess.CalledProcessError as exc:
            return Failure(f"rsvg-convert failed for {svg_path}: {exc.stderr}")

    return Success(tuple(written))


def main() -> int:
    """Imperative shell entry point."""
    repo_root = Path(__file__).resolve().parents[2]
    targets = (
        repo_root / "article" / "figures" / "dataset_overview" / "figure_bulk_cohorts_ppt_slide.svg",
        repo_root / "output" / "reports" / "figure_bulk_cohorts_ppt_slide.svg",
    )

    res = validate_and_export(render_svg(), targets)
    if isinstance(res, Failure):
        sys.stderr.write(f"[ERROR] {res.failure()}\n")
        return 1
    paths = res.unwrap()
    for p in paths:
        size_bytes = p.stat().st_size
        sys.stdout.write(f"[OK] Generated {p.suffix.upper().lstrip('.')}: {p.relative_to(repo_root)} ({size_bytes:,} bytes)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
