#!/usr/bin/env python3
"""
Publication-Grade 2x2 Master Figure Assembler across 180mm Total Journal Width.

Combines four analyses into a unified 2x2 publication vector figure:
- Row 1:
  - Panel a (Top-Left, 130mm): 3x3 Single-Cell Landscape UMAP Grid (continuous Milo log2FC)
  - Panel b (Top-Right, 50mm): Cross-Cohort Differential Abundance Diverging Bar Chart
- Row 2:
  - Panel c (Bottom-Left, 130mm): 3x3 Empirical Deconvolution Concordance Scatter Grid
  - Panel d (Bottom-Right, 50mm): Concordance Summary Correlation Bar Plot

Adheres strictly to Nature Methods minimalist wireframe standard:
- Pure white canvas (#FFFFFF), zero drop shadows, hairline borders (#CBD5E1)
- Inkscape-compliant native layer architecture (Background, Row1, Row2)
- Minimum font size >= 10pt (13.5px at 96 CSS DPI)
- 300 DPI rasterization via rsvg-convert
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
from typing import Final
import xml.etree.ElementTree as ET

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# ------------------------------------------------------------------------------
# Dimensional & Geometry Constants
# ------------------------------------------------------------------------------
TOTAL_WIDTH_MM: Final[float] = 180.0
LEFT_COL_WIDTH_MM: Final[float] = 130.0
RIGHT_COL_WIDTH_MM: Final[float] = 50.0

ROW_HEIGHT_MM: Final[float] = 123.31
GUTTER_HEIGHT_MM: Final[float] = 6.4
TOTAL_HEIGHT_MM: Final[float] = ROW_HEIGHT_MM * 2.0 + GUTTER_HEIGHT_MM  # ~253.02 mm

# ViewBox Coordinates (Matching 486 px left + 188 px right, 461 px rows + 24 px gutter)
VIEWBOX_LEFT_WIDTH: Final[int] = 486
VIEWBOX_RIGHT_WIDTH: Final[int] = 188
VIEWBOX_TOTAL_WIDTH: Final[int] = VIEWBOX_LEFT_WIDTH + VIEWBOX_RIGHT_WIDTH  # 674 px

VIEWBOX_ROW_HEIGHT: Final[int] = 461
VIEWBOX_GUTTER_HEIGHT: Final[int] = 24
VIEWBOX_TOTAL_HEIGHT: Final[int] = VIEWBOX_ROW_HEIGHT * 2 + VIEWBOX_GUTTER_HEIGHT  # 946 px


class AssemblerConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    panel_a_svg: Path = Path(
        "output/manuscript_pipeline/4a3_plot_cross_cohort_umaps_grid/nine_cohorts_umaps_logfc_grid.svg"
    )
    panel_b_svg: Path = Path(
        "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.svg"
    )
    panel_c_svg: Path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
    )
    panel_d_svg: Path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.svg"
    )

    row1_composite_svg: Path = Path(
        "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/nine_cohorts_single_cell_composite.svg"
    )
    row1_composite_png: Path = Path(
        "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/nine_cohorts_single_cell_composite.png"
    )

    master_svg: Path = Path(
        "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/immunotherapy_response_2x2_master.svg"
    )
    master_png: Path = Path(
        "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/immunotherapy_response_2x2_master.png"
    )
    article_master_svg: Path = Path("article/figures/milo-deconv/immunotherapy_response_2x2_master.svg")
    article_master_png: Path = Path("article/figures/milo-deconv/immunotherapy_response_2x2_master.png")


def extract_inner_svg(svg_str: str) -> Result[str, str]:
    """Extracts elements inside the outermost <svg> tag while stripping root tags."""
    clean = svg_str.strip()
    # Strip XML declaration if present
    clean = re.sub(r"<\?xml[^>]*\?>", "", clean).strip()
    # Strip doctype if present
    clean = re.sub(r"<!DOCTYPE[^>]*>", "", clean).strip()

    m = re.search(r"^<svg[^>]*>(.*)</svg>$", clean, re.DOTALL)
    if not m:
        return Failure("Could not parse root <svg> wrapper")
    return Success(m.group(1).strip())


def stitch_row1_composite(
    panel_a_str: str,
    panel_b_str: str,
) -> Result[str, str]:
    """Stitches Panel a (130mm) and Panel b (50mm) into a 180mm x 154.07mm Row 1 composite."""
    inner_a = extract_inner_svg(panel_a_str)
    inner_b = extract_inner_svg(panel_b_str)

    if isinstance(inner_a, Failure):
        return inner_a
    if isinstance(inner_b, Failure):
        return inner_b

    svg_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     version="1.1"
     width="{TOTAL_WIDTH_MM:.1f}mm"
     height="{ROW_HEIGHT_MM:.2f}mm"
     viewBox="0 0 {VIEWBOX_TOTAL_WIDTH} {VIEWBOX_ROW_HEIGHT}">
  <defs/>
  <rect id="bg-rect" width="{VIEWBOX_TOTAL_WIDTH}" height="{VIEWBOX_ROW_HEIGHT}" fill="#FFFFFF"/>
  <g inkscape:groupmode="layer" id="layer-panel-a" inkscape:label="Panel_A_Landscape" transform="translate(0, 0)">
    {inner_a.unwrap()}
  </g>
  <g inkscape:groupmode="layer" id="layer-panel-b" inkscape:label="Panel_B_DifferentialAbundance" transform="translate({VIEWBOX_LEFT_WIDTH}, 0)">
    {inner_b.unwrap()}
  </g>
</svg>"""
    return Success(svg_content)


def stitch_master_2x2_figure(
    panel_a_str: str,
    panel_b_str: str,
    panel_c_str: str,
    panel_d_str: str,
) -> Result[str, str]:
    """Stitches all four panels into the master 2x2 publication vector figure."""
    inner_a = extract_inner_svg(panel_a_str)
    inner_b = extract_inner_svg(panel_b_str)
    inner_c = extract_inner_svg(panel_c_str)
    inner_d = extract_inner_svg(panel_d_str)

    for name, res in [("a", inner_a), ("b", inner_b), ("c", inner_c), ("d", inner_d)]:
        if isinstance(res, Failure):
            return Failure(f"Failed to extract inner SVG for Panel {name}: {res}")

    y_row2 = VIEWBOX_ROW_HEIGHT + VIEWBOX_GUTTER_HEIGHT

    svg_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     version="1.1"
     width="{TOTAL_WIDTH_MM:.1f}mm"
     height="{TOTAL_HEIGHT_MM:.2f}mm"
     viewBox="0 0 {VIEWBOX_TOTAL_WIDTH} {VIEWBOX_TOTAL_HEIGHT}">
  <defs/>
  <!-- Layer 0: Background Canvas -->
  <g inkscape:groupmode="layer" id="layer-background" inkscape:label="Background">
    <rect id="master-bg-canvas" width="{VIEWBOX_TOTAL_WIDTH}" height="{VIEWBOX_TOTAL_HEIGHT}" fill="#FFFFFF"/>
  </g>

  <!-- Layer 1: Row 1 Single-Cell Differential Abundance Testing (Panels a & b) -->
  <g inkscape:groupmode="layer" id="layer-row1-single-cell" inkscape:label="Row1_SingleCell_DA" transform="translate(0, 0)">
    <!-- Panel a: 3x3 UMAP Landscape Grid (130 mm) -->
    <g id="panel-a-umap-grid" transform="translate(0, 0)">
      {inner_a.unwrap()}
    </g>
    <!-- Panel b: Differential Abundance Diverging Bars (50 mm) -->
    <g id="panel-b-da-bars" transform="translate({VIEWBOX_LEFT_WIDTH}, 0)">
      {inner_b.unwrap()}
    </g>
  </g>

  <!-- Layer 2: Row 2 Deconvolution Concordance Benchmark (Panels c & d) -->
  <g inkscape:groupmode="layer" id="layer-row2-deconvolution" inkscape:label="Row2_Deconvolution_Benchmark" transform="translate(0, {y_row2})">
    <!-- Panel c: 3x3 Deconvolution Scatter Grid (130 mm) -->
    <g id="panel-c-scatter-grid" transform="translate(0, 0)">
      {inner_c.unwrap()}
    </g>
    <!-- Panel d: Concordance Summary Bar Plot (50 mm) -->
    <g id="panel-d-correlation-bar" transform="translate({VIEWBOX_LEFT_WIDTH}, 0)">
      {inner_d.unwrap()}
    </g>
  </g>
</svg>"""
    return Success(svg_content)


def rasterize_svg_with_rsvg(svg_path: Path, png_path: Path, width_px: int = 2126) -> Result[Path, str]:
    """Rasterizes SVG to high-resolution 300 DPI PNG via rsvg-convert."""
    cmd = [
        "rsvg-convert",
        "-w", str(width_px),
        "-o", str(png_path),
        str(svg_path),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if not png_path.exists() or png_path.stat().st_size == 0:
            return Failure(f"rsvg-convert produced empty PNG at {png_path}")
        return Success(png_path)
    except subprocess.CalledProcessError as exc:
        return Failure(f"rsvg-convert failed: {exc.stderr}")
    except FileNotFoundError:
        return Failure("rsvg-convert executable not found in PATH")


def execute_assembly(config: AssemblerConfig) -> Result[Path, str]:
    """Orchestrates ingestion, assembly, validation, and rasterization."""
    print("=== Assembling Master 2x2 Immunotherapy Publication Figure (180 mm) ===")

    # Verify input SVGs exist
    for p_name, p_path in [
        ("Panel a (UMAP grid)", config.panel_a_svg),
        ("Panel b (DA bars)", config.panel_b_svg),
        ("Panel c (Scatter grid)", config.panel_c_svg),
        ("Panel d (Correlation bar)", config.panel_d_svg),
    ]:
        if not p_path.exists():
            return Failure(f"Missing input SVG for {p_name} at {p_path}")
        print(f"  Found {p_name:28s}: {p_path}")

    a_str = config.panel_a_svg.read_text(encoding="utf-8")
    b_str = config.panel_b_svg.read_text(encoding="utf-8")
    c_str = config.panel_c_svg.read_text(encoding="utf-8")
    d_str = config.panel_d_svg.read_text(encoding="utf-8")

    # 1. Row 1 composite
    print(f"\n1. Generating Row 1 Composite SVG (180 mm x {ROW_HEIGHT_MM:.2f} mm)...")
    config.row1_composite_svg.parent.mkdir(parents=True, exist_ok=True)
    match stitch_row1_composite(a_str, b_str):
        case Failure(err):
            return Failure(f"Failed to stitch Row 1 composite: {err}")
        case Success(row1_svg):
            config.row1_composite_svg.write_text(row1_svg, encoding="utf-8")
            rasterize_svg_with_rsvg(config.row1_composite_svg, config.row1_composite_png, width_px=2126)
            print(f"  Saved Row 1 SVG -> {config.row1_composite_svg}")
            print(f"  Saved Row 1 PNG -> {config.row1_composite_png}")

    # 2. Master 2x2 figure
    print(f"\n2. Generating Master 2x2 Figure SVG (180 mm x {TOTAL_HEIGHT_MM:.2f} mm)...")
    config.master_svg.parent.mkdir(parents=True, exist_ok=True)
    config.article_master_svg.parent.mkdir(parents=True, exist_ok=True)

    match stitch_master_2x2_figure(a_str, b_str, c_str, d_str):
        case Failure(err):
            return Failure(f"Failed to stitch master 2x2 figure: {err}")
        case Success(master_svg):
            config.master_svg.write_text(master_svg, encoding="utf-8")
            config.article_master_svg.write_text(master_svg, encoding="utf-8")

            # Rasterize master figure at 300 DPI
            print("  Rasterizing 300 DPI PNG via rsvg-convert (2126 px width)...")
            rasterize_svg_with_rsvg(config.master_svg, config.master_png, width_px=2126)
            config.article_master_png.write_bytes(config.master_png.read_bytes())

            print(f"  Saved Master SVG -> {config.master_svg}")
            print(f"  Saved Master PNG -> {config.master_png}")
            print(f"  Saved Article Master SVG -> {config.article_master_svg}")
            print(f"  Saved Article Master PNG -> {config.article_master_png}")

    # Validate XML well-formedness
    try:
        ET.parse(config.master_svg)
        print("\nXML validation PASSED: Master SVG is well-formed.")
    except Exception as exc:
        return Failure(f"XML validation failed for master SVG: {exc}")

    return Success(config.master_svg)


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble 2x2 master immunotherapy figure.")
    parser.add_argument(
        "--panel-a",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a3_plot_cross_cohort_umaps_grid/nine_cohorts_umaps_logfc_grid.svg"
        ),
    )
    parser.add_argument(
        "--panel-b",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.svg"
        ),
    )
    parser.add_argument(
        "--panel-c",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
        ),
    )
    parser.add_argument(
        "--panel-d",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.svg"
        ),
    )
    parser.add_argument(
        "--out-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/immunotherapy_response_2x2_master.svg"
        ),
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure/immunotherapy_response_2x2_master.png"
        ),
    )
    args = parser.parse_args()

    config = AssemblerConfig(
        panel_a_svg=args.panel_a,
        panel_b_svg=args.panel_b,
        panel_c_svg=args.panel_c,
        panel_d_svg=args.panel_d,
        master_svg=args.out_svg,
        master_png=args.out_png,
    )

    match execute_assembly(config):
        case Failure(err):
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"\nAssembly Complete: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
