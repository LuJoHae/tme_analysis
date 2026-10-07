#!/usr/bin/env python3
"""
Combine All 12 Stability Selection Paths Figures into a Publication-Grade Composite.

Creates a unified 2-column x 6-row Nature Methods wireframe figure displaying
stability selection probability paths across all 12 curated immunotherapy cohorts:
  - 3 pooled cohorts (Pan-Cancer n=1,015, Melanoma n=338, RCC n=263)
  - 9 individual clinical trial cohorts (Rosenberg, McDermott, Liu, Riaz, Gide,
    Padron, Anders, Hugo, Choueiri)

Conforms strictly to:
  - Nature Methods Minimalist Wireframe standard via plotting_utils tokens.
  - Native Inkscape layer structure (namespaces, layer groups, semantic identifiers).
  - High-contrast publication typography and Okabe-Ito palette.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Final, Sequence
import xml.etree.ElementTree as ET
from returns.result import Failure, Result, Success

from plotting_utils import (
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_CARD_BG,
    COLOR_DIVIDER_RULE,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_HAIRLINE,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
    FONT_SERIF_MATH,
    OKABE_BLUISH_GREEN,
)
from plotting_utils.svg import verify_and_rasterize_svg, wrap_inkscape_svg

CANVAS_WIDTH: Final[int] = 1950
CANVAS_HEIGHT: Final[int] = 2880

CARD_WIDTH: Final[float] = 930.0
CARD_HEIGHT: Final[float] = 425.0
MARGIN_LEFT: Final[float] = 30.0
MARGIN_TOP: Final[float] = 135.0
COL_GUTTER: Final[float] = 30.0
ROW_GUTTER: Final[float] = 25.0


@dataclass(frozen=True)
class PanelConfig:
    """Metadata configuration for an individual cohort panel in the composite figure."""

    panel_id: str
    cohort_id: str
    display_name: str
    n_samples: int
    n_responders: int
    n_non_responders: int
    disease_context: str
    svg_filename: str
    selected_genes: tuple[str, ...]
    grid_row: int
    grid_col: int


PANEL_CONFIGS: Final[tuple[PanelConfig, ...]] = (
    # Row 0: Pooled Master Cohorts
    PanelConfig(
        panel_id="a",
        cohort_id="pancancer",
        display_name="Pan-Cancer Combined",
        n_samples=1015,
        n_responders=319,
        n_non_responders=696,
        disease_context="Curated Immunotherapy Panel • 12 Clinical Trials",
        svg_filename="paths_comparison_pancancer.svg",
        selected_genes=("HLA-C", "LAG3", "TBX21"),
        grid_row=0,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="b",
        cohort_id="melanoma",
        display_name="Melanoma Combined",
        n_samples=338,
        n_responders=132,
        n_non_responders=206,
        disease_context="Cutaneous Melanoma Pooled • Gide, Hugo, Liu, Riaz",
        svg_filename="paths_comparison_melanoma.svg",
        selected_genes=("HLA-A", "IKZF2", "TNFSF18", "TNFSF9"),
        grid_row=0,
        grid_col=1,
    ),
    # Row 1: Pooled RCC & Major Trial
    PanelConfig(
        panel_id="c",
        cohort_id="rcc",
        display_name="RCC Combined",
        n_samples=263,
        n_responders=75,
        n_non_responders=188,
        disease_context="Renal Cell Carcinoma Pooled • McDermott, Choueiri",
        svg_filename="paths_comparison_rcc.svg",
        selected_genes=("PTGS2",),
        grid_row=1,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="d",
        cohort_id="Rosenberg-iAtlas",
        display_name="Rosenberg et al.",
        n_samples=298,
        n_responders=68,
        n_non_responders=230,
        disease_context="Bladder Urothelial Carcinoma • Atezolizumab (anti-PD-L1)",
        svg_filename="paths_comparison_rosenberg_iatlas.svg",
        selected_genes=("IFNG", "TGFB1"),
        grid_row=1,
        grid_col=1,
    ),
    # Row 2: Renal & Melanoma Trials
    PanelConfig(
        panel_id="e",
        cohort_id="McDermott-iAtlas",
        display_name="McDermott et al.",
        n_samples=247,
        n_responders=72,
        n_non_responders=175,
        disease_context="Renal Cell Carcinoma • IMmotion150 (Atezolizumab +/- Bevacizumab)",
        svg_filename="paths_comparison_mcdermott_iatlas.svg",
        selected_genes=("PTGS2",),
        grid_row=2,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="f",
        cohort_id="Liu-iAtlas",
        display_name="Liu et al.",
        n_samples=122,
        n_responders=48,
        n_non_responders=74,
        disease_context="Cutaneous Melanoma • Nivolumab / Pembrolizumab (anti-PD-1)",
        svg_filename="paths_comparison_liu_iatlas.svg",
        selected_genes=(),
        grid_row=2,
        grid_col=1,
    ),
    # Row 3: Melanoma Trials
    PanelConfig(
        panel_id="g",
        cohort_id="Riaz-iAtlas",
        display_name="Riaz et al.",
        n_samples=98,
        n_responders=20,
        n_non_responders=78,
        disease_context="Cutaneous Melanoma • Nivolumab (anti-PD-1)",
        svg_filename="paths_comparison_riaz_iatlas.svg",
        selected_genes=("TNFSF9",),
        grid_row=3,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="h",
        cohort_id="Gide-iAtlas",
        display_name="Gide et al.",
        n_samples=91,
        n_responders=49,
        n_non_responders=42,
        disease_context="Cutaneous Melanoma • anti-PD-1 monotherapy / + Ipilimumab",
        svg_filename="paths_comparison_gide_iatlas.svg",
        selected_genes=("HLA-A",),
        grid_row=3,
        grid_col=1,
    ),
    # Row 4: Pancreatic & Smaller Bladder Trial
    PanelConfig(
        panel_id="i",
        cohort_id="Padron-iAtlas",
        display_name="Padron et al.",
        n_samples=85,
        n_responders=38,
        n_non_responders=47,
        disease_context="Pancreatic Ductal Adenocarcinoma (PDAC) • Nivolumab + Gem/Nab-Pac",
        svg_filename="paths_comparison_padron_iatlas.svg",
        selected_genes=(),
        grid_row=4,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="j",
        cohort_id="Anders-iAtlas",
        display_name="Anders et al.",
        n_samples=31,
        n_responders=7,
        n_non_responders=24,
        disease_context="Bladder Urothelial Carcinoma • anti-PD-1",
        svg_filename="paths_comparison_anders_iatlas.svg",
        selected_genes=(),
        grid_row=4,
        grid_col=1,
    ),
    # Row 5: Targeted Melanoma & RCC Cohorts
    PanelConfig(
        panel_id="k",
        cohort_id="Hugo-iAtlas",
        display_name="Hugo et al.",
        n_samples=27,
        n_responders=14,
        n_non_responders=13,
        disease_context="Cutaneous Melanoma • Pembrolizumab (anti-PD-1)",
        svg_filename="paths_comparison_hugo_iatlas.svg",
        selected_genes=(),
        grid_row=5,
        grid_col=0,
    ),
    PanelConfig(
        panel_id="l",
        cohort_id="Choueiri-iAtlas",
        display_name="Choueiri et al.",
        n_samples=16,
        n_responders=3,
        n_non_responders=13,
        disease_context="Renal Cell Carcinoma • Nivolumab (anti-PD-1)",
        svg_filename="paths_comparison_choueiri_iatlas.svg",
        selected_genes=(),
        grid_row=5,
        grid_col=1,
    ),
)


def clean_nested_svg(svg_path: Path) -> Result[str, str]:
    """Load an individual Vega-Lite SVG file and serialize back to clean SVG XML."""
    if not svg_path.exists():
        return Failure(f"Source SVG does not exist: {svg_path}")

    try:
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

        tree = ET.parse(svg_path)
        root = tree.getroot()

        # Remove the white rect at the background of the vega chart so card background shows cleanly
        for rect in list(root.findall("{http://www.w3.org/2000/svg}rect")):
            if rect.attrib.get("fill") == "#FFFFFF":
                root.remove(rect)

        # Build parent map to remove elements safely
        parent_map = {c: p for p in root.iter() for c in p}

        # Remove the redundant top-level role-title element
        for g in list(root.iter("{http://www.w3.org/2000/svg}g")):
            if g.attrib.get("class") == "mark-group role-title":
                texts = [
                    t.text
                    for t in g.iter("{http://www.w3.org/2000/svg}text")
                    if t.text and "Stability Selection Paths" in t.text
                ]
                if texts and g in parent_map:
                    parent_map[g].remove(g)

        return Success(ET.tostring(root, encoding="unicode"))
    except Exception as exc:
        return Failure(f"Failed to process SVG {svg_path.name}: {exc}")


def build_panel_svg(
    config: PanelConfig,
    cleaned_svg_content: str,
    card_x: float,
    card_y: float,
    card_w: float,
    card_h: float,
) -> tuple[str, str, str]:
    """Build an Inkscape-compliant layer tuple (layer_id, layer_label, inner_xml) for a single cohort panel."""
    layer_id = f"layer-{config.panel_id}-panel"
    layer_label = f"Panel_{config.panel_id.upper()}_{config.cohort_id}"

    # Summary pill
    if config.selected_genes:
        genes_str = ", ".join(config.selected_genes)
        pill_fill = "#ECFDF5"
        pill_stroke = "#A7F3D0"
        pill_text_color = OKABE_BLUISH_GREEN
        pill_text = f"{len(config.selected_genes)} Stable: {genes_str}"
    else:
        pill_fill = "#F8FAFC"
        pill_stroke = "#E2E8F0"
        pill_text_color = COLOR_TEXT_MUTED
        pill_text = "0 features reach π_thr = 0.75"

    pill_w = max(160.0, len(pill_text) * 5.8 + 18.0)
    pill_x = card_w - pill_w - 18.0

    inner_xml = f"""  <g id="grp-panel-{config.panel_id}" transform="translate({card_x}, {card_y})">
    <!-- Panel Background Card -->
    <rect id="card-bg-{config.panel_id}" x="0" y="0" width="{card_w}" height="{card_h}"
          fill="{COLOR_CARD_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />

    <!-- Card Header Strip -->
    <rect id="card-hdr-{config.panel_id}" x="0" y="0" width="{card_w}" height="42"
          fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" rx="0" />

    <!-- Panel Letter Badge (Minimalist Nature Style) -->
    <circle cx="24" cy="21" r="11" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="24" y="25.5" font-family="{FONT_SANS}" font-size="12" font-weight="700"
          fill="#FFFFFF" text-anchor="middle">{config.panel_id}</text>

    <!-- Cohort Title & Stratification Info -->
    <text x="44" y="19" font-family="{FONT_SANS}" font-size="13" font-weight="700"
          fill="{COLOR_TEXT_PRIMARY}">{config.display_name}</text>
    <text x="44" y="32" font-family="{FONT_SANS}" font-size="9.5" font-weight="500"
          fill="{COLOR_TEXT_MUTED}">n = {config.n_samples:,} patients ({config.n_responders} R / {config.n_non_responders} NR) • {config.disease_context}</text>

    <!-- Selection Status Pill Badge -->
    <rect x="{pill_x}" y="10" width="{pill_w}" height="22" rx="0"
          fill="{pill_fill}" stroke="{pill_stroke}" stroke-width="0.75" />
    <text x="{pill_x + pill_w / 2.0}" y="24.5" font-family="{FONT_SANS}" font-size="9.5" font-weight="600"
          fill="{pill_text_color}" text-anchor="middle">{pill_text}</text>

    <!-- Embedded Scaled Vega-Lite Content -->
    <g id="vega-content-{config.panel_id}" transform="translate(18, 52) scale(0.965)">
{cleaned_svg_content}
    </g>
  </g>"""

    return layer_id, layer_label, inner_xml


def render_figure_header() -> tuple[str, str, str]:
    """Generate Layer 01: Figure title, subtitle, and controlled error bounding banner."""
    inner_xml = f"""  <!-- Main Title -->
  <text x="30" y="44" font-family="{FONT_SANS}" font-size="22" font-weight="700"
        fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.4">
    Comprehensive Finite-Sample Stability Selection Paths Across 12 Clinical Immunotherapy Cohorts
  </text>

  <!-- Subtitle -->
  <text x="30" y="68" font-family="{FONT_SANS}" font-size="12.5" font-weight="400"
        fill="{COLOR_TEXT_MUTED}">
    Pan-cancer benchmarking of 1,015 pre-treatment patients comparing Meinshausen &amp; Bühlmann (2010) and Shah &amp; Samworth (2013) Complementary Pairs Selection
  </text>

  <!-- Statistical Control Chip -->
  <g id="grp-header-chip" transform="translate(1250, 20)">
    <rect x="0" y="0" width="670" height="54" fill="{COLOR_SUBTLE_FILL}"
          stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
    <text x="335" y="22" font-family="{FONT_SANS}" font-size="9.5" font-weight="700"
          fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" letter-spacing="0.4">
      FINITE-SAMPLE ERROR BOUNDS &amp; REGULARIZATION BUDGETING
    </text>
    <text x="335" y="40" font-family="{FONT_SERIF_MATH}" font-size="11" font-style="italic"
          fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
      Cutoff π<tspan font-size="8.5" baseline-shift="sub">thr</tspan> = 0.75   |   Budget q<tspan font-size="8.5" baseline-shift="sub">target</tspan> ≤ 20.0   |   PFER ≤ 0.97   |   2B = 100 Complementary Subsamples
    </text>
  </g>

  <!-- Top Divider Rule -->
  <line x1="30" y1="92" x2="{CANVAS_WIDTH - 30}" y2="92"
        stroke="{COLOR_DIVIDER_RULE}" stroke-width="1.0" />"""
    return "layer-01-header", "01_Figure_Header", inner_xml


def render_figure_footer() -> tuple[str, str, str]:
    """Generate Layer 03: Bottom metadata rule, methodology citation, and explanatory footnotes."""
    fy = CANVAS_HEIGHT - 32
    inner_xml = f"""  <line x1="30" y1="{fy}" x2="{CANVAS_WIDTH - 30}" y2="{fy}"
        stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
  <text x="30" y="{fy + 18}" font-family="{FONT_SANS}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_HAIRLINE}">
    Panels a–l represent independent stability paths fitted across 100 complementary subsample pairs (size ⌊n/2⌋). Solid colored paths denote stable biomarkers exceeding π_thr = 0.75.
  </text>
  <text x="{CANVAS_WIDTH - 30}" y="{fy + 18}" font-family="{FONT_SANS}" font-size="9.5" font-weight="500"
        fill="{COLOR_TEXT_MUTED}" text-anchor="end">
    Nature Methods Minimalist Wireframe Standard • Scaled at 300 DPI Vector Precision
  </text>"""
    return "layer-99-footer", "99_Figure_Footer", inner_xml


def assemble_composite_svg(configs: Sequence[PanelConfig], input_dir: Path) -> Result[str, str]:
    """Assemble all 12 panels into a publication-grade composite SVG document."""
    layers: list[tuple[str, str, str]] = [render_figure_header()]

    for cfg in configs:
        panel_x = MARGIN_LEFT + cfg.grid_col * (CARD_WIDTH + COL_GUTTER)
        panel_y = MARGIN_TOP + cfg.grid_row * (CARD_HEIGHT + ROW_GUTTER)

        svg_path = input_dir / cfg.svg_filename
        match clean_nested_svg(svg_path):
            case Failure(err):
                return Failure(f"Error on panel {cfg.panel_id} ({cfg.cohort_id}): {err}")
            case Success(cleaned_svg):
                panel_layer = build_panel_svg(
                    cfg, cleaned_svg, panel_x, panel_y, CARD_WIDTH, CARD_HEIGHT
                )
                layers.append(panel_layer)

    layers.append(render_figure_footer())

    full_svg = wrap_inkscape_svg(
        layers=layers,
        width=CANVAS_WIDTH,
        height=CANVAS_HEIGHT,
        title="Composite Stability Paths (12 Cohorts)",
        bg_color=COLOR_CANVAS_BG,
    )
    return Success(full_svg)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble all 12 cohort stability paths into a single publication-grade composite SVG."
    )
    default_fig_dir = (
        Path("output/stability_selection_fitters_benchmark/figures")
        if (Path("output/stability_selection_fitters_benchmark/figures") / "paths_comparison_pancancer.svg").exists()
        else Path("output/stability_selection_iatlas_immunotherapy/figures")
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_fig_dir,
        help="Directory containing the individual paths_comparison_*.svg files",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=default_fig_dir / "all_cohorts_stability_paths_composite.svg",
        help="Destination path for composite SVG",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=default_fig_dir / "all_cohorts_stability_paths_composite.png",
        help="Destination path for 300 DPI rasterized PNG",
    )
    parser.add_argument(
        "--rsvg-path",
        type=str,
        default="/opt/homebrew/bin/rsvg-convert",
        help="Path to rsvg-convert executable",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    print("==================================================================")
    print(" Stability Selection Composite Figure Assembler (Nature Methods)")
    print("==================================================================")
    print(f" Input Directory:  {args.input_dir}")
    print(f" Output SVG Path:  {args.output_svg}")
    print(f" Output PNG Path:  {args.output_png}")
    print(f" Panels to Merge:  {len(PANEL_CONFIGS)} cohorts (2 cols x 6 rows)")

    match assemble_composite_svg(PANEL_CONFIGS, args.input_dir):
        case Failure(err):
            print(f"[-] Assembly failed: {err}")
            return 1
        case Success(svg_content):
            res = verify_and_rasterize_svg(
                svg_content_or_path=svg_content,
                output_svg_path=args.output_svg,
                output_png_path=args.output_png,
                width=CANVAS_WIDTH,
                height=CANVAS_HEIGHT,
                scale=2.0,
                min_layers=len(PANEL_CONFIGS) + 2,
                rsvg_bin=args.rsvg_path,
            )
            match res:
                case Failure(err):
                    print(f"[-] Verification/Rasterization failed: {err}")
                    return 1
                case Success(svg_p):
                    print(f"[+] Successfully wrote composite SVG to: {svg_p}")
                    if args.output_png.exists():
                        print(f"[+] Successfully rasterized 300 DPI PNG to: {args.output_png}")
                    print("[+] Composite figure pipeline finished successfully.")
                    return 0


if __name__ == "__main__":
    sys.exit(main())
