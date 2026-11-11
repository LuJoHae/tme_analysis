#!/usr/bin/env python3
"""Generate the 70mm x 150mm TME & ICB Benchmark Pipeline Graphical Abstract SVG.

Constructs a publication-grade vertical graphical abstract (width=70mm, height=150mm,
viewBox="0 0 462 990", minimum font-size >= 10.0) adhering to Nature Methods minimalist
wireframe standards:
  - Panel 1: Tumor Microenvironment Ecosystem (embeds assets/tme_tumor_microenvironment_ecosystem.svg)
             with 4 annotated cellular compartments.
  - Panel 2: Deconvolution Inference Pipeline from scRNA-seq reference, Bayesian Deconvolution
             (BayesPrism optical prism illustration), and Response classification (Logistic
             Regression & Random Forest).
  - Panel 3: 9 Pre-Treatment Clinical ICB Cohorts (5 cancer types with sample sizes, N = 1,015)
             and 27 published RNA-seq baseline predictors grouped into 3 categories.
  - Panel 4: Complete-matrix PR-AUC benchmark evaluation across classifiers and cohorts.

Also embeds the generated 70mm x 150mm graphical abstract into the top-left slot of
Page 1 in article/figures/inkscape/figure_01_overview.svg.
"""

from __future__ import annotations

import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Final, Sequence

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# Ensure scripts/figures is in sys.path for shared style tokens
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style_config import (  # noqa: E402
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_PANEL_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
    FONT_SERIF_MATH,
    OKABE_BLUISH_GREEN,
    OKABE_BLUE,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
)

# Physical dimensions and coordinate grid (6.6 units = 1 mm, 462x990)
FIG_WIDTH_MM: Final[str] = "70mm"
FIG_HEIGHT_MM: Final[str] = "150mm"
VIEWBOX_W: Final[int] = 462
VIEWBOX_H: Final[int] = 990
MIN_FONT_SIZE_PT: Final[float] = 10.0


class CompartmentCardSpec(BaseModel):
    """Immutable specification for a simplified TME compartment callout card."""

    model_config = ConfigDict(frozen=True)
    badge_num: str
    title: str
    subtitle: str
    accent_color: str
    y_pos: float
    target_x: float
    target_y: float


class CancerCohortSpec(BaseModel):
    """Immutable specification for a cancer indication row in Panel 3."""

    model_config = ConfigDict(frozen=True)
    label: str
    accent_color: str
    text_color: str
    y_pos: float


class PredictorCategorySpec(BaseModel):
    """Immutable specification for a baseline predictor category card in Panel 3."""

    model_config = ConfigDict(frozen=True)
    title: str
    accent_color: str
    y_pos: float


class HeatmapRowSpec(BaseModel):
    """Immutable specification for a row in the miniature PR-AUC matrix schematic."""

    model_config = ConfigDict(frozen=True)
    label: str
    values: tuple[str, ...]
    scores: tuple[float, ...]


COMPARTMENT_CARDS: Final[tuple[CompartmentCardSpec, ...]] = (
    CompartmentCardSpec(
        badge_num="1",
        title="1. Malignant Cells",
        subtitle="Proliferating &amp; hypoxic tumor",
        accent_color="#C2410C",
        y_pos=40.0,
        target_x=110.0,
        target_y=88.0,
    ),
    CompartmentCardSpec(
        badge_num="2",
        title="2. Effector Immune",
        subtitle="CD8+ T, NK &amp; B lymphocytes",
        accent_color=OKABE_BLUE,
        y_pos=97.0,
        target_x=148.0,
        target_y=125.0,
    ),
    CompartmentCardSpec(
        badge_num="3",
        title="3. Suppressive Immune",
        subtitle="Tregs, M2 TAMs &amp; MDSCs",
        accent_color="#6C32A6",
        y_pos=154.0,
        target_x=130.0,
        target_y=168.0,
    ),
    CompartmentCardSpec(
        badge_num="4",
        title="4. Stroma &amp; Vessels",
        subtitle="CAFs, collagen ECM &amp; vessels",
        accent_color="#8C6239",
        y_pos=211.0,
        target_x=96.0,
        target_y=225.0,
    ),
)


CANCER_COHORTS: Final[tuple[CancerCohortSpec, ...]] = (
    CancerCohortSpec(
        label="Clear Cell Renal (n = 263)",
        accent_color=OKABE_BLUISH_GREEN,
        text_color=OKABE_BLUISH_GREEN,
        y_pos=30.0,
    ),
    CancerCohortSpec(
        label="Melanoma (n = 338)",
        accent_color=OKABE_ORANGE,
        text_color=OKABE_ORANGE,
        y_pos=58.0,
    ),
    CancerCohortSpec(
        label="Pancreatic Ductal (n = 85)",
        accent_color=OKABE_REDDISH_PURPLE,
        text_color=OKABE_REDDISH_PURPLE,
        y_pos=86.0,
    ),
    CancerCohortSpec(
        label="Triple-Neg. Breast (n = 31)",
        accent_color=OKABE_VERMILION,
        text_color=OKABE_VERMILION,
        y_pos=114.0,
    ),
    CancerCohortSpec(
        label="Urothelial Bladder (n = 298)",
        accent_color=OKABE_SKY_BLUE,
        text_color="#0284C7",
        y_pos=142.0,
    ),
)


PREDICTOR_CATEGORIES: Final[tuple[PredictorCategorySpec, ...]] = (
    PredictorCategorySpec(
        title="• 16 Multi-Gene Signatures",
        accent_color=OKABE_BLUISH_GREEN,
        y_pos=33.0,
    ),
    PredictorCategorySpec(
        title="• 5 Single-Gene Checkpoints",
        accent_color=OKABE_BLUE,
        y_pos=86.5,
    ),
    PredictorCategorySpec(
        title="• 6 Systems &amp; Stromal Models",
        accent_color=OKABE_REDDISH_PURPLE,
        y_pos=140.0,
    ),
)


HEATMAP_ROWS: Final[tuple[HeatmapRowSpec, ...]] = (
    HeatmapRowSpec(
        label="Classifier 1",
        values=(".52**", ".78**", ".52**", ".72**", ".62**", ".58**"),
        scores=(0.52, 0.78, 0.52, 0.72, 0.62, 0.58),
    ),
    HeatmapRowSpec(
        label="Classifier 2",
        values=(".26", ".48", ".37", ".86**", ".63*", ".54**"),
        scores=(0.26, 0.48, 0.37, 0.86, 0.63, 0.54),
    ),
    HeatmapRowSpec(
        label="Classifier 3",
        values=(".26", ".50", ".38", ".84**", ".64*", ".55**"),
        scores=(0.26, 0.50, 0.38, 0.84, 0.64, 0.55),
    ),
    HeatmapRowSpec(
        label="Classifier 4",
        values=(".24", ".50", ".31", ".84**", ".58", ".42**"),
        scores=(0.24, 0.50, 0.31, 0.84, 0.58, 0.42),
    ),
    HeatmapRowSpec(
        label="Classifier 5",
        values=(".12", ".71", ".16", ".50", ".21", ".33"),
        scores=(0.12, 0.71, 0.16, 0.50, 0.21, 0.33),
    ),
)


def pr_auc_to_color(score: float) -> tuple[str, str]:
    """Map PR-AUC score in [0.10, 0.85] to diverging Red-White-Blue fill & text color."""
    match score:
        case s if s >= 0.70:
            return ("#1B5E91", "#FFFFFF")
        case s if s >= 0.58:
            return ("#4386B8", "#FFFFFF")
        case s if s >= 0.50:
            return ("#87BBD9", "#0F172A")
        case s if s >= 0.42:
            return ("#CFE4F0", "#0F172A")
        case s if s >= 0.36:
            return ("#F1F5F9", "#334155")
        case s if s >= 0.28:
            return ("#F9D5C4", "#0F172A")
        case s if s >= 0.20:
            return ("#E6866A", "#FFFFFF")
        case _:
            return ("#B82E2E", "#FFFFFF")


def extract_tme_asset_parts(asset_path: Path, prefix: str = "tme_eco_") -> Result[tuple[str, str], str]:
    """Extract and namespace <defs> children and artwork <g> layers from the TME SVG asset."""
    if not asset_path.exists():
        return Failure(f"TME asset not found at {asset_path}")

    raw_svg = asset_path.read_text(encoding="utf-8")

    # Collect all id="..." declarations inside the TME asset and prefix them to prevent collisions
    ids_found = set(re.findall(r'\bid="([^"]+)"', raw_svg))
    namespaced_svg = raw_svg
    for old_id in sorted(ids_found, key=len, reverse=True):
        new_id = f"{prefix}{old_id}"
        namespaced_svg = namespaced_svg.replace(f'id="{old_id}"', f'id="{new_id}"')
        namespaced_svg = namespaced_svg.replace(f"url(#{old_id})", f"url(#{new_id})")
        namespaced_svg = namespaced_svg.replace(f'href="#{old_id}"', f'href="#{new_id}"')

    # Extract <defs> inner content
    defs_match = re.search(r"<defs>(.*?)</defs>", namespaced_svg, flags=re.DOTALL)
    defs_inner = defs_match.group(1).strip() if defs_match else ""

    # Extract everything after </defs> up to </svg>, stripping the 800x740 background layer
    body_match = re.search(r"</defs>(.*)</svg>", namespaced_svg, flags=re.DOTALL)
    if not body_match:
        return Failure("Could not extract SVG body from TME ecosystem asset.")

    body_content = body_match.group(1)
    # Remove the 00_Canvas_Background layer so its white rectangle does not overflow
    body_clean = re.sub(
        r"<!-- 00_Canvas_Background -->.*?<g[^>]*id=\"[^\"]*layer-00-background\"[^>]*>.*?</g>",
        "",
        body_content,
        flags=re.DOTALL,
    )
    # Convert inkscape:groupmode="layer" on nested asset groups to normal groups so top-level layers stay clean
    body_clean = body_clean.replace('inkscape:groupmode="layer"', "")
    return Success((defs_inner, body_clean.strip()))


def render_compartment_swatch(badge_num: str, cx: float, cy: float) -> str:
    """Render miniature biological cell icon swatch for each of the 4 TME compartments."""
    match badge_num:
        case "1":
            return f"""
            <g transform="translate({cx}, {cy}) scale(1.75)">
              <circle cx="0" cy="0" r="5.6" fill="#F5AAAA" stroke="#C2410C" stroke-width="0.7" />
              <circle cx="-0.5" cy="-0.3" r="2.7" fill="#9C3B3B" />
            </g>"""
        case "2":
            return f"""
            <g transform="translate({cx}, {cy}) scale(1.75)">
              <circle cx="-1.5" cy="-1.2" r="4.5" fill="#3FA4D8" stroke="#0E5482" stroke-width="0.7" />
              <circle cx="-1.5" cy="-1.2" r="2.5" fill="#063250" />
              <circle cx="2.8" cy="2.2" r="3.6" fill="#A4DEC4" stroke="#328761" stroke-width="0.6" />
              <circle cx="2.8" cy="2.2" r="1.9" fill="#328761" />
            </g>"""
        case "3":
            return f"""
            <g transform="translate({cx}, {cy}) scale(1.75)">
              <polygon points="0,-5.8 1.8,-2.2 5.5,-2.8 2.8,0.4 4.6,4.1 0.8,2.6 -2.2,5.1 -2.0,1.2 -5.5,-0.8 -1.8,-2.2"
                       fill="#8F59C7" stroke="#4F1C7E" stroke-width="0.55" />
              <circle cx="0" cy="0" r="2.0" fill="#3B0F60" />
            </g>"""
        case _:
            return f"""
            <g transform="translate({cx}, {cy}) scale(1.75)">
              <path d="M -5.8 2.6 Q -1.5 -2.6 5.8 -2.6 Q 1.5 2.6 -5.8 2.6 Z"
                    fill="#B89B7D" stroke="#6E5034" stroke-width="0.6" />
              <circle cx="0" cy="0" r="1.4" fill="#5A3E24" />
              <ellipse cx="2.2" cy="3.1" rx="3.4" ry="1.6" fill="#801C1C" stroke="#DE6464" stroke-width="0.5" />
            </g>"""


def render_patient_mini(x: float, y: float, color: str, scale: float = 0.96) -> str:
    """Render a compact vector patient silhouette icon."""
    return f"""
    <g transform="translate({x}, {y}) scale({scale})">
      <circle cx="5" cy="3.2" r="2.5" fill="{color}" />
      <path d="M 1.5 7.2 C 1.5 5.8, 3.0 5.4, 5.0 5.4 C 7.0 5.4, 8.5 5.8, 8.5 7.2
               L 8.0 11.5 L 6.8 11.5 L 6.8 16.0 L 5.5 16.0 L 5.5 11.8
               L 4.5 11.8 L 4.5 16.0 L 3.2 16.0 L 3.2 11.5 L 2.0 11.5 Z" fill="{color}" />
    </g>"""


def generate_graphical_abstract_svg(tme_defs: str, tme_body: str) -> str:
    """Construct the complete 70mm x 150mm SVG document with native Inkscape layers and font-size >= 10."""
    parts: list[str] = []

    # XML & SVG Root
    parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{FIG_WIDTH_MM}" height="{FIG_HEIGHT_MM}"
     viewBox="0 0 {VIEWBOX_W} {VIEWBOX_H}"
     version="1.1" id="fig-tme-pipeline-graphical-abstract"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  <sodipodi:namedview
     id="namedview-tme-ga"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:document-units="mm" />
  <defs id="defs-main">
    <marker id="ga-arr-slate" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
    </marker>
    <marker id="ga-arr-blue" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
    <marker id="ga-arr-green" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUISH_GREEN}" />
    </marker>
    <clipPath id="ga-clip-tme-asset">
      <rect x="13.0" y="41.0" width="196.0" height="220.0" />
    </clipPath>
    {tme_defs}
  </defs>
""")

    # =========================================================================
    # LAYER 00: Canvas Background & Wireframe Containers
    # =========================================================================
    parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND & PANEL CONTAINERS -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{VIEWBOX_W}" height="{VIEWBOX_H}" fill="{COLOR_CANVAS_BG}" />
    <rect id="card-panel-a" x="5" y="5" width="452" height="266" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.2" />
    <rect id="card-panel-b" x="5" y="285" width="452" height="214" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.2" />
    <rect id="card-panel-c" x="5" y="513" width="452" height="232" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.2" />
    <rect id="card-panel-d" x="5" y="759" width="452" height="226" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.2" />
  </g>
""")

    # =========================================================================
    # LAYER 01: Panel Section Banners (Uniform neutral styling)
    # =========================================================================
    parts.append(f"""
  <!-- LAYER 01: PANEL BANNERS -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Panel_Banners">
    <!-- Panel 1 Banner -->
    <g id="banner-panel-a">
      <rect x="5" y="5" width="452" height="27.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <text x="12" y="23.5" font-size="13.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        TUMOR MICROENVIRONMENT ECOSYSTEM
      </text>
    </g>

    <!-- Panel 2 Banner -->
    <g id="banner-panel-b">
      <rect x="5" y="285" width="452" height="27.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <text x="12" y="303.5" font-size="13.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        DECONVOLUTION INFERENCE PIPELINE
      </text>
    </g>

    <!-- Panel 3 Banner -->
    <g id="banner-panel-c">
      <rect x="5" y="513" width="452" height="27.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <text x="12" y="531.5" font-size="12.6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        9 PRE-TREATMENT ICB COHORTS &amp; 27 LITERATURE BIOMARKERS
      </text>
    </g>

    <!-- Panel 4 Banner -->
    <g id="banner-panel-d">
      <rect x="5" y="759" width="452" height="27.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <text x="12" y="777.5" font-size="13.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        PR-AUC MATRIX BENCHMARK
      </text>
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 02: PANEL 1 — TME Ecosystem Asset & Simplified Compartments
    # =========================================================================
    compartment_cards_svg: list[str] = []
    callout_pins_svg: list[str] = []

    for spec in COMPARTMENT_CARDS:
        swatch_svg = render_compartment_swatch(spec.badge_num, 18.0, 25.5)
        compartment_cards_svg.append(f"""
      <g id="grp-tme-card-{spec.badge_num}" transform="translate(218, {spec.y_pos})">
        <rect x="0" y="0" width="232" height="51.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
        <rect x="0" y="0" width="4.5" height="51.0" fill="{spec.accent_color}" />
        {swatch_svg}
        <text x="34" y="21.5" font-size="12.2" font-weight="700" fill="{spec.accent_color}">{spec.title}</text>
        <text x="34" y="39.5" font-size="10.5" fill="{COLOR_TEXT_SECONDARY}">{spec.subtitle}</text>
      </g>""")

        card_mid_y = spec.y_pos + 25.5
        callout_pins_svg.append(f"""
      <g id="grp-tme-pin-{spec.badge_num}">
        <path d="M {spec.target_x} {spec.target_y} L 204 {card_mid_y} L 218 {card_mid_y}"
              stroke="{spec.accent_color}" stroke-width="1.1" stroke-dasharray="2.5,1.8" fill="none" />
        <circle cx="{spec.target_x}" cy="{spec.target_y}" r="8.8" fill="{spec.accent_color}" stroke="#FFFFFF" stroke-width="1.1" />
        <text x="{spec.target_x}" y="{spec.target_y + 3.6}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{spec.badge_num}</text>
      </g>""")

    parts.append(f"""
  <!-- LAYER 02: PANEL 1 — TUMOR MICROENVIRONMENT ECOSYSTEM -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_1_TME_Ecosystem">
    <!-- Left Frame: Embedded TME Vector Asset -->
    <rect id="frame-tme-asset" x="12" y="40.0" width="198" height="222.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
    <g id="grp-tme-artwork-clipped" clip-path="url(#ga-clip-tme-asset)">
      <g id="grp-tme-artwork" inkscape:label="Embedded_TME_Ecosystem_Asset" transform="translate(13.0, 48.0) scale(0.245)">
        {tme_body}
      </g>
    </g>
    <rect x="15.0" y="43.0" width="128" height="20.0" fill="#FFFFFF" opacity="0.92" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
    <text x="79.0" y="57.0" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
      Heterogeneous TME
    </text>

    <!-- Callout Leader Pins & 4 Compartment Cards -->
    {"".join(callout_pins_svg)}
    {"".join(compartment_cards_svg)}
  </g>
""")

    # =========================================================================
    # LAYER 03: PANEL 2 — Deconvolution Inference Pipeline
    # =========================================================================
    parts.append(f"""
  <!-- LAYER 03: PANEL 2 — DECONVOLUTION INFERENCE PIPELINE -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_2_Deconv_Pipeline">
    <!-- Module B1: Single-Cell Atlas Reference -->
    <g id="grp-module-b1" transform="translate(12, 318.0)">
      <rect x="0" y="0" width="139" height="174.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <rect x="0" y="0" width="139" height="23.0" fill="#E2E8F0" />
      <text x="69.5" y="16.2" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        1. scRNA Reference
      </text>

      <!-- Mini UMAP + Reference Matrix Schematic -->
      <g transform="translate(6, 28.0)">
        <rect x="0" y="0" width="127" height="60.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
        <!-- Mini UMAP cluster dots -->
        <circle cx="18" cy="17.0" r="3.8" fill="{OKABE_BLUE}" />
        <circle cx="25" cy="13.5" r="3.2" fill="{OKABE_BLUE}" />
        <circle cx="22" cy="24.5" r="3.5" fill="{OKABE_BLUE}" />
        <circle cx="39" cy="15.5" r="3.6" fill="{OKABE_BLUISH_GREEN}" />
        <circle cx="44" cy="22.5" r="3.0" fill="{OKABE_BLUISH_GREEN}" />
        <circle cx="29" cy="35.0" r="3.8" fill="{OKABE_ORANGE}" />
        <circle cx="36" cy="38.5" r="3.2" fill="{OKABE_ORANGE}" />
        <circle cx="15" cy="36.5" r="3.5" fill="{OKABE_VERMILION}" />
        <text x="29" y="54.5" font-size="10.0" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Leiden</text>

        <path d="M 52 26.5 L 65 26.5" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.0" marker-end="url(#ga-arr-slate)" />

        <!-- Mini Signature Matrix Phi -->
        <g transform="translate(69, 5.5)">
          <rect x="0" y="0" width="50" height="32.0" fill="#F8FAFC" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
          <rect x="3.5" y="3.5" width="9.0" height="6.5" fill="{OKABE_BLUE}" />
          <rect x="14.5" y="3.5" width="9.0" height="6.5" fill="#CFE4F0" />
          <rect x="25.5" y="3.5" width="9.0" height="6.5" fill="#F1F5F9" />
          <rect x="36.5" y="3.5" width="10.0" height="6.5" fill="#F9D5C4" />
          <rect x="3.5" y="12.5" width="9.0" height="6.5" fill="#F1F5F9" />
          <rect x="14.5" y="12.5" width="9.0" height="6.5" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="25.5" y="12.5" width="9.0" height="6.5" fill="#CFE4F0" />
          <rect x="36.5" y="12.5" width="10.0" height="6.5" fill="#F1F5F9" />
          <rect x="3.5" y="21.5" width="9.0" height="6.5" fill="#F9D5C4" />
          <rect x="14.5" y="21.5" width="9.0" height="6.5" fill="#F1F5F9" />
          <rect x="25.5" y="21.5" width="9.0" height="6.5" fill="{OKABE_ORANGE}" />
          <rect x="36.5" y="21.5" width="10.0" height="6.5" fill="{OKABE_VERMILION}" />
          <text x="25" y="49.0" font-size="10.2" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
            Φ (G×K)
          </text>
        </g>
      </g>

      <!-- Bullet Steps -->
      <text x="5.5" y="109.0" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">• Solid-Tumor Atlases</text>
      <text x="9.0" y="124.5" font-size="10.0" fill="{COLOR_TEXT_SECONDARY}">Matched .h5ad cohorts</text>
      <text x="5.5" y="146.5" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">• Harmony &amp; Leiden</text>
      <text x="9.0" y="162.0" font-size="10.0" fill="{COLOR_TEXT_SECONDARY}">Cell state clusters K</text>
    </g>

    <!-- Arrow B1 -> B2 -->
    <path d="M 152.0 405.0 L 160.0 405.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />

    <!-- Module B2: Bayesian Deconvolution (BayesPrism Optical Prism Illustration) -->
    <g id="grp-module-b2" transform="translate(161.5, 318.0)">
      <rect x="0" y="0" width="139" height="174.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <rect x="0" y="0" width="139" height="23.0" fill="#E2E8F0" />
      <text x="69.5" y="16.2" font-size="10.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" letter-spacing="-0.3">
        2. Bayesian Deconvolution
      </text>

      <!-- BayesPrism Dispersion Prism & Stacked Fractions Visual -->
      <g transform="translate(6, 28.0)">
        <rect x="0" y="0" width="127" height="60.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
        <!-- Incoming mixed bulk beam into prism left face -->
        <line x1="4.5" y1="28.0" x2="20.0" y2="25.0" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="2.2" />
        <line x1="4.5" y1="28.0" x2="20.0" y2="25.0" stroke="{OKABE_ORANGE}" stroke-width="0.9" stroke-dasharray="2.0,1.5" />
        <!-- BayesPrism Optical Prism -->
        <polygon points="28,7.0 12,43.0 44,43.0" fill="#EFF6FF" stroke="{OKABE_BLUE}" stroke-width="1.2" />
        <line x1="20.0" y1="25.0" x2="34.0" y2="21.0" stroke="#93C5FD" stroke-width="0.9" />
        <line x1="20.0" y1="25.0" x2="36.5" y2="28.0" stroke="#93C5FD" stroke-width="0.9" />
        <line x1="20.0" y1="25.0" x2="39.5" y2="35.0" stroke="#93C5FD" stroke-width="0.9" />
        <!-- Dispersed Cell-Lineage Rays Exiting Right Face of Prism -->
        <line x1="33.0" y1="18.5" x2="64.0" y2="10.0" stroke="{OKABE_BLUE}" stroke-width="1.5" />
        <line x1="35.5" y1="24.0" x2="64.0" y2="18.0" stroke="{OKABE_BLUISH_GREEN}" stroke-width="1.5" />
        <line x1="38.0" y1="29.5" x2="64.0" y2="25.5" stroke="{OKABE_ORANGE}" stroke-width="1.5" />
        <line x1="40.5" y1="35.0" x2="64.0" y2="33.5" stroke="{OKABE_VERMILION}" stroke-width="1.5" />

        <!-- Stacked Cell State Fractions Bar -->
        <g transform="translate(68, 5.5)">
          <rect x="0" y="0" width="49" height="8.5" fill="{OKABE_BLUE}" />
          <rect x="0" y="8.5" width="49" height="7.5" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="0" y="16.0" width="49" height="7.5" fill="{OKABE_ORANGE}" />
          <rect x="0" y="23.5" width="49" height="8.5" fill="{OKABE_VERMILION}" />
          <rect x="0" y="0" width="49" height="32.0" fill="none" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.7" />
          <text x="24.5" y="49.0" font-size="10.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Fractions θ</text>
        </g>
      </g>

      <!-- Bullet Steps -->
      <text x="5.5" y="109.0" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">• BayesPrism / Insta</text>
      <text x="9.0" y="124.5" font-size="10.0" fill="{COLOR_TEXT_SECONDARY}">Posterior abundance θ</text>
      <text x="5.5" y="146.5" font-size="10.0" font-weight="700" fill="{OKABE_BLUISH_GREEN}">• Regularizes sibling</text>
      <text x="9.0" y="162.0" font-size="10.0" font-weight="700" fill="{OKABE_BLUISH_GREEN}">cell state cancelation</text>
    </g>

    <!-- Arrow B2 -> B3 -->
    <path d="M 301.5 405.0 L 309.5 405.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />

    <!-- Module B3: Compositional CLR & Response Classifier (Logistic Regression & Random Forest) -->
    <g id="grp-module-b3" transform="translate(311, 318.0)">
      <rect x="0" y="0" width="139" height="174.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <rect x="0" y="0" width="139" height="23.0" fill="#E2E8F0" />
      <text x="69.5" y="16.2" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        3. Response
      </text>

      <!-- CLR + Logistic Regression & Random Forest Visual -->
      <g transform="translate(6, 28.0)">
        <rect x="0" y="0" width="127" height="60.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
        <!-- Diverging CLR feature bars -->
        <line x1="20" y1="6.0" x2="20" y2="38.5" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.9" />
        <rect x="20" y="8.5" width="14" height="5.8" fill="{OKABE_BLUE}" />
        <rect x="20" y="17.5" width="10" height="5.8" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="8" y="26.5" width="12" height="5.8" fill="{OKABE_VERMILION}" />
        <text x="21" y="54.5" font-size="10.0" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">CLR(θ)</text>

        <path d="M 37 25.5 L 47 25.5" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.0" marker-end="url(#ga-arr-slate)" />

        <!-- Mini Icon 1: Logistic Regression Sigmoid Curve -->
        <g id="grp-icon-logreg" transform="translate(49, 5.5)">
          <rect x="0" y="0" width="35" height="48.0" fill="#F8FAFC" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
          <line x1="4.5" y1="6.0" x2="4.5" y2="41.0" stroke="{COLOR_TEXT_MUTED}" stroke-width="0.7" />
          <line x1="4.5" y1="41.0" x2="31.5" y2="41.0" stroke="{COLOR_TEXT_MUTED}" stroke-width="0.7" />
          <line x1="4.5" y1="23.5" x2="31.5" y2="23.5" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" stroke-dasharray="1.8,1.4" />
          <circle cx="8.5" cy="38.0" r="1.6" fill="{OKABE_VERMILION}" />
          <circle cx="12.5" cy="38.0" r="1.6" fill="{OKABE_VERMILION}" />
          <circle cx="23.5" cy="9.0" r="1.6" fill="{OKABE_BLUISH_GREEN}" />
          <circle cx="28.0" cy="9.0" r="1.6" fill="{OKABE_BLUISH_GREEN}" />
          <path d="M 5.5 38.0 C 15.0 38.0, 16.0 9.0, 30.5 9.0" fill="none" stroke="{OKABE_BLUE}" stroke-width="1.6" />
        </g>

        <!-- Mini Icon 2: Random Forest Decision Trees -->
        <g id="grp-icon-rf" transform="translate(87, 5.5)">
          <rect x="0" y="0" width="36" height="48.0" fill="#F8FAFC" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
          <line x1="18.0" y1="8.0" x2="9.5" y2="21.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.0" />
          <line x1="18.0" y1="8.0" x2="26.5" y2="21.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.0" />
          <line x1="9.5" y1="21.0" x2="5.0" y2="35.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.9" />
          <line x1="9.5" y1="21.0" x2="14.0" y2="35.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.9" />
          <line x1="26.5" y1="21.0" x2="22.0" y2="35.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.9" />
          <line x1="26.5" y1="21.0" x2="31.0" y2="35.0" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.9" />
          <circle cx="18.0" cy="8.0" r="2.6" fill="{OKABE_BLUE}" />
          <circle cx="9.5" cy="21.0" r="2.3" fill="{OKABE_ORANGE}" />
          <circle cx="26.5" cy="21.0" r="2.3" fill="{OKABE_ORANGE}" />
          <rect x="2.8" y="34.0" width="4.4" height="4.4" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="11.8" y="34.0" width="4.4" height="4.4" fill="{OKABE_VERMILION}" />
          <rect x="19.8" y="34.0" width="4.4" height="4.4" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="28.8" y="34.0" width="4.4" height="4.4" fill="{OKABE_BLUISH_GREEN}" />
        </g>
      </g>

      <!-- Bullet Steps -->
      <text x="5.5" y="109.0" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">• Ensemble Classifier</text>
      <text x="9.0" y="124.5" font-size="10.0" fill="{COLOR_TEXT_SECONDARY}">ExtraTrees / Logistic</text>
      <text x="5.5" y="148.5" font-size="10.0" font-weight="700" fill="#000000">• Pooled Out-of-Fold CV</text>
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 04: PANEL 3 — 9 Clinical ICB Cohorts & 27 Literature Predictors
    # =========================================================================
    cohort_rows_svg: list[str] = []
    for idx, c_spec in enumerate(CANCER_COHORTS):
        cohort_rows_svg.append(f"""
      <g id="grp-cancer-cohort-{idx + 1}" transform="translate(6.0, {c_spec.y_pos})">
        <rect x="0" y="0" width="202" height="24.5" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
        <rect x="0" y="0" width="4.2" height="24.5" fill="{c_spec.accent_color}" />
        {render_patient_mini(7.0, 3.2, c_spec.accent_color, 0.96)}
        <text x="20.0" y="16.2" font-size="10.8" font-weight="700" fill="{c_spec.text_color}">{c_spec.label}</text>
      </g>""")

    predictor_cat_svg: list[str] = []
    for idx, p_spec in enumerate(PREDICTOR_CATEGORIES):
        predictor_cat_svg.append(f"""
      <g id="grp-pred-cat-{idx + 1}" transform="translate(6.0, {p_spec.y_pos})">
        <rect x="0" y="0" width="204" height="45.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.7" />
        <rect x="0" y="0" width="4.2" height="45.0" fill="{p_spec.accent_color}" />
        <text x="11.0" y="26.8" font-size="11.0" font-weight="700" fill="{p_spec.accent_color}">
          {p_spec.title}
        </text>
      </g>""")

    parts.append(f"""
  <!-- LAYER 04: PANEL 3 — 9 CLINICAL COHORTS & 27 LITERATURE PREDICTORS -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_3_Cohorts_Predictors">
    <!-- Left Sub-Card: 9 Pre-Treatment Clinical ICB Cohorts (N = 1,015) -->
    <g id="grp-panel-c-cohorts" transform="translate(12, 546.0)">
      <rect x="0" y="0" width="214" height="192.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <rect x="0" y="0" width="214" height="24.0" fill="#E2E8F0" />
      <text x="107" y="16.8" font-size="11.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        9 Cohorts (N = 1,015)
      </text>

      {"".join(cohort_rows_svg)}

      <!-- Binary Endpoint Footer -->
      <text x="107" y="184.0" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        RECIST: <tspan fill="{OKABE_BLUISH_GREEN}">319 Resp</tspan> vs. <tspan fill="{OKABE_VERMILION}">696 Non-Resp</tspan>
      </text>
    </g>

    <!-- Right Sub-Card: 27 Published Transcriptomic Predictors (3 Categories Only) -->
    <g id="grp-panel-c-predictors" transform="translate(234, 546.0)">
      <rect x="0" y="0" width="216" height="192.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
      <rect x="0" y="0" width="216" height="24.0" fill="#E2E8F0" />
      <text x="108" y="16.8" font-size="11.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        27 Baseline Predictors
      </text>

      {"".join(predictor_cat_svg)}
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 05: PANEL 4 — Complete-Matrix PR-AUC Benchmark
    # =========================================================================
    col_headers = ("Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4", "Cohort 5", "Cohort 6")
    col_x_positions: list[float] = [88.0 + c_idx * 57.0 for c_idx in range(len(col_headers))]

    header_cells_svg: list[str] = []
    for c_idx, cname in enumerate(col_headers):
        cx = col_x_positions[c_idx] + 27.5
        header_cells_svg.append(
            f'<text x="{cx}" y="40.0" font-size="10.2" font-weight="700" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">{cname}</text>'
        )

    matrix_rows_svg: list[str] = []
    for r_idx, r_spec in enumerate(HEATMAP_ROWS):
        ry = 47.0 + r_idx * 25.0
        row_cells: list[str] = [
            f'<text x="83.0" y="{ry + 15.8}" font-size="10.6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{r_spec.label}</text>'
        ]
        for c_idx, (val_str, num_score) in enumerate(zip(r_spec.values, r_spec.scores)):
            cell_x = col_x_positions[c_idx]
            bg_col, fg_col = pr_auc_to_color(num_score)
            row_cells.append(
                f'<rect x="{cell_x}" y="{ry}" width="55.0" height="22.5" fill="{bg_col}" stroke="#FFFFFF" stroke-width="0.6" />'
                f'<text x="{cell_x + 27.5}" y="{ry + 15.5}" font-size="10.6" font-weight="700" fill="{fg_col}" text-anchor="middle">{val_str}</text>'
            )
        matrix_rows_svg.append(f'<g id="grp-hm-row-{r_idx}">{"".join(row_cells)}</g>')

    parts.append(f"""
  <!-- LAYER 05: PANEL 4 — PR-AUC BENCHMARK MATRIX -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_4_Pooling_PRAUC">
    <!-- Full-Panel Complete-Matrix PR-AUC Heatmap Schematic -->
    <g id="grp-panel-d-heatmap" transform="translate(12, 793.0)">
      <rect x="0" y="0" width="438" height="184.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />

      <!-- Title -->
      <text x="8.0" y="20.0" font-size="11.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Pre-Treatment Standard ICB Response (Complete Matrix, PR-AUC)
      </text>

      {"".join(header_cells_svg)}
      {"".join(matrix_rows_svg)}
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 06: Gutter-Constrained Inter-Panel Flow Connectors
    # =========================================================================
    parts.append(f"""
  <!-- LAYER 06: INTER-PANEL FLOW CONNECTORS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Flow_Connectors">
    <!-- Panel 1 -> Panel 2 -->
    <path id="conn-a-to-b" d="M 231 271 L 231 283.5" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />
    <!-- Panel 2 -> Panel 3 -->
    <path id="conn-b-to-c" d="M 231 499 L 231 511.5" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />
    <!-- Panel 3 -> Panel 4 (Left Cohort Stream & Right Predictor Stream) -->
    <path id="conn-c-to-d-left" d="M 119 745 L 119 757.5" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />
    <path id="conn-c-to-d-right" d="M 342 745 L 342 757.5" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.3" marker-end="url(#ga-arr-slate)" />
  </g>
</svg>
""")

    return "".join(parts)


def validate_svg(svg_text: str) -> Result[int, str]:
    """Validate XML well-formedness, Inkscape layers, uniqueness of IDs, and minimum font size >= 10."""
    try:
        root = ET.fromstring(svg_text)
    except ET.ParseError as exc:
        return Failure(f"XML parse error: {exc}")

    if root.attrib.get("width") != FIG_WIDTH_MM or root.attrib.get("height") != FIG_HEIGHT_MM:
        return Failure(
            f"Dimension mismatch: width={root.attrib.get('width')}, height={root.attrib.get('height')}"
        )

    inkscape_ns = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
    layers = [el for el in root.iter() if el.attrib.get(inkscape_ns) == "layer"]
    if len(layers) < 6:
        return Failure(f"Expected at least 6 Inkscape layers, found {len(layers)}")

    seen_ids: set[str] = set()
    dup_ids: list[str] = []
    small_fonts: list[str] = []
    for el in root.iter():
        eid = el.attrib.get("id", "")
        if eid:
            if eid in seen_ids:
                dup_ids.append(eid)
            seen_ids.add(eid)
        fs_raw = el.attrib.get("font-size", "")
        if fs_raw:
            fs_val = float(re.sub(r"[a-zA-Z%]+$", "", fs_raw))
            if fs_val < MIN_FONT_SIZE_PT:
                small_fonts.append(f"{el.tag}(font-size={fs_raw})")

    if dup_ids:
        return Failure(f"Duplicate XML IDs detected: {dup_ids[:10]}")
    if small_fonts:
        return Failure(f"Font sizes smaller than {MIN_FONT_SIZE_PT} detected: {small_fonts[:10]}")

    return Success(len(layers))


def embed_into_figure_01_overview(overview_path: Path, standalone_svg: str) -> Result[Path, str]:
    """Embed the 70mm x 150mm graphical abstract into the top-left slot of Page 1 in figure_01_overview.svg."""
    if not overview_path.exists():
        return Failure(f"Overview SVG not found at {overview_path}")

    # Prefix all IDs in standalone_svg with 'p1_ga_' so there is zero collision inside figure_01_overview.svg
    all_ids = set(re.findall(r'\bid="([^"]+)"', standalone_svg))
    prefixed_svg = standalone_svg
    for old_id in sorted(all_ids, key=len, reverse=True):
        new_id = f"p1_ga_{old_id}"
        prefixed_svg = prefixed_svg.replace(f'id="{old_id}"', f'id="{new_id}"')
        prefixed_svg = prefixed_svg.replace(f"url(#{old_id})", f"url(#{new_id})")
        prefixed_svg = prefixed_svg.replace(f'href="#{old_id}"', f'href="#{new_id}"')

    # Extract <defs> and body from prefixed_svg
    defs_match = re.search(r"<defs[^>]*>(.*?)</defs>", prefixed_svg, flags=re.DOTALL)
    body_match = re.search(r"</defs>(.*)</svg>", prefixed_svg, flags=re.DOTALL)
    if not defs_match or not body_match:
        return Failure("Failed to extract defs/body for figure_01_overview embedding.")

    defs_block = f'<g id="p1_ga_defs_container" style="display:none"><defs id="p1_ga_defs">{defs_match.group(1)}</defs></g>'
    body_block = body_match.group(1).replace('inkscape:groupmode="layer"', "")

    scale_in_g180 = (70.0 / float(VIEWBOX_W)) / 0.98714269
    wrapper_group = (
        f'<g id="g-page01-tme-pipeline-abstract" inkscape:label="Page01_TME_Pipeline_Graphical_Abstract" '
        f'transform="matrix({scale_in_g180:.8f},0,0,{scale_in_g180:.8f},-13.193823,-3.488709)">'
        f"{defs_block}{body_block}</g>"
    )

    overview_text = overview_path.read_text(encoding="utf-8")

    # Remove any prior embedded version if re-running idempotently
    overview_text = re.sub(
        r'<g id="g-page01-tme-pipeline-abstract".*?</g>\s*(?=<g\s+id="g186"|<g id="g186")',
        "",
        overview_text,
        flags=re.DOTALL,
    )

    target_marker = '<g\n         id="g186"'
    if target_marker not in overview_text:
        target_marker = 'id="g186"'
        idx = overview_text.find(target_marker)
        if idx == -1:
            return Failure("Could not locate g186 in figure_01_overview.svg")
        g_start = overview_text.rfind("<g", 0, idx)
        updated_overview = overview_text[:g_start] + wrapper_group + "\n      " + overview_text[g_start:]
    else:
        updated_overview = overview_text.replace(target_marker, wrapper_group + "\n      " + target_marker, 1)

    overview_path.write_text(updated_overview, encoding="utf-8")
    return Success(overview_path)


def render_png_preview(svg_path: Path, png_path: Path, width_px: int, height_px: int) -> Result[Path, str]:
    """Render 300 DPI PNG raster using rsvg-convert or headless vl_convert fallback."""
    from plotting_utils import rasterize_svg_to_png

    return rasterize_svg_to_png(
        svg_path=svg_path,
        png_path=png_path,
        width_px=width_px,
        height_px=height_px,
        background=None,
    )


def main() -> int:
    """Imperative shell for generating, validating, and exporting the graphical abstract SVG."""
    repo_root = Path(__file__).resolve().parent.parent.parent
    tme_asset_path = repo_root / "assets" / "tme_tumor_microenvironment_ecosystem.svg"

    out_svg_paths: Sequence[Path] = (
        repo_root / "assets" / "tme_pipeline_graphical_abstract.svg",
        repo_root / "article" / "figures" / "deconvolution" / "tme_pipeline_graphical_abstract.svg",
    )
    out_png_paths: Sequence[Path] = (
        repo_root / "assets" / "rendered_png" / "tme_pipeline_graphical_abstract.png",
        repo_root / "article" / "figures" / "deconvolution" / "tme_pipeline_graphical_abstract.png",
    )
    overview_svg_path = repo_root / "article" / "figures" / "inkscape" / "figure_01_overview.svg"

    match extract_tme_asset_parts(tme_asset_path):
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            return 1
        case Success((tme_defs, tme_body)):
            svg_content = generate_graphical_abstract_svg(tme_defs, tme_body)

    match validate_svg(svg_content):
        case Failure(err):
            print(f"[VALIDATION ERROR] {err}", file=sys.stderr)
            return 1
        case Success(n_layers):
            print(
                f"[OK] Validated SVG XML ({FIG_WIDTH_MM} x {FIG_HEIGHT_MM}) with {n_layers} Inkscape layers "
                f"and minimum font-size >= {MIN_FONT_SIZE_PT} pt."
            )

    for svg_p in out_svg_paths:
        svg_p.parent.mkdir(parents=True, exist_ok=True)
        svg_p.write_text(svg_content, encoding="utf-8")
        print(f"[SAVED SVG] {svg_p}")

    for png_p in out_png_paths:
        match render_png_preview(out_svg_paths[0], png_p, width_px=1120, height_px=2400):
            case Failure(err):
                print(f"[WARN] {err}", file=sys.stderr)
            case Success(p):
                print(f"[SAVED PNG] {p}")

    match embed_into_figure_01_overview(overview_svg_path, svg_content):
        case Failure(err):
            print(f"[WARN] Could not embed into figure_01_overview.svg: {err}", file=sys.stderr)
        case Success(ov_p):
            print(f"[UPDATED OVERVIEW] Embedded 70mm x 150mm panel into {ov_p}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
