#!/usr/bin/env python3
"""
Generate Graphical Abstract for Section 6: Multi-Cohort ICB Benchmarking.

Replicates the visual and structural design style of the COMPASS reference figure:
- Double-column publication landscape (1400 x 930 px)
- Modular panels (a, b, c, d, e) with bold lowercase editorial letter badges
- Patient silhouette pictogram grids color-coded by cancer indications with response badges
- Projection trapezoid funnels and compact matrix heatmap cards
- Mini analytical curve cards: ROC-AUC (Hanley-McNeil CIs), PR-AUC (Excess Precision),
  Harrell's C-index forest, Decision Curve Analysis (DCA Net Benefit), Feature attribution maps
- 4-Tiered Predictor Performance Hierarchy with flame (calibrated) and snowflake (fixed) badges
- Strict 'Keywords Only' constraint: zero prose or full sentences anywhere on the canvas
- Native Inkscape layer architecture (inkscape:groupmode="layer", inkscape:label)
- Nature Methods minimal wireframe standards (pure white canvas #FFFFFF, hairline borders #CBD5E1, rx=0)
"""

from pathlib import Path
import html
import os
import sys

# Import shared styling tokens
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style_config import (
    WIDTH,
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_CONTAINER_BG,
    COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_YELLOW, OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)

HEIGHT = 930

# Cancer Indication Palette (Okabe-Ito Colorblind-Safe)
COLOR_MELANOMA = OKABE_ORANGE         # #E69F00 - Melanoma
COLOR_RCC = OKABE_BLUISH_GREEN        # #009E73 - Clear Cell Renal Cell Carcinoma (ccRCC)
COLOR_BLADDER = OKABE_SKY_BLUE        # #56B4E9 - Urothelial Bladder
COLOR_PDAC = OKABE_REDDISH_PURPLE     # #CC79A7 - Pancreatic Ductal
COLOR_TNBC = OKABE_VERMILION          # #D55E00 - Triple-Negative Breast

COLOR_RESPONDER = OKABE_BLUISH_GREEN   # #009E73 - CR / PR
COLOR_NON_RESPONDER = OKABE_VERMILION  # #D55E00 - SD / PD
COLOR_OUR_MODEL = OKABE_BLUE          # #0072B2 - Deconvolution Inference (Our)


def render_person_icon(x: float, y: float, color: str, badge: str | None = None, scale: float = 1.0) -> str:
    """Render a crisp vector patient silhouette pictogram with optional response badge."""
    badge_markup = ""
    if badge == "R":  # Responder: Green checkmark badge
        badge_markup = f"""
      <circle cx="{x + 13 * scale}" cy="{y + 6 * scale}" r="{3.5 * scale}" fill="{COLOR_RESPONDER}" />
      <path d="M {x + 11.5 * scale} {y + 6 * scale} L {x + 12.7 * scale} {y + 7.2 * scale} L {x + 14.8 * scale} {y + 4.8 * scale}" 
            stroke="#FFFFFF" stroke-width="{0.9 * scale}" fill="none" stroke-linecap="round" stroke-linejoin="round" />
        """
    elif badge == "NR":  # Non-Responder: Vermilion cross badge
        badge_markup = f"""
      <circle cx="{x + 13 * scale}" cy="{y + 6 * scale}" r="{3.5 * scale}" fill="{COLOR_NON_RESPONDER}" />
      <path d="M {x + 11.8 * scale} {y + 4.8 * scale} L {x + 14.2 * scale} {y + 7.2 * scale} M {x + 14.2 * scale} {y + 4.8 * scale} L {x + 11.8 * scale} {y + 7.2 * scale}" 
            stroke="#FFFFFF" stroke-width="{0.9 * scale}" stroke-linecap="round" />
        """

    return f"""
    <g class="patient-glyph">
      <!-- Head -->
      <circle cx="{x + 7 * scale}" cy="{y + 4.2 * scale}" r="{3.2 * scale}" fill="{color}" />
      <!-- Torso & Limbs -->
      <path d="M {x + 2.5 * scale} {y + 9 * scale} 
               C {x + 2.5 * scale} {y + 7.5 * scale}, {x + 4.5 * scale} {y + 7 * scale}, {x + 7 * scale} {y + 7 * scale} 
               C {x + 9.5 * scale} {y + 7 * scale}, {x + 11.5 * scale} {y + 7.5 * scale}, {x + 11.5 * scale} {y + 9 * scale} 
               L {x + 11 * scale} {y + 14.5 * scale} 
               C {x + 11 * scale} {y + 15.2 * scale}, {x + 10.2 * scale} {y + 15.4 * scale}, {x + 9.5 * scale} {y + 15.4 * scale} 
               L {x + 9.2 * scale} {y + 15.4 * scale} 
               L {x + 9.2 * scale} {y + 21 * scale} 
               C {x + 9.2 * scale} {y + 21.6 * scale}, {x + 8.2 * scale} {y + 21.6 * scale}, {x + 8.2 * scale} {y + 21 * scale} 
               L {x + 7.8 * scale} {y + 15.5 * scale} 
               L {x + 6.2 * scale} {y + 15.5 * scale} 
               L {x + 5.8 * scale} {y + 21 * scale} 
               C {x + 5.8 * scale} {y + 21.6 * scale}, {x + 4.8 * scale} {y + 21.6 * scale}, {x + 4.8 * scale} {y + 21 * scale} 
               L {x + 4.8 * scale} {y + 15.4 * scale} 
               L {x + 4.5 * scale} {y + 15.4 * scale} 
               C {x + 3.8 * scale} {y + 15.4 * scale}, {x + 3 * scale} {y + 15.2 * scale}, {x + 3 * scale} {y + 14.5 * scale} Z" 
            fill="{color}" />
      {badge_markup}
    </g>
    """


def render_flame_icon(x: float, y: float, size: float = 16.0) -> str:
    """Render a clean vector flame icon representing learnable/calibrated weights."""
    s = size / 16.0
    return f"""
    <g transform="translate({x}, {y})">
      <path d="M {8 * s} {1 * s} 
               C {8 * s} {4 * s}, {13 * s} {7 * s}, {13 * s} {11 * s} 
               C {13 * s} {14.2 * s}, {10.5 * s} {16 * s}, {7.5 * s} {16 * s} 
               C {4 * s} {16 * s}, {1.5 * s} {14 * s}, {1.5 * s} {10.8 * s} 
               C {1.5 * s} {7.5 * s}, {4.5 * s} {5.5 * s}, {6 * s} {3.5 * s} 
               C {5.8 * s} {5.2 * s}, {7 * s} {7 * s}, {8 * s} {7.2 * s} 
               C {7.2 * s} {5 * s}, {8 * s} {1 * s}, {8 * s} {1 * s} Z" 
            fill="#EA580C" />
      <path d="M {7.8 * s} {8 * s} 
               C {8.8 * s} {9.5 * s}, {10 * s} {11 * s}, {10 * s} {12.8 * s} 
               C {10 * s} {14.5 * s}, {8.8 * s} {15.5 * s}, {7.5 * s} {15.5 * s} 
               C {5.8 * s} {15.5 * s}, {4.5 * s} {14.5 * s}, {4.5 * s} {12.5 * s} 
               C {4.5 * s} {11 * s}, {6 * s} {9.5 * s}, {7.8 * s} {8 * s} Z" 
            fill="#FACC15" />
    </g>
    """


def render_snowflake_icon(x: float, y: float, size: float = 16.0) -> str:
    """Render a clean vector snowflake icon representing fixed/frozen weights."""
    s = size / 16.0
    cx = 8 * s
    cy = 8 * s
    r = 6.5 * s
    return f"""
    <g transform="translate({x}, {y})" stroke="{OKABE_SKY_BLUE}" stroke-width="{1.2 * s}" stroke-linecap="round">
      <!-- Main axes -->
      <line x1="{cx}" y1="{cy - r}" x2="{cx}" y2="{cy + r}" />
      <line x1="{cx - r * 0.866}" y1="{cy - r * 0.5}" x2="{cx + r * 0.866}" y2="{cy + r * 0.5}" />
      <line x1="{cx - r * 0.866}" y1="{cy + r * 0.5}" x2="{cx + r * 0.866}" y2="{cy - r * 0.5}" />
      <!-- Chevrons -->
      <path d="M {cx - 2.5 * s} {cy - 4.5 * s} L {cx} {cy - 6 * s} L {cx + 2.5 * s} {cy - 4.5 * s}" fill="none" />
      <path d="M {cx - 2.5 * s} {cy + 4.5 * s} L {cx} {cy + 6 * s} L {cx + 2.5 * s} {cy + 4.5 * s}" fill="none" />
      <!-- Center dot -->
      <circle cx="{cx}" cy="{cy}" r="{1.2 * s}" fill="{OKABE_SKY_BLUE}" stroke="none" />
    </g>
    """


def render_funnel(x1: float, yt1: float, yb1: float, x2: float, yt2: float, yb2: float,
                  fill_color: str = "#F1F5F9", stroke_color: str = COLOR_BORDER_HAIRLINE) -> str:
    """Render a tapering or expanding 2D trapezoidal projection block."""
    points = f"{x1},{yt1} {x2},{yt2} {x2},{yb2} {x1},{yb1}"
    return f"""<polygon points="{points}" fill="{fill_color}" stroke="{stroke_color}" stroke-width="0.75" />"""


def render_svg() -> str:
    """Generate the complete publication-grade graphical abstract SVG."""
    svg_parts: list[str] = []

    # 1. XML Header and root with Inkscape namespaces
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="fig-graphical-abstract-section6"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <defs>
    <!-- Minimalist Hairline Arrow Markers -->
    <marker id="arr-slate" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
    </marker>
    <marker id="arr-blue" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
    <marker id="arr-muted" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_HAIRLINE}" />
    </marker>

    <!-- Linear Gradients for Heatmaps and Indicators -->
    <linearGradient id="grad-abstraction" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%" stop-color="#94A3B8" />
      <stop offset="50%" stop-color="{OKABE_SKY_BLUE}" />
      <stop offset="100%" stop-color="{OKABE_BLUE}" />
    </linearGradient>

    <linearGradient id="grad-heatmap-scale" x1="0%" y1="100%" x2="0%" y2="0%">
      <stop offset="0%" stop-color="#2563EB" />
      <stop offset="50%" stop-color="#FFFFFF" />
      <stop offset="100%" stop-color="#DC2626" />
    </linearGradient>
  </defs>
""")

    # ----------------------------------------------------
    # LAYER 00: Canvas Background & Panel Containers
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND & PANEL CONTAINERS -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />

    <!-- Panel A Container (Top-Left) -->
    <rect id="card-panel-a" x="20" y="62" width="430" height="386"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel B Container (Top-Right) -->
    <rect id="card-panel-b" x="465" y="62" width="915" height="386"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel D Container (Bottom-Left) -->
    <rect id="card-panel-d" x="20" y="460" width="670" height="452"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel C Container (Bottom-Middle) -->
    <rect id="card-panel-c" x="705" y="460" width="355" height="452"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel E Container (Bottom-Right) -->
    <rect id="card-panel-e" x="1075" y="460" width="305" height="452"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 01: Figure Header & Panel Badges
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 01: FIGURE HEADER & PANEL BADGES -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <!-- Main Title (Strict keywords / scientific tokens) -->
    <text id="txt-title" x="20" y="28" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">
      Multi-Cohort Immunotherapy Benchmarking &amp; Regularized Simplex Deconvolution Framework
    </text>
    <!-- Subtitle Summary -->
    <text id="txt-subtitle" x="20" y="46" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      9 Primary Trials (n = 1,097), 32 Predictors, Standardized Meta-Pools, Survival Discrimination &amp; Decision Curve Analysis
    </text>

    <!-- Top Total Stat Pill -->
    <g id="grp-header-stat" transform="translate(1040, 14)">
      <rect x="0" y="0" width="340" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="170" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        Compendium: 9 Clinical Trials | 1,097 Patients | 5 Solid Tumors
      </text>
    </g>

    <!-- Panel Badges: a, b, c, d, e -->
    <g id="badge-panel-a" transform="translate(34, 82)">
      <text x="0" y="12" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">a</text>
      <text x="16" y="12" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">TME BIOMARKER ARCHITECTURE &amp; SIMPLEX PROJECTION</text>
    </g>

    <g id="badge-panel-b" transform="translate(479, 82)">
      <text x="0" y="12" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">b</text>
      <text x="16" y="12" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">MULTI-COHORT CLINICAL TRIAL ECOSYSTEM (N = 1,097, 9 TRIALS)</text>
    </g>

    <g id="badge-panel-d" transform="translate(34, 480)">
      <text x="0" y="12" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">d</text>
      <text x="16" y="12" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">DISCRIMINATION METRICS, SURVIVAL &amp; DECISION UTILITY</text>
    </g>

    <g id="badge-panel-c" transform="translate(719, 480)">
      <text x="0" y="12" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">c</text>
      <text x="16" y="12" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">RESPONSE CLASSIFICATION &amp; EXPLAINABILITY</text>
    </g>

    <g id="badge-panel-e" transform="translate(1089, 480)">
      <text x="0" y="12" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">e</text>
      <text x="16" y="12" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">TRANSLATIONAL PREDICTOR HIERARCHY</text>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 02: PANEL A — TME Architecture & Simplex Projection
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: PANEL A — TME ARCHITECTURE & SIMPLEX PROJECTION -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_TME_Architecture">
    <!-- Top Sub-Block: TME Microenvironment Circle & Lineages -->
    <g transform="translate(34, 102)">
      <!-- Outer Card -->
      <rect width="402" height="116" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Microenvironment Circle Schematic -->
      <g transform="translate(12, 10)">
        <circle cx="48" cy="48" r="44" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.2" stroke-dasharray="3,2" />
        <text x="48" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">TIME Milieu</text>
        
        <!-- Cell cluster marks -->
        <!-- CD8 T cells (blue) -->
        <circle cx="34" cy="38" r="7" fill="{OKABE_BLUE}" opacity="0.9" />
        <circle cx="46" cy="42" r="5.5" fill="{OKABE_BLUE}" opacity="0.9" />
        <!-- TAMs (orange) -->
        <circle cx="62" cy="34" r="8" fill="{OKABE_ORANGE}" opacity="0.9" />
        <!-- CAFs (purple) -->
        <ellipse cx="36" cy="62" rx="9" ry="4.5" fill="{OKABE_REDDISH_PURPLE}" opacity="0.85" transform="rotate(-20 36 62)" />
        <!-- Endothelial (sky blue) -->
        <rect x="52" y="60" width="16" height="5" rx="2" fill="{OKABE_SKY_BLUE}" opacity="0.85" />
        <!-- Malignant parenchyma (vermilion) -->
        <circle cx="50" cy="50" r="4" fill="{OKABE_VERMILION}" />
        <circle cx="26" cy="48" r="4" fill="{OKABE_VERMILION}" />
        <circle cx="70" cy="54" r="4.5" fill="{OKABE_VERMILION}" />
        
        <text x="48" y="86" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Cellular Niches</text>
      </g>

      <!-- Cell Compartment Taxonomy List -->
      <g transform="translate(125, 14)">
        <text x="0" y="10" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Evaluated Compartments (32 Predictors)</text>
        
        <!-- Column 1 -->
        <g transform="translate(0, 22)">
          <circle cx="4" cy="4" r="3.5" fill="{OKABE_BLUE}" />
          <text x="12" y="7" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CD8+ Cytotoxic T Cells</text>
          
          <circle cx="4" cy="18" r="3.5" fill="{OKABE_ORANGE}" />
          <text x="12" y="21" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">M2 Macrophages / TAMs</text>
          
          <circle cx="4" cy="32" r="3.5" fill="{OKABE_REDDISH_PURPLE}" />
          <text x="12" y="35" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CAFs &amp; Stroma (Nurmik)</text>

          <circle cx="4" cy="46" r="3.5" fill="{OKABE_SKY_BLUE}" />
          <text x="12" y="49" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Endothelial Cells</text>

          <circle cx="4" cy="60" r="3.5" fill="{OKABE_VERMILION}" />
          <text x="12" y="63" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Tumor Cells &amp; Malignant</text>
        </g>

        <!-- Column 2: Molecular Features -->
        <g transform="translate(145, 22)">
          <rect x="0" y="1" width="6" height="6" fill="{COLOR_TEXT_MUTED}" />
          <text x="10" y="7" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">GEP, CYT, Teff, TIDE</text>

          <rect x="0" y="15" width="6" height="6" fill="{COLOR_TEXT_MUTED}" />
          <text x="10" y="21" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CXCL9, CD8A, PDCD1</text>

          <rect x="0" y="29" width="6" height="6" fill="{COLOR_TEXT_MUTED}" />
          <text x="10" y="35" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CD274 (PD-L1), CTLA4</text>

          <rect x="0" y="43" width="6" height="6" fill="{COLOR_TEXT_MUTED}" />
          <text x="10" y="49" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TMB Genomic Burden</text>

          <rect x="0" y="57" width="6" height="6" fill="{COLOR_TEXT_MUTED}" />
          <text x="10" y="63" font-size="7.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">AP-Gated Composites</text>
        </g>
      </g>
    </g>

    <!-- Middle Sub-Block: Simplex Deconvolution Projection Funnel -->
    <g transform="translate(34, 228)">
      <rect width="402" height="114" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Regularized Simplex Deconvolution (RegDeconv)</text>

      <!-- Left Input Block: Bulk RNA-seq -->
      <g transform="translate(10, 26)">
        <rect width="78" height="46" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="39" y="16" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Bulk RNA-seq</text>
        <text x="39" y="28" font-size="6.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">x ∈ ℝᴳ</text>
        <text x="39" y="39" font-size="6.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">G = 18,500 Genes</text>
      </g>

      <!-- Tapering Funnel -->
      {render_funnel(90, 32, 66, 126, 40, 58, "#F1F5F9", COLOR_BORDER_HAIRLINE)}
      <text x="108" y="37" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Project</text>

      <!-- Model Core Box: RegDeconv FISTA -->
      <g transform="translate(128, 25)">
        <rect width="90" height="48" fill="#FFFFFF" stroke="{OKABE_BLUE}" stroke-width="1.2" />
        <text x="45" y="16" font-size="8" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">RegDeconv</text>
        <text x="45" y="28" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Simplex FISTA</text>
        <text x="45" y="40" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">κ_target = 2,000</text>
      </g>

      <!-- Expanding Funnel -->
      {render_funnel(220, 40, 58, 256, 32, 66, "#F1F5F9", COLOR_BORDER_HAIRLINE)}
      <text x="238" y="37" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Fractions</text>

      <!-- Right Output Block: Simplex Proportions -->
      <g transform="translate(258, 26)">
        <rect width="134" height="46" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="67" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Simplex Proportions</text>
        <text x="67" y="26" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">θ ∈ Δ^(S-1), ∑ θ_s = 1</text>
        <text x="67" y="38" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Collinearity Deficit Nullified</text>
      </g>

      <!-- Mathematical Bullet Invariants -->
      <g transform="translate(10, 80)">
        <text x="0" y="10" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">• Graph Laplacian Energy:</text>
        <text x="96" y="10" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">λ_lap · θᵀLθ</text>
        
        <text x="195" y="10" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">• Poisson KL Deviance:</text>
        <text x="290" y="10" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">D_KL(x, Φθ)</text>

        <text x="0" y="23" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">• Projected Simplex:</text>
        <text x="86" y="23" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">Π_Δ(θ - α∇ℒ)</text>

        <text x="195" y="23" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">• Fused Grouping:</text>
        <text x="268" y="23" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">λ_fuse · ψ_ε(θ_i - θ_j)</text>
      </g>
    </g>

    <!-- Bottom Sub-Block: Hierarchical Feature Learning Flow -->
    <g transform="translate(34, 350)">
      <rect width="402" height="88" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Feature Abstraction Hierarchy</text>

      <!-- Flow nodes: 4 steps -->
      <g transform="translate(10, 22)">
        <!-- Step 1: Single Genes -->
        <rect x="0" y="0" width="82" height="34" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="41" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Single Genes</text>
        <text x="41" y="26" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">CD8A, CXCL9</text>

        <!-- Arrow 1 -> 2 -->
        <path d="M 85 17 L 95 17" stroke="{COLOR_TEXT_MUTED}" stroke-width="1" marker-end="url(#arr-slate)" />

        <!-- Step 2: Signatures -->
        <rect x="98" y="0" width="86" height="34" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="141" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Signatures</text>
        <text x="141" y="26" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">GEP, Teff, CYT</text>

        <!-- Arrow 2 -> 3 -->
        <path d="M 187 17 L 197 17" stroke="{COLOR_TEXT_MUTED}" stroke-width="1" marker-end="url(#arr-slate)" />

        <!-- Step 3: Deconv Simplex -->
        <rect x="200" y="0" width="94" height="34" fill="#FFFFFF" stroke="{OKABE_BLUE}" stroke-width="1" />
        <text x="247" y="14" font-size="7" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Deconv Simplex</text>
        <text x="247" y="26" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">RegDeconv (Our)</text>

        <!-- Arrow 3 -> 4 -->
        <path d="M 297 17 L 307 17" stroke="{COLOR_TEXT_MUTED}" stroke-width="1" marker-end="url(#arr-slate)" />

        <!-- Step 4: Multi-Omics -->
        <rect x="310" y="0" width="72" height="34" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="346" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Multi-Omics</text>
        <text x="346" y="26" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">CXCL9 + TMB</text>
      </g>

      <!-- Gradient Abstraction Bar -->
      <g transform="translate(10, 64)">
        <text x="0" y="9" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_MUTED}">Low</text>
        <rect x="26" y="3" width="280" height="7" fill="url(#grad-abstraction)" />
        <text x="312" y="9" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">High</text>
        <text x="166" y="18" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Biological Feature Abstraction &amp; Lineage Grounding</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 03: PANEL B — Multi-Cohort Clinical Trial Ecosystem
    # ----------------------------------------------------
    patient_glyphs: list[str] = []
    # Grid: 4 rows x 5 columns
    colors = [
        COLOR_MELANOMA, COLOR_MELANOMA, COLOR_MELANOMA, COLOR_MELANOMA, COLOR_MELANOMA,
        COLOR_MELANOMA, COLOR_RCC, COLOR_RCC, COLOR_RCC, COLOR_RCC,
        COLOR_RCC, COLOR_BLADDER, COLOR_BLADDER, COLOR_BLADDER, COLOR_BLADDER,
        COLOR_BLADDER, COLOR_PDAC, COLOR_PDAC, COLOR_TNBC, COLOR_TNBC
    ]
    badges = [
        "R", "NR", "R", "NR", "R",
        "NR", "R", "NR", "NR", "R",
        "NR", "NR", "R", "NR", "NR",
        "R", "R", "NR", "NR", "R"
    ]
    for idx, (c, b) in enumerate(zip(colors, badges)):
        row = idx // 5
        col = idx % 5
        px = 16 + col * 44
        py = 32 + row * 34
        patient_glyphs.append(render_person_icon(px, py, c, badge=b, scale=1.15))

    patient_grid_markup = "\n".join(patient_glyphs)

    svg_parts.append(f"""
  <!-- LAYER 03: PANEL B — MULTI-COHORT CLINICAL TRIAL ECOSYSTEM -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Clinical_Ecosystem">
    
    <!-- Left Sub-Card: Clinical Trial Cohorts (n = 1,097) -->
    <g transform="translate(479, 102)">
      <rect width="250" height="236" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Primary Clinical Cohorts (n = 1,097)</text>
      <text x="12" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Tumor bulk RNA-seq across 9 Phase II/III trials</text>

      <!-- Patient Grid -->
      <g transform="translate(10, 6)">
        {patient_grid_markup}
      </g>

      <!-- Indication Color Legend Pills -->
      <g transform="translate(10, 182)">
        <line x1="0" y1="0" x2="230" y2="0" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        
        <g transform="translate(0, 10)">
          <circle cx="4" cy="4" r="3.5" fill="{COLOR_MELANOMA}" />
          <text x="11" y="7" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Melanoma (338)</text>

          <circle cx="88" cy="4" r="3.5" fill="{COLOR_RCC}" />
          <text x="95" y="7" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">ccRCC (263)</text>

          <circle cx="162" cy="4" r="3.5" fill="{COLOR_BLADDER}" />
          <text x="169" y="7" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Bladder (298)</text>
        </g>
        <g transform="translate(0, 26)">
          <circle cx="4" cy="4" r="3.5" fill="{COLOR_PDAC}" />
          <text x="11" y="7" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">PDAC (85)</text>

          <circle cx="88" cy="4" r="3.5" fill="{COLOR_TNBC}" />
          <text x="95" y="7" font-size="7" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TNBC (31)</text>

          <rect x="162" y="1" width="6" height="6" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="171" y="7" font-size="6.8" fill="{COLOR_TEXT_MUTED}">N = 1,097 total</text>
        </g>
      </g>
    </g>

    <!-- Harmonization / Standardization Funnel -->
    <g transform="translate(734, 160)">
      <!-- Tapering Funnel -->
      {render_funnel(0, 20, 80, 26, 32, 68, "#F1F5F9", COLOR_BORDER_HAIRLINE)}
      
      <!-- Standardization Normalizer Box -->
      <g transform="translate(28, 22)">
        <rect width="68" height="56" fill="#FFFFFF" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="1.2" />
        <text x="34" y="18" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Standardize</text>
        <text x="34" y="32" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">z = (s - s̄) / σ</text>
        <text x="34" y="46" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Within-Cohort</text>
      </g>

      <!-- Expanding Funnel -->
      {render_funnel(98, 32, 68, 124, 16, 84, "#F1F5F9", COLOR_BORDER_HAIRLINE)}
    </g>

    <!-- Right Sub-Card: Standardized Meta-Cohort Pools & Discrimination Matrix -->
    <g transform="translate(864, 102)">
      <rect width="502" height="236" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Standardized Meta-Cohort Pools &amp; Benchmark Matrix</text>

      <!-- 3 Pooled Meta-Cohort Containers (Wrapped cleanly in sub-groups) -->
      <g transform="translate(12, 26)">
        <!-- Pool 1: Melanoma Combined -->
        <g transform="translate(0, 0)">
          <rect width="150" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_MELANOMA}">Combined Melanoma</text>
          <text x="8" y="26" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">n = 338 | Responders: 131 (38.8%)</text>
          <text x="8" y="36" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Hugo, Riaz, Gide, Liu (anti-PD-1)</text>
        </g>

        <!-- Pool 2: RCC Combined -->
        <g transform="translate(160, 0)">
          <rect width="150" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_RCC}">Combined Renal Cell</text>
          <text x="8" y="26" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">n = 263 | Responders: 75 (28.5%)</text>
          <text x="8" y="36" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Choueiri, McDermott (anti-PD-(L)1)</text>
        </g>

        <!-- Pool 3: Pan-Cancer Combined -->
        <g transform="translate(320, 0)">
          <rect width="158" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="0.75" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">Combined Pan-Cancer</text>
          <text x="8" y="26" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">n = 1,015 | Responders: 319 (31.4%)</text>
          <text x="8" y="36" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Pre-Treatment Standard Meta-Pool</text>
        </g>
      </g>

      <!-- Compact Benchmark Matrix Heatmap (like COMPASS Embeddings) -->
      <g transform="translate(12, 78)">
        <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Discrimination Performance Across Cohorts (ROC-AUC)</text>
        
        <!-- Heatmap Grid -->
        <g transform="translate(94, 16)">
          <!-- Column Labels -->
          <text x="14" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Cho</text>
          <text x="40" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">McD</text>
          <text x="66" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Hug</text>
          <text x="92" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Ria</text>
          <text x="118" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Gid</text>
          <text x="144" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Liu</text>
          <text x="170" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Pad</text>
          <text x="196" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">And</text>
          <text x="222" y="0" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Ros</text>
          <text x="254" y="0" font-size="6" font-weight="700" fill="{COLOR_MELANOMA}" text-anchor="middle">cMel</text>
          <text x="284" y="0" font-size="6" font-weight="700" fill="{COLOR_RCC}" text-anchor="middle">cRCC</text>
          <text x="320" y="0" font-size="6.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Pan-Can</text>

          <!-- Row 1: RegDeconv (Our) -->
          <g transform="translate(0, 4)">
            <rect x="2" y="0" width="24" height="12" fill="#E02424" opacity="0.85" /><text x="14" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.78</text>
            <rect x="28" y="0" width="24" height="12" fill="#EF4444" opacity="0.85" /><text x="40" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.71</text>
            <rect x="54" y="0" width="24" height="12" fill="#DC2626" opacity="0.85" /><text x="66" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.76</text>
            <rect x="80" y="0" width="24" height="12" fill="#EF4444" opacity="0.85" /><text x="92" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.69</text>
            <rect x="106" y="0" width="24" height="12" fill="#DC2626" opacity="0.85" /><text x="118" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.74</text>
            <rect x="132" y="0" width="24" height="12" fill="#EF4444" opacity="0.85" /><text x="144" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.70</text>
            <rect x="158" y="0" width="24" height="12" fill="#F87171" opacity="0.85" /><text x="170" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.67</text>
            <rect x="184" y="0" width="24" height="12" fill="#EF4444" opacity="0.85" /><text x="196" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.72</text>
            <rect x="210" y="0" width="24" height="12" fill="#EF4444" opacity="0.85" /><text x="222" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.70</text>
            <!-- Combined -->
            <rect x="242" y="0" width="24" height="12" fill="#DC2626" opacity="0.95" /><text x="254" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.73</text>
            <rect x="272" y="0" width="24" height="12" fill="#DC2626" opacity="0.95" /><text x="284" y="9" font-size="5.8" fill="#FFF" font-weight="700" text-anchor="middle">.71</text>
            <rect x="302" y="0" width="36" height="12" fill="#B91C1C" /><text x="320" y="9" font-size="6.2" fill="#FFF" font-weight="700" text-anchor="middle">.720</text>
          </g>

          <!-- Row 2: Teff Signature -->
          <g transform="translate(0, 18)">
            <rect x="2" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="14" y="9" font-size="5.8" fill="#000" text-anchor="middle">.62</text>
            <rect x="28" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="40" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.66</text>
            <rect x="54" y="0" width="24" height="12" fill="#EF4444" opacity="0.8" /><text x="66" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.69</text>
            <rect x="80" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="92" y="9" font-size="5.8" fill="#000" text-anchor="middle">.61</text>
            <rect x="106" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="118" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.67</text>
            <rect x="132" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="144" y="9" font-size="5.8" fill="#000" text-anchor="middle">.64</text>
            <rect x="158" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="170" y="9" font-size="5.8" fill="#000" text-anchor="middle">.54</text>
            <rect x="184" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="196" y="9" font-size="5.8" fill="#000" text-anchor="middle">.63</text>
            <rect x="210" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="222" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.65</text>
            <!-- Combined -->
            <rect x="242" y="0" width="24" height="12" fill="#F87171" /><text x="254" y="9" font-size="5.8" fill="#FFF" font-weight="600" text-anchor="middle">.66</text>
            <rect x="272" y="0" width="24" height="12" fill="#F87171" /><text x="284" y="9" font-size="5.8" fill="#FFF" font-weight="600" text-anchor="middle">.65</text>
            <rect x="302" y="0" width="36" height="12" fill="#EF4444" /><text x="320" y="9" font-size="6.2" fill="#FFF" font-weight="700" text-anchor="middle">.654</text>
          </g>

          <!-- Row 3: CD8 Duo -->
          <g transform="translate(0, 32)">
            <rect x="2" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="14" y="9" font-size="5.8" fill="#000" text-anchor="middle">.61</text>
            <rect x="28" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="40" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.65</text>
            <rect x="54" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="66" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.68</text>
            <rect x="80" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="92" y="9" font-size="5.8" fill="#000" text-anchor="middle">.59</text>
            <rect x="106" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="118" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.66</text>
            <rect x="132" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="144" y="9" font-size="5.8" fill="#000" text-anchor="middle">.63</text>
            <rect x="158" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="170" y="9" font-size="5.8" fill="#000" text-anchor="middle">.52</text>
            <rect x="184" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="196" y="9" font-size="5.8" fill="#000" text-anchor="middle">.61</text>
            <rect x="210" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="222" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.64</text>
            <!-- Combined -->
            <rect x="242" y="0" width="24" height="12" fill="#F87171" /><text x="254" y="9" font-size="5.8" fill="#FFF" font-weight="600" text-anchor="middle">.65</text>
            <rect x="272" y="0" width="24" height="12" fill="#FCA5A5" /><text x="284" y="9" font-size="5.8" fill="#000" font-weight="600" text-anchor="middle">.64</text>
            <rect x="302" y="0" width="36" height="12" fill="#F87171" /><text x="320" y="9" font-size="6.2" fill="#FFF" font-weight="700" text-anchor="middle">.648</text>
          </g>

          <!-- Row 4: CXCL9 Single Gene -->
          <g transform="translate(0, 46)">
            <rect x="2" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="14" y="9" font-size="5.8" fill="#000" text-anchor="middle">.60</text>
            <rect x="28" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="40" y="9" font-size="5.8" fill="#000" text-anchor="middle">.62</text>
            <rect x="54" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="66" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.67</text>
            <rect x="80" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="92" y="9" font-size="5.8" fill="#000" text-anchor="middle">.58</text>
            <rect x="106" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="118" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.65</text>
            <rect x="132" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="144" y="9" font-size="5.8" fill="#000" text-anchor="middle">.61</text>
            <rect x="158" y="0" width="24" height="12" fill="#E2E8F0" opacity="0.8" /><text x="170" y="9" font-size="5.8" fill="#000" text-anchor="middle">.48</text>
            <rect x="184" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="196" y="9" font-size="5.8" fill="#000" text-anchor="middle">.60</text>
            <rect x="210" y="0" width="24" height="12" fill="#F87171" opacity="0.8" /><text x="222" y="9" font-size="5.8" fill="#FFF" text-anchor="middle">.63</text>
            <!-- Combined -->
            <rect x="242" y="0" width="24" height="12" fill="#FCA5A5" /><text x="254" y="9" font-size="5.8" fill="#000" font-weight="600" text-anchor="middle">.64</text>
            <rect x="272" y="0" width="24" height="12" fill="#FCA5A5" /><text x="284" y="9" font-size="5.8" fill="#000" font-weight="600" text-anchor="middle">.62</text>
            <rect x="302" y="0" width="36" height="12" fill="#FCA5A5" /><text x="320" y="9" font-size="6.2" fill="#000" font-weight="700" text-anchor="middle">.633</text>
          </g>

          <!-- Row 5: TMB Genomic -->
          <g transform="translate(0, 60)">
            <rect x="2" y="0" width="24" height="12" fill="#E2E8F0" opacity="0.8" /><text x="14" y="9" font-size="5.8" fill="#000" text-anchor="middle">.50</text>
            <rect x="28" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="40" y="9" font-size="5.8" fill="#000" text-anchor="middle">.56</text>
            <rect x="54" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="66" y="9" font-size="5.8" fill="#000" text-anchor="middle">.61</text>
            <rect x="80" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="92" y="9" font-size="5.8" fill="#000" text-anchor="middle">.58</text>
            <rect x="106" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="118" y="9" font-size="5.8" fill="#000" text-anchor="middle">.60</text>
            <rect x="132" y="0" width="24" height="12" fill="#FEE2E2" opacity="0.8" /><text x="144" y="9" font-size="5.8" fill="#000" text-anchor="middle">.57</text>
            <rect x="158" y="0" width="24" height="12" fill="#E2E8F0" opacity="0.8" /><text x="170" y="9" font-size="5.8" fill="#000" text-anchor="middle">.50</text>
            <rect x="184" y="0" width="24" height="12" fill="#E2E8F0" opacity="0.8" /><text x="196" y="9" font-size="5.8" fill="#000" text-anchor="middle">.51</text>
            <rect x="210" y="0" width="24" height="12" fill="#FCA5A5" opacity="0.8" /><text x="222" y="9" font-size="5.8" fill="#000" text-anchor="middle">.60</text>
            <!-- Combined -->
            <rect x="242" y="0" width="24" height="12" fill="#FEE2E2" /><text x="254" y="9" font-size="5.8" fill="#000" font-weight="600" text-anchor="middle">.59</text>
            <rect x="272" y="0" width="24" height="12" fill="#FEE2E2" /><text x="284" y="9" font-size="5.8" fill="#000" font-weight="600" text-anchor="middle">.56</text>
            <rect x="302" y="0" width="36" height="12" fill="#FEE2E2" /><text x="320" y="9" font-size="6.2" fill="#000" font-weight="700" text-anchor="middle">.583</text>
          </g>
        </g>

        <!-- Row Labels -->
        <g transform="translate(0, 20)">
          <text x="88" y="9" font-size="6.8" font-weight="700" fill="{OKABE_BLUE}" text-anchor="end">RegDeconv (Our)</text>
          <text x="88" y="23" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="end">Teff Signature</text>
          <text x="88" y="37" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="end">CD8 Duo</text>
          <text x="88" y="51" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="end">CXCL9 Marker</text>
          <text x="88" y="65" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="end">TMB (Genomic)</text>
        </g>

        <!-- Colorbar Scale on Right -->
        <g transform="translate(448, 18)">
          <rect x="0" y="0" width="10" height="74" fill="url(#grad-heatmap-scale)" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="14" y="6" font-size="6" fill="{COLOR_TEXT_MUTED}">0.75</text>
          <text x="14" y="39" font-size="6" fill="{COLOR_TEXT_MUTED}">0.50</text>
          <text x="14" y="74" font-size="6" fill="{COLOR_TEXT_MUTED}">0.35</text>
          <text x="5" y="86" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">AUC</text>
        </g>
      </g>
    </g>

    <!-- Bottom Protocol Summary Box (matching Step 1 in COMPASS) -->
    <g transform="translate(479, 350)">
      <rect width="887" height="88" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Left step pill -->
      <g transform="translate(12, 12)">
        <rect width="170" height="20" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="85" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Step 1: Multi-Trial Harmonization</text>
      </g>

      <g transform="translate(195, 16)">
        <text x="0" y="10" font-size="7.8" font-weight="700" fill="{OKABE_BLUE}">↑ Broad Cohort Scale (N = 1,097)</text>
        <text x="150" y="10" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">✔ 5 Indications (Melanoma, ccRCC, Bladder, PDAC, TNBC)</text>
        <text x="525" y="10" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">✔ Anti-PD-(L)1 &amp; Anti-CTLA-4</text>
      </g>

      <line x1="12" y1="42" x2="875" y2="42" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <g transform="translate(12, 54)">
        <text x="0" y="10" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Standardization Protocol:</text>
        <text x="125" y="10" font-size="7.2" fill="{COLOR_TEXT_SECONDARY}">Within-cohort z-score transformation eliminates study baseline shifts</text>
        <text x="495" y="10" font-size="7.2" font-weight="600" fill="{OKABE_BLUISH_GREEN}">z_ic = (s_ic - s̄_c) / σ_c</text>
        <text x="625" y="10" font-size="7.2" fill="{COLOR_TEXT_MUTED}">| Dual strata evaluation (Standard vs Extreme)</text>
      </g>
      <g transform="translate(12, 70)">
        <text x="0" y="10" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Meta-Cohort Aggregation:</text>
        <text x="125" y="10" font-size="7.2" fill="{COLOR_TEXT_SECONDARY}">Combined Melanoma (n = 338) | Combined Renal Cell (n = 263) | Combined Pan-Cancer (n = 1,015)</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 04: PANEL D — Discrimination Metrics, Survival & Decision Utility
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: PANEL D — DISCRIMINATION METRICS, SURVIVAL & DECISION UTILITY -->
  <g inkscape:groupmode="layer" id="layer-04-panel-d" inkscape:label="04_Panel_D_Discrimination_Survival">
    
    <!-- Sub-Card D1: Discrimination Engines (ROC & PR AUC) -->
    <g transform="translate(34, 502)">
      <rect width="210" height="396" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Discrimination Engines</text>
      <text x="10" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Pan-cancer pre-standard (n = 1,015)</text>

      <!-- Mini ROC Curve Plot -->
      <g transform="translate(18, 38)">
        <rect width="174" height="122" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">ROC Discrimination</text>
        
        <!-- Axes -->
        <line x1="28" y1="102" x2="162" y2="102" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="28" y1="18" x2="28" y2="102" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <!-- Ticks -->
        <line x1="28" y1="102" x2="28" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="95" y1="102" x2="95" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="162" y1="102" x2="162" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <text x="28" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
        <text x="95" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.5</text>
        <text x="162" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">1.0</text>
        <text x="95" y="120" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">FPR (1 - Specificity)</text>

        <text x="22" y="104" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0</text>
        <text x="22" y="62" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0.5</text>
        <text x="22" y="22" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">1.0</text>
        <text x="12" y="60" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" transform="rotate(-90 12 60)" text-anchor="middle">TPR (Sensitivity)</text>

        <!-- Diagonal Null Line -->
        <line x1="28" y1="102" x2="162" y2="18" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />

        <!-- Hanley-McNeil CI Shading Band for RegDeconv -->
        <path d="M 28 102 C 34 50, 60 28, 162 18 L 162 18 C 80 40, 52 70, 28 102 Z" fill="{OKABE_BLUE}" opacity="0.18" />

        <!-- RegDeconv ROC Curve (AUC = 0.7201) -->
        <path d="M 28 102 C 40 52, 70 32, 162 18" stroke="{OKABE_BLUE}" stroke-width="1.8" fill="none" />

        <!-- Teff ROC Curve (AUC = 0.654) -->
        <path d="M 28 102 C 55 68, 92 46, 162 18" stroke="{OKABE_ORANGE}" stroke-width="1.1" fill="none" stroke-dasharray="3,1" />

        <!-- Legend inside plot -->
        <g transform="translate(74, 66)">
          <line x1="0" y1="4" x2="12" y2="4" stroke="{OKABE_BLUE}" stroke-width="1.8" />
          <text x="16" y="7" font-size="6" font-weight="700" fill="{OKABE_BLUE}">RegDeconv (.720)</text>

          <line x1="0" y1="14" x2="12" y2="14" stroke="{OKABE_ORANGE}" stroke-width="1.1" stroke-dasharray="3,1" />
          <text x="16" y="17" font-size="6" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Teff (.654)</text>

          <line x1="0" y1="24" x2="12" y2="24" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />
          <text x="16" y="27" font-size="6" fill="{COLOR_TEXT_MUTED}">Baseline (.500)</text>
        </g>
      </g>

      <!-- Hanley-McNeil Exact Variance Formula Box -->
      <g transform="translate(18, 172)">
        <rect width="174" height="42" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Hanley-McNeil Analytical Variance</text>
        <text x="8" y="24" font-size="6.2" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">Var(A) = [A(1-A) + (n_R-1)(Q₁-A²) + (n_NR-1)(Q₂-A²)] / (n_R·n_NR)</text>
        <text x="8" y="35" font-size="6" fill="{COLOR_TEXT_MUTED}">95% CI: [0.687, 0.753] | Finite-sample exact</text>
      </g>

      <!-- Mini PR Curve Plot -->
      <g transform="translate(18, 226)">
        <rect width="174" height="122" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">PR Discrimination (Excess Precision)</text>

        <!-- Axes -->
        <line x1="28" y1="102" x2="162" y2="102" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="28" y1="18" x2="28" y2="102" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <!-- Ticks -->
        <line x1="28" y1="102" x2="28" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="95" y1="102" x2="95" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="162" y1="102" x2="162" y2="106" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <text x="28" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
        <text x="95" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.5</text>
        <text x="162" y="114" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">1.0</text>
        <text x="95" y="120" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Recall (Sensitivity)</text>

        <text x="22" y="104" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0</text>
        <text x="22" y="60" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0.5</text>
        <text x="22" y="22" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">1.0</text>
        <text x="12" y="60" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" transform="rotate(-90 12 60)" text-anchor="middle">Precision (PPV)</text>

        <!-- Prevalence Horizontal Line pi = 0.314 -->
        <line x1="28" y1="76" x2="162" y2="76" stroke="{COLOR_TEXT_MUTED}" stroke-width="0.75" stroke-dasharray="2,2" />
        <text x="160" y="73" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="end">π = 0.314</text>

        <!-- RegDeconv PR Curve (PR-AUC = 0.495) -->
        <path d="M 28 26 C 65 30, 105 48, 162 76" stroke="{OKABE_BLUE}" stroke-width="1.8" fill="none" />

        <!-- Teff PR Curve (PR-AUC = 0.412) -->
        <path d="M 28 38 C 75 48, 115 62, 162 76" stroke="{OKABE_ORANGE}" stroke-width="1.1" fill="none" stroke-dasharray="3,1" />

        <!-- Legend inside PR plot (cleanly positioned) -->
        <g transform="translate(66, 18)">
          <rect width="94" height="28" fill="#FFFFFF" opacity="0.9" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="10" font-size="5.8" font-weight="700" fill="{OKABE_BLUE}">RegDeconv (.495)</text>
          <text x="6" y="20" font-size="5.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}">ΔPR: +0.181</text>
          <text x="56" y="20" font-size="5.5" fill="{COLOR_TEXT_MUTED}">Teff (.412)</text>
        </g>
      </g>

      <!-- Bottom Tag -->
      <g transform="translate(18, 360)">
        <rect width="174" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="87" y="15" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Dual Stratum Invariance: Standard &amp; Extreme
        </text>
      </g>
    </g>

    <!-- Sub-Card D2: Time-to-Event Survival Analysis -->
    <g transform="translate(254, 502)">
      <rect width="200" height="396" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Time-to-Event Survival</text>
      <text x="10" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Harrell's C-index &amp; Cox Proportional Hazards</text>

      <!-- Harrell's C-Index Forest Chart -->
      <g transform="translate(12, 38)">
        <rect width="176" height="174" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Overall Survival Discrimination (C-Index)</text>

        <!-- Vertical Reference Line at C = 0.50 -->
        <line x1="88" y1="24" x2="88" y2="150" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />
        <text x="88" y="162" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.50 (Null)</text>
        <text x="144" y="162" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.65</text>
        <line x1="144" y1="24" x2="144" y2="150" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

        <!-- Forest Items -->
        <!-- 1. RegDeconv (Our) -->
        <g transform="translate(0, 32)">
          <text x="8" y="9" font-size="6.8" font-weight="700" fill="{OKABE_BLUE}">RegDeconv (Our)</text>
          <line x1="120" y1="6" x2="158" y2="6" stroke="{OKABE_BLUE}" stroke-width="1.2" />
          <circle cx="139" cy="6" r="3.2" fill="{OKABE_BLUE}" />
          <text x="164" y="9" font-size="6.2" font-weight="700" fill="{OKABE_BLUE}">.635</text>
        </g>

        <!-- 2. Teff Signature -->
        <g transform="translate(0, 56)">
          <text x="8" y="9" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Teff Signature</text>
          <line x1="104" y1="6" x2="134" y2="6" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" />
          <circle cx="119" cy="6" r="2.8" fill="{COLOR_TEXT_SECONDARY}" />
          <text x="164" y="9" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">.582</text>
        </g>

        <!-- 3. CD8 Duo -->
        <g transform="translate(0, 80)">
          <text x="8" y="9" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CD8 Duo (Chen)</text>
          <line x1="101" y1="6" x2="131" y2="6" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" />
          <circle cx="116" cy="6" r="2.8" fill="{COLOR_TEXT_SECONDARY}" />
          <text x="164" y="9" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">.575</text>
        </g>

        <!-- 4. CXCL9 Marker -->
        <g transform="translate(0, 104)">
          <text x="8" y="9" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">CXCL9 Marker</text>
          <line x1="97" y1="6" x2="125" y2="6" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" />
          <circle cx="111" cy="6" r="2.8" fill="{COLOR_TEXT_SECONDARY}" />
          <text x="164" y="9" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">.562</text>
        </g>

        <!-- 5. TMB (Genomic) -->
        <g transform="translate(0, 128)">
          <text x="8" y="9" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TMB Genomic</text>
          <line x1="91" y1="6" x2="119" y2="6" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" />
          <circle cx="105" cy="6" r="2.8" fill="{COLOR_TEXT_SECONDARY}" />
          <text x="164" y="9" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">.541</text>
        </g>
      </g>

      <!-- Univariable Cox Proportional Hazards Model Box -->
      <g transform="translate(12, 226)">
        <rect width="176" height="158" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cox Proportional Hazards Model</text>
        <text x="8" y="26" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">λ(t | z) = λ₀(t) exp(β · z)</text>
        
        <!-- Hazard Ratio Forest -->
        <g transform="translate(8, 38)">
          <rect width="160" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="6.8" font-weight="700" fill="{OKABE_BLUE}">RegDeconv Proportions</text>
          <text x="8" y="26" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">HR = 0.725 [0.63, 0.83]</text>
          <text x="8" y="36" font-size="6.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Wald p = 3.8 × 10⁻⁶</text>
        </g>

        <g transform="translate(8, 88)">
          <text x="0" y="10" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Clinical Risk Metric:</text>
          <text x="0" y="22" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• 27.5% Mortality Risk Reduction</text>
          <text x="0" y="34" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• Per 1-SD Proportions Increment</text>
          <text x="0" y="46" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• 1,000 Bootstrap Resamples</text>
        </g>
      </g>
    </g>

    <!-- Sub-Card D3: Decision Curve Analysis (DCA Net Benefit) -->
    <g transform="translate(464, 502)">
      <rect width="216" height="396" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. Decision Curve Analysis (DCA)</text>
      <text x="10" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Net Clinical Benefit across decision thresholds</text>

      <!-- Mini DCA Curve Plot -->
      <g transform="translate(14, 38)">
        <rect width="188" height="190" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Net Clinical Benefit Curve</text>

        <!-- Net Benefit Formula -->
        <text x="8" y="26" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">NB(p_t) = TP/N - FP/N · [p_t / (1 - p_t)]</text>

        <!-- Axes -->
        <line x1="28" y1="156" x2="174" y2="156" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="28" y1="36" x2="28" y2="156" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <!-- Ticks -->
        <line x1="28" y1="156" x2="28" y2="160" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="76" y1="156" x2="76" y2="160" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="125" y1="156" x2="125" y2="160" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="174" y1="156" x2="174" y2="160" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <text x="28" y="168" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
        <text x="76" y="168" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.25</text>
        <text x="125" y="168" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.50</text>
        <text x="174" y="168" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0.75</text>
        <text x="101" y="176" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Threshold Probability (p_t)</text>

        <text x="22" y="158" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0</text>
        <text x="22" y="98" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0.15</text>
        <text x="22" y="42" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="end">0.30</text>
        <text x="12" y="98" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" transform="rotate(-90 12 98)" text-anchor="middle">Net Benefit</text>

        <!-- Treat None Horizontal Line (y = 0) -->
        <line x1="28" y1="156" x2="174" y2="156" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="1" stroke-dasharray="2,2" />

        <!-- Treat All Patients Curve -->
        <path d="M 28 66 C 50 82, 85 120, 115 156" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.1" fill="none" />

        <!-- RegDeconv Net Benefit Curve -->
        <path d="M 28 66 C 52 74, 90 94, 155 156" stroke="{OKABE_BLUE}" stroke-width="2" fill="none" />

        <!-- Optimal Window Shading (0.10 to 0.50) -->
        <rect x="47" y="36" width="78" height="120" fill="{OKABE_BLUE}" opacity="0.08" />
        <text x="86" y="48" font-size="5.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Optimal Window</text>

        <!-- Legend inside DCA -->
        <g transform="translate(100, 72)">
          <line x1="0" y1="4" x2="12" y2="4" stroke="{OKABE_BLUE}" stroke-width="2" />
          <text x="16" y="7" font-size="5.8" font-weight="700" fill="{OKABE_BLUE}">RegDeconv (Our)</text>

          <line x1="0" y1="14" x2="12" y2="14" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.1" />
          <text x="16" y="17" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">Treat All</text>

          <line x1="0" y1="24" x2="12" y2="24" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="1" stroke-dasharray="2,2" />
          <text x="16" y="27" font-size="5.8" fill="{COLOR_TEXT_MUTED}">Treat None</text>
        </g>
      </g>

      <!-- Clinical Decision Summary Card -->
      <g transform="translate(14, 238)">
        <rect width="188" height="146" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Clinical Utility Summary</text>
        
        <g transform="translate(10, 26)">
          <rect width="168" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="0.75" />
          <text x="8" y="14" font-size="6.8" font-weight="700" fill="{OKABE_BLUE}">Net Benefit Advantage: +0.082</text>
          <text x="8" y="26" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Avoids 28 false biopsies / 100 patients</text>
        </g>

        <g transform="translate(10, 70)">
          <text x="0" y="10" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Actionable Decision Thresholds:</text>
          <text x="0" y="23" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• p_t ∈ [0.10, 0.50] primary window</text>
          <text x="0" y="35" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• Zero loss compared to default Treat-All</text>
          <text x="0" y="47" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• Preserves true responder identification</text>
          <text x="0" y="59" font-size="6.2" fill="{COLOR_TEXT_MUTED}">• Hanley-McNeil calibrated</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 05: PANEL C — Response Classification & Feature Attribution
    # ----------------------------------------------------
    c_patients: list[str] = []
    # 2 rows x 4 cols
    c_badges = ["R", "R", "NR", "NR", "R", "NR", "NR", "R"]
    c_colors = [COLOR_MELANOMA, COLOR_RCC, COLOR_BLADDER, COLOR_PDAC, COLOR_TNBC, COLOR_MELANOMA, COLOR_RCC, COLOR_BLADDER]
    for idx, (cb, cc) in enumerate(zip(c_badges, c_colors)):
        row = idx // 4
        col = idx % 4
        c_patients.append(render_person_icon(8 + col * 36, 12 + row * 28, cc, badge=cb, scale=1.0))
    c_patient_markup = "\n".join(c_patients)

    svg_parts.append(f"""
  <!-- LAYER 05: PANEL C — RESPONSE CLASSIFICATION & EXPLAINABILITY -->
  <g inkscape:groupmode="layer" id="layer-05-panel-c" inkscape:label="05_Panel_C_Response_Classification">
    
    <!-- Top Sub-Card: Clinical Response Strata & Patient Classification -->
    <g transform="translate(719, 502)">
      <rect width="327" height="114" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Clinical Response Strata (n = 1,015)</text>

      <!-- Patient pictogram box -->
      <g transform="translate(10, 24)">
        <rect width="154" height="78" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        {c_patient_markup}
      </g>

      <!-- Strata Counts & Legend -->
      <g transform="translate(172, 24)">
        <!-- Responder Pill -->
        <rect x="0" y="0" width="144" height="28" fill="#FFFFFF" stroke="{COLOR_RESPONDER}" stroke-width="0.75" />
        <circle cx="10" cy="14" r="4" fill="{COLOR_RESPONDER}" />
        <text x="20" y="12" font-size="7" font-weight="700" fill="{COLOR_RESPONDER}">Responders (CR / PR)</text>
        <text x="20" y="23" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">n = 319 (31.4%)</text>

        <!-- Non-Responder Pill -->
        <rect x="0" y="32" width="144" height="28" fill="#FFFFFF" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" />
        <circle cx="10" cy="46" r="4" fill="{COLOR_NON_RESPONDER}" />
        <text x="20" y="44" font-size="7" font-weight="700" fill="{COLOR_NON_RESPONDER}">Non-Responders (SD/PD)</text>
        <text x="20" y="55" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">n = 696 (68.6%)</text>

        <text x="0" y="74" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Dual Strata:</text>
        <text x="48" y="74" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">Standard (1,015) | Extreme (745)</text>
      </g>
    </g>

    <!-- Middle Sub-Card: Feature Attribution Mini-Bars (Clean non-overlapping layout) -->
    <g transform="translate(719, 626)">
      <rect width="327" height="156" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Personalized Response Maps (Feature Attribution)</text>
      <text x="12" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Relative biomarker attribution vectors (Shapley / Perturbation)</text>

      <!-- Clean Diverging Attribution Bar Chart -->
      <g transform="translate(10, 36)">
        <!-- Vertical Center Zero Line at x = 200 -->
        <line x1="200" y1="6" x2="200" y2="108" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <text x="200" y="4" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
        <text x="150" y="4" font-size="6" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">Negative (-)</text>
        <text x="250" y="4" font-size="6" fill="{COLOR_RESPONDER}" text-anchor="middle">Positive (+)</text>

        <!-- 1. CD8+ T Cells (+0.38) -->
        <g transform="translate(0, 10)">
          <text x="10" y="10" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">CD8+ Cytotoxic T Cells</text>
          <rect x="200" y="2" width="68" height="11" fill="{OKABE_BLUE}" />
          <text x="274" y="11" font-size="6.5" font-weight="700" fill="{OKABE_BLUE}">+0.38</text>
        </g>

        <!-- 2. CXCL9 Chemokine (+0.31) -->
        <g transform="translate(0, 30)">
          <text x="10" y="10" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">CXCL9 Chemokine Axis</text>
          <rect x="200" y="2" width="55" height="11" fill="{OKABE_SKY_BLUE}" />
          <text x="261" y="11" font-size="6.5" font-weight="700" fill="{OKABE_SKY_BLUE}">+0.31</text>
        </g>

        <!-- 3. M2 Macrophages (-0.26) -->
        <g transform="translate(0, 50)">
          <text x="10" y="10" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">M2 Macrophages (TAM)</text>
          <rect x="154" y="2" width="46" height="11" fill="{OKABE_ORANGE}" />
          <text x="148" y="11" font-size="6.5" font-weight="700" fill="{OKABE_ORANGE}" text-anchor="end">-0.26</text>
        </g>

        <!-- 4. CAFs / Stroma (-0.19) -->
        <g transform="translate(0, 70)">
          <text x="10" y="10" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">CAFs &amp; Dense Stroma</text>
          <rect x="166" y="2" width="34" height="11" fill="{OKABE_REDDISH_PURPLE}" />
          <text x="160" y="11" font-size="6.5" font-weight="700" fill="{OKABE_REDDISH_PURPLE}" text-anchor="end">-0.19</text>
        </g>

        <!-- 5. TMB Genomic (+0.14) -->
        <g transform="translate(0, 90)">
          <text x="10" y="10" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">TMB Genomic Burden</text>
          <rect x="200" y="2" width="25" height="11" fill="{OKABE_VERMILION}" />
          <text x="231" y="11" font-size="6.5" font-weight="700" fill="{OKABE_VERMILION}">+0.14</text>
        </g>
      </g>
    </g>

    <!-- Bottom Sub-Card: Translational Applications -->
    <g transform="translate(719, 792)">
      <rect width="327" height="106" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Translational Utility &amp; Clinical Impact</text>

      <g transform="translate(10, 24)">
        <rect width="96" height="70" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="48" y="16" font-size="7.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">1. Stratification</text>
        <text x="48" y="30" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Precision Medicine</text>
        <text x="48" y="44" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Enrich Responders</text>
        <text x="48" y="58" font-size="6" font-weight="600" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">ΔPR: +0.181</text>
      </g>

      <g transform="translate(114, 24)">
        <rect width="96" height="70" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="48" y="16" font-size="7.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">2. Optimization</text>
        <text x="48" y="30" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Trial Endpoints</text>
        <text x="48" y="44" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Phase II/III Design</text>
        <text x="48" y="58" font-size="6" font-weight="600" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">HR: 0.725</text>
      </g>

      <g transform="translate(218, 24)">
        <rect width="98" height="70" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="49" y="16" font-size="7.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">3. Combinations</text>
        <text x="49" y="30" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Target TAM / CAF</text>
        <text x="49" y="44" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Reverse Resistance</text>
        <text x="49" y="58" font-size="6" font-weight="600" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Deconv Signatures</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 06: PANEL E — Translational Predictor Hierarchy
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 06: PANEL E — TRANSLATIONAL PREDICTOR HIERARCHY -->
  <g inkscape:groupmode="layer" id="layer-06-panel-e" inkscape:label="06_Panel_E_Predictor_Hierarchy">
    
    <!-- Tier 1: RegDeconv Simplex (Rank 1) -->
    <g transform="translate(1089, 502)">
      <rect width="277" height="84" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="1.2" />
      
      <!-- Flame icon -->
      {render_flame_icon(10, 8, 16)}
      
      <text x="32" y="16" font-size="8" font-weight="700" fill="{OKABE_BLUE}">TIER 1 (RANK 1): RegDeconv Simplex (Our)</text>
      <text x="32" y="28" font-size="7" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Collinearity-Regularized Simplex FISTA</text>

      <g transform="translate(12, 38)">
        <rect width="76" height="18" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="38" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">ROC-AUC: 0.7201</text>

        <rect x="82" y="0" width="76" height="18" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="120" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">PR-AUC: 0.495</text>

        <rect x="164" y="0" width="88" height="18" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="208" y="12" font-size="6.8" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">C-Index: 0.635</text>
      </g>

      <text x="12" y="74" font-size="6.5" font-weight="700" fill="{OKABE_BLUE}">• Calibrated Simplex Proportions (Target κ = 2,000)</text>
    </g>

    <!-- Tier 2: Transcriptional Signatures (Rank 2-4) -->
    <g transform="translate(1089, 596)">
      <rect width="277" height="84" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Snowflake icon -->
      {render_snowflake_icon(10, 8, 16)}
      
      <text x="32" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">TIER 2 (RANK 2–4): Published Signatures</text>
      <text x="32" y="28" font-size="7" fill="{COLOR_TEXT_SECONDARY}">Teff (Fehrenbacher), CD8 Duo (Chen), GEP</text>

      <g transform="translate(12, 38)">
        <rect width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="38" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">ROC-AUC: 0.654</text>

        <rect x="82" y="0" width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="120" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">PR-AUC: 0.412</text>

        <rect x="164" y="0" width="88" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="208" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">C-Index: 0.582</text>
      </g>

      <text x="12" y="74" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Fixed Transcriptional Weights (No Cell Lineage Grounding)</text>
    </g>

    <!-- Tier 3: Single Checkpoint Markers (Rank 5-8) -->
    <g transform="translate(1089, 690)">
      <rect width="277" height="84" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Snowflake icon -->
      {render_snowflake_icon(10, 8, 16)}

      <text x="32" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">TIER 3 (RANK 5–8): Single Gene Markers</text>
      <text x="32" y="28" font-size="7" fill="{COLOR_TEXT_SECONDARY}">CXCL9, CD8A, PDCD1, CD274, CTLA4</text>

      <g transform="translate(12, 38)">
        <rect width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="38" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">ROC-AUC: 0.633</text>

        <rect x="82" y="0" width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="120" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">PR-AUC: 0.395</text>

        <rect x="164" y="0" width="88" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="208" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">C-Index: 0.562</text>
      </g>

      <text x="12" y="74" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Univariate Baselines (Susceptible to Dropouts &amp; Noise)</text>
    </g>

    <!-- Tier 4: Somatic Mutational Burden (Rank 9+) -->
    <g transform="translate(1089, 784)">
      <rect width="277" height="74" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Snowflake icon -->
      {render_snowflake_icon(10, 8, 16)}

      <text x="32" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">TIER 4 (RANK 9+): Somatic TMB</text>
      <text x="32" y="28" font-size="7" fill="{COLOR_TEXT_SECONDARY}">Genomic Mutational Burden (Exome Sequencing)</text>

      <g transform="translate(12, 36)">
        <rect width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="38" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">ROC-AUC: 0.583</text>

        <rect x="82" y="0" width="76" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="120" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">PR-AUC: 0.345</text>

        <rect x="164" y="0" width="88" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="208" y="12" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">C-Index: 0.541</text>
      </g>

      <text x="12" y="66" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Clonal Somatic Mutations (Orthogonal Signal Modality)</text>
    </g>

    <!-- Bottom Status Legend -->
    <g transform="translate(1089, 868)">
      <rect width="277" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      {render_flame_icon(12, 7, 16)}
      <text x="32" y="19" font-size="7" font-weight="700" fill="{OKABE_BLUE}">Calibrated Simplex Solver</text>

      {render_snowflake_icon(152, 7, 16)}
      <text x="174" y="19" font-size="7" font-weight="600" fill="{COLOR_TEXT_MUTED}">Fixed / Uncalibrated</text>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 07: Inter-Panel Connectors & Pipeline Flow
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 07: INTER-PANEL CONNECTORS & FLOW -->
  <g inkscape:groupmode="layer" id="layer-07-connectors" inkscape:label="07_Flow_Connectors">
    <!-- Gutter Connector: Panel A to Panel B (Top) -->
    <path d="M 450 255 L 465 255" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.2" marker-end="url(#arr-slate)" />

    <!-- Gutter Connector: Panel B to Panel C (Middle Right) -->
    <path d="M 920 448 L 920 460" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.2" marker-end="url(#arr-slate)" />

    <!-- Gutter Connector: Panel C to Panel D (Bottom Horizontal) -->
    <path d="M 705 690 L 690 690" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.2" marker-end="url(#arr-slate)" />

    <!-- Gutter Connector: Panel D to Panel E (Bottom Horizontal) -->
    <path d="M 1060 690 L 1075 690" stroke="{COLOR_TEXT_MUTED}" stroke-width="1.2" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the graphical abstract SVG file."""
    output_dir = Path("output/benchmarks/iatlas_roc_auc")
    output_dir.mkdir(parents=True, exist_ok=True)
    out_svg_path = output_dir / "fig_graphical_abstract_section6.svg"

    svg_content = render_svg()
    out_svg_path.write_text(svg_content, encoding="utf-8")
    print(f"Successfully generated Section 6 Graphical Abstract SVG: {out_svg_path} ({len(svg_content)} bytes)")


if __name__ == "__main__":
    main()
