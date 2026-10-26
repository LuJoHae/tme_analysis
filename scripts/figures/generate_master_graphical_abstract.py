#!/usr/bin/env python3
"""
Generate Master Publication-Grade Graphical Abstract SVG (Restructured & Technical Version).

Design Specifications:
- DIN A4 Textbody Dimensions: width="175mm", height="170mm", viewBox="0 0 1400 1360".
- Single unified figure: zero panel segregation boxes, zero lowercase letter badges (a, b, c, d).
- Nature Methods minimalist wireframe standard: pure white canvas (#FFFFFF), 0.75 pt hairlines (#CBD5E1), rx=0, zero drop shadows.
- Okabe-Ito Colorblind-Safe palette strictly on data marks.
- Formal LaTeX math symbols and serif typography hierarchy (Times New Roman / Georgia).
- RESTRUCTURED LAYOUT:
    Left Column (x=30, y=78, w=510, h=538):
        Bulk RNA-seq ICB Clinical Cohorts (Discovery & Harmonization)
        Spans downwards alongside both Single-Cell Database and Pipeline!
        Bulk Matrix B drops directly down into BayesPrism below it!
    Right Column Top (x=560, y=78, w=810, h=152):
        Solid Tumor Single-Cell Database (.h5ad Discovery & Repository)
    Right Column Bottom (x=560, y=246, w=810, h=370):
        Automated Single-Cell Processing Engine (Collapsed into right half, 2 stacked rows):
        Row 2A: 1. QC & Filtering -> 2. log(1+p) & Sampling -> 3. Technical Truncated SVD PCA -> 4. Dense kNN Graph
        Row 2B: 5. Technical 4-Phase Harmony -> 6. Technical 3-Phase Leiden -> 7. Cell State Matrix Phi
- Native Inkscape layers:
    00_Canvas_Background
    01_Figure_Header
    02_Data_Ingestion_Stream
    03_Automated_Pipeline_Stream
    04_Deconvolution_Prediction_Stream
    05_Triangulation_Concordance_Stream
    06_Flow_Conduits
"""

from __future__ import annotations

from pathlib import Path
import html
import sys

# Ensure scripts/figures is in sys.path for nature_style_config
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_style_config import (
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_CONTAINER_BG,
    COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_YELLOW, OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)

# Semantic Colors
COLOR_RESPONDER = OKABE_BLUISH_GREEN    # #009E73 (CR / PR)
COLOR_NON_RESPONDER = OKABE_VERMILION   # #D55E00 (SD / PD)
COLOR_MELANOMA = OKABE_ORANGE           # #E69F00
COLOR_RCC = OKABE_BLUISH_GREEN          # #009E73
COLOR_BLADDER = OKABE_SKY_BLUE          # #56B4E9
COLOR_NSCLC = OKABE_BLUE                # #0072B2
COLOR_CRC = OKABE_REDDISH_PURPLE        # #CC79A7
COLOR_TNBC = OKABE_VERMILION            # #D55E00
COLOR_HNSCC = "#DDAA33"                 # Amber / Ochre
COLOR_CD8 = OKABE_BLUE                  # #0072B2
COLOR_MCMC = OKABE_REDDISH_PURPLE       # #CC79A7

CANVAS_WIDTH = 1400
CANVAS_HEIGHT = 1360


def svg_defs() -> str:
    """Return standard SVG defs containing hairline arrow markers."""
    return f"""
  <defs>
    <marker id="arr-slate" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
    </marker>
    <marker id="arr-blue" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
    <marker id="arr-green" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_RESPONDER}" />
    </marker>
    <marker id="arr-vermilion" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_NON_RESPONDER}" />
    </marker>
    <marker id="arr-purple" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_MCMC}" />
    </marker>
  </defs>
"""


def render_patient_glyph(x: float, y: float, color: str, badge: str | None = None, scale: float = 1.0) -> str:
    """Render a crisp vector patient silhouette pictogram with optional response badge."""
    badge_markup = ""
    if badge == "R":
        badge_markup = f"""
      <circle cx="{x + 12.5 * scale}" cy="{y + 5.5 * scale}" r="{3.5 * scale}" fill="{COLOR_RESPONDER}" />
      <path d="M {x + 11.2 * scale} {y + 5.5 * scale} L {x + 12.2 * scale} {y + 6.7 * scale} L {x + 14.2 * scale} {y + 4.3 * scale}" 
            stroke="#FFFFFF" stroke-width="{0.9 * scale}" fill="none" stroke-linecap="round" stroke-linejoin="round" />
        """
    elif badge == "NR":
        badge_markup = f"""
      <circle cx="{x + 12.5 * scale}" cy="{y + 5.5 * scale}" r="{3.5 * scale}" fill="{COLOR_NON_RESPONDER}" />
      <path d="M {x + 11.2 * scale} {y + 4.2 * scale} L {x + 13.8 * scale} {y + 6.8 * scale} M {x + 13.8 * scale} {y + 4.2 * scale} L {x + 11.2 * scale} {y + 6.8 * scale}" 
            stroke="#FFFFFF" stroke-width="{0.9 * scale}" stroke-linecap="round" />
        """

    return f"""
    <g class="patient-glyph">
      <!-- Head -->
      <circle cx="{x + 6.5 * scale}" cy="{y + 4 * scale}" r="{3 * scale}" fill="{color}" />
      <!-- Torso & Limbs -->
      <path d="M {x + 2.5 * scale} {y + 8.5 * scale} 
               C {x + 2.5 * scale} {y + 7.2 * scale}, {x + 4.2 * scale} {y + 6.8 * scale}, {x + 6.5 * scale} {y + 6.8 * scale} 
               C {x + 8.8 * scale} {y + 6.8 * scale}, {x + 10.5 * scale} {y + 7.2 * scale}, {x + 10.5 * scale} {y + 8.5 * scale} 
               L {x + 10 * scale} {y + 13.5 * scale} 
               L {x + 8.5 * scale} {y + 13.5 * scale} 
               L {x + 8.5 * scale} {y + 19.5 * scale} 
               L {x + 7.2 * scale} {y + 19.5 * scale} 
               L {x + 7.2 * scale} {y + 14 * scale} 
               L {x + 5.8 * scale} {y + 14 * scale} 
               L {x + 5.8 * scale} {y + 19.5 * scale} 
               L {x + 4.5 * scale} {y + 19.5 * scale} 
               L {x + 4.5 * scale} {y + 13.5 * scale} 
               L {x + 3 * scale} {y + 13.5 * scale} Z" 
            fill="{color}" />
      {badge_markup}
    </g>
    """


def render_master_graphical_abstract() -> str:
    """Render the master publication-grade Graphical Abstract SVG."""
    parts: list[str] = []

    # SVG Root Header
    parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="175mm" height="170mm" viewBox="0 0 {CANVAS_WIDTH} {CANVAS_HEIGHT}"
     version="1.1" id="figure-graphical-abstract-master"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
""")

    # --------------------------------------------------------------------------
    # LAYER 00: Canvas Background
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{CANVAS_WIDTH}" height="{CANVAS_HEIGHT}" fill="{COLOR_CANVAS_BG}" />
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 01: Figure Header Banner
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 01: FIGURE HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <g id="grp-header" transform="translate(30, 18)">
      <!-- Main Title -->
      <text x="0" y="20" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">
        Decoupling Tumor Microenvironment Heterogeneity: Automated Single-Cell Processing, Hierarchical Deconvolution &amp; In Silico Triangulation
      </text>
      <!-- Subtitle -->
      <text x="0" y="38" font-size="10" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        End-to-end framework integrating solid tumor single-cell atlases, BayesPrism deconvolution, stability-selected immunotherapy response modeling, and Milo graph calibration.
      </text>

      <!-- Top Right Badges -->
      <g transform="translate(1140, 4)">
        <rect x="0" y="0" width="200" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <circle cx="16" cy="14" r="4" fill="{OKABE_BLUE}" />
        <text x="26" y="18" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Graphical Abstract</text>
        <text x="145" y="18" font-size="7.5" fill="{COLOR_TEXT_MUTED}">DIN A4 Textbody</text>
      </g>

      <!-- Dividing Rule -->
      <line x1="0" y1="50" x2="1340" y2="50" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
    </g>
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 02: Data Ingestion & Storage Stream (y = 78 to 616, h = 538)
    # Reorganized: Left Column (Bulk) & Right Column Top (Single-Cell Database)
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 02: DATA INGESTION STREAM -->
  <g inkscape:groupmode="layer" id="layer-02-data-ingestion" inkscape:label="02_Data_Ingestion_Stream">
    <!-- ====================================================================== -->
    <!-- LEFT COLUMN: Bulk RNA-seq ICB Cohorts (x=30, y=78, w=510, h=538)       -->
    <!-- Spans downwards alongside both Single-Cell Database and Pipeline!      -->
    <!-- ====================================================================== -->
    <g id="stream-1a-bulk-icb" transform="translate(30, 78)">
      <rect width="510" height="538" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Banner -->
      <rect x="0" y="0" width="510" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Bulk RNA-seq ICB Clinical Cohorts (Discovery &amp; Harmonization)
      </text>
      <text x="390" y="16" font-size="7.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">N = 1,097 Patients</text>

      <!-- Block 1: All 7 Cancer Indications with Patient Pictograms (y = 30 to 105) -->
      <g transform="translate(14, 30)">
        <!-- Row 1: Melanoma, ccRCC, Bladder, NSCLC -->
        <g transform="translate(0, 0)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_MELANOMA}">Melanoma (N=485)</text>
          {render_patient_glyph(0, 11, COLOR_MELANOMA, badge="R", scale=0.7)}
          {render_patient_glyph(13, 11, COLOR_MELANOMA, badge="NR", scale=0.7)}
          {render_patient_glyph(26, 11, COLOR_MELANOMA, badge="R", scale=0.7)}
        </g>

        <g transform="translate(125, 0)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_RCC}">ccRCC (N=263)</text>
          {render_patient_glyph(0, 11, COLOR_RCC, badge="R", scale=0.7)}
          {render_patient_glyph(13, 11, COLOR_RCC, badge="NR", scale=0.7)}
          {render_patient_glyph(26, 11, COLOR_RCC, badge="R", scale=0.7)}
        </g>

        <g transform="translate(245, 0)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_BLADDER}">mUC / Bladder (N=349)</text>
          {render_patient_glyph(0, 11, COLOR_BLADDER, badge="R", scale=0.7)}
          {render_patient_glyph(13, 11, COLOR_BLADDER, badge="NR", scale=0.7)}
          {render_patient_glyph(26, 11, COLOR_BLADDER, badge="NR", scale=0.7)}
        </g>

        <g transform="translate(375, 0)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_NSCLC}">NSCLC (Lung)</text>
          {render_patient_glyph(0, 11, COLOR_NSCLC, badge="R", scale=0.7)}
          {render_patient_glyph(13, 11, COLOR_NSCLC, badge="NR", scale=0.7)}
        </g>

        <!-- Row 2: CRC, TNBC, HNSCC + Response Legend -->
        <g transform="translate(0, 38)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_CRC}">CRC</text>
          {render_patient_glyph(0, 11, COLOR_CRC, badge="R", scale=0.7)}
          {render_patient_glyph(13, 11, COLOR_CRC, badge="NR", scale=0.7)}
        </g>

        <g transform="translate(125, 38)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_TNBC}">TNBC</text>
          {render_patient_glyph(0, 11, COLOR_TNBC, badge="R", scale=0.7)}
        </g>

        <g transform="translate(245, 38)">
          <text x="0" y="8" font-size="6.5" font-weight="700" fill="{COLOR_HNSCC}">HNSCC</text>
          {render_patient_glyph(0, 11, COLOR_HNSCC, badge="NR", scale=0.7)}
        </g>

        <!-- Response Status Legend (2-row vertical card) -->
        <g transform="translate(365, 36)">
          <rect width="115" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="10" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Response Status:</text>
          <!-- Row 1: Responder -->
          <circle cx="12" cy="20" r="3" fill="{COLOR_RESPONDER}" />
          <text x="18" y="22.5" font-size="5.8" font-weight="600" fill="{COLOR_RESPONDER}">Responder (CR/PR)</text>
          <!-- Row 2: Non-Responder -->
          <circle cx="12" cy="30" r="3" fill="{COLOR_NON_RESPONDER}" />
          <text x="18" y="32.5" font-size="5.8" font-weight="600" fill="{COLOR_NON_RESPONDER}">Non-Resp (SD/PD)</text>
        </g>
      </g>

      <!-- Dividing Rule -->
      <line x1="14" y1="110" x2="496" y2="110" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Block 2: Biopsy -> 8-Row Bulk Matrix B -> RECIST Labels (y = 118 to 275) -->
      <g transform="translate(14, 118)">
        <!-- Biopsy Tube -->
        <g transform="translate(0, 42)">
          <rect x="0" y="0" width="16" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
          <rect x="2" y="16" width="12" height="24" fill="{COLOR_MELANOMA}" opacity="0.4" />
          <text x="8" y="52" font-size="5.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Biopsy</text>
        </g>

        <path d="M 20 62 L 30 62" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- 8-Row Bulk Expression Matrix B -->
        <g transform="translate(36, 0)">
          <rect width="144" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="72" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
            Bulk Matrix <tspan font-family="{FONT_SERIF_MATH}" font-weight="bold">B</tspan> ∈ ℝ<tspan baseline-shift="super" font-size="70%">G × N</tspan>
          </text>
          
          <!-- 8 Rows of Matrix Cells -->
          <g transform="translate(10, 18)">
            <!-- Row 1: CD8A -->
            <rect x="0" y="0" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.9" />
            <rect x="18" y="0" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.3" />
            <rect x="36" y="0" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.85" />
            <rect x="54" y="0" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.4" />
            <rect x="72" y="0" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.1" />
            <text x="96" y="7" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD8A</text>

            <!-- Row 2: IFNG -->
            <rect x="0" y="11" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.75" />
            <rect x="18" y="11" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.2" />
            <rect x="36" y="11" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.8" />
            <rect x="54" y="11" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.3" />
            <rect x="72" y="11" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.15" />
            <text x="96" y="18" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">IFNG</text>

            <!-- Row 3: CXCL9 -->
            <rect x="0" y="22" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.95" />
            <rect x="18" y="22" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.1" />
            <rect x="36" y="22" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.9" />
            <rect x="54" y="22" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.5" />
            <rect x="72" y="22" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.05" />
            <text x="96" y="29" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CXCL9</text>

            <!-- Row 4: PDCD1 -->
            <rect x="0" y="33" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.6" />
            <rect x="18" y="33" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.8" />
            <rect x="36" y="33" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.4" />
            <rect x="54" y="33" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.7" />
            <rect x="72" y="33" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.2" />
            <text x="96" y="40" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">PDCD1</text>

            <!-- Row 5: CTLA4 -->
            <rect x="0" y="44" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.5" />
            <rect x="18" y="44" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.6" />
            <rect x="36" y="44" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.3" />
            <rect x="54" y="44" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.65" />
            <rect x="72" y="44" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.3" />
            <text x="96" y="51" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CTLA4</text>

            <!-- Row 6: TCF7 -->
            <rect x="0" y="55" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.85" />
            <rect x="18" y="55" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.1" />
            <rect x="36" y="55" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.7" />
            <rect x="54" y="55" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.2" />
            <rect x="72" y="55" width="16" height="8.5" fill="{OKABE_BLUE}" opacity="0.1" />
            <text x="96" y="62" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">TCF7</text>

            <!-- Row 7: CD163 -->
            <rect x="0" y="66" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.1" />
            <rect x="18" y="66" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.9" />
            <rect x="36" y="66" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.2" />
            <rect x="54" y="66" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.8" />
            <rect x="72" y="66" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.7" />
            <text x="96" y="73" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD163</text>

            <!-- Row 8: MLANA -->
            <rect x="0" y="77" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.2" />
            <rect x="18" y="77" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.85" />
            <rect x="36" y="77" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.1" />
            <rect x="54" y="77" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.7" />
            <rect x="72" y="77" width="16" height="8.5" fill="{COLOR_NON_RESPONDER}" opacity="0.9" />
            <text x="96" y="84" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">MLANA</text>
          </g>
          <text x="72" y="138" font-size="5.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">G signature genes × N bulk biopsies</text>
        </g>

        <!-- RECIST Response Vector -->
        <g transform="translate(186, 0)">
          <rect width="134" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="67" y="14" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">RECIST Labels</text>
          <text x="67" y="27" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
            <tspan font-style="italic">y</tspan> ∈ {{0, 1}}<tspan baseline-shift="super" font-size="70%">N</tspan>
          </text>
          <g transform="translate(12, 34)">
            <rect width="110" height="16" fill="{COLOR_RESPONDER}" />
            <text x="55" y="11" font-size="6.2" font-weight="700" fill="#FFFFFF" text-anchor="middle">R (1): CR / PR</text>
            
            <rect y="21" width="110" height="16" fill="{COLOR_NON_RESPONDER}" />
            <text x="55" y="32" font-size="6.2" font-weight="700" fill="#FFFFFF" text-anchor="middle">NR (0): SD / PD</text>
          </g>
          <text x="67" y="88" font-size="5.6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Clinical Endpoints:</text>
          <text x="67" y="100" font-size="5.5" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Objective Radiographic Response</text>
          <text x="67" y="112" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Progression-Free Survival (PFS)</text>
          <text x="67" y="124" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Overall Survival (OS)</text>
          <text x="67" y="138" font-size="5.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Harmonized RECIST 1.1 Criteria</text>
        </g>

        <!-- Study Design Overview Card -->
        <g transform="translate(326, 0)">
          <rect width="154" height="150" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Clinical Multi-Trial Design:</text>
          
          <rect x="10" y="24" width="134" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="16" y="35" font-size="5.8" font-weight="600" fill="{OKABE_BLUE}">7 Solid Tumor Types</text>
          <text x="16" y="44" font-size="5.2" fill="{COLOR_TEXT_MUTED}">Melanoma, ccRCC, Bladder, NSCLC...</text>

          <rect x="10" y="52" width="134" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="16" y="63" font-size="5.8" font-weight="600" fill="{COLOR_RESPONDER}">9 Multicenter Trials</text>
          <text x="16" y="72" font-size="5.2" fill="{COLOR_TEXT_MUTED}">anti-PD-1, anti-PD-L1, anti-CTLA-4</text>

          <rect x="10" y="80" width="134" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="16" y="91" font-size="5.8" font-weight="600" fill="{OKABE_ORANGE}">Pre-Treatment Biopsies</text>
          <text x="16" y="100" font-size="5.2" fill="{COLOR_TEXT_MUTED}">Standardized RNA-seq library profiles</text>

          <text x="10" y="122" font-size="5.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Cross-Cohort Calibration:</text>
          <text x="10" y="134" font-size="5.2" fill="{COLOR_TEXT_SECONDARY}">Enables generalizable response modeling</text>
          <text x="10" y="144" font-size="5.2" fill="{COLOR_TEXT_MUTED}">across independent clinical centres</text>
        </g>
      </g>

      <!-- Dividing Rule -->
      <line x1="14" y1="278" x2="496" y2="278" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Block 3: All 9 Standardized ICB Clinical Trials Table (y = 286 to 486, h = 200) -->
      <g transform="translate(14, 286)">
        <rect width="482" height="198" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          All 9 Standardized ICB Clinical Trials (Pre-Treatment Biopsies):
        </text>

        <!-- Table Column Headers -->
        <rect x="10" y="20" width="462" height="15" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="16" y="30" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Trial Cohort</text>
        <text x="135" y="30" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cancer Indication</text>
        <text x="260" y="30" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Immunotherapy Regimen</text>
        <text x="415" y="30" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">Sample Size (N)</text>

        <!-- Table Rows (9 Standardized Trials, each 16px row) -->
        <g transform="translate(10, 37)">
          <!-- Row 1: Riaz -->
          <text x="6" y="10" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Riaz et al. (2017)</text>
          <text x="125" y="10" font-size="5.8" fill="{COLOR_MELANOMA}">Melanoma</text>
          <text x="250" y="10" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 (Nivolumab)</text>
          <text x="405" y="10" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 98</text>
          <line x1="0" y1="14" x2="462" y2="14" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 2: Gide -->
          <text x="6" y="25" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Gide et al. (2019)</text>
          <text x="125" y="25" font-size="5.8" fill="{COLOR_MELANOMA}">Melanoma</text>
          <text x="250" y="25" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 / anti-CTLA-4</text>
          <text x="405" y="25" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 73</text>
          <line x1="0" y1="29" x2="462" y2="29" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 3: Hugo -->
          <text x="6" y="40" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Hugo et al. (2016)</text>
          <text x="125" y="40" font-size="5.8" fill="{COLOR_MELANOMA}">Melanoma</text>
          <text x="250" y="40" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 (Pembrolizumab)</text>
          <text x="405" y="40" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 27</text>
          <line x1="0" y1="44" x2="462" y2="44" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 4: VanAllen -->
          <text x="6" y="55" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">VanAllen et al. (2015)</text>
          <text x="125" y="55" font-size="5.8" fill="{COLOR_MELANOMA}">Melanoma</text>
          <text x="250" y="55" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-CTLA-4 (Ipilimumab)</text>
          <text x="405" y="55" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 42</text>
          <line x1="0" y1="59" x2="462" y2="59" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 5: Liu -->
          <text x="6" y="70" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Liu et al. (2019)</text>
          <text x="125" y="70" font-size="5.8" fill="{COLOR_MELANOMA}">Melanoma</text>
          <text x="250" y="70" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 (Nivolumab / Pembro)</text>
          <text x="405" y="70" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 121</text>
          <line x1="0" y1="74" x2="462" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 6: McDermott -->
          <text x="6" y="85" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">McDermott et al. (2018)</text>
          <text x="125" y="85" font-size="5.8" fill="{COLOR_RCC}">ccRCC</text>
          <text x="250" y="85" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-L1 (Atezolizumab)</text>
          <text x="405" y="85" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 247</text>
          <line x1="0" y1="89" x2="462" y2="89" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 7: Braun -->
          <text x="6" y="100" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Braun et al. (2020)</text>
          <text x="125" y="100" font-size="5.8" fill="{COLOR_RCC}">ccRCC</text>
          <text x="250" y="100" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 (Nivolumab)</text>
          <text x="405" y="100" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 181</text>
          <line x1="0" y1="104" x2="462" y2="104" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 8: Choueiri -->
          <text x="6" y="115" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Choueiri et al. (2016)</text>
          <text x="125" y="115" font-size="5.8" fill="{COLOR_RCC}">ccRCC</text>
          <text x="250" y="115" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-1 (Nivolumab)</text>
          <text x="405" y="115" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 16</text>
          <line x1="0" y1="119" x2="462" y2="119" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Row 9: Mariathasan -->
          <text x="6" y="130" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Mariathasan et al. (2018)</text>
          <text x="125" y="130" font-size="5.8" fill="{COLOR_BLADDER}">mUC (Bladder)</text>
          <text x="250" y="130" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">anti-PD-L1 (Atezolizumab)</text>
          <text x="405" y="130" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="end">N = 348</text>
          <line x1="0" y1="134" x2="462" y2="134" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.3" />

          <!-- Summary Footnote -->
          <text x="6" y="148" font-size="5.2" fill="{COLOR_TEXT_MUTED}">
            Total: 9 Standardized Cohorts | 1,097 Patients with Complete Pre-Treatment Transcriptomes &amp; RECIST Response
          </text>
        </g>
      </g>

      <!-- Block 4: Bottom Harmonization Invariant Rule (y = 494) -->
      <g transform="translate(14, 494)">
        <line x1="0" y1="0" x2="482" y2="0" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        <text x="0" y="12" font-size="6.2" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Harmonization Invariant:</text>
        <text x="105" y="12" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">
          Standardized Ensembl 111 (GRCh38) gene mapping with cross-cohort RECIST radiographic outcome harmonization.
        </text>
        <text x="0" y="24" font-size="5.2" fill="{COLOR_TEXT_MUTED}">
          Feeds directly into BayesPrism deconvolution engine below with zero inter-panel cross-talk.
        </text>
      </g>
    </g>

    <!-- ====================================================================== -->
    <!-- RIGHT COLUMN TOP: Solid Tumor Single-Cell Database (x=560, y=78, w=810, h=152) -->
    <!-- ====================================================================== -->
    <g id="stream-1b-sc-database" transform="translate(560, 78)">
      <rect width="810" height="152" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Banner -->
      <rect x="0" y="0" width="810" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Solid Tumor Single-Cell Database (<tspan font-family="{FONT_SANS}">.h5ad</tspan> Discovery &amp; Repository)
      </text>
      <text x="590" y="16" font-size="7.5" fill="{COLOR_TEXT_MUTED}">GEO, CELLxGENE, &amp; Zenodo Integration</text>

      <!-- Cylinder Database Graphic -->
      <g transform="translate(16, 32)">
        <ellipse cx="26" cy="10" rx="24" ry="6" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <rect x="2" y="10" width="48" height="68" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <ellipse cx="26" cy="78" rx="24" ry="6" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <path d="M 2 32 A 24 6 0 0 0 50 32" fill="none" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <path d="M 2 55 A 24 6 0 0 0 50 55" fill="none" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <text x="26" y="24" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">.h5ad</text>
        <text x="26" y="46" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">AnnData</text>
        <text x="26" y="69" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">HDF5</text>
      </g>

      <!-- Detailed AnnData Data Specification -->
      <g transform="translate(76, 30)">
        <rect width="520" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="13" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Standardized Dual-Layer AnnData Specification:</text>
        
        <!-- .X slot -->
        <rect x="12" y="18" width="496" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="18" y="29.5" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">adata.X</text>
        <text x="80" y="29.5" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Raw integer UMI count matrix (CSR compressed sparse row)</text>

        <!-- .layers slot -->
        <rect x="12" y="37" width="496" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="18" y="48.5" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_MCMC}">layers["log1p_norm"]</text>
        <text x="135" y="48.5" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">CP10K Library-normalized: ln( 1 + 10<tspan baseline-shift="super" font-size="70%">4</tspan>·counts/depth )</text>

        <!-- .obs slot -->
        <rect x="12" y="56" width="496" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="18" y="67.5" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">adata.obs</text>
        <text x="80" y="67.5" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">patient_id, sample_id, clinical_response, cell_type, pct_counts_mt</text>

        <!-- .var slot -->
        <rect x="12" y="75" width="496" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="18" y="86.5" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_ORANGE}">adata.var</text>
        <text x="80" y="86.5" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Ensembl Release 111 mapped HGNC gene symbols (GRCh38)</text>
      </g>

      <!-- Metadata Audit Stats Card -->
      <g transform="translate(606, 30)">
        <rect width="190" height="96" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="95" y="14" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Metadata Audit</text>
        
        <g transform="translate(20, 22)">
          <text x="0" y="16" font-size="16" font-weight="700" fill="{OKABE_BLUE}">350+</text>
          <text x="0" y="27" font-size="6" fill="{COLOR_TEXT_MUTED}">Solid Tumor Cohorts</text>
        </g>

        <g transform="translate(105, 22)">
          <text x="0" y="16" font-size="16" font-weight="700" fill="{COLOR_RESPONDER}">9</text>
          <text x="0" y="27" font-size="6" fill="{COLOR_TEXT_MUTED}">Verified ICB Trials</text>
        </g>

        <line x1="15" y1="56" x2="175" y2="56" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        <text x="95" y="68" font-size="6" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">100% Patient Response Audit</text>
        <text x="95" y="78" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Standardized Ensembl 111 Invariant</text>
        <text x="95" y="88" font-size="5.2" fill="{OKABE_BLUE}" text-anchor="middle">Direct feed to Automated Pipeline</text>
      </g>

      <!-- Bottom Harmonization Standard -->
      <g transform="translate(16, 134)">
        <line x1="0" y1="0" x2="780" y2="0" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        <text x="0" y="11" font-size="6.2" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Universal Harmonization Standard:</text>
        <text x="135" y="11" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">
          Ensembl Release 111 (GRCh38) gene nomenclature with strictly verified RECIST radiographic criteria.
        </text>
      </g>
    </g>
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 03: Automated Pipeline Stream (Right Column Bottom: x=560, y=246, w=810, h=370)
    # Ordered in TWO ROWS: Row 2A (Steps 1, 2, 3, 4) & Row 2B (Steps 5, 6, 7)
    # Technical Vector Representations of Truncated SVD PCA, Harmony, and Leiden
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 03: AUTOMATED PIPELINE STREAM (RIGHT HALF) -->
  <g inkscape:groupmode="layer" id="layer-03-pipeline" inkscape:label="03_Automated_Pipeline_Stream">
    <g id="stream-02-sc-pipeline" transform="translate(560, 246)">
      <rect width="810" height="370" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Section Label Banner -->
      <rect x="0" y="0" width="810" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Automated Single-Cell Processing Engine &amp; Cell State Profile Generation (Standardized Pipeline)
      </text>
      <text x="590" y="16" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Scripts/01-03 Pipeline Engine</text>

      <!-- ================= ROW 2A: Steps 1, 2, 3, 4 (y = 30 to 185, h = 155) ================= -->
      <g transform="translate(12, 30)">
        <!-- Step 1: QC & Filtering (w = 188, h = 150) -->
        <g transform="translate(0, 0)">
          <rect width="188" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="188" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="94" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">1. Quality Control &amp; Filtering</text>

          <g transform="translate(10, 22)">
            <text x="0" y="12" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              n<tspan baseline-shift="sub" font-size="70%">counts</tspan> ≥ 500 <tspan font-family="{FONT_SANS}" font-size="5.5" fill="{COLOR_TEXT_MUTED}">(UMI depth)</tspan>
            </text>
            <text x="0" y="27" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
              200 ≤ n<tspan baseline-shift="sub" font-size="70%">genes</tspan> ≤ 8,000 <tspan font-family="{FONT_SANS}" font-size="5.5" fill="{COLOR_TEXT_MUTED}">(detection)</tspan>
            </text>
            <text x="0" y="42" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_NON_RESPONDER}">
              f<tspan baseline-shift="sub" font-size="70%">mito</tspan> &lt; 15% <tspan font-family="{FONT_SANS}" font-size="5.5" fill="{COLOR_TEXT_MUTED}">(debris &amp; apoptotic)</tspan>
            </text>
            <text x="0" y="57" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">
              f<tspan baseline-shift="sub" font-size="70%">ribo</tspan> &gt; 3% <tspan font-family="{FONT_SANS}" font-size="5.5" fill="{COLOR_TEXT_MUTED}">(translational capacity)</tspan>
            </text>
          </g>

          <!-- QC Viability Curve Graphic -->
          <g transform="translate(14, 88)">
            <rect width="160" height="38" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <path d="M 10 30 Q 35 6, 80 18 T 150 28" fill="none" stroke="{OKABE_BLUE}" stroke-width="1.2" />
            <line x1="50" y1="5" x2="50" y2="35" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />
            <text x="45" y="10" font-size="5" fill="{COLOR_NON_RESPONDER}" text-anchor="end">Cutoff</text>
            <text x="100" y="14" font-size="5.2" font-weight="600" fill="{COLOR_RESPONDER}">Retained Viable</text>
          </g>

          <text x="94" y="138" font-size="5.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Discards apoptotic debris &amp; multiplets</text>
          <text x="94" y="146" font-size="5.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">prior to dimensionality reduction</text>
        </g>

        <!-- Arrow 1 -> 2 -->
        <path d="M 191 75 L 199 75" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Step 2: Log(1p) Scaling & Sampling (w = 188, h = 150) -->
        <g transform="translate(202, 0)">
          <rect width="188" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="188" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="94" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">2. log(1+p) Scaling &amp; Sampling</text>

          <g transform="translate(10, 24)">
            <text x="0" y="10" font-size="6.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">CP10K Library Scaling:</text>
            <text x="0" y="24" font-size="7.2" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              x̃<tspan baseline-shift="sub" font-size="70%">c,g</tspan> = 10<tspan baseline-shift="super" font-size="70%">4</tspan> · x<tspan baseline-shift="sub" font-size="70%">c,g</tspan> / ∑<tspan baseline-shift="sub" font-size="65%">g</tspan> x<tspan baseline-shift="sub" font-size="70%">c,g</tspan>
            </text>

            <text x="0" y="40" font-size="6.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Variance Stabilization:</text>
            <text x="0" y="54" font-size="7.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_MCMC}">
              ln( 1 + x̃<tspan baseline-shift="sub" font-size="70%">c,g</tspan> )
            </text>

            <text x="0" y="70" font-size="6.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Sampling Engine:</text>
            <text x="0" y="82" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">Stratified quota subsampling</text>
            <text x="0" y="92" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">N ≤ 1,200 per lineage</text>
          </g>

          <g transform="translate(14, 122)">
            <rect width="160" height="20" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="80" y="13" font-size="5.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Equalizes lineage representation</text>
          </g>
        </g>

        <!-- Arrow 2 -> 3 -->
        <path d="M 393 75 L 401 75" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Step 3: Technical PCA Latent Space Projection (Truncated SVD) (w = 196, h = 150) -->
        <g transform="translate(404, 0)">
          <rect width="196" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="196" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="98" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">3. PCA Latent Space Projection</text>

          <g transform="translate(8, 22)">
            <text x="0" y="9" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              X̃ = UΣV<tspan baseline-shift="super" font-size="70%">⊤</tspan>,   z<tspan baseline-shift="sub" font-size="70%">i</tspan> = V<tspan baseline-shift="sub" font-size="70%">30</tspan><tspan baseline-shift="super" font-size="70%">⊤</tspan>(x<tspan baseline-shift="sub" font-size="70%">i</tspan> - μ) ∈ ℝ<tspan baseline-shift="super" font-size="70%">30</tspan>
            </text>
          </g>

          <!-- Detailed Technical 2D PCA Eigenspace Diagram -->
          <g transform="translate(10, 36)">
            <rect width="176" height="95" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

            <!-- Eigenspace Axes PC1 (v1) and PC2 (v2) -->
            <!-- Principal Axis PC1 (v1) -->
            <line x1="20" y1="78" x2="162" y2="22" stroke="{OKABE_BLUE}" stroke-width="1.2" marker-end="url(#arr-blue)" />
            <text x="166" y="24" font-size="5.8" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_BLUE}">v<tspan baseline-shift="sub" font-size="70%">1</tspan></text>

            <!-- Orthogonal Axis PC2 (v2) -->
            <line x1="85" y1="70" x2="65" y2="18" stroke="{OKABE_ORANGE}" stroke-width="1.0" marker-end="url(#arr-slate)" />
            <text x="65" y="14" font-size="5.8" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_ORANGE}">v<tspan baseline-shift="sub" font-size="70%">2</tspan></text>

            <!-- Orthogonal symbol at intersection -->
            <path d="M 82 58 L 86 56 L 88 62" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.6" fill="none" />

            <!-- Data Point Cloud & Orthogonal Projections -->
            <!-- Point 1 -->
            <circle cx="45" cy="58" r="2.2" fill="{COLOR_RESPONDER}" />
            <line x1="45" y1="58" x2="52" y2="65" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

            <!-- Point 2 -->
            <circle cx="65" cy="42" r="2.2" fill="{COLOR_RESPONDER}" />
            <line x1="65" y1="42" x2="72" y2="57" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

            <!-- Point 3 -->
            <circle cx="95" cy="55" r="2.2" fill="{COLOR_RESPONDER}" />
            <line x1="95" y1="55" x2="90" y2="50" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

            <!-- Point 4 -->
            <circle cx="115" cy="35" r="2.2" fill="{COLOR_NON_RESPONDER}" />
            <line x1="115" y1="35" x2="120" y2="38" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

            <!-- Point 5 -->
            <circle cx="135" cy="26" r="2.2" fill="{COLOR_NON_RESPONDER}" />
            <line x1="135" y1="26" x2="138" y2="31" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

            <!-- Point 6 -->
            <circle cx="78" cy="32" r="2.2" fill="{OKABE_ORANGE}" />

            <!-- Point 7 -->
            <circle cx="108" cy="48" r="2.2" fill="{OKABE_ORANGE}" />

            <!-- Scree Inset: Eigenvalue Variance Decay Curve (Top 30 PCs) -->
            <g transform="translate(118, 52)">
              <rect width="52" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="26" y="8" font-size="4.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Scree Spectrum</text>
              <path d="M 6 32 Q 12 14, 24 24 T 48 30" fill="none" stroke="{OKABE_BLUE}" stroke-width="0.9" />
              <line x1="28" y1="12" x2="28" y2="34" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.6" stroke-dasharray="1,1" />
              <text x="28" y="37" font-size="4" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">d=30</text>
            </g>
          </g>

          <text x="98" y="142" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Orthogonal eigenvectors capture biological covariance</text>
        </g>

        <!-- Arrow 3 -> 4 -->
        <path d="M 603 75 L 611 75" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Step 4: Manifold kNN Graph Construction (w = 170, h = 150) -->
        <g transform="translate(614, 0)">
          <rect width="170" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="170" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="85" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">4. kNN Graph Construction</text>

          <g transform="translate(8, 22)">
            <text x="0" y="9" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              k = 30 Nearest Neighbors in ℝ<tspan baseline-shift="super" font-size="70%">30</tspan>
            </text>
          </g>

          <!-- Dense 25-Node kNN Manifold Graphic -->
          <g transform="translate(8, 36)">
            <rect width="154" height="95" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

            <!-- Cluster A (Left, Blue) -->
            <line x1="20" y1="25" x2="45" y2="20" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="20" y1="25" x2="30" y2="55" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="45" y1="20" x2="30" y2="55" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="45" y1="20" x2="65" y2="35" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="30" y1="55" x2="65" y2="35" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="30" y1="55" x2="50" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="65" y1="35" x2="50" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

            <!-- Bridges across clusters -->
            <line x1="65" y1="35" x2="90" y2="40" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />
            <line x1="50" y1="70" x2="95" y2="65" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />

            <!-- Cluster B (Right, Orange) -->
            <line x1="90" y1="40" x2="115" y2="25" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="90" y1="40" x2="95" y2="65" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="115" y1="25" x2="135" y2="35" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="90" y1="40" x2="135" y2="35" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="95" y1="65" x2="130" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <line x1="135" y1="35" x2="130" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

            <!-- Vertices Cluster A -->
            <circle cx="20" cy="25" r="2.2" fill="{OKABE_BLUE}" />
            <circle cx="45" cy="20" r="2.5" fill="{OKABE_BLUE}" />
            <circle cx="30" cy="55" r="2.2" fill="{OKABE_BLUE}" />
            <circle cx="65" cy="35" r="3.2" fill="{OKABE_BLUE}" stroke="#0F172A" stroke-width="0.5" /> <!-- Hub -->
            <circle cx="50" cy="70" r="2.5" fill="{OKABE_BLUE}" />

            <!-- Vertices Cluster B -->
            <circle cx="90" cy="40" r="3.0" fill="{OKABE_ORANGE}" stroke="#0F172A" stroke-width="0.5" />
            <circle cx="115" cy="25" r="2.2" fill="{OKABE_ORANGE}" />
            <circle cx="95" cy="65" r="2.5" fill="{OKABE_ORANGE}" />
            <circle cx="135" cy="35" r="3.2" fill="{OKABE_ORANGE}" stroke="#0F172A" stroke-width="0.5" /> <!-- Hub -->
            <circle cx="130" cy="70" r="2.2" fill="{OKABE_ORANGE}" />
          </g>

          <text x="85" y="142" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Manifold topology G=(V,E) in ℝ³⁰</text>
        </g>
      </g>

      <!-- Connecting Flow Pipe: Down from Row 2A (Step 4) to Row 2B (Step 5) -->
      <path d="M 700 180 L 700 190 L 138 190 L 138 200" stroke="{OKABE_BLUE}" stroke-width="1.0" fill="none" marker-end="url(#arr-blue)" />

      <!-- ================= ROW 2B: Steps 5, 6, 7 (y = 200 to 355, h = 155) ================= -->
      <g transform="translate(12, 200)">
        <!-- Step 5: Technical Harmony Batch Correction (harmonypy) (w = 254, h = 155) -->
        <g transform="translate(0, 0)">
          <rect width="254" height="155" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="254" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="127" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">5. Harmony Batch Correction (harmonypy)</text>

          <!-- Mathematical Formula Header -->
          <g transform="translate(8, 20)">
            <text x="0" y="8" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              Soft K-Means: R<tspan baseline-shift="sub" font-size="70%">ki</tspan> ∝ exp(-||z<tspan baseline-shift="sub" font-size="70%">i</tspan> - μ<tspan baseline-shift="sub" font-size="70%">k</tspan>||<tspan baseline-shift="super" font-size="70%">2</tspan>/2σ<tspan baseline-shift="super" font-size="70%">2</tspan>) · (O<tspan baseline-shift="sub" font-size="70%">kb</tspan>/E<tspan baseline-shift="sub" font-size="70%">kb</tspan>)<tspan baseline-shift="super" font-size="70%">-θ</tspan>
            </text>
            <text x="0" y="18" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
              Linear Shift: ẑ<tspan baseline-shift="sub" font-size="70%">i</tspan> = z<tspan baseline-shift="sub" font-size="70%">i</tspan> - ∑<tspan baseline-shift="sub" font-size="65%">k</tspan> R<tspan baseline-shift="sub" font-size="70%">ki</tspan> w<tspan baseline-shift="sub" font-size="70%">k,b(i)</tspan>
            </text>
          </g>

          <!-- 3-Stage Technical Transformation Illustration -->
          <g transform="translate(8, 42)">
            <!-- Stage 1: Uncorrected Batched (Left) -->
            <g transform="translate(0, 0)">
              <rect width="74" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="37" y="10" font-size="5.2" font-weight="700" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">1. Uncorrected</text>

              <!-- Batch 1 (Orange Circles) & Batch 2 (Blue Squares) with Batch Offset delta -->
              <!-- T Cell Cluster Batched -->
              <circle cx="20" cy="26" r="2.2" fill="{OKABE_ORANGE}" />
              <circle cx="28" cy="22" r="2.2" fill="{OKABE_ORANGE}" />
              <rect x="36" y="32" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
              <rect x="44" y="28" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
              <!-- Batch offset vector delta -->
              <line x1="24" y1="24" x2="40" y2="34" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="1.5,1.5" marker-end="url(#arr-slate)" />
              <text x="36" y="44" font-size="4.5" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">δ<tspan baseline-shift="sub" font-size="70%">batch</tspan></text>

              <!-- Myeloid Cluster Batched -->
              <circle cx="22" cy="66" r="2.2" fill="{OKABE_ORANGE}" />
              <circle cx="30" cy="62" r="2.2" fill="{OKABE_ORANGE}" />
              <rect x="46" y="74" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
              <rect x="54" y="70" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
              <line x1="26" y1="64" x2="50" y2="72" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="1.5,1.5" marker-end="url(#arr-slate)" />

              <text x="37" y="90" font-size="4.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Cohort Shifts</text>
            </g>

            <!-- Arrow 1 -> 2 -->
            <path d="M 76 48 L 81 48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" marker-end="url(#arr-slate)" />

            <!-- Stage 2: Soft Clustering & Shift Vectors (Middle) -->
            <g transform="translate(83, 0)">
              <rect width="76" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="38" y="10" font-size="5.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">2. Diversity Shifts</text>

              <!-- Cluster Centroid mu_k (Crosshair) -->
              <circle cx="38" cy="30" r="14" fill="{OKABE_BLUE}" opacity="0.08" stroke="{OKABE_BLUE}" stroke-width="0.5" stroke-dasharray="2,2" />
              <path d="M 34 30 L 42 30 M 38 26 L 38 34" stroke="{OKABE_BLUE}" stroke-width="1.0" />
              <text x="44" y="27" font-size="5" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_BLUE}">μ<tspan baseline-shift="sub" font-size="70%">k</tspan></text>

              <!-- Shift Vectors w_kb pointing inward -->
              <circle cx="22" cy="25" r="2.0" fill="{OKABE_ORANGE}" />
              <line x1="22" y1="25" x2="35" y2="29" stroke="{COLOR_RESPONDER}" stroke-width="0.8" marker-end="url(#arr-green)" />
              <text x="26" y="20" font-size="4.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">w<tspan baseline-shift="sub" font-size="70%">k,1</tspan></text>

              <rect x="52" y="32" width="4" height="4" fill="{OKABE_BLUE}" />
              <line x1="52" y1="34" x2="41" y2="31" stroke="{COLOR_RESPONDER}" stroke-width="0.8" marker-end="url(#arr-green)" />
              <text x="50" y="44" font-size="4.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">w<tspan baseline-shift="sub" font-size="70%">k,2</tspan></text>

              <!-- Entropy penalty label -->
              <rect x="6" y="58" width="64" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.4" />
              <text x="38" y="67" font-size="4.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Batch Diversity (O/E)<tspan baseline-shift="super" font-size="70%">-θ</tspan></text>
              <text x="38" y="73.5" font-size="4.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Penalizes monotypic clusters</text>

              <text x="38" y="90" font-size="4.5" fill="{COLOR_RESPONDER}" text-anchor="middle">Ridge Optimization</text>
            </g>

            <!-- Arrow 2 -> 3 -->
            <path d="M 161 48 L 166 48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" marker-end="url(#arr-slate)" />

            <!-- Stage 3: Harmonized Integrated (Right) -->
            <g transform="translate(168, 0)">
              <rect width="70" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="35" y="10" font-size="5.2" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">3. Harmonized</text>

              <!-- Co-localized Integrated T Cells -->
              <circle cx="28" cy="27" r="2.0" fill="{OKABE_ORANGE}" />
              <rect x="36" y="24" width="4" height="4" fill="{OKABE_BLUE}" />
              <circle cx="34" cy="33" r="2.0" fill="{OKABE_ORANGE}" />
              <rect x="26" y="34" width="4" height="4" fill="{OKABE_BLUE}" />
              <circle cx="41" cy="29" r="2.0" fill="{OKABE_ORANGE}" />
              <text x="35" y="44" font-size="4.8" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">Integrated T</text>

              <!-- Co-localized Integrated Myeloid -->
              <circle cx="28" cy="62" r="2.0" fill="{OKABE_ORANGE}" />
              <rect x="36" y="59" width="4" height="4" fill="{OKABE_BLUE}" />
              <circle cx="33" cy="68" r="2.0" fill="{OKABE_ORANGE}" />
              <rect x="42" y="65" width="4" height="4" fill="{OKABE_BLUE}" />
              <text x="35" y="78" font-size="4.8" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">Integrated Myeloid</text>

              <text x="35" y="90" font-size="4.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Preserves Lineage</text>
            </g>
          </g>

          <text x="127" y="148" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Iterative clustering removes platform batch offsets without count distortion</text>
        </g>

        <!-- Arrow 5 -> 6 -->
        <path d="M 257 75 L 263 75" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Step 6: Technical Leiden Community Detection (w = 254, h = 155) -->
        <g transform="translate(266, 0)">
          <rect width="254" height="155" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="254" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="127" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">6. Leiden Community Detection</text>

          <!-- Mathematical Quality & Refinement Formulation -->
          <g transform="translate(8, 20)">
            <text x="0" y="8" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              CPM Quality: ℋ<tspan baseline-shift="sub" font-size="70%">CPM</tspan>(𝒫) = ∑<tspan baseline-shift="sub" font-size="65%">c</tspan> [ e(c) - γ(n<tspan baseline-shift="sub" font-size="70%">c</tspan> choose 2) ]
            </text>
            <text x="0" y="18" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">
              Refinement Guarantee: e(v, C'\\{{v\\}}) ≥ γ|v|(|C| - |v|)
            </text>
          </g>

          <!-- 3-Phase Algorithmic Loop Diagram -->
          <g transform="translate(8, 42)">
            <!-- Phase 1: Local Move -->
            <g transform="translate(0, 0)">
              <rect width="74" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="37" y="10" font-size="5.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Phase 1: Move</text>

              <!-- Node movement into neighbor cluster -->
              <ellipse cx="26" cy="35" rx="18" ry="14" fill="{OKABE_BLUE}" opacity="0.12" stroke="{OKABE_BLUE}" stroke-width="0.5" />
              <circle cx="18" cy="32" r="2.2" fill="{OKABE_BLUE}" />
              <circle cx="28" cy="28" r="2.2" fill="{OKABE_BLUE}" />
              <circle cx="24" cy="42" r="2.2" fill="{OKABE_BLUE}" />

              <!-- Moving Node v -->
              <circle cx="56" cy="48" r="2.8" fill="{OKABE_ORANGE}" stroke="#0F172A" stroke-width="0.6" />
              <path d="M 52 46 L 36 38" stroke="{OKABE_ORANGE}" stroke-width="0.8" marker-end="url(#arr-slate)" />
              <text x="58" y="58" font-size="4.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_ORANGE}">v</text>

              <text x="37" y="75" font-size="4.5" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">Greedy node move</text>
              <text x="37" y="83" font-size="4.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">Δℋ &gt; 0</text>
              <text x="37" y="91" font-size="4.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Fast modularity step</text>
            </g>

            <!-- Arrow Phase 1 -> Phase 2 -->
            <path d="M 76 48 L 81 48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" marker-end="url(#arr-slate)" />

            <!-- Phase 2: Refinement (Resolves Louvain Disconnection Defect!) -->
            <g transform="translate(83, 0)">
              <rect width="76" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="38" y="10" font-size="5.2" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">Phase 2: Refine</text>

              <!-- Louvain flaw resolution: splitting poorly connected sub-clusters -->
              <ellipse cx="24" cy="32" rx="14" ry="12" fill="{COLOR_RESPONDER}" opacity="0.12" stroke="{COLOR_RESPONDER}" stroke-width="0.6" stroke-dasharray="1.5,1.5" />
              <circle cx="18" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
              <circle cx="26" cy="34" r="2.0" fill="{COLOR_RESPONDER}" />
              <text x="24" y="48" font-size="4.5" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">C'<tspan baseline-shift="sub" font-size="70%">1</tspan></text>

              <!-- Weak bridge line -->
              <line x1="34" y1="32" x2="48" y2="32" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.6" stroke-dasharray="1.5,1.5" />
              <text x="41" y="28" font-size="4" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">Cut</text>

              <ellipse cx="58" cy="32" rx="14" ry="12" fill="{COLOR_RESPONDER}" opacity="0.12" stroke="{COLOR_RESPONDER}" stroke-width="0.6" stroke-dasharray="1.5,1.5" />
              <circle cx="54" cy="30" r="2.0" fill="{COLOR_RESPONDER}" />
              <circle cx="62" cy="36" r="2.0" fill="{COLOR_RESPONDER}" />
              <text x="58" y="48" font-size="4.5" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">C'<tspan baseline-shift="sub" font-size="70%">2</tspan></text>

              <rect x="4" y="58" width="68" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.4" />
              <text x="38" y="67" font-size="4.5" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">Well-Connectedness</text>
              <text x="38" y="74" font-size="4.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Guarantees connected</text>
              <text x="38" y="80" font-size="4.2" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">Eliminates Louvain defect</text>

              <text x="38" y="91" font-size="4.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Sub-community split</text>
            </g>

            <!-- Arrow Phase 2 -> Phase 3 -->
            <path d="M 161 48 L 166 48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" marker-end="url(#arr-slate)" />

            <!-- Phase 3: Aggregation & Final Partition -->
            <g transform="translate(168, 0)">
              <rect width="70" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
              <text x="35" y="10" font-size="5.2" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Phase 3: Output</text>

              <!-- Final 3 Guaranteed Connected Lineage Communities -->
              <!-- C1: CD8 T -->
              <ellipse cx="24" cy="28" rx="14" ry="10" fill="{OKABE_BLUE}" opacity="0.15" stroke="{OKABE_BLUE}" stroke-width="0.6" />
              <circle cx="20" cy="26" r="1.8" fill="{OKABE_BLUE}" />
              <circle cx="26" cy="30" r="1.8" fill="{OKABE_BLUE}" />
              <text x="24" y="39" font-size="4.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">CD8+ T</text>

              <!-- C2: CD4 T -->
              <ellipse cx="50" cy="28" rx="14" ry="10" fill="{OKABE_SKY_BLUE}" opacity="0.2" stroke="{OKABE_SKY_BLUE}" stroke-width="0.6" />
              <circle cx="48" cy="26" r="1.8" fill="{OKABE_SKY_BLUE}" />
              <circle cx="54" cy="30" r="1.8" fill="{OKABE_SKY_BLUE}" />
              <text x="50" y="39" font-size="4.5" font-weight="700" fill="{OKABE_SKY_BLUE}" text-anchor="middle">CD4+ T</text>

              <!-- Inter-community bridge -->
              <line x1="37" y1="28" x2="38" y2="48" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

              <!-- C3: Myeloid -->
              <ellipse cx="37" cy="60" rx="20" ry="12" fill="{OKABE_ORANGE}" opacity="0.15" stroke="{OKABE_ORANGE}" stroke-width="0.6" />
              <circle cx="28" cy="58" r="1.8" fill="{OKABE_ORANGE}" />
              <circle cx="37" cy="62" r="1.8" fill="{OKABE_ORANGE}" />
              <circle cx="46" cy="58" r="1.8" fill="{OKABE_ORANGE}" />
              <text x="37" y="73" font-size="4.5" font-weight="700" fill="{OKABE_ORANGE}" text-anchor="middle">Myeloid / Mac</text>

              <text x="35" y="85" font-size="4.2" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">Connected: 100%</text>
              <text x="35" y="91" font-size="4.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Optimal Partition 𝒫*</text>
            </g>
          </g>

          <text x="127" y="148" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Guarantees strictly connected, locally optimal modular communities</text>
        </g>

        <!-- Arrow 6 -> 7 -->
        <path d="M 522 75 L 528 75" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Step 7: Cell State Reference Signature Matrix (Phi) (w = 256, h = 155) -->
        <g transform="translate(532, 0)">
          <rect width="256" height="155" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect width="256" height="16" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="128" y="11" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">7. Cell State Matrix (Φ ∈ ℝ<tspan baseline-shift="super" font-size="70%">G×K</tspan>)</text>

          <g transform="translate(10, 20)">
            <text x="0" y="8" font-size="6.2" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
              Centroids: ϕ<tspan baseline-shift="sub" font-size="70%">k</tspan> = <tspan font-style="normal">1</tspan>/|C<tspan baseline-shift="sub" font-size="70%">k</tspan>| ∑<tspan baseline-shift="sub" font-size="65%">c∈C<tspan baseline-shift="sub" font-size="65%">k</tspan></tspan> x<tspan baseline-shift="sub" font-size="70%">c</tspan>
            </text>
          </g>

          <!-- Heatmap Representation -->
          <g transform="translate(8, 34)">
            <!-- Lineage Headers -->
            <rect x="0" y="0" width="36" height="7.5" fill="{OKABE_BLUE}" />
            <text x="18" y="5.5" font-size="4.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">CD8+ T</text>

            <rect x="38" y="0" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" />
            <text x="56" y="5.5" font-size="4.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">CD4+ T</text>

            <rect x="76" y="0" width="36" height="7.5" fill="{COLOR_RESPONDER}" />
            <text x="94" y="5.5" font-size="4.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">B Cells</text>

            <rect x="114" y="0" width="36" height="7.5" fill="{COLOR_MELANOMA}" />
            <text x="132" y="5.5" font-size="4.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Myeloid</text>

            <rect x="152" y="0" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" />
            <text x="170" y="5.5" font-size="4.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Malignant</text>

            <!-- Row 1: CD8A -->
            <rect x="0" y="10" width="36" height="7.5" fill="{OKABE_BLUE}" opacity="0.95" />
            <rect x="38" y="10" width="36" height="7.5" fill="{OKABE_BLUE}" opacity="0.1" />
            <rect x="76" y="10" width="36" height="7.5" fill="{OKABE_BLUE}" opacity="0.05" />
            <rect x="114" y="10" width="36" height="7.5" fill="{OKABE_BLUE}" opacity="0.05" />
            <rect x="152" y="10" width="36" height="7.5" fill="{OKABE_BLUE}" opacity="0.05" />
            <text x="195" y="16" font-size="5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD8A / CD8B</text>

            <!-- Row 2: FOXP3 -->
            <rect x="0" y="20" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" opacity="0.1" />
            <rect x="38" y="20" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" opacity="0.85" />
            <rect x="76" y="20" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
            <rect x="114" y="20" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
            <rect x="152" y="20" width="36" height="7.5" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
            <text x="195" y="26" font-size="5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD4 / FOXP3</text>

            <!-- Row 3: MS4A1 -->
            <rect x="0" y="30" width="36" height="7.5" fill="{COLOR_RESPONDER}" opacity="0.05" />
            <rect x="38" y="30" width="36" height="7.5" fill="{COLOR_RESPONDER}" opacity="0.05" />
            <rect x="76" y="30" width="36" height="7.5" fill="{COLOR_RESPONDER}" opacity="0.9" />
            <rect x="114" y="30" width="36" height="7.5" fill="{COLOR_RESPONDER}" opacity="0.05" />
            <rect x="152" y="30" width="36" height="7.5" fill="{COLOR_RESPONDER}" opacity="0.05" />
            <text x="195" y="36" font-size="5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">MS4A1 (CD20)</text>

            <!-- Row 4: CD68 -->
            <rect x="0" y="40" width="36" height="7.5" fill="{COLOR_MELANOMA}" opacity="0.05" />
            <rect x="38" y="40" width="36" height="7.5" fill="{COLOR_MELANOMA}" opacity="0.05" />
            <rect x="76" y="40" width="36" height="7.5" fill="{COLOR_MELANOMA}" opacity="0.05" />
            <rect x="114" y="40" width="36" height="7.5" fill="{COLOR_MELANOMA}" opacity="0.9" />
            <rect x="152" y="40" width="36" height="7.5" fill="{COLOR_MELANOMA}" opacity="0.05" />
            <text x="195" y="46" font-size="5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD68 / CD163</text>

            <!-- Row 5: MLANA -->
            <rect x="0" y="50" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
            <rect x="38" y="50" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
            <rect x="76" y="50" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
            <rect x="114" y="50" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
            <rect x="152" y="50" width="36" height="7.5" fill="{COLOR_NON_RESPONDER}" opacity="0.95" />
            <text x="195" y="56" font-size="5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">MLANA / PMEL</text>
          </g>

          <g transform="translate(10, 108)">
            <rect width="236" height="38" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="118" y="13" font-size="5.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">
              Simplex Invariant: ∑<tspan baseline-shift="sub" font-size="65%">g=1</tspan><tspan baseline-shift="super" font-size="65%">G</tspan> ϕ<tspan baseline-shift="sub" font-size="65%">k,g</tspan> = 1,  ϕ<tspan baseline-shift="sub" font-size="65%">k,g</tspan> ≥ 0
            </text>
            <text x="118" y="24" font-size="5.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Reference Matrix Φ (G signature genes × K states)</text>
            <text x="118" y="33" font-size="5.2" fill="{COLOR_RESPONDER}" text-anchor="middle">Ready for Hierarchical Deconvolution</text>
          </g>
        </g>
      </g>
    </g>
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 04: Bulk Deconvolution & Response Prediction (Row 3: y = 635 to 885, h = 250)
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 04: DECONVOLUTION & RESPONSE PREDICTION STREAM -->
  <g inkscape:groupmode="layer" id="layer-04-deconv-prediction" inkscape:label="04_Deconvolution_Prediction_Stream">
    <!-- Stream 3A: BayesPrism Deconvolution Engine (x=30, y=635, w=650, h=250) -->
    <g id="stream-03a-bayesprism" transform="translate(30, 635)">
      <rect width="650" height="250" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Section Label Banner -->
      <rect x="0" y="0" width="650" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        BayesPrism Deconvolution Engine (Hierarchical Collinearity Resolution)
      </text>
      <text x="500" y="16" font-size="7.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">b ≈ Φθ,  θ ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan></text>

      <!-- Upper Theory & Formulation Cards -->
      <g transform="translate(14, 34)">
        <!-- Pseudobulk Model Card -->
        <rect width="295" height="66" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Pseudobulk Mixture Model:</text>
        <text x="10" y="32" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
          <tspan font-weight="bold">b</tspan><tspan baseline-shift="sub" font-size="70%">n</tspan> ≈ <tspan font-weight="bold">Φ</tspan> · θ<tspan baseline-shift="sub" font-size="70%">n</tspan>,   where  θ<tspan baseline-shift="sub" font-size="70%">n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan>
        </text>
        <text x="10" y="48" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
          ∑<tspan baseline-shift="sub" font-size="70%">k=1</tspan><tspan baseline-shift="super" font-size="70%">K</tspan> θ<tspan baseline-shift="sub" font-size="70%">k,n</tspan> = 1,   θ<tspan baseline-shift="sub" font-size="70%">k,n</tspan> ≥ 0  (∀ n ∈ {{1,...,N}})
        </text>
        <text x="10" y="60" font-size="6" fill="{COLOR_TEXT_MUTED}">Simplex-constrained physical cell state proportions</text>

        <!-- Null-Space Variance Cancellation Card -->
        <g transform="translate(305, 0)">
          <rect width="315" height="66" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_MCMC}">Null-Space Variance Cancellation:</text>
          <text x="10" y="30" font-size="8.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            <tspan font-weight="bold">M</tspan><tspan font-weight="bold">v</tspan> = <tspan font-weight="bold">0</tspan>  ⟹  <tspan font-weight="bold">v</tspan> ∈ Null(<tspan font-weight="bold">M</tspan>)
          </text>
          <text x="10" y="44" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
            Var(θ̂<tspan baseline-shift="sub" font-size="70%">n,t</tspan><tspan baseline-shift="super" font-size="70%">(0)</tspan>) = ∑ Var(θ̂<tspan baseline-shift="sub" font-size="70%">s</tspan>) + 2∑ Cov(θ̂<tspan baseline-shift="sub" font-size="70%">s</tspan>, θ̂<tspan baseline-shift="sub" font-size="70%">s'</tspan>) = 𝒪(1)
          </text>
          <text x="10" y="58" font-size="6" fill="{COLOR_TEXT_MUTED}">Cancels collinear variance explosion (Gram condition κ(Φ) → ∞)</text>
        </g>
      </g>

      <!-- Middle Flow: Hierarchy -> MCMC Gibbs -> Output Fractions -->
      <g transform="translate(14, 110)">
        <!-- Two-Tier Hierarchy -->
        <g transform="translate(0, 0)">
          <rect width="180" height="82" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Two-Tier Partition Mapping:</text>
          <rect x="10" y="22" width="160" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="90" y="36" font-size="7" font-weight="600" fill="{OKABE_BLUE}" text-anchor="middle">
            Broad Types 𝒯  (T &lt; S)
          </text>

          <path d="M 90 44 L 90 52" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

          <rect x="10" y="52" width="160" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="90" y="66" font-size="7" font-weight="600" fill="{COLOR_MCMC}" text-anchor="middle">
            Fine States 𝒮  (S states)
          </text>
          <text x="10" y="80" font-size="5" fill="{COLOR_TEXT_MUTED}">Aggregation operator M ∈ {{0, 1}}<tspan baseline-shift="super" font-size="70%">T × S</tspan></text>
        </g>

        <!-- Arrow Hierarchy -> Gibbs -->
        <path d="M 186 41 L 202 41" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- MCMC Gibbs Sampling Card -->
        <g transform="translate(208, 0)">
          <rect width="210" height="82" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_MCMC}">MCMC Gibbs Sampling Cycle:</text>
          
          <text x="10" y="30" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            1. Latent counts: Z<tspan baseline-shift="sub" font-size="70%">n,g,s</tspan> ~ Multinomial( X<tspan baseline-shift="sub" font-size="70%">n,g</tspan>, p )
          </text>
          <text x="10" y="44" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            2. State fractions: θ<tspan baseline-shift="sub" font-size="70%">n,s</tspan> ~ Dirichlet( α + Z )
          </text>
          <text x="10" y="58" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            3. Expression profile refinement: ψ<tspan baseline-shift="sub" font-size="70%">t</tspan>
          </text>
          <text x="10" y="74" font-size="5.8" fill="{COLOR_TEXT_MUTED}">
            Shrinkage prior: γ<tspan baseline-shift="sub" font-size="70%">t,g</tspan> ~ 𝒩( 0, σ<tspan baseline-shift="super" font-size="70%">2</tspan> ) across 2,000 steps
          </text>
        </g>

        <!-- Arrow Gibbs -> Output -->
        <path d="M 419 41 L 435 41" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

        <!-- Output Cell Fractions Card -->
        <g transform="translate(440, 0)">
          <rect width="180" height="82" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_RESPONDER}">Estimated Cell Fractions:</text>
          
          <text x="10" y="30" font-size="7.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            θ̂<tspan baseline-shift="sub" font-size="70%">n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan>
          </text>

          <!-- Stacked Proportion Bar -->
          <g transform="translate(10, 38)">
            <rect x="0" y="0" width="45" height="14" fill="{OKABE_BLUE}" />
            <rect x="45" y="0" width="35" height="14" fill="{OKABE_SKY_BLUE}" />
            <rect x="80" y="0" width="25" height="14" fill="{COLOR_RESPONDER}" />
            <rect x="105" y="0" width="35" height="14" fill="{COLOR_MELANOMA}" />
            <rect x="140" y="0" width="20" height="14" fill="{COLOR_NON_RESPONDER}" />
          </g>
          <text x="10" y="66" font-size="5.8" fill="{COLOR_TEXT_MUTED}">CD8 | CD4 | B | Myeloid | Tumor</text>
          <text x="10" y="76" font-size="6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Patient-specific cellular composition</text>
        </g>
      </g>

      <!-- Bottom Benchmark Footer Badges -->
      <g transform="translate(14, 206)">
        <text x="0" y="10" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Benchmarked Deconvolution Suite:</text>
        <g transform="translate(125, 0)">
          <rect x="0" y="0" width="68" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="34" y="9.5" font-size="5.8" font-weight="600" fill="{OKABE_BLUE}" text-anchor="middle">BayesPrism (Gibbs)</text>

          <rect x="74" y="0" width="65" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="106" y="9.5" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">InstaPrism (EM)</text>

          <rect x="145" y="0" width="72" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="181" y="9.5" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">RegDeconv (FISTA)</text>

          <rect x="223" y="0" width="58" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="252" y="9.5" font-size="5.8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">CIBERSORTx</text>
        </g>
      </g>
    </g>

    <!-- Stream 3B: Stability Selection & Machine Learning Classifier (x=700, y=635, w=670, h=250) -->
    <g id="stream-03b-stability-selection" transform="translate(700, 635)">
      <rect width="670" height="250" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Section Label Banner -->
      <rect x="0" y="0" width="670" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Stability Selection &amp; Random Forest Response Classification
      </text>
      <text x="490" y="16" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Finite-Sample Error Controlled Biomarker Discovery</text>

      <!-- Upper Tier: SS-CPSS Theory & Elaborated Random Forest Ensemble -->
      <g transform="translate(14, 34)">
        <!-- Stability Selection Formulation Card -->
        <g transform="translate(0, 0)">
          <rect width="360" height="82" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Complementary Pairs Stability Selection (SS-CPSS):
          </text>
          
          <text x="10" y="30" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
            Π<tspan baseline-shift="sub" font-size="70%">k</tspan>(λ) = 1/(2B) ∑<tspan baseline-shift="sub" font-size="65%">b=1</tspan><tspan baseline-shift="super" font-size="65%">2B</tspan> 𝕀( k ∈ Ŝ<tspan baseline-shift="super" font-size="70%">(b)</tspan>(λ) )
          </text>
          <text x="10" y="44" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Ŝ<tspan baseline-shift="sub" font-size="70%">stable</tspan> = {{ k : max<tspan baseline-shift="sub" font-size="65%">λ</tspan> Π<tspan baseline-shift="sub" font-size="70%">k</tspan>(λ) ≥ π<tspan baseline-shift="sub" font-size="70%">thr</tspan> }}  (π<tspan baseline-shift="sub" font-size="70%">thr</tspan> = 0.60)
          </text>

          <text x="10" y="58" font-size="6.8" font-weight="600" fill="{COLOR_NON_RESPONDER}">
            Shah-Samworth Finite-Sample PFER Bound:
          </text>
          <text x="10" y="72" font-size="7.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            𝔼[ V ] ≤ 1/(2π<tspan baseline-shift="sub" font-size="70%">thr</tspan> - 1) · q<tspan baseline-shift="sub" font-size="70%">Λ</tspan><tspan baseline-shift="super" font-size="70%">2</tspan> / p
          </text>
        </g>

        <!-- Elaborated Multi-Depth Random Forest Ensemble -->
        <g transform="translate(372, 0)">
          <rect width="270" height="82" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Random Forest Ensemble (Decision Trees):</text>

          <!-- Tree 1 -->
          <g transform="translate(18, 24)">
            <!-- Root -->
            <circle cx="20" cy="4" r="2.5" fill="{OKABE_BLUE}" />
            <line x1="20" y1="4" x2="10" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <line x1="20" y1="4" x2="30" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <!-- Depth 1 -->
            <circle cx="10" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <circle cx="30" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <line x1="10" y1="15" x2="5" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="10" y1="15" x2="15" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="25" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="35" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <!-- Leaves (R vs NR) -->
            <rect x="3.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <rect x="13.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <rect x="23.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <rect x="33.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <text x="20" y="35" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree 1</text>
          </g>

          <!-- Tree 2 -->
          <g transform="translate(70, 24)">
            <circle cx="20" cy="4" r="2.5" fill="{OKABE_BLUE}" />
            <line x1="20" y1="4" x2="10" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <line x1="20" y1="4" x2="30" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <circle cx="10" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <circle cx="30" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <line x1="10" y1="15" x2="5" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="10" y1="15" x2="15" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="25" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="35" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <rect x="3.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <rect x="13.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <rect x="23.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <rect x="33.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <text x="20" y="35" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree 2</text>
          </g>

          <!-- Tree B -->
          <g transform="translate(122, 24)">
            <circle cx="20" cy="4" r="2.5" fill="{OKABE_BLUE}" />
            <line x1="20" y1="4" x2="10" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <line x1="20" y1="4" x2="30" y2="15" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
            <circle cx="10" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <circle cx="30" cy="15" r="2.0" fill="{OKABE_BLUE}" />
            <line x1="10" y1="15" x2="5" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="10" y1="15" x2="15" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="25" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <line x1="30" y1="15" x2="35" y2="25" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.6" />
            <rect x="3.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <rect x="13.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <rect x="23.5" y="24" width="3" height="3" fill="{COLOR_RESPONDER}" />
            <rect x="33.5" y="24" width="3" height="3" fill="{COLOR_NON_RESPONDER}" />
            <text x="20" y="35" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree B</text>
          </g>

          <!-- Voting Aggregator -->
          <g transform="translate(170, 20)">
            <text x="0" y="8" font-size="5.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">P̂( Y=1 | θ̂ ) = 1/B ∑ T<tspan baseline-shift="sub" font-size="70%">b</tspan>(θ̂)</text>
            <rect x="0" y="14" width="46" height="13" fill="{COLOR_RESPONDER}" />
            <text x="23" y="23" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Responder (R)</text>

            <rect x="50" y="14" width="46" height="13" fill="{COLOR_NON_RESPONDER}" />
            <text x="73" y="23" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Non-Resp (NR)</text>
          </g>

          <text x="135" y="74" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Ensemble voting on stable fractions</text>
        </g>
      </g>

      <!-- Lower Tier: ROC Performance & Top Biomarkers -->
      <g transform="translate(14, 126)">
        <!-- ROC Curve Mini-Canvas -->
        <g transform="translate(0, 0)">
          <rect width="84" height="74" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <line x1="12" y1="62" x2="74" y2="62" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.5" />
          <line x1="12" y1="62" x2="12" y2="10" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.5" />
          
          <!-- Diagonal Chance Line -->
          <line x1="12" y1="62" x2="74" y2="10" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1.5,1.5" />
          
          <!-- Classifier ROC Curve (AUC = 0.812) -->
          <path d="M 12 62 C 16 32, 28 16, 74 10" fill="none" stroke="{OKABE_BLUE}" stroke-width="1.2" />
          
          <text x="43" y="70" font-size="4.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">1 - Specificity</text>
          <text x="6" y="36" font-size="4.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90, 6, 36)">Sensitivity</text>
        </g>

        <!-- Performance Summary & Top Biomarkers -->
        <g transform="translate(94, 0)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Cross-Cohort Discrimination Benchmark (9 ICB Trials):
          </text>
          <text x="0" y="24" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">
            Pooled ROC-AUC = 0.812 (95% CI: 0.76 - 0.86, Hanley-McNeil test)
          </text>

          <g transform="translate(0, 36)">
            <text x="0" y="0" font-size="6" fill="{COLOR_TEXT_SECONDARY}">
              • <tspan font-weight="700" fill="{COLOR_RESPONDER}">CD8+ Effector Memory T cells (Tem)</tspan>: Selected in 100% of CPSS subsamples (OR = 2.45, p &lt; 0.001)
            </text>
            <text x="0" y="12" font-size="6" fill="{COLOR_TEXT_SECONDARY}">
              • <tspan font-weight="700" fill="{COLOR_RESPONDER}">TCF7+ Stem-like Progenitor Exhausted</tspan>: Strong positive predictor of durable response
            </text>
            <text x="0" y="24" font-size="6" fill="{COLOR_TEXT_SECONDARY}">
              • <tspan font-weight="700" fill="{COLOR_NON_RESPONDER}">CD163+ / MRC1+ M2 Macrophages</tspan>: Stable negative predictor of response (OR = 0.42, p &lt; 0.005)
            </text>
          </g>
        </g>
      </g>
    </g>
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 05: Single-Cell Response Ground Truth & Triangulation Stream (y = 910 to 1326, h = 416)
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 05: TRIANGULATION CONCORDANCE STREAM -->
  <g inkscape:groupmode="layer" id="layer-05-triangulation" inkscape:label="05_Triangulation_Concordance_Stream">
    <g id="stream-04-triangulation" transform="translate(30, 910)">
      <rect width="1340" height="416" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Section Label Banner -->
      <rect x="0" y="0" width="1340" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Single-Cell Response Ground Truth &amp; In Silico Triangulation Concordance
      </text>
      <text x="890" y="16" font-size="7.5" fill="{COLOR_TEXT_MUTED}">
        Cross-Checking Pseudobulk Deconvolution vs. Milo Single-Cell Neighborhood Testing
      </text>

      <!-- Sub-Branch 4A: Milo Continuous Graph Testing (x=14, y=34, w=440, h=370) -->
      <g transform="translate(14, 34)">
        <rect width="440" height="370" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <rect width="440" height="20" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Branch A: Single-Cell Milo Continuous Differential Abundance (Ground Truth)
        </text>

        <!-- Milo Formulation Box -->
        <g transform="translate(12, 28)">
          <rect width="416" height="74" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Negative Binomial Neighborhood GLM:</text>
          <text x="10" y="30" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
            log( 𝔼[ N<tspan baseline-shift="sub" font-size="70%">v,j</tspan> ] ) = β<tspan baseline-shift="sub" font-size="70%">0,v</tspan> + β<tspan baseline-shift="sub" font-size="70%">v</tspan><tspan baseline-shift="super" font-size="70%">milo</tspan> · Y<tspan baseline-shift="sub" font-size="70%">j</tspan> + log( TotalCells<tspan baseline-shift="sub" font-size="70%">j</tspan> )
          </text>
          <text x="10" y="44" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
            Spatial-FDR(v) = ∑<tspan baseline-shift="sub" font-size="65%">u</tspan> w<tspan baseline-shift="sub" font-size="65%">uv</tspan> 𝕀(p<tspan baseline-shift="sub" font-size="65%">u</tspan> ≤ p<tspan baseline-shift="sub" font-size="65%">v</tspan>) / ∑<tspan baseline-shift="sub" font-size="65%">u</tspan> w<tspan baseline-shift="sub" font-size="65%">uv</tspan>  ·  p<tspan baseline-shift="sub" font-size="65%">v</tspan> / F̂(p<tspan baseline-shift="sub" font-size="65%">v</tspan>)
          </text>
          <text x="10" y="56" font-size="5.8" fill="{COLOR_TEXT_MUTED}">
            w_uv: Jaccard overlap index between overlapping neighborhoods 𝒩(u) and 𝒩(v)
          </text>
          <text x="10" y="68" font-size="6" font-weight="600" fill="{COLOR_RESPONDER}">
            Attributes response effect β_milo directly to individual single cells without clustering bias
          </text>
        </g>

        <!-- Continuous Graph Neighborhoods Illustration -->
        <g transform="translate(12, 114)">
          <rect width="416" height="154" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Single-Cell Manifold kNN Topology with Differential Neighborhoods:
          </text>

          <g transform="translate(30, 20)">
            <line x1="40" y1="35" x2="110" y2="25" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="40" y1="35" x2="70" y2="85" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="110" y1="25" x2="190" y2="40" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="70" y1="85" x2="150" y2="95" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="110" y1="25" x2="150" y2="95" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="190" y1="40" x2="260" y2="35" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="190" y1="40" x2="230" y2="85" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="260" y1="35" x2="330" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="230" y1="85" x2="330" y2="70" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
            <line x1="150" y1="95" x2="230" y2="85" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

            <circle cx="70" cy="50" r="38" fill="{COLOR_RESPONDER}" opacity="0.14" stroke="{COLOR_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />
            <circle cx="170" cy="60" r="38" fill="{COLOR_RESPONDER}" opacity="0.10" stroke="{COLOR_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />
            <circle cx="270" cy="60" r="38" fill="{COLOR_NON_RESPONDER}" opacity="0.14" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />

            <circle cx="40" cy="35" r="4" fill="{COLOR_RESPONDER}" />
            <circle cx="70" cy="85" r="5" fill="{COLOR_RESPONDER}" stroke="#0F172A" stroke-width="0.8" />
            <circle cx="110" cy="25" r="4" fill="{COLOR_RESPONDER}" />
            <circle cx="150" cy="95" r="4" fill="{COLOR_RESPONDER}" />
            <circle cx="190" cy="40" r="5" fill="{COLOR_RESPONDER}" stroke="#0F172A" stroke-width="0.8" />
            <circle cx="230" cy="85" r="4" fill="{COLOR_NON_RESPONDER}" />
            <circle cx="260" cy="35" r="5" fill="{COLOR_NON_RESPONDER}" stroke="#0F172A" stroke-width="0.8" />
            <circle cx="330" cy="70" r="4" fill="{COLOR_NON_RESPONDER}" />

            <text x="70" y="105" font-size="6" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">
              𝒩(v<tspan baseline-shift="sub" font-size="70%">R1</tspan>): CD8+ Tem (β<tspan baseline-shift="super" font-size="70%">milo</tspan> &gt; 0)
            </text>
            <text x="270" y="105" font-size="6" font-weight="700" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">
              𝒩(v<tspan baseline-shift="sub" font-size="70%">NR</tspan>): Term Tex (β<tspan baseline-shift="super" font-size="70%">milo</tspan> &lt; 0)
            </text>
          </g>
        </g>

        <!-- Bottom Attribution Dossier -->
        <g transform="translate(12, 276)">
          <rect width="416" height="84" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Single-Cell Ground Truth Attribution Summary:
          </text>
          <text x="10" y="27" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
            • <tspan font-weight="600" fill="{COLOR_RESPONDER}">Effector Memory (Tem)</tspan>: β<tspan baseline-shift="super" font-size="70%">milo</tspan> = +1.82 (Spatial FDR &lt; 0.005, highly responder-enriched)
          </text>
          <text x="10" y="40" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
            • <tspan font-weight="600" fill="{COLOR_RESPONDER}">TCF7+ Progenitor Tex</tspan>: β<tspan baseline-shift="super" font-size="70%">milo</tspan> = +1.34 (Spatial FDR &lt; 0.010, stem-like reserve)
          </text>
          <text x="10" y="53" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
            • <tspan font-weight="600" fill="{COLOR_NON_RESPONDER}">Terminally Exhausted Tex</tspan>: β<tspan baseline-shift="super" font-size="70%">milo</tspan> = -1.65 (Spatial FDR &lt; 0.005, dysfunction)
          </text>
          <text x="10" y="66" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
            • <tspan font-weight="600" fill="{COLOR_NON_RESPONDER}">CD163+ M2 Macrophages</tspan>: β<tspan baseline-shift="super" font-size="70%">milo</tspan> = -1.48 (Spatial FDR &lt; 0.010, immunosuppressive)
          </text>
        </g>
      </g>

      <!-- Sub-Branch 4B: In Silico Pseudo-bulk & Deconv (x=468, y=34, w=380, h=370) -->
      <g transform="translate(468, 34)">
        <rect width="380" height="370" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <rect width="380" height="20" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Branch B: In Silico Pseudo-bulk Aggregation &amp; Deconvolution
        </text>

        <!-- Pooling Formula Card -->
        <g transform="translate(12, 28)">
          <rect width="356" height="78" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Per-Patient Single-Cell Pooling:</text>
          <text x="10" y="28" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
            <tspan font-weight="bold">b</tspan><tspan baseline-shift="sub" font-size="70%">j</tspan><tspan baseline-shift="super" font-size="70%">pseudo</tspan> = ∑<tspan baseline-shift="sub" font-size="65%">c ∈ Patient<tspan baseline-shift="sub" font-size="65%">j</tspan></tspan> <tspan font-weight="bold">x</tspan><tspan baseline-shift="sub" font-size="70%">c</tspan>
          </text>
          <text x="10" y="42" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">
            Preserves physical transcript yield &amp; cellular mass proportions per patient
          </text>
          <text x="10" y="55" font-size="6.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Deconvolve with reference <tspan font-weight="bold">Φ</tspan>  ⟹  θ̂<tspan baseline-shift="sub" font-size="70%">j</tspan><tspan baseline-shift="super" font-size="70%">pseudo</tspan> ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan>
          </text>
          <text x="10" y="68" font-size="5.8" fill="{COLOR_TEXT_MUTED}">Eliminates batch discrepancy between reference and mixture</text>
        </g>

        <!-- Logistic Regression on Pseudo-bulk Fractions -->
        <g transform="translate(12, 114)">
          <rect width="356" height="154" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Logistic Regression on Deconvoluted Fractions:
          </text>

          <text x="10" y="32" font-size="7.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            logit( P( Y<tspan baseline-shift="sub" font-size="70%">j</tspan>=1 ) ) = β<tspan baseline-shift="sub" font-size="70%">0</tspan> + ∑<tspan baseline-shift="sub" font-size="65%">k</tspan> β<tspan baseline-shift="sub" font-size="70%">k</tspan><tspan baseline-shift="super" font-size="70%">deconv</tspan> · θ̂<tspan baseline-shift="sub" font-size="70%">k,j</tspan><tspan baseline-shift="super" font-size="70%">pseudo</tspan>
          </text>

          <!-- Deconvolution Fractions Table -->
          <g transform="translate(10, 42)">
            <rect width="336" height="100" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="10" y="13" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
              Estimated Coefficients (β_deconv):
            </text>

            <text x="10" y="30" font-size="6" fill="{COLOR_TEXT_SECONDARY}">• CD8+ Effector Memory:</text>
            <text x="180" y="30" font-size="6" font-weight="700" fill="{COLOR_RESPONDER}">β<tspan baseline-shift="super" font-size="70%">deconv</tspan> = +2.15 (p &lt; 0.001)</text>

            <text x="10" y="46" font-size="6" fill="{COLOR_TEXT_SECONDARY}">• TCF7+ Progenitor Tex:</text>
            <text x="180" y="46" font-size="6" font-weight="700" fill="{COLOR_RESPONDER}">β<tspan baseline-shift="super" font-size="70%">deconv</tspan> = +1.42 (p = 0.008)</text>

            <text x="10" y="62" font-size="6" fill="{COLOR_TEXT_SECONDARY}">• CD163+ M2 Macrophages:</text>
            <text x="180" y="62" font-size="6" font-weight="700" fill="{COLOR_NON_RESPONDER}">β<tspan baseline-shift="super" font-size="70%">deconv</tspan> = -1.68 (p = 0.003)</text>

            <text x="10" y="78" font-size="6" fill="{COLOR_TEXT_SECONDARY}">• Terminally Exhausted Tex:</text>
            <text x="180" y="78" font-size="6" font-weight="700" fill="{COLOR_NON_RESPONDER}">β<tspan baseline-shift="super" font-size="70%">deconv</tspan> = -1.52 (p = 0.006)</text>

            <text x="10" y="93" font-size="5.2" fill="{COLOR_TEXT_MUTED}">Simplex regularized estimation with adaptive condition deficit</text>
          </g>
        </g>

        <!-- Self-Consistency Check Metric -->
        <g transform="translate(12, 276)">
          <rect width="356" height="84" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">In Silico Self-Consistency Metric:</text>
          <text x="10" y="27" font-size="6.2" font-weight="600" fill="{OKABE_BLUE}">
            • Reconstruction Fidelity: Pearson r = 0.94 (true vs. deconvoluted)
          </text>
          <text x="10" y="40" font-size="6" fill="{COLOR_TEXT_SECONDARY}">
            • Proves that deconvolution faithfully separates known mixtures
          </text>
          <text x="10" y="53" font-size="6" fill="{COLOR_TEXT_SECONDARY}">
            • Directly transferable to empirical patient cohorts
          </text>
          <text x="10" y="68" font-size="5.5" fill="{COLOR_TEXT_MUTED}">
            scripts/concordance_calibration/02_pseudobulk_self_deconvolution.py
          </text>
        </g>
      </g>

      <!-- Sub-Branch 4C: Triangulation Concordance Scatter (x=862, y=34, w=464, h=370) -->
      <g transform="translate(862, 34)">
        <rect width="464" height="370" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <rect width="464" height="20" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Branch C: Triangulation Concordance Scatter (β_deconv vs. β_milo)
        </text>

        <!-- 4-Quadrant Scatter Plot Canvas -->
        <g transform="translate(24, 30)">
          <rect width="416" height="240" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

          <!-- Axes & Quadrants -->
          <g transform="translate(48, 16)">
            <!-- Axes Lines -->
            <line x1="160" y1="10" x2="160" y2="190" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" /> <!-- Y-axis -->
            <line x1="10" y1="100" x2="310" y2="100" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" /> <!-- X-axis -->

            <!-- Quadrant 1: Top Right (Concordant Responder) -->
            <rect x="160" y="10" width="150" height="90" fill="{COLOR_RESPONDER}" opacity="0.08" />
            <text x="235" y="25" font-size="6.8" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">
              Q1: Concordant Responder
            </text>

            <!-- Quadrant 3: Bottom Left (Concordant Non-Responder) -->
            <rect x="10" y="100" width="150" height="90" fill="{COLOR_NON_RESPONDER}" opacity="0.08" />
            <text x="85" y="180" font-size="6.8" font-weight="700" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">
              Q3: Concordant Non-Responder
            </text>

            <!-- Quadrants 2 & 4: Discordant Off-Diagonal -->
            <text x="85" y="25" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Q2: Discordant (SC-, Bulk+)</text>
            <text x="235" y="180" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Q4: Discordant (SC+, Bulk-)</text>

            <!-- Trend Line -->
            <line x1="30" y1="174" x2="290" y2="26" stroke="{OKABE_BLUE}" stroke-width="1.2" stroke-dasharray="3,2" />

            <!-- Scatter Data Points -->
            <!-- Q1: Concordant Responders -->
            <circle cx="265" cy="40" r="4.5" fill="{OKABE_BLUE}" /> <!-- CD8 Tem -->
            <text x="273" y="38" font-size="6.5" font-weight="700" fill="{OKABE_BLUE}">CD8+ Tem</text>

            <circle cx="225" cy="55" r="4.0" fill="{COLOR_RESPONDER}" /> <!-- CXCL9+ Mac -->
            <text x="200" y="50" font-size="6" font-style="italic" fill="{COLOR_RESPONDER}">CXCL9+ Mac</text>

            <circle cx="240" cy="65" r="4.0" fill="{OKABE_SKY_BLUE}" /> <!-- TCF7+ Tex -->
            <text x="247" y="69" font-size="6" font-style="italic" fill="{OKABE_SKY_BLUE}">TCF7+ Tex</text>

            <!-- Q3: Concordant Non-Responders -->
            <circle cx="65" cy="145" r="4.0" fill="{COLOR_NON_RESPONDER}" /> <!-- CD163+ M2 -->
            <text x="40" y="152" font-size="6" font-style="italic" fill="{COLOR_NON_RESPONDER}">CD163+ M2</text>

            <circle cx="95" cy="130" r="4.0" fill="{COLOR_MELANOMA}" /> <!-- Term Tex -->
            <text x="80" y="122" font-size="6" font-style="italic" fill="{COLOR_MELANOMA}">Term Tex</text>

            <circle cx="130" cy="165" r="4.0" fill="{OKABE_ORANGE}" /> <!-- Malignant -->
            <text x="140" y="168" font-size="6" font-style="italic" fill="{OKABE_ORANGE}">Malignant</text>

            <!-- Axis Labels -->
            <text x="160" y="202" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
              Single-Cell Milo Response Effect (β<tspan baseline-shift="super" font-size="70%">milo</tspan>) ⟶
            </text>
            <text x="-100" y="-12" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" transform="rotate(-90)">
              Deconvolution Response Effect (β<tspan baseline-shift="super" font-size="70%">deconv</tspan>) ⟶
            </text>
          </g>
        </g>

        <!-- Triangulation Summary Dossier -->
        <g transform="translate(24, 276)">
          <rect width="416" height="84" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Triangulation Concordance Summary:
          </text>
          <text x="10" y="28" font-size="6.2" font-weight="600" fill="{OKABE_BLUE}">
            • Rank Correlation: Spearman ρ = 0.84 (p = 2.4 × 10⁻⁵, strong concordance)
          </text>
          <text x="10" y="42" font-size="6.2" font-weight="600" fill="{COLOR_RESPONDER}">
            • Directional Agreement: 87.5% in concordant quadrants (Q1 &amp; Q3)
          </text>
          <text x="10" y="56" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">
            • Discrepancy Attribution: Resolves collinear leakage &amp; cell size mRNA asymmetry (S_k)
          </text>
          <text x="10" y="70" font-size="5.5" fill="{COLOR_TEXT_MUTED}">
            Validates that deconvolution reconstructs genuine single-cell biology across ICB trials
          </text>
        </g>
      </g>
    </g>
  </g>
""")

    # --------------------------------------------------------------------------
    # LAYER 06: Flow Conduits & Connectors
    # RESTRUCTURED: Bulk Matrix B drops straight down into BayesPrism!
    # Database drops into Pipeline; Pipeline feeds BayesPrism and Milo/Pseudo-bulk!
    # --------------------------------------------------------------------------
    parts.append(f"""
  <!-- LAYER 06: FLOW CONDUITS -->
  <g inkscape:groupmode="layer" id="layer-06-conduits" inkscape:label="06_Flow_Conduits">
    <!-- Stream 1B (Database) to Stream 2 (Pipeline): direct vertical drop in horizontal gutter (y=230 to 246) -->
    <path d="M 965 230 L 965 246" stroke="{OKABE_BLUE}" stroke-width="1.2" marker-end="url(#arr-blue)" />

    <!-- Stream 1A (Bulk Matrix B) to Stream 3A (BayesPrism): DIRECT VERTICAL DROP at x=150 (y=616 to 635)! -->
    <!-- Zero margin snaking needed now that Left Column sits directly above BayesPrism! -->
    <path d="M 150 616 L 150 635" stroke="{OKABE_BLUE}" stroke-width="1.2" marker-end="url(#arr-blue)" />
    <rect x="130" y="621" width="40" height="10" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="150" y="628.5" font-size="5.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">Bulk B</text>

    <!-- Stream 2 (Reference Phi) down to Stream 3A (BayesPrism): routes in horizontal gutter (y=626) with bridge gap over x=690 -->
    <path d="M 1218 616 L 1218 626 L 696 626 M 684 626 L 420 626 L 420 635" 
          stroke="{OKABE_BLUE}" stroke-width="1.2" fill="none" marker-end="url(#arr-blue)" />
    <rect x="400" y="621" width="40" height="10" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="420" y="628.5" font-size="5.5" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">Ref Φ</text>

    <!-- Stream 3A (Deconv Fractions) across to Stream 3B (Stability Selection): horizontal gutter (x=680 to 700) -->
    <path d="M 680 755 L 700 755" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" marker-end="url(#arr-slate)" />
    <rect x="677" y="743" width="26" height="9" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="690" y="749.5" font-size="5.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">θ̂</text>

    <!-- Stream 2 (scRNA-seq Pipeline Step 6) down to Stream 4 (Triangulation): drops through vertical gutter at x=690, splits at y=898 -->
    <path d="M 952 616 L 952 622 L 690 622 L 690 898 M 690 898 L 245 898 L 245 910 M 690 898 L 658 898 L 658 910" 
          stroke="{OKABE_BLUE}" stroke-width="1.2" fill="none" marker-end="url(#arr-blue)" />
    <rect x="625" y="892" width="130" height="12" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="690" y="900.5" font-size="5.2" font-weight="600" fill="{OKABE_BLUE}" text-anchor="middle">scRNA-seq + ICB Response Labels</text>

    <!-- Triangulation Branch B (Pseudo-bulk) across to Branch C (Scatter Plot): gutter at x=878 to 892 -->
    <path d="M 878 1095 L 891 1095" stroke="{COLOR_RESPONDER}" stroke-width="1.0" marker-end="url(#arr-green)" />
    <rect x="874" y="1083" width="22" height="9" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="885" y="1089.5" font-size="5.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}" text-anchor="middle">β̂<tspan baseline-shift="super" font-size="70%">deconv</tspan></text>

    <!-- Triangulation Branch A (Milo ground truth) to Branch C: routes via bottom gutter (y=1322) -->
    <path d="M 245 1314 L 245 1322 L 1094 1322 L 1094 1314" 
          stroke="{OKABE_BLUE}" stroke-width="1.0" fill="none" marker-end="url(#arr-blue)" />
    <rect x="605" y="1317" width="170" height="11" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="690" y="1325" font-size="5.2" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">Single-Cell Milo Response Effect β̂<tspan baseline-shift="super" font-size="70%">milo</tspan></text>
  </g>
</svg>
""")

    return "".join(parts)


def main() -> None:
    output_dir = Path("assets")
    output_dir.mkdir(parents=True, exist_ok=True)

    master_svg_path = output_dir / "graphical_abstract.svg"
    svg_content = render_master_graphical_abstract()
    master_svg_path.write_text(svg_content, encoding="utf-8")
    print(f"Master Graphical Abstract successfully generated at: {master_svg_path}")


if __name__ == "__main__":
    main()
