#!/usr/bin/env python3
"""
Generate 8 Modular Publication-Grade Vector SVG Assets for Graphical Abstract (Elaborated Version).

Follows:
- Nature Methods minimalist wireframe standard (pure white background, 0.75 pt hairlines, rx=0, zero drop shadows).
- Okabe-Ito Colorblind-Safe scientific palette strictly on data marks.
- Formal LaTeX math symbols and serif typography hierarchy (Times New Roman / Georgia).
- Native Inkscape layers and semantic groups.
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


# ==============================================================================
# ASSET 1: Bulk ICB Cohorts (Elaborated)
# ==============================================================================
def generate_asset_01_bulk_icb_cohorts() -> str:
    """Generate Asset 1: Bulk RNA-seq ICB Patient Cohorts (full left column stream)."""
    w, h = 510, 538
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-01-bulk-icb"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-01-bulk-icb" inkscape:label="01_Bulk_ICB_Cohorts">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Banner -->
    <rect x="0" y="0" width="{w}" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
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

      <!-- Response Status Legend: Placed in Row 2! -->
      <g transform="translate(355, 38)">
        <rect width="125" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <circle cx="10" cy="10" r="3" fill="{COLOR_RESPONDER}" />
        <text x="18" y="12.5" font-size="5.5" font-weight="600" fill="{COLOR_RESPONDER}">Responder (CR/PR)</text>
        <circle cx="10" cy="21" r="3" fill="{COLOR_NON_RESPONDER}" />
        <text x="18" y="23.5" font-size="5.5" font-weight="600" fill="{COLOR_NON_RESPONDER}">Non-Resp (SD/PD)</text>
      </g>
    </g>

    <!-- Dividing Rule -->
    <line x1="14" y1="114" x2="{w-14}" y2="114" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <!-- Block 2: Biopsy Specimen & 8-Row Bulk Matrix Representation (y = 120 to 270) -->
    <g transform="translate(14, 120)">
      <!-- Biopsy Tube -->
      <rect x="0" y="8" width="18" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <rect x="3" y="22" width="12" height="25" fill="{COLOR_MELANOMA}" opacity="0.4" />
      <text x="9" y="58" font-size="5.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Biopsy</text>

      <path d="M 22 28 L 30 28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

      <!-- 8-Row Bulk Expression Matrix B -->
      <g transform="translate(34, 0)">
        <rect width="144" height="150" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="72" y="14" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Bulk Matrix <tspan font-family="{FONT_SERIF_MATH}" font-weight="bold">B</tspan> ∈ ℝ<tspan baseline-shift="super" font-size="70%">G × N</tspan>
        </text>
        
        <!-- 8 Rows of Matrix Cells -->
        <g transform="translate(10, 20)">
          <!-- Row 1: CD8A -->
          <rect x="0" y="0" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.9" />
          <rect x="18" y="0" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.3" />
          <rect x="36" y="0" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="54" y="0" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.4" />
          <rect x="72" y="0" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.1" />
          <text x="96" y="7" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD8A</text>

          <!-- Row 2: IFNG -->
          <rect x="0" y="11" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.75" />
          <rect x="18" y="11" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.2" />
          <rect x="36" y="11" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.8" />
          <rect x="54" y="11" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.3" />
          <rect x="72" y="11" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.15" />
          <text x="96" y="18" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">IFNG</text>

          <!-- Row 3: CXCL9 -->
          <rect x="0" y="22" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.95" />
          <rect x="18" y="22" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.1" />
          <rect x="36" y="22" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.9" />
          <rect x="54" y="22" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.5" />
          <rect x="72" y="22" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.05" />
          <text x="96" y="29" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CXCL9</text>

          <!-- Row 4: PDCD1 -->
          <rect x="0" y="33" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.6" />
          <rect x="18" y="33" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.8" />
          <rect x="36" y="33" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.4" />
          <rect x="54" y="33" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.7" />
          <rect x="72" y="33" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.2" />
          <text x="96" y="40" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">PDCD1</text>

          <!-- Row 5: CTLA4 -->
          <rect x="0" y="44" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.5" />
          <rect x="18" y="44" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.6" />
          <rect x="36" y="44" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.3" />
          <rect x="54" y="44" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.65" />
          <rect x="72" y="44" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.3" />
          <text x="96" y="51" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CTLA4</text>

          <!-- Row 6: TCF7 -->
          <rect x="0" y="55" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="18" y="55" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.1" />
          <rect x="36" y="55" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.7" />
          <rect x="54" y="55" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.2" />
          <rect x="72" y="55" width="16" height="9" fill="{OKABE_BLUE}" opacity="0.1" />
          <text x="96" y="62" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">TCF7</text>

          <!-- Row 7: CD163 -->
          <rect x="0" y="66" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.1" />
          <rect x="18" y="66" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.9" />
          <rect x="36" y="66" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.2" />
          <rect x="54" y="66" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.8" />
          <rect x="72" y="66" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.7" />
          <text x="96" y="73" font-size="5.2" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD163</text>

          <!-- Row 8: MLANA -->
          <rect x="0" y="77" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.2" />
          <rect x="18" y="77" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.85" />
          <rect x="36" y="77" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.1" />
          <rect x="54" y="77" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.7" />
          <rect x="72" y="77" width="16" height="9" fill="{COLOR_NON_RESPONDER}" opacity="0.9" />
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
    <line x1="14" y1="278" x2="{w-14}" y2="278" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

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
      <line x1="0" y1="0" x2="{w-28}" y2="0" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
      <text x="0" y="12" font-size="6.2" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Harmonization Invariant:</text>
      <text x="105" y="12" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">
        Standardized Ensembl 111 (GRCh38) gene mapping with cross-cohort RECIST radiographic outcome harmonization.
      </text>
      <text x="0" y="24" font-size="5.2" fill="{COLOR_TEXT_MUTED}">
        Feeds directly into BayesPrism deconvolution engine below with zero inter-panel cross-talk.
      </text>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 2: Single-Cell Tumor Database (.h5ad) (Elaborated & Cleaned)
# ==============================================================================
def generate_asset_02_sc_database_h5ad() -> str:
    """Generate Asset 2: Solid Tumor scRNA-seq Retrieval & Local Database (.h5ad)."""
    w, h = 810, 152
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-02-sc-database"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-02-sc-database" inkscape:label="02_SC_Database_H5AD">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Banner -->
    <rect x="0" y="0" width="{w}" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
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
      <line x1="0" y1="0" x2="{w-32}" y2="0" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
      <text x="0" y="11" font-size="6.2" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Universal Harmonization Standard:</text>
      <text x="135" y="11" font-size="5.8" fill="{COLOR_TEXT_SECONDARY}">
        Ensembl Release 111 (GRCh38) gene nomenclature with strictly verified RECIST radiographic criteria.
      </text>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 3: Automated Single-Cell Pipeline (Ordered in Two Rows!)
# ==============================================================================
def generate_asset_03_automated_sc_pipeline() -> str:
    """Generate Asset 3: Highly Automated Single-Cell Processing Engine ordered in two rows."""
    w, h = 810, 370
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-03-sc-pipeline"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-03-sc-pipeline" inkscape:label="03_Automated_SC_Pipeline">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Section Label Banner -->
    <rect x="0" y="0" width="{w}" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
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

          <!-- Principal Axis PC1 (v1) -->
          <line x1="20" y1="78" x2="162" y2="22" stroke="{OKABE_BLUE}" stroke-width="1.2" marker-end="url(#arr-blue)" />
          <text x="166" y="24" font-size="5.8" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_BLUE}">v<tspan baseline-shift="sub" font-size="70%">1</tspan></text>

          <!-- Orthogonal Axis PC2 (v2) -->
          <line x1="85" y1="70" x2="65" y2="18" stroke="{OKABE_ORANGE}" stroke-width="1.0" marker-end="url(#arr-slate)" />
          <text x="65" y="14" font-size="5.8" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_ORANGE}">v<tspan baseline-shift="sub" font-size="70%">2</tspan></text>

          <!-- Orthogonal symbol at intersection -->
          <path d="M 82 58 L 86 56 L 88 62" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.6" fill="none" />

          <!-- Data Point Cloud & Orthogonal Projections -->
          <circle cx="45" cy="58" r="2.2" fill="{COLOR_RESPONDER}" />
          <line x1="45" y1="58" x2="52" y2="65" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

          <circle cx="65" cy="42" r="2.2" fill="{COLOR_RESPONDER}" />
          <line x1="65" y1="42" x2="72" y2="57" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

          <circle cx="95" cy="55" r="2.2" fill="{COLOR_RESPONDER}" />
          <line x1="95" y1="55" x2="90" y2="50" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

          <circle cx="115" cy="35" r="2.2" fill="{COLOR_NON_RESPONDER}" />
          <line x1="115" y1="35" x2="120" y2="38" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

          <circle cx="135" cy="26" r="2.2" fill="{COLOR_NON_RESPONDER}" />
          <line x1="135" y1="26" x2="138" y2="31" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" stroke-dasharray="1,1" />

          <circle cx="78" cy="32" r="2.2" fill="{OKABE_ORANGE}" />
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
          <circle cx="65" cy="35" r="3.2" fill="{OKABE_BLUE}" stroke="#0F172A" stroke-width="0.5" />
          <circle cx="50" cy="70" r="2.5" fill="{OKABE_BLUE}" />

          <!-- Vertices Cluster B -->
          <circle cx="90" cy="40" r="3.0" fill="{OKABE_ORANGE}" stroke="#0F172A" stroke-width="0.5" />
          <circle cx="115" cy="25" r="2.2" fill="{OKABE_ORANGE}" />
          <circle cx="95" cy="65" r="2.5" fill="{OKABE_ORANGE}" />
          <circle cx="135" cy="35" r="3.2" fill="{OKABE_ORANGE}" stroke="#0F172A" stroke-width="0.5" />
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

            <circle cx="20" cy="26" r="2.2" fill="{OKABE_ORANGE}" />
            <circle cx="28" cy="22" r="2.2" fill="{OKABE_ORANGE}" />
            <rect x="36" y="32" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
            <rect x="44" y="28" width="4.4" height="4.4" fill="{OKABE_BLUE}" />
            <line x1="24" y1="24" x2="40" y2="34" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="1.5,1.5" marker-end="url(#arr-slate)" />
            <text x="36" y="44" font-size="4.5" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">δ<tspan baseline-shift="sub" font-size="70%">batch</tspan></text>

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

            <circle cx="38" cy="30" r="14" fill="{OKABE_BLUE}" opacity="0.08" stroke="{OKABE_BLUE}" stroke-width="0.5" stroke-dasharray="2,2" />
            <path d="M 34 30 L 42 30 M 38 26 L 38 34" stroke="{OKABE_BLUE}" stroke-width="1.0" />
            <text x="44" y="27" font-size="5" font-family="{FONT_SERIF_MATH}" font-weight="bold" fill="{OKABE_BLUE}">μ<tspan baseline-shift="sub" font-size="70%">k</tspan></text>

            <circle cx="22" cy="25" r="2.0" fill="{OKABE_ORANGE}" />
            <line x1="22" y1="25" x2="35" y2="29" stroke="{COLOR_RESPONDER}" stroke-width="0.8" marker-end="url(#arr-green)" />
            <text x="26" y="20" font-size="4.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">w<tspan baseline-shift="sub" font-size="70%">k,1</tspan></text>

            <rect x="52" y="32" width="4" height="4" fill="{OKABE_BLUE}" />
            <line x1="52" y1="34" x2="41" y2="31" stroke="{COLOR_RESPONDER}" stroke-width="0.8" marker-end="url(#arr-green)" />
            <text x="50" y="44" font-size="4.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">w<tspan baseline-shift="sub" font-size="70%">k,2</tspan></text>

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

            <circle cx="28" cy="27" r="2.0" fill="{OKABE_ORANGE}" />
            <rect x="36" y="24" width="4" height="4" fill="{OKABE_BLUE}" />
            <circle cx="34" cy="33" r="2.0" fill="{OKABE_ORANGE}" />
            <rect x="26" y="34" width="4" height="4" fill="{OKABE_BLUE}" />
            <circle cx="41" cy="29" r="2.0" fill="{OKABE_ORANGE}" />
            <text x="35" y="44" font-size="4.8" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">Integrated T</text>

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
            Refinement Guarantee: e(v, C'\{{v\}}) ≥ γ|v|(|C| - |v|)
          </text>
        </g>

        <!-- 3-Phase Algorithmic Loop Diagram -->
        <g transform="translate(8, 42)">
          <!-- Phase 1: Local Move -->
          <g transform="translate(0, 0)">
            <rect width="74" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="37" y="10" font-size="5.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Phase 1: Move</text>

            <ellipse cx="26" cy="35" rx="18" ry="14" fill="{OKABE_BLUE}" opacity="0.12" stroke="{OKABE_BLUE}" stroke-width="0.5" />
            <circle cx="18" cy="32" r="2.2" fill="{OKABE_BLUE}" />
            <circle cx="28" cy="28" r="2.2" fill="{OKABE_BLUE}" />
            <circle cx="24" cy="42" r="2.2" fill="{OKABE_BLUE}" />

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

            <ellipse cx="24" cy="32" rx="14" ry="12" fill="{COLOR_RESPONDER}" opacity="0.12" stroke="{COLOR_RESPONDER}" stroke-width="0.6" stroke-dasharray="1.5,1.5" />
            <circle cx="18" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
            <circle cx="26" cy="34" r="2.0" fill="{COLOR_RESPONDER}" />
            <text x="24" y="48" font-size="4.5" font-weight="600" fill="{COLOR_RESPONDER}" text-anchor="middle">C'<tspan baseline-shift="sub" font-size="70%">1</tspan></text>

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

            <ellipse cx="24" cy="28" rx="14" ry="10" fill="{OKABE_BLUE}" opacity="0.15" stroke="{OKABE_BLUE}" stroke-width="0.6" />
            <circle cx="20" cy="26" r="1.8" fill="{OKABE_BLUE}" />
            <circle cx="26" cy="30" r="1.8" fill="{OKABE_BLUE}" />
            <text x="24" y="39" font-size="4.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">CD8+ T</text>

            <ellipse cx="50" cy="28" rx="14" ry="10" fill="{OKABE_SKY_BLUE}" opacity="0.2" stroke="{OKABE_SKY_BLUE}" stroke-width="0.6" />
            <circle cx="48" cy="26" r="1.8" fill="{OKABE_SKY_BLUE}" />
            <circle cx="54" cy="30" r="1.8" fill="{OKABE_SKY_BLUE}" />
            <text x="50" y="39" font-size="4.5" font-weight="700" fill="{OKABE_SKY_BLUE}" text-anchor="middle">CD4+ T</text>

            <line x1="37" y1="28" x2="38" y2="48" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />

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
</svg>
"""
    return svg


# ==============================================================================
# ASSET 4: Cell State Profiles (Kept synchronized)
# ==============================================================================
def generate_asset_04_cell_state_profiles() -> str:
    """Generate Asset 4: Cell State Profile Calculation (Centroid Matrix Phi)."""
    return generate_asset_04_cell_state_profiles_impl()


def generate_asset_04_cell_state_profiles_impl() -> str:
    w, h = 380, 200
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-04-cell-profiles"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-04-cell-profiles" inkscape:label="04_Cell_State_Profiles">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <text x="14" y="20" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Cell State Profile Matrix (<tspan font-family="{FONT_SERIF_MATH}" font-weight="bold">Φ</tspan>)
    </text>
    <text x="14" y="33" font-size="8" fill="{COLOR_TEXT_MUTED}">
      Cluster-Wise Centroid Expression Averaging Across G Marker Genes
    </text>
    <line x1="14" y1="41" x2="{w-14}" y2="41" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <g transform="translate(14, 48)">
      <rect x="0" y="0" width="160" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="80" y="14" font-size="7" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Centroid Expression Formula</text>
      <text x="80" y="34" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">
        ϕ<tspan baseline-shift="sub" font-size="70%">k</tspan> = <tspan font-style="normal">1</tspan>/|C<tspan baseline-shift="sub" font-size="70%">k</tspan>| ∑<tspan baseline-shift="sub" font-size="65%">c∈C<tspan baseline-shift="sub" font-size="65%">k</tspan></tspan> x<tspan baseline-shift="sub" font-size="70%">c</tspan>
      </text>
      <text x="80" y="46" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Average profile of cluster C_k</text>
    </g>

    <g transform="translate(186, 48)">
      <rect x="0" y="0" width="180" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="90" y="14" font-size="7" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Reference Signature Matrix</text>
      <text x="90" y="32" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">
        <tspan font-weight="bold">Φ</tspan> = [ ϕ<tspan baseline-shift="sub" font-size="70%">1</tspan>, ϕ<tspan baseline-shift="sub" font-size="70%">2</tspan>, ..., ϕ<tspan baseline-shift="sub" font-size="70%">K</tspan> ] ∈ ℝ<tspan baseline-shift="super" font-size="70%">G × K</tspan>
      </text>
      <text x="90" y="46" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">G marker genes × K cell states</text>
    </g>

    <g transform="translate(14, 108)">
      <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Lineage Centroid Heatmap (Φ):</text>
      
      <g transform="translate(0, 16)">
        <rect x="0" y="0" width="50" height="8" fill="{OKABE_BLUE}" />
        <text x="25" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">CD8+ T</text>
        <rect x="52" y="0" width="50" height="8" fill="{OKABE_SKY_BLUE}" />
        <text x="77" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">CD4+ T</text>
        <rect x="104" y="0" width="50" height="8" fill="{COLOR_RESPONDER}" />
        <text x="129" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">B Cells</text>
        <rect x="156" y="0" width="50" height="8" fill="{COLOR_MELANOMA}" />
        <text x="181" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Myeloid</text>
        <rect x="208" y="0" width="50" height="8" fill="{COLOR_MCMC}" />
        <text x="233" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Stromal</text>
        <rect x="260" y="0" width="50" height="8" fill="{COLOR_NON_RESPONDER}" />
        <text x="285" y="6.5" font-size="5.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Tumor</text>
      </g>

      <g transform="translate(0, 26)">
        <rect x="0" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.9" />
        <rect x="52" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.1" />
        <rect x="104" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.05" />
        <rect x="156" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.05" />
        <rect x="208" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.05" />
        <rect x="260" y="0" width="50" height="8" fill="{OKABE_BLUE}" opacity="0.05" />
        <text x="316" y="6.5" font-size="6" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD8A/B</text>

        <rect x="0" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.1" />
        <rect x="52" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.85" />
        <rect x="104" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
        <rect x="156" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
        <rect x="208" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
        <rect x="260" y="9" width="50" height="8" fill="{OKABE_SKY_BLUE}" opacity="0.05" />
        <text x="316" y="15.5" font-size="6" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD4/FOXP3</text>

        <rect x="0" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.05" />
        <rect x="52" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.05" />
        <rect x="104" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.9" />
        <rect x="156" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.05" />
        <rect x="208" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.05" />
        <rect x="260" y="18" width="50" height="8" fill="{COLOR_RESPONDER}" opacity="0.05" />
        <text x="316" y="24.5" font-size="6" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">MS4A1</text>

        <rect x="0" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.05" />
        <rect x="52" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.05" />
        <rect x="104" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.05" />
        <rect x="156" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.85" />
        <rect x="208" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.1" />
        <rect x="260" y="27" width="50" height="8" fill="{COLOR_MELANOMA}" opacity="0.05" />
        <text x="316" y="33.5" font-size="6" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">CD68/CD163</text>

        <rect x="0" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
        <rect x="52" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
        <rect x="104" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
        <rect x="156" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
        <rect x="208" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.05" />
        <rect x="260" y="36" width="50" height="8" fill="{COLOR_NON_RESPONDER}" opacity="0.9" />
        <text x="316" y="42.5" font-size="6" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">MLANA/PMEL</text>
      </g>
    </g>

    <g transform="translate(14, 182)">
      <text x="0" y="10" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
        Simplex Signature: ∑<tspan baseline-shift="sub" font-size="70%">g=1</tspan><tspan baseline-shift="super" font-size="70%">G</tspan> ϕ<tspan baseline-shift="sub" font-size="70%">k,g</tspan> = 1,   ϕ<tspan baseline-shift="sub" font-size="70%">k,g</tspan> ≥ 0  (∀ k ∈ {{1,...,K}})
      </text>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 5: BayesPrism Deconvolution Engine
# ==============================================================================
def generate_asset_05_bayesprism_deconv() -> str:
    """Generate Asset 5: BayesPrism Hierarchical Deconvolution Engine."""
    w, h = 460, 220
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-05-bayesprism"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-05-bayesprism" inkscape:label="05_BayesPrism_Deconvolution">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <text x="14" y="20" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      BayesPrism Deconvolution Engine
    </text>
    <text x="14" y="33" font-size="8" fill="{COLOR_TEXT_MUTED}">
      Hierarchical MCMC Gibbs Sampling &amp; Collinearity Null-Space Cancellation
    </text>
    <line x1="14" y1="41" x2="{w-14}" y2="41" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <g transform="translate(14, 48)">
      <rect x="0" y="0" width="180" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="90" y="14" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Pseudobulk Mixture Model</text>
      <text x="90" y="32" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">
        <tspan font-weight="bold">b</tspan><tspan baseline-shift="sub" font-size="70%">n</tspan> ≈ <tspan font-weight="bold">Φ</tspan> · θ<tspan baseline-shift="sub" font-size="70%">n</tspan>,   θ<tspan baseline-shift="sub" font-size="70%">n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan>
      </text>
      <text x="90" y="45" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
        ∑<tspan baseline-shift="sub" font-size="70%">k=1</tspan><tspan baseline-shift="super" font-size="70%">K</tspan> θ<tspan baseline-shift="sub" font-size="70%">k,n</tspan> = 1,   θ<tspan baseline-shift="sub" font-size="70%">k,n</tspan> ≥ 0
      </text>
    </g>

    <g transform="translate(204, 48)">
      <rect x="0" y="0" width="242" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="121" y="14" font-size="7" font-weight="700" fill="{COLOR_MCMC}" text-anchor="middle">Hierarchical Null-Space Cancellation</text>
      <text x="121" y="30" font-size="8.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        <tspan font-weight="bold">M</tspan><tspan font-weight="bold">v</tspan> = <tspan font-weight="bold">0</tspan>  ⟹  <tspan font-weight="bold">v</tspan> ∈ Null(<tspan font-weight="bold">M</tspan>)
      </text>
      <text x="121" y="44" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">
        Cancels divergent variance (κ(Φ) → ∞) across broad types
      </text>
    </g>

    <g transform="translate(14, 110)">
      <g transform="translate(0, 0)">
        <rect x="0" y="0" width="124" height="66" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Two-Tier Hierarchy:</text>
        <rect x="8" y="20" width="108" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="62" y="31" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Broad Types <tspan font-weight="bold">𝒯</tspan> (T ≪ S)
        </text>
        <rect x="8" y="40" width="108" height="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="62" y="51" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}" text-anchor="middle">
          Fine States <tspan font-weight="bold">𝒮</tspan> (S states)
        </text>
      </g>

      <path d="M 128 33 L 146 33" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

      <g transform="translate(150, 0)">
        <rect x="0" y="0" width="134" height="66" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="13" font-size="7" font-weight="700" fill="{COLOR_MCMC}">MCMC Gibbs Sampling:</text>
        <text x="8" y="26" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
          1. Z<tspan baseline-shift="sub" font-size="70%">n,g,s</tspan> ~ Multinom(X<tspan baseline-shift="sub" font-size="70%">n,g</tspan>, p)
        </text>
        <text x="8" y="38" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
          2. θ<tspan baseline-shift="sub" font-size="70%">n,s</tspan> ~ Dirichlet(α + Z)
        </text>
        <text x="8" y="50" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
          3. ψ<tspan baseline-shift="sub" font-size="70%">t</tspan> Profile Refinement
        </text>
        <text x="8" y="60" font-size="6" fill="{COLOR_TEXT_MUTED}">Convergence across 2,000 Gibbs steps</text>
      </g>

      <path d="M 288 33 L 306 33" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />

      <g transform="translate(310, 0)">
        <rect x="0" y="0" width="122" height="66" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="13" font-size="7" font-weight="700" fill="{COLOR_RESPONDER}">Estimated Cell Fractions:</text>
        <text x="8" y="24" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          θ̂<tspan baseline-shift="sub" font-size="70%">n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="70%">K-1</tspan>
        </text>
        <g transform="translate(8, 30)">
          <rect x="0" y="0" width="106" height="12" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect x="0" y="0" width="35" height="12" fill="{OKABE_BLUE}" />
          <rect x="35" y="0" width="22" height="12" fill="{OKABE_SKY_BLUE}" />
          <rect x="57" y="0" width="18" height="12" fill="{COLOR_RESPONDER}" />
          <rect x="75" y="0" width="20" height="12" fill="{COLOR_MELANOMA}" />
          <rect x="95" y="0" width="11" height="12" fill="{COLOR_NON_RESPONDER}" />
        </g>
        <text x="8" y="52" font-size="6" fill="{COLOR_TEXT_MUTED}">CD8 | CD4 | B | Mye | Tum</text>
        <text x="8" y="60" font-size="6" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Patient-specific proportions</text>
      </g>
    </g>

    <g transform="translate(14, 186)">
      <text x="0" y="12" font-size="7" font-weight="600" fill="{COLOR_TEXT_MUTED}">Multi-Tool Deconvolution Suite:</text>
      <rect x="140" y="3" width="70" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="175" y="13" font-size="6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">BayesPrism</text>
      <rect x="216" y="3" width="70" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="251" y="13" font-size="6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">InstaPrism (EM)</text>
      <rect x="292" y="3" width="80" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="332" y="13" font-size="6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">RegDeconv (FISTA)</text>
      <rect x="378" y="3" width="60" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="408" y="13" font-size="6" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">CIBERSORTx</text>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 6: Stability Selection & Elaborate Random Forest Classifier
# ==============================================================================
def generate_asset_06_stability_selection_classifier() -> str:
    """Generate Asset 6: Stability Selection & Elaborate Multi-Depth Random Forest Classifier."""
    w, h = 500, 230
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-06-stability-selection"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-06-stability-selection" inkscape:label="06_Stability_Selection_Classifier">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <text x="14" y="20" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Stability Selection &amp; Random Forest Classification
    </text>
    <text x="14" y="33" font-size="8" fill="{COLOR_TEXT_MUTED}">
      Finite-Sample Error Controlled Feature Filtering &amp; Binary Clinical Response Prediction
    </text>
    <line x1="14" y1="41" x2="{w-14}" y2="41" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <!-- Stability Selection Theory Box -->
    <g transform="translate(14, 48)">
      <rect x="0" y="0" width="230" height="96" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Finite-Sample Stability Selection (SS-CPSS)
      </text>

      <text x="10" y="28" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
        Π̂<tspan baseline-shift="sub" font-size="70%">k</tspan>(λ) = <tspan font-style="normal">1</tspan>/(2B) ∑<tspan baseline-shift="sub" font-size="65%">b=1</tspan><tspan baseline-shift="super" font-size="65%">2B</tspan> 𝕀( k ∈ Ŝ<tspan baseline-shift="super" font-size="65%">λ</tspan>(I<tspan baseline-shift="sub" font-size="65%">b</tspan>) )
      </text>

      <text x="10" y="42" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
        Ŝ<tspan baseline-shift="sub" font-size="70%">stable</tspan> = {{ k : max<tspan baseline-shift="sub" font-size="70%">λ∈Λ</tspan> Π̂<tspan baseline-shift="sub" font-size="70%">k</tspan>(λ) ≥ π<tspan baseline-shift="sub" font-size="70%">thr</tspan> }}  (π<tspan baseline-shift="sub" font-size="70%">thr</tspan> = 0.60)
      </text>

      <!-- PFER Error Bound Box -->
      <rect x="8" y="50" width="214" height="38" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="62" font-size="6.8" font-weight="700" fill="{COLOR_NON_RESPONDER}">Shah-Samworth Finite-Sample PFER Bound:</text>
      <text x="12" y="76" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
        E[ V ] ≤ <tspan font-style="normal">1</tspan>/(2π<tspan baseline-shift="sub" font-size="70%">thr</tspan> - 1) · q<tspan baseline-shift="sub" font-size="70%">Λ</tspan><tspan baseline-shift="super" font-size="70%">2</tspan> / p
      </text>
    </g>

    <!-- Elaborated Random Forest Classifier (Multi-Depth Decision Trees) -->
    <g transform="translate(254, 48)">
      <rect x="0" y="0" width="232" height="96" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        Random Forest Ensemble (Decision Trees)
      </text>

      <!-- 3 Elaborate Decision Trees with 3 Depth Levels -->
      <g transform="translate(8, 20)">
        <!-- Tree 1 -->
        <g transform="translate(0, 0)">
          <!-- Root -->
          <circle cx="26" cy="6" r="2.8" fill="{OKABE_BLUE}" />
          <!-- Level 1 Splits -->
          <line x1="26" y1="9" x2="14" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="26" y1="9" x2="38" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="14" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <circle cx="38" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <!-- Level 2 Splits to Leaves -->
          <line x1="14" y1="20" x2="8" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="14" y1="20" x2="20" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="32" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="44" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="8" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="20" cy="28" r="2.0" fill="{COLOR_NON_RESPONDER}" />
          <circle cx="32" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="44" cy="28" r="2.0" fill="{COLOR_NON_RESPONDER}" />
          <text x="26" y="38" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree 1</text>
        </g>

        <!-- Tree 2 -->
        <g transform="translate(56, 0)">
          <circle cx="26" cy="6" r="2.8" fill="{OKABE_BLUE}" />
          <line x1="26" y1="9" x2="14" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="26" y1="9" x2="38" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="14" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <circle cx="38" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <line x1="14" y1="20" x2="8" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="14" y1="20" x2="20" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="32" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="44" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="8" cy="28" r="2.0" fill="{COLOR_NON_RESPONDER}" />
          <circle cx="20" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="32" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="44" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <text x="26" y="38" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree 2</text>
        </g>

        <!-- Tree B -->
        <g transform="translate(112, 0)">
          <circle cx="26" cy="6" r="2.8" fill="{OKABE_BLUE}" />
          <line x1="26" y1="9" x2="14" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="26" y1="9" x2="38" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="14" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <circle cx="38" cy="18" r="2.2" fill="{COLOR_TEXT_PRIMARY}" />
          <line x1="14" y1="20" x2="8" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="14" y1="20" x2="20" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="32" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="38" y1="20" x2="44" y2="28" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <circle cx="8" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="20" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <circle cx="32" cy="28" r="2.0" fill="{COLOR_NON_RESPONDER}" />
          <circle cx="44" cy="28" r="2.0" fill="{COLOR_RESPONDER}" />
          <text x="26" y="38" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Tree B</text>
        </g>
      </g>

      <!-- Prediction Equation & Output Badges -->
      <g transform="translate(10, 62)">
        <text x="0" y="8" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          Ensemble Vote: P̂( Y=1 | θ̂ ) = <tspan font-style="normal">1</tspan>/B ∑ T<tspan baseline-shift="sub" font-size="65%">b</tspan>( θ̂ )
        </text>
        <g transform="translate(0, 14)">
          <rect width="102" height="15" fill="{COLOR_RESPONDER}" />
          <text x="51" y="10" font-size="6.2" font-weight="700" fill="#FFFFFF" text-anchor="middle">Responder (R)</text>

          <rect x="108" width="104" height="15" fill="{COLOR_NON_RESPONDER}" />
          <text x="160" y="10" font-size="6.2" font-weight="700" fill="#FFFFFF" text-anchor="middle">Non-Responder (NR)</text>
        </g>
      </g>
    </g>

    <!-- Bottom ROC-AUC & Multi-Cohort Performance -->
    <g transform="translate(14, 150)">
      <rect width="472" height="68" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Mini ROC Curve Graphic -->
      <g transform="translate(14, 8)">
        <line x1="0" y1="46" x2="60" y2="46" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="0" y1="0" x2="0" y2="46" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="0" y1="46" x2="60" y2="0" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" stroke-dasharray="2,2" />
        <path d="M 0 46 C 8 14, 22 6, 60 0" fill="none" stroke="{OKABE_BLUE}" stroke-width="1.8" />
        <text x="30" y="55" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">1 - Specificity</text>
      </g>

      <g transform="translate(90, 10)">
        <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cross-Cohort Discrimination (9 Clinical ICB Trials):</text>
        <text x="0" y="24" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
          Pooled ROC-AUC = 0.812  (95% CI: 0.76 - 0.86, Hanley-McNeil test)
        </text>
        <text x="0" y="38" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
          • <tspan font-weight="600" fill="{COLOR_RESPONDER}">CD8+ Tem &amp; TCF7+ Progenitor Tex</tspan>: Selected in 100% of CPSS subsamples (OR = 2.45)
        </text>
        <text x="0" y="50" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
          • <tspan font-weight="600" fill="{COLOR_NON_RESPONDER}">CD163+ / MRC1+ M2 Macrophages</tspan>: Stable negative predictor of response (OR = 0.42)
        </text>
      </g>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 7: Single-Cell Milo Differential Abundance
# ==============================================================================
def generate_asset_07_sc_milo_differential_abundance() -> str:
    """Generate Asset 7: Single-Cell Continuous Neighborhood Testing (Milo)."""
    w, h = 440, 240
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-07-milo-da"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-07-milo-da" inkscape:label="07_Milo_Differential_Abundance">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <text x="14" y="20" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Single-Cell Milo Continuous Differential Abundance
    </text>
    <text x="14" y="33" font-size="8" fill="{COLOR_TEXT_MUTED}">
      Direct Attribution of Clinical Response Labels to Individual Cells via kNN Graphs
    </text>
    <line x1="14" y1="41" x2="{w-14}" y2="41" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <g transform="translate(14, 48)">
      <rect x="0" y="0" width="220" height="96" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Negative Binomial Neighborhood GLM</text>
      <text x="10" y="30" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
        log( E[ N<tspan baseline-shift="sub" font-size="70%">v,j</tspan> ] ) = β<tspan baseline-shift="sub" font-size="70%">0,v</tspan> + β<tspan baseline-shift="sub" font-size="70%">v</tspan><tspan baseline-shift="super" font-size="70%">milo</tspan> · Y<tspan baseline-shift="sub" font-size="70%">j</tspan> + log( TotalCells<tspan baseline-shift="sub" font-size="70%">j</tspan> )
      </text>
      <text x="10" y="46" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
        M<tspan baseline-shift="sub" font-size="70%">c,v</tspan> = 1 if c ∈ 𝒩(v)  [Incidence M ∈ {{0, 1}}<tspan baseline-shift="super" font-size="70%">C × V</tspan>]
      </text>
      <rect x="8" y="54" width="204" height="34" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="65" font-size="6.8" font-weight="700" fill="{COLOR_MCMC}">Spatial False Discovery Rate (FDR):</text>
      <text x="14" y="78" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
        Spatial-FDR(v) = ∑<tspan baseline-shift="sub" font-size="65%">u</tspan> w<tspan baseline-shift="sub" font-size="65%">uv</tspan> 𝕀(p<tspan baseline-shift="sub" font-size="65%">u</tspan> ≤ p<tspan baseline-shift="sub" font-size="65%">v</tspan>) / ∑<tspan baseline-shift="sub" font-size="65%">u</tspan> w<tspan baseline-shift="sub" font-size="65%">uv</tspan> · p<tspan baseline-shift="sub" font-size="65%">v</tspan> / F̂(p<tspan baseline-shift="sub" font-size="65%">v</tspan>)
      </text>
    </g>

    <g transform="translate(244, 48)">
      <rect x="0" y="0" width="182" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Graph Neighborhoods 𝒩(v)</text>
      <g transform="translate(18, 20)">
        <line x1="20" y1="20" x2="60" y2="15" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="20" y1="20" x2="35" y2="50" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="60" y1="15" x2="105" y2="22" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="35" y1="50" x2="75" y2="55" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="60" y1="15" x2="75" y2="55" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="105" y1="22" x2="130" y2="45" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="75" y1="55" x2="130" y2="45" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

        <circle cx="35" cy="28" r="22" fill="{COLOR_RESPONDER}" opacity="0.15" stroke="{COLOR_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />
        <circle cx="95" cy="35" r="22" fill="{COLOR_NON_RESPONDER}" opacity="0.15" stroke="{COLOR_NON_RESPONDER}" stroke-width="0.75" stroke-dasharray="2,2" />

        <circle cx="20" cy="20" r="3.5" fill="{COLOR_RESPONDER}" />
        <circle cx="35" cy="50" r="4.5" fill="{COLOR_RESPONDER}" stroke="#0F172A" stroke-width="0.75" />
        <circle cx="60" cy="15" r="3.5" fill="{COLOR_RESPONDER}" />
        <circle cx="75" cy="55" r="3.5" fill="{COLOR_NON_RESPONDER}" />
        <circle cx="105" cy="22" r="4.5" fill="{COLOR_NON_RESPONDER}" stroke="#0F172A" stroke-width="0.75" />
        <circle cx="130" cy="45" r="3.5" fill="{COLOR_NON_RESPONDER}" />

        <text x="35" y="70" font-size="5.5" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">𝒩(v<tspan baseline-shift="sub" font-size="70%">R</tspan>): β<tspan baseline-shift="super" font-size="70%">milo</tspan> &gt; 0</text>
        <text x="110" y="70" font-size="5.5" font-weight="700" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">𝒩(v<tspan baseline-shift="sub" font-size="70%">NR</tspan>): β<tspan baseline-shift="super" font-size="70%">milo</tspan> &lt; 0</text>
      </g>
    </g>

    <g transform="translate(14, 152)">
      <rect x="0" y="0" width="412" height="74" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Attributed Single-Cell Response Ground Truth:</text>
      <g transform="translate(10, 24)">
        <rect x="0" y="0" width="190" height="42" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="7" font-weight="700" fill="{COLOR_RESPONDER}">Responder-Enriched Cells (Spatial FDR &lt; 0.10):</text>
        <text x="8" y="24" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• CD8+ T effector memory (GZMB+, PRF1+)</text>
        <text x="8" y="34" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• TCF7+ stem-like progenitor exhausted T cells</text>

        <rect x="200" y="0" width="192" height="42" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="7" font-weight="700" fill="{COLOR_NON_RESPONDER}">Non-Responder-Enriched Cells (Spatial FDR &lt; 0.10):</text>
        <text x="8" y="24" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• Terminally exhausted T cells (HAVCR2+, LAG3+)</text>
        <text x="8" y="34" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• M2-polarized suppressive macrophages (CD163+)</text>
      </g>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# ASSET 8: In Silico Pseudo-bulk & Triangulation Concordance
# ==============================================================================
def generate_asset_08_pseudobulk_triangulation_concordance() -> str:
    """Generate Asset 8: In Silico Pseudo-bulk & Triangulation Concordance."""
    w, h = 480, 240
    svg = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="{w}" height="{h}" viewBox="0 0 {w} {h}" version="1.1" id="asset-08-triangulation"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  {svg_defs()}
  <g inkscape:groupmode="layer" id="layer-08-triangulation" inkscape:label="08_Triangulation_Concordance">
    <rect x="0.5" y="0.5" width="{w-1}" height="{h-1}" fill="{COLOR_CONTAINER_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <text x="14" y="20" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      In Silico Pseudo-bulk &amp; Triangulation Concordance
    </text>
    <text x="14" y="33" font-size="8" fill="{COLOR_TEXT_MUTED}">
      Cross-Checking Bulk Deconvolution Fractions Against Single-Cell Milo Ground Truth
    </text>
    <line x1="14" y1="41" x2="{w-14}" y2="41" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

    <g transform="translate(14, 48)">
      <rect x="0" y="0" width="220" height="96" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">In Silico Pseudo-bulk Aggregation</text>
      <text x="10" y="30" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
        <tspan font-weight="bold">b</tspan><tspan baseline-shift="sub" font-size="70%">j</tspan><tspan baseline-shift="super" font-size="70%">pseudo</tspan> = ∑<tspan baseline-shift="sub" font-size="65%">c ∈ Patient<tspan baseline-shift="sub" font-size="65%">j</tspan></tspan> <tspan font-weight="bold">x</tspan><tspan baseline-shift="sub" font-size="70%">c</tspan>
      </text>
      <text x="10" y="42" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Pools single cells per patient</text>

      <rect x="8" y="50" width="204" height="38" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="62" font-size="6.8" font-weight="700" fill="{COLOR_MCMC}">Deconv &amp; Logistic Regression on Fractions:</text>
      <text x="12" y="74" font-size="7" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
        logit( P( Y<tspan baseline-shift="sub" font-size="70%">j</tspan>=1 ) ) = β<tspan baseline-shift="sub" font-size="70%">0</tspan> + ∑<tspan baseline-shift="sub" font-size="65%">k</tspan> β<tspan baseline-shift="sub" font-size="70%">k</tspan><tspan baseline-shift="super" font-size="70%">deconv</tspan> · θ̂<tspan baseline-shift="sub" font-size="70%">k,j</tspan><tspan baseline-shift="super" font-size="70%">pseudo</tspan>
      </text>
      <text x="12" y="84" font-size="6" fill="{COLOR_TEXT_MUTED}">Yields cell state effect size β_deconv</text>
    </g>

    <g transform="translate(244, 48)">
      <rect x="0" y="0" width="222" height="178" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Triangulation Concordance Scatter Plot</text>
      <g transform="translate(30, 24)">
        <line x1="80" y1="10" x2="80" y2="130" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <line x1="10" y1="70" x2="150" y2="70" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />

        <rect x="80" y="10" width="70" height="60" fill="{COLOR_RESPONDER}" opacity="0.08" />
        <text x="115" y="24" font-size="6" font-weight="700" fill="{COLOR_RESPONDER}" text-anchor="middle">Concordant R</text>

        <rect x="10" y="70" width="70" height="60" fill="{COLOR_NON_RESPONDER}" opacity="0.08" />
        <text x="45" y="122" font-size="6" font-weight="700" fill="{COLOR_NON_RESPONDER}" text-anchor="middle">Concordant NR</text>

        <text x="45" y="24" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Discordant</text>
        <text x="115" y="122" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Discordant</text>

        <circle cx="120" cy="35" r="3.5" fill="{OKABE_BLUE}" />
        <text x="125" y="33" font-size="5.5" font-weight="600" fill="{OKABE_BLUE}">CD8 Tem</text>

        <circle cx="105" cy="48" r="3" fill="{OKABE_BLUE}" />
        <text x="111" y="47" font-size="5" fill="{COLOR_TEXT_SECONDARY}">TCF7+ Tex</text>

        <circle cx="95" cy="38" r="3" fill="{COLOR_RESPONDER}" />

        <circle cx="45" cy="100" r="3.5" fill="{COLOR_NON_RESPONDER}" />
        <text x="15" y="99" font-size="5.5" font-weight="600" fill="{COLOR_NON_RESPONDER}">M2 Macro</text>

        <circle cx="55" cy="88" r="3" fill="{COLOR_NON_RESPONDER}" />
        <text x="25" y="86" font-size="5" fill="{COLOR_TEXT_SECONDARY}">Term Tex</text>

        <circle cx="65" cy="110" r="2.5" fill="{COLOR_MELANOMA}" />

        <line x1="20" y1="120" x2="140" y2="20" stroke="{OKABE_BLUE}" stroke-width="1.2" stroke-dasharray="3,2" />

        <text x="80" y="142" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Single-Cell Milo Effect (β<tspan baseline-shift="super" font-size="70%">milo</tspan>) →
        </text>
        <text x="-70" y="-14" font-size="6.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" transform="rotate(-90)">
          Deconvolution Effect (β<tspan baseline-shift="super" font-size="70%">deconv</tspan>) →
        </text>
      </g>
    </g>

    <g transform="translate(14, 152)">
      <rect x="0" y="0" width="220" height="74" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Triangulation Validation Metrics:</text>
      <text x="10" y="28" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{COLOR_RESPONDER}">
        • Spearman Rank Correlation: ρ = 0.84  (p &lt; 10<tspan baseline-shift="super" font-size="70%">-4</tspan>)
      </text>
      <text x="10" y="40" font-size="7" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
        • Directional Concordance: 87.5% in diagonal Q1 &amp; Q3
      </text>
      <text x="10" y="52" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">
        • Discordance Diagnosis: Collinearity leakage / cell size
      </text>
      <text x="10" y="64" font-size="6" fill="{COLOR_TEXT_MUTED}">
        scripts/concordance_calibration/05_plot_triangulation_figures.py
      </text>
    </g>
  </g>
</svg>
"""
    return svg


# ==============================================================================
# MAIN EXECUTOR
# ==============================================================================
def main() -> None:
    output_dir = Path("assets")
    output_dir.mkdir(parents=True, exist_ok=True)

    generators = [
        ("01_bulk_icb_cohorts.svg", generate_asset_01_bulk_icb_cohorts),
        ("02_sc_database_h5ad.svg", generate_asset_02_sc_database_h5ad),
        ("03_automated_sc_pipeline.svg", generate_asset_03_automated_sc_pipeline),
        ("04_cell_state_profiles.svg", generate_asset_04_cell_state_profiles),
        ("05_bayesprism_deconv.svg", generate_asset_05_bayesprism_deconv),
        ("06_stability_selection_classifier.svg", generate_asset_06_stability_selection_classifier),
        ("07_sc_milo_differential_abundance.svg", generate_asset_07_sc_milo_differential_abundance),
        ("08_pseudobulk_triangulation_concordance.svg", generate_asset_08_pseudobulk_triangulation_concordance),
    ]

    for filename, gen_func in generators:
        filepath = output_dir / filename
        svg_content = gen_func()
        filepath.write_text(svg_content, encoding="utf-8")
        print(f"Generated asset: {filepath}")


if __name__ == "__main__":
    main()
