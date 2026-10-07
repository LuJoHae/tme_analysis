#!/usr/bin/env python3
"""
Generate Pure Dataset Compendium Overview Figure (Data-Only).

Implements the authentic Nature Methods / Nature editorial art style:
- Double-column landscape (1400 px width, 900 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to data marks (bars, tiles, dots)
- Authentic editorial letter badges (a, b, c, d)
- Exclusively summarizes datasets, clinical stratifications, assay completeness, and lineages
- Native Inkscape layers and semantic groups
"""

from pathlib import Path
import html
from nature_style_config import (
    WIDTH,
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS
)

HEIGHT = 900


def render_svg() -> str:
    """Declaratively render the complete Nature Methods pure dataset compendium SVG."""
    svg_parts: list[str] = []

    # 1. XML Header and SVG root with Inkscape namespaces
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-dataset-compendium-nature"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
""")

    # ----------------------------------------------------
    # LAYER 0: Canvas Background & Panel Containers
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND & PANEL CONTAINERS -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />

    <!-- Panel A Container Card (Top-Left) -->
    <rect id="card-panel-a" x="20" y="68" width="665" height="402"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel B Container Card (Top-Right) -->
    <rect id="card-panel-b" x="715" y="68" width="665" height="402"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel C Container Card (Bottom-Left) -->
    <rect id="card-panel-c" x="20" y="484" width="665" height="396"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel D Container Card (Bottom-Right) -->
    <rect id="card-panel-d" x="715" y="484" width="665" height="396"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 1: Figure Header, Subtitle & Panel Badges
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 01: FIGURE HEADER & PANEL BADGES -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <!-- Main Title -->
    <text id="txt-title" x="20" y="30" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.1">
      Study Dataset Compendium: Multi-Omic Immunotherapy Trials, Single-Cell References, and Primary Baseline Cohorts
    </text>
    <!-- Subtitle Summary -->
    <text id="txt-subtitle" x="20" y="48" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Multi-modal data matrix summarizing 9 clinical trials (n = 1,097), 5 TCGA baseline cohorts (n = 2,932), and pan-cancer single-cell references (40,002 cells)
    </text>

    <!-- Top Total Stat Box -->
    <g id="grp-header-stat" transform="translate(1085, 14)">
      <rect x="0" y="0" width="295" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="147" y="18" font-size="9" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        Compendium: &gt;4,080 Biopsies + 40,002 Single Cells
      </text>
    </g>

    <!-- Panel Badges [a], [b], [c], [d] -->
    <g id="badge-panel-a" transform="translate(34, 86)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">a</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Multi-Cohort Dataset Ecosystem &amp; Specimen Hierarchy</text>
    </g>

    <g id="badge-panel-b" transform="translate(729, 86)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">b</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Clinical Trial Cohorts: Sample Size &amp; RECIST Response Rates (n = 1,097)</text>
    </g>

    <g id="badge-panel-c" transform="translate(34, 502)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">c</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Multi-Omic Assay Completeness Matrix Across Cohorts</text>
    </g>

    <g id="badge-panel-d" transform="translate(729, 502)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">d</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cellular Lineage Composition &amp; Genomic Characteristics</text>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: PANEL A — Multi-Cohort Dataset Ecosystem
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: PANEL A — DATASET ECOSYSTEM & SPECIMEN HIERARCHY -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Dataset_Ecosystem">

    <!-- Card 1: Single-Cell Reference & Clinical Validation -->
    <g id="grp-card-scrna" transform="translate(34, 114)">
      <rect width="198" height="342" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Single-Cell Atlases (40k cells)</text>
      <line x1="10" y1="26" x2="188" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Pan-cancer reference details -->
      <g transform="translate(8, 34)">
        <rect width="182" height="120" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Pan-Cancer Atlas (40,002 c)</text>
        <text x="10" y="32" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 6 Lineages, 22 Phenotypes</text>
        <text x="10" y="46" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 58 Silhouette K-means States</text>
        <text x="10" y="60" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Leiden Clustering res 0.1–1.5</text>
        <text x="10" y="74" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Truncated Selective Inf. p &lt; 10⁻¹⁶</text>
        
        <!-- Lineage marks -->
        <g transform="translate(10, 84)">
          <rect x="0" y="0" width="24" height="14" fill="{OKABE_BLUE}" /><text x="12" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">CD8</text>
          <rect x="28" y="0" width="24" height="14" fill="{OKABE_SKY_BLUE}" /><text x="40" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">CD4</text>
          <rect x="56" y="0" width="24" height="14" fill="{OKABE_ORANGE}" /><text x="68" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Mye</text>
          <rect x="84" y="0" width="24" height="14" fill="{OKABE_BLUISH_GREEN}" /><text x="96" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">B</text>
          <rect x="112" y="0" width="24" height="14" fill="{COLOR_TEXT_MUTED}" /><text x="124" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Str</text>
          <rect x="140" y="0" width="20" height="14" fill="{OKABE_VERMILION}" /><text x="150" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Mal</text>
        </g>
      </g>

      <!-- Matched clinical cohort: Sade-Feldman -->
      <g transform="translate(8, 164)">
        <rect width="182" height="102" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Sade-Feldman (Cell 2018)</text>
        <text x="10" y="32" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Matched Melanoma Biopsies</text>
        <text x="10" y="46" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• n = 51 biopsies (32 patients)</text>
        <text x="10" y="60" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 16,291 CD45+ cells (12 states)</text>
        <text x="10" y="74" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Smart-seq2: 20 Pre- / 31 Post-ICB</text>
        <rect x="10" y="82" width="162" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="91" y="92.5" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Matched Response Ground Truth</text>
      </g>

      <!-- Additional References -->
      <g transform="translate(8, 274)">
        <rect width="182" height="58" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Supporting References</text>
        <text x="10" y="32" font-size="8" fill="{COLOR_TEXT_MUTED}">• Jerby-Arnon 2018: 7,186 c (GSE115978)</text>
        <text x="10" y="46" font-size="8" fill="{COLOR_TEXT_MUTED}">• Li 2019 (Melanoma) | Ma 2019 (Liver)</text>
      </g>
    </g>

    <!-- Card 2: Clinical Immunotherapy Trial Cohorts (n = 1,097) -->
    <g id="grp-card-trials" transform="translate(242, 114)">
      <rect width="240" height="342" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Clinical ICI Trial Cohorts (n = 1,097)</text>
      <line x1="10" y1="26" x2="230" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Melanoma -->
      <g transform="translate(8, 34)">
        <rect width="224" height="60" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_ORANGE}">Melanoma (n = 347 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Hugo (n=27, Pembro) | Riaz (n=107, Nivo)</text>
        <text x="8" y="46" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Liu (n=122, Anti-PD1) | Gide (n=91, Combo)</text>
      </g>

      <!-- Bladder -->
      <g transform="translate(8, 102)">
        <rect width="224" height="48" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_SKY_BLUE}">Bladder (n = 347 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Rosenberg / IMvigor210 (Anti-PD-L1 Atezo)</text>
      </g>

      <!-- Renal -->
      <g transform="translate(8, 158)">
        <rect width="224" height="48" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Renal Cell (n = 279 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• McDermott / IMmotion150 (263) | Choueiri (16)</text>
      </g>

      <!-- Pancreatic & Breast -->
      <g transform="translate(8, 214)">
        <rect width="224" height="56" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_REDDISH_PURPLE}">PDAC / BC (n = 124 patients)</text>
        <text x="8" y="30" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Padron (n=93, PDAC, Nivo + Sotigalimab)</text>
        <text x="8" y="44" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Anders (n=31, TNBC, Atezo + Chemo)</text>
      </g>

      <!-- Extended Compendium -->
      <g transform="translate(8, 276)">
        <rect width="224" height="56" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">+ Extended Validation Cohorts</text>
        <text x="8" y="30" font-size="7.2" fill="{COLOR_TEXT_MUTED}">VanAllen (110) | Snyder (25) | Rose (50) | Ravi (44)</text>
        <text x="8" y="44" font-size="7.2" fill="{COLOR_TEXT_MUTED}">Prat (65) | Freeman (104) | Auslander | Lauss</text>
      </g>
    </g>

    <!-- Card 3: TCGA Baseline Reference Cohorts (n = 2,932) -->
    <g id="grp-card-tcga" transform="translate(492, 114)">
      <rect width="182" height="342" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. TCGA Baselines (n = 2,932)</text>
      <line x1="10" y1="26" x2="172" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- 5 Primary Types Stack -->
      <g transform="translate(8, 34)">
        <rect width="166" height="162" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">5 Untreated Primary Types</text>
        
        <g transform="translate(10, 30)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-SKCM</text>
          <rect x="70" y="2" width="52" height="11" fill="{OKABE_ORANGE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">471</text>
        </g>
        <g transform="translate(10, 48)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-BLCA</text>
          <rect x="70" y="2" width="46" height="11" fill="{OKABE_SKY_BLUE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">414</text>
        </g>
        <g transform="translate(10, 66)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-BRCA</text>
          <rect x="70" y="2" width="85" height="11" fill="{OKABE_REDDISH_PURPLE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1,102</text>
        </g>
        <g transform="translate(10, 84)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-KIRC</text>
          <rect x="70" y="2" width="60" height="11" fill="{OKABE_BLUISH_GREEN}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">538</text>
        </g>
        <g transform="translate(10, 102)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-PAAD</text>
          <rect x="70" y="2" width="20" height="11" fill="{COLOR_TEXT_MUTED}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">178</text>
        </g>
        <text x="10" y="146" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Total Untreated: 2,703–2,932</text>
      </g>

      <!-- Baseline Microenvironment summary -->
      <g transform="translate(8, 206)">
        <rect width="166" height="126" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="18" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Baseline Microenvironment</text>
        <text x="10" y="34" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• K = 100 Reference Centroids</text>
        <text x="10" y="48" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 58 Deconvolution States</text>
        <text x="10" y="62" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Primary Organ State Space</text>
        <rect x="10" y="72" width="146" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="14" y="87" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Unperturbed TME Baseline</text>
        <text x="14" y="100" font-size="7" fill="{COLOR_TEXT_MUTED}">Anchors metastatic remodeling</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: PANEL B — Clinical Trial Stratification & Response
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 03: PANEL B — CLINICAL TRIAL STRATIFICATION & RECIST -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Trial_Stratification">

    <!-- Container Table & Bar Grid -->
    <g id="grp-trial-bars" transform="translate(729, 114)">
      <!-- Table Header Bar -->
      <rect width="637" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Trial Cohort</text>
      <text x="110" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cancer Type</text>
      <text x="185" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Therapeutic Agent</text>
      <text x="310" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">N</text>
      <text x="345" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">RECIST Response Rate (CR/PR vs. SD/PD)</text>
      <text x="590" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">ORR %</text>

      <!-- Row 1: Hugo -->
      <g transform="translate(0, 30)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Hugo 2016</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_ORANGE}">Melanoma</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Pembrolizumab (PD-1)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">27</text>
        <!-- Stacked Bar (Width 220px total) -->
        <rect x="345" y="2" width="114" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="459" y="2" width="106" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">14 R</text>
        <text x="459" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">13 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">51.9%</text>
      </g>

      <!-- Row 2: Riaz -->
      <g transform="translate(0, 56)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Riaz 2017</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_ORANGE}">Melanoma</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Nivolumab (PD-1)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">107</text>
        <rect x="345" y="2" width="70" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="415" y="2" width="150" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="8">34 R</text>
        <text x="415" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">73 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">31.8%</text>
      </g>

      <!-- Row 3: Liu -->
      <g transform="translate(0, 82)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Liu 2019</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_ORANGE}">Melanoma</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Nivo / Pembro (PD-1)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">122</text>
        <rect x="345" y="2" width="86" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="431" y="2" width="134" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="8">48 R</text>
        <text x="431" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">74 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">39.3%</text>
      </g>

      <!-- Row 4: Gide -->
      <g transform="translate(0, 108)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Gide 2019</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_ORANGE}">Melanoma</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Anti-PD-1 +/- Ipi (Combo)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">91</text>
        <rect x="345" y="2" width="152" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="497" y="2" width="68" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">63 R</text>
        <text x="497" y="13" font-size="7" font-weight="700" fill="#FFF" dx="8">28 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">69.2%</text>
      </g>

      <!-- Row 5: Rosenberg -->
      <g transform="translate(0, 134)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Rosenberg 2016</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_SKY_BLUE}">Bladder</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Atezolizumab (PD-L1)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">347</text>
        <rect x="345" y="2" width="43" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="388" y="2" width="177" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="5">68 R</text>
        <text x="388" y="13" font-size="7" font-weight="700" fill="#FFF" dx="12">279 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">19.6%</text>
      </g>

      <!-- Row 6: McDermott -->
      <g transform="translate(0, 160)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">McDermott 2018</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Renal Cell</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Atezo +/- Bevacizumab</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">263</text>
        <rect x="345" y="2" width="81" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="426" y="2" width="139" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="8">97 R</text>
        <text x="426" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">166 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">36.9%</text>
      </g>

      <!-- Row 7: Choueiri -->
      <g transform="translate(0, 186)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Choueiri 2016</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Renal Cell</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Nivolumab (PD-1)</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">16</text>
        <rect x="345" y="2" width="82" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="427" y="2" width="138" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="8">6 R</text>
        <text x="427" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">10 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">37.5%</text>
      </g>

      <!-- Row 8: Padron -->
      <g transform="translate(0, 212)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Padron 2022</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">PDAC</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Nivo + Sotigalimab + Chemo</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">93</text>
        <rect x="345" y="2" width="38" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="383" y="2" width="182" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="4">16</text>
        <text x="383" y="13" font-size="7" font-weight="700" fill="#FFF" dx="12">77 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">17.2%</text>
      </g>

      <!-- Row 9: Anders -->
      <g transform="translate(0, 238)">
        <text x="10" y="14" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Anders 2021</text>
        <text x="110" y="14" font-size="7.5" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">TNBC</text>
        <text x="185" y="14" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Atezo + Nab-paclitaxel</text>
        <text x="310" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">31</text>
        <rect x="345" y="2" width="78" height="15" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="423" y="2" width="142" height="15" fill="{OKABE_VERMILION}" />
        <text x="345" y="13" font-size="7" font-weight="700" fill="#FFF" dx="6">11 R</text>
        <text x="423" y="13" font-size="7" font-weight="700" fill="#FFF" dx="10">20 NR</text>
        <text x="590" y="14" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">35.5%</text>
      </g>

      <!-- Overall Summary Row -->
      <g transform="translate(0, 268)">
        <rect width="637" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="21" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Overall Trial Compendium:</text>
        <text x="150" y="21" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">N = 1,097 Total</text>
        
        <!-- Legend Indicator -->
        <rect x="270" y="11" width="12" height="12" fill="{OKABE_BLUISH_GREEN}" />
        <text x="286" y="20.5" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Responders: 357 (32.5%)</text>

        <rect x="430" y="11" width="12" height="12" fill="{OKABE_VERMILION}" />
        <text x="446" y="20.5" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Non-Responders: 740 (67.5%)</text>
      </g>

      <!-- Biopsy Timing Breakdown Strip -->
      <g transform="translate(0, 310)">
        <rect width="310" height="26" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="17" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Pre-Treatment Baseline Biopsies: n = 852 (77.7%)</text>

        <rect x="327" y="0" width="310" height="26" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="339" y="17" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">On-Treatment / Post-Progression: n = 245 (22.3%)</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: PANEL C — Multi-Omic Assay Completeness Matrix
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: PANEL C — MULTI-OMIC ASSAY MATRIX -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Assay_Matrix">

    <g id="grp-assay-matrix" transform="translate(34, 528)">
      <!-- Matrix Table Wireframe -->
      <rect width="637" height="340" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
""")

    modalities = [
        ("Bulk RNA-seq", "TPM / raw counts", OKABE_BLUE),
        ("Single-Cell RNA-seq", "Smart-seq2 / 10x Chromium", OKABE_SKY_BLUE),
        ("WES / Somatic Mutations", "Binary somatic alteration calls", OKABE_ORANGE),
        ("Subclonal VAF Spectrum", "Allele frequencies (0.01 – 0.50)", OKABE_ORANGE),
        ("Tumor Mutational Burden", "Exome-wide TMB quantification", OKABE_ORANGE),
        ("RECIST 1.1 Response", "Annotated CR/PR vs SD/PD", OKABE_BLUISH_GREEN),
        ("Longitudinal Sampling", "Paired Pre- vs On/Post-treatment", OKABE_REDDISH_PURPLE),
        ("Survival Endpoints", "Annotated PFS / OS intervals", OKABE_BLUE),
    ]

    cohorts = [
        ("Hugo", True, False, True, True, True, True, False, True),
        ("Riaz", True, False, True, True, True, True, True, True),
        ("Liu", True, False, True, True, True, True, False, True),
        ("Gide", True, False, True, True, True, True, False, True),
        ("Rosen", True, False, True, True, True, True, False, True),
        ("McDer", True, False, True, True, True, True, False, True),
        ("Choue", True, False, True, True, True, True, False, True),
        ("Padro", True, False, True, True, True, True, False, True),
        ("Ander", True, False, True, True, True, True, False, True),
        ("Sade", True, True, False, False, False, True, True, False),
        ("SKCM", True, False, True, False, True, False, False, True),
        ("BLCA", True, False, True, False, True, False, False, True),
        ("BRCA", True, False, True, False, True, False, False, True),
        ("KIRC", True, False, True, False, True, False, False, True),
        ("PAAD", True, False, True, False, True, False, False, True),
    ]

    # Generate matrix header
    svg_parts.append("""      <!-- Cohort Column Headers -->
      <g transform="translate(150, 22)">
""")
    for i, (c_name, *_) in enumerate(cohorts):
        cx = i * 32 + 16
        is_trial = i < 9
        is_sc = i == 9
        col_text_color = COLOR_TEXT_PRIMARY if is_trial else (OKABE_BLUE if is_sc else COLOR_TEXT_MUTED)
        svg_parts.append(f"""        <text x="{cx}" y="0" font-size="8" font-weight="700" fill="{col_text_color}" text-anchor="middle">{c_name}</text>
""")
    svg_parts.append(f"""      </g>
      <!-- Sub-headers for Cohort groups -->
      <g transform="translate(150, 30)">
        <line x1="2" y1="0" x2="286" y2="0" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="1" />
        <line x1="290" y1="0" x2="318" y2="0" stroke="{OKABE_BLUE}" stroke-width="1" />
        <line x1="322" y1="0" x2="478" y2="0" stroke="{COLOR_TEXT_MUTED}" stroke-width="1" />
        <text x="144" y="12" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">9 Clinical Trial Cohorts</text>
        <text x="304" y="12" font-size="7" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">scRNA</text>
        <text x="400" y="12" font-size="7" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">5 TCGA Baselines</text>
      </g>
""")

    # Generate matrix rows
    for row_idx, (m_name, m_desc, m_color) in enumerate(modalities):
        ry = 56 + row_idx * 33
        bg_fill = COLOR_SUBTLE_FILL if row_idx % 2 == 0 else COLOR_PANEL_BG
        svg_parts.append(f"""
      <!-- Row {row_idx}: {m_name} -->
      <g transform="translate(6, {ry})">
        <rect width="625" height="28" fill="{bg_fill}" />
        <text x="14" y="12" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">{m_name}</text>
        <text x="14" y="23" font-size="6.5" fill="{COLOR_TEXT_MUTED}">{m_desc}</text>
        
        <!-- Cohort Assay Dots -->
        <g transform="translate(144, 0)">
""")
        for c_idx, (_, *flags) in enumerate(cohorts):
            has_assay = flags[row_idx]
            dot_x = c_idx * 32 + 16
            dot_y = 14
            if has_assay:
                svg_parts.append(f"""          <circle cx="{dot_x}" cy="{dot_y}" r="5" fill="{m_color}" />
""")
            else:
                svg_parts.append(f"""          <circle cx="{dot_x}" cy="{dot_y}" r="1.5" fill="{COLOR_BORDER_HAIRLINE}" />
""")
        svg_parts.append("""        </g>
      </g>
""")

    # Bottom Legend inside Panel C
    svg_parts.append(f"""
      <g transform="translate(12, 322)">
        <circle cx="8" cy="8" r="4.5" fill="{OKABE_BLUE}" />
        <text x="18" y="11.5" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Assay Available &amp; Analyzed</text>
        <circle cx="160" cy="8" r="2" fill="{COLOR_BORDER_HAIRLINE}" />
        <text x="170" y="11.5" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Not Available / Excluded</text>
        <text x="390" y="11.5" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Matrix spans all 15 core multi-omic cohort assets</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: PANEL D — Cellular & Genomic Landscape
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 05: PANEL D — CELLULAR & GENOMIC DATA LANDSCAPE -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Cellular_Genomic_Landscape">

    <!-- Sub-block 1: Single-Cell Lineage Hierarchy (Top) -->
    <g id="grp-lineage-breakdown" transform="translate(729, 528)">
      <rect width="637" height="162" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Single-Cell Reference Cellular Compartments (40,002 Cells)</text>
      <text x="540" y="18" font-size="8" font-weight="700" fill="{OKABE_BLUE}">6 Major Lineages</text>
      <line x1="10" y1="24" x2="627" y2="24" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Horizontal Stacked Lineage Bar -->
      <g transform="translate(12, 34)">
        <rect x="0" y="0" width="192" height="20" fill="{OKABE_BLUE}" />
        <rect x="194" y="0" width="140" height="20" fill="{OKABE_SKY_BLUE}" />
        <rect x="336" y="0" width="114" height="20" fill="{OKABE_ORANGE}" />
        <rect x="452" y="0" width="74" height="20" fill="{OKABE_BLUISH_GREEN}" />
        <rect x="528" y="0" width="48" height="20" fill="{COLOR_TEXT_HAIRLINE}" />
        <rect x="578" y="0" width="35" height="20" fill="{OKABE_VERMILION}" />

        <text x="96" y="14" font-size="8" font-weight="700" fill="#FFF" text-anchor="middle">CD8+ T (31.4%)</text>
        <text x="264" y="14" font-size="8" font-weight="700" fill="#FFF" text-anchor="middle">CD4+ T (22.8%)</text>
        <text x="393" y="14" font-size="8" font-weight="700" fill="#FFF" text-anchor="middle">Myeloid (18.6%)</text>
        <text x="489" y="14" font-size="7.5" font-weight="700" fill="#FFF" text-anchor="middle">B (12.1%)</text>
        <text x="552" y="14" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Str</text>
        <text x="595" y="14" font-size="6.5" font-weight="700" fill="#FFF" text-anchor="middle">Mal</text>
      </g>

      <!-- Detailed Cell Sub-types -->
      <g transform="translate(12, 62)">
        <g transform="translate(0, 0)">
          <rect width="192" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="13" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">CD8+ T Cell Lineage (12,560 c)</text>
          <text x="6" y="24" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Naive, Memory / Stem-like (TCF7, IL7R)</text>
          <text x="6" y="34" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Effector, Exhausted (PDCD1, HAVCR2)</text>
        </g>

        <g transform="translate(196, 0)">
          <rect width="138" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="13" font-size="7.5" font-weight="700" fill="{OKABE_SKY_BLUE}">CD4+ T Cells (9,120 c)</text>
          <text x="6" y="24" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Helper (RPS/RPL active)</text>
          <text x="6" y="34" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Regulatory Treg (FOXP3)</text>
        </g>

        <g transform="translate(338, 0)">
          <rect width="134" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="13" font-size="7.5" font-weight="700" fill="{OKABE_ORANGE}">Myeloid Lineage (7,440 c)</text>
          <text x="6" y="24" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• M1 Pro-inflammatory</text>
          <text x="6" y="34" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• M2/MDSC Suppressive</text>
        </g>

        <g transform="translate(476, 0)">
          <rect width="137" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="13" font-size="7.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}">B / Plasma / Stromal</text>
          <text x="6" y="24" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Naive B (MS4A1, CD19)</text>
          <text x="6" y="34" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Fibroblasts, Endothelial</text>
        </g>
      </g>

      <!-- Sade-Feldman CD45+ Immune Distribution -->
      <g transform="translate(12, 108)">
        <rect width="613" height="42" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Matched Clinical Single-Cell Cohort (Sade-Feldman, n = 51 biopsies, 16,291 CD45+ cells):</text>
        <text x="8" y="26" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• 12 Validated Immune States: 5 CD8+ T (Tem/Trm clusters), CD4+ T, Treg, Naive B, Plasma cells, Macrophages, pDC, Monocytes</text>
        <text x="8" y="36" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Longitudinal Biopsy Partition: 20 Baseline Pre-treatment biopsies vs. 31 On-treatment / Post-progression biopsies</text>
      </g>
    </g>

    <!-- Sub-block 2: Genomic Burden & Mutation Spectrum (Bottom) -->
    <g id="grp-genomic-landscape" transform="translate(729, 698)">
      <rect width="637" height="170" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Somatic Genomic Landscape, Subclonal VAF &amp; TMB Distributions</text>
      <text x="540" y="18" font-size="8" font-weight="700" fill="{OKABE_ORANGE}">DNA &amp; Exome</text>
      <line x1="10" y1="24" x2="627" y2="24" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Indication TMB Stack -->
      <g transform="translate(12, 34)">
        <rect width="298" height="126" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Median DNA-TMB by Tumor Indication</text>
        
        <!-- Melanoma -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Melanoma (SKCM)</text>
          <rect x="105" y="2" width="130" height="9" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="105" y="2" width="98" height="9" fill="{OKABE_ORANGE}" />
          <text x="245" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">15.3 mut/Mb</text>
        </g>
        <!-- Bladder -->
        <g transform="translate(10, 42)">
          <text x="0" y="10" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Bladder (BLCA)</text>
          <rect x="105" y="2" width="130" height="9" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="105" y="2" width="62" height="9" fill="{OKABE_SKY_BLUE}" />
          <text x="245" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">9.8 mut/Mb</text>
        </g>
        <!-- Renal -->
        <g transform="translate(10, 60)">
          <text x="0" y="10" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Renal Cell (KIRC)</text>
          <rect x="105" y="2" width="130" height="9" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="105" y="2" width="26" height="9" fill="{OKABE_BLUISH_GREEN}" />
          <text x="245" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">4.1 mut/Mb</text>
        </g>
        <!-- Breast -->
        <g transform="translate(10, 78)">
          <text x="0" y="10" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Breast (TNBC)</text>
          <rect x="105" y="2" width="130" height="9" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="105" y="2" width="24" height="9" fill="{OKABE_REDDISH_PURPLE}" />
          <text x="245" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3.8 mut/Mb</text>
        </g>
        <!-- Pancreatic -->
        <g transform="translate(10, 96)">
          <text x="0" y="10" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Pancreatic (PDAC)</text>
          <rect x="105" y="2" width="130" height="9" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="105" y="2" width="16" height="9" fill="{COLOR_TEXT_MUTED}" />
          <text x="245" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2.6 mut/Mb</text>
        </g>
      </g>

      <!-- VAF & Driver Alterations Sub-card -->
      <g transform="translate(322, 34)">
        <rect width="303" height="126" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        
        <!-- VAF Spectrum Box -->
        <g transform="translate(10, 10)">
          <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Subclonal VAF Spectrum Architecture</text>
          <text x="0" y="22" font-size="7.2" fill="{COLOR_TEXT_MUTED}">• Allele Frequency Evaluation Range: 0.01 – 0.50</text>
          <text x="0" y="34" font-size="7.2" fill="{COLOR_TEXT_MUTED}">• Clonal Trunk Peaks: VAF ~0.25 – 0.45</text>
          <text x="0" y="46" font-size="7.2" fill="{COLOR_TEXT_MUTED}">• Subclonal Branch Tail: VAF &lt; 0.10 (immune editing)</text>
        </g>

        <!-- Driver Alteration Highlight: TGM6 -->
        <g transform="translate(10, 68)">
          <rect width="283" height="48" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_VERMILION}">TGM6 Alteration (OR = 27.38, p = 8.72 × 10⁻⁶)</text>
          <text x="8" y="30" font-size="7.2" fill="{COLOR_TEXT_SECONDARY}">• Mutant Cohort Response Rate: 92.9% (13/14 Responders)</text>
          <text x="8" y="42" font-size="7.2" fill="{COLOR_TEXT_MUTED}">• Wildtype Response Rate: 32.2% (66/205) | Melanoma Freq: 6.4%</text>
        </g>
      </g>
    </g>
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate pure dataset compendium SVG."""
    output_dir = Path("article/figures/dataset_overview")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "figure_dataset_only_overview.svg"
    svg_content = render_svg()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Nature Methods dataset compendium SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
