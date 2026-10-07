#!/usr/bin/env python3
"""
Generate Figure 1: Comprehensive Multi-Omic Immunotherapy Dataset Compendium,
Single-Cell References, and Analytical Deconvolution Framework.

Implements the authentic Nature Methods / Nature editorial art style:
- Double-column landscape (1400 px width, 880 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional elements
- Authentic editorial letter badges (a, b, c, d)
- Formal mathematical formulas in serif typography
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
    FONT_SANS, FONT_SERIF_MATH
)

HEIGHT = 880


def render_svg() -> str:
    """Declaratively render the complete Nature Methods Figure 1 SVG DOM."""
    svg_parts: list[str] = []

    # 1. XML Header and SVG root with Inkscape namespaces
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure1-study-overview-nature"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <defs>
    <!-- Minimalist Hairline Arrow Markers -->
    <marker id="arr-slate" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
    </marker>
    <marker id="arr-blue" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
  </defs>
""")

    # ----------------------------------------------------
    # LAYER 0: Canvas Background & Panel Containers
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND & PANEL CONTAINERS -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />

    <!-- Panel A Container Card (Top-Left) -->
    <rect id="card-panel-a" x="20" y="68" width="665" height="392"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel B Container Card (Top-Right) -->
    <rect id="card-panel-b" x="715" y="68" width="665" height="392"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel C Container Card (Bottom-Left) -->
    <rect id="card-panel-c" x="20" y="474" width="665" height="388"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel D Container Card (Bottom-Right) -->
    <rect id="card-panel-d" x="715" y="474" width="665" height="388"
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
      Figure 1: Multi-Omic Immunotherapy Dataset Compendium, Single-Cell References, and Analytical Deconvolution Framework
    </text>
    <!-- Subtitle Summary -->
    <text id="txt-subtitle" x="20" y="48" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Pan-cancer single-cell atlas (40,002 cells), 9 clinical trials (n = 1,097), 5 TCGA baseline cohorts (n = 2,932), and matched clinical deconvolution benchmarking
    </text>

    <!-- Top Total Stat Box -->
    <g id="grp-header-stat" transform="translate(1085, 14)">
      <rect x="0" y="0" width="295" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="147" y="18" font-size="9" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        Total Compendium: &gt;4,080 Multi-Omic Biopsies
      </text>
    </g>

    <!-- Panel Badges [a], [b], [c], [d] -->
    <g id="badge-panel-a" transform="translate(34, 86)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">a</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Multi-Cohort Dataset Ecosystem &amp; Discovery Compendium</text>
    </g>

    <g id="badge-panel-b" transform="translate(729, 86)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">b</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Multi-Omic Profiling Modalities &amp; Clinical Annotations</text>
    </g>

    <g id="badge-panel-c" transform="translate(34, 492)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">c</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Unified Analytical Framework &amp; Simulation Suite</text>
    </g>

    <g id="badge-panel-d" transform="translate(729, 492)">
      <text x="0" y="14" font-size="14" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">d</text>
      <text x="18" y="14" font-size="11" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Biomarker Concordance, Vulnerabilities &amp; Generalization</text>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: PANEL A — Multi-Cohort Dataset Ecosystem
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: PANEL A — MULTI-COHORT DATASET ECOSYSTEM -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Dataset_Ecosystem">

    <!-- Card 1: Single-Cell Reference & Clinical Discovery -->
    <g id="grp-card-scrna" transform="translate(34, 114)">
      <rect width="198" height="332" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Single-Cell References (40k cells)</text>
      <line x1="10" y1="26" x2="188" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Pan-cancer reference details -->
      <g transform="translate(8, 34)">
        <rect width="182" height="114" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Pan-Cancer Reference</text>
        <text x="10" y="30" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Total: 40,002 single cells</text>
        <text x="10" y="44" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 6 Major Lineages, 22 Types</text>
        <text x="10" y="58" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 58 Silhouette K-means States</text>
        
        <!-- Lineage glyphs -->
        <g transform="translate(10, 68)">
          <rect x="0" y="0" width="24" height="14" fill="{OKABE_BLUE}" /><text x="12" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">CD8</text>
          <rect x="27" y="0" width="24" height="14" fill="{OKABE_SKY_BLUE}" /><text x="39" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">CD4</text>
          <rect x="54" y="0" width="24" height="14" fill="{OKABE_ORANGE}" /><text x="66" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Mye</text>
          <rect x="81" y="0" width="24" height="14" fill="{OKABE_BLUISH_GREEN}" /><text x="93" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">B</text>
          <rect x="108" y="0" width="24" height="14" fill="{COLOR_TEXT_MUTED}" /><text x="120" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Str</text>
          <rect x="135" y="0" width="23" height="14" fill="{OKABE_VERMILION}" /><text x="146.5" y="10.5" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Mal</text>
          <text x="75" y="27" font-size="6.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Selective Inf. p &lt; 10⁻¹⁶</text>
        </g>
      </g>

      <!-- Matched clinical cohort: Sade-Feldman -->
      <g transform="translate(8, 158)">
        <rect width="182" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Sade-Feldman (Cell 2018)</text>
        <text x="10" y="30" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Matched Melanoma Validation</text>
        <text x="10" y="44" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• n = 51 biopsies (32 patients)</text>
        <text x="10" y="58" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 16,291 CD45+ cells (12 states)</text>
        <text x="10" y="72" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Smart-seq2: 20 Pre- / 31 Post-ICB</text>
        <rect x="10" y="78" width="162" height="14" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="91" y="88.5" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Ground Truth Milo DA Benchmark</text>
      </g>

      <!-- Additional References -->
      <g transform="translate(8, 262)">
        <rect width="182" height="60" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Supporting Single-Cell Atlases</text>
        <text x="10" y="31" font-size="8" fill="{COLOR_TEXT_MUTED}">• Jerby-Arnon 2018: 7,186 c (GSE115978)</text>
        <text x="10" y="45" font-size="8" fill="{COLOR_TEXT_MUTED}">• Li 2019 (GSE123139) | Ma 2019 (Liver)</text>
      </g>
    </g>

    <!-- Card 2: Clinical Immunotherapy Trial Cohorts (n = 1,097) -->
    <g id="grp-card-trials" transform="translate(242, 114)">
      <rect width="242" height="332" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Clinical ICI Trial Cohorts (n = 1,097)</text>
      <line x1="10" y1="26" x2="232" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Melanoma -->
      <g transform="translate(8, 34)">
        <rect width="226" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_ORANGE}">Melanoma (n = 347 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Hugo (n=27, Pembro) | Riaz (n=107, Nivo)</text>
        <text x="8" y="46" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Liu (n=122, Anti-PD1) | Gide (n=91, Combo)</text>
      </g>

      <!-- Bladder -->
      <g transform="translate(8, 98)">
        <rect width="226" height="46" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_SKY_BLUE}">Bladder (n = 347 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Rosenberg / IMvigor210 (Anti-PD-L1 Atezo)</text>
      </g>

      <!-- Renal -->
      <g transform="translate(8, 150)">
        <rect width="226" height="46" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Renal Cell (n = 279 patients)</text>
        <text x="8" y="32" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• McDermott / IMmotion150 (263) | Choueiri (16)</text>
      </g>

      <!-- Pancreatic & Breast -->
      <g transform="translate(8, 202)">
        <rect width="226" height="56" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{OKABE_REDDISH_PURPLE}">PDAC / BC (n = 124 patients)</text>
        <text x="8" y="30" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Padron (n=93, PDAC, Nivo + Sotigalimab)</text>
        <text x="8" y="44" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">• Anders (n=31, TNBC, Atezo + Chemo)</text>
      </g>

      <!-- Extended Compendium -->
      <g transform="translate(8, 264)">
        <rect width="226" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">+ Extended Validation Cohorts</text>
        <text x="8" y="30" font-size="7.2" fill="{COLOR_TEXT_MUTED}">VanAllen (110) | Snyder (25) | Rose (50) | Ravi (44)</text>
        <text x="8" y="44" font-size="7.2" fill="{COLOR_TEXT_MUTED}">Prat (65) | Freeman (104) | Auslander | Lauss</text>
      </g>
    </g>

    <!-- Card 3: TCGA Baseline Reference Cohorts (n = 2,932) -->
    <g id="grp-card-tcga" transform="translate(492, 114)">
      <rect width="182" height="332" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. TCGA Baselines (n = 2,932)</text>
      <line x1="10" y1="26" x2="172" y2="26" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- 5 Primary Types Stack -->
      <g transform="translate(8, 34)">
        <rect width="166" height="154" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">5 Untreated Primary Types</text>
        
        <g transform="translate(10, 28)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-SKCM</text>
          <rect x="70" y="2" width="52" height="11" fill="{OKABE_ORANGE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">471</text>
        </g>
        <g transform="translate(10, 46)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-BLCA</text>
          <rect x="70" y="2" width="46" height="11" fill="{OKABE_SKY_BLUE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">414</text>
        </g>
        <g transform="translate(10, 64)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-BRCA</text>
          <rect x="70" y="2" width="85" height="11" fill="{OKABE_REDDISH_PURPLE}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1,102</text>
        </g>
        <g transform="translate(10, 82)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-KIRC</text>
          <rect x="70" y="2" width="60" height="11" fill="{OKABE_BLUISH_GREEN}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">538</text>
        </g>
        <g transform="translate(10, 100)">
          <text x="0" y="11" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">TCGA-PAAD</text>
          <rect x="70" y="2" width="20" height="11" fill="{COLOR_TEXT_MUTED}" opacity="0.8" />
          <text x="128" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">178</text>
        </g>
        <text x="10" y="138" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Total: 2,703–2,932</text>
      </g>

      <!-- Microenvironment baseline -->
      <g transform="translate(8, 196)">
        <rect width="166" height="126" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="18" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Baseline Microenvironment</text>
        <text x="10" y="34" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• K = 100 Centroids</text>
        <text x="10" y="48" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• 58 Deconv States</text>
        <text x="10" y="62" font-size="8" fill="{COLOR_TEXT_SECONDARY}">• Primary Organ Space</text>
        <rect x="10" y="72" width="146" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="14" y="87" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Unperturbed Baseline</text>
        <text x="14" y="100" font-size="7" fill="{COLOR_TEXT_MUTED}">Anchors metastatic remodeling</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: PANEL B — Multi-Omic Profiling Modalities
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 03: PANEL B — MULTI-OMIC PROFILING MODALITIES -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Profiling_Modalities">

    <!-- Track 1: RNA -->
    <g id="grp-track-rna" transform="translate(729, 114)">
      <rect width="637" height="74" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Transcriptomic Modalities &amp; Expression Matrices</text>
      
      <g transform="translate(14, 30)">
        <rect x="0" y="0" width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Single-Cell RNA-seq (Smart-seq2 &amp; 10x)</text>
        <text x="8" y="26" font-size="7" fill="{COLOR_TEXT_MUTED}">UMI counts, kNN manifold, cell-state DGE signatures</text>

        <rect x="312" y="0" width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="320" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Bulk Tumor RNA-seq (n = 1,097 Trials + 2,932 TCGA)</text>
        <text x="320" y="26" font-size="7" fill="{COLOR_TEXT_MUTED}">RSEM TPM / FPKM matrices decompounded into cell fractions</text>
      </g>
    </g>

    <!-- Track 2: DNA / Somatic Genomics -->
    <g id="grp-track-genomics" transform="translate(729, 196)">
      <rect width="637" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Somatic Genomics, Subclonal VAF Spectrum &amp; Predictive Drivers</text>
      
      <g transform="translate(14, 30)">
        <!-- VAF Spectrum Card -->
        <g transform="translate(0, 0)">
          <rect x="0" y="0" width="195" height="54" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Subclonal VAF Spectrum</text>
          <text x="8" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">• Allele frequencies: 0.01 – 0.50</text>
          <text x="8" y="39" font-size="7" fill="{COLOR_TEXT_MUTED}">• Unsupervised KDE Valley detection</text>
          <text x="8" y="49" font-size="7" fill="{COLOR_TEXT_MUTED}">• TMB Reliability Score (TRS)</text>
        </g>

        <!-- Driver Alteration: TGM6 -->
        <g transform="translate(205, 0)">
          <rect x="0" y="0" width="200" height="54" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{OKABE_VERMILION}">TGM6 Alteration (OR = 27.38)</text>
          <text x="8" y="28" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Responders: 92.9% (13/14 mutated)</text>
          <text x="8" y="42" font-size="7" fill="{COLOR_TEXT_MUTED}">• Wildtype: 32.2% (p = 8.72 × 10⁻⁶)</text>
        </g>

        <!-- Predictive Modeling -->
        <g transform="translate(415, 0)">
          <rect x="0" y="0" width="192" height="54" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Predictive Modeling</text>
          <text x="8" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">• Random Forest, Adaline, MLP</text>
          <text x="8" y="39" font-size="7" fill="{COLOR_TEXT_MUTED}">• Stratified 5-Fold Cross-Validation</text>
          <text x="8" y="49" font-size="7" fill="{COLOR_TEXT_MUTED}">• ANOVA Feature Selection (AUC .74)</text>
        </g>
      </g>
    </g>

    <!-- Track 3: Standardized Clinical Endpoints -->
    <g id="grp-track-clinical" transform="translate(729, 300)">
      <rect width="637" height="66" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. Harmonized Clinical Response Endpoints &amp; Criteria</text>
      
      <g transform="translate(14, 28)">
        <g transform="translate(0, 0)">
          <rect x="0" y="0" width="144" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="18" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Responders: CR / PR</text>
        </g>

        <g transform="translate(154, 0)">
          <rect x="0" y="0" width="150" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="18" font-size="8" font-weight="700" fill="{OKABE_VERMILION}">Non-Responders: SD / PD</text>
        </g>

        <g transform="translate(314, 0)">
          <rect x="0" y="0" width="293" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="18" font-size="7.5" fill="{COLOR_TEXT_SECONDARY}">Progression-Free Survival (PFS) &amp; Overall Survival (OS)</text>
        </g>
      </g>
    </g>

    <!-- Track 4: Longitudinal Timing & Regimens -->
    <g id="grp-track-timing" transform="translate(729, 374)">
      <rect width="637" height="72" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">4. Biopsy Sampling Timing &amp; Therapeutic Regimens</text>
      
      <g transform="translate(14, 30)">
        <g transform="translate(0, 0)">
          <rect x="0" y="0" width="144" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Pre-Treatment Baseline</text>
          <text x="8" y="24" font-size="7" fill="{COLOR_TEXT_MUTED}">Untreated / Ipi-Naive (n = 852)</text>
        </g>

        <g transform="translate(154, 0)">
          <rect x="0" y="0" width="150" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">On-Treatment / Post-Prog.</text>
          <text x="8" y="24" font-size="7" fill="{COLOR_TEXT_MUTED}">Pre-treated / Resistant (n = 245)</text>
        </g>

        <g transform="translate(314, 0)">
          <rect x="0" y="0" width="293" height="30" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="18" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Regimens: Anti-PD-1/L1, Anti-CTLA-4, Nivo+Ipi, Combo+Chemo</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: PANEL C — Unified Analytical Framework
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: PANEL C — UNIFIED ANALYTICAL FRAMEWORK -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Analytical_Pipeline">

    <!-- Stream 1: Single-Cell Milo Graph DA -->
    <g id="grp-pipeline-milo" transform="translate(34, 520)">
      <rect width="637" height="74" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stream 1: Single-Cell Graph Differential Abundance (Milo)</text>
      
      <g transform="translate(14, 28)">
        <g transform="translate(0, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Continuous kNN Neighborhood Graph (k = 15)</text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Non-parametric testing across phenotypic manifolds without clustering</text>
        </g>

        <g transform="translate(312, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">Negative Binomial GLM: log(𝔼[N<tspan baseline-shift="sub" font-size="75%">v,j</tspan>]) = β<tspan baseline-shift="sub" font-size="75%">0</tspan> + β<tspan baseline-shift="sub" font-size="75%">v</tspan> Y<tspan baseline-shift="sub" font-size="75%">j</tspan></text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Cell-state neighborhood effect size: logFC<tspan baseline-shift="sub" font-size="75%">k</tspan> = median(β<tspan baseline-shift="sub" font-size="75%">v</tspan>)</text>
        </g>
      </g>
    </g>

    <!-- Stream 2: Bulk Deconvolution (InstaPrism / BayesPrism) -->
    <g id="grp-pipeline-deconv" transform="translate(34, 602)">
      <rect width="637" height="74" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stream 2: High-Resolution Bulk Deconvolution (InstaPrism / BayesPrism)</text>
      
      <g transform="translate(14, 28)">
        <g transform="translate(0, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">Simplex Inversion: b<tspan baseline-shift="sub" font-size="75%">j</tspan> ≈ Φ · f<tspan baseline-shift="sub" font-size="75%">j</tspan><tspan baseline-shift="super" font-size="75%">mRNA</tspan> (58 states)</text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Reference signature matrix Φ constructed from scRNA-seq centroids</text>
        </g>

        <g transform="translate(312, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">Point-Biserial Logit Effect: β̂<tspan baseline-shift="sub" font-size="75%">k</tspan> = 2 r<tspan baseline-shift="sub" font-size="75%">pb</tspan> / √(1 - r<tspan baseline-shift="sub" font-size="75%">pb</tspan>²)</text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Standardized association with RECIST clinical response outcome</text>
        </g>
      </g>
    </g>

    <!-- Stream 3: Multi-Replicate Distortion Simulation Suite -->
    <g id="grp-pipeline-stress" transform="translate(34, 684)">
      <rect width="637" height="84" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stream 3: Multi-Replicate Perturbation Suite (N = 5,485 Simulations)</text>
      
      <g transform="translate(14, 28)">
        <!-- Sub-card 1: 8 Distortion Modes -->
        <g transform="translate(0, 0)">
          <rect width="195" height="44" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">8 Distortion Modes (240 runs)</text>
          <text x="8" y="26" font-size="7" fill="{COLOR_TEXT_MUTED}">• Cell Size Asymmetry (1× – 50×)</text>
          <text x="8" y="37" font-size="7" fill="{COLOR_TEXT_MUTED}">• Marker Dysregulation, Ambient Soup</text>
        </g>

        <!-- Sub-card 2: 4 Compound Regimes -->
        <g transform="translate(205, 0)">
          <rect width="200" height="44" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">4 Compound Regimes (120 runs)</text>
          <text x="8" y="26" font-size="7" fill="{COLOR_TEXT_MUTED}">• Core Needle Biopsy (Size + Sparsity)</text>
          <text x="8" y="37" font-size="7" fill="{COLOR_TEXT_MUTED}">• Inflamed TME, Triple Breakdown</text>
        </g>

        <!-- Sub-card 3: 2D Interaction Surface -->
        <g transform="translate(415, 0)">
          <rect width="192" height="44" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2D Surface (125 runs)</text>
          <text x="8" y="26" font-size="7" fill="{COLOR_TEXT_MUTED}">• 5×5 Factorial: Size × Activation</text>
          <text x="8" y="37" font-size="7" fill="{COLOR_TEXT_MUTED}">• Directional Sign Inversion Maps</text>
        </g>
      </g>
    </g>

    <!-- Stream 4: Subclonal VAF Dynamics & Joint Cutoff Optimization -->
    <g id="grp-pipeline-vaf" transform="translate(34, 776)">
      <rect width="637" height="74" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stream 4: Subclonal VAF &amp; Joint TMB Threshold Optimization</text>
      
      <g transform="translate(14, 28)">
        <g transform="translate(0, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Unsupervised KDE Valley Detection</text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Automated inflection point discovery without outcome supervision</text>
        </g>

        <g transform="translate(312, 0)">
          <rect width="295" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">Regularized Joint Cutoff: f<tspan baseline-shift="sub" font-size="75%">reg</tspan> = f - λ (Δlog t<tspan baseline-shift="sub" font-size="75%">vaf</tspan>² + Δlog t<tspan baseline-shift="sub" font-size="75%">tmb</tspan>²)</text>
          <text x="8" y="25" font-size="7" fill="{COLOR_TEXT_MUTED}">Penalizes excessive baseline divergence (λ = 0.20, t<tspan baseline-shift="sub" font-size="75%">vaf</tspan>=0.05, t<tspan baseline-shift="sub" font-size="75%">tmb</tspan>=10)</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: PANEL D — Biomarker Synthesis & Benchmarks
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 05: PANEL D — BIOMARKER SYNTHESIS & BENCHMARKS -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Biomarker_Synthesis">

    <!-- Sub-Panel 1: Diagnostic Concordance Quadrants (Sade-Feldman n=51) -->
    <g id="grp-diagnostic-quadrants" transform="translate(729, 520)">
      <rect width="300" height="330" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Milo vs. Deconv. Concordance (n = 51)</text>
      <text x="240" y="18" font-size="8" font-weight="700" fill="{OKABE_BLUE}">ρ = 0.853</text>
      <line x1="12" y1="24" x2="288" y2="24" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Quadrant Coordinate System -->
      <g transform="translate(30, 36)">
        <rect x="0" y="0" width="240" height="224" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <line x1="120" y1="0" x2="120" y2="224" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" stroke-dasharray="3,3" />
        <line x1="0" y1="112" x2="240" y2="112" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" stroke-dasharray="3,3" />

        <!-- Axis Labels -->
        <text x="120" y="238" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Milo Neighborhood logFC →</text>
        <text x="-112" y="-12" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Deconv. β (Response) →</text>

        <!-- Quadrant Watermarks -->
        <text x="180" y="75" font-size="6.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" opacity="0.3" text-anchor="middle">CONCORDANT RESPONDER</text>
        <text x="60" y="145" font-size="6.5" font-weight="700" fill="{OKABE_BLUE}" opacity="0.3" text-anchor="middle">CONCORDANT NON-RESP.</text>
        <text x="180" y="200" font-size="6.5" font-weight="700" fill="{OKABE_VERMILION}" opacity="0.3" text-anchor="middle">DISCORDANT</text>

        <!-- Plotted Biomarker Cell States -->
        <!-- 04_Naive B cells -->
        <circle cx="205" cy="35" r="6" fill="{OKABE_BLUISH_GREEN}" />
        <text x="205" y="23" font-size="7.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">04_Naive B</text>
        <text x="205" y="48" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">β = +2.07</text>

        <!-- 08_Treg -->
        <circle cx="175" cy="106" r="4.5" fill="{COLOR_TEXT_HAIRLINE}" />
        <text x="175" y="98" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">08_Treg</text>

        <!-- 05_Macrophages -->
        <circle cx="45" cy="180" r="6" fill="{OKABE_BLUE}" />
        <text x="45" y="168" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">05_Macro</text>
        <text x="45" y="194" font-size="6.5" font-weight="600" fill="{OKABE_BLUE}" text-anchor="middle">β = -1.56</text>

        <!-- 07_Tem/Trm Cytotoxic T -->
        <circle cx="50" cy="132" r="5" fill="{OKABE_BLUE}" />
        <text x="62" y="135" font-size="6.5" fill="{OKABE_BLUE}">07_Tem/Trm (β = -1.05)</text>

        <!-- 10_pDC -->
        <circle cx="70" cy="116" r="4" fill="{OKABE_BLUE}" />
        <text x="82" y="119" font-size="6.5" fill="{OKABE_BLUE}">10_pDC</text>

        <!-- 01_Tem/Trm (Discordant) -->
        <circle cx="185" cy="142" r="6" fill="{OKABE_VERMILION}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
        <text x="185" y="157" font-size="7" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">01_Tem (Discordant)</text>
      </g>

      <!-- Bottom Concordance Stat Summary -->
      <g transform="translate(14, 284)">
        <rect width="272" height="36" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Concordance: Spearman ρ = 0.853 (p = 4.2 × 10⁻⁴)</text>
        <text x="8" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">Sign agreement: 75.0% | Pre-ICB: ρ = 0.818 | Post-ICB: ρ = 0.783</text>
      </g>
    </g>

    <!-- Sub-Panel 2: Deconvolution Failure Modes & Vulnerabilities -->
    <g id="grp-vulnerability-gauges" transform="translate(1045, 520)">
      <rect width="321" height="158" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Primary Deconvolution Failure Modes</text>
      <line x1="12" y1="24" x2="309" y2="24" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Failure Mode 1: Marker Dysregulation -->
      <g transform="translate(12, 34)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Patient Marker Dysregulation (Most Lethal)</text>
        <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_MUTED}">Bulk marker suppression decouples from physical cells</text>
        <rect x="0" y="25" width="297" height="12" fill="{COLOR_BORDER_HAIRLINE}" />
        <rect x="0" y="25" width="280" height="12" fill="{OKABE_VERMILION}" />
        <text x="140" y="34.5" font-size="7" font-weight="700" fill="#FFFFFF" text-anchor="middle">Spearman ρ: +0.97 → -0.04 (Sign Agr. 47.5%)</text>
      </g>

      <!-- Failure Mode 2: Cell Size Asymmetry -->
      <g transform="translate(12, 78)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Cell Size / mRNA Asymmetry (1× – 50×)</text>
        <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_MUTED}">Disproportionate mRNA yield drives rank collapse</text>
        <rect x="0" y="25" width="297" height="12" fill="{COLOR_BORDER_HAIRLINE}" />
        <rect x="0" y="25" width="210" height="12" fill="{OKABE_ORANGE}" />
        <text x="105" y="34.5" font-size="7" font-weight="700" fill="#FFFFFF" text-anchor="middle">Spearman ρ: +0.82 → +0.25 (Δρ = 0.57)</text>
      </g>

      <!-- Failure Mode 3: Compound Clinical Core Needle Biopsy -->
      <g transform="translate(12, 122)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. Compound Core Needle Biopsy (Most Vulnerable Regime)</text>
        <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Size disparity + Biopsy sparsity + Unmodeled ghost tumor: ρ = 0.27</text>
      </g>
    </g>

    <!-- Sub-Panel 3: TCGA Baseline Generalization Gap -->
    <g id="grp-generalization-map" transform="translate(1045, 692)">
      <rect width="321" height="158" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="12" y="18" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Metastatic Generalization Gap (K = 100 Centroids)</text>
      <text x="240" y="18" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}">Acc: 26.3%</text>
      <line x1="12" y1="24" x2="309" y2="24" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Performance Table Grid -->
      <g transform="translate(12, 34)">
        <rect width="297" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="12.5" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Primary Type</text>
        <text x="90" y="12.5" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Precision</text>
        <text x="160" y="12.5" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Recall</text>
        <text x="230" y="12.5" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">F1-Score</text>

        <text x="10" y="32" font-size="7.5" fill="{COLOR_TEXT_PRIMARY}">SKCM (Melanoma)</text>
        <text x="90" y="32" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">0.89</text>
        <text x="160" y="32" font-size="7.5" fill="{COLOR_TEXT_MUTED}">0.52</text>
        <text x="230" y="32" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">0.66</text>

        <text x="10" y="48" font-size="7.5" fill="{COLOR_TEXT_PRIMARY}">KIRC (Renal)</text>
        <text x="90" y="48" font-size="7.5" fill="{COLOR_TEXT_MUTED}">0.46</text>
        <text x="160" y="48" font-size="7.5" fill="{COLOR_TEXT_MUTED}">0.28</text>
        <text x="230" y="48" font-size="7.5" fill="{COLOR_TEXT_MUTED}">0.34</text>

        <text x="10" y="64" font-size="7.5" fill="{COLOR_TEXT_PRIMARY}">BLCA (Bladder)</text>
        <text x="90" y="64" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>
        <text x="160" y="64" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>
        <text x="230" y="64" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>

        <text x="10" y="80" font-size="7.5" fill="{COLOR_TEXT_PRIMARY}">PAAD (Pancreas)</text>
        <text x="90" y="80" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>
        <text x="160" y="80" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>
        <text x="230" y="80" font-size="7.5" fill="{OKABE_VERMILION}">0.00</text>
      </g>

      <!-- Explanation Callout -->
      <g transform="translate(12, 120)">
        <rect width="297" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="12" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">High Rank Order (ρ = 0.72–0.86) vs. Remodeling Gap</text>
        <text x="8" y="22" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Bladder/PDAC zero recall caused by stromal overgrowth &amp; chemo remodeling</text>
      </g>
    </g>
  </g>

  <!-- LAYER 06: FLOW CONNECTORS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Flow_Connectors">
    <!-- Connector: Panel A to Panel B -->
    <path d="M 685 150 L 710 150" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 685 240 L 710 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 685 330 L 710 330" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector: Panel C to Panel D -->
    <path d="M 685 556 L 725 556" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1030 600 L 1042 600" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1030 770 L 1042 770" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate Figure 1 vector SVG file."""
    output_dir = Path("article/figures/dataset_overview")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "figure1_study_cohorts_overview.svg"
    svg_content = render_svg()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Figure 1 vector SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
