#!/usr/bin/env python3
"""
Generate Nature Methods Reference Vector Figure: How BayesPrism Works.

Implements the authentic Nature Methods / Nature editorial art style:
- Double-column landscape (180 mm / 1400 px width, 630 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional elements
- Formal display equations with serif mathematics (Georgia/Times) and formal numbering (1), (2), (3)
- Formal Bayesian plate diagram with double circles for observed variables
- Real coordinate axes with outward tick marks on trace plots
"""

from pathlib import Path
from nature_style_config import (
    WIDTH, HEIGHT,
    COLOR_CANVAS_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)


def render_nature_methods_bayesprism() -> str:
    """Render the BayesPrism workflow in pure Nature Methods editorial style."""
    svg_parts: list[str] = []

    # 1. XML Header and Root SVG
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-bayesprism-nature-methods"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <defs>
    <!-- Minimalist Hairline Arrow Markers -->
    <marker id="arr-slate" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
    </marker>
    <marker id="arr-blue" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
    <marker id="arr-purple" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_REDDISH_PURPLE}" />
    </marker>
    <marker id="arr-green" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUISH_GREEN}" />
    </marker>
  </defs>
""")

    # ----------------------------------------------------
    # LAYER 0: Canvas Background & Outer Wireframe
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND & OUTER WIREFRAME -->
  <g inkscape:groupmode="layer" id="layer-00-canvas" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />
    <!-- Outer Hairline Frame (No drop shadow, sharp corners rx=0) -->
    <rect id="frame-outer" x="18" y="14" width="1364" height="602" fill="{COLOR_CANVAS_BG}"
          stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 1: Header, Subtitle & Reference Tag
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 01: FIGURE HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <g id="grp-header" transform="translate(36, 24)">
      <!-- Main Title -->
      <text id="txt-title" x="0" y="18" font-size="13.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.1">
        BayesPrism: Probabilistic Generative Deconvolution of Cell Types &amp; Lineage Expression
      </text>
      <!-- Subtitle -->
      <text id="txt-subtitle" x="0" y="34" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        Joint Bayesian inversion of bulk RNA-seq mixtures into cell fractions (θ) and cell-type-specific read counts (Z)
      </text>

      <!-- Academic Reference Tag on Right -->
      <g id="grp-ref-tag" transform="translate(1085, 4)">
        <rect x="0" y="0" width="240" height="24" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="120" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
          Chu et al., Nat. Cancer 3, 482–498 (2022)
        </text>
      </g>
    </g>

    <!-- Header Divider Line -->
    <line x1="18" y1="68" x2="1382" y2="68" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: Column 1 — Data Inputs & Prior Specification
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: COLUMN 1 — DATA INPUTS & PRIORS -->
  <g inkscape:groupmode="layer" id="layer-02-stage1" inkscape:label="02_Stage_1_Inputs_and_Priors">
    <g id="grp-col-1" transform="translate(36, 82)">
      <!-- Column 1 Frame (width=305, height=518, rx=0) -->
      <rect width="305" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Column Header Banner -->
      <line x1="0" y1="28" x2="305" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="305" y2="0" stroke="{OKABE_BLUE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        1. INPUT DATA &amp; PRIOR SPECIFICATION
      </text>

      <!-- Block 1: Single-Cell RNA-seq Reference -->
      <g id="col1-box-sc" transform="translate(10, 36)">
        <rect width="285" height="136" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Single-cell reference matrix X ∈ ℕ<tspan baseline-shift="super" font-size="75%">C × G</tspan></text>

        <!-- Mini Matrix Schematic -->
        <g transform="translate(10, 24)">
          <rect width="70" height="60" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <line x1="23" y1="0" x2="23" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="46" y1="0" x2="46" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="0" y1="30" x2="70" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <rect x="2" y="2" width="19" height="26" fill="{OKABE_BLUE}" opacity="0.4" />
          <rect x="25" y="32" width="19" height="26" fill="{OKABE_ORANGE}" opacity="0.4" />
          <rect x="48" y="12" width="20" height="38" fill="{OKABE_BLUISH_GREEN}" opacity="0.4" />
          <text x="35" y="70" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
          <text x="-30" y="-3" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Cells (C) →</text>
        </g>

        <!-- Annotations -->
        <g transform="translate(90, 24)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_SECONDARY}">Two-tier hierarchy:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Cell types t ∈ {{1, …, T}} (coarse)</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Cell states s ∈ {{1, …, S}} (fine)</text>
          
          <rect x="0" y="40" width="185" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="52" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Dirichlet base prior:</text>
          <text x="6" y="64" font-size="6.8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">α<tspan baseline-shift="sub" font-size="75%">t, g</tspan> = α<tspan baseline-shift="sub" font-size="75%">0</tspan> + ∑<tspan baseline-shift="sub" font-size="75%">c ∈ C_t</tspan> X<tspan baseline-shift="sub" font-size="75%">c, g</tspan></text>
          <text x="6" y="74" font-size="6" fill="{COLOR_TEXT_MUTED}">Static profile: φ<tspan baseline-shift="sub" font-size="75%">t, ·</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> = α<tspan baseline-shift="sub" font-size="75%">t, ·</tspan> / ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> α<tspan baseline-shift="sub" font-size="75%">t, g</tspan></text>
        </g>
        <text x="10" y="128" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Captures cell-type baseline expression profiles on simplex Δ<tspan baseline-shift="super" font-size="75%">G-1</tspan>.</text>
      </g>

      <!-- Block 2: Bulk RNA-seq Raw Counts -->
      <g id="col1-box-bulk" transform="translate(10, 180)">
        <rect width="285" height="136" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Bulk RNA-seq count matrix Y ∈ ℕ<tspan baseline-shift="super" font-size="75%">G × N</tspan></text>

        <!-- Mini Matrix Bulk Bars -->
        <g transform="translate(10, 24)">
          <rect width="70" height="60" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <line x1="23" y1="0" x2="23" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="46" y1="0" x2="46" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <rect x="4" y="12" width="15" height="46" fill="{OKABE_BLUE}" opacity="0.6" />
          <rect x="27" y="24" width="15" height="34" fill="{OKABE_BLUE}" opacity="0.8" />
          <rect x="50" y="6" width="15" height="52" fill="{OKABE_BLUE}" opacity="0.4" />
          <text x="35" y="70" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
          <text x="-30" y="-3" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Tumors (N) →</text>
        </g>

        <!-- Annotations -->
        <g transform="translate(90, 24)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_SECONDARY}">Input properties:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Raw integer read counts Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan></text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Library size: N<tspan baseline-shift="sub" font-size="75%">n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan></text>
          
          <rect x="0" y="40" width="185" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="52" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Malignant compartment:</text>
          <text x="6" y="64" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Accommodates tumor expression drift</text>
          <text x="6" y="74" font-size="6" fill="{COLOR_TEXT_MUTED}">TME reference treated as static baseline</text>
        </g>
        <text x="10" y="128" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Observed across N patient clinical trial biopsy cohorts.</text>
      </g>

      <!-- Block 3: Pre-deconvolution Gene Filtering -->
      <g id="col1-box-filtering" transform="translate(10, 324)">
        <rect width="285" height="182" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Pre-deconvolution signature selection:</text>

        <!-- Step 1 -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">a. Outlier gene filtering</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Removes tumor markers, ribosomal and mitochondrial RNAs</text>
          <text x="0" y="32" font-size="6.5" font-weight="600" fill="{OKABE_VERMILION}">Prevents malignant read leakage into non-malignant TME</text>
          <line x1="0" y1="38" x2="265" y2="38" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        </g>

        <!-- Step 2 -->
        <g transform="translate(10, 68)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">b. Cell-state sub-clustering</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Estimates state priors α<tspan baseline-shift="sub" font-size="75%">s, g</tspan> within coarse lineage groups</text>
          <text x="0" y="32" font-size="6.5" font-weight="600" fill="{OKABE_ORANGE}">Reduces reference condition number κ(Φ) to prevent sign flips</text>
          <line x1="0" y1="38" x2="265" y2="38" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        </g>

        <!-- Step 3 -->
        <g transform="translate(10, 112)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">c. Baseline profile normalization</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Normalizes reference centroid vector to unit simplex Δ<tspan baseline-shift="super" font-size="75%">G-1</tspan></text>
          <text x="0" y="32" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Provides initial condition for Gibbs sampler</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: Column 2 — Probabilistic Generative Model
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 03: COLUMN 2 — GENERATIVE MIXTURE MODEL -->
  <g inkscape:groupmode="layer" id="layer-03-stage2" inkscape:label="03_Stage_2_Generative_Model">
    <g id="grp-col-2" transform="translate(366, 82)">
      <!-- Column 2 Frame (width=325, height=518, rx=0) -->
      <rect width="325" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Column Header Banner -->
      <line x1="0" y1="28" x2="325" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="325" y2="0" stroke="{OKABE_ORANGE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        2. GENERATIVE MIXTURE MODEL
      </text>

      <!-- Block 1: Bulk Generative Likelihood -->
      <g id="col2-box-likelihood" transform="translate(10, 36)">
        <rect width="305" height="108" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="15" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Multinomial mixture likelihood:</text>

        <!-- Formal Equation (1) -->
        <rect x="8" y="24" width="289" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="135" y="41" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" font-family="{FONT_SERIF_MATH}">
          Y<tspan baseline-shift="sub" font-size="75%">·, n</tspan> ~ Multinomial( N<tspan baseline-shift="sub" font-size="75%">n</tspan>,  ψ<tspan baseline-shift="sub" font-size="75%">n</tspan> )
        </text>
        <text x="135" y="56" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle" font-family="{FONT_SERIF_MATH}">
          where  ψ<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>
        </text>
        <text x="282" y="48" font-size="8.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(1)</text>

        <g transform="translate(10, 74)">
          <text x="0" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Cell fraction simplex: θ<tspan baseline-shift="sub" font-size="75%">·, n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="75%">T-1</tspan> (∑<tspan baseline-shift="sub" font-size="75%">t</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> = 1)</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Lineage expression: φ<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="75%">G-1</tspan> (∑<tspan baseline-shift="sub" font-size="75%">g</tspan> φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> = 1)</text>
          <text x="0" y="30" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• N<tspan baseline-shift="sub" font-size="75%">n</tspan> = total sequenced reads in patient tumor sample n</text>
        </g>
      </g>

      <!-- Block 2: Formal Bayesian Plate Diagram -->
      <g id="col2-box-plate" transform="translate(10, 152)">
        <rect width="305" height="210" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="15" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Bayesian graphical plate diagram:</text>

        <!-- Outer Plate Rectangle: G Genes x N Samples -->
        <rect x="20" y="26" width="265" height="154" fill="none" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.75" />
        <text x="276" y="174" font-size="8" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{COLOR_TEXT_SECONDARY}" text-anchor="end">N, G, T</text>

        <!-- DAG Nodes -->
        <!-- Prior Node: alpha_tg -->
        <g transform="translate(42, 40)">
          <circle cx="18" cy="18" r="16" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" />
          <text x="18" y="22" font-size="9" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">α<tspan baseline-shift="sub" font-size="75%">tg</tspan></text>
          <text x="18" y="44" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Prior</text>
        </g>

        <!-- Arrow alpha -> phi -->
        <path d="M 78 58 L 104 58" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

        <!-- Latent Node: phi_tgn -->
        <g transform="translate(108, 40)">
          <circle cx="18" cy="18" r="16" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_ORANGE}" stroke-width="1.4" />
          <text x="18" y="22" font-size="9" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{OKABE_ORANGE}" text-anchor="middle">φ<tspan baseline-shift="sub" font-size="75%">tgn</tspan></text>
          <text x="18" y="44" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Profile</text>
        </g>

        <!-- Arrow phi -> Z -->
        <path d="M 144 58 L 170 58" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

        <!-- Latent Node: Z_tgn (Latent Counts) -->
        <g transform="translate(174, 40)">
          <circle cx="18" cy="18" r="16" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_VERMILION}" stroke-width="1.6" />
          <text x="18" y="22" font-size="9" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Z<tspan baseline-shift="sub" font-size="75%">tgn</tspan></text>
          <text x="18" y="44" font-size="6" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Latent</text>
        </g>

        <!-- Arrow Z -> Y -->
        <path d="M 210 58 L 234 58" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.4" fill="none" marker-end="url(#arr-slate)" />

        <!-- Observed Node: Y_gn (Formal Double Circle for Observed) -->
        <g transform="translate(238, 40)">
          <circle cx="18" cy="18" r="18" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="1.4" />
          <circle cx="18" cy="18" r="15" fill="none" stroke="{OKABE_BLUE}" stroke-width="0.8" />
          <text x="18" y="22" font-size="9.5" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Y<tspan baseline-shift="sub" font-size="75%">gn</tspan></text>
          <text x="18" y="46" font-size="6" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Observed</text>
        </g>

        <!-- Latent Node: theta_tn -->
        <g transform="translate(108, 106)">
          <circle cx="18" cy="18" r="16" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_BLUISH_GREEN}" stroke-width="1.4" />
          <text x="18" y="22" font-size="9" font-family="{FONT_SERIF_MATH}" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">θ<tspan baseline-shift="sub" font-size="75%">tn</tspan></text>
          <text x="18" y="44" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Fraction</text>
        </g>

        <!-- Diagonal Arrow theta -> Z -->
        <path d="M 142 116 L 176 76" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

        <text x="152" y="196" font-size="6.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Double circle indicates observed variable; single circles denote latent states</text>
      </g>

      <!-- Block 3: Latent Count Allocation Principle -->
      <g id="col2-box-allocation" transform="translate(10, 370)">
        <rect width="305" height="136" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="15" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Latent read allocation principle:</text>

        <!-- Formal Equation (2) -->
        <rect x="8" y="24" width="289" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="135" y="45" font-size="10.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" font-family="{FONT_SERIF_MATH}">
          Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>,    Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ∈ ℕ<tspan baseline-shift="sub" font-size="75%">0</tspan>
        </text>
        <text x="282" y="45" font-size="8.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(2)</text>

        <g transform="translate(10, 68)">
          <text x="0" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Preserves raw integer read counts across all genes</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Inverts bulk mixture into exact lineage read counts</text>
          <text x="0" y="32" font-size="6.8" font-weight="700" fill="{OKABE_VERMILION}">• Jointly infers both cellular composition (θ) and expression (Z)</text>
          <text x="0" y="43" font-size="6.5" fill="{COLOR_TEXT_MUTED}">• Eliminates distortion caused by ad-hoc non-linear log scaling</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: Column 3 — Two-Stage MCMC Gibbs Inference
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: COLUMN 3 — MCMC GIBBS INFERENCE -->
  <g inkscape:groupmode="layer" id="layer-04-stage3" inkscape:label="03_Stage_3_Gibbs_Inference">
    <g id="grp-col-3" transform="translate(716, 82)">
      <!-- Column 3 Frame (width=325, height=518, rx=0) -->
      <rect width="325" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Column Header Banner -->
      <line x1="0" y1="28" x2="325" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="325" y2="0" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        3. TWO-STAGE MCMC GIBBS INFERENCE
      </text>

      <!-- Block 1: Stage 1 Initialization -->
      <g id="col3-box-stage1" transform="translate(10, 36)">
        <rect width="305" height="88" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stage 1: Coarse cell fraction initialization</text>

        <text x="10" y="28" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Estimates baseline fractions θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan><tspan baseline-shift="super" font-size="75%">(0)</tspan> using static single-cell profiles φ<tspan baseline-shift="sub" font-size="75%">t, ·</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan></text>
        <text x="10" y="39" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Identifies tumor-specific and highly variable non-linear genes</text>
        <text x="10" y="50" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Excludes genes with severe baseline deviation from expected prior</text>
        
        <rect x="8" y="58" width="289" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="14" y="72" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Yields robust starting fractions θ<tspan baseline-shift="sub" font-size="75%">0</tspan> &amp; filtered signature subset G<tspan baseline-shift="sub" font-size="75%">sub</tspan></text>
      </g>

      <!-- Block 2: Stage 2 Joint Gibbs Sampling Loop -->
      <g id="col3-box-stage2" transform="translate(10, 132)">
        <rect width="305" height="210" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Stage 2: Joint iterative Gibbs sampling (Z &amp; φ):</text>

        <!-- Step A: Sample Z -->
        <g transform="translate(8, 24)">
          <rect width="289" height="66" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{OKABE_VERMILION}">Step A: Sample Latent Read Counts Z</text>
          <text x="135" y="30" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" font-family="{FONT_SERIF_MATH}">
            Z<tspan baseline-shift="sub" font-size="75%">·, g, n</tspan> | Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>, θ, φ ~ Multinomial( Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>,  π<tspan baseline-shift="sub" font-size="75%">g, n</tspan> )
          </text>
          <text x="275" y="30" font-size="8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(3)</text>
          <text x="8" y="44" font-size="6.8" fill="{COLOR_TEXT_MUTED}">where allocation vector π<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ∝ θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan></text>
          <text x="8" y="56" font-size="6.5" font-weight="600" fill="{OKABE_VERMILION}">Allocates bulk reads into lineages based on fraction &amp; profile</text>
        </g>

        <!-- Step B: Sample Phi -->
        <g transform="translate(8, 96)">
          <rect width="289" height="66" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{OKABE_ORANGE}">Step B: Update Cell-Type Profile φ</text>
          <text x="135" y="30" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" font-family="{FONT_SERIF_MATH}">
            φ<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> | Z, α ~ Dirichlet( α<tspan baseline-shift="sub" font-size="75%">t, ·</tspan> + Z<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> )
          </text>
          <text x="275" y="30" font-size="8" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(4)</text>
          <text x="8" y="44" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Exact closed-form update via Dirichlet-Multinomial conjugacy</text>
          <text x="8" y="56" font-size="6.5" font-weight="600" fill="{OKABE_ORANGE}">Adapts profile to patient-specific transcriptional activation</text>
        </g>

        <!-- Loop Indicator Banner -->
        <rect x="8" y="170" width="289" height="30" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="144" y="184" font-size="7" font-weight="700" fill="{OKABE_REDDISH_PURPLE}" text-anchor="middle">
          Iterate Step A ⇄ Step B for M Gibbs Cycles
        </text>
        <text x="144" y="194" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Guarantees convergence to stationary posterior distribution</text>
      </g>

      <!-- Block 3: MCMC Convergence & Formal Axis Plot -->
      <g id="col3-box-convergence" transform="translate(10, 350)">
        <rect width="305" height="156" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Convergence &amp; posterior marginalization:</text>

        <!-- Formal Coordinate Axes Trace Plot with Tick Marks -->
        <g transform="translate(12, 24)">
          <rect width="115" height="60" fill="{COLOR_CANVAS_BG}" />
          <!-- Axes -->
          <line x1="16" y1="8" x2="16" y2="48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="16" y1="48" x2="110" y2="48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <!-- Ticks Y -->
          <line x1="13" y1="12" x2="16" y2="12" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="13" y1="30" x2="16" y2="30" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="13" y1="48" x2="16" y2="48" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <!-- Ticks X -->
          <line x1="45" y1="48" x2="45" y2="51" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="75" y1="48" x2="75" y2="51" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="105" y1="48" x2="105" y2="51" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />

          <!-- Axis Labels -->
          <text x="6" y="32" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90 6,32)">θ<tspan baseline-shift="sub" font-size="75%">tn</tspan></text>
          <text x="65" y="58" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">MCMC iteration (m)</text>

          <!-- Gibbs Trace Line -->
          <path d="M 18 38 L 26 18 L 34 32 L 42 22 L 50 26 L 58 24 L 68 25 L 78 24 L 88 25 L 98 24 L 108 24" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="1.2" fill="none" />
          <line x1="48" y1="8" x2="48" y2="48" stroke="{OKABE_VERMILION}" stroke-width="0.75" stroke-dasharray="1.5,1.5" />
          <text x="32" y="14" font-size="5" fill="{COLOR_TEXT_MUTED}">Burn-in</text>
          <text x="75" y="14" font-size="5" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">Stationary</text>
        </g>

        <!-- Marginalization Equations -->
        <g transform="translate(136, 24)">
          <rect width="158" height="60" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Posterior expectation:</text>
          <text x="6" y="28" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" font-family="{FONT_SERIF_MATH}">
            θ̂<tspan baseline-shift="sub" font-size="75%">t, n</tspan> = (1/M) ∑ θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>
          </text>
          <text x="6" y="38" font-size="6" fill="{COLOR_TEXT_MUTED}">→ Posterior cell fraction θ*</text>
          <text x="6" y="48" font-size="7.5" font-weight="700" fill="{OKABE_VERMILION}" font-family="{FONT_SERIF_MATH}">
            Ẑ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> = (1/M) ∑ Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>
          </text>
          <text x="6" y="57" font-size="6" fill="{OKABE_VERMILION}">→ Deconvolved expression Z*</text>
        </g>

        <g transform="translate(10, 94)">
          <text x="0" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Evaluated across M stationary draws after discarding burn-in iterations</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Accurately infers patient-specific somatic malignant CNA programs</text>
          <text x="0" y="32" font-size="6.8" font-weight="700" fill="{OKABE_REDDISH_PURPLE}">• Eliminates requirement for patient-matched tumor scRNA-seq references</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: Column 4 — Deconvolved Dual Outputs
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 05: COLUMN 4 — DECONVOLVED OUTPUTS -->
  <g inkscape:groupmode="layer" id="layer-05-stage4" inkscape:label="04_Stage_4_Dual_Outputs">
    <g id="grp-col-4" transform="translate(1066, 82)">
      <!-- Column 4 Frame (width=298, height=518, rx=0) -->
      <rect width="298" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Column Header Banner -->
      <line x1="0" y1="28" x2="298" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="298" y2="0" stroke="{OKABE_BLUISH_GREEN}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        4. DECONVOLVED DUAL OUTPUTS
      </text>

      <!-- Block 1: Output 1 — Fractions Theta -->
      <g id="col4-box-out1" transform="translate(10, 36)">
        <rect width="278" height="116" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Output 1: Cell Fractions θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> (Simplex Δ<tspan baseline-shift="super" font-size="75%">T-1</tspan>)</text>

        <!-- Minimalist Stacked Bar (Nature style) -->
        <g transform="translate(8, 24)">
          <rect width="262" height="18" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <rect x="0" y="0" width="80" height="18" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="80" y="0" width="60" height="18" fill="{OKABE_ORANGE}" opacity="0.85" />
          <rect x="140" y="0" width="44" height="18" fill="{OKABE_BLUISH_GREEN}" opacity="0.85" />
          <rect x="184" y="0" width="78" height="18" fill="{OKABE_VERMILION}" opacity="0.85" />
          <text x="40" y="12" font-size="6.5" font-weight="700" fill="#FFF" text-anchor="middle">CD8 T (31%)</text>
          <text x="110" y="12" font-size="6.5" font-weight="700" fill="#FFF" text-anchor="middle">Mye (23%)</text>
          <text x="162" y="12" font-size="6.5" font-weight="700" fill="#FFF" text-anchor="middle">B (17%)</text>
          <text x="223" y="12" font-size="6.5" font-weight="700" fill="#FFF" text-anchor="middle">Tumor (29%)</text>
        </g>

        <g transform="translate(8, 54)">
          <text x="0" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Patient-level lineage fractions: ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> = 1</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Sub-state resolution: T<tspan baseline-shift="sub" font-size="75%">reg</tspan>, Exhausted CD8, M1/M2 macrophages</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Correlates directly with RECIST immunotherapy response</text>
          <text x="0" y="44" font-size="6.8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Standard cellular composition endpoint</text>
        </g>
      </g>

      <!-- Block 2: Output 2 — Expression Tensor Z -->
      <g id="col4-box-out2" transform="translate(10, 160)">
        <rect width="278" height="136" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Output 2: Expression Tensor Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> (T × G × N)</text>

        <!-- Clean Isometric Slices -->
        <g transform="translate(8, 24)">
          <rect width="262" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect x="6" y="4" width="44" height="18" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_BLUE}" stroke-width="0.8" />
          <rect x="11" y="8" width="44" height="18" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_ORANGE}" stroke-width="0.8" />
          <rect x="16" y="12" width="44" height="18" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_BLUISH_GREEN}" stroke-width="0.8" />
          <text x="38" y="24" font-size="6.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">CD8+ T</text>

          <g transform="translate(68, 4)">
            <text x="0" y="11" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3D Count Tensor</text>
            <text x="0" y="22" font-size="6.5" fill="{COLOR_TEXT_MUTED}">T lineages × G genes × N tumors</text>
          </g>
          <rect x="194" y="8" width="60" height="16" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_VERMILION}" stroke-width="0.8" />
          <text x="224" y="19" font-size="6.5" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Exact Sum</text>
        </g>

        <g transform="translate(8, 70)">
          <text x="0" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Lineage-specific read count matrix per patient</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Removes confounding from cell type abundance changes</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Distinguishes true transcriptional activation from cell proliferation</text>
          <text x="0" y="44" font-size="6.8" font-weight="700" fill="{OKABE_VERMILION}">Enables lineage-specific differential expression (DESeq2)</text>
        </g>
      </g>

      <!-- Block 3: Translational Capabilities -->
      <g id="col4-box-translational" transform="translate(10, 304)">
        <rect width="278" height="202" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Downstream translational capabilities:</text>

        <!-- Cap 1 -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Lineage-Specific DGE</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Run DESeq2 / edgeR directly on Ẑ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan></text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Isolates true gene induction per immune cell type</text>
          <text x="0" y="42" font-size="6.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Solves cellularity confounding</text>
          <line x1="0" y1="48" x2="258" y2="48" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        </g>

        <!-- Cap 2 -->
        <g transform="translate(10, 82)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Tumor Heterogeneity &amp; CNA</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Deconvolves patient malignant profiles</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Infers patient-specific somatic driver alterations</text>
          <text x="0" y="42" font-size="6.5" font-weight="700" fill="{OKABE_VERMILION}">No patient-matched scRNA required</text>
          <line x1="0" y1="48" x2="258" y2="48" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        </g>

        <!-- Cap 3 -->
        <g transform="translate(10, 140)">
          <text x="0" y="10" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. High-Throughput Trials</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Scales across 1,000s of bulk trial samples</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Connects bulk survival (OS/PFS) to single-cell</text>
          <text x="0" y="42" font-size="6.5" font-weight="700" fill="{OKABE_REDDISH_PURPLE}">Scalable trial biomarker discovery</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 6: Gutter Connecting Arrows
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 06: GUTTER CONNECTING ARROWS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Pipeline_Connectors">
    <!-- Connector 1 -> 2 (x=341 to x=366) -->
    <path d="M 342 240 L 364 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 342 410 L 364 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 2 -> 3 (x=691 to x=716) -->
    <path d="M 692 240 L 714 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 692 410 L 714 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 3 -> 4 (x=1041 to x=1066) -->
    <path d="M 1042 240 L 1064 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1042 410 L 1064 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the Nature Methods BayesPrism SVG."""
    output_dir = Path("article/figures/deconvolution")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "how_bayesprism_works.svg"
    svg_content = render_nature_methods_bayesprism()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Nature Methods BayesPrism SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
