#!/usr/bin/env python3
"""
Generate Publication Visual Overview: BayesPrism Deconvolution Framework.

Adheres strictly to Nature Methods / Nature editorial standards:
- Double-column landscape (1400 px width, 630 px height)
- Minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional marks
- High data-ink ratio: concise, professional scientific nomenclature with zero redundant text
- Clear graphical storytelling (heterogeneous tissue, reference signatures, read partitioning, MCMC Gibbs cycle, dual outputs)
- Native Inkscape layers and semantic groups with coordinate isolation
"""

from pathlib import Path
from nature_style_config import (
    WIDTH, HEIGHT,
    COLOR_CANVAS_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS
)


def render_visual_bayesprism() -> str:
    """Render the simplified, professional visual BayesPrism workflow SVG."""
    svg_parts: list[str] = []

    # 1. XML Header and Root SVG
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-bayesprism-visual-overview"
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
    <marker id="arr-orange" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_ORANGE}" />
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
        BayesPrism: Probabilistic Deconvolution of Cell Composition and Lineage Transcriptomes
      </text>
      <!-- Subtitle -->
      <text id="txt-subtitle" x="0" y="34" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        Partitioning composite bulk RNA-seq read counts into cell fractions (θ) and lineage-specific transcript counts (Z)
      </text>

      <!-- Academic Reference Tag -->
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
    # LAYER 2: Column 1 — Input Modalities
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: COLUMN 1 — INPUT MODALITIES -->
  <g inkscape:groupmode="layer" id="layer-02-col1" inkscape:label="02_Column_1_Inputs">
    <g id="grp-col-1" transform="translate(36, 82)">
      <!-- Frame -->
      <rect width="305" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Banner -->
      <line x1="0" y1="28" x2="305" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="305" y2="0" stroke="{OKABE_BLUE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        1. INPUT PROFILING MODALITIES
      </text>

      <!-- Box 1A: Patient Bulk RNA-seq (Spacious & Clean, Objective Box Removed) -->
      <g id="box-bulk-input" transform="translate(10, 36)">
        <rect width="285" height="134" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">A. Tumor Tissue RNA-seq (Bulk Mixture)</text>
        <text x="10" y="27" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Intact clinical biopsy; composite mixture of cellular transcriptomes</text>

        <!-- Compact Schematic: Tissue -> RNA-seq -> Read Vector Y -->
        <g transform="translate(10, 36)">
          <!-- Left: Microenvironment Specimen Core -->
          <g transform="translate(0, 0)">
            <rect width="112" height="84" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="56" y="13" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Tumor Biopsy Core</text>

            <!-- Histology schematic dots -->
            <g transform="translate(14, 18)">
              <!-- Malignant Cluster Boundary -->
              <rect x="28" y="6" width="46" height="34" rx="3" fill="{OKABE_VERMILION}" opacity="0.12" stroke="{OKABE_VERMILION}" stroke-width="0.5" stroke-dasharray="2,2" />
              <circle cx="38" cy="16" r="4" fill="{OKABE_VERMILION}" />
              <circle cx="52" cy="22" r="4.5" fill="{OKABE_VERMILION}" />
              <circle cx="62" cy="30" r="3.5" fill="{OKABE_VERMILION}" />
              
              <!-- Infiltrating Immune Cells -->
              <circle cx="8" cy="10" r="3" fill="{OKABE_BLUE}" /><!-- CD8 -->
              <circle cx="14" cy="28" r="3" fill="{OKABE_BLUE}" />
              <circle cx="22" cy="42" r="3" fill="{OKABE_BLUE}" />
              <circle cx="10" cy="40" r="3" fill="{OKABE_ORANGE}" /><!-- Mye -->
              <circle cx="70" cy="10" r="3" fill="{OKABE_ORANGE}" />
              <circle cx="78" cy="36" r="3" fill="{OKABE_BLUISH_GREEN}" /><!-- B -->
              <circle cx="74" cy="44" r="3" fill="{OKABE_BLUISH_GREEN}" />
            </g>

            <!-- Lineage Micro-Legend -->
            <text x="56" y="74" font-size="5.6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">CD8+ T • Mye • B • Malignant</text>
          </g>

          <!-- Middle: Centered Arrow & RNA-seq Label (Fixed Alignment in Gap x=112..136, y=42) -->
          <g transform="translate(112, 0)">
            <text x="12" y="35" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">RNA-seq</text>
            <line x1="3" y1="42" x2="21" y2="42" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" marker-end="url(#arr-slate)" />
          </g>

          <!-- Right: Composite Counts Vector Y -->
          <g transform="translate(136, 0)">
            <rect width="129" height="84" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="64" y="13" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Composite Counts (Y)</text>

            <!-- Slim Bar Profile (Genes g1 ... g8) -->
            <g transform="translate(14, 18)">
              <!-- Hairline axes -->
              <line x1="0" y1="0" x2="0" y2="44" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
              <line x1="0" y1="44" x2="98" y2="44" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
              
              <!-- Subtle Gene Bars -->
              <rect x="5" y="16" width="7" height="28" fill="{OKABE_BLUE}" opacity="0.85" />
              <rect x="17" y="6" width="7" height="38" fill="{OKABE_VERMILION}" opacity="0.85" />
              <rect x="29" y="22" width="7" height="22" fill="{OKABE_ORANGE}" opacity="0.85" />
              <rect x="41" y="8" width="7" height="36" fill="{OKABE_BLUE}" opacity="0.85" />
              <rect x="53" y="26" width="7" height="18" fill="{OKABE_BLUISH_GREEN}" opacity="0.85" />
              <rect x="65" y="4" width="7" height="40" fill="{OKABE_VERMILION}" opacity="0.85" />
              <rect x="77" y="18" width="7" height="26" fill="{OKABE_ORANGE}" opacity="0.85" />
              <rect x="89" y="12" width="7" height="32" fill="{OKABE_BLUE}" opacity="0.85" />
              
              <text x="49" y="54" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (g = 1 … G) →</text>
            </g>
          </g>
        </g>
      </g>

      <!-- Box 1B: Single-Cell Reference Signatures (Footer Removed, Generous Height) -->
      <g id="box-sc-input" transform="translate(10, 180)">
        <rect width="285" height="326" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="18" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">B. Single-Cell Reference Signatures</text>
        <text x="10" y="30" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Baseline expression profiles across curated lineages</text>

        <!-- Pure Cell Catalog Grid -->
        <g transform="translate(10, 42)">
          <!-- Lineage 1: CD8 T -->
          <g transform="translate(0, 0)">
            <rect width="265" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="20" cy="29" r="8" fill="{OKABE_BLUE}" /><circle cx="20" cy="29" r="3.5" fill="#FFF" />
            <text x="38" y="22" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">CD8+ T Lymphocytes</text>
            <text x="38" y="36" font-size="6.8" fill="{COLOR_TEXT_MUTED}">CD8A, CD8B, GZMB, IFNG, PRF1</text>
            <text x="38" y="48" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Cytotoxic effector &amp; exhaustion markers</text>
          </g>

          <!-- Lineage 2: Myeloid -->
          <g transform="translate(0, 68)">
            <rect width="265" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="20" cy="29" r="8" fill="{OKABE_ORANGE}" /><circle cx="20" cy="29" r="3.5" fill="#FFF" />
            <text x="38" y="22" font-size="7.5" font-weight="700" fill="{OKABE_ORANGE}">Myeloid Compartment</text>
            <text x="38" y="36" font-size="6.8" fill="{COLOR_TEXT_MUTED}">CD14, CD68, CD163, MRC1, APOE</text>
            <text x="38" y="48" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">M1 inflammatory vs. M2 immunosuppressive</text>
          </g>

          <!-- Lineage 3: B Cells -->
          <g transform="translate(0, 136)">
            <rect width="265" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="20" cy="29" r="8" fill="{OKABE_BLUISH_GREEN}" /><circle cx="20" cy="29" r="3.5" fill="#FFF" />
            <text x="38" y="22" font-size="7.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}">B &amp; Plasma Cells</text>
            <text x="38" y="36" font-size="6.8" fill="{COLOR_TEXT_MUTED}">MS4A1 (CD20), CD19, MZB1, JCHAIN</text>
            <text x="38" y="48" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Humoral response &amp; antibody secretion</text>
          </g>

          <!-- Lineage 4: Malignant -->
          <g transform="translate(0, 204)">
            <rect width="265" height="58" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="20" cy="29" r="8" fill="{OKABE_VERMILION}" /><circle cx="20" cy="29" r="3.5" fill="#FFF" />
            <text x="38" y="22" font-size="7.5" font-weight="700" fill="{OKABE_VERMILION}">Malignant Compartment</text>
            <text x="38" y="36" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Patient-independent baseline prior (α)</text>
            <text x="38" y="48" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">Heterogeneous somatic expression states</text>
          </g>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: Column 2 — Probabilistic Read Allocation
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 03: COLUMN 2 — READ ALLOCATION PRINCIPLE -->
  <g inkscape:groupmode="layer" id="layer-03-col2" inkscape:label="03_Column_2_Read_Allocation">
    <g id="grp-col-2" transform="translate(366, 82)">
      <!-- Frame -->
      <rect width="325" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Banner -->
      <line x1="0" y1="28" x2="325" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="325" y2="0" stroke="{OKABE_ORANGE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        2. PROBABILISTIC READ ALLOCATION
      </text>

      <!-- Section: BayesPrism Joint Inversion (Conventional Deconvolution Box Removed) -->
      <g id="box-concept-compare" transform="translate(10, 36)">
        <rect width="305" height="74" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Methodological Paradigm</text>

        <!-- BayesPrism Dynamic Inversion -->
        <g transform="translate(10, 24)">
          <rect width="285" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="0.75" />
          <text x="8" y="13" font-size="7.2" font-weight="700" fill="{OKABE_BLUE}">BayesPrism (Joint Generative Allocation)</text>
          <text x="8" y="25" font-size="6.8" fill="{COLOR_TEXT_PRIMARY}">• Partitions observed reads Y into lineage-specific counts Z</text>
          <text x="8" y="35" font-size="6.8" fill="{OKABE_BLUISH_GREEN}">• Jointly infers cellular composition θ and state-specific expression φ</text>
        </g>
      </g>

      <!-- Section: Optical Dispersion Metaphor (Equilateral Triangle, Bottom Text Removed) -->
      <g id="box-prism-graphic" transform="translate(10, 118)">
        <rect width="305" height="388" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Deconvolution Mechanism: Partitioning Counts</text>

        <!-- Optical Dispersion Graphic -->
        <g transform="translate(10, 30)">
          <!-- Mixed Read Stream (Left) -->
          <g transform="translate(0, 78)">
            <rect width="70" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
            <text x="35" y="12" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Bulk Counts</text>
            <text x="35" y="22" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Y ∈ ℕ<tspan baseline-shift="super" font-size="75%">G × N</tspan></text>
            <!-- Ray into Prism Face (from x=70 to x=115 at y=14, which touches left face of equilateral prism) -->
            <line x1="70" y1="14" x2="114" y2="14" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" marker-end="url(#arr-slate)" />
          </g>

          <!-- Central Optical Prism (Equilateral Geometry: Base=100, Height=86, Ratio Balanced) -->
          <g transform="translate(88, 48)">
            <polygon points="50,0 100,86 0,86" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="1.5" />
            <text x="50" y="46" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Bayes</text>
            <text x="50" y="58" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">Prism</text>
            <text x="50" y="71" font-size="6.2" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Inversion</text>
          </g>

          <!-- Lineage Read Rays (Emanate from Right Face of Equilateral Prism into Output Boxes) -->
          <!-- CD8 T: ray from (148, 64) -> (205, 24) -->
          <path d="M 148 64 L 205 24" stroke="{OKABE_BLUE}" stroke-width="1.6" fill="none" marker-end="url(#arr-blue)" />
          <g transform="translate(210, 10)">
            <rect width="75" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="0.75" />
            <text x="6" y="12" font-size="7" font-weight="700" fill="{OKABE_BLUE}">Z<tspan baseline-shift="sub" font-size="75%">CD8</tspan> Counts</text>
            <text x="6" y="22" font-size="6" fill="{COLOR_TEXT_MUTED}">Cytotoxic T cells</text>
          </g>

          <!-- Myeloid: ray from (158, 80) -> (205, 60) -->
          <path d="M 158 80 L 205 60" stroke="{OKABE_ORANGE}" stroke-width="1.6" fill="none" marker-end="url(#arr-orange)" />
          <g transform="translate(210, 46)">
            <rect width="75" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_ORANGE}" stroke-width="0.75" />
            <text x="6" y="12" font-size="7" font-weight="700" fill="{OKABE_ORANGE}">Z<tspan baseline-shift="sub" font-size="75%">Mye</tspan> Counts</text>
            <text x="6" y="22" font-size="6" fill="{COLOR_TEXT_MUTED}">Myeloid lineage</text>
          </g>

          <!-- B Cells: ray from (168, 97) -> (205, 96) -->
          <path d="M 168 97 L 205 96" stroke="{OKABE_BLUISH_GREEN}" stroke-width="1.6" fill="none" marker-end="url(#arr-green)" />
          <g transform="translate(210, 82)">
            <rect width="75" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUISH_GREEN}" stroke-width="0.75" />
            <text x="6" y="12" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Z<tspan baseline-shift="sub" font-size="75%">B</tspan> Counts</text>
            <text x="6" y="22" font-size="6" fill="{COLOR_TEXT_MUTED}">B &amp; plasma cells</text>
          </g>

          <!-- Malignant: ray from (178, 114) -> (205, 132) -->
          <path d="M 178 114 L 205 132" stroke="{OKABE_VERMILION}" stroke-width="1.6" fill="none" marker-end="url(#arr-slate)" />
          <g transform="translate(210, 118)">
            <rect width="75" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_VERMILION}" stroke-width="0.75" />
            <text x="6" y="12" font-size="7" font-weight="700" fill="{OKABE_VERMILION}">Z<tspan baseline-shift="sub" font-size="75%">Mal</tspan> Counts</text>
            <text x="6" y="22" font-size="6" fill="{COLOR_TEXT_MUTED}">Malignant cells</text>
          </g>
        </g>

        <!-- Conservation Law Card (Clean, Prominent Mathematical Invariance) -->
        <g transform="translate(10, 210)">
          <rect width="285" height="74" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="10" y="16" font-size="7.5" font-weight="700" fill="{OKABE_BLACK}">Read Conservation Constraint:</text>
          <text x="10" y="32" font-size="7.5" font-weight="600" fill="{OKABE_VERMILION}">∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>,   Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ∈ ℕ<tspan baseline-shift="sub" font-size="75%">0</tspan></text>
          <text x="10" y="48" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">Strict mathematical invariance preserves raw read count scale.</text>
          <text x="10" y="60" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Avoids artifacts induced by non-linear library normalization.</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: Column 3 — Joint MCMC Gibbs Inference
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: COLUMN 3 — MCMC GIBBS INFERENCE -->
  <g inkscape:groupmode="layer" id="layer-04-col3" inkscape:label="03_Column_3_Iterative_Loop">
    <g id="grp-col-3" transform="translate(716, 82)">
      <!-- Frame -->
      <rect width="325" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Banner -->
      <line x1="0" y1="28" x2="325" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="325" y2="0" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        3. JOINT MCMC GIBBS INFERENCE
      </text>

      <!-- Latent Parameter Coupling -->
      <g id="box-coupling" transform="translate(10, 36)">
        <rect width="305" height="54" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Coupled Parameter Inversion:</text>
        <text x="10" y="27" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Cell fractions θ and lineage profiles φ are mutually dependent.</text>
        <text x="10" y="40" font-size="6.8" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">Resolved via cyclical two-stage Gibbs sampling.</text>
      </g>

      <!-- 2-Step Gibbs Sampling Engine -->
      <g id="box-cycle-diagram" transform="translate(10, 98)">
        <rect width="305" height="236" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cyclical Posterior Sampling Engine</text>

        <!-- Step 1: Estimate Fractions -->
        <g transform="translate(16, 28)">
          <rect width="215" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="1" />
          <circle cx="16" cy="18" r="8" fill="{OKABE_BLUE}" /><text x="16" y="21.5" font-size="7.5" font-weight="700" fill="#FFF" text-anchor="middle">1</text>
          <text x="32" y="16" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">Step 1: Sample Cell Proportions (θ)</text>
          <text x="32" y="28" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">Conditioned on current lineage expression φ</text>
          <text x="32" y="39" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Estimates lineage abundance vector θ ∈ Δ<tspan baseline-shift="super" font-size="75%">T-1</tspan></text>
        </g>

        <!-- Transition Down Arrow -->
        <g transform="translate(123, 82)">
          <line x1="0" y1="0" x2="0" y2="18" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" marker-end="url(#arr-slate)" />
        </g>

        <!-- Step 2: Update Expression -->
        <g transform="translate(16, 104)">
          <rect width="215" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_ORANGE}" stroke-width="1" />
          <circle cx="16" cy="18" r="8" fill="{OKABE_ORANGE}" /><text x="16" y="21.5" font-size="7.5" font-weight="700" fill="#FFF" text-anchor="middle">2</text>
          <text x="32" y="16" font-size="7.5" font-weight="700" fill="{OKABE_ORANGE}">Step 2: Sample Expression (Z &amp; φ)</text>
          <text x="32" y="28" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">Conditioned on estimated fractions θ and priors α</text>
          <text x="32" y="39" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Adapts to patient mutations &amp; activation states</text>
        </g>

        <!-- Clean Loopback Arc (Clear of border and cards) -->
        <path d="M 233 130 C 265 130, 265 54, 233 54" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="1.3" fill="none" marker-end="url(#arr-purple)" />
        
        <!-- Clean Horizontal Badge in Right Gutter -->
        <g transform="translate(244, 76)">
          <rect x="0" y="0" width="50" height="32" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="0.75" />
          <text x="25" y="13" font-size="6.2" font-weight="700" fill="{OKABE_REDDISH_PURPLE}" text-anchor="middle">Gibbs</text>
          <text x="25" y="23" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">m = 1…M</text>
        </g>

        <!-- Convergence Property -->
        <g transform="translate(16, 168)">
          <rect width="273" height="48" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Ergodic Posterior Convergence:</text>
          <text x="8" y="27" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Markov chain reaches unique stationary joint posterior</text>
          <text x="8" y="38" font-size="6.8" fill="{COLOR_TEXT_MUTED}">distribution p(θ, Z, φ | Y, α).</text>
        </g>
      </g>

      <!-- MCMC Trace Plot (Footer Text Removed, Full Posterior Label Removed) -->
      <g id="box-convergence-trace" transform="translate(10, 342)">
        <rect width="305" height="164" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Posterior Sampling Convergence</text>

        <g transform="translate(20, 36)">
          <!-- Coordinate Axes -->
          <line x1="20" y1="10" x2="20" y2="70" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          <line x1="20" y1="70" x2="160" y2="70" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="0.75" />
          
          <text x="8" y="44" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90 8,44)">Parameter</text>
          <text x="90" y="82" font-size="6.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">MCMC Iteration (m) →</text>

          <!-- Trace -->
          <path d="M 22 55 L 34 22 L 46 48 L 58 30 L 70 38 L 82 35 L 96 36 L 112 35 L 128 36 L 144 35 L 158 35"
                stroke="{OKABE_REDDISH_PURPLE}" stroke-width="1.3" fill="none" />
          
          <line x1="75" y1="10" x2="75" y2="70" stroke="{OKABE_VERMILION}" stroke-width="0.75" stroke-dasharray="2,2" />
          <text x="46" y="18" font-size="5.8" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Burn-in</text>
          <text x="115" y="18" font-size="5.8" font-weight="700" fill="{OKABE_REDDISH_PURPLE}" text-anchor="middle">Equilibrium</text>
        </g>

        <!-- Expectation Card (Clean Mathematical Point Estimates) -->
        <g transform="translate(185, 36)">
          <rect width="110" height="70" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="16" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Expectation:</text>
          <text x="8" y="34" font-size="6.8" fill="{COLOR_TEXT_MUTED}">θ̂ = (1/M) ∑ θ<tspan baseline-shift="super" font-size="75%">(m)</tspan></text>
          <text x="8" y="52" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Ẑ = (1/M) ∑ Z<tspan baseline-shift="super" font-size="75%">(m)</tspan></text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: Column 4 — Dual Deconvolved Outputs
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 05: COLUMN 4 — DUAL OUTPUTS (UNCLUTTERED & SPACIOUS) -->
  <g inkscape:groupmode="layer" id="layer-05-col4" inkscape:label="04_Column_4_Dual_Outputs">
    <g id="grp-col-4" transform="translate(1066, 82)">
      <!-- Frame -->
      <rect width="298" height="518" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      
      <!-- Banner -->
      <line x1="0" y1="28" x2="298" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <line x1="0" y1="0" x2="298" y2="0" stroke="{OKABE_BLUISH_GREEN}" stroke-width="2.5" />
      <text x="12" y="19" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        4. DUAL DECONVOLVED OUTPUTS
      </text>

      <!-- Output 1: Cellular Composition (θ) (Clinical Endpoint Box Removed) -->
      <g id="box-out-proportions" transform="translate(10, 36)">
        <rect width="278" height="174" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Output 1: Cellular Composition (θ ∈ Δ<tspan baseline-shift="super" font-size="75%">T-1</tspan>)</text>
        <text x="10" y="27" font-size="6.8" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Quantification of patient-level immune infiltration</text>

        <!-- Composition Stacked Bar -->
        <g transform="translate(8, 38)">
          <rect width="262" height="26" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <rect x="0" y="0" width="82" height="26" fill="{OKABE_BLUE}" opacity="0.9" />
          <rect x="82" y="0" width="60" height="26" fill="{OKABE_ORANGE}" opacity="0.9" />
          <rect x="142" y="0" width="44" height="26" fill="{OKABE_BLUISH_GREEN}" opacity="0.9" />
          <rect x="186" y="0" width="76" height="26" fill="{OKABE_VERMILION}" opacity="0.9" />
          
          <text x="41" y="16.5" font-size="6.8" font-weight="700" fill="#FFF" text-anchor="middle">CD8 T (31%)</text>
          <text x="112" y="16.5" font-size="6.8" font-weight="700" fill="#FFF" text-anchor="middle">Mye (23%)</text>
          <text x="164" y="16.5" font-size="6.8" font-weight="700" fill="#FFF" text-anchor="middle">B (17%)</text>
          <text x="224" y="16.5" font-size="6.8" font-weight="700" fill="#FFF" text-anchor="middle">Tumor (29%)</text>
        </g>

        <!-- Lineage Sub-State Specificity Details -->
        <g transform="translate(8, 76)">
          <rect width="262" height="82" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="15" font-size="7.2" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Sub-State Abundance Resolution:</text>
          <text x="8" y="30" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">• CD8+ T cells: cytotoxic effector vs. exhausted vs. memory</text>
          <text x="8" y="46" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">• Myeloid: pro-inflammatory M1 vs. immunosuppressive M2</text>
          <text x="8" y="62" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">• Lymphocytes: regulatory T (Treg), naive B, and plasma cells</text>
        </g>
      </g>

      <!-- Output 2: Lineage-Specific Transcriptomes (Z) (Confounder Callout Removed, Generous Spacing) -->
      <g id="box-out-expression" transform="translate(10, 220)">
        <rect width="278" height="286" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Output 2: Lineage-Specific Transcriptomes (Z)</text>
        <text x="10" y="28" font-size="6.8" font-weight="600" fill="{OKABE_VERMILION}">Isolated gene expression matrices (T × G × N)</text>

        <!-- Lineage Expression Cards -->
        <g transform="translate(8, 38)">
          <!-- CD8 T Transcriptome -->
          <g transform="translate(0, 0)">
            <rect width="262" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="0.75" />
            <text x="8" y="15" font-size="7.2" font-weight="700" fill="{OKABE_BLUE}">CD8+ T-Cell Transcriptome</text>
            <text x="8" y="28" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">Cytolytic effector programs &amp; exhaustion markers (GZMB, IFNG)</text>
            <text x="8" y="42" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Resolves activation shifts independent of cell frequency</text>
          </g>

          <!-- Myeloid Transcriptome -->
          <g transform="translate(0, 60)">
            <rect width="262" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_ORANGE}" stroke-width="0.75" />
            <text x="8" y="15" font-size="7.2" font-weight="700" fill="{OKABE_ORANGE}">Myeloid Lineage Transcriptome</text>
            <text x="8" y="28" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">Immunosuppressive M2 vs. pro-inflammatory M1 polarization</text>
            <text x="8" y="42" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Phagocytic activity and antigen presentation capacity</text>
          </g>

          <!-- B & Plasma Transcriptome -->
          <g transform="translate(0, 120)">
            <rect width="262" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUISH_GREEN}" stroke-width="0.75" />
            <text x="8" y="15" font-size="7.2" font-weight="700" fill="{OKABE_BLUISH_GREEN}">B-Cell &amp; Plasma Transcriptome</text>
            <text x="8" y="28" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">Immunoglobulin repertoire and antibody-mediated immunity</text>
            <text x="8" y="42" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Tertiary lymphoid structure (TLS) signature evaluation</text>
          </g>

          <!-- Malignant Transcriptome -->
          <g transform="translate(0, 180)">
            <rect width="262" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_VERMILION}" stroke-width="0.75" />
            <text x="8" y="15" font-size="7.2" font-weight="700" fill="{OKABE_VERMILION}">Malignant Cell Transcriptome</text>
            <text x="8" y="28" font-size="6.6" fill="{COLOR_TEXT_SECONDARY}">Patient somatic driver alterations &amp; immune evasion pathways</text>
            <text x="8" y="42" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Captures patient-specific aberrant gene expression states</text>
          </g>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 6: Connecting Gutter Arrows
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 06: GUTTER CONNECTING ARROWS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Pipeline_Connectors">
    <!-- Col 1 -> Col 2 -->
    <path d="M 342 200 L 364 200" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 342 390 L 364 390" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />

    <!-- Col 2 -> Col 3 -->
    <path d="M 692 200 L 714 200" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 692 390 L 714 390" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />

    <!-- Col 3 -> Col 4 -->
    <path d="M 1042 200 L 1064 200" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1042 390 L 1064 390" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1.2" fill="none" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the simplified visual BayesPrism SVG."""
    output_dir = Path("article/figures/deconvolution")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "how_bayesprism_works_visual.svg"
    svg_content = render_visual_bayesprism()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Visual BayesPrism SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
