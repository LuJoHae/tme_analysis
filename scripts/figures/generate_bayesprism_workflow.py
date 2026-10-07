#!/usr/bin/env python3
"""
Generate Publication-Grade Vector Figure: How BayesPrism Works.

Illustrates the Bayesian statistical framework for joint deconvolution of
cell type fractions and cell-type-specific gene expression from bulk RNA-seq
using a single-cell reference (Chu et al., Nature Cancer 2022).
Structured with native Inkscape layers, semantic groups, and Nature editorial styling.
"""

from pathlib import Path
import html

# Canvas Dimensions (Double-column landscape, 16:10 aspect ratio ~180 mm journal width)
WIDTH = 1400
HEIGHT = 880

# Nature / Nature Cancer Editorial Palette
COLOR_CANVAS_BG = "#F8FAFC"
COLOR_CARD_BG = "#FFFFFF"
COLOR_CARD_BORDER = "#E2E8F0"

COLOR_TEXT_DARK = "#0F172A"       # Slate 900
COLOR_TEXT_MUTED = "#475569"      # Slate 600
COLOR_TEXT_LIGHT = "#94A3B8"      # Slate 400

COLOR_NAVY = "#1D3557"            # Deep Navy (Headers, structural frames)
COLOR_STEEL = "#457B9D"           # Steel Blue (Bulk RNA-seq)
COLOR_TEAL = "#2A9D8F"            # Deep Teal (Single-cell reference)
COLOR_EMERALD = "#059669"         # Emerald Green (Cell fractions, outputs)
COLOR_PURPLE = "#7C3AED"          # Violet / Purple (Bayesian MCMC inference)
COLOR_AMBER = "#D97706"           # Warm Amber (Priors & likelihood model)
COLOR_CRIMSON = "#DC2626"         # Crimson Red (Malignant cells & allocation)


def escape_txt(text: str) -> str:
    """Safely escape text for XML/SVG rendering."""
    return html.escape(text)


def render_svg() -> str:
    """Declaratively render the complete publication-grade SVG DOM."""
    svg_parts: list[str] = []

    # 1. XML Header and SVG root with Inkscape namespaces
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-bayesprism-mechanism"
     style="background-color: {COLOR_CANVAS_BG}; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;">

  <defs>
    <!-- Card Dropshadow -->
    <filter id="card-shadow" x="-3%" y="-2%" width="106%" height="106%" filterUnits="userSpaceOnUse">
      <feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="#0F172A" flood-opacity="0.04" />
    </filter>

    <!-- Arrow Markers -->
    <marker id="arrow-navy" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR_NAVY}" />
    </marker>
    <marker id="arrow-teal" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR_TEAL}" />
    </marker>
    <marker id="arrow-purple" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR_PURPLE}" />
    </marker>
    <marker id="arrow-emerald" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR_EMERALD}" />
    </marker>
    <marker id="arrow-amber" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR_AMBER}" />
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
    <rect id="card-panel-a" x="20" y="68" width="665" height="392" rx="10" ry="10"
          fill="{COLOR_CARD_BG}" stroke="{COLOR_CARD_BORDER}" stroke-width="1.2" filter="url(#card-shadow)" />

    <!-- Panel B Container Card (Top-Right) -->
    <rect id="card-panel-b" x="715" y="68" width="665" height="392" rx="10" ry="10"
          fill="{COLOR_CARD_BG}" stroke="{COLOR_CARD_BORDER}" stroke-width="1.2" filter="url(#card-shadow)" />

    <!-- Panel C Container Card (Bottom-Left) -->
    <rect id="card-panel-c" x="20" y="474" width="665" height="388" rx="10" ry="10"
          fill="{COLOR_CARD_BG}" stroke="{COLOR_CARD_BORDER}" stroke-width="1.2" filter="url(#card-shadow)" />

    <!-- Panel D Container Card (Bottom-Right) -->
    <rect id="card-panel-d" x="715" y="474" width="665" height="388" rx="10" ry="10"
          fill="{COLOR_CARD_BG}" stroke="{COLOR_CARD_BORDER}" stroke-width="1.2" filter="url(#card-shadow)" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 1: Header, Subtitle & Panel Badges
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 01: FIGURE HEADER & PANEL BADGES -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <!-- Main Title -->
    <text id="txt-title" x="25" y="32" font-size="16.5" font-weight="700" fill="{COLOR_NAVY}" letter-spacing="-0.2">
      BayesPrism Statistical Mechanism | Joint Bayesian Deconvolution of Cell Types &amp; Expression Profiles
    </text>
    <!-- Subtitle Summary -->
    <text id="txt-subtitle" x="25" y="52" font-size="11.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Probabilistic generative mixture modeling inverting bulk RNA-seq mixtures into cell fractions (θ) and cell-type-specific read counts (Z)
    </text>

    <!-- Top Method Citation Pill -->
    <g id="grp-header-stat-pill" transform="translate(1075, 16)">
      <rect x="0" y="0" width="305" height="34" rx="17" ry="17" fill="#F5F3FF" stroke="#DDD6FE" stroke-width="1" />
      <circle cx="18" cy="17" r="7" fill="{COLOR_PURPLE}" />
      <text x="18" y="21" font-size="8.5" font-weight="700" fill="#FFF" text-anchor="middle">BP</text>
      <text x="32" y="21.5" font-size="10" font-weight="700" fill="#5B21B6">Chu et al., Nature Cancer 2022</text>
    </g>

    <!-- Panel Badges [a], [b], [c], [d] -->
    <g id="badge-panel-a" transform="translate(34, 82)">
      <circle cx="11" cy="11" r="11" fill="{COLOR_NAVY}" />
      <text x="11" y="15" font-size="12" font-weight="700" fill="#FFFFFF" text-anchor="middle">a</text>
      <text x="30" y="15" font-size="13" font-weight="700" fill="{COLOR_NAVY}">Data Inputs &amp; Prior Construction (Single-Cell &amp; Bulk)</text>
    </g>

    <g id="badge-panel-b" transform="translate(729, 82)">
      <circle cx="11" cy="11" r="11" fill="{COLOR_NAVY}" />
      <text x="11" y="15" font-size="12" font-weight="700" fill="#FFFFFF" text-anchor="middle">b</text>
      <text x="30" y="15" font-size="13" font-weight="700" fill="{COLOR_NAVY}">Probabilistic Generative Mixture Model (Likelihood &amp; Prior)</text>
    </g>

    <g id="badge-panel-c" transform="translate(34, 488)">
      <circle cx="11" cy="11" r="11" fill="{COLOR_NAVY}" />
      <text x="11" y="15" font-size="12" font-weight="700" fill="#FFFFFF" text-anchor="middle">c</text>
      <text x="30" y="15" font-size="13" font-weight="700" fill="{COLOR_NAVY}">Two-Stage Gibbs Sampling &amp; Latent Count Allocation (MCMC)</text>
    </g>

    <g id="badge-panel-d" transform="translate(729, 488)">
      <circle cx="11" cy="11" r="11" fill="{COLOR_NAVY}" />
      <text x="11" y="15" font-size="12" font-weight="700" fill="#FFFFFF" text-anchor="middle">d</text>
      <text x="30" y="15" font-size="13" font-weight="700" fill="{COLOR_NAVY}">Deconvolved Dual Outputs: Fractions (θ) &amp; Expression (Z)</text>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: PANEL A — Data Inputs & Prior Construction
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 02: PANEL A — DATA INPUTS & PRIOR CONSTRUCTION -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Inputs_and_Priors">

    <!-- Card 1: Single-Cell RNA-seq Reference Input -->
    <g id="grp-input-scrna" transform="translate(34, 114)">
      <rect width="310" height="190" rx="8" ry="8" fill="#F0FDFA" stroke="#99F6E4" stroke-width="1.2" />
      <path d="M 0 8 C 0 3.6 3.6 0 8 0 L 302 0 C 306.4 0 310 3.6 310 8 L 310 28 L 0 28 Z" fill="#CCFBF1" />
      <circle cx="14" cy="14" r="5" fill="{COLOR_TEAL}" />
      <text x="26" y="18" font-size="10.5" font-weight="700" fill="#115E59">Input 1: Single-Cell RNA-seq Reference (X)</text>
      <rect x="250" y="5" width="52" height="18" rx="9" fill="{COLOR_TEAL}" />
      <text x="276" y="17.5" font-size="8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Reference</text>

      <!-- Matrix Diagram -->
      <g transform="translate(12, 38)">
        <!-- Cells x Genes Matrix -->
        <rect width="110" height="96" rx="4" fill="#FFFFFF" stroke="#99F6E4" stroke-width="1" />
        <!-- Matrix Grid lines -->
        <line x1="22" y1="0" x2="22" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="44" y1="0" x2="44" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="66" y1="0" x2="66" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="88" y1="0" x2="88" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="0" y1="32" x2="110" y2="32" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="0" y1="64" x2="110" y2="64" stroke="#E2E8F0" stroke-width="0.8" />
        <!-- Heatmap pills -->
        <rect x="2" y="2" width="18" height="28" fill="#3B82F6" opacity="0.7" rx="2" />
        <rect x="24" y="34" width="18" height="28" fill="#F59E0B" opacity="0.7" rx="2" />
        <rect x="68" y="66" width="18" height="28" fill="#10B981" opacity="0.7" rx="2" />

        <text x="55" y="112" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
        <text x="-48" y="-4" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Cells (C) →</text>

        <!-- Cell Types & State annotations -->
        <g transform="translate(125, 4)">
          <text x="0" y="12" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_DARK}">Two-Tier Annotation Hierarchy:</text>
          <text x="0" y="26" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Cell Types t ∈ {{1, ..., T}} (coarse)</text>
          <text x="0" y="38" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Cell States s ∈ {{1, ..., S}} (fine sub-clusters)</text>
          <rect x="0" y="46" width="158" height="42" rx="4" fill="#FFFFFF" stroke="#CCFBF1" stroke-width="0.8" />
          <text x="6" y="58" font-size="7.5" font-weight="700" fill="#115E59">Dirichlet Prior: α_tg</text>
          <text x="6" y="70" font-size="7" fill="{COLOR_TEXT_MUTED}">α_tg = pseudo_count + Σ x_cg</text>
          <text x="6" y="81" font-size="7" fill="{COLOR_TEXT_MUTED}">Mean profile: φ_tg^ref = α_tg / Σ α_tg'</text>
        </g>
      </g>
      <text x="12" y="174" font-size="7.5" fill="{COLOR_TEXT_MUTED}">X_cg = raw count of gene g in cell c of cell type t</text>
    </g>

    <!-- Card 2: Bulk RNA-seq Count Input -->
    <g id="grp-input-bulk" transform="translate(360, 114)">
      <rect width="310" height="190" rx="8" ry="8" fill="#EFF6FF" stroke="#BFDBFE" stroke-width="1.2" />
      <path d="M 0 8 C 0 3.6 3.6 0 8 0 L 302 0 C 306.4 0 310 3.6 310 8 L 310 28 L 0 28 Z" fill="#DBEAFE" />
      <circle cx="14" cy="14" r="5" fill="{COLOR_STEEL}" />
      <text x="26" y="18" font-size="10.5" font-weight="700" fill="{COLOR_NAVY}">Input 2: Bulk RNA-seq Raw Counts (Y)</text>
      <rect x="250" y="5" width="52" height="18" rx="9" fill="{COLOR_STEEL}" />
      <text x="276" y="17.5" font-size="8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Bulk Mixture</text>

      <g transform="translate(12, 38)">
        <!-- Samples x Genes Matrix -->
        <rect width="110" height="96" rx="4" fill="#FFFFFF" stroke="#BFDBFE" stroke-width="1" />
        <line x1="28" y1="0" x2="28" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="56" y1="0" x2="56" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="84" y1="0" x2="84" y2="96" stroke="#E2E8F0" stroke-width="0.8" />
        <line x1="0" y1="48" x2="110" y2="48" stroke="#E2E8F0" stroke-width="0.8" />
        <!-- Bulk read bars -->
        <rect x="4" y="6" width="20" height="36" fill="{COLOR_STEEL}" opacity="0.6" rx="2" />
        <rect x="32" y="16" width="20" height="26" fill="{COLOR_STEEL}" opacity="0.8" rx="2" />
        <rect x="60" y="8" width="20" height="34" fill="{COLOR_STEEL}" opacity="0.4" rx="2" />

        <text x="55" y="112" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
        <text x="-48" y="-4" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Samples (N) →</text>

        <g transform="translate(125, 4)">
          <text x="0" y="12" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_DARK}">Key Input Properties:</text>
          <text x="0" y="26" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Must be integer raw read counts Y_gn</text>
          <text x="0" y="38" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Total reads in sample n: N_n = Σ Y_gn</text>
          
          <rect x="0" y="46" width="158" height="42" rx="4" fill="#FFFFFF" stroke="#BFDBFE" stroke-width="0.8" />
          <text x="6" y="58" font-size="7.5" font-weight="700" fill="{COLOR_NAVY}">Malignant Cell Handling</text>
          <text x="6" y="70" font-size="7" fill="{COLOR_TEXT_MUTED}">Allows tumor-specific drift;</text>
          <text x="6" y="81" font-size="7" fill="{COLOR_TEXT_MUTED}">TME reference treated as static base.</text>
        </g>
      </g>
      <text x="12" y="174" font-size="7.5" fill="{COLOR_TEXT_MUTED}">Y_gn = observed bulk read counts across patient tumors</text>
    </g>

    <!-- Bottom Sub-Card: Pre-processing & Gene Filtering -->
    <g id="grp-preprocessing-banner" transform="translate(34, 316)">
      <rect width="636" height="130" rx="8" ry="8" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="1.2" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_NAVY}">Pre-deconvolution Gene Filtering &amp; Signature Selection:</text>
      
      <g transform="translate(14, 30)">
        <rect width="195" height="86" rx="5" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="0.8" />
        <rect x="6" y="6" width="18" height="18" rx="9" fill="#FEE2E2" />
        <text x="15" y="18" font-size="8" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle">1</text>
        <text x="30" y="18" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_DARK}">Outlier Gene Filtering</text>
        <text x="8" y="36" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Removes tumor marker genes</text>
        <text x="8" y="49" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Ribosomal &amp; mitochondrial RNAs</text>
        <text x="8" y="62" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Genes with extreme patient drift</text>
        <text x="8" y="75" font-size="7" fill="{COLOR_CRIMSON}">Prevents tumor leakage into TME</text>
      </g>

      <g transform="translate(220, 30)">
        <rect width="195" height="86" rx="5" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="0.8" />
        <rect x="6" y="6" width="18" height="18" rx="9" fill="#FEF3C7" />
        <text x="15" y="18" font-size="8" font-weight="700" fill="{COLOR_AMBER}" text-anchor="middle">2</text>
        <text x="30" y="18" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_DARK}">Cell-State Sub-Clustering</text>
        <text x="8" y="36" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Resolves sub-states (e.g. CD8 T)</text>
        <text x="8" y="49" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Computes state-level priors α_sg</text>
        <text x="8" y="62" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Reduces condition number κ(Φ)</text>
        <text x="8" y="75" font-size="7" fill="{COLOR_AMBER}">Mitigates collinear sign flipping</text>
      </g>

      <g transform="translate(426, 30)">
        <rect width="195" height="86" rx="5" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="0.8" />
        <rect x="6" y="6" width="18" height="18" rx="9" fill="#DCFCE7" />
        <text x="15" y="18" font-size="8" font-weight="700" fill="{COLOR_EMERALD}" text-anchor="middle">3</text>
        <text x="30" y="18" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_DARK}">Construction of φ_tg^ref</text>
        <text x="8" y="36" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Normalizes to simplex Δ^(G-1)</text>
        <text x="8" y="49" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Linear probability vectors</text>
        <text x="8" y="62" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Establishes baseline expectation</text>
        <text x="8" y="75" font-size="7" fill="{COLOR_EMERALD}">Initialization for MCMC sampling</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: PANEL B — Probabilistic Generative Model
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 03: PANEL B — PROBABILISTIC GENERATIVE MODEL -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Generative_Model">

    <!-- Card: The Multinomial-Dirichlet Mixture Likelihood -->
    <g id="grp-generative-model" transform="translate(729, 114)">
      <rect width="637" height="332" rx="8" ry="8" fill="#FFFBEB" stroke="#FDE68A" stroke-width="1.2" />
      <path d="M 0 8 C 0 3.6 3.6 0 8 0 L 629 0 C 633.4 0 637 3.6 637 8 L 637 28 L 0 28 Z" fill="#FEF3C7" />
      <circle cx="14" cy="14" r="5" fill="{COLOR_AMBER}" />
      <text x="26" y="18" font-size="10.5" font-weight="700" fill="#92400E">Generative Probabilistic Mixture Model (Multinomial Likelihood)</text>
      <rect x="520" y="5" width="105" height="18" rx="9" fill="{COLOR_AMBER}" />
      <text x="572.5" y="17.5" font-size="8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Bayesian Mixture</text>

      <!-- Central Formula Display Box -->
      <g transform="translate(14, 38)">
        <rect width="609" height="76" rx="6" fill="#FFFFFF" stroke="#FDE68A" stroke-width="1" />
        
        <text x="16" y="24" font-size="11" font-weight="700" fill="{COLOR_NAVY}">Bulk Generative Likelihood:</text>
        <text x="16" y="44" font-size="13" font-weight="700" fill="#B45309" font-family="Georgia, serif">
          Y<tspan baseline-shift="sub" font-size="75%">·, n</tspan> ~ Multinomial( N<tspan baseline-shift="sub" font-size="75%">n</tspan>,  ψ<tspan baseline-shift="sub" font-size="75%">n</tspan> )
        </text>
        <text x="16" y="64" font-size="11.5" font-weight="600" fill="{COLOR_TEXT_DARK}" font-family="Georgia, serif">
          where  ψ<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t</tspan>  θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>
        </text>

        <!-- Variable Legend Pills on Right -->
        <g transform="translate(360, 10)">
          <rect width="235" height="56" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_NAVY}">Key Latent Variables:</text>
          <text x="8" y="27" font-size="7" fill="{COLOR_TEXT_MUTED}">• θ_tn ∈ Δ^(T-1) : Fraction of cell type t in sample n</text>
          <text x="8" y="39" font-size="7" fill="{COLOR_TEXT_MUTED}">• φ_tgn ∈ Δ^(G-1) : Expression of gene g in cell type t</text>
          <text x="8" y="50" font-size="7" fill="{COLOR_TEXT_MUTED}">• N_n = ∑ Y_gn : Total library size of bulk sample n</text>
        </g>
      </g>

      <!-- Graphical Model (Bayesian Directed Acyclic Graph - DAG) -->
      <g transform="translate(14, 126)">
        <rect width="609" height="192" rx="6" fill="#FFFFFF" stroke="#FDE68A" stroke-width="1" />
        <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_NAVY}">Probabilistic Graphical Model (Plate Diagram):</text>

        <!-- Plate Diagram Nodes -->
        <!-- Prior Node: alpha_tg -->
        <g transform="translate(60, 40)">
          <circle cx="25" cy="25" r="22" fill="#F0FDFA" stroke="{COLOR_TEAL}" stroke-width="1.8" />
          <text x="25" y="29" font-size="11" font-weight="700" fill="{COLOR_TEAL}" text-anchor="middle" font-family="Georgia, serif">α_tg</text>
          <text x="25" y="60" font-size="7" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">scRNA Prior</text>
        </g>

        <!-- Arrow alpha -> phi -->
        <path d="M 107 65 L 140 65" stroke="{COLOR_NAVY}" stroke-width="1.5" fill="none" marker-end="url(#arrow-navy)" />

        <!-- Latent Node: phi_tgn -->
        <g transform="translate(145, 40)">
          <circle cx="25" cy="25" r="22" fill="#FFFBEB" stroke="{COLOR_AMBER}" stroke-width="1.8" />
          <text x="25" y="29" font-size="11" font-weight="700" fill="{COLOR_AMBER}" text-anchor="middle" font-family="Georgia, serif">φ_tgn</text>
          <text x="25" y="60" font-size="7" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Cell-Type Profile</text>
        </g>

        <!-- Arrow phi -> Z -->
        <path d="M 192 65 L 225 65" stroke="{COLOR_NAVY}" stroke-width="1.5" fill="none" marker-end="url(#arrow-navy)" />

        <!-- Latent Node: Z_tgn (Latent Count Allocation) -->
        <g transform="translate(230, 40)">
          <circle cx="25" cy="25" r="22" fill="#FEE2E2" stroke="{COLOR_CRIMSON}" stroke-width="2" />
          <text x="25" y="29" font-size="12" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle" font-family="Georgia, serif">Z_tgn</text>
          <text x="25" y="60" font-size="7.5" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle">Latent Counts</text>
        </g>

        <!-- Arrow Z -> Y -->
        <path d="M 277 65 L 310 65" stroke="{COLOR_NAVY}" stroke-width="2" fill="none" marker-end="url(#arrow-navy)" />

        <!-- Observed Node: Y_gn (Bulk Observed Counts) -->
        <g transform="translate(315, 40)">
          <circle cx="25" cy="25" r="23" fill="#DBEAFE" stroke="{COLOR_STEEL}" stroke-width="2.5" />
          <text x="25" y="30" font-size="13" font-weight="700" fill="{COLOR_NAVY}" text-anchor="middle" font-family="Georgia, serif">Y_gn</text>
          <text x="25" y="60" font-size="8" font-weight="700" fill="{COLOR_NAVY}" text-anchor="middle">Observed Bulk</text>
        </g>

        <!-- Latent Node: theta_tn (Cell Proportions) -->
        <g transform="translate(145, 120)">
          <circle cx="25" cy="25" r="22" fill="#ECFDF5" stroke="{COLOR_EMERALD}" stroke-width="1.8" />
          <text x="25" y="29" font-size="11" font-weight="700" fill="{COLOR_EMERALD}" text-anchor="middle" font-family="Georgia, serif">θ_tn</text>
          <text x="25" y="58" font-size="7" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Cell Proportions</text>
        </g>

        <!-- Arrow theta -> Z -->
        <path d="M 192 135 L 235 85" stroke="{COLOR_NAVY}" stroke-width="1.5" fill="none" marker-end="url(#arrow-navy)" />

        <!-- Plate Annotations / Explanatory Text Box on Right -->
        <g transform="translate(390, 20)">
          <rect width="205" height="156" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="10" y="18" font-size="8" font-weight="700" fill="{COLOR_NAVY}">The Latent Read Allocation Principle:</text>
          <text x="10" y="34" font-size="7" fill="{COLOR_TEXT_MUTED}">Bulk count Y_gn is the sum of reads</text>
          <text x="10" y="46" font-size="7" fill="{COLOR_TEXT_MUTED}">originating from all cell types t:</text>
          
          <rect x="8" y="54" width="189" height="24" rx="3" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="0.8" />
          <text x="102" y="70" font-size="9.5" font-weight="700" fill="{COLOR_NAVY}" text-anchor="middle" font-family="Georgia, serif">
            Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t</tspan>  Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>
          </text>

          <text x="10" y="94" font-size="7" fill="{COLOR_TEXT_MUTED}">• Linear simplex: Σ_t θ_tn = 1</text>
          <text x="10" y="106" font-size="7" fill="{COLOR_TEXT_MUTED}">• Profile simplex: Σ_g φ_tgn = 1</text>
          <text x="10" y="118" font-size="7" fill="{COLOR_TEXT_MUTED}">• Plates: G genes, T cell types,</text>
          <text x="10" y="130" font-size="7" fill="{COLOR_TEXT_MUTED}">  N clinical bulk samples</text>
          <text x="10" y="146" font-size="7" font-weight="700" fill="{COLOR_EMERALD}">Jointly infers both θ and Z</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: PANEL C — Two-Stage Gibbs Sampling Engine
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 04: PANEL C — TWO-STAGE GIBBS SAMPLING ENGINE -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Gibbs_Sampling">

    <!-- Stage 1 & Stage 2 MCMC Flow -->
    <g id="grp-mcmc-sampling" transform="translate(34, 520)">
      <rect width="637" height="330" rx="8" ry="8" fill="#F5F3FF" stroke="#DDD6FE" stroke-width="1.2" />
      <path d="M 0 8 C 0 3.6 3.6 0 8 0 L 629 0 C 633.4 0 637 3.6 637 8 L 637 28 L 0 28 Z" fill="#EDE9FE" />
      <circle cx="14" cy="14" r="5" fill="{COLOR_PURPLE}" />
      <text x="26" y="18" font-size="10.5" font-weight="700" fill="#5B21B6">Two-Stage MCMC Gibbs Sampling &amp; Latent Inference Workflow</text>
      <rect x="525" y="5" width="100" height="18" rx="9" fill="{COLOR_PURPLE}" />
      <text x="575" y="17.5" font-size="8" font-weight="700" fill="#FFFFFF" text-anchor="middle">MCMC Engine</text>

      <!-- Stage 1 Container -->
      <g transform="translate(14, 38)">
        <rect width="295" height="154" rx="6" fill="#FFFFFF" stroke="#DDD6FE" stroke-width="1" />
        <rect x="8" y="8" width="60" height="16" rx="3" fill="#EDE9FE" />
        <text x="38" y="19" font-size="7.5" font-weight="700" fill="#5B21B6" text-anchor="middle">STAGE 1</text>
        <text x="74" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_DARK}">Coarse Cell Fraction Initialization</text>

        <text x="10" y="38" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Initial estimate of cell fractions θ_tn</text>
        <text x="10" y="50" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Uses static single-cell mean profile φ_tg^ref</text>
        <text x="10" y="62" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Identifies tumor-specific / highly variable genes</text>
        <text x="10" y="74" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Filters genes where bulk expression strongly</text>
        <text x="10" y="86" font-size="7.5" fill="{COLOR_TEXT_MUTED}">  deviates from non-malignant prior expectation</text>

        <rect x="8" y="96" width="279" height="48" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
        <text x="14" y="110" font-size="7.5" font-weight="700" fill="{COLOR_NAVY}">Stage 1 Output:</text>
        <text x="14" y="123" font-size="7" fill="{COLOR_TEXT_MUTED}">Initial θ_tn baseline &amp; robust gene subset G_sub.</text>
        <text x="14" y="135" font-size="7" fill="{COLOR_TEXT_MUTED}">Sets starting point for joint gene expression sampling.</text>
      </g>

      <!-- Stage 2 Container -->
      <g transform="translate(325, 38)">
        <rect width="298" height="154" rx="6" fill="#FFFFFF" stroke="#DDD6FE" stroke-width="1" />
        <rect x="8" y="8" width="60" height="16" rx="3" fill="#EDE9FE" />
        <text x="38" y="19" font-size="7.5" font-weight="700" fill="#5B21B6" text-anchor="middle">STAGE 2</text>
        <text x="74" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_DARK}">Joint Gibbs Sampling (Z &amp; φ)</text>

        <!-- Iterative Gibbs Sampling Loop -->
        <g transform="translate(8, 30)">
          <!-- Step 2A: Sample Z | Y, theta, phi -->
          <rect width="282" height="52" rx="4" fill="#FEF2F2" stroke="#FECACA" stroke-width="0.8" />
          <text x="8" y="14" font-size="8" font-weight="700" fill="{COLOR_CRIMSON}">Step A: Sample Latent Read Counts Z</text>
          <text x="8" y="26" font-size="8.5" font-weight="700" fill="#991B1B" font-family="Georgia, serif">
            Z<tspan baseline-shift="sub" font-size="75%">·, g, n</tspan> | Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>, θ, φ ~ Multinomial( Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>,  π<tspan baseline-shift="sub" font-size="75%">g, n</tspan> )
          </text>
          <text x="8" y="42" font-size="7" fill="{COLOR_TEXT_MUTED}">where allocation probability π<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ∝ θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan></text>

          <!-- Step 2B: Sample phi | alpha, Z -->
          <g transform="translate(0, 58)">
            <rect width="282" height="52" rx="4" fill="#FFFBEB" stroke="#FDE68A" stroke-width="0.8" />
            <text x="8" y="14" font-size="8" font-weight="700" fill="#B45309">Step B: Update Cell-Type Profile φ</text>
            <text x="8" y="26" font-size="8.5" font-weight="700" fill="#B45309" font-family="Georgia, serif">
              φ<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> | Z, α ~ Dirichlet( α<tspan baseline-shift="sub" font-size="75%">t, ·</tspan> + Z<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> )
            </text>
            <text x="8" y="42" font-size="7" fill="{COLOR_TEXT_MUTED}">Dirichlet conjugacy enables closed-form exact posterior update</text>
          </g>
        </g>
      </g>

      <!-- Sampling Convergence Banner -->
      <g transform="translate(14, 204)">
        <rect width="609" height="114" rx="6" fill="#FFFFFF" stroke="#DDD6FE" stroke-width="1" />
        <text x="12" y="18" font-size="9" font-weight="700" fill="{COLOR_NAVY}">MCMC Convergence &amp; Posterior Marginalization:</text>

        <!-- Gibbs Iteration Graph Graphic -->
        <g transform="translate(12, 28)">
          <!-- Trace line schematic -->
          <rect width="180" height="66" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <path d="M 10 48 L 30 20 L 50 38 L 70 26 L 90 32 L 110 28 L 130 30 L 150 29 L 170 30" stroke="{COLOR_PURPLE}" stroke-width="1.8" fill="none" />
          <line x1="70" y1="6" x2="70" y2="60" stroke="#EF4444" stroke-width="1" stroke-dasharray="2,2" />
          <text x="40" y="58" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Burn-in</text>
          <text x="120" y="58" font-size="6.5" fill="{COLOR_PURPLE}">Posterior Sampling</text>
          <text x="90" y="78" font-size="7" font-weight="600" fill="{COLOR_TEXT_DARK}" text-anchor="middle">Gibbs Trace Plot (Convergence)</text>
        </g>

        <!-- Marginalization Equations -->
        <g transform="translate(205, 26)">
          <text x="0" y="14" font-size="8" font-weight="700" fill="{COLOR_TEXT_DARK}">Posterior Expectation across Post-Burn-in Iterations (M):</text>
          
          <rect x="0" y="22" width="390" height="50" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="12" y="38" font-size="8.5" font-weight="700" fill="{COLOR_NAVY}" font-family="Georgia, serif">
            E[θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan>] = (1/M) ∑<tspan baseline-shift="sub" font-size="75%">m</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>   (Posterior Cell Fraction θ*)
          </text>
          <text x="12" y="56" font-size="8.5" font-weight="700" fill="{COLOR_CRIMSON}" font-family="Georgia, serif">
            E[Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>] = (1/M) ∑<tspan baseline-shift="sub" font-size="75%">m</tspan> Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>   (Deconvolved Read Count Matrix Z*)
          </text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: PANEL D — Deconvolved Dual Outputs
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 05: PANEL D — DUAL DECONVOLVED OUTPUTS -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Dual_Outputs">

    <g id="grp-deconvolved-outputs" transform="translate(729, 520)">
      <rect width="637" height="330" rx="8" ry="8" fill="#ECFDF5" stroke="#A7F3D0" stroke-width="1.2" />
      <path d="M 0 8 C 0 3.6 3.6 0 8 0 L 629 0 C 633.4 0 637 3.6 637 8 L 637 28 L 0 28 Z" fill="#D1FAE5" />
      <circle cx="14" cy="14" r="5" fill="{COLOR_EMERALD}" />
      <text x="26" y="18" font-size="10.5" font-weight="700" fill="#065F46">Dual Deconvolved Outputs &amp; Downstream Applications</text>
      <rect x="525" y="5" width="100" height="18" rx="9" fill="{COLOR_EMERALD}" />
      <text x="575" y="17.5" font-size="8" font-weight="700" fill="#FFFFFF" text-anchor="middle">Outputs (θ, Z)</text>

      <!-- Output 1: Cell Type Proportions (Theta) -->
      <g transform="translate(14, 38)">
        <rect width="295" height="140" rx="6" fill="#FFFFFF" stroke="#A7F3D0" stroke-width="1" />
        <rect x="8" y="8" width="68" height="16" rx="3" fill="#ECFDF5" />
        <text x="42" y="19" font-size="7.5" font-weight="700" fill="#047857" text-anchor="middle">OUTPUT 1</text>
        <text x="82" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_DARK}">Cell Type Fractions (θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan>)</text>

        <!-- Proportion Bar Graphic -->
        <g transform="translate(10, 32)">
          <rect width="275" height="20" rx="3" fill="#E2E8F0" />
          <rect x="0" y="0" width="85" height="20" rx="3" fill="#3B82F6" />
          <rect x="87" y="0" width="60" height="20" rx="3" fill="#F59E0B" />
          <rect x="149" y="0" width="45" height="20" rx="3" fill="#10B981" />
          <rect x="196" y="0" width="79" height="20" rx="3" fill="#EF4444" />
          <text x="42" y="13" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">CD8 T (31%)</text>
          <text x="117" y="13" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Mye (22%)</text>
          <text x="171" y="13" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">B (16%)</text>
          <text x="235" y="13" font-size="7" font-weight="700" fill="#FFF" text-anchor="middle">Tumor (29%)</text>
        </g>

        <text x="10" y="70" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Cell type proportions per patient sample n</text>
        <text x="10" y="82" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Sum-to-one constraint: ∑<tspan baseline-shift="sub" font-size="75%">t</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> = 1</text>
        <text x="10" y="94" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Fine state fractions: T_reg, M1, M2, Naive B</text>
        <text x="10" y="106" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Correlates with RECIST immunotherapy response</text>
        <rect x="8" y="114" width="279" height="18" rx="3" fill="#ECFDF5" />
        <text x="147" y="126" font-size="7" font-weight="700" fill="#047857" text-anchor="middle">Standard deconvolution endpoint</text>
      </g>

      <!-- Output 2: Cell-Type-Specific Expression (Z) -->
      <g transform="translate(325, 38)">
        <rect width="298" height="140" rx="6" fill="#FFFFFF" stroke="#A7F3D0" stroke-width="1" />
        <rect x="8" y="8" width="68" height="16" rx="3" fill="#FEE2E2" />
        <text x="42" y="19" font-size="7.5" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle">OUTPUT 2</text>
        <text x="82" y="20" font-size="9" font-weight="700" fill="{COLOR_TEXT_DARK}">Expression Tensor (Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>)</text>

        <!-- Compact 3D Tensor Graphic & Header Bar (y=28 to y=62) -->
        <g transform="translate(10, 28)">
          <rect width="278" height="34" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <!-- Mini Isometric Tensor Slices -->
          <rect x="8" y="5" width="46" height="18" rx="2" fill="#DBEAFE" stroke="#93C5FD" opacity="0.6" />
          <rect x="13" y="9" width="46" height="18" rx="2" fill="#FEF3C7" stroke="#FDE68A" opacity="0.8" />
          <rect x="18" y="13" width="46" height="18" rx="2" fill="#DCFCE7" stroke="#86EFAC" />
          <text x="41" y="25" font-size="6.5" font-weight="700" fill="{COLOR_NAVY}" text-anchor="middle">CD8+ T</text>

          <g transform="translate(74, 4)">
            <text x="0" y="11" font-size="7.5" font-weight="700" fill="{COLOR_NAVY}">3D Count Tensor Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>:</text>
            <text x="0" y="23" font-size="7" fill="{COLOR_TEXT_MUTED}">T cell types × G genes × N samples</text>
          </g>
          <rect x="204" y="8" width="66" height="18" rx="9" fill="#FEE2E2" />
          <text x="237" y="20" font-size="6.5" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle">Unique to BP</text>
        </g>

        <text x="10" y="74" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Deconvolved integer count matrix for each lineage</text>
        <text x="10" y="86" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Removes confounding from cell type abundance</text>
        <text x="10" y="98" font-size="7.5" fill="{COLOR_TEXT_MUTED}">• Recovers patient-specific transcriptional activation</text>
        <rect x="8" y="114" width="282" height="18" rx="3" fill="#FEE2E2" />
        <text x="149" y="126" font-size="7" font-weight="700" fill="{COLOR_CRIMSON}" text-anchor="middle">Enables cell-type-specific differential expression</text>
      </g>

      <!-- Downstream Translational Applications Banner -->
      <g transform="translate(14, 190)">
        <rect width="609" height="128" rx="6" fill="#FFFFFF" stroke="#A7F3D0" stroke-width="1" />
        <text x="12" y="18" font-size="9.5" font-weight="700" fill="{COLOR_NAVY}">Downstream Analytical &amp; Translational Capabilities:</text>

        <!-- 3 Feature Cards -->
        <g transform="translate(10, 28)">
          <rect width="190" height="88" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="8" y="15" font-size="8" font-weight="700" fill="{COLOR_NAVY}">1. Cell-Type Specific DGE</text>
          <text x="8" y="28" font-size="7" fill="{COLOR_TEXT_MUTED}">• Run DESeq2 / edgeR directly</text>
          <text x="8" y="40" font-size="7" fill="{COLOR_TEXT_MUTED}">  on Ẑ_tgn for a target cell type</text>
          <text x="8" y="52" font-size="7" fill="{COLOR_TEXT_MUTED}">• Distinguishes true gene induction</text>
          <text x="8" y="64" font-size="7" fill="{COLOR_TEXT_MUTED}">  from cell number expansion</text>
          <text x="8" y="78" font-size="6.5" font-weight="700" fill="{COLOR_EMERALD}">Solves activation confounding</text>
        </g>

        <g transform="translate(208, 28)">
          <rect width="190" height="88" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="8" y="15" font-size="8" font-weight="700" fill="{COLOR_NAVY}">2. Tumor Heterogeneity</text>
          <text x="8" y="28" font-size="7" fill="{COLOR_TEXT_MUTED}">• Deconvolves malignant profile</text>
          <text x="8" y="40" font-size="7" fill="{COLOR_TEXT_MUTED}">  without requiring patient-matched</text>
          <text x="8" y="52" font-size="7" fill="{COLOR_TEXT_MUTED}">  tumor scRNA-seq references</text>
          <text x="8" y="64" font-size="7" fill="{COLOR_TEXT_MUTED}">• Isolates somatic driver programs</text>
          <text x="8" y="78" font-size="6.5" font-weight="700" fill="{COLOR_CRIMSON}">Captures patient-specific CNA</text>
        </g>

        <g transform="translate(406, 28)">
          <rect width="193" height="88" rx="4" fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.8" />
          <text x="8" y="15" font-size="8" font-weight="700" fill="{COLOR_NAVY}">3. High-Throughput Trials</text>
          <text x="8" y="28" font-size="7" fill="{COLOR_TEXT_MUTED}">• Scales to thousands of bulk</text>
          <text x="8" y="40" font-size="7" fill="{COLOR_TEXT_MUTED}">  immunotherapy trial samples</text>
          <text x="8" y="52" font-size="7" fill="{COLOR_TEXT_MUTED}">• Connects bulk survival (OS/PFS)</text>
          <text x="8" y="64" font-size="7" fill="{COLOR_TEXT_MUTED}">  to single-cell biology</text>
          <text x="8" y="78" font-size="6.5" font-weight="700" fill="{COLOR_PURPLE}">Scalable biomarker discovery</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 6: Flow Connectors & Inter-Panel Arrows
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 06: FLOW CONNECTORS & ARROWS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Flow_Connectors_and_Legends">

    <!-- Top Gutter Horizontal Arrow: Panel A to Panel B -->
    <g id="conn-a-to-b">
      <path d="M 685 200 L 710 200" stroke="{COLOR_AMBER}" stroke-width="2.2" fill="none" marker-end="url(#arrow-amber)" />
      <path d="M 685 380 L 710 380" stroke="{COLOR_NAVY}" stroke-width="2" fill="none" marker-end="url(#arrow-navy)" />
    </g>

    <!-- Bottom Gutter Horizontal Arrow: Panel C to Panel D -->
    <g id="conn-c-to-d">
      <path d="M 685 600 L 725 600" stroke="{COLOR_PURPLE}" stroke-width="2.2" fill="none" marker-end="url(#arrow-purple)" />
      <path d="M 685 750 L 725 750" stroke="{COLOR_EMERALD}" stroke-width="2.2" fill="none" marker-end="url(#arrow-emerald)" />
    </g>
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate BayesPrism mechanism SVG."""
    output_dir = Path("article/figures/deconvolution")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "how_bayesprism_works.svg"
    svg_content = render_svg()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated BayesPrism mechanism SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
