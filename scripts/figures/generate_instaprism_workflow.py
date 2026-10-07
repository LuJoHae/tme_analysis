#!/usr/bin/env python3
"""
Generate Nature Methods Reference Vector Figure: How InstaPrism Works.

Implements the authentic Nature Methods / Nature editorial art style:
- Double-column landscape (180 mm / 1400 px width, 630 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional elements
- Formal display equations with serif mathematics (Georgia/Times) and formal numbering (1) to (6)
- Real coordinate axes with outward tick marks on convergence curves
- Native Inkscape layers and semantic groups
"""

from pathlib import Path
from nature_style_config import (
    WIDTH, HEIGHT,
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)


def render_nature_methods_instaprism() -> str:
    """Render the InstaPrism workflow in pure Nature Methods editorial style."""
    svg_parts: list[str] = []

    # 1. XML Header and Root SVG
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-instaprism-nature-methods"
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
  </defs>
""")

    # ----------------------------------------------------
    # LAYER 0: Canvas Background
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 00: CANVAS BACKGROUND -->
  <g inkscape:groupmode="layer" id="layer-00-canvas" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 1: Editorial Header
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 01: EDITORIAL HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Editorial_Header">
    <g id="grp-header" transform="translate(36, 24)">
      <!-- Main Title -->
      <text id="txt-title" x="0" y="20" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.1">
        InstaPrism: Fast Deterministic Fixed-Point Bayesian Deconvolution
      </text>
      <!-- Subtitle -->
      <text id="txt-subtitle" x="0" y="38" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        Vectorized Expectation-Maximization inverting bulk RNA-seq mixtures into cell fractions (θ) and expression matrices (Z) in seconds
      </text>

      <!-- Academic Reference Box -->
      <g id="grp-reference" transform="translate(1044, 2)">
        <rect x="0" y="0" width="320" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <text x="160" y="18" font-size="9" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
          Korsunsky et al., Nat. Commun. 13, 2022 (Adapted)
        </text>
      </g>
    </g>

    <!-- Hairline Separator -->
    <line x1="36" y1="74" x2="1364" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: Stage 1 — Input Data & Initialization
    # ----------------------------------------------------
    # Column 1: x=36, y=88, width=300, height=518
    svg_parts.append(f"""
  <!-- LAYER 02: STAGE 1 — INPUT DATA & INITIALIZATION -->
  <g inkscape:groupmode="layer" id="layer-02-stage1" inkscape:label="02_Stage1_Inputs_and_Init">
    <g id="col-stage1" transform="translate(36, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="300" height="518" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        1. INPUT DATA &amp; INITIALIZATION
      </text>
      <line x1="14" y1="30" x2="286" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 1A: Single-cell Reference Centroids Matrix -->
      <g id="panel-ref-matrix" transform="translate(14, 40)">
        <rect width="272" height="136" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Reference Centroids Matrix Φ<tspan font-family="{FONT_SERIF_MATH}" font-style="italic">ref</tspan> ∈ [0, 1]<tspan font-size="75%" baseline-shift="super">S × G</tspan>
        </text>

        <!-- Mini Matrix Heatmap -->
        <g transform="translate(10, 26)">
          <rect width="70" height="60" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <line x1="23" y1="0" x2="23" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="46" y1="0" x2="46" y2="60" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="0" y1="30" x2="70" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <rect x="3" y="3" width="17" height="24" fill="{OKABE_SKY_BLUE}" opacity="0.85" />
          <rect x="26" y="33" width="17" height="24" fill="{OKABE_ORANGE}" opacity="0.85" />
          <rect x="49" y="8" width="18" height="49" fill="{OKABE_BLUE}" opacity="0.85" />
          <text x="35" y="70" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
          <text x="-30" y="-3" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">States (S) →</text>
        </g>

        <!-- Specs -->
        <g transform="translate(88, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Reference Properties:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• S cell types/states (e.g. 58 states)</text>
          <text x="0" y="31" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• G signature genes</text>
          <rect x="0" y="38" width="174" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="49" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Simplex Normalized:</text>
          <text x="6" y="58" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> = 1  (Probability column)</text>
        </g>
        <text x="10" y="126" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Single-cell cluster centroids or reference signature matrix.</text>
      </g>

      <!-- Panel 1B: Bulk RNA-seq count vector -->
      <g id="panel-bulk-counts" transform="translate(14, 186)">
        <rect width="272" height="136" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Bulk RNA-seq Expression Vector Y ∈ ℕ₀<tspan font-size="75%" baseline-shift="super">G</tspan>
        </text>

        <!-- Mini Matrix Bulk Bars -->
        <g transform="translate(10, 26)">
          <rect width="70" height="60" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <rect x="5" y="14" width="16" height="46" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="27" y="24" width="16" height="36" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="49" y="8" width="16" height="52" fill="{OKABE_BLUE}" opacity="0.85" />
          <text x="35" y="70" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Genes (G) →</text>
          <text x="-30" y="-3" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Bulk Vector</text>
        </g>

        <!-- Specs -->
        <g transform="translate(88, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Input Properties:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Raw count vector Y<tspan baseline-shift="sub" font-size="75%">g</tspan></text>
          <text x="0" y="31" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Sequencing depth: N = ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> Y<tspan baseline-shift="sub" font-size="75%">g</tspan></text>
          <rect x="0" y="38" width="174" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="49" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Independent Sample Inversion:</text>
          <text x="6" y="58" font-size="6" fill="{COLOR_TEXT_SECONDARY}">Deconvolves each bulk biopsy in milliseconds</text>
        </g>
        <text x="10" y="126" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Enables massive scaling across thousands of clinical trial samples.</text>
      </g>

      <!-- Panel 1C: Initialization & Memory Buffer Allocation -->
      <g id="panel-initialization" transform="translate(14, 332)">
        <rect width="272" height="172" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Vectorized Initialization &amp; Guards:
        </text>

        <!-- Step a: Uniform starting fractions -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">a. Unbiased Simplex Initialization</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">θ<tspan baseline-shift="sub" font-size="75%">s</tspan><tspan baseline-shift="super" font-size="75%">(0)</tspan> = 1 / S for each cell state s ∈ {{1, ..., S}}</text>
          <text x="0" y="30" font-size="6.5" font-weight="600" fill="{OKABE_BLUE}">Unbiased centroid start on simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan></text>
        </g>

        <!-- Step b: Pre-allocated buffers -->
        <g transform="translate(10, 68)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">b. In-Place Pre-allocated Buffers</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">Pre-allocates P (G × S) and Z (G × S) memory</text>
          <text x="0" y="30" font-size="6.5" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">Zero memory allocation inside EM inner loop</text>
        </g>

        <!-- Step c: Numerical guards -->
        <g transform="translate(10, 112)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">c. Numerical Stability Guards</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">Adds machine epsilon (ε = 10⁻¹²) to prevent zero division</text>
          <text x="0" y="30" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Guarantees smooth convex convergence</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: Stage 2 — Probabilistic Model & EM Formulation
    # ----------------------------------------------------
    # Column 2: x=356, y=88, width=320, height=518
    svg_parts.append(f"""
  <!-- LAYER 03: STAGE 2 — EM FORMULATION -->
  <g inkscape:groupmode="layer" id="layer-03-stage2" inkscape:label="03_Stage2_EM_Formulation">
    <g id="col-stage2" transform="translate(356, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="320" height="518" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        2. EM OPTIMIZATION FORMULATION
      </text>
      <line x1="14" y1="30" x2="306" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 2A: Multinomial Mixture Likelihood -->
      <g id="panel-mixture-likelihood" transform="translate(14, 40)">
        <rect width="292" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Bulk Generative Mixture Likelihood:
        </text>

        <!-- Display Equation (1) -->
        <g id="eq-mixture" transform="translate(10, 24)">
          <rect width="272" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="70" y="25" font-size="11" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Y ~ Multinomial( N,  ψ )
          </text>
          <text x="18" y="37" font-size="9" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
            where  ψ<tspan baseline-shift="sub" font-size="75%">g</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan>  Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s</tspan>
          </text>
          <!-- Formal Equation Number (1) -->
          <text x="256" y="26" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(1)</text>
        </g>

        <text x="10" y="80" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Cell fraction simplex: θ<tspan baseline-shift="sub" font-size="75%">s</tspan> ∈ Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan>  (∑<tspan baseline-shift="sub" font-size="75%">s</tspan> θ<tspan baseline-shift="sub" font-size="75%">s</tspan> = 1, θ<tspan baseline-shift="sub" font-size="75%">s</tspan> ≥ 0)</text>
        <text x="10" y="90" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Sequencing library size: N = ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> Y<tspan baseline-shift="sub" font-size="75%">g</tspan></text>
      </g>

      <!-- Panel 2B: Maximum Likelihood Objective & Latent Variable -->
      <g id="panel-em-objective" transform="translate(14, 146)">
        <rect width="292" height="196" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Log-Likelihood Objective &amp; Latent Variable:
        </text>

        <!-- Display Equation (2) -->
        <g id="eq-log-likelihood" transform="translate(10, 24)">
          <rect width="272" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="16" y="26" font-size="10" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            max<tspan baseline-shift="sub" font-size="75%">θ ∈ Δ</tspan>  ℓ(θ) = ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Y<tspan baseline-shift="sub" font-size="75%">g</tspan> log( ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> θ<tspan baseline-shift="sub" font-size="75%">s</tspan> )
          </text>
          <!-- Formal Equation Number (2) -->
          <text x="256" y="26" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(2)</text>
        </g>
        <text x="12" y="76" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Concave optimization problem bounded on probability simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan>.</text>

        <!-- Latent Counts Principle -->
        <g transform="translate(10, 88)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Latent Count Decomposition Principle:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">Bulk reads Y<tspan baseline-shift="sub" font-size="75%">g</tspan> represent the sum of unobserved cell-state reads Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan>:</text>

          <!-- Display Equation (3) -->
          <rect x="0" y="28" width="272" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="56" y="46" font-size="11" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
            Y<tspan baseline-shift="sub" font-size="75%">g</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan>  Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan>,   Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> ∈ ℕ₀
          </text>
          <text x="256" y="46" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(3)</text>

          <text x="0" y="68" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• BayesPrism: Z is sampled stochastically via Gibbs MCMC.</text>
          <text x="0" y="78" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• InstaPrism: Z is replaced by its conditional expectation 𝔼[Z | Y, θ].</text>
          <text x="0" y="88" font-size="6.8" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">Yields deterministic, monotonically converging EM updates.</text>
        </g>
      </g>

      <!-- Panel 2C: Algorithm Comparison (BayesPrism vs InstaPrism) -->
      <g id="panel-comparison" transform="translate(14, 352)">
        <rect width="292" height="152" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Method Paradigm Comparison:
        </text>

        <!-- BayesPrism Block -->
        <g transform="translate(10, 26)">
          <rect width="272" height="52" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{OKABE_REDDISH_PURPLE}">BayesPrism: Stochastic Gibbs MCMC</text>
          <text x="8" y="26" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Draws random samples from Multinomial &amp; Dirichlet posteriors</text>
          <text x="8" y="37" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Requires 1,000s of iterations + burn-in (~30–120s / sample)</text>
          <text x="8" y="47" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_MUTED}">Stochastic convergence; computationally intensive</text>
        </g>

        <!-- InstaPrism Block -->
        <g transform="translate(10, 86)">
          <rect width="272" height="56" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">InstaPrism: Deterministic Fixed-Point EM</text>
          <text x="8" y="26" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Analytic vectorized matrix multiplications strictly in-place</text>
          <text x="8" y="37" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Converges in 50–200 iterations (~0.05s / sample, &gt;500x faster)</text>
          <text x="8" y="48" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Guaranteed monotonic ascent to exact MLE</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: Stage 3 — 3-Step Fixed-Point Iteration Engine
    # ----------------------------------------------------
    # Column 3: x=696, y=88, width=320, height=518
    svg_parts.append(f"""
  <!-- LAYER 04: STAGE 3 — FIXED-POINT EM ENGINE -->
  <g inkscape:groupmode="layer" id="layer-04-stage3" inkscape:label="04_Stage3_Fixed_Point_Engine">
    <g id="col-stage3" transform="translate(696, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="320" height="518" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        3. FIXED-POINT EM ITERATION ENGINE
      </text>
      <line x1="14" y1="30" x2="306" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 3A: 3 In-Place Update Steps per Iteration (t) -->
      <g id="panel-loop-steps" transform="translate(14, 40)">
        <rect width="292" height="306" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          The 3 In-Place Update Steps per Iteration (t):
        </text>

        <!-- STEP 1: Probability Matrix P (E-step) -->
        <g id="step1-p" transform="translate(8, 26)">
          <rect width="276" height="72" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Step 1 (E-Step): Update Posterior Matrix P
          </text>
          
          <!-- Display Equation (4) -->
          <g transform="translate(6, 20)">
            <text x="4" y="16" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
              P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = ( Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> ) / ∑<tspan baseline-shift="sub" font-size="75%">s'</tspan> ( Φ<tspan baseline-shift="sub" font-size="75%">g, s'</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s'</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> )
            </text>
            <text x="250" y="16" font-size="9" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(4)</text>
          </g>
          <text x="8" y="52" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Posterior probability that bulk read of gene g originates from cell state s.</text>
          <text x="8" y="62" font-size="6.5" font-weight="600" fill="{OKABE_BLUE}">Row-stochastic normalization: ∑<tspan baseline-shift="sub" font-size="75%">s</tspan> P<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = 1</text>
        </g>

        <!-- Down connector arrow -->
        <line x1="146" y1="102" x2="146" y2="112" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" marker-end="url(#arr-slate)" />

        <!-- STEP 2: Latent Read Allocation Z -->
        <g id="step2-z" transform="translate(8, 116)">
          <rect width="276" height="72" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Step 2 (Expectation): Allocate Latent Counts Z
          </text>
          
          <!-- Display Equation (5) -->
          <g transform="translate(6, 20)">
            <text x="24" y="16" font-size="10" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
              Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan> · P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan>
            </text>
            <text x="250" y="16" font-size="9" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(5)</text>
          </g>
          <text x="8" y="52" font-size="6.5" fill="{COLOR_TEXT_MUTED}">In-place broadcast elementwise product: np.multiply(bulk, P)</text>
          <text x="8" y="62" font-size="6.5" font-weight="600" fill="{OKABE_VERMILION}">Read conservation invariant strictly holds: ∑<tspan baseline-shift="sub" font-size="75%">s</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan></text>
        </g>

        <!-- Down connector arrow -->
        <line x1="146" y1="192" x2="146" y2="202" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" marker-end="url(#arr-slate)" />

        <!-- STEP 3: M-Step Update Cell Fractions Theta -->
        <g id="step3-theta" transform="translate(8, 206)">
          <rect width="276" height="72" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="14" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
            Step 3 (M-Step): Update Cell Fractions θ
          </text>
          
          <!-- Display Equation (6) -->
          <g transform="translate(6, 20)">
            <text x="6" y="16" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUISH_GREEN}">
              θ<tspan baseline-shift="sub" font-size="75%">s</tspan><tspan baseline-shift="super" font-size="75%">(t+1)</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> / ∑<tspan baseline-shift="sub" font-size="75%">s'</tspan> ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s'</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan>
            </text>
            <text x="250" y="16" font-size="9" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(6)</text>
          </g>
          <text x="8" y="52" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Sum allocated counts across genes for each cell state, then normalize.</text>
          <text x="8" y="62" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Fixed-point update strictly stays on simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan></text>
        </g>

        <!-- Loop Banner -->
        <rect x="8" y="282" width="276" height="18" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="146" y="294" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Iterate Step 1 ⇄ 2 ⇄ 3 until ||θ<tspan baseline-shift="super" font-size="75%">(t+1)</tspan> - θ<tspan baseline-shift="super" font-size="75%">(t)</tspan>||₁ &lt; 10⁻⁷
        </text>
      </g>

      <!-- Panel 3B: Monotonic Convergence & Theoretical Guarantees -->
      <g id="panel-convergence" transform="translate(14, 352)">
        <rect width="292" height="152" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Monotonic EM Convergence Properties:
        </text>

        <!-- Real coordinate axes with outward ticks -->
        <g id="grp-convergence-axis" transform="translate(16, 32)">
          <!-- Y axis line -->
          <line x1="30" y1="10" x2="30" y2="52" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <!-- X axis line -->
          <line x1="30" y1="52" x2="114" y2="52" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />

          <!-- Outward Ticks on Y-axis -->
          <line x1="27" y1="12" x2="30" y2="12" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <line x1="27" y1="32" x2="30" y2="32" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <line x1="27" y1="52" x2="30" y2="52" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <text x="24" y="15" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}" text-anchor="end">max</text>
          <text x="24" y="55" font-size="6" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}" text-anchor="end">ℓ₀</text>

          <!-- Outward Ticks on X-axis -->
          <line x1="30" y1="52" x2="30" y2="55" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <line x1="72" y1="52" x2="72" y2="55" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <line x1="114" y1="52" x2="114" y2="55" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <text x="30" y="62" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
          <text x="72" y="62" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">50</text>
          <text x="114" y="62" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">100</text>

          <!-- Axis Labels -->
          <text x="72" y="70" font-size="6.5" font-weight="500" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">EM Iteration (t)</text>
          <text x="-32" y="14" font-size="6.5" font-weight="500" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle" transform="rotate(-90)">ℓ(θ)</text>

          <!-- Smooth monotonic convergence curve -->
          <path d="M 30 50 Q 45 18 114 12" stroke="{OKABE_BLUE}" stroke-width="1.6" fill="none" />
          <circle cx="114" cy="12" r="2.5" fill="{OKABE_BLUE}" />
        </g>

        <!-- Theoretical Guarantees Text -->
        <g transform="translate(138, 30)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Theoretical Guarantees:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Guaranteed monotonic ascent:</text>
          <text x="8" y="32" font-size="7.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">ℓ(θ<tspan baseline-shift="super" font-size="75%">(t+1)</tspan>) ≥ ℓ(θ<tspan baseline-shift="super" font-size="75%">(t)</tspan>)</text>
          <text x="0" y="44" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Strictly zero Monte Carlo variance</text>
          <text x="0" y="55" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Convex geometry on simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan></text>
        </g>

        <!-- Equivalence Explanation Box -->
        <rect x="10" y="106" width="272" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="16" y="118" font-size="6.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Why InstaPrism Reaches the Same Optimum as BayesPrism:</text>
        <text x="16" y="128" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• The fixed point of this EM engine is the exact Maximum Likelihood Estimator (MLE)</text>
        <text x="16" y="137" font-size="6.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Matches BayesPrism posterior mode without sampling overhead or burn-in.</text>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: Stage 4 — Dual Deconvolved Outputs & Impact
    # ----------------------------------------------------
    # Column 4: x=1036, y=88, width=328, height=518
    svg_parts.append(f"""
  <!-- LAYER 05: STAGE 4 — DUAL OUTPUTS & IMPACT -->
  <g inkscape:groupmode="layer" id="layer-05-stage4" inkscape:label="05_Stage4_Outputs_and_Impact">
    <g id="col-stage4" transform="translate(1036, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="328" height="518" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        4. DECONVOLVED DUAL OUTPUTS &amp; SPEED
      </text>
      <line x1="14" y1="30" x2="314" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 4A: Cell Fractions (Theta) -->
      <g id="panel-output-fractions" transform="translate(14, 40)">
        <rect width="300" height="116" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Output 1: Cell Fractions θ<tspan baseline-shift="sub" font-size="75%">s</tspan> (Simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan>)
        </text>

        <!-- Proportions Stack Bar (Okabe-Ito Colors) -->
        <g transform="translate(10, 26)">
          <rect width="280" height="18" fill="{COLOR_BORDER_HAIRLINE}" />
          <rect x="0" y="0" width="86" height="18" fill="{OKABE_BLUE}" />
          <rect x="86" y="0" width="64" height="18" fill="{OKABE_ORANGE}" />
          <rect x="150" y="0" width="48" height="18" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="198" y="0" width="82" height="18" fill="{OKABE_VERMILION}" />
          <text x="43" y="12" font-size="6.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">CD8 T (31%)</text>
          <text x="118" y="12" font-size="6.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Mye (23%)</text>
          <text x="174" y="12" font-size="6.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">B (17%)</text>
          <text x="239" y="12" font-size="6.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">Tumor (29%)</text>
        </g>

        <!-- Specs -->
        <text x="10" y="58" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Patient-level lineage fractions: ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> θ<tspan baseline-shift="sub" font-size="75%">s</tspan> = 1</text>
        <text x="10" y="70" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Resolves 58+ fine immune, stromal, and malignant sub-states</text>
        <text x="10" y="82" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• High concordance with BayesPrism MCMC posterior mean (ρ &gt; 0.98)</text>
        <text x="10" y="96" font-size="7" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Validated against single-cell Milo neighborhood DA</text>
      </g>

      <!-- Panel 4B: Allocated Expression Matrix (Z) -->
      <g id="panel-output-expression" transform="translate(14, 164)">
        <rect width="300" height="138" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Output 2: Expression Matrix Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> (G × S)
        </text>

        <!-- 2D matrix wireframe graphic -->
        <g transform="translate(10, 26)">
          <rect width="280" height="34" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <rect x="6" y="5" width="46" height="24" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_TEXT_HAIRLINE}" stroke-width="0.5" />
          <line x1="21" y1="5" x2="21" y2="29" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="36" y1="5" x2="36" y2="29" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <rect x="8" y="7" width="11" height="8" fill="{OKABE_VERMILION}" opacity="0.8" />
          <rect x="23" y="18" width="11" height="8" fill="{OKABE_BLUE}" opacity="0.8" />

          <g transform="translate(62, 5)">
            <text x="0" y="11" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Allocated Counts Z (G × S)</text>
            <text x="0" y="21" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Read counts per gene &amp; cell state</text>
          </g>

          <rect x="210" y="8" width="60" height="16" fill="{COLOR_CANVAS_BG}" stroke="{OKABE_VERMILION}" stroke-width="0.75" />
          <text x="240" y="19" font-size="6" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Exact Sum</text>
        </g>

        <!-- Specs -->
        <text x="10" y="74" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Cell-state gene expression: Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan> · P<tspan baseline-shift="sub" font-size="75%">g, s</tspan></text>
        <text x="10" y="86" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Read conservation invariant strictly holds: ∑<tspan baseline-shift="sub" font-size="75%">s</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan></text>
        <text x="10" y="98" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Disentangles cell abundance from cell activation</text>
        <text x="10" y="112" font-size="7" font-weight="600" fill="{OKABE_VERMILION}">Enables lineage-specific differential expression (DESeq2)</text>
      </g>

      <!-- Panel 4C: Computational & Clinical Impact -->
      <g id="panel-impact" transform="translate(14, 310)">
        <rect width="300" height="194" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Downstream Translational Impact:
        </text>

        <!-- 1. Speedup -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. 500x Speedup over MCMC</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• ~0.05 seconds per bulk tumor sample</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Runs in pure vectorized Python / NumPy</text>
          <text x="0" y="41" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Enables real-time cohort re-analysis</text>
        </g>

        <!-- 2. Multi-Cohort Clinical Trials -->
        <g transform="translate(10, 78)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Multi-Cohort Clinical Trials</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Deconvolutes 14 cohorts (9 trials + 5 TCGA)</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Thousands of patients processed in minutes</text>
          <text x="0" y="41" font-size="6.8" font-weight="600" fill="{OKABE_REDDISH_PURPLE}">Scalable biomarker discovery pipeline</text>
        </g>

        <!-- 3. Milo Concordance -->
        <g transform="translate(10, 132)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. Benchmarked Against Milo DA</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Preserves rank stability across inflamed TME (ρ = 0.967)</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Accurately isolates response-associated sub-states</text>
          <text x="0" y="41" font-size="6.8" font-weight="600" fill="{OKABE_BLUISH_GREEN}">High fidelity with single-cell ground truth</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 6: Pipeline Connectors
    # ----------------------------------------------------
    svg_parts.append(f"""
  <!-- LAYER 06: PIPELINE CONNECTORS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Pipeline_Connectors">
    <!-- Connector 1 -> 2 -->
    <path d="M 336 240 L 354 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 336 410 L 354 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 2 -> 3 -->
    <path d="M 676 240 L 694 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 676 410 L 694 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 3 -> 4 -->
    <path d="M 1016 240 L 1034 240" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1016 410 L 1034 410" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the Nature Methods InstaPrism SVG."""
    output_dir = Path("article/figures/deconvolution")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "how_instaprism_works.svg"
    svg_content = render_nature_methods_instaprism()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Nature Methods InstaPrism SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
