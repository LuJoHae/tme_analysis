#!/usr/bin/env python3
"""
Generate Nature Methods Reference Vector Figure: Mathematical Formulations Compendium.

Implements the authentic Nature Methods / Nature editorial art style:
- Full-page landscape (180 mm / 1400 px width, 960 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional mathematical elements
- Formal display equations with serif mathematics (Georgia/Times) and formal numbering (1) to (18)
- Cross-modality mathematical bridge mapping physical single-cell counts to bulk mRNA mass fractions
- Native Inkscape layers and semantic groups
"""

from pathlib import Path
from nature_style_config import (
    WIDTH,
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)

HEIGHT = 960


def render_nature_methods_math() -> str:
    """Render the mathematical formulations comparison SVG in pure Nature Methods editorial style."""
    svg_parts: list[str] = []

    # 1. XML Header and Root SVG
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-mathematical-formulations-nature"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <defs>
    <!-- Minimalist Hairline Arrow Markers -->
    <marker id="arr-slate" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{COLOR_TEXT_SECONDARY}" />
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
  <!-- LAYER 01: FIGURE HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <g id="grp-header" transform="translate(36, 24)">
      <!-- Main Title -->
      <text id="txt-title" x="0" y="20" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.1">
        Mathematical Foundations: BayesPrism, InstaPrism, and Milo / milopy
      </text>
      <!-- Subtitle -->
      <text id="txt-subtitle" x="0" y="38" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        Rigorous statistical formulations comparing Bayesian MCMC joint deconvolution, deterministic fixed-point EM, and graph Negative Binomial GLMs
      </text>

      <!-- Reference Box -->
      <g id="grp-reference" transform="translate(1040, 2)">
        <rect x="0" y="0" width="324" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <text x="162" y="18" font-size="9" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
          Mathematical Compendium &amp; Mapping Bridge
        </text>
      </g>
    </g>

    <!-- Hairline Separator -->
    <line x1="36" y1="74" x2="1364" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
  </g>
""")

    # ----------------------------------------------------
    # LAYER 2: Column 1 — BayesPrism Mathematical Formulation
    # ----------------------------------------------------
    # Column 1: x=36, y=88, width=426, height=696
    svg_parts.append(f"""
  <!-- LAYER 02: COLUMN 1 — BAYESPRISM MATHEMATICAL FORMULATION -->
  <g inkscape:groupmode="layer" id="layer-02-bayesprism" inkscape:label="02_BayesPrism_Formulation">
    <g id="grp-col-bayesprism" transform="translate(36, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="426" height="696" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        BAYESPRISM: STOCHASTIC MCMC GIBBS
      </text>
      <line x1="14" y1="30" x2="412" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Section 1: Likelihood & Simplex Invariants -->
      <g id="bp-sec1-likelihood" transform="translate(14, 40)">
        <rect width="398" height="144" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          1. Multinomial Mixture Likelihood &amp; Simplex Invariants
        </text>

        <!-- Display Equation (1) -->
        <g id="eq-bp-1" transform="translate(10, 24)">
          <rect width="378" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="16" y="24" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Y<tspan baseline-shift="sub" font-size="75%">·, n</tspan> ~ Multinomial( N<tspan baseline-shift="sub" font-size="75%">n</tspan>,  ψ<tspan baseline-shift="sub" font-size="75%">n</tspan> ),   where  ψ<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>
          </text>
          <text x="360" y="24" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(1)</text>
        </g>

        <!-- Simplex constraints -->
        <g transform="translate(10, 70)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Dual Simplex Constraints:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Cell Fractions: θ<tspan baseline-shift="sub" font-size="75%">·, n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="75%">T-1</tspan>  ⟹  ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> = 1,   θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> ≥ 0</text>
          <text x="0" y="32" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Cell Expression: φ<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> ∈ Δ<tspan baseline-shift="super" font-size="75%">G-1</tspan>  ⟹  ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> = 1,   φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ≥ 0</text>
          <text x="0" y="43" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Bulk Sequencing Depth: N<tspan baseline-shift="sub" font-size="75%">n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan> (observed total integer reads)</text>
          <text x="0" y="55" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• T lineages, G signature genes, N bulk tumor samples</text>
        </g>
      </g>

      <!-- Section 2: Prior & Latent Count Allocation -->
      <g id="bp-sec2-prior" transform="translate(14, 192)">
        <rect width="398" height="136" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          2. Dirichlet Prior &amp; Latent Count Allocation
        </text>

        <!-- Prior Formula (2) -->
        <g id="eq-bp-2" transform="translate(10, 24)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            α<tspan baseline-shift="sub" font-size="75%">t, g</tspan> = α<tspan baseline-shift="sub" font-size="75%">0</tspan> + ∑<tspan baseline-shift="sub" font-size="75%">c ∈ C_t</tspan> X<tspan baseline-shift="sub" font-size="75%">c, g</tspan>,    φ<tspan baseline-shift="sub" font-size="75%">t, ·</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> = α<tspan baseline-shift="sub" font-size="75%">t, ·</tspan> / ∑<tspan baseline-shift="sub" font-size="75%">g</tspan> α<tspan baseline-shift="sub" font-size="75%">t, g</tspan>
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(2)</text>
        </g>

        <!-- Latent Sum Formula (3) -->
        <g id="eq-bp-3" transform="translate(10, 60)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="24" y="20" font-size="10" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
            Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">t=1</tspan><tspan baseline-shift="super" font-size="75%">T</tspan> Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>,    Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ∈ ℕ<tspan baseline-shift="sub" font-size="75%">0</tspan>
          </text>
          <text x="360" y="20" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(3)</text>
        </g>

        <g transform="translate(10, 98)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Preserves raw integer counts without artificial non-linear transforms</text>
          <text x="0" y="22" font-size="7" font-weight="600" fill="{OKABE_VERMILION}">• Yields 3D gene expression tensor Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> (T lineages × G genes × N tumors)</text>
        </g>
      </g>

      <!-- Section 3: Joint Gibbs Sampling Engine -->
      <g id="bp-sec3-gibbs" transform="translate(14, 336)">
        <rect width="398" height="210" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          3. Joint Gibbs Sampling Engine (Two-Stage MCMC)
        </text>

        <!-- Joint posterior probability -->
        <g transform="translate(10, 24)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Joint Posterior Distribution:</text>
          <rect x="0" y="16" width="378" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="18" y="33" font-size="8.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            p( Z, θ, Φ | Y, α ) ∝ p( Y | Z ) · p( Z | θ, Φ ) · p( Φ | α ) · p( θ )
          </text>
          <text x="360" y="33" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(4)</text>
        </g>

        <!-- Step A: Sample Z -->
        <g transform="translate(10, 76)">
          <rect width="378" height="52" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Step A: Sample Latent Reads (Z | Y, θ, Φ)</text>
          <text x="16" y="28" font-size="8.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
            Z<tspan baseline-shift="sub" font-size="75%">·, g, n</tspan> | Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>, θ, Φ ~ Multinomial( Y<tspan baseline-shift="sub" font-size="75%">g, n</tspan>,  π<tspan baseline-shift="sub" font-size="75%">g, n</tspan> )
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(5)</text>
          <text x="16" y="42" font-size="7.2" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
            where allocation vector π<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> = ( θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan> ) / ∑<tspan baseline-shift="sub" font-size="75%">t'</tspan> ( θ<tspan baseline-shift="sub" font-size="75%">t', n</tspan> · φ<tspan baseline-shift="sub" font-size="75%">t', g, n</tspan> )
          </text>
        </g>

        <!-- Step B: Sample Phi -->
        <g transform="translate(10, 136)">
          <rect width="378" height="46" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Step B: Update Cell Profile via Conjugacy (Φ | Z, α)</text>
          <text x="16" y="28" font-size="8.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            φ<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> | Z, α ~ Dirichlet( α<tspan baseline-shift="sub" font-size="75%">t, ·</tspan> + Z<tspan baseline-shift="sub" font-size="75%">t, ·, n</tspan> )
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(6)</text>
          <text x="8" y="40" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Dirichlet-Multinomial conjugacy yields exact analytical conditional.</text>
        </g>

        <text x="10" y="198" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Iterate Step A ⇄ Step B for M Gibbs cycles after burn-in.</text>
      </g>

      <!-- Section 4: Posterior Marginalization -->
      <g id="bp-sec4-marginalization" transform="translate(14, 554)">
        <rect width="398" height="128" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          4. Posterior Monte Carlo Marginalization
        </text>

        <rect x="10" y="24" width="378" height="48" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="40" font-size="8.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          𝔼[θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan>] = (1 / M) ∑<tspan baseline-shift="sub" font-size="75%">m=1</tspan><tspan baseline-shift="super" font-size="75%">M</tspan> θ<tspan baseline-shift="sub" font-size="75%">t, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>   (Posterior Cell Fraction θ*)
        </text>
        <text x="12" y="58" font-size="8.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
          𝔼[Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan>] = (1 / M) ∑<tspan baseline-shift="sub" font-size="75%">m=1</tspan><tspan baseline-shift="super" font-size="75%">M</tspan> Z<tspan baseline-shift="sub" font-size="75%">t, g, n</tspan><tspan baseline-shift="super" font-size="75%">(m)</tspan>   (Deconvolved Read Matrix Z*)
        </text>

        <g transform="translate(10, 80)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Evaluated across M stationary MCMC cycles after discarding burn-in</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Runtime: ~30 – 120s per bulk sample due to O(M · G · T) random draws</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• Accurately captures patient-specific somatic malignant CNA profiles</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 3: Column 2 — InstaPrism Mathematical Formulation
    # ----------------------------------------------------
    # Column 2: x=487, y=88, width=426, height=696
    svg_parts.append(f"""
  <!-- LAYER 03: COLUMN 2 — INSTAPRISM MATHEMATICAL FORMULATION -->
  <g inkscape:groupmode="layer" id="layer-03-instaprism" inkscape:label="03_InstaPrism_Formulation">
    <g id="grp-col-instaprism" transform="translate(487, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="426" height="696" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        INSTAPRISM: DETERMINISTIC FIXED-POINT EM
      </text>
      <line x1="14" y1="30" x2="412" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Section 1: Objective & Concave Maximization -->
      <g id="ip-sec1-objective" transform="translate(14, 40)">
        <rect width="398" height="144" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          1. Mixture Log-Likelihood &amp; Simplex Optimization
        </text>

        <!-- Equation 1.1 (7) -->
        <g id="eq-ip-1" transform="translate(10, 24)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Y<tspan baseline-shift="sub" font-size="75%">g</tspan> ~ Multinomial( N,  ψ<tspan baseline-shift="sub" font-size="75%">g</tspan> ),   where  ψ<tspan baseline-shift="sub" font-size="75%">g</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s</tspan>
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(7)</text>
        </g>

        <!-- Concave Objective Formula (8) -->
        <g id="eq-ip-2" transform="translate(10, 60)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="20" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            max<tspan baseline-shift="sub" font-size="75%">θ ∈ Δ^(S-1)</tspan>  ℓ(θ) = ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Y<tspan baseline-shift="sub" font-size="75%">g</tspan> · log( ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s</tspan> )
          </text>
          <text x="360" y="20" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(8)</text>
        </g>

        <g transform="translate(10, 98)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Strictly concave log-likelihood maximized over unit probability simplex</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Replaces random sampling with analytic fixed-point EM map</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• Closed-form expectation 𝔼[Z|Y, θ] guarantees monotonic ascent</text>
        </g>
      </g>

      <!-- Section 2: Conditional Latent Read Allocation -->
      <g id="ip-sec2-allocation" transform="translate(14, 192)">
        <rect width="398" height="136" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          2. Conditional Expectation &amp; Read Conservation
        </text>

        <!-- Read Conservation Invariant (9) -->
        <g id="eq-ip-3" transform="translate(10, 24)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="20" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
            Y<tspan baseline-shift="sub" font-size="75%">g</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">s=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan>,   where  Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = 𝔼[ Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan> | Y<tspan baseline-shift="sub" font-size="75%">g</tspan>, θ ]
          </text>
          <text x="360" y="20" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(9)</text>
        </g>

        <!-- Vectorized Hadamard (10) -->
        <g id="eq-ip-4" transform="translate(10, 60)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="8.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan> · P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan>  ⟹  Z<tspan baseline-shift="super" font-size="75%">(t)</tspan> = Y ⊙ P<tspan baseline-shift="super" font-size="75%">(t)</tspan>  (Vectorized Hadamard)
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(10)</text>
        </g>

        <g transform="translate(10, 98)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Broadcasting in-place product: np.multiply(bulk[:, None], P, out=Z)</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Eliminates sampling noise while matching BayesPrism count output</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_VERMILION}">• Outputs full G × S expression matrix for cell-state differential analysis</text>
        </g>
      </g>

      <!-- Section 3: 3-Step Fixed-Point EM Engine -->
      <g id="ip-sec3-engine" transform="translate(14, 336)">
        <rect width="398" height="210" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          3. 3-Step In-Place Fixed-Point Iteration Engine
        </text>

        <!-- Step 1: E-Step (11) -->
        <g transform="translate(10, 24)">
          <rect width="378" height="46" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Step 1 (E-Step): Posterior Probability Matrix (P)</text>
          <text x="12" y="28" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = ( Φ<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> ) / ∑<tspan baseline-shift="sub" font-size="75%">s'=1</tspan><tspan baseline-shift="super" font-size="75%">S</tspan> ( Φ<tspan baseline-shift="sub" font-size="75%">g, s'</tspan><tspan baseline-shift="super" font-size="75%">ref</tspan> · θ<tspan baseline-shift="sub" font-size="75%">s'</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> )
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(11)</text>
          <text x="8" y="40" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Row-stochastic normalization: ∑<tspan baseline-shift="sub" font-size="75%">s</tspan> P<tspan baseline-shift="sub" font-size="75%">g, s</tspan> = 1</text>
        </g>

        <!-- Step 2: Expectation -->
        <g transform="translate(10, 76)">
          <rect width="378" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Step 2 (Expectation): Allocate Reads to Cell States (Z)</text>
          <text x="12" y="28" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_VERMILION}">
            Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = Y<tspan baseline-shift="sub" font-size="75%">g</tspan> · P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan>
          </text>
          <text x="8" y="38" font-size="6.5" fill="{COLOR_TEXT_MUTED}">In-place vector multiplication without intermediate array allocations</text>
        </g>

        <!-- Step 3: M-Step (12) -->
        <g transform="translate(10, 124)">
          <rect width="378" height="48" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Step 3 (M-Step): Simplex Fraction Update (θ)</text>
          <text x="12" y="28" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUISH_GREEN}">
            θ<tspan baseline-shift="sub" font-size="75%">s</tspan><tspan baseline-shift="super" font-size="75%">(t+1)</tspan> = (1 / N) ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Z<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan> = (1 / N) ∑<tspan baseline-shift="sub" font-size="75%">g=1</tspan><tspan baseline-shift="super" font-size="75%">G</tspan> Y<tspan baseline-shift="sub" font-size="75%">g</tspan> P<tspan baseline-shift="sub" font-size="75%">g, s</tspan><tspan baseline-shift="super" font-size="75%">(t)</tspan>
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(12)</text>
          <text x="8" y="42" font-size="6.5" fill="{COLOR_TEXT_MUTED}">Column sum over genes normalized to simplex Δ<tspan baseline-shift="super" font-size="75%">S-1</tspan></text>
        </g>

        <!-- Repeat text -->
        <rect x="10" y="178" width="378" height="20" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="199" y="191" font-size="6.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
          Iterate Step 1 ⇄ 2 ⇄ 3 until ||θ<tspan baseline-shift="super" font-size="75%">(t+1)</tspan> - θ<tspan baseline-shift="super" font-size="75%">(t)</tspan>||₁ &lt; 10⁻⁷
        </text>
      </g>

      <!-- Section 4: Monotonic Convergence & Speed -->
      <g id="ip-sec4-convergence" transform="translate(14, 554)">
        <rect width="398" height="128" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          4. Monotonic Convergence Guarantee &amp; Scalability
        </text>

        <rect x="10" y="24" width="378" height="48" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="40" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
          ℓ( θ<tspan baseline-shift="super" font-size="75%">(t+1)</tspan> ) ≥ ℓ( θ<tspan baseline-shift="super" font-size="75%">(t)</tspan> ),   with lim<tspan baseline-shift="sub" font-size="75%">t→∞</tspan> θ<tspan baseline-shift="super" font-size="75%">(t)</tspan> = θ<tspan baseline-shift="sub" font-size="75%">MLE</tspan> = argmax<tspan baseline-shift="sub" font-size="75%">θ ∈ Δ</tspan> ℓ(θ)
        </text>
        <text x="12" y="58" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">
          Concordance with BayesPrism MCMC: Spearman ρ &gt; 0.98 across all cell states
        </text>

        <g transform="translate(10, 80)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Monotonic ascent guaranteed by Baum-Welch / EM fixed-point theorem</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Runtime: ~0.05s per bulk sample (&gt;500x speedup over MCMC sampling)</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• Enables deconvolution of 10,000+ patient clinical trial cohorts in minutes</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 4: Column 3 — Milo / milopy Mathematical Formulation
    # ----------------------------------------------------
    # Column 3: x=938, y=88, width=426, height=696
    svg_parts.append(f"""
  <!-- LAYER 04: COLUMN 3 — MILO / MILOPY MATHEMATICAL FORMULATION -->
  <g inkscape:groupmode="layer" id="layer-04-milopy" inkscape:label="04_Milo_Formulation">
    <g id="grp-col-milopy" transform="translate(938, 88)">
      <!-- Wireframe Outer Frame -->
      <rect width="426" height="696" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        MILO / MILOPY: GRAPH NEGATIVE BINOMIAL GLM
      </text>
      <line x1="14" y1="30" x2="412" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Section 1: kNN Graph & Neighborhood Incidence -->
      <g id="milo-sec1-graph" transform="translate(14, 40)">
        <rect width="398" height="144" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          1. kNN Graph Manifold &amp; Medoid Refinement
        </text>

        <!-- Graph Definition (13) -->
        <g id="eq-milo-1" transform="translate(10, 24)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            G = (V, E),   E = {{ (u, v) | v ∈ kNN(u) in ℝ<tspan baseline-shift="super" font-size="75%">d</tspan> }},   d = 30,  k = 30
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(13)</text>
        </g>

        <!-- Medoid Refinement (14) -->
        <g id="eq-milo-2" transform="translate(10, 60)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            v* = argmin<tspan baseline-shift="sub" font-size="75%">u ∈ N(v)</tspan> || x<tspan baseline-shift="sub" font-size="75%">u</tspan> - median({{ x<tspan baseline-shift="sub" font-size="75%">w</tspan> | w ∈ N(v) }}) ||<tspan baseline-shift="sub" font-size="75%">2</tspan>
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(14)</text>
        </g>

        <g transform="translate(10, 98)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Neighborhood incidence: M<tspan baseline-shift="sub" font-size="75%">c, v</tspan> = 𝕀( c ∈ N(v) ),   M ∈ {{0, 1}}<tspan baseline-shift="super" font-size="75%">C × V</tspan></text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Shifts random sample vertices to medoids for uniform coverage</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUISH_GREEN}">• Overlapping continuous ensembles avoid discrete clustering resolution bias</text>
        </g>
      </g>

      <!-- Section 2: Patient Cell Counting -->
      <g id="milo-sec2-counting" transform="translate(14, 192)">
        <rect width="398" height="136" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          2. Patient Cell Counting &amp; Biological Replication
        </text>

        <!-- Cell Counting Formula (15) -->
        <g id="eq-milo-3" transform="translate(10, 24)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            N<tspan baseline-shift="sub" font-size="75%">v, j</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">c ∈ sample j</tspan>  M<tspan baseline-shift="sub" font-size="75%">c, v</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">c ∈ sample j</tspan>  𝕀( c ∈ N(v) )
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(15)</text>
        </g>

        <!-- Library Size Offset (16) -->
        <g id="eq-milo-4" transform="translate(10, 60)">
          <rect width="378" height="30" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="12" y="19" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            offset<tspan baseline-shift="sub" font-size="75%">j</tspan> = log( TotalCells<tspan baseline-shift="sub" font-size="75%">j</tspan> ) + log( TMM<tspan baseline-shift="sub" font-size="75%">j</tspan> )
          </text>
          <text x="360" y="19" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(16)</text>
        </g>

        <g transform="translate(10, 98)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Eliminates pseudo-replication by treating patient samples (J) as true replicates</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Effective sequencing library size normalized via Trimmed Mean of M-values</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• Accommodates complex design formulas: `~ condition + batch + age`</text>
        </g>
      </g>

      <!-- Section 3: Quasi-Likelihood Negative Binomial GLM -->
      <g id="milo-sec3-glm" transform="translate(14, 336)">
        <rect width="398" height="210" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          3. Negative Binomial Quasi-Likelihood GLM
        </text>

        <!-- GLM Mean Formula (17) -->
        <g transform="translate(10, 24)">
          <rect width="378" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">GLM Log-Link Model:</text>
          <text x="12" y="28" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            log( 𝔼[N<tspan baseline-shift="sub" font-size="75%">v, j</tspan>] ) = β<tspan baseline-shift="sub" font-size="75%">0, v</tspan> + β<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">milo</tspan> · Y<tspan baseline-shift="sub" font-size="75%">j</tspan> + log( TotalCells<tspan baseline-shift="sub" font-size="75%">j</tspan> )
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(17)</text>
        </g>

        <!-- Variance & Dispersion -->
        <g transform="translate(10, 72)">
          <rect width="378" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Variance &amp; Overdispersion Structure:</text>
          <text x="12" y="28" font-size="8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_SECONDARY}">
            Var( N<tspan baseline-shift="sub" font-size="75%">v, j</tspan> ) = μ<tspan baseline-shift="sub" font-size="75%">v, j</tspan> + φ<tspan baseline-shift="sub" font-size="75%">v</tspan> · μ<tspan baseline-shift="sub" font-size="75%">v, j</tspan><tspan baseline-shift="super" font-size="75%">2</tspan>,   d<tspan baseline-shift="sub" font-size="75%">v, j</tspan> = σ<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">2</tspan> · ( μ<tspan baseline-shift="sub" font-size="75%">v, j</tspan> + φ<tspan baseline-shift="sub" font-size="75%">v</tspan> μ<tspan baseline-shift="sub" font-size="75%">v, j</tspan><tspan baseline-shift="super" font-size="75%">2</tspan> )
          </text>
        </g>

        <!-- QL F-Test (18) -->
        <g transform="translate(10, 120)">
          <rect width="378" height="54" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Quasi-Likelihood F-Test &amp; Spatial FDR:</text>
          <text x="12" y="28" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            F<tspan baseline-shift="sub" font-size="75%">v</tspan> = [ ( Dev<tspan baseline-shift="sub" font-size="75%">reduced, v</tspan> - Dev<tspan baseline-shift="sub" font-size="75%">full, v</tspan> ) / Δdf ] / σ̂<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">2</tspan>
          </text>
          <text x="360" y="28" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(18)</text>
          <text x="12" y="44" font-size="7.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
            w<tspan baseline-shift="sub" font-size="75%">v</tspan> ∝ 1 / connectivity_density(v),   p<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">spatial</tspan> = Weighted-BH( p<tspan baseline-shift="sub" font-size="75%">v</tspan>, w<tspan baseline-shift="sub" font-size="75%">v</tspan> )
          </text>
        </g>

        <text x="10" y="196" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Weighted Benjamini-Hochberg corrects for correlated tests across overlapping nhoods.</text>
      </g>

      <!-- Section 4: Effect Size & Cluster Summary -->
      <g id="milo-sec4-effect" transform="translate(14, 554)">
        <rect width="398" height="128" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          4. Neighborhood &amp; Cluster-Level Effect Sizes
        </text>

        <rect x="10" y="24" width="378" height="48" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="40" font-size="8.8" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          logFC<tspan baseline-shift="sub" font-size="75%">v</tspan> = β<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">milo</tspan> / log(2),   logFC<tspan baseline-shift="sub" font-size="75%">k</tspan> = median<tspan baseline-shift="sub" font-size="75%">v ∈ C_k</tspan>( logFC<tspan baseline-shift="sub" font-size="75%">v</tspan> )
        </text>
        <text x="12" y="58" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">
          Validated in Sade-Feldman: Naive B cells (+2.066), M2 Macrophages (-1.562)
        </text>

        <g transform="translate(10, 80)">
          <text x="0" y="10" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Evaluates physical cell numbers, inherently independent of cell size S<tspan baseline-shift="sub" font-size="75%">k</tspan></text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Provides single-cell ground truth for calibrating bulk deconvolution</text>
          <text x="0" y="32" font-size="7" font-weight="600" fill="{OKABE_BLUE}">• High rank stability across inflamed and perturbed TME regimes (ρ = 0.967)</text>
        </g>
      </g>
    </g>
  </g>
""")

    # ----------------------------------------------------
    # LAYER 5: Bottom Banner — Cross-Modality Mathematical Bridge
    # ----------------------------------------------------
    # Width: 1328px (x=36, y=798, height=144)
    svg_parts.append(f"""
  <!-- LAYER 05: CROSS-MODALITY MATHEMATICAL BRIDGE -->
  <g inkscape:groupmode="layer" id="layer-05-bridge" inkscape:label="05_Cross_Modality_Bridge">
    <g id="grp-cross-bridge" transform="translate(36, 798)">
      <!-- Banner Container -->
      <rect width="1328" height="144" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="14" y="20" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        CROSS-MODALITY MATHEMATICAL BRIDGE: SINGLE-CELL MILO DA (COUNTS) ⇄ BULK DECONVOLUTION (mRNA FRACTIONS)
      </text>
      <line x1="14" y1="28" x2="1314" y2="28" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- 3 Bridge Cards -->
      <!-- Card 1: mRNA Mass Decoupling Equation -->
      <g transform="translate(14, 38)">
        <rect width="420" height="92" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. mRNA Mass Decoupling (Cell Size Bias S<tspan baseline-shift="sub" font-size="75%">k</tspan>):</text>
        <rect x="8" y="24" width="404" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="20" y="42" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          f<tspan baseline-shift="sub" font-size="75%">k, j</tspan><tspan baseline-shift="super" font-size="75%">mRNA</tspan> = ( p<tspan baseline-shift="sub" font-size="75%">k, j</tspan> · S<tspan baseline-shift="sub" font-size="75%">k</tspan> ) / ∑<tspan baseline-shift="sub" font-size="75%">m=1</tspan><tspan baseline-shift="super" font-size="75%">K</tspan> ( p<tspan baseline-shift="sub" font-size="75%">m, j</tspan> · S<tspan baseline-shift="sub" font-size="75%">m</tspan> )
        </text>
        <text x="10" y="66" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• p<tspan baseline-shift="sub" font-size="75%">k, j</tspan> = physical cell fraction (tested by Milo); S<tspan baseline-shift="sub" font-size="75%">k</tspan> = per-cell mRNA yield.</text>
        <text x="10" y="78" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">• Large cells (plasma, macrophages, 10–50x mRNA) decouple f<tspan baseline-shift="super" font-size="75%">mRNA</tspan> from p<tspan baseline-shift="sub" font-size="75%">k</tspan>.</text>
      </g>

      <!-- Card 2: Standardized Point-Biserial Effect Size -->
      <g transform="translate(454, 38)">
        <rect width="420" height="92" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Standardized Point-Biserial Deconvolution Effect Size:</text>
        <rect x="8" y="24" width="404" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="16" y="42" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
          r<tspan baseline-shift="sub" font-size="75%">pb, k</tspan> = Corr( Y<tspan baseline-shift="sub" font-size="75%">j</tspan>, f̂<tspan baseline-shift="sub" font-size="75%">k, j</tspan> ),    β̂<tspan baseline-shift="sub" font-size="75%">k</tspan><tspan baseline-shift="super" font-size="75%">deconv</tspan> = ( 2 · r<tspan baseline-shift="sub" font-size="75%">pb, k</tspan> ) / √( 1 - r<tspan baseline-shift="sub" font-size="75%">pb, k</tspan><tspan baseline-shift="super" font-size="75%">2</tspan> + ε )
        </text>
        <text x="10" y="66" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Standardizes patient deconvolution fractions f̂<tspan baseline-shift="sub" font-size="75%">k, j</tspan> against binary response Y<tspan baseline-shift="sub" font-size="75%">j</tspan>.</text>
        <text x="10" y="78" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">• Places bulk deconvolution on identical numerical scale as Milo logFC<tspan baseline-shift="sub" font-size="75%">k</tspan>.</text>
      </g>

      <!-- Card 3: Directional Concordance Quadrants -->
      <g transform="translate(894, 38)">
        <rect width="420" height="92" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">3. Four-Quadrant Concordance Classification Criterion:</text>
        <rect x="8" y="24" width="404" height="28" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="12" y="37" font-size="7.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}">
          Concordant Responder:  β̂<tspan baseline-shift="sub" font-size="75%">k</tspan><tspan baseline-shift="super" font-size="75%">deconv</tspan> &gt; +0.1  ∧  logFC<tspan baseline-shift="sub" font-size="75%">k</tspan> &gt; +0.1
        </text>
        <text x="12" y="47" font-size="7.5" font-weight="600" fill="{OKABE_VERMILION}">
          Concordant Non-Responder:  β̂<tspan baseline-shift="sub" font-size="75%">k</tspan><tspan baseline-shift="super" font-size="75%">deconv</tspan> &lt; -0.1  ∧  logFC<tspan baseline-shift="sub" font-size="75%">k</tspan> &lt; -0.1
        </text>
        <text x="10" y="66" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Discordant: opposite sign agreement due to state activation or ghost cells.</text>
        <text x="10" y="78" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">• High global concordance in melanoma: Spearman ρ = 0.853 (p = 4.18 × 10⁻⁴).</text>
      </g>
    </g>
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the Nature Methods mathematical formulations SVG."""
    output_dir = Path("article/figures/deconvolution")
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_file = output_dir / "mathematical_formulations_comparison.svg"
    svg_content = render_nature_methods_math()

    with open(svg_file, "w", encoding="utf-8") as f:
        f.write(svg_content)

    print(f"Successfully generated Nature Methods math formulations SVG at: {svg_file}")
    print(f"File size: {len(svg_content):,} bytes")


if __name__ == "__main__":
    main()
