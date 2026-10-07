#!/usr/bin/env python3
"""
Generate Nature Methods Reference Vector Figure: How Milo / milopy Works.

Implements the authentic Nature Methods / Nature editorial art style:
- Double-column landscape (180 mm / 1400 px width, 680 px height)
- Pure minimalist wireframe (rx=0, pure white canvas, zero drop shadows, 0.5-0.75 pt hairlines)
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to functional elements
- Formal display equations with serif mathematics (Times/Georgia) and formal numbering (1) to (4)
- Rich visual diagrams in all 4 stages:
  1. kNN graph & index sampling with direct visual contrast: Discrete Clusters vs Continuous Milo
  2. Overlapping neighborhoods, sparse incidence matrix, and cross-cohort purity/composition (C = B^T A)
  3. Patient count table, TMM normalization, and Pseudo-Replication Fallacy vs Hierarchical GLM
  4. Quasi-Likelihood GLM, Spatial FDR weighting, and dual outputs (Mini Volcano + UMAP logFC)
- Native Inkscape layers and semantic groups
- Automatic 300 DPI PNG rasterization via vl-convert or rsvg-convert
"""

from __future__ import annotations

import html
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

try:
    import vl_convert as vlc  # type: ignore
except ImportError:
    vlc = None

from nature_style_config import (
    COLOR_CANVAS_BG, COLOR_PANEL_BG, COLOR_BORDER_HAIRLINE, COLOR_DIVIDER_RULE, COLOR_SUBTLE_FILL,
    COLOR_TEXT_PRIMARY, COLOR_TEXT_SECONDARY, COLOR_TEXT_MUTED, COLOR_TEXT_HAIRLINE,
    OKABE_BLACK, OKABE_ORANGE, OKABE_SKY_BLUE, OKABE_BLUISH_GREEN,
    OKABE_BLUE, OKABE_VERMILION, OKABE_REDDISH_PURPLE,
    FONT_SANS, FONT_SERIF_MATH
)

# Canvas Dimensions (Double-column landscape, 180 mm / 1400 px width, 680 px height)
WIDTH = 1400
HEIGHT = 680


def render_nature_methods_milopy() -> str:
    """Render the Milo / milopy workflow in pure Nature Methods editorial style."""
    svg_parts: list[str] = []

    # 1. XML Header and Root SVG
    svg_parts.append(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://www.sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-milopy-nature-methods"
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
    <g id="grp-header" transform="translate(36, 20)">
      <!-- Main Title -->
      <text id="txt-title" x="0" y="20" font-size="15" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.1">
        Milo / milopy: Continuous Differential Abundance on Single-Cell k-Nearest Neighbor Graphs
      </text>
      <!-- Subtitle -->
      <text id="txt-subtitle" x="0" y="38" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
        Non-parametric neighborhood testing bypassing discrete clustering boundaries using Negative Binomial quasi-likelihood GLMs
      </text>

      <!-- Academic Reference Box -->
      <g id="grp-reference" transform="translate(1040, 6)">
        <rect x="0" y="0" width="324" height="28" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
        <text x="162" y="18" font-size="9" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
          Dann et al., Nat. Biotechnol. 40, 245–253 (2022)
        </text>
      </g>
    </g>

    <!-- Hairline Separator -->
    <line x1="36" y1="66" x2="1364" y2="66" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
  </g>
""")

    # =========================================================================
    # LAYER 2: Stage 1 — kNN Graph & Index Vertex Sampling
    # Column 1: x=36, y=78, width=316, height=586
    # =========================================================================
    svg_parts.append(f"""
  <!-- LAYER 02: STAGE 1 — KNN GRAPH & INDEX SAMPLING -->
  <g inkscape:groupmode="layer" id="layer-02-stage1" inkscape:label="02_Stage1_Graph_and_Sampling">
    <g id="col-stage1" transform="translate(36, 78)">
      <!-- Wireframe Outer Frame -->
      <rect width="316" height="586" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        1. kNN GRAPH &amp; INDEX SAMPLING
      </text>
      <line x1="14" y1="30" x2="302" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 1A: kNN Graph Construction -->
      <g id="panel-knn-graph" transform="translate(14, 38)">
        <rect width="288" height="152" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          k-Nearest Neighbor Graph G = (V, E)
        </text>

        <!-- Mini kNN Graph Diagram -->
        <g transform="translate(10, 26)">
          <rect width="78" height="74" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <!-- Graph Edges -->
          <line x1="16" y1="18" x2="38" y2="12" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="38" y1="12" x2="62" y2="22" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="16" y1="18" x2="24" y2="48" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="38" y1="12" x2="46" y2="40" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="62" y1="22" x2="58" y2="52" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="24" y1="48" x2="46" y2="40" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="46" y1="40" x2="58" y2="52" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="24" y1="48" x2="32" y2="66" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="46" y1="40" x2="52" y2="64" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />
          <line x1="58" y1="52" x2="52" y2="64" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1" />

          <!-- Graph Nodes (Color coded by phenotypic states) -->
          <circle cx="16" cy="18" r="3.5" fill="{OKABE_SKY_BLUE}" />
          <circle cx="38" cy="12" r="4.5" fill="{OKABE_BLUE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <circle cx="62" cy="22" r="3.5" fill="{OKABE_SKY_BLUE}" />
          <circle cx="24" cy="48" r="3.5" fill="{OKABE_BLUISH_GREEN}" />
          <circle cx="46" cy="40" r="4.5" fill="{OKABE_BLUE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <circle cx="58" cy="52" r="3.5" fill="{OKABE_VERMILION}" />
          <circle cx="32" cy="66" r="3.5" fill="{OKABE_BLUISH_GREEN}" />
          <circle cx="52" cy="64" r="4.5" fill="{OKABE_VERMILION}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />

          <text x="39" y="82" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">PCA Latent Manifold</text>
        </g>

        <!-- Hyperparameter Specs -->
        <g transform="translate(96, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Graph Hyperparameters:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Latent PCs: d = 30 dimensions</text>
          <text x="0" y="31" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Nearest neighbors: k = 30</text>
          <text x="0" y="41" font-size="7" fill="{COLOR_TEXT_SECONDARY}">• Metric: Euclidean on PCs</text>
          <rect x="0" y="48" width="180" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="59" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Trajectory Preservation:</text>
          <text x="6" y="68" font-size="6" fill="{COLOR_TEXT_SECONDARY}">Preserves continuous phenotypic transitions</text>
        </g>
        <text x="10" y="142" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Constructed via `milopy.core.build_graph(adata, k=30, d=30)`.</text>
      </g>

      <!-- Panel 1B: Discrete Clusters vs Continuous Milo -->
      <g id="panel-clusters-vs-milo" transform="translate(14, 198)">
        <rect width="288" height="190" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Discrete Clustering vs. Continuous Milo:
        </text>

        <!-- Graphic: Comparison Schematics -->
        <g transform="translate(10, 24)">
          <!-- Sub-panel Left: Discrete Hard Clusters -->
          <g transform="translate(0, 0)">
            <rect width="130" height="66" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="65" y="11" font-size="6.5" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Discrete Clusters (Leiden)</text>
            <!-- Hard dividing line -->
            <line x1="65" y1="16" x2="65" y2="60" stroke="{OKABE_VERMILION}" stroke-width="1" stroke-dasharray="2,2" />
            <circle cx="35" cy="38" r="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="95" cy="38" r="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <!-- Intermediate state cut in half -->
            <circle cx="65" cy="38" r="4.5" fill="{OKABE_ORANGE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
            <text x="35" y="41" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Clust 1</text>
            <text x="95" y="41" font-size="6" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Clust 2</text>
            <text x="65" y="64" font-size="5.5" font-weight="600" fill="{OKABE_VERMILION}" text-anchor="middle">Dilutes boundary shifts</text>
          </g>

          <!-- Sub-panel Right: Continuous Milo Neighborhoods -->
          <g transform="translate(138, 0)">
            <rect width="130" height="66" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="65" y="11" font-size="6.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Continuous Milo Nhoods</text>
            <!-- Overlapping circles spanning the boundary -->
            <circle cx="45" cy="38" r="18" fill="none" stroke="{OKABE_BLUE}" stroke-width="1" stroke-dasharray="2,2" />
            <circle cx="85" cy="38" r="18" fill="none" stroke="{OKABE_BLUISH_GREEN}" stroke-width="1" stroke-dasharray="2,2" />
            <circle cx="65" cy="38" r="16" fill="none" stroke="{OKABE_ORANGE}" stroke-width="1.2" />
            <circle cx="65" cy="38" r="4" fill="{OKABE_ORANGE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
            <text x="65" y="64" font-size="5.5" font-weight="600" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Direct local detection</text>
          </g>
        </g>

        <!-- Descriptive comparative text -->
        <g transform="translate(10, 102)">
          <text x="0" y="8" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Fundamental Advantages of Milo:</text>
          <text x="0" y="20" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">1. <tspan font-weight="700">No Hard Boundaries:</tspan> Captures cell state continua (e.g. CD8 exhaustion)</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">2. <tspan font-weight="700">Resolution-Invariant:</tspan> Overcomes Leiden clustering resolution sensitivity</text>
          <text x="0" y="44" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">3. <tspan font-weight="700">High Statistical Power:</tspan> Avoids diluting effect size across broad clusters</text>
          <text x="0" y="58" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Detects subtle immune remodelings missed by discrete clustering</text>
        </g>
      </g>

      <!-- Panel 1C: Index Vertex Sampling & Medoid Refinement -->
      <g id="panel-sampling" transform="translate(14, 396)">
        <rect width="288" height="176" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Index Sampling &amp; Medoid Refinement
        </text>

        <!-- Diagram: Random vs Refined -->
        <g transform="translate(10, 26)">
          <rect width="78" height="66" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <!-- Target cloud -->
          <circle cx="39" cy="32" r="22" fill="{COLOR_DIVIDER_RULE}" opacity="0.4" />
          <circle cx="26" cy="22" r="2.5" fill="{COLOR_TEXT_HAIRLINE}" />
          <circle cx="50" cy="24" r="2.5" fill="{COLOR_TEXT_HAIRLINE}" />
          <circle cx="32" cy="44" r="2.5" fill="{COLOR_TEXT_HAIRLINE}" />
          <circle cx="48" cy="42" r="2.5" fill="{COLOR_TEXT_HAIRLINE}" />
          <!-- Initial Random vs Refined Node -->
          <circle cx="22" cy="36" r="3.5" fill="{OKABE_ORANGE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
          <line x1="26" y1="35" x2="34" y2="32" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" marker-end="url(#arr-slate)" />
          <circle cx="39" cy="32" r="4.5" fill="{OKABE_BLUE}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="1" />
          <text x="39" y="60" font-size="6" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Medoid Shift</text>
        </g>

        <!-- Specs -->
        <g transform="translate(96, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Two-Step Sampling Protocol:</text>
          <text x="0" y="21" font-size="7" fill="{COLOR_TEXT_SECONDARY}">1. Random draw: prop = 0.1 (10% cells)</text>
          <text x="0" y="31" font-size="7" fill="{COLOR_TEXT_SECONDARY}">2. Refine: shift to median of kNN profile</text>
          <rect x="0" y="38" width="180" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="49" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Manifold Invariance:</text>
          <text x="6" y="58" font-size="6" fill="{COLOR_TEXT_SECONDARY}">Uniform sampling across dense &amp; sparse states</text>
        </g>
        <text x="10" y="112" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Eliminates sampling density bias across diverse cell populations.</text>
        <text x="10" y="124" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Selects representative index vertices v ∈ V<tspan baseline-shift="sub" font-size="75%">indices</tspan>.</text>
        <text x="10" y="138" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Executed via `milopy.core.make_nhoods(adata, prop=0.1, refined=True)`.</text>
      </g>
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 3: Stage 2 — Neighborhood Definition, Topology & Composition
    # Column 2: x=372, y=78, width=316, height=586
    # =========================================================================
    svg_parts.append(f"""
  <!-- LAYER 03: STAGE 2 — NEIGHBORHOODS & COMPOSITION -->
  <g inkscape:groupmode="layer" id="layer-03-stage2" inkscape:label="03_Stage2_Neighborhoods">
    <g id="col-stage2" transform="translate(372, 78)">
      <!-- Wireframe Outer Frame -->
      <rect width="316" height="586" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        2. NEIGHBORHOODS &amp; TOPOLOGY
      </text>
      <line x1="14" y1="30" x2="302" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 2A: Overlapping Neighborhood Ensemble -->
      <g id="panel-nhood-def" transform="translate(14, 38)">
        <rect width="288" height="152" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Overlapping Neighborhood Ensemble:
        </text>

        <!-- Diagram: 3 Overlapping Neighborhoods Wireframe -->
        <g transform="translate(12, 26)">
          <rect width="82" height="74" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <!-- Nhood 1 (Blue) -->
          <circle cx="32" cy="30" r="22" fill="none" stroke="{OKABE_BLUE}" stroke-width="1.2" stroke-dasharray="3,2" />
          <circle cx="32" cy="30" r="3.5" fill="{OKABE_BLUE}" />
          <text x="32" y="23" font-size="6.5" font-weight="700" fill="{OKABE_BLUE}" text-anchor="middle">v<tspan baseline-shift="sub" font-size="75%">1</tspan></text>

          <!-- Nhood 2 (Orange) -->
          <circle cx="56" cy="26" r="20" fill="none" stroke="{OKABE_ORANGE}" stroke-width="1.2" stroke-dasharray="3,2" />
          <circle cx="56" cy="26" r="3.5" fill="{OKABE_ORANGE}" />
          <text x="56" y="19" font-size="6.5" font-weight="700" fill="{OKABE_ORANGE}" text-anchor="middle">v<tspan baseline-shift="sub" font-size="75%">2</tspan></text>

          <!-- Nhood 3 (Purple) -->
          <circle cx="44" cy="50" r="18" fill="none" stroke="{OKABE_REDDISH_PURPLE}" stroke-width="1.2" stroke-dasharray="3,2" />
          <circle cx="44" cy="50" r="3.5" fill="{OKABE_REDDISH_PURPLE}" />
          <text x="44" y="63" font-size="6.5" font-weight="700" fill="{OKABE_REDDISH_PURPLE}" text-anchor="middle">v<tspan baseline-shift="sub" font-size="75%">3</tspan></text>

          <!-- Shared intersection cells -->
          <circle cx="44" cy="32" r="3" fill="{OKABE_VERMILION}" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.5" />
          <text x="41" y="80" font-size="5.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Shared Boundary Cells</text>
        </g>

        <!-- Specs -->
        <g transform="translate(102, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Neighborhood Formulation:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• N(v) = {{v}} ∪ Neighbors(v) in G</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Index vertices: V nhoods across C cells</text>
          <rect x="0" y="38" width="174" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="49" font-size="6.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Overlapping Property:</text>
          <text x="6" y="58" font-size="6" fill="{COLOR_TEXT_SECONDARY}">Cells belong to multiple neighborhoods</text>
        </g>
        <text x="10" y="142" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Enables continuous statistical resolution across transition states.</text>
      </g>

      <!-- Panel 2B: Sparse Binary Incidence Matrix -->
      <g id="panel-incidence-matrix" transform="translate(14, 198)">
        <rect width="288" height="175" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Sparse Binary Incidence Matrix (adata.obsm['nhoods']):
        </text>

        <!-- Matrix graphic -->
        <g transform="translate(14, 26)">
          <rect width="104" height="74" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <line x1="26" y1="0" x2="26" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="52" y1="0" x2="52" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <line x1="78" y1="0" x2="78" y2="74" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <!-- 1s pattern -->
          <rect x="3" y="6" width="20" height="14" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="29" y="6" width="20" height="14" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="29" y="27" width="20" height="14" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="55" y="27" width="20" height="14" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="55" y="50" width="20" height="18" fill="{OKABE_BLUE}" opacity="0.85" />
          <rect x="81" y="50" width="20" height="18" fill="{OKABE_BLUE}" opacity="0.85" />

          <text x="52" y="86" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">Nhoods (V) →</text>
          <text x="-37" y="-4" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90)">Cells (C) →</text>
        </g>

        <g transform="translate(130, 26)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Matrix Properties:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Binary sparse CSC matrix</text>
          <text x="0" y="31" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Rows = Single cells (C)</text>
          <text x="0" y="41" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Columns = Nhoods (V)</text>
          <rect x="0" y="48" width="144" height="26" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="59" font-size="6.5" font-weight="700" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">M<tspan baseline-shift="sub" font-size="75%">c, v</tspan> = 1 if c ∈ N(v)</text>
          <text x="6" y="68" font-size="6" fill="{COLOR_TEXT_SECONDARY}">Ultra-fast sparse arithmetic</text>
        </g>
        <text x="10" y="130" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• <tspan font-weight="700">Vectorized Cell Scoring:</tspan> c<tspan baseline-shift="sub" font-size="75%">resp</tspan> = A · v<tspan baseline-shift="sub" font-size="75%">resp</tspan> in &lt; 5 ms.</text>
        <text x="10" y="142" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Avoids slow iterative neighborhood unwrapping loops.</text>
        <text x="10" y="156" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Enables instant single-cell DA mapping across 50k+ cells.</text>
      </g>

      <!-- Panel 2C: Neighborhood Distance & Composition -->
      <g id="panel-nhood-distance" transform="translate(14, 381)">
        <rect width="288" height="191" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Neighborhood Distance &amp; Composition:
        </text>

        <!-- Display Equation (1) -->
        <g id="eq-nhood-dist" transform="translate(10, 24)">
          <rect width="268" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="14" y="21" font-size="9" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            D<tspan baseline-shift="sub" font-size="75%">u, v</tspan> = || median(N<tspan baseline-shift="sub" font-size="75%">u</tspan>) - median(N<tspan baseline-shift="sub" font-size="75%">v</tspan>) ||<tspan baseline-shift="sub" font-size="75%">2</tspan>
          </text>
          <text x="248" y="21" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(1)</text>
        </g>

        <!-- Cross-Cohort Composition Formulation Box -->
        <g transform="translate(10, 64)">
          <rect width="268" height="42" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cross-Cohort Neighborhood Composition:</text>
          <text x="8" y="27" font-size="9" font-weight="700" font-family="{FONT_SERIF_MATH}" fill="{OKABE_BLUE}">
            C = B<tspan baseline-shift="super" font-size="75%">T</tspan> · A  ∈  ℕ₀<tspan baseline-shift="super" font-size="75%">K × V</tspan>
          </text>
          <text x="8" y="38" font-size="6.2" fill="{COLOR_TEXT_SECONDARY}">B = one-hot dataset matrix; A = incidence matrix</text>
        </g>

        <g transform="translate(10, 114)">
          <text x="0" y="8" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Purity Metric: Π<tspan baseline-shift="sub" font-size="75%">v</tspan> = max<tspan baseline-shift="sub" font-size="75%">k</tspan> C<tspan baseline-shift="sub" font-size="75%">k, v</tspan> / ∑<tspan baseline-shift="sub" font-size="75%">k</tspan> C<tspan baseline-shift="sub" font-size="75%">k, v</tspan></text>
          <text x="0" y="20" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Normalized Shannon Entropy: H<tspan baseline-shift="sub" font-size="75%">v</tspan> ∈ [0, 1]</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Distinguishes true conserved responses from batch artefacts</text>
          <text x="0" y="46" font-size="6.8" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Guarantees neighborhood cross-dataset biological validity</text>
        </g>
      </g>
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 4: Stage 3 — Patient Cell Counting Matrix & Error Control
    # Column 3: x=708, y=78, width=316, height=586
    # =========================================================================
    svg_parts.append(f"""
  <!-- LAYER 04: STAGE 3 — PATIENT CELL COUNTING -->
  <g inkscape:groupmode="layer" id="layer-04-stage3" inkscape:label="04_Stage3_Cell_Counting">
    <g id="col-stage3" transform="translate(708, 78)">
      <!-- Wireframe Outer Frame -->
      <rect width="316" height="586" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        3. PATIENT CELL COUNTING MATRIX
      </text>
      <line x1="14" y1="30" x2="302" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 3A: Counting Formulation -->
      <g id="panel-counting-def" transform="translate(14, 38)">
        <rect width="288" height="96" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Count Cells per Neighborhood &amp; Patient Sample:
        </text>

        <!-- Display Equation (2) -->
        <g id="eq-counting" transform="translate(10, 24)">
          <rect width="268" height="38" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="20" y="24" font-size="10" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            N<tspan baseline-shift="sub" font-size="75%">v, j</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">c ∈ sample j</tspan>  𝕀( c ∈ N(v) )
          </text>
          <text x="248" y="24" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(2)</text>
        </g>
        <text x="10" y="78" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Produces count table: V neighborhoods × J patient biological samples</text>
        <text x="10" y="88" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Maintains true biological replication structure (no pseudo-replication)</text>
      </g>

      <!-- Panel 3B: Neighborhood Counts Table Display -->
      <g id="panel-counts-table" transform="translate(14, 142)">
        <rect width="288" height="200" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Neighborhood Count Table (V nhoods × J samples):
        </text>

        <!-- Table schematic -->
        <g transform="translate(10, 26)">
          <rect width="268" height="92" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          
          <!-- Column headers: Patients -->
          <text x="75" y="14" font-size="7" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="middle">Responders (Y=1)</text>
          <text x="195" y="14" font-size="7" font-weight="700" fill="{OKABE_VERMILION}" text-anchor="middle">Non-Responders (Y=0)</text>
          <line x1="134" y1="4" x2="134" y2="88" stroke="{COLOR_DIVIDER_RULE}" stroke-width="1" stroke-dasharray="2,2" />

          <!-- Table rows -->
          <g transform="translate(8, 22)">
            <!-- Header labels -->
            <text x="0" y="10" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Nhood v1:</text>
            <text x="0" y="25" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Nhood v2:</text>
            <text x="0" y="40" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Nhood v3:</text>
            <text x="0" y="55" font-size="6.5" font-weight="600" fill="{COLOR_TEXT_MUTED}">Nhood v4:</text>

            <!-- Numbers (Responder enriched in v1, NR enriched in v3) -->
            <text x="50" y="10" font-size="6.8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">42</text>
            <text x="80" y="10" font-size="6.8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">38</text>
            <text x="110" y="10" font-size="6.8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">51</text>
            <text x="165" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">2</text>
            <text x="195" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">0</text>
            <text x="225" y="10" font-size="6.8" fill="{COLOR_TEXT_MUTED}">4</text>

            <text x="50" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">12</text>
            <text x="80" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">15</text>
            <text x="110" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">10</text>
            <text x="165" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">14</text>
            <text x="195" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">11</text>
            <text x="225" y="25" font-size="6.8" fill="{COLOR_TEXT_MUTED}">13</text>

            <text x="50" y="40" font-size="6.8" fill="{COLOR_TEXT_MUTED}">1</text>
            <text x="80" y="40" font-size="6.8" fill="{COLOR_TEXT_MUTED}">3</text>
            <text x="110" y="40" font-size="6.8" fill="{COLOR_TEXT_MUTED}">0</text>
            <text x="165" y="40" font-size="6.8" font-weight="700" fill="{OKABE_VERMILION}">64</text>
            <text x="195" y="40" font-size="6.8" font-weight="700" fill="{OKABE_VERMILION}">72</text>
            <text x="225" y="40" font-size="6.8" font-weight="700" fill="{OKABE_VERMILION}">58</text>

            <text x="50" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">25</text>
            <text x="80" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">19</text>
            <text x="110" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">22</text>
            <text x="165" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">21</text>
            <text x="195" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">28</text>
            <text x="225" y="55" font-size="6.8" fill="{COLOR_TEXT_MUTED}">24</text>
          </g>
        </g>

        <!-- Library Size Normalization Box -->
        <g transform="translate(10, 128)">
          <rect width="268" height="44" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="13" font-size="7" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Library Size Offset &amp; Normalization:</text>
          <text x="6" y="25" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• Sample cellularity: TotalCells<tspan baseline-shift="sub" font-size="75%">j</tspan> = ∑<tspan baseline-shift="sub" font-size="75%">c</tspan> 𝕀( c ∈ sample j )</text>
          <text x="6" y="36" font-size="6.5" fill="{COLOR_TEXT_SECONDARY}">• TMM (Trimmed Mean of M-values) calculates effective library size</text>
        </g>
        <text x="10" y="188" font-size="6.8" fill="{COLOR_TEXT_MUTED}">Stored in `adata.uns['nhood_counts']` ready for GLM testing.</text>
      </g>

      <!-- Panel 3C: Eliminating Pseudo-Replication -->
      <g id="panel-why-replicates" transform="translate(14, 350)">
        <rect width="288" height="222" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Eliminating Single-Cell Pseudo-Replication:
        </text>

        <!-- Visual Schematic of Error Inflation -->
        <g transform="translate(10, 24)">
          <rect width="268" height="74" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{OKABE_VERMILION}">The Pseudo-Replication Fallacy (Direct Cell Tests):</text>
          <text x="8" y="25" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Pooling 50,000 cells directly treats cells as independent observations.</text>
          <text x="8" y="37" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Severe variance underestimation: p-values artificially crash to p &lt; 10⁻⁵⁰.</text>
          <text x="8" y="49" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Single dominant patient completely biases the entire result.</text>
          <rect x="8" y="55" width="252" height="14" fill="#FEE2E2" stroke="#EF4444" stroke-width="0.5" />
          <text x="134" y="65" font-size="6" font-weight="700" fill="#991B1B" text-anchor="middle">Catastrophic Type-I Error Explosion: False discoveries dominate</text>
        </g>

        <!-- Milo Statistical Solution -->
        <g transform="translate(10, 106)">
          <rect width="268" height="88" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="8" y="13" font-size="7.5" font-weight="700" fill="{OKABE_BLUE}">The Milo Hierarchical Solution (Patient GLM):</text>
          <text x="8" y="25" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• True biological units: Evaluates variation across J patient biopsies.</text>
          <text x="8" y="37" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Models patient overdispersion with Negative Binomial parameter φ<tspan baseline-shift="sub" font-size="75%">v</tspan>.</text>
          <text x="8" y="49" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Accommodates clinical covariates &amp; batch: ~ condition + batch.</text>
          <text x="8" y="61" font-size="6.8" font-weight="600" fill="{OKABE_BLUISH_GREEN}">• Calibrated error control: Robust against extreme sample imbalance.</text>
          <text x="8" y="75" font-size="6.2" fill="{COLOR_TEXT_MUTED}">Solves the PDAC 15 R vs 2 NR dilemma via empirical shrinkage.</text>
        </g>
      </g>
    </g>
  </g>
""")

    # =========================================================================
    # LAYER 5: Stage 4 — Negative Binomial GLM & Spatial DA
    # Column 4: x=1044, y=78, width=320, height=586
    # =========================================================================
    svg_parts.append(f"""
  <!-- LAYER 05: STAGE 4 — GLM & SPATIAL DA -->
  <g inkscape:groupmode="layer" id="layer-05-stage4" inkscape:label="05_Stage4_GLM_and_DA">
    <g id="col-stage4" transform="translate(1044, 78)">
      <!-- Wireframe Outer Frame -->
      <rect width="320" height="586" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      
      <!-- Stage Section Header -->
      <text x="14" y="22" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="0.3">
        4. NEGATIVE BINOMIAL GLM &amp; SPATIAL DA
      </text>
      <line x1="14" y1="30" x2="306" y2="30" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />

      <!-- Panel 4A: Negative Binomial QL-GLM -->
      <g id="panel-glm-model" transform="translate(14, 38)">
        <rect width="292" height="130" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Negative Binomial QL-GLM:
        </text>

        <!-- Display Equation (3) -->
        <g id="eq-glm" transform="translate(10, 24)">
          <rect width="272" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="6" y="21" font-size="8.2" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            log( 𝔼[N<tspan baseline-shift="sub" font-size="75%">v, j</tspan>] ) = β<tspan baseline-shift="sub" font-size="75%">0, v</tspan> + β<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">milo</tspan> · Y<tspan baseline-shift="sub" font-size="75%">j</tspan> + log( TotalCells<tspan baseline-shift="sub" font-size="75%">j</tspan> )
          </text>
          <text x="250" y="21" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(3)</text>
        </g>
        <text x="10" y="72" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Log-link with library size offset; Quasi-Likelihood F-test</text>
        <text x="10" y="84" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Y<tspan baseline-shift="sub" font-size="75%">j</tspan> ∈ {{0, 1}} : Binary clinical response (Res vs NR)</text>
        <text x="10" y="96" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• β<tspan baseline-shift="sub" font-size="75%">v</tspan><tspan baseline-shift="super" font-size="75%">milo</tspan> : Differential abundance effect size (logFC = β / ln 2)</text>
        <text x="10" y="110" font-size="6.8" font-weight="600" fill="{OKABE_BLUE}">Empirical Bayes dispersion shrinkage squeezes variance to trend</text>
      </g>

      <!-- Panel 4B: Spatial FDR Weighted Benjamini-Hochberg -->
      <g id="panel-spatial-fdr" transform="translate(14, 176)">
        <rect width="292" height="130" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Spatial FDR (Weighted Benjamini-Hochberg):
        </text>

        <!-- Display Equation (4) -->
        <g id="eq-spatial-fdr" transform="translate(10, 24)">
          <rect width="272" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
          <text x="24" y="21" font-size="9.5" font-weight="600" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_PRIMARY}">
            weight(v) ∝ 1 / connectivity_density(v)
          </text>
          <text x="250" y="21" font-size="9.5" font-family="{FONT_SERIF_MATH}" fill="{COLOR_TEXT_MUTED}">(4)</text>
        </g>
        <text x="10" y="72" font-size="6.8" fill="{COLOR_TEXT_MUTED}">• Weight is inversely proportional to kNN overlap degree</text>
        <text x="10" y="84" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Overlapping neighborhoods induce correlated test statistics</text>
        <text x="10" y="96" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Prevents dense manifold regions from dominating error rate</text>
        <text x="10" y="110" font-size="6.8" font-weight="600" fill="{OKABE_BLUISH_GREEN}">Rigorous false discovery control across graph tests (FDR &lt; 0.10)</text>
      </g>

      <!-- Panel 4C: Dual Outputs & Biological Discovery -->
      <g id="panel-impact" transform="translate(14, 314)">
        <rect width="292" height="258" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="10" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          Dual Analytic Outputs &amp; Biomarkers:
        </text>

        <!-- Mini Visual: Volcano Plot & UMAP logFC side by side -->
        <g transform="translate(10, 24)">
          <!-- Mini Volcano -->
          <g transform="translate(0, 0)">
            <rect width="130" height="82" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="65" y="11" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Cohort Volcano Grid</text>
            <!-- Axes -->
            <line x1="16" y1="68" x2="120" y2="68" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
            <line x1="68" y1="16" x2="68" y2="68" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" stroke-dasharray="2,2" />
            <line x1="16" y1="16" x2="16" y2="68" stroke="{COLOR_TEXT_PRIMARY}" stroke-width="0.75" />
            <!-- Significance threshold line -->
            <line x1="16" y1="44" x2="120" y2="44" stroke="#EF4444" stroke-width="0.75" stroke-dasharray="2,2" />
            <text x="96" y="41" font-size="5" font-weight="600" fill="#EF4444">FDR = 0.10</text>
            <!-- Significant Dots (Responders > 0) -->
            <circle cx="95" cy="25" r="2.5" fill="{OKABE_BLUISH_GREEN}" />
            <circle cx="104" cy="32" r="2.5" fill="{OKABE_BLUISH_GREEN}" />
            <circle cx="88" cy="36" r="2" fill="{OKABE_BLUISH_GREEN}" />
            <circle cx="112" cy="22" r="2" fill="{OKABE_BLUISH_GREEN}" />
            <circle cx="82" cy="40" r="1.8" fill="{OKABE_BLUISH_GREEN}" />
            <!-- Significant Dots (Non-Responders < 0) -->
            <circle cx="36" cy="28" r="2.5" fill="{OKABE_VERMILION}" />
            <circle cx="44" cy="34" r="2" fill="{OKABE_VERMILION}" />
            <circle cx="28" cy="24" r="2" fill="{OKABE_VERMILION}" />
            <circle cx="50" cy="41" r="1.8" fill="{OKABE_VERMILION}" />
            <!-- Non-significant points -->
            <circle cx="62" cy="54" r="1.5" fill="{COLOR_TEXT_HAIRLINE}" />
            <circle cx="74" cy="56" r="1.5" fill="{COLOR_TEXT_HAIRLINE}" />
            <circle cx="58" cy="60" r="1.5" fill="{COLOR_TEXT_HAIRLINE}" />
            <circle cx="78" cy="62" r="1.5" fill="{COLOR_TEXT_HAIRLINE}" />
            <circle cx="68" cy="58" r="1.5" fill="{COLOR_TEXT_HAIRLINE}" />
            <text x="68" y="77" font-size="5.5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">log2FC</text>
            <text x="12" y="42" font-size="5" fill="{COLOR_TEXT_MUTED}" text-anchor="middle" transform="rotate(-90 12 42)">-log10 FDR</text>
          </g>

          <!-- Mini UMAP logFC -->
          <g transform="translate(138, 0)">
            <rect width="130" height="82" fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <text x="65" y="11" font-size="6" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">UMAP log2FC Gradient</text>
            <!-- Cloud of points -->
            <circle cx="45" cy="38" r="16" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <circle cx="85" cy="44" r="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
            <!-- Intermediate background cells -->
            <circle cx="44" cy="48" r="2" fill="{COLOR_TEXT_HAIRLINE}" opacity="0.6" />
            <circle cx="65" cy="40" r="2" fill="{COLOR_TEXT_HAIRLINE}" opacity="0.6" />
            <circle cx="70" cy="48" r="2" fill="{COLOR_TEXT_HAIRLINE}" opacity="0.6" />
            <!-- Enriched non-responder zone -->
            <circle cx="38" cy="32" r="4.5" fill="{OKABE_VERMILION}" opacity="0.85" />
            <circle cx="48" cy="36" r="3.5" fill="{OKABE_VERMILION}" opacity="0.85" />
            <circle cx="32" cy="40" r="3" fill="{OKABE_VERMILION}" opacity="0.85" />
            <!-- Enriched responder zone -->
            <circle cx="84" cy="40" r="4.5" fill="{OKABE_BLUISH_GREEN}" opacity="0.85" />
            <circle cx="92" cy="46" r="5" fill="{OKABE_BLUISH_GREEN}" opacity="0.85" />
            <circle cx="78" cy="50" r="3.5" fill="{OKABE_BLUISH_GREEN}" opacity="0.85" />
            <!-- Clamped colorbar -->
            <g transform="translate(15, 68)">
              <rect width="100" height="4" fill="{COLOR_BORDER_HAIRLINE}" />
              <rect x="0" y="0" width="35" height="4" fill="{OKABE_VERMILION}" />
              <rect x="35" y="0" width="30" height="4" fill="#CBD5E1" />
              <rect x="65" y="0" width="35" height="4" fill="{OKABE_BLUISH_GREEN}" />
              <text x="0" y="10" font-size="5" font-weight="600" fill="{COLOR_TEXT_MUTED}">-5</text>
              <text x="50" y="10" font-size="5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="middle">0</text>
              <text x="100" y="10" font-size="5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="end">+5</text>
            </g>
          </g>
        </g>

        <!-- Biological Findings & Concordance Bullets -->
        <g transform="translate(10, 118)">
          <text x="0" y="10" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">1. Validated Immunotherapy Biomarkers:</text>
          <text x="0" y="21" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• <tspan font-weight="600" fill="{OKABE_BLUISH_GREEN}">Naive B cells</tspan> enriched in Responders (β = +2.07, p = 3e-4)</text>
          <text x="0" y="32" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• <tspan font-weight="600" fill="{OKABE_VERMILION}">M2 Macrophages &amp; Tex</tspan> enriched in Non-Responders (β = -1.56)</text>
          
          <text x="0" y="48" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">2. Milo DA vs. Pseudobulk Deconvolution:</text>
          <text x="0" y="59" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Concordant in metastatic melanoma (Spearman ρ = 0.853)</text>
          <text x="0" y="70" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Disentangles physical cell counts from mRNA mass fractions</text>
          <text x="0" y="82" font-size="6.8" fill="{COLOR_TEXT_SECONDARY}">• Solves state-activation confounding in cytokine-induced shifts</text>
          <text x="0" y="98" font-size="6.8" font-weight="700" fill="{OKABE_BLUE}">Gold-standard continuous single-cell ground truth</text>
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
    <path d="M 352 230 L 370 230" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 352 420 L 370 420" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 2 -> 3 -->
    <path d="M 688 230 L 706 230" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 688 420 L 706 420" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />

    <!-- Connector 3 -> 4 -->
    <path d="M 1024 230 L 1042 230" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
    <path d="M 1024 420 L 1042 420" stroke="{COLOR_TEXT_SECONDARY}" stroke-width="1" fill="none" marker-end="url(#arr-slate)" />
  </g>
</svg>
""")

    return "".join(svg_parts)


def main() -> None:
    """Generate and write the Nature Methods milopy SVG and 300 DPI PNG."""
    print("=" * 70)
    print("GENERATING NATURE METHODS MILOPY METHODOLOGY FIGURE")
    print("=" * 70)

    svg_content = render_nature_methods_milopy()

    # Destination directories
    deconv_dir = Path("article/figures/deconvolution")
    reports_dir = Path("output/reports")
    deconv_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    svg_path_article = deconv_dir / "how_milopy_works.svg"
    png_path_article = deconv_dir / "how_milopy_works.png"
    svg_path_reports = reports_dir / "how_milopy_works.svg"
    png_path_reports = reports_dir / "how_milopy_works.png"

    # Write SVG files
    svg_path_article.write_text(svg_content, encoding="utf-8")
    svg_path_reports.write_text(svg_content, encoding="utf-8")
    print(f"[EXPORTED] SVG -> {svg_path_article} ({svg_path_article.stat().st_size / 1e3:.1f} KB)")
    print(f"[EXPORTED] SVG -> {svg_path_reports} ({svg_path_reports.stat().st_size / 1e3:.1f} KB)")

    # Validate XML well-formedness
    try:
        tree = ET.parse(svg_path_article)
        root = tree.getroot()
        layers = [
            elem for elem in root.iter()
            if elem.attrib.get("{http://www.inkscape.org/namespaces/inkscape}groupmode") == "layer"
        ]
        print(f"[VERIFIED] SVG is 100% valid XML ({root.tag}) with width={root.attrib.get('width')}, height={root.attrib.get('height')}")
        print(f"[VERIFIED] Inkscape native layers detected: {len(layers)}")
    except Exception as exc:
        print(f"[FATAL XML ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    # Render 300 DPI high-resolution PNG rasters via vl-convert or rsvg-convert
    converted = False
    if vlc is not None:
        try:
            png_bytes = vlc.svg_to_png(svg_content, scale=3.125)  # 300 DPI (300 / 96 = 3.125)
            png_path_article.write_bytes(png_bytes)
            png_path_reports.write_bytes(png_bytes)
            print(f"[EXPORTED] 300 DPI PNG via vl-convert -> {png_path_article} ({len(png_bytes) / 1e3:.1f} KB)")
            print(f"[EXPORTED] 300 DPI PNG via vl-convert -> {png_path_reports} ({len(png_bytes) / 1e3:.1f} KB)")
            converted = True
        except Exception as exc:
            print(f"[vl-convert warning] {exc}")

    if not converted:
        try:
            subprocess.run(
                ["rsvg-convert", "-d", "300", "-p", "300", str(svg_path_article), "-o", str(png_path_article)],
                check=True,
            )
            subprocess.run(
                ["rsvg-convert", "-d", "300", "-p", "300", str(svg_path_reports), "-o", str(png_path_reports)],
                check=True,
            )
            print(f"[EXPORTED] 300 DPI PNG via rsvg-convert -> {png_path_article} ({png_path_article.stat().st_size / 1e3:.1f} KB)")
            print(f"[EXPORTED] 300 DPI PNG via rsvg-convert -> {png_path_reports} ({png_path_reports.stat().st_size / 1e3:.1f} KB)")
            converted = True
        except Exception as exc:
            print(f"[PNG CONVERSION WARNING] {exc}", file=sys.stderr)

    print("\n" + "=" * 70)
    print("MILOPY METHOD FIGURE GENERATION COMPLETED SUCCESSFULLY.")
    print("=" * 70)


if __name__ == "__main__":
    main()
