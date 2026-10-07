---
name: publication-vector-figures
description: >-
  Actionable guide, architecture patterns, and code templates for generating publication-grade,
  top-tier academic journal vector figures (SVG) with native Inkscape layers, clean semantic grouping,
  and 300 DPI verification via rsvg-convert.
---

# Publication Vector Figures & Inkscape Standards Guide

This skill provides step-by-step procedures and code patterns for constructing publication-grade vector graphics (`SVG`) suitable for top-tier journals (*Nature*, *Cell*, *Science*, *Nature Medicine*) and fully editable in Inkscape.

---

## 1. Figure Scope Classification: Data Compendium vs. Study Overview

Always determine the exact conceptual scope before laying out elements:

| Scope Type | Purpose | Included Elements | Strictly Excluded |
| :--- | :--- | :--- | :--- |
| **Pure Dataset Compendium ("Data-Only")** | Illustrate exclusively the data assets assembled across the study. | Patient cohorts, sample sizes ($N$), RECIST response rates, clinical stratifications, multi-omic assay completeness matrix, cellular compartment breakdown, genomic burden / TMB / VAF ranges. | **Zero** analytical algorithms, mathematical formulas, deconvolution equations, simulation sweeps, ROC curves, or biomarker outcome predictions. |
| **Study Design / Graphical Abstract** | End-to-end translation from specimens to biological discovery. | Specimen ecosystem, wet-lab / sequencing workflows, computational pipeline streams, mathematical equations, benchmark results, and biomarker discovery plots. | Redundant multi-page assay tables (keep to high-level synthesis). |

---

## 2. Inkscape XML Namespace & Layer Architecture

Inkscape relies on standard SVG with special namespaced attributes. Without these attributes, Inkscape treats all elements as a single flat canvas.

### Required XML Header & Namespaces
```xml
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="1400" height="900" viewBox="0 0 1400 900"
     version="1.1" id="figure-id">
```

### Native Layer Structure
Declare top-level `<g>` elements with `inkscape:groupmode="layer"` and human-readable `inkscape:label`:

```xml
<!-- 00 Canvas Background -->
<g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
  <rect id="bg-canvas" ... />
  <rect id="card-panel-a" ... />
</g>

<!-- 01 Figure Header -->
<g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
  <text id="txt-title" ...>Figure Title</text>
  <g id="badge-panel-a"> ... </g>
</g>

<!-- 02 Panel A -->
<g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Title">
  <g id="grp-cohort-cards" inkscape:label="Cohort Cards"> ... </g>
</g>

<!-- 03 Inter-Panel Connectors & Legends -->
<g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Flow_Connectors">
  <path id="arrow-a-to-b" ... />
</g>
```

---

## 3. Formatting & Design Invariants

1. **Editable Typography**:
   - Use native `<text>` tags with font fallbacks: `font-family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"`.
   - Never convert text to outlines/paths unless custom typography is strictly required.
   - **Avoid combining unicode marks**: Do not write combining circumflexes like `β̂` (`\u0302`), as cairo/rsvg-convert will render them with offset or missing glyphs. Write standalone unicode symbols (`β`, `Δ`, `ρ`, `OR`, `p-value`).

2. **Preventing Coordinate Clashes in Grid Cards**:
   - When placing multiple sub-cards in a horizontal track, wrap each sub-card in a `<g transform="translate(X, Y)">` and write local coordinates relative to `(0, 0)`:
   ```xml
   <!-- Card 1 -->
   <g transform="translate(0, 0)">
     <rect width="180" height="42" ... />
     <text x="8" y="14">Card 1 Title</text>
   </g>
   <!-- Card 2 -->
   <g transform="translate(190, 0)">
     <rect width="180" height="42" ... />
     <text x="8" y="14">Card 2 Title</text>
   </g>
   ```

3. **Gutter-Constrained Connectors**:
   - Connectors and flow arrows must stay strictly within gutters (margins between panel cards).
   - Never draw a connector that cuts across an unrelated container card, table, or scatter plot.

4. **Nature Methods Minimalist Wireframe & Okabe-Ito Scientific Palette**:
   - **Canvas Background**: `#FFFFFF` (Pure white).
   - **Container Cards**: `#FFFFFF` with hairline stroke `#CBD5E1` (`0.5 pt` to `0.75 pt`), `rx=0` (sharp square corners).
   - **Zero Drop Shadows**: `feDropShadow` and filter glows are strictly prohibited.
   - **Colorblind-Safe Palette (Okabe-Ito)**: Applied strictly to data marks, curves, DAG nodes, and scatter points—**never** as decorative card background fills:
     - Slate / Dark: `#1E293B` (Structural text, primary labels)
     - Sky Blue: `#56B4E9` (Bladder, CD4+, secondary modalities)
     - Bluish Green: `#009E73` (Renal, Responders, B-cells)
     - Orange: `#E69F00` (Melanoma, Myeloid, discordant highlights)
     - Blue: `#0072B2` (CD8+, primary deconvolution, correlations)
     - Vermilion: `#D55E00` (Malignant, Non-responders, key alterations)
     - Reddish Purple: `#CC79A7` (Breast, PDAC, secondary markers)

5. **Formal Mathematical Typesetting**:
   - Math font: `"Times New Roman", Times, Georgia, serif`.
   - Center display equations; use italic styling for scalar variables ($x, y, \beta$) and upright Roman for operators/functions ($\log, \Pr, \mathbb{E}$).
   - Number formal equations with right-aligned tags `(1)`, `(2)`, etc.
   - Use native `<tspan baseline-shift="sub">` and `<tspan baseline-shift="super">` for indices, and standalone clean unicode symbols (`β`, `Δ`, `ρ`, `θ̂`, `∑`, `Φ`).
   - Draw horizontal division rules with `<line stroke="#1E293B" stroke-width="0.75" />` instead of diagonal slashes in formal display fractions.

6. **Probabilistic Plate Notation & Coordinate Axes**:
   - **Bayesian Plate Notation**: Double circles (`<circle ... stroke-width="1.5" /><circle ... stroke-width="0.75" />`) for observed data ($Y$), single circles for latent variables/parameters ($\theta, Z$), and labeled bounding rectangles for plates with indexing labels ($g \in \{1,\dots,G\}$, $n \in \{1,\dots,N\}$).
   - **Coordinate Axes**: Line plots, trace plots, and scatter plots must include real continuous axes with outward tick marks (`stroke="#1E293B"`, `stroke-width="0.75"`), explicit axis labels with units, and legible numeric tick labels.

7. **Shared Design Tokens (`nature_style_config.py`)**:
   - All figure scripts should import shared styling tokens from `scripts/figures/nature_style_config.py` to ensure visual uniformity across the manuscript.

---

## 4. Procedural Generation & Verification Workflow

Always write a declarative, pure Python script to render the SVG, then verify:

```bash
# 1. Generate the SVG
/usr/bin/python3 scripts/figures/generate_figure.py

# 2. Validate XML syntax and Inkscape layers
/usr/bin/python3 -c "
import xml.etree.ElementTree as ET
tree = ET.parse('output.svg')
layers = [e for e in tree.getroot().iter() if e.attrib.get('{http://www.inkscape.org/namespaces/inkscape}groupmode') == 'layer']
print(f'Valid XML with {len(layers)} Inkscape layers')
"

# 3. Render 300 DPI PNG with rsvg-convert (2x scale for 1400px = 2800px width)
/opt/homebrew/bin/rsvg-convert -w 2800 -h 1800 output.svg -o output.png

# 4. Compile document
/opt/homebrew/bin/typst compile article/main.typ article/main.pdf
```
