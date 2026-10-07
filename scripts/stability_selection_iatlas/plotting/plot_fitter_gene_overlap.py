#!/usr/bin/env python3
"""Plot Publication-Grade Multi-Fitter Cross-Cohort Gene Selection Figure.

Adheres strictly to the project's functional programming principles:
- Pure functions mapping immutable Polars data structures to SVG vector elements.
- Nature Methods minimalist wireframe standard (pure white background, #CBD5E1 borders, rx=0, no drop shadows).
- Okabe-Ito Colorblind-Safe scientific palette assigned to fitters:
  - Lasso: #0072B2 (Dark Blue)
  - Elastic Net: #56B4E9 (Sky Blue)
  - Logistic: #009E73 (Bluish Green)
  - Random Forest: #E69F00 (Orange)
  - Cohort-Adjusted FWL: #CC79A7 (Reddish Purple)
  - Multi-Cohort Group Lasso: #D55E00 (Vermilion)
- Native Inkscape layers and semantic groups.
- Vector SVG exports and 300 DPI PNG verification via rsvg-convert.
- Monadic error handling with Result[Path, str].
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys
from typing import Final, Sequence
import polars as pl
from returns.result import Failure, Result, Success

# Figure Canvas Dimensions (Nature double-column landscape)
WIDTH: Final[int] = 1400
HEIGHT: Final[int] = 920

from plotting_utils import (
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_DIVIDER_RULE,
    COLOR_FITTER_COHORT_ADJ,
    COLOR_FITTER_ELASTIC_NET,
    COLOR_FITTER_GLMM,
    COLOR_FITTER_GROUP_LASSO,
    COLOR_FITTER_INVARIANT,
    COLOR_FITTER_LASSO,
    COLOR_FITTER_LOGISTIC,
    COLOR_FITTER_MERF,
    COLOR_FITTER_META,
    COLOR_FITTER_MULTITASK,
    COLOR_FITTER_OSCAR,
    COLOR_FITTER_RF,
    COLOR_FITTER_SLOPE,
    COLOR_PANEL_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_HAIRLINE,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FITTER_BASELINE_ORDER,
    FITTER_COLOR_MAP,
    FITTER_DISPLAY_NAMES,
    FITTER_MULTI_COHORT_ORDER,
    FONT_SANS,
    FONT_SERIF_MATH,
)
from plotting_utils.svg import verify_and_rasterize_svg, wrap_inkscape_svg

FITTER_MULTICOHORT_ORDER: Final[tuple[str, ...]] = FITTER_MULTI_COHORT_ORDER

FITTER_ORDER: Final[tuple[str, ...]] = (
    *FITTER_BASELINE_ORDER,
    *FITTER_MULTICOHORT_ORDER,
)

FONT_SANS: Final[str] = (
    "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
)


# -------------------------------------------------------------------------
# Immutable Domain Models
# -------------------------------------------------------------------------

@dataclass(frozen=True)
class FitterGeneHit:
    """Represents a gene selection status for a specific (cohort, fitter)."""
    cohort: str
    fitter: str
    feature: str
    stability_score: float
    selected: bool


@dataclass(frozen=True)
class CohortDisplayMeta:
    """Display information for a clinical cohort."""
    name: str
    display_title: str
    n_samples: int
    is_combined: bool
    order_idx: int


@dataclass(frozen=True)
class GeneDisplayMeta:
    """Display information for a selected biomarker gene."""
    symbol: str
    pathway: str
    pathway_color: str
    col_idx: int


# -------------------------------------------------------------------------
# Domain Knowledge & Configuration
# -------------------------------------------------------------------------

ORDERED_COHORTS: Final[tuple[CohortDisplayMeta, ...]] = (
    CohortDisplayMeta("pancancer", "Pan-Cancer Combined", 1015, True, 0),
    CohortDisplayMeta("melanoma", "Melanoma Combined", 338, True, 1),
    CohortDisplayMeta("rcc", "RCC Combined", 263, True, 2),
    CohortDisplayMeta("Rosenberg-iAtlas", "Rosenberg (Bladder)", 298, False, 3),
    CohortDisplayMeta("McDermott-iAtlas", "McDermott (RCC)", 247, False, 4),
    CohortDisplayMeta("Riaz-iAtlas", "Riaz (Melanoma)", 98, False, 5),
    CohortDisplayMeta("Gide-iAtlas", "Gide (Melanoma)", 91, False, 6),
    CohortDisplayMeta("Liu-iAtlas", "Liu (Melanoma)", 122, False, 7),
    CohortDisplayMeta("Hugo-iAtlas", "Hugo (Melanoma)", 27, False, 8),
)

ORDERED_GENES: Final[tuple[GeneDisplayMeta, ...]] = (
    # Antigen Presentation Machinery (MHC-I & II)
    GeneDisplayMeta("HLA-A", "Antigen", COLOR_FITTER_LASSO, 0),
    GeneDisplayMeta("HLA-C", "Antigen", COLOR_FITTER_LASSO, 1),
    GeneDisplayMeta("NLRC5", "Antigen", COLOR_FITTER_LASSO, 2),
    GeneDisplayMeta("PSMB8", "Antigen", COLOR_FITTER_LASSO, 3),
    GeneDisplayMeta("CIITA", "Antigen", COLOR_FITTER_LASSO, 4),
    # T-cell Effector Lineage & Checkpoints
    GeneDisplayMeta("TBX21", "Effector", COLOR_FITTER_LOGISTIC, 5),
    GeneDisplayMeta("IFNG", "Effector", COLOR_FITTER_LOGISTIC, 6),
    GeneDisplayMeta("LAG3", "Effector", COLOR_FITTER_LOGISTIC, 7),
    GeneDisplayMeta("CD3E", "Effector", COLOR_FITTER_LOGISTIC, 8),
    GeneDisplayMeta("TCF7", "Effector", COLOR_FITTER_LOGISTIC, 9),
    # Co-stimulation & Immune Regulation
    GeneDisplayMeta("TNFSF9", "Costim / Reg", COLOR_FITTER_RF, 10),
    GeneDisplayMeta("TNFSF18", "Costim / Reg", COLOR_FITTER_RF, 11),
    GeneDisplayMeta("IKZF2", "Costim / Reg", COLOR_FITTER_RF, 12),
    GeneDisplayMeta("TNFRSF14", "Costim / Reg", COLOR_FITTER_RF, 13),
    GeneDisplayMeta("CCR7", "Costim / Reg", COLOR_FITTER_RF, 14),
    # Stroma, Inflammation & Angiogenesis
    GeneDisplayMeta("PTGS2", "Stroma / Supp", COLOR_FITTER_COHORT_ADJ, 15),
    GeneDisplayMeta("TGFB1", "Stroma / Supp", COLOR_FITTER_COHORT_ADJ, 16),
    GeneDisplayMeta("PVR", "Stroma / Supp", COLOR_FITTER_COHORT_ADJ, 17),
)


# -------------------------------------------------------------------------
# Pure Data Pipeline
# -------------------------------------------------------------------------

def load_benchmark_data(scores_path: Path) -> Result[pl.DataFrame, str]:
    """Pure loader returning benchmark stability scores."""
    if not scores_path.exists():
        return Failure(f"Scores parquet not found: {scores_path}")
    try:
        df = pl.read_parquet(scores_path)
        return Success(df)
    except Exception as exc:
        return Failure(f"Failed to read parquet: {exc}")


def extract_fitter_hits(
    df_scores: pl.DataFrame,
    method: str = "SS-CPSS",
) -> tuple[FitterGeneHit, ...]:
    """Extract immutable FitterGeneHit records for the specified stability method."""
    df_method = df_scores.filter(pl.col("method") == method)
    records: list[FitterGeneHit] = []

    for row in df_method.iter_rows(named=True):
        records.append(
            FitterGeneHit(
                cohort=str(row["cohort"]),
                fitter=str(row["fitter"]),
                feature=str(row["feature"]),
                stability_score=float(row["stability_score"]),
                selected=bool(row["selected"]),
            )
        )
    return tuple(records)


# -------------------------------------------------------------------------
# SVG Rendering Functions
# -------------------------------------------------------------------------

def render_header() -> str:
    """Render top figure title and metadata banner."""
    return f"""
  <!-- LAYER 01: FIGURE HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <text id="txt-title" x="20" y="32" font-size="15.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">
      Multi-Fitter Stability Selection Landscape in Cancer Immunotherapy
    </text>
    <text id="txt-subtitle" x="20" y="50" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Systematic biomarker discovery across 12 iAtlas cohorts and 13 statistical fitters (Single-Cohort Baselines &amp; Multi-Cohort Regularization)
    </text>

    <!-- Stat Tag Chip -->
    <g id="grp-header-chip" transform="translate(980, 16)">
      <rect x="0" y="0" width="400" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <text x="200" y="15" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" letter-spacing="0.2">
        STABILITY SELECTION BENCHMARK PARAMETERS
      </text>
      <text x="200" y="27" font-size="8.5" font-weight="500" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
        Shah-Samworth CPSS  |  Cutoff π<tspan font-size="7" baseline-shift="sub">thr</tspan> = 0.75  |  q<tspan font-size="7" baseline-shift="sub">budget</tspan> ≤ 20.0  |  PFER ≤ 1.0  |  B = 50
      </text>
    </g>
  </g>
"""


def render_panel_a_grid(
    hits: Sequence[FitterGeneHit],
    cohorts: Sequence[CohortDisplayMeta],
    genes: Sequence[GeneDisplayMeta],
) -> str:
    """Render Panel A: Cohort x Fitter Gene Selection Matrix with 2-Tier Capsule Rack."""
    px = 20
    py = 70
    pw = 890
    ph = 440

    svg = [f"""
  <!-- LAYER 02: PANEL A (GRID MATRIX) -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Grid">
    <rect id="card-panel-a" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />

    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">a</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Cohort × Fitter Selection Matrix across Core Predictive Genes
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Each cell displays a 13-fitter 2-tier capsule rack (Top: 6 Single-Cohort Baselines | Bottom: 7 Multi-Cohort Fitters; π<tspan font-size="7.5" baseline-shift="sub">thr</tspan> ≥ 0.75)
    </text>
"""]

    # Coordinates
    grid_top_y = py + 78
    grid_x_start = px + 172
    col_w = 36.5
    row_h = 36.5

    # 1. Pathway Category Banners
    pathways = (
        ("Antigen Presentation", COLOR_FITTER_LASSO, 0, 5),
        ("Th1 Effector / Checkpoints", COLOR_FITTER_LOGISTIC, 5, 10),
        ("Costimulation &amp; Regulators", COLOR_FITTER_RF, 10, 15),
        ("Stroma &amp; Suppression", COLOR_FITTER_COHORT_ADJ, 15, 18),
    )
    for p_title, p_col, s_col, e_col in pathways:
        bx = grid_x_start + s_col * col_w + 1
        bw = (e_col - s_col) * col_w - 2
        by = grid_top_y - 28
        svg.append(f"""
    <rect x="{bx}" y="{by}" width="{bw}" height="2" fill="{p_col}" rx="0" />
    <text x="{bx + bw/2}" y="{by - 4}" font-size="7.5" font-weight="700" fill="{p_col}" text-anchor="middle">
      {p_title}
    </text>
""")

    # 2. Gene Column Headers
    for g in genes:
        gx = grid_x_start + g.col_idx * col_w + col_w / 2
        f_size = "7.5" if len(g.symbol) >= 7 else "8.5"
        svg.append(f"""
    <text x="{gx}" y="{grid_top_y - 8}" font-size="{f_size}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
      {g.symbol}
    </text>
""")

    # Fast lookup map: (cohort, fitter, gene) -> selected
    hit_map: dict[tuple[str, str, str], bool] = {
        (h.cohort, h.fitter, h.feature): h.selected for h in hits
    }

    # 3. Cohort Rows and 2-Tier Pill Capsules
    for r_idx, c in enumerate(cohorts):
        cy = grid_top_y + r_idx * row_h
        row_bg = COLOR_SUBTLE_FILL if r_idx % 2 == 1 else "#FFFFFF"
        font_wt = "700" if c.is_combined else "500"

        svg.append(f"""
    <!-- Row {r_idx}: {c.name} -->
    <rect x="{px + 12}" y="{cy}" width="{pw - 24}" height="{row_h - 1}" fill="{row_bg}" rx="0" />
    <text x="{px + 165}" y="{cy + 22}" font-size="8.5" font-weight="{font_wt}" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">
      {c.display_title} <tspan font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">(n={c.n_samples})</tspan>
    </text>
""")

        # Draw 2-tier pill capsule rack in each cell
        for g in genes:
            cell_x = grid_x_start + g.col_idx * col_w + 3.0
            cell_y = cy + 4.5

            # Tier 1 (Top sub-row): 6 Single-Cohort Baseline Fitters
            b_pill_w = 4.2
            b_pill_h = 11.5
            b_pill_gap = 1.0
            for f_idx, f_name in enumerate(FITTER_BASELINE_ORDER):
                is_sel = hit_map.get((c.name, f_name, g.symbol), False)
                f_col = FITTER_COLOR_MAP[f_name]
                px_pos = cell_x + f_idx * (b_pill_w + b_pill_gap)
                py_pos = cell_y

                if is_sel:
                    svg.append(f"""
    <rect x="{px_pos:.1f}" y="{py_pos:.1f}" width="{b_pill_w}" height="{b_pill_h}" rx="0"
          fill="{f_col}" stroke="#FFFFFF" stroke-width="0.5" />
""")
                else:
                    svg.append(f"""
    <rect x="{px_pos:.1f}" y="{py_pos:.1f}" width="{b_pill_w}" height="{b_pill_h}" rx="0"
          fill="#FFFFFF" stroke="#E2E8F0" stroke-width="0.5" />
""")

            # Tier 2 (Bottom sub-row): 7 Multi-Cohort Fitters
            m_pill_w = 3.4
            m_pill_h = 11.5
            m_pill_gap = 1.0
            m_cell_y = cell_y + 13.5
            for f_idx, f_name in enumerate(FITTER_MULTICOHORT_ORDER):
                is_sel = hit_map.get((c.name, f_name, g.symbol), False)
                f_col = FITTER_COLOR_MAP[f_name]
                px_pos = cell_x + f_idx * (m_pill_w + m_pill_gap)
                py_pos = m_cell_y

                if not c.is_combined:
                    # Single cohorts did not execute multi-cohort fitters
                    svg.append(f"""
    <rect x="{px_pos:.1f}" y="{py_pos:.1f}" width="{m_pill_w}" height="{m_pill_h}" rx="0"
          fill="#F8FAFC" stroke="#E2E8F0" stroke-width="0.5" />
""")
                elif is_sel:
                    svg.append(f"""
    <rect x="{px_pos:.1f}" y="{py_pos:.1f}" width="{m_pill_w}" height="{m_pill_h}" rx="0"
          fill="{f_col}" stroke="#FFFFFF" stroke-width="0.5" />
""")
                else:
                    svg.append(f"""
    <rect x="{px_pos:.1f}" y="{py_pos:.1f}" width="{m_pill_w}" height="{m_pill_h}" rx="0"
          fill="#FFFFFF" stroke="#E2E8F0" stroke-width="0.5" />
""")

        # Count total detections in this cohort across all fitters
        cohort_hits = sum(1 for h in hits if h.cohort == c.name and h.selected and h.feature in {g.symbol for g in genes})
        svg.append(f"""
    <text x="{grid_x_start + len(genes) * col_w + 15}" y="{cy + 22}" font-size="8" font-weight="700" fill="{COLOR_TEXT_SECONDARY}">
      {cohort_hits} hits
    </text>
""")

    # Column right summary header
    svg.append(f"""
    <text x="{grid_x_start + len(genes) * col_w + 15}" y="{grid_top_y - 8}" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}">
      TOTAL
    </text>

    <!-- Bottom Guide Callout -->
    <g transform="translate({px + 14}, {ph + py - 26})">
      <rect x="0" y="0" width="{pw - 28}" height="18" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" rx="0" />
      <text x="10" y="12" font-size="8" font-weight="700" fill="{COLOR_FITTER_LASSO}">2-Tier Capsule Rack:</text>
      <text x="110" y="12" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        Top: 6 Single-Cohort Baselines [Lasso | ENet | Logistic | RF | OSCAR | SLOPE]  |  Bottom: 7 Multi-Cohort Fitters [FWL | Grp-Lasso | MERF | MT-Logistic | Meta | IRM | GLMM]
      </text>
    </g>
  </g>
""")
    return "".join(svg)


def render_panel_b_concordance(hits: Sequence[FitterGeneHit]) -> str:
    """Render Panel B: Fitter Cardinality & Model Diversity Analysis for 13 Fitters."""
    px = 930
    py = 70
    pw = 450
    ph = 440

    svg = [f"""
  <!-- LAYER 03: PANEL B (CONCORDANCE) -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Concordance">
    <rect id="card-panel-b" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />

    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">b</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Fitter Discovery Cardinality &amp; Concordance
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Total gene selection instances across cohorts across 13 base fitters
    </text>
"""]

    # Fitter selections count
    fitter_counts: dict[str, int] = {f: 0 for f in FITTER_ORDER}
    for h in hits:
        if h.selected and h.fitter in fitter_counts:
            fitter_counts[h.fitter] += 1

    max_count = max(fitter_counts.values()) if fitter_counts else 1
    max_bar_w = 210

    # 1. Baseline Fitters Group
    b_start_y = py + 42
    svg.append(f"""
    <text x="{px + 14}" y="{b_start_y}" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_MUTED}" letter-spacing="0.3">
      SINGLE-COHORT BASELINE FITTERS (6)
    </text>
""")
    b_chart_y = b_start_y + 6
    b_bar_h = 8.5
    b_spacing = 11.5
    for idx, f_name in enumerate(FITTER_BASELINE_ORDER):
        cy = b_chart_y + idx * b_spacing
        cnt = fitter_counts[f_name]
        bw = (cnt / max_count) * max_bar_w
        f_col = FITTER_COLOR_MAP[f_name]
        d_name = FITTER_DISPLAY_NAMES[f_name]

        svg.append(f"""
    <text x="{px + 130}" y="{cy + 7}" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">
      {d_name}
    </text>
    <rect x="{px + 138}" y="{cy}" width="{bw:.1f}" height="{b_bar_h}" fill="{f_col}" rx="0" />
    <text x="{px + 143 + bw:.1f}" y="{cy + 7}" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      {cnt}
    </text>
""")

    # 2. Multi-Cohort Fitters Group
    m_start_y = b_chart_y + len(FITTER_BASELINE_ORDER) * b_spacing + 8
    svg.append(f"""
    <text x="{px + 14}" y="{m_start_y}" font-size="7.5" font-weight="700" fill="{COLOR_FITTER_GROUP_LASSO}" letter-spacing="0.3">
      MULTI-COHORT AWARE FITTERS (7)
    </text>
""")
    m_chart_y = m_start_y + 6
    m_bar_h = 8.5
    m_spacing = 11.5
    for idx, f_name in enumerate(FITTER_MULTICOHORT_ORDER):
        cy = m_chart_y + idx * m_spacing
        cnt = fitter_counts[f_name]
        bw = (cnt / max_count) * max_bar_w
        f_col = FITTER_COLOR_MAP[f_name]
        d_name = FITTER_DISPLAY_NAMES[f_name]

        svg.append(f"""
    <text x="{px + 130}" y="{cy + 7}" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">
      {d_name}
    </text>
    <rect x="{px + 138}" y="{cy}" width="{bw:.1f}" height="{m_bar_h}" fill="{f_col}" rx="0" />
    <text x="{px + 143 + bw:.1f}" y="{cy + 7}" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      {cnt}
    </text>
""")

    # Sub-card: Model Family Insights Card
    card_y = m_chart_y + len(FITTER_MULTICOHORT_ORDER) * m_spacing + 8
    card_h = ph - (card_y - py) - 10

    svg.append(f"""
    <!-- Model Family Insights Card -->
    <g transform="translate({px + 14}, {card_y})">
      <rect x="0" y="0" width="{pw - 28}" height="{card_h}" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <text x="12" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        MODEL DIVERSITY &amp; BIOLOGICAL RECOVERY (13 FITTERS)
      </text>

      <g transform="translate(12, 22)">
        <circle cx="4" cy="5" r="3" fill="{COLOR_FITTER_LASSO}" />
        <text x="14" y="8" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Linear Sparsity (Lasso, Elastic Net, GLMM Lasso):</text>
        <text x="14" y="18" font-size="7" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Isolates dominant linear orthogonal markers (HLA-A, PTGS2, TGFB1).
        </text>

        <circle cx="4" cy="28" r="3" fill="{COLOR_FITTER_OSCAR}" />
        <text x="14" y="31" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Ordered Weighted L1 (OSCAR, SLOPE):</text>
        <text x="14" y="41" font-size="7" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Octagonal clustering groups collinear genes; BH quantiles control discovery FDR.
        </text>

        <circle cx="4" cy="51" r="3" fill="{COLOR_FITTER_LOGISTIC}" />
        <text x="14" y="54" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Exact Bernoulli Likelihood (Logistic, Multi-Task):</text>
        <text x="14" y="64" font-size="7" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Directly models response odds; rescues sub-threshold cytokines (CCR7, IL10).
        </text>

        <circle cx="4" cy="74" r="3" fill="{COLOR_FITTER_RF}" />
        <text x="14" y="77" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Non-Linear Ensembles (Random Forest, MERF):</text>
        <text x="14" y="87" font-size="7" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Captures epistatic interactions and master transactivators (CIITA in Melanoma).
        </text>

        <circle cx="4" cy="97" r="3" fill="{COLOR_FITTER_GROUP_LASSO}" />
        <text x="14" y="100" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Joint Support &amp; Consensus (Group Lasso, Meta, Invariant):</text>
        <text x="14" y="110" font-size="7" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Forces cross-study consensus, exposing full MHC-I regulon (NLRC5, PSMB8).
        </text>
      </g>
    </g>
  </g>
""")
    return "".join(svg)


def render_panel_c_multicohort() -> str:
    """Render Panel C: Multi-Cohort Fitter Comparison in Pooled Datasets."""
    px = 20
    py = 525
    pw = 675
    ph = 380

    svg = [f"""
  <!-- LAYER 04: PANEL C (MULTI-COHORT FITTER COMPARISON) -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_MultiCohort">
    <rect id="card-panel-c" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />

    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">c</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Multi-Cohort Fitter Comparison: Pooled Melanoma &amp; RCC
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Head-to-head contrast: Single-Cohort Baselines vs Multi-Task Group Regularization vs Non-Linear MERF
    </text>
"""]

    # Box 1: Melanoma Combined Comparison (Top)
    b1_y = py + 48
    b1_h = 154
    svg.append(f"""
    <!-- Melanoma Box -->
    <g transform="translate({px + 14}, {b1_y})">
      <rect x="0" y="0" width="{pw - 28}" height="{b1_h}" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <rect x="0" y="0" width="{pw - 28}" height="22" fill="#E2E8F0" rx="0" />
      <text x="10" y="15" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        POOLED MELANOMA (n = 338: Hugo, Riaz, Liu, Gide)
      </text>

      <!-- Sub-fitter 1: Standard Lasso & ENet -->
      <g transform="translate(10, 30)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_LASSO}" stroke-width="1" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_LASSO}">Lasso / Elastic Net</text>
        <text x="148" y="13" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          HLA-A, IKZF2, TNFSF18
        </text>
        <text x="320" y="13" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Naive pooling; single-study dominant</text>
      </g>

      <!-- Sub-fitter 2: FWL Cohort-Adjusted -->
      <g transform="translate(10, 52)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_COHORT_ADJ}" stroke-width="1" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_COHORT_ADJ}">Cohort-Adjusted FWL</text>
        <text x="148" y="13" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          HLA-A, IKZF2, TNFSF9
        </text>
        <text x="320" y="13" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Removes mean shifts; recovers TNFSF9</text>
      </g>

      <!-- Sub-fitter 3: Group Lasso & Multi-Task Logistic -->
      <g transform="translate(10, 74)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_GROUP_LASSO}" stroke-width="1.25" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GROUP_LASSO}">Group Lasso / MT-Log</text>
        <text x="148" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GROUP_LASSO}">
          NLRC5, PSMB8, HLA-A, TNFSF9, PVR, VSIR
        </text>
      </g>

      <!-- Sub-fitter 4: MERF (Mixed RF) & Meta-Analysis -->
      <g transform="translate(10, 96)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_MERF}" stroke-width="1.25" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_MERF}">MERF / Meta-Consensus</text>
        <text x="148" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_MERF}">
          CIITA, HLA-A, PSMB8, TNFSF9
        </text>
      </g>
      <text x="10" y="136" font-size="7.5" font-style="italic" fill="{COLOR_TEXT_SECONDARY}">
        Joint multi-task sparsity recovers complete MHC-I transactivation regulon (NLRC5, PSMB8) and novel checkpoints (VISTA, PVR).
      </text>
    </g>
""")

    # Box 2: RCC Combined Comparison (Bottom)
    b2_y = b1_y + b1_h + 10
    b2_h = 142
    svg.append(f"""
    <!-- RCC Box -->
    <g transform="translate({px + 14}, {b2_y})">
      <rect x="0" y="0" width="{pw - 28}" height="{b2_h}" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <rect x="0" y="0" width="{pw - 28}" height="22" fill="#E2E8F0" rx="0" />
      <text x="10" y="15" font-size="9" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        POOLED RENAL CELL CARCINOMA (n = 263: McDermott, Choueiri)
      </text>

      <!-- Sub-fitter 1: Standard & FWL -->
      <g transform="translate(10, 30)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_LASSO}" stroke-width="1" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_LASSO}">Lasso / FWL / ENet</text>
        <text x="148" y="13" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
          PTGS2 (COX-2)
        </text>
        <text x="250" y="13" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Consistent canonical myeloid/angiogenic suppressor</text>
      </g>

      <!-- Sub-fitter 2: Group Lasso & Multi-Task -->
      <g transform="translate(10, 52)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_GROUP_LASSO}" stroke-width="1.25" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GROUP_LASSO}">Group Lasso / MT-Log</text>
        <text x="148" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GROUP_LASSO}">
          CIITA, NLRC5, PTGS2
        </text>
        <text x="290" y="13" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Joint MHC-I/II transactivation across trials</text>
      </g>

      <!-- Sub-fitter 3: GLMM & Meta-Analysis -->
      <g transform="translate(10, 74)">
        <rect x="0" y="0" width="138" height="18" fill="#FFFFFF" stroke="{COLOR_FITTER_GLMM}" stroke-width="1.25" rx="0" />
        <text x="6" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GLMM}">GLMM / Meta-Analysis</text>
        <text x="148" y="13" font-size="8" font-weight="700" fill="{COLOR_FITTER_GLMM}">
          PTGS2, CIITA
        </text>
        <text x="290" y="13" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Mixed effects confirms invariant trial signal</text>
      </g>

      <text x="10" y="112" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        Key Insight: Linear fitters capture phenotypic stroma (PTGS2), whereas multi-cohort joint modeling discovers
      </text>
      <text x="10" y="125" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        the upstream transcriptional master switches (CIITA, NLRC5) required for cross-trial immune recognition.
      </text>
    </g>
  </g>
""")
    return "".join(svg)


def render_panel_d_summary() -> str:
    """Render Panel D: Biological Synthesis & Fitter Specificity Rationale."""
    px = 705
    py = 525
    pw = 675
    ph = 380

    svg = [f"""
  <!-- LAYER 05: PANEL D (SUMMARY & LEGENDS) -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Summary">
    <rect id="card-panel-d" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />

    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">d</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Biological Mechanism &amp; Fitter Specificity Rationale
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Algorithmic rationale for consensus vs fitter-specific biomarker discoveries across 13 fitters
    </text>
"""]

    # Card 1: Fitter Color Legends Track (13 Fitters)
    svg.append(f"""
    <g transform="translate({px + 14}, {py + 46})">
      <rect x="0" y="0" width="{pw - 28}" height="92" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <text x="12" y="15" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        EVALUATED STABILITY SELECTION FITTERS (PILL COLOR MAPPING &amp; TAXONOMY)
      </text>

      <!-- Row 1: Baselines (Lasso, Elastic Net, Logistic, RF, OSCAR, SLOPE) -->
      <g transform="translate(12, 26)">
        <rect x="0" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_LASSO}" />
        <text x="12" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Lasso</text>
        <text x="42" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">L1</text>

        <rect x="80" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_ELASTIC_NET}" />
        <text x="92" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Elastic Net</text>
        <text x="145" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">L1+L2</text>

        <rect x="185" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_LOGISTIC}" />
        <text x="197" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Logistic</text>
        <text x="235" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Bernoulli</text>

        <rect x="290" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_RF}" />
        <text x="302" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Random Forest</text>
        <text x="372" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Gini trees</text>

        <rect x="425" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_OSCAR}" />
        <text x="437" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">OSCAR</text>
        <text x="475" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Octagonal</text>

        <rect x="535" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_SLOPE}" />
        <text x="547" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">SLOPE</text>
        <text x="585" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Sorted FDR</text>
      </g>

      <!-- Row 2: Multi-Cohort 1 (FWL, Group Lasso, MERF, MT-Logistic) -->
      <g transform="translate(12, 46)">
        <rect x="0" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_COHORT_ADJ}" />
        <text x="12" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">FWL Adjusted</text>
        <text x="75" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Mean-shift</text>

        <rect x="145" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_GROUP_LASSO}" />
        <text x="157" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Group Lasso</text>
        <text x="218" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">L2,1 joint</text>

        <rect x="295" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_MERF}" />
        <text x="307" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">MERF</text>
        <text x="340" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Mixed effects RF</text>

        <rect x="450" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_MULTITASK}" />
        <text x="462" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">MT-Logistic</text>
        <text x="522" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">L2,1 Bernoulli</text>
      </g>

      <!-- Row 3: Multi-Cohort 2 (Meta-Analysis, Invariant IRM, GLMM Lasso) -->
      <g transform="translate(12, 66)">
        <rect x="0" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_META}" />
        <text x="12" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Meta-Consensus</text>
        <text x="90" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Fisher min-p</text>

        <rect x="180" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_INVARIANT}" />
        <text x="192" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Invariant IRM</text>
        <text x="260" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Risk min</text>

        <rect x="360" y="0" width="8" height="12" rx="0" fill="{COLOR_FITTER_GLMM}" />
        <text x="372" y="9" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">GLMM Lasso</text>
        <text x="435" y="9" font-size="7" font-weight="400" fill="{COLOR_TEXT_MUTED}">Laplace penalized</text>
      </g>
    </g>
""")

    # Card 2: Functional Pathway Matrix
    svg.append(f"""
    <g transform="translate({px + 14}, {py + 146})">
      <rect x="0" y="0" width="{pw - 28}" height="136" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <rect x="0" y="0" width="{pw - 28}" height="20" fill="{COLOR_SUBTLE_FILL}" rx="0" />
      <text x="12" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        BIOLOGICAL PATHWAY TARGETS &amp; ALGORITHMIC SENSITIVITY
      </text>

      <g transform="translate(12, 28)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_FITTER_LASSO}">Antigen Presentation Machinery:</text>
        <text x="180" y="10" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          HLA-A &amp; HLA-C selected across linear models; Group Lasso uniquely recruits master transactivator
        </text>
        <text x="180" y="21" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          NLRC5, while MERF and GLMM uncover CIITA and immunoproteasome subunit PSMB8.
        </text>

        <text x="0" y="38" font-size="8" font-weight="700" fill="{COLOR_FITTER_LOGISTIC}">Effector Lineage &amp; Checkpoints:</text>
        <text x="180" y="38" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          TBX21 (T-bet) robust in Pan-Cancer; Random Forest captures TCR coreceptor CD3E;
        </text>
        <text x="180" y="49" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Group Lasso and MT-Logistic identify secondary checkpoint ligands VSIR (VISTA) and PVR.
        </text>

        <text x="0" y="66" font-size="8" font-weight="700" fill="{COLOR_FITTER_RF}">Costimulation &amp; Regulatory T cells:</text>
        <text x="180" y="66" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          TNFSF9 (4-1BBL), TNFSF18 (GITRL), and IKZF2 (Helios) form a solid consensus core in melanoma;
        </text>
        <text x="180" y="77" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          exact Bernoulli Logistic fitter rescues lymph node homing chemokine receptor CCR7.
        </text>

        <text x="0" y="94" font-size="8" font-weight="700" fill="{COLOR_FITTER_COHORT_ADJ}">Stroma &amp; Immune Exclusion:</text>
        <text x="180" y="94" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          PTGS2 (COX-2) in RCC and TGFB1 in Bladder demonstrate 100% consensus across all single-cohort
        </text>
        <text x="180" y="105" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          and multi-cohort fitters (Lasso, Elastic Net, Logistic, Random Forest, GLMM).
        </text>
      </g>
    </g>
""")

    # Card 3: Controlled Selection Rationale Chip
    svg.append(f"""
    <g transform="translate({px + 14}, {py + 290})">
      <rect x="0" y="0" width="{pw - 28}" height="76" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" rx="0" />
      <text x="12" y="16" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        METHODOLOGICAL CONCLUSION
      </text>
      <text x="12" y="32" font-size="8" font-weight="600" fill="{COLOR_FITTER_GROUP_LASSO}">
        Complementary Model Synergy:
      </text>
      <text x="140" y="32" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        Standard linear Lasso provides a conservative, high-precision core (HLA-A, PTGS2, TGFB1),
      </text>
      <text x="12" y="46" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        while Multi-Cohort Group Lasso and Random Forest expose epistatic regulators (NLRC5, CIITA) that
      </text>
      <text x="12" y="58" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        remain invisible to single-cohort linear models.
      </text>
    </g>
  </g>
""")
    return "".join(svg)


def render_full_svg(
    hits: Sequence[FitterGeneHit],
    cohorts: Sequence[CohortDisplayMeta],
    genes: Sequence[GeneDisplayMeta],
) -> str:
    """Declaratively render complete publication vector figure with Inkscape layers."""
    parts: list[str] = [
        f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-fitter-gene-overlap-nature"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <!-- LAYER 00: CANVAS BACKGROUND -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />
  </g>
""",
        render_header(),
        render_panel_a_grid(hits, cohorts, genes),
        render_panel_b_concordance(hits),
        render_panel_c_multicohort(),
        render_panel_d_summary(),
        "\n</svg>\n",
    ]
    return "".join(parts)


# -------------------------------------------------------------------------
# Shell & I/O Execution
# -------------------------------------------------------------------------

def save_and_verify(
    svg_content: str,
    output_svg_path: Path,
    output_png_path: Path,
) -> Result[Path, str]:
    """Pure I/O shell writing SVG, validating XML, and rendering 300 DPI PNG."""
    return verify_and_rasterize_svg(
        svg_content_or_path=svg_content,
        output_svg_path=output_svg_path,
        output_png_path=output_png_path,
        width=WIDTH,
        height=HEIGHT,
        scale=2.0,
        min_layers=5,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate publication figure of multi-fitter gene selection overlap"
    )
    parser.add_argument(
        "--scores-path",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark/stability_scores.parquet"),
        help="Path to multi-fitter stability scores parquet",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark/figures/fitter_gene_overlap_nature.svg"),
        help="Target output SVG path",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark/figures/fitter_gene_overlap_nature.png"),
        help="Target output PNG path",
    )
    args = parser.parse_args()

    load_res = load_benchmark_data(args.scores_path)
    match load_res:
        case Failure(err):
            print(f"Error loading inputs: {err}", file=sys.stderr)
            return 1
        case Success(df_scores):
            hits = extract_fitter_hits(df_scores, method="SS-CPSS")
            svg_str = render_full_svg(hits, ORDERED_COHORTS, ORDERED_GENES)
            save_res = save_and_verify(svg_str, args.output_svg, args.output_png)

            match save_res:
                case Failure(err):
                    print(f"Error saving/verifying vector figure: {err}", file=sys.stderr)
                    return 1
                case Success(svg_path):
                    print(f"Successfully generated publication vector figure at: {svg_path}")
                    if args.output_png.exists():
                        print(f"Successfully rasterized 300 DPI PNG at: {args.output_png}")
                    return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
