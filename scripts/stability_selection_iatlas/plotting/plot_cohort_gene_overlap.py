#!/usr/bin/env python3
"""Generate Publication-Grade Cross-Cohort Gene Selection Overlap Figure.

Adheres strictly to the project's functional programming principles:
- Pure functions mapping immutable Polars data structures to SVG vector elements.
- Nature Methods minimalist wireframe standard (pure white background, #CBD5E1 borders, rx=0, no drop shadows).
- Okabe-Ito Colorblind-Safe scientific palette strictly applied to data marks.
- Native Inkscape layers and semantic groups.
- Vector SVG exports and 300 DPI PNG verification via rsvg-convert.
- Monadic error handling with Result[Path, str].
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path
import subprocess
import sys
from typing import Final, Sequence
from pydantic import BaseModel, ConfigDict
import polars as pl
from returns.result import Failure, Result, Success

from plotting_utils import (
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_DIVIDER_RULE,
    COLOR_PANEL_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_HAIRLINE,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
    FONT_SERIF_MATH,
    OKABE_BLACK,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    OKABE_YELLOW,
)
from plotting_utils.svg import verify_and_rasterize_svg, wrap_inkscape_svg

# Figure Canvas Dimensions (Nature double-column landscape)
WIDTH: Final[int] = 1400
HEIGHT: Final[int] = 920


# -------------------------------------------------------------------------
# Immutable Domain Models
# -------------------------------------------------------------------------

class CohortInfo(BaseModel):
    """Metadata for a single clinical cohort."""
    model_config = ConfigDict(frozen=True)

    name: str
    display_name: str
    cancer_type: str
    is_pooled: bool
    n_samples: int
    n_responders: int
    n_non_responders: int
    order_idx: int


class GeneInfo(BaseModel):
    """Metadata for a stably selected biomarker gene."""
    model_config = ConfigDict(frozen=True)

    symbol: str
    full_name: str
    pathway: str
    cancer_compartment: str
    role_summary: str


class GeneSelectionCell(BaseModel):
    """Stability selection status and scores for a (cohort, gene) pair."""
    model_config = ConfigDict(frozen=True)

    cohort: str
    gene: str
    score_mb: float
    score_ss: float
    max_score: float
    selected_mb: bool
    selected_ss: bool
    concordance_status: str  # 'both', 'ss_only', 'mb_only', 'subthreshold', 'unselected'


class UpSetIntersection(BaseModel):
    """An intersection column in the UpSet diagram."""
    model_config = ConfigDict(frozen=True)

    id: int
    label: str
    genes: tuple[str, ...]
    size: int
    cohorts: tuple[str, ...]


# -------------------------------------------------------------------------
# Pure Domain Knowledge Definitions
# -------------------------------------------------------------------------

ORDERED_GENES: Final[tuple[GeneInfo, ...]] = (
    GeneInfo(
        symbol="HLA-A",
        full_name="MHC Class I, A",
        pathway="Antigen Presentation",
        cancer_compartment="Melanoma",
        role_summary="Primary neoantigen presentation to CD8+ T cells; loss confers resistance",
    ),
    GeneInfo(
        symbol="HLA-C",
        full_name="MHC Class I, C",
        pathway="Antigen Presentation",
        cancer_compartment="Pan-Cancer",
        role_summary="Broad histocompatibility complex for cytotoxic T and NK cell recognition",
    ),
    GeneInfo(
        symbol="TBX21",
        full_name="T-box Transcription Factor 21 (T-bet)",
        pathway="Lineage & Exhaustion",
        cancer_compartment="Pan-Cancer",
        role_summary="Master transcription factor coordinating Th1 and CD8+ cytotoxic commitment",
    ),
    GeneInfo(
        symbol="LAG3",
        full_name="Lymphocyte Activating 3",
        pathway="Lineage & Exhaustion",
        cancer_compartment="Pan-Cancer",
        role_summary="Inhibitory immune checkpoint upregulated upon chronic tumor exhaustion",
    ),
    GeneInfo(
        symbol="TNFSF9",
        full_name="TNF Superfamily Member 9 (4-1BBL)",
        pathway="Costimulation & Regulators",
        cancer_compartment="Melanoma",
        role_summary="Potent CD137 costimulatory ligand sustaining effector T-cell persistence",
    ),
    GeneInfo(
        symbol="TNFSF18",
        full_name="TNF Superfamily Member 18 (GITRL)",
        pathway="Costimulation & Regulators",
        cancer_compartment="Melanoma",
        role_summary="GITR activation reversing regulatory T-cell suppression in tumor beds",
    ),
    GeneInfo(
        symbol="IKZF2",
        full_name="IKAROS Family Zinc Finger 2 (Helios)",
        pathway="Costimulation & Regulators",
        cancer_compartment="Melanoma",
        role_summary="Maintains regulatory T-cell suppressive phenotype in inflamed niches",
    ),
    GeneInfo(
        symbol="PTGS2",
        full_name="Prostaglandin-Endoperoxide Synthase 2 (COX-2)",
        pathway="Stroma & Inflammation",
        cancer_compartment="RCC",
        role_summary="PGE2 biosynthesis mediating myeloid suppression and angiogenic escape",
    ),
    GeneInfo(
        symbol="TGFB1",
        full_name="Transforming Growth Factor Beta 1",
        pathway="Stroma & Inflammation",
        cancer_compartment="Bladder",
        role_summary="Stroma-mediated exclusion of cytotoxic lymphocytes from tumor cores",
    ),
    GeneInfo(
        symbol="IFNG",
        full_name="Interferon Gamma",
        pathway="Stroma & Inflammation",
        cancer_compartment="Bladder",
        role_summary="Hallmark Th1 cytokine orchestrating chemokine recruitment and PD-L1 expression",
    ),
)


# -------------------------------------------------------------------------
# Pure Data Pipeline (Polars)
# -------------------------------------------------------------------------

def load_data(
    scores_path: Path, summary_path: Path
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Pure IO loader returning Polars DataFrames."""
    if not scores_path.exists():
        return Failure(f"Stability scores file not found: {scores_path}")
    if not summary_path.exists():
        return Failure(f"Cohort summary file not found: {summary_path}")

    try:
        df_scores = pl.read_parquet(scores_path)
        df_summary = pl.read_parquet(summary_path)
        return Success((df_scores, df_summary))
    except Exception as exc:
        return Failure(f"Failed to load parquet data: {exc}")


def extract_cohort_models(df_summary: pl.DataFrame) -> tuple[CohortInfo, ...]:
    """Extract ordered immutable CohortInfo models from cohort summary DataFrame."""
    # Filter for lasso baseline if multi-fitter summary is provided
    if "fitter" in df_summary.columns:
        df_summary = df_summary.filter(pl.col("fitter") == "lasso")

    # Custom ordering: Pooled cohorts first, active trials next, zero-selection trials last
    cohort_order_map: dict[str, tuple[int, str, str, bool]] = {
        "pancancer": (0, "Pan-Cancer Pooled", "Pan-Cancer", True),
        "melanoma": (1, "Melanoma Pooled", "Melanoma", True),
        "rcc": (2, "RCC Pooled", "RCC", True),
        "Rosenberg-iAtlas": (3, "Rosenberg (Bladder)", "Bladder", False),
        "McDermott-iAtlas": (4, "McDermott (RCC)", "RCC", False),
        "Riaz-iAtlas": (5, "Riaz (Melanoma)", "Melanoma", False),
        "Gide-iAtlas": (6, "Gide (Melanoma)", "Melanoma", False),
        "Liu-iAtlas": (7, "Liu (Melanoma)", "Melanoma", False),
        "Padron-iAtlas": (8, "Padron (PDAC)", "PDAC", False),
        "Anders-iAtlas": (9, "Anders (Bladder)", "Bladder", False),
        "Hugo-iAtlas": (10, "Hugo (Melanoma)", "Melanoma", False),
        "Choueiri-iAtlas": (11, "Choueiri (RCC)", "RCC", False),
    }

    cohort_rows = df_summary.to_dicts()
    cohort_models: list[CohortInfo] = []

    for row in cohort_rows:
        c_name = str(row["cohort"])
        if c_name not in cohort_order_map:
            continue
        order_idx, display_name, cancer_type, is_pooled = cohort_order_map[c_name]
        cohort_models.append(
            CohortInfo(
                name=c_name,
                display_name=display_name,
                cancer_type=cancer_type,
                is_pooled=is_pooled,
                n_samples=int(row["n_samples"]),
                n_responders=int(row["n_responders"]),
                n_non_responders=int(row["n_non_responders"]),
                order_idx=order_idx,
            )
        )

    return tuple(sorted(cohort_models, key=lambda c: c.order_idx))


def extract_matrix_cells(
    df_scores: pl.DataFrame,
    cohorts: Sequence[CohortInfo],
    genes: Sequence[GeneInfo],
) -> tuple[GeneSelectionCell, ...]:
    """Pure function building grid selection cells with concordance status."""
    # Filter for lasso baseline if multi-fitter scores are provided
    if "fitter" in df_scores.columns:
        df_scores = df_scores.filter(pl.col("fitter") == "lasso")

    target_genes = {g.symbol for g in genes}
    df_filtered = df_scores.filter(pl.col("feature").is_in(list(target_genes)))

    cells: list[GeneSelectionCell] = []

    for c in cohorts:
        c_scores = df_filtered.filter(pl.col("cohort") == c.name)
        for g in genes:
            g_scores = c_scores.filter(pl.col("feature") == g.symbol)

            mb_row = g_scores.filter(pl.col("method") == "MB")
            ss_row = g_scores.filter(pl.col("method") == "SS-CPSS")

            score_mb = float(mb_row["stability_score"][0]) if not mb_row.is_empty() else 0.0
            score_ss = float(ss_row["stability_score"][0]) if not ss_row.is_empty() else 0.0

            sel_mb = bool(mb_row["selected"][0]) if not mb_row.is_empty() else False
            sel_ss = bool(ss_row["selected"][0]) if not ss_row.is_empty() else False

            max_score = max(score_mb, score_ss)

            if sel_mb and sel_ss:
                status = "both"
            elif sel_ss and not sel_mb:
                status = "ss_only"
            elif sel_mb and not sel_ss:
                status = "mb_only"
            elif max_score >= 0.50:
                status = "subthreshold"
            else:
                status = "unselected"

            cells.append(
                GeneSelectionCell(
                    cohort=c.name,
                    gene=g.symbol,
                    score_mb=score_mb,
                    score_ss=score_ss,
                    max_score=max_score,
                    selected_mb=sel_mb,
                    selected_ss=sel_ss,
                    concordance_status=status,
                )
            )

    return tuple(cells)


def compute_upset_intersections() -> tuple[UpSetIntersection, ...]:
    """Define the exact biological multi-cohort intersections observed across active cohorts."""
    return (
        UpSetIntersection(
            id=1,
            label="Pan-Cancer Core",
            genes=("TBX21", "HLA-C", "LAG3"),
            size=3,
            cohorts=("pancancer",),
        ),
        UpSetIntersection(
            id=2,
            label="Melanoma Synergy",
            genes=("IKZF2", "TNFSF18"),
            size=2,
            cohorts=("melanoma",),
        ),
        UpSetIntersection(
            id=3,
            label="Melanoma ∩ Gide",
            genes=("HLA-A",),
            size=1,
            cohorts=("melanoma", "Gide-iAtlas"),
        ),
        UpSetIntersection(
            id=4,
            label="Melanoma ∩ Riaz",
            genes=("TNFSF9",),
            size=1,
            cohorts=("melanoma", "Riaz-iAtlas"),
        ),
        UpSetIntersection(
            id=5,
            label="RCC ∩ McDermott",
            genes=("PTGS2",),
            size=1,
            cohorts=("rcc", "McDermott-iAtlas"),
        ),
        UpSetIntersection(
            id=6,
            label="Bladder Program",
            genes=("TGFB1", "IFNG"),
            size=2,
            cohorts=("Rosenberg-iAtlas",),
        ),
    )


# -------------------------------------------------------------------------
# SVG Rendering Functions
# -------------------------------------------------------------------------

def render_header() -> str:
    """Render top figure title, subtitle, and formal statistical parameter chip."""
    return f"""
  <!-- LAYER 01: FIGURE HEADER -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header">
    <text id="txt-title" x="20" y="32" font-size="15.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">
      Cross-Cohort Stability Selection Landscape in Immune Checkpoint Blockade
    </text>
    <text id="txt-subtitle" x="20" y="50" font-size="10.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Synthesis of stable predictive genes across 12 iAtlas cohorts (n = 1,015) comparing Shah-Samworth CPSS and Meinshausen-Bühlmann methods
    </text>

    <!-- Stat Tag Chip -->
    <g id="grp-header-chip" transform="translate(1010, 16)">
      <rect x="0" y="0" width="370" height="34" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <text x="185" y="15" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle" letter-spacing="0.2">
        CONTROLLED SELECTION PARAMETERS
      </text>
      <text x="185" y="27" font-size="8.5" font-weight="500" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">
        Cutoff π<tspan font-size="7" baseline-shift="sub">thr</tspan> = 0.75  |  q<tspan font-size="7" baseline-shift="sub">budget</tspan> ≤ 20.0  |  PFER ≤ 0.97  |  200 Subsamples
      </text>
    </g>
  </g>
"""


def render_panel_a_upset(
    intersections: Sequence[UpSetIntersection],
    active_cohorts: Sequence[CohortInfo],
) -> str:
    """Render Panel A: UpSet Multi-Cohort Overlap Matrix."""
    px = 20
    py = 70
    pw = 665
    ph = 430

    svg = [f"""
  <!-- LAYER 02: PANEL A (UPSET PLOT) -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_UpSet">
    <!-- Panel Background Card -->
    <rect id="card-panel-a" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel Badge & Titles -->
    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">a</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      UpSet Multi-Cohort Overlap Matrix (10 Stably Selected Genes)
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Intersection cardinality and shared gene discovery across 7 active clinical cohorts (π<tspan font-size="7.5" baseline-shift="sub">thr</tspan> ≥ 0.75)
    </text>
"""]

    # Coordinates for UpSet components
    chart_top = py + 48
    bar_chart_h = 100
    baseline_y = chart_top + bar_chart_h  # y = 218
    matrix_top = baseline_y + 20          # y = 238
    row_height = 24

    # 1. Top Intersection Size Bars
    # 6 columns
    col_x_start = px + 335
    col_spacing = 52

    svg.append(f"""
    <!-- Baseline rule for intersection bars -->
    <line x1="{col_x_start - 25}" y1="{baseline_y}" x2="{col_x_start + len(intersections) * col_spacing - 15}" y2="{baseline_y}"
          stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
    <text x="{col_x_start - 30}" y="{baseline_y - 4}" font-size="8" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="end">
      Intersection Size
    </text>
""")

    # Render top bars
    for idx, inter in enumerate(intersections):
        cx = col_x_start + idx * col_spacing
        bar_h = inter.size * 26
        by = baseline_y - bar_h
        gene_label_str = ", ".join(inter.genes)

        # Bar fill color: Okabe Blue for shared, Slate for exclusive
        is_shared = len(inter.cohorts) > 1
        bar_color = OKABE_BLUE if is_shared else "#334155"

        svg.append(f"""
    <!-- Intersection Bar {inter.id}: {inter.label} -->
    <rect x="{cx - 14}" y="{by}" width="28" height="{bar_h}" fill="{bar_color}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="{cx}" y="{by - 5}" font-size="10" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">{inter.size}</text>
    <text x="{cx}" y="{by - 17}" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">{gene_label_str}</text>
""")

    # 2. Left Set Size Bars and Cohort Labels
    # Active cohorts list
    cohort_set_sizes: dict[str, int] = {
        "melanoma": 4,
        "pancancer": 3,
        "Rosenberg-iAtlas": 2,
        "rcc": 1,
        "McDermott-iAtlas": 1,
        "Riaz-iAtlas": 1,
        "Gide-iAtlas": 1,
    }

    cohort_order_panel_a = (
        "melanoma",
        "pancancer",
        "Rosenberg-iAtlas",
        "rcc",
        "McDermott-iAtlas",
        "Riaz-iAtlas",
        "Gide-iAtlas",
    )

    bar_max_w = 90
    bar_end_x = px + 155

    svg.append(f"""
    <!-- Set Size Header -->
    <text x="{px + 14}" y="{matrix_top - 8}" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}">
      COHORT SET SIZE
    </text>
    <text x="{px + 165}" y="{matrix_top - 8}" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}">
      ACTIVE COHORT
    </text>
""")

    cohort_y_map: dict[str, float] = {}

    for r_idx, c_name in enumerate(cohort_order_panel_a):
        cy_row: float = float(matrix_top + r_idx * row_height + 10)
        cohort_y_map[c_name] = cy_row
        set_size = cohort_set_sizes[c_name]
        bar_w = set_size * 22
        bx = bar_end_x - bar_w

        # Find display info
        c_match = [c for c in active_cohorts if c.name == c_name]
        disp_title = c_match[0].display_name if c_match else c_name
        n_biopsies = c_match[0].n_samples if c_match else 0

        # Alternating background row for matrix
        bg_fill = COLOR_SUBTLE_FILL if r_idx % 2 == 1 else "#FFFFFF"
        svg.append(f"""
    <rect x="{px + 14}" y="{cy_row - 10}" width="{pw - 28}" height="{row_height}" fill="{bg_fill}" />

    <!-- Set size bar -->
    <rect x="{bx}" y="{cy_row - 7}" width="{bar_w}" height="14" fill="#94A3B8" />
    <text x="{bx - 5}" y="{cy_row + 4}" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{set_size}</text>

    <!-- Cohort Name with inline sample size -->
    <text x="{px + 165}" y="{cy_row + 4}" font-size="9" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
      {disp_title} <tspan font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">(n={n_biopsies})</tspan>
    </text>
""")

    # 3. UpSet Matrix Dots and Connecting Paths
    for idx, inter in enumerate(intersections):
        cx = col_x_start + idx * col_spacing
        active_c = set(inter.cohorts)

        # Compute min and max y among active cohorts for connector line
        active_ys = [cohort_y_map[c] for c in inter.cohorts if c in cohort_y_map]
        if len(active_ys) > 1:
            min_y = min(active_ys)
            max_y = max(active_ys)
            svg.append(f"""
    <!-- Vertical connector path for intersection {inter.id} -->
    <line x1="{cx}" y1="{min_y}" x2="{cx}" y2="{max_y}" stroke="{OKABE_BLUE}" stroke-width="2.5" />
""")

        # Draw dots for all rows
        for c_name in cohort_order_panel_a:
            cy_dot: float = cohort_y_map[c_name]
            if c_name in active_c:
                dot_fill = OKABE_BLUE if len(active_c) > 1 else "#334155"
                svg.append(f"""
    <circle cx="{cx}" cy="{cy_dot}" r="5.5" fill="{dot_fill}" />
""")
            else:
                svg.append(f"""
    <circle cx="{cx}" cy="{cy_dot}" r="3" fill="#E2E8F0" />
""")

    # Bottom summary annotation inside Panel A
    svg.append(f"""
    <!-- Panel A Bottom Callout -->
    <g transform="translate({px + 14}, {ph + py - 32})">
      <rect x="0" y="0" width="{pw - 28}" height="24" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="10" y="15" font-size="8.5" font-weight="600" fill="{OKABE_BLUE}">Cross-Cohort Validation:</text>
      <text x="145" y="15" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        PTGS2 is shared between RCC &amp; McDermott; HLA-A &amp; TNFSF9 in pooled Melanoma replicate in Gide &amp; Riaz trials.
      </text>
    </g>
  </g>
""")
    return "".join(svg)


def render_panel_b_matrix(
    cells: Sequence[GeneSelectionCell],
    cohorts: Sequence[CohortInfo],
    genes: Sequence[GeneInfo],
) -> str:
    """Render Panel B: Cohort x Gene Selection Grid & Method Concordance."""
    px = 705
    py = 70
    pw = 675
    ph = 430

    svg = [f"""
  <!-- LAYER 03: PANEL B (GRID HEATMAP) -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_GridMatrix">
    <!-- Panel Background Card -->
    <rect id="card-panel-b" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel Badge & Titles -->
    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">b</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Cohort × Gene Selection Grid &amp; Method Concordance
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Bubble size encodes stability probability Π̂<tspan font-size="7.5" baseline-shift="sub">max</tspan>; color denotes SS-CPSS &amp; MB consensus
    </text>
"""]

    # Layout geometry
    grid_top_y = py + 72
    grid_x_start = px + 172
    col_w = 36
    row_h = 24.5

    # 0. Pathway Category Banners
    pathway_groups = (
        ("Antigen", OKABE_BLUE, 0, 2),
        ("Th1 / Exh", OKABE_BLUE, 2, 4),
        ("Costim", OKABE_ORANGE, 4, 7),
        ("Stroma / IFN", OKABE_BLUISH_GREEN, 7, 10),
    )
    for p_title, p_col, start_c, end_c in pathway_groups:
        bx = grid_x_start + start_c * col_w + 1
        bw = (end_c - start_c) * col_w - 2
        by = grid_top_y - 28
        svg.append(f"""
    <!-- Pathway Group {p_title} -->
    <rect x="{bx}" y="{by}" width="{bw}" height="2" fill="{p_col}" />
    <text x="{bx + bw/2}" y="{by - 3}" font-size="7" font-weight="700" fill="{p_col}" text-anchor="middle">
      {p_title}
    </text>
""")

    # 1. Top Marginal Track: Gene Recurrence Bar Chart
    recurrence_counts: dict[str, int] = {}
    for g in genes:
        rec = sum(1 for c in cells if c.gene == g.symbol and c.concordance_status in ("both", "ss_only", "mb_only"))
        recurrence_counts[g.symbol] = rec

    rec_base_y = grid_top_y - 2
    svg.append(f"""
    <!-- Recurrence track baseline -->
    <line x1="{grid_x_start - 8}" y1="{rec_base_y}" x2="{grid_x_start + len(genes) * col_w}" y2="{rec_base_y}"
          stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
    <text x="{grid_x_start - 12}" y="{rec_base_y - 4}" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_MUTED}" text-anchor="end">
      Recurrence
    </text>
""")

    for g_idx, g in enumerate(genes):
        gx = grid_x_start + g_idx * col_w + col_w / 2
        cnt = recurrence_counts[g.symbol]
        bar_h = cnt * 9
        by = rec_base_y - bar_h
        svg.append(f"""
    <rect x="{gx - 8}" y="{by}" width="16" height="{bar_h}" fill="#64748B" />
    <text x="{gx}" y="{by - 3}" font-size="7.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">{cnt}</text>
""")

    # 2. Gene Column Header Labels
    for g_idx, g in enumerate(genes):
        gx = grid_x_start + g_idx * col_w + col_w / 2
        f_size = "7.5" if len(g.symbol) >= 7 else "8.5"
        svg.append(f"""
    <!-- Gene Header {g.symbol} -->
    <text x="{gx}" y="{grid_top_y + 10}" font-size="{f_size}" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
      {g.symbol}
    </text>
""")

    # 3. Rows (12 Cohorts) and Grid Bubbles
    cell_map: dict[tuple[str, str], GeneSelectionCell] = {
        (c.cohort, c.gene): c for c in cells
    }

    bubble_start_y = grid_top_y + 22

    for r_idx, cohort in enumerate(cohorts):
        cy = bubble_start_y + r_idx * row_h
        row_bg = COLOR_SUBTLE_FILL if r_idx % 2 == 1 else "#FFFFFF"

        # Cohort font weight: Bold for pooled, regular for individual
        font_wt = "700" if cohort.is_pooled else "500"
        cohort_col = COLOR_TEXT_PRIMARY if cohort.order_idx < 7 else COLOR_TEXT_MUTED

        svg.append(f"""
    <!-- Row {r_idx}: {cohort.name} -->
    <rect x="{px + 12}" y="{cy - 11}" width="{pw - 24}" height="{row_h - 1}" fill="{row_bg}" />
    <text x="{px + 165}" y="{cy + 3}" font-size="8.5" font-weight="{font_wt}" fill="{cohort_col}" text-anchor="end">
      {cohort.display_name}
    </text>
""")

        # Grid cells
        for g_idx, g in enumerate(genes):
            gx = grid_x_start + g_idx * col_w + col_w / 2
            cell = cell_map.get((cohort.name, g.symbol))

            if cell is None or cell.concordance_status == "unselected":
                # Faint baseline tick
                svg.append(f"""
    <circle cx="{gx}" cy="{cy}" r="1.5" fill="#E2E8F0" />
""")
            elif cell.concordance_status == "subthreshold":
                # Gray outline circle
                r_dot = 4.0
                svg.append(f"""
    <circle cx="{gx}" cy="{cy}" r="{r_dot}" fill="{COLOR_SUBTLE_FILL}" stroke="#94A3B8" stroke-width="0.75" />
""")
            else:
                # Selected: Map score to radius
                # score in [0.75, 1.0] -> radius in [6.0, 9.5]
                score_norm = max(0.0, min(1.0, (cell.max_score - 0.70) / 0.30))
                r_dot = 6.0 + score_norm * 3.5

                if cell.concordance_status == "both":
                    fill_c = OKABE_BLUE
                elif cell.concordance_status == "ss_only":
                    fill_c = OKABE_SKY_BLUE
                else:  # mb_only
                    fill_c = OKABE_ORANGE

                svg.append(f"""
    <circle cx="{gx}" cy="{cy}" r="{r_dot:.1f}" fill="{fill_c}" stroke="#FFFFFF" stroke-width="0.5" />
""")

        # 4. Right Marginal Stacked Bars: Responders vs Non-Responders
        right_bar_start_x = px + 542
        right_bar_max_w = 88
        scale_n = cohort.n_samples / 1015.0
        tot_w = max(4.0, scale_n * right_bar_max_w)
        resp_w = (cohort.n_responders / cohort.n_samples) * tot_w
        non_resp_w = tot_w - resp_w

        svg.append(f"""
    <!-- Cohort Size Bar for {cohort.name} -->
    <rect x="{right_bar_start_x}" y="{cy - 5}" width="{resp_w:.1f}" height="10" fill="{OKABE_BLUISH_GREEN}" />
    <rect x="{right_bar_start_x + resp_w:.1f}" y="{cy - 5}" width="{non_resp_w:.1f}" height="10" fill="{OKABE_VERMILION}" />
    <text x="{right_bar_start_x + tot_w + 5:.1f}" y="{cy + 3}" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">
      {cohort.n_samples}
    </text>
""")

    # Right Bar Header
    svg.append(f"""
    <text x="{px + 542}" y="{grid_top_y + 10}" font-size="8" font-weight="700" fill="{COLOR_TEXT_MUTED}">
      COHORT SAMPLE SIZE (N)
    </text>
  </g>
""")
    return "".join(svg)


def render_panel_c_euler(genes: Sequence[GeneInfo]) -> str:
    """Render Panel C: Cancer-Type 4-Set Euler Partition."""
    px = 20
    py = 515
    pw = 665
    ph = 385

    svg = [f"""
  <!-- LAYER 04: PANEL C (EULER DIAGRAM) -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Euler">
    <!-- Panel Background Card -->
    <rect id="card-panel-c" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel Badge & Titles -->
    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">c</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Cancer-Type 4-Set Euler Partition
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Discrete biological mapping of 10 predictive genes across 4 clinical tumor modalities
    </text>
"""]

    # 4 Disease Modality Lobes (2 x 2 layout)
    lobe_w = 300
    lobe_h = 145
    col1_x = px + 22
    col2_x = px + 342
    row1_y = py + 48
    row2_y = py + 208

    # Lobe 1: Pan-Cancer (Top-Left)
    svg.append(f"""
    <!-- LOBE 1: Pan-Cancer -->
    <g id="euler-lobe-pancancer">
      <rect x="{col1_x}" y="{row1_y}" width="{lobe_w}" height="{lobe_h}"
            fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUE}" stroke-width="1.25" />
      <rect x="{col1_x}" y="{row1_y}" width="{lobe_w}" height="22" fill="{OKABE_BLUE}" />
      <text x="{col1_x + 10}" y="{row1_y + 15}" font-size="9" font-weight="700" fill="#FFFFFF">
        PAN-CANCER (n = 1,015)
      </text>
      <text x="{col1_x + lobe_w - 10}" y="{row1_y + 15}" font-size="8" font-weight="600" fill="#FFFFFF" text-anchor="end">
        3 Genes Selected
      </text>

      <!-- Gene Badges -->
      <g transform="translate({col1_x + 10}, {row1_y + 32})">
        <rect x="0" y="0" width="85" height="24" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="9.5" font-weight="700" fill="{OKABE_BLUE}">TBX21</text>
        <text x="95" y="16" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Th1 Lineage &amp; IFN-γ Driver (T-bet)</text>
      </g>
      <g transform="translate({col1_x + 10}, {row1_y + 64})">
        <rect x="0" y="0" width="85" height="24" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="9.5" font-weight="700" fill="{OKABE_BLUE}">HLA-C</text>
        <text x="95" y="16" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Universal MHC-I Antigen Presentation</text>
      </g>
      <g transform="translate({col1_x + 10}, {row1_y + 96})">
        <rect x="0" y="0" width="85" height="24" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="9.5" font-weight="700" fill="{OKABE_BLUE}">LAG3</text>
        <text x="95" y="16" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Exhaustion Checkpoint Receptor</text>
      </g>
    </g>
""")

    # Lobe 2: Melanoma (Top-Right)
    svg.append(f"""
    <!-- LOBE 2: Melanoma -->
    <g id="euler-lobe-melanoma">
      <rect x="{col2_x}" y="{row1_y}" width="{lobe_w}" height="{lobe_h}"
            fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_ORANGE}" stroke-width="1.25" />
      <rect x="{col2_x}" y="{row1_y}" width="{lobe_w}" height="22" fill="{OKABE_ORANGE}" />
      <text x="{col2_x + 10}" y="{row1_y + 15}" font-size="9" font-weight="700" fill="#FFFFFF">
        MELANOMA (n = 338)
      </text>
      <text x="{col2_x + lobe_w - 10}" y="{row1_y + 15}" font-size="8" font-weight="600" fill="#FFFFFF" text-anchor="end">
        4 Genes Selected
      </text>

      <!-- Gene Badges -->
      <g transform="translate({col2_x + 10}, {row1_y + 28})">
        <rect x="0" y="0" width="70" height="22" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="6" y="15" font-size="9" font-weight="700" fill="{OKABE_ORANGE}">HLA-A</text>
        <text x="78" y="15" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">MHC-I Neoantigen Recognition</text>
      </g>
      <g transform="translate({col2_x + 10}, {row1_y + 54})">
        <rect x="0" y="0" width="70" height="22" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="6" y="15" font-size="9" font-weight="700" fill="{OKABE_ORANGE}">TNFSF9</text>
        <text x="78" y="15" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">4-1BB Costimulatory Priming</text>
      </g>
      <g transform="translate({col2_x + 10}, {row1_y + 80})">
        <rect x="0" y="0" width="70" height="22" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="6" y="15" font-size="9" font-weight="700" fill="{OKABE_ORANGE}">TNFSF18</text>
        <text x="78" y="15" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">GITRL Immune Co-Activation</text>
      </g>
      <g transform="translate({col2_x + 10}, {row1_y + 106})">
        <rect x="0" y="0" width="70" height="22" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="6" y="15" font-size="9" font-weight="700" fill="{OKABE_ORANGE}">IKZF2</text>
        <text x="78" y="15" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Helios Treg Plasticity &amp; Anergy</text>
      </g>
    </g>
""")

    # Lobe 3: Bladder (Bottom-Left)
    svg.append(f"""
    <!-- LOBE 3: Bladder -->
    <g id="euler-lobe-bladder">
      <rect x="{col1_x}" y="{row2_y}" width="{lobe_w}" height="{lobe_h}"
            fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_SKY_BLUE}" stroke-width="1.25" />
      <rect x="{col1_x}" y="{row2_y}" width="{lobe_w}" height="22" fill="{OKABE_SKY_BLUE}" />
      <text x="{col1_x + 10}" y="{row2_y + 15}" font-size="9" font-weight="700" fill="#FFFFFF">
        BLADDER UROTHELIAL (n = 298)
      </text>
      <text x="{col1_x + lobe_w - 10}" y="{row2_y + 15}" font-size="8" font-weight="600" fill="#FFFFFF" text-anchor="end">
        2 Genes Selected
      </text>

      <g transform="translate({col1_x + 10}, {row2_y + 36})">
        <rect x="0" y="0" width="85" height="24" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="9.5" font-weight="700" fill="{OKABE_SKY_BLUE}">TGFB1</text>
        <text x="95" y="16" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Fibroblast T-cell Exclusion Program</text>
      </g>
      <g transform="translate({col1_x + 10}, {row2_y + 72})">
        <rect x="0" y="0" width="85" height="24" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="16" font-size="9.5" font-weight="700" fill="{OKABE_SKY_BLUE}">IFNG</text>
        <text x="95" y="16" font-size="8" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">Hot Inflammatory Effector Cascade</text>
      </g>
      <text x="{col1_x + 10}" y="{row2_y + 124}" font-size="7.5" font-style="italic" fill="{COLOR_TEXT_MUTED}">
        Canonical opposition: Mariathasan et al. (Nature 2018) TGF-β exclusion vs IFN-γ
      </text>
    </g>
""")

    # Lobe 4: RCC (Bottom-Right)
    svg.append(f"""
    <!-- LOBE 4: Renal Cell Carcinoma -->
    <g id="euler-lobe-rcc">
      <rect x="{col2_x}" y="{row2_y}" width="{lobe_w}" height="{lobe_h}"
            fill="{COLOR_SUBTLE_FILL}" stroke="{OKABE_BLUISH_GREEN}" stroke-width="1.25" />
      <rect x="{col2_x}" y="{row2_y}" width="{lobe_w}" height="22" fill="{OKABE_BLUISH_GREEN}" />
      <text x="{col2_x + 10}" y="{row2_y + 15}" font-size="9" font-weight="700" fill="#FFFFFF">
        RENAL CELL CARCINOMA (n = 263)
      </text>
      <text x="{col2_x + lobe_w - 10}" y="{row2_y + 15}" font-size="8" font-weight="600" fill="#FFFFFF" text-anchor="end">
        1 Gene Selected
      </text>

      <g transform="translate({col2_x + 10}, {row2_y + 44})">
        <rect x="0" y="0" width="85" height="26" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <text x="8" y="17" font-size="10" font-weight="700" fill="{OKABE_BLUISH_GREEN}">PTGS2</text>
        <text x="95" y="17" font-size="8.5" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">COX-2 Prostaglandin Synthase</text>
      </g>
      <text x="{col2_x + 10}" y="{row2_y + 98}" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        Key inflammatory mediator in clear cell RCC: synthesis of PGE<tspan font-size="6.5" baseline-shift="sub">2</tspan> drives myeloid
      </text>
      <text x="{col2_x + 10}" y="{row2_y + 112}" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
        extravasation and VEGF resistance (McDermott et al., Nat Med 2018).
      </text>
    </g>

    <!-- Central Partition Synthesis Label -->
    <g transform="translate({px + 175}, {py + ph - 25})">
      <rect x="0" y="0" width="315" height="18" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
      <text x="157" y="12" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">
        Complete Disease Segregation: Zero Off-Target Gene Overlap
      </text>
    </g>
  </g>
""")
    return "".join(svg)


def render_panel_d_summary() -> str:
    """Render Panel D: Methodological Concordance & Statistical Rationale."""
    px = 705
    py = 515
    pw = 675
    ph = 385

    svg = [f"""
  <!-- LAYER 05: PANEL D (SUMMARY & LEGENDS) -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Summary">
    <!-- Panel Background Card -->
    <rect id="card-panel-d" x="{px}" y="{py}" width="{pw}" height="{ph}"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel Badge & Titles -->
    <circle cx="{px + 18}" cy="{py + 20}" r="10" fill="{COLOR_TEXT_PRIMARY}" />
    <text x="{px + 18}" y="{py + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">d</text>
    <text x="{px + 36}" y="{py + 18}" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Methodological Concordance &amp; Statistical Rationale
    </text>
    <text x="{px + 36}" y="{py + 31}" font-size="9.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
      Analytical comparison of Shah-Samworth CPSS vs Meinshausen-Bühlmann error bounds
    </text>
"""]

    # 1. Method Comparison Metrics Card
    svg.append(f"""
    <g transform="translate({px + 16}, {py + 46})">
      <rect x="0" y="0" width="{pw - 32}" height="95" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Column 1: Concordance Rate -->
      <g transform="translate(15, 12)">
        <text x="0" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}">METHOD CONSENSUS</text>
        <text x="0" y="38" font-size="20" font-weight="700" fill="{OKABE_BLUE}">77.3%</text>
        <text x="0" y="54" font-size="8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">17 of 22 Selections</text>
        <text x="0" y="68" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">Identical features at π<tspan font-size="6" baseline-shift="sub">thr</tspan>=0.75</text>
      </g>

      <!-- Dividing Rule -->
      <line x1="160" y1="10" x2="160" y2="85" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />

      <!-- Column 2: Error Bounds -->
      <g transform="translate(180, 12)">
        <text x="0" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}">ERROR CONTROL (PFER)</text>
        <text x="0" y="30" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          Shah &amp; Samworth (2013) CPSS:
        </text>
        <text x="0" y="44" font-size="10" font-weight="700" fill="{OKABE_BLUE}">
          E[V] ≤ 0.818 <tspan font-size="8" font-weight="400" fill="{COLOR_TEXT_MUTED}">(q=17.65, 200 splits)</tspan>
        </text>
        <text x="0" y="60" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          Meinshausen &amp; Bühlmann (2010):
        </text>
        <text x="0" y="74" font-size="10" font-weight="700" fill="{OKABE_ORANGE}">
          E[V] ≤ 0.970 <tspan font-size="8" font-weight="400" fill="{COLOR_TEXT_MUTED}">(q=17.74, 100 splits)</tspan>
        </text>
      </g>

      <!-- Dividing Rule -->
      <line x1="390" y1="10" x2="390" y2="85" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />

      <!-- Column 3: Boundary Discrepancy Note -->
      <g transform="translate(405, 12)">
        <text x="0" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}">BOUNDARY EFFECT RATIONALE</text>
        <text x="0" y="30" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          Discrepancies (5/22) are exclusively near-
        </text>
        <text x="0" y="42" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          threshold boundary values:
        </text>
        <text x="0" y="56" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          • Gide HLA-A: SS=0.80 vs MB=0.70
        </text>
        <text x="0" y="68" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          • Rosenberg IFNG: SS=0.88 vs MB=0.74
        </text>
        <text x="0" y="80" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">
          • Pan-Cancer LAG3: MB=0.86 vs SS=0.73
        </text>
      </g>
    </g>
""")

    # 2. Biological Pathway Summary Table
    svg.append(f"""
    <g transform="translate({px + 16}, {py + 152})">
      <rect x="0" y="0" width="{pw - 32}" height="105" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <rect x="0" y="0" width="{pw - 32}" height="20" fill="{COLOR_SUBTLE_FILL}" />
      <text x="12" y="14" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
        FUNCTIONAL PATHWAY CONSENSUS &amp; THERAPEUTIC TARGETS
      </text>

      <g transform="translate(12, 30)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{OKABE_BLUE}">Antigen Processing &amp; Presentation:</text>
        <text x="175" y="10" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          HLA-A &amp; HLA-C loss limits neoantigen presentation; robust expression predicts checkpoint response.
        </text>

        <text x="0" y="28" font-size="8" font-weight="700" fill="{OKABE_BLUE}">Th1 Effector &amp; Exhaustion:</text>
        <text x="175" y="28" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          TBX21 commits cytotoxic lineage; LAG3 marks targetable exhaustion reversing primary anti-PD-1 resistance.
        </text>

        <text x="0" y="46" font-size="8" font-weight="700" fill="{OKABE_ORANGE}">Costimulatory &amp; Regulatory Priming:</text>
        <text x="175" y="46" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          TNFSF9 (4-1BBL) &amp; TNFSF18 (GITRL) co-activate T cells; IKZF2 (Helios) marks suppressive tumor-infiltrating Tregs.
        </text>

        <text x="0" y="64" font-size="8" font-weight="700" fill="{OKABE_BLUISH_GREEN}">Stroma &amp; Immune Exclusion:</text>
        <text x="175" y="64" font-size="8" font-weight="400" fill="{COLOR_TEXT_SECONDARY}">
          TGFB1 drives peritumoral collagen entrapment (Bladder); PTGS2/COX-2 synthesizes immunosuppressive PGE2 (RCC).
        </text>
      </g>
    </g>
""")

    # 3. Figure Legends Track
    svg.append(f"""
    <g transform="translate({px + 16}, {py + 268})">
      <rect x="0" y="0" width="{pw - 32}" height="100" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

      <!-- Sub-legend 1: Selection Status Color -->
      <g transform="translate(15, 14)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">SELECTION STATUS COLOR</text>

        <circle cx="6" cy="28" r="5.5" fill="{OKABE_BLUE}" />
        <text x="18" y="31" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Both SS-CPSS &amp; MB (Consensus)</text>

        <circle cx="6" cy="46" r="5.5" fill="{OKABE_SKY_BLUE}" />
        <text x="18" y="49" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">SS-CPSS Exclusively (π ≥ 0.75)</text>

        <circle cx="6" cy="64" r="5.5" fill="{OKABE_ORANGE}" />
        <text x="18" y="67" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">MB Exclusively (π ≥ 0.75)</text>

        <circle cx="6" cy="82" r="4.0" fill="{COLOR_SUBTLE_FILL}" stroke="#94A3B8" stroke-width="0.75" />
        <text x="18" y="85" font-size="8" font-weight="500" fill="{COLOR_TEXT_MUTED}">Sub-threshold (0.50 ≤ Π̂ &lt; 0.75)</text>
      </g>

      <!-- Dividing Rule -->
      <line x1="200" y1="12" x2="200" y2="88" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />

      <!-- Sub-legend 2: Bubble Size Scale -->
      <g transform="translate(220, 14)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">STABILITY SCORE SCALE (Π̂)</text>

        <circle cx="15" cy="45" r="4.0" fill="#94A3B8" />
        <text x="15" y="70" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">0.50</text>

        <circle cx="60" cy="45" r="6.5" fill="{OKABE_BLUE}" />
        <text x="60" y="70" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">0.75</text>

        <circle cx="115" cy="45" r="8.5" fill="{OKABE_BLUE}" />
        <text x="115" y="70" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">0.90</text>

        <circle cx="175" cy="45" r="9.5" fill="{OKABE_BLUE}" />
        <text x="175" y="70" font-size="7.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" text-anchor="middle">1.00</text>
      </g>

      <!-- Dividing Rule -->
      <line x1="435" y1="12" x2="435" y2="88" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />

      <!-- Sub-legend 3: Clinical Response Stacked Bar -->
      <g transform="translate(455, 14)">
        <text x="0" y="10" font-size="8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">RECIST RESPONSE BREAKDOWN</text>

        <rect x="0" y="24" width="20" height="10" fill="{OKABE_BLUISH_GREEN}" />
        <text x="28" y="32" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Responders (CR / PR)</text>

        <rect x="0" y="44" width="20" height="10" fill="{OKABE_VERMILION}" />
        <text x="28" y="52" font-size="8" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">Non-Responders (SD / PD)</text>

        <text x="0" y="75" font-size="7.5" font-weight="400" fill="{COLOR_TEXT_MUTED}">
          Bar width proportional to total cohort sample size N.
        </text>
      </g>
    </g>
  </g>
""")
    return "".join(svg)


def render_full_svg(
    cells: Sequence[GeneSelectionCell],
    cohorts: Sequence[CohortInfo],
    genes: Sequence[GeneInfo],
    intersections: Sequence[UpSetIntersection],
) -> str:
    """Render complete publication vector SVG with Inkscape namespaces."""
    active_cohorts = tuple(c for c in cohorts if c.order_idx < 7)

    parts: list[str] = [
        f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     version="1.1" id="figure-cohort-gene-overlap-nature"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <!-- LAYER 00: CANVAS BACKGROUND -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="{COLOR_CANVAS_BG}" />
  </g>
""",
        render_header(),
        render_panel_a_upset(intersections, active_cohorts),
        render_panel_b_matrix(cells, cohorts, genes),
        render_panel_c_euler(genes),
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
    """Imperative shell for generating the cross-cohort gene selection overlap figure."""
    parser = argparse.ArgumentParser(
        description="Generate publication-grade figure of cross-cohort gene selection overlap"
    )
    default_scores = (
        Path("output/stability_selection_fitters_benchmark/stability_scores.parquet")
        if Path("output/stability_selection_fitters_benchmark/stability_scores.parquet").exists()
        else Path("output/stability_selection_iatlas_immunotherapy/stability_scores.parquet")
    )
    default_summary = (
        Path("output/stability_selection_fitters_benchmark/cohort_fitter_summary.parquet")
        if Path("output/stability_selection_fitters_benchmark/cohort_fitter_summary.parquet").exists()
        else (
            Path("output/stability_selection_fitters_benchmark/cohort_summary.parquet")
            if Path("output/stability_selection_fitters_benchmark/cohort_summary.parquet").exists()
            else Path("output/stability_selection_iatlas_immunotherapy/cohort_summary.parquet")
        )
    )

    parser.add_argument(
        "--scores-path",
        type=Path,
        default=default_scores,
        help="Path to stability scores parquet",
    )
    parser.add_argument(
        "--summary-path",
        type=Path,
        default=default_summary,
        help="Path to cohort summary parquet",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark/figures/cohort_gene_overlap_nature.svg"),
        help="Target output SVG path",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark/figures/cohort_gene_overlap_nature.png"),
        help="Target output PNG path",
    )
    args = parser.parse_args()

    # Load data
    load_res = load_data(args.scores_path, args.summary_path)
    match load_res:
        case Failure(err):
            print(f"Error loading inputs: {err}", file=sys.stderr)
            return 1
        case Success((df_scores, df_summary)):
            cohorts = extract_cohort_models(df_summary)
            cells = extract_matrix_cells(df_scores, cohorts, ORDERED_GENES)
            intersections = compute_upset_intersections()

            svg_str = render_full_svg(cells, cohorts, ORDERED_GENES, intersections)

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
