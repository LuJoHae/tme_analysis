#!/usr/bin/env python3
"""
Publication-Grade 3x3 Single-Cell Landscape UMAP Grid across 9 Clinical Cohorts.
Generates Panel a of the master 2x2 immunotherapy figure.
Adheres strictly to Nature Methods minimalist wireframe standard, Okabe-Ito palette,
and >=10pt (13.5px) minimum typography.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Final

import altair as alt
import numpy as np
import pandas as pd
import polars as pl
import vl_convert as vlc
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# Disable Altair row limit
alt.data_transformers.disable_max_rows()

# ------------------------------------------------------------------------------
# Typography & Palette Constants
# ------------------------------------------------------------------------------
MIN_FONT_PT: Final[float] = 10.0
MIN_FONT_PX: Final[float] = 13.5  # 10 pt * (96 / 72) = 13.33 -> 13.5 px
TITLE_FONT_PX: Final[float] = 14.0

TARGET_WIDTH_MM: Final[float] = 130.0
TARGET_HEIGHT_MM: Final[float] = 123.31

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_VERMILION: Final[str] = "#D55E00"
COLOR_NEUTRAL_GREY: Final[str] = "#CBD5E1"
COLOR_HAIRLINE: Final[str] = "#CBD5E1"

# ------------------------------------------------------------------------------
# Cohort Metadata & Ordering Matrix (Exact match with deconvolution Panel c)
# ------------------------------------------------------------------------------
COHORT_ORDER_AND_METADATA: Final[tuple[dict[str, str], ...]] = (
    # Row 1: Melanoma, Pan-Cancer, and Colorectal
    {"cohort": "GSE120575", "letter": "a", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018"},
    {"cohort": "CELLxGENE_7b20c613", "letter": "", "cancer_type": "Pan-Cancer", "author": "Gondal et al.", "year": "2025"},
    {"cohort": "CELLxGENE_05a8c945", "letter": "", "cancer_type": "Colorectal", "author": "Marteau et al.", "year": "2026"},
    # Row 2: Breast and Non-Small Cell Lung Cancer
    {"cohort": "CELLxGENE_6f9de485", "letter": "", "cancer_type": "Breast", "author": "Yan et al.", "year": "2026"},
    {"cohort": "GSE207422", "letter": "", "cancer_type": "NSCLC", "author": "Hu et al.", "year": "2023"},
    {"cohort": "GSE243013", "letter": "", "cancer_type": "NSCLC", "author": "Liu et al.", "year": "2025"},
    # Row 3: NSCLC, Head & Neck, and Pancreatic
    {"cohort": "GSE233203", "letter": "", "cancer_type": "NSCLC", "author": "Hong et al.", "year": "2025"},
    {"cohort": "GSE200996", "letter": "", "cancer_type": "HNSCC", "author": "Luoma et al.", "year": "2022"},
    {"cohort": "GSE316195", "letter": "", "cancer_type": "PDAC", "author": "Raufi et al.", "year": "2026"},
)


class UmapGridConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    results_dir: Path = Path("output/results/milopy_response_analysis")
    out_svg: Path = Path("results/milopy_analysis/nine_cohorts_umaps_logfc_grid.svg")
    out_png: Path = Path("results/milopy_analysis/nine_cohorts_umaps_logfc_grid.png")
    article_png: Path = Path("article/figures/single-cell/nine_cohorts_umaps_logfc_grid.png")
    panel_width: int = 142
    panel_height: int = 79
    max_cells_per_panel: int = 1800
    point_size: int = 14
    seed: int = 42


def load_and_subsample_cells(
    cohort_dir: Path,
    config: UmapGridConfig,
) -> Result[pd.DataFrame, str]:
    """Loads cell scores and intelligently subsamples while preserving significant cells."""
    cell_file = cohort_dir / "cell_level_scores.parquet"
    if not cell_file.exists():
        return Failure(f"Missing {cell_file}")

    try:
        df = pl.read_parquet(cell_file)
        if "UMAP1" not in df.columns or "UMAP2" not in df.columns or "cell_logfc" not in df.columns:
            return Failure(f"Missing essential columns in {cell_file}")

        pdf = df.select(["UMAP1", "UMAP2", "cell_logfc", "da_status", "cell_type"]).to_pandas()
        pdf["abs_logfc"] = np.abs(pdf["cell_logfc"])
        pdf = pdf.sort_values(by="abs_logfc", ascending=True)

        # Smart sampling: guarantee significant cells are retained
        target = config.max_cells_per_panel
        if len(pdf) > target:
            if "da_status" in pdf.columns:
                sig_mask = pdf["da_status"] != "Not Significant"
            else:
                sig_mask = pdf["abs_logfc"] > 0.5

            sig_pdf = pdf[sig_mask]
            non_sig_pdf = pdf[~sig_mask]

            if len(sig_pdf) > 0 and len(non_sig_pdf) > 0:
                n_sig = min(len(sig_pdf), int(target * 0.65))
                n_nonsig = target - n_sig
                sample_sig = sig_pdf.sample(n=n_sig, random_state=config.seed)
                sample_nonsig = non_sig_pdf.sample(n=n_nonsig, random_state=config.seed)
                pdf = pd.concat([sample_nonsig, sample_sig], axis=0).sort_values(by="abs_logfc", ascending=True)
            else:
                pdf = pdf.sample(n=target, random_state=config.seed).sort_values(by="abs_logfc", ascending=True)

        return Success(pdf)
    except Exception as exc:
        return Failure(f"Error loading {cell_file}: {exc}")


def build_umap_single_panel(
    pdf: pd.DataFrame,
    meta: dict[str, str],
    config: UmapGridConfig,
    show_x_title: bool,
    show_y_title: bool,
) -> alt.Chart:
    """Builds a single UMAP panel with >=10pt typography and clean coordinate frame."""
    x_title = "UMAP 1" if show_x_title else None
    y_title = "UMAP 2" if show_y_title else None

    # Title lines: letter prefix on the top-left panel only
    letter_prefix = f"({meta['letter']})  " if meta["letter"] else ""
    title_lines = [
        f"{letter_prefix}{meta['cancer_type']} ({meta['year']})",
        f"{meta['author']}",
    ]

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=config.point_size, opacity=0.85, stroke=COLOR_HAIRLINE, strokeWidth=0.5)
        .encode(
            x=alt.X(
                "UMAP1:Q",
                title=x_title,
                axis=alt.Axis(
                    grid=False,
                    ticks=False,
                    labels=False,
                    domain=True,
                    domainColor=COLOR_HAIRLINE,
                    domainWidth=0.75,
                    titleFontSize=MIN_FONT_PX,
                    titleFont="sans-serif",
                    titlePadding=4,
                ),
            ),
            y=alt.Y(
                "UMAP2:Q",
                title=y_title,
                axis=alt.Axis(
                    grid=False,
                    ticks=False,
                    labels=False,
                    domain=True,
                    domainColor=COLOR_HAIRLINE,
                    domainWidth=0.75,
                    titleFontSize=MIN_FONT_PX,
                    titleFont="sans-serif",
                    titlePadding=4,
                ),
            ),
            color=alt.Color(
                "cell_logfc:Q",
                scale=alt.Scale(
                    domain=[-5.0, 0.0, 5.0],
                    range=[OKABE_BLUE, "#FFFFFF", OKABE_VERMILION],
                    interpolate="rgb",
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    direction="horizontal",
                    title="Milo log2FC (Responder vs Non-Responder; clamped at ±5)",
                    titleFontSize=MIN_FONT_PX,
                    labelFontSize=MIN_FONT_PX,
                    titleLimit=0,
                    gradientLength=260,
                    values=[-5.0, -2.5, 0.0, 2.5, 5.0],
                ),
            ),
        )
        .properties(
            width=config.panel_width,
            height=config.panel_height,
            title=alt.TitleParams(
                text=title_lines,
                fontSize=TITLE_FONT_PX,
                font="sans-serif",
                anchor="start",
                offset=4,
            ),
        )
    )
    return chart


def assemble_3x3_umap_grid(panels: list[alt.Chart]) -> alt.VConcatChart:
    """Assembles 9 UMAP panels into a 3x3 composite grid with shared colorbar."""
    row1 = alt.hconcat(panels[0], panels[1], panels[2], spacing=12)
    row2 = alt.hconcat(panels[3], panels[4], panels[5], spacing=12)
    row3 = alt.hconcat(panels[6], panels[7], panels[8], spacing=12)

    grid = (
        alt.vconcat(row1, row2, row3, spacing=8)
        .resolve_scale(color="shared")
        .properties(
            padding={"left": 4, "right": 4, "top": 4, "bottom": 5},
        )
        .configure_view(strokeWidth=0.75, stroke=COLOR_HAIRLINE)
    )
    return grid


def enforce_svg_physical_dimensions_mm(svg_text: str, width_mm: float = 130.0, height_mm: float = 154.07) -> str:
    """Enforces exact mm width and height on SVG tag while keeping viewBox intact."""
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    if not m:
        return svg_text

    svg_updated = re.sub(
        r'(<svg[^>]*?)\s+width="[^"]+"\s+height="[^"]+"',
        rf'\1 width="{width_mm:.1f}mm" height="{height_mm:.2f}mm"',
        svg_text,
        count=1,
    )
    return svg_updated


def generate_nine_cohorts_umaps_grid(config: UmapGridConfig) -> Result[Path, str]:
    """Generates Panel a 3x3 UMAP grid vector SVG and 300 DPI raster PNG."""
    config.out_svg.parent.mkdir(parents=True, exist_ok=True)
    config.article_png.parent.mkdir(parents=True, exist_ok=True)

    print("=== Constructing 3x3 Single-Cell Landscape UMAP Grid across 9 Cohorts ===")
    panels: list[alt.Chart] = []

    for idx, meta in enumerate(COHORT_ORDER_AND_METADATA):
        cohort_id = meta["cohort"]
        cohort_dir = config.results_dir / cohort_id
        match load_and_subsample_cells(cohort_dir, config):
            case Failure(err):
                return Failure(f"Failed to load cells for {cohort_id}: {err}")
            case Success(pdf):
                show_x = (idx >= 6)  # Bottom row gets X title
                show_y = (idx % 3 == 0)  # Left column gets Y title
                chart = build_umap_single_panel(pdf, meta, config, show_x, show_y)
                panels.append(chart)
                print(f"  Panel {meta['letter'] or ' '}: {cohort_id:18s} | Cells={len(pdf):,}")

    grid = assemble_3x3_umap_grid(panels)
    grid_dict = grid.to_dict()

    print(f"\nExporting publication vector SVG (width={TARGET_WIDTH_MM:.0f}mm, height={TARGET_HEIGHT_MM:.2f}mm) to {config.out_svg}...")
    raw_svg = vlc.vegalite_to_svg(grid_dict)
    svg_str = enforce_svg_physical_dimensions_mm(raw_svg, width_mm=TARGET_WIDTH_MM, height_mm=TARGET_HEIGHT_MM)

    with open(config.out_svg, "w", encoding="utf-8") as f:
        f.write(svg_str)

    print(f"Exporting 300 DPI raster PNG to {config.out_png}...")
    png_bytes = vlc.vegalite_to_png(grid_dict, scale=3.125)
    with open(config.out_png, "wb") as f:
        f.write(png_bytes)

    with open(config.article_png, "wb") as f:
        f.write(png_bytes)

    print("\nSuccessfully exported Panel a UMAP grid assets!")
    return Success(config.out_svg)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot 3x3 composite UMAP grid across 9 cohorts.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/results/milopy_response_analysis"),
        help="Root milopy analysis directory",
    )
    parser.add_argument(
        "--out-svg",
        type=Path,
        default=Path("results/milopy_analysis/nine_cohorts_umaps_logfc_grid.svg"),
        help="Output SVG path",
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=Path("results/milopy_analysis/nine_cohorts_umaps_logfc_grid.png"),
        help="Output PNG path",
    )
    args = parser.parse_args()

    config = UmapGridConfig(
        results_dir=args.results_dir,
        out_svg=args.out_svg,
        out_png=args.out_png,
    )

    match generate_nine_cohorts_umaps_grid(config):
        case Failure(err):
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"Complete: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
