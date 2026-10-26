#!/usr/bin/env python3
"""
Cross-Cohort Deconvolution Scatter Grid (Panel A Composite).

Visualizes the baseline empirical concordance between single-cell graph differential
abundance (Milopy Mean logFC) and patient-level pseudobulk deconvolution effect (beta)
across all nine single-cell immunotherapy response cohorts.

Design Specifications:
- 3x3 multi-panel composite grid (Panels A through I).
- Evaluates 5 independent baseline sweeps per cohort to demonstrate sampling stability.
- Strictly single uniform mark color (Okabe-Ito blue #1F77B4 with 1px black outline) with NO separate quadrant coloring.
- Strictly NO in-plot text annotations over data points (clean vector marks, cell states available via SVG hover tooltips).
- Minimalist Nature Methods wireframe styling: white canvas, hairline borders, dashed regression lines, and zero crosshairs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import altair as alt
import polars as pl
import vl_convert as vlc
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# ------------------------------------------------------------------------------
# Cohort Metadata & Ordering Matrix
# ------------------------------------------------------------------------------

COHORT_ORDER_AND_METADATA: Final[tuple[dict[str, str], ...]] = (
    # Row 1: Melanoma, Pan-Cancer, and Colorectal (Single panel letter 'c' for the whole 3x3 grid)
    {"cohort": "GSE120575", "letter": "c", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018"},
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


# ------------------------------------------------------------------------------
# Configuration & Data Models
# ------------------------------------------------------------------------------

# 130 mm at 96 CSS px/inch (25.4 mm/inch) = 491.34 px; 10 pt = 10 * (96/72) = 13.33 px
TARGET_WIDTH_MM: Final[float] = 130.0
TARGET_WIDTH_PX: Final[float] = TARGET_WIDTH_MM * 96.0 / 25.4  # 491.34 px
MIN_FONT_PT: Final[float] = 10.0
MIN_FONT_PX: Final[float] = 13.5  # 10.125 pt at 96 DPI (strictly >= 10pt)


class PlotConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    benchmark_dir: Path = Path(
        "output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark"
    )
    out_svg: Path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
    )
    out_png: Path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.png"
    )
    article_png: Path = Path("article/figures/milo-deconv/nine_cohorts_deconv_scatter_grid.png")
    panel_width: int = 118
    panel_height: int = 78
    point_size: int = 42
    point_color: str = "#1F77B4"  # Single uniform Okabe-Ito blue (no quadrant coloring)


class CohortPanelData(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    cohort_id: str
    panel_letter: str
    cancer_type: str
    author: str
    year: str
    points_df: pl.DataFrame
    rho_mean: float
    rho_std: float
    sign_mean: float
    sign_std: float
    n_states: int
    n_sweeps: int


# ------------------------------------------------------------------------------
# Pure Functional Data Ingestion
# ------------------------------------------------------------------------------

def load_cohort_panel_data(
    benchmark_dir: Path,
    meta: dict[str, str],
) -> Result[CohortPanelData, str]:
    """Pure function loading 5 baseline sweeps points and summary statistics for a given cohort."""
    cohort_id = meta["cohort"]
    cohort_dir = benchmark_dir / cohort_id

    pts_file = cohort_dir / f"{cohort_id}_baseline_sweeps_points.parquet"
    sum_file = cohort_dir / f"{cohort_id}_baseline_sweeps_summary.parquet"

    # Fallback to single self-concordance if sweeps parquet not found
    if not pts_file.exists():
        pts_file = cohort_dir / f"{cohort_id}_self_concordance.parquet"
    if not sum_file.exists():
        sum_file = cohort_dir / f"{cohort_id}_concordance_summary.parquet"

    if not pts_file.exists():
        return Failure(f"Points parquet not found for {cohort_id} at {pts_file}")
    if not sum_file.exists():
        return Failure(f"Summary parquet not found for {cohort_id} at {sum_file}")

    try:
        raw_pts = pl.read_parquet(pts_file)
        raw_sum = pl.read_parquet(sum_file)

        # Filter to Combined condition
        pts = raw_pts.filter(pl.col("condition") == "Combined")
        sums = raw_sum.filter(pl.col("condition") == "Combined")

        if pts.height == 0:
            return Failure(f"No points found under 'Combined' condition for {cohort_id}")

        n_sweeps = int(pts["sweep_idx"].n_unique()) if "sweep_idx" in pts.columns else 1
        n_states = int(pts["cell_state"].n_unique())

        if "spearman_rho_milo_deconv" in sums.columns:
            rho_col = "spearman_rho_milo_deconv"
        elif "spearman_rho" in sums.columns:
            rho_col = "spearman_rho"
        else:
            return Failure(f"Missing spearman rho column in summary for {cohort_id}")

        rho_mean_val = sums[rho_col].mean()
        rho_std_val = sums[rho_col].std() if sums.height > 1 else 0.0
        sign_mean_val = sums["sign_concordance_pct"].mean()
        sign_std_val = sums["sign_concordance_pct"].std() if sums.height > 1 else 0.0

        rho_mean = float(rho_mean_val) if rho_mean_val is not None else 0.0
        rho_std = float(rho_std_val) if rho_std_val is not None and not (rho_std_val != rho_std_val) else 0.0
        sign_mean = float(sign_mean_val) if sign_mean_val is not None else 0.0
        sign_std = float(sign_std_val) if sign_std_val is not None and not (sign_std_val != sign_std_val) else 0.0

        return Success(
            CohortPanelData(
                cohort_id=cohort_id,
                panel_letter=meta["letter"],
                cancer_type=meta["cancer_type"],
                author=meta["author"],
                year=meta["year"],
                points_df=pts,
                rho_mean=rho_mean,
                rho_std=rho_std,
                sign_mean=sign_mean,
                sign_std=sign_std,
                n_states=n_states,
                n_sweeps=n_sweeps,
            )
        )
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to load concordance data for {cohort_id}: {exc}")


# ------------------------------------------------------------------------------
# Pure Functional Altair Chart Construction
# ------------------------------------------------------------------------------

def build_cohort_scatter_panel(
    panel_data: CohortPanelData,
    config: PlotConfig,
) -> alt.Chart:
    """
    Construct a single Panel A scatter chart sized for a 130mm 3x3 grid with >=10pt typography:
    - Single uniform mark color (#1F77B4) with black outline.
    - Zero separate quadrant coloring.
    - Zero text annotations over data marks.
    - Dashed linear regression trend line and zero reference crosshairs.
    """
    df_pd = panel_data.points_df.to_pandas()

    is_left_col = panel_data.cohort_id in {"GSE120575", "CELLxGENE_6f9de485", "GSE233203"}
    is_bottom_row = panel_data.cohort_id in {"GSE233203", "GSE200996", "GSE316195"}

    x_title = "Milopy logFC" if is_bottom_row else None
    y_title = "Deconv. β" if is_left_col else None

    c_base = (
        alt.Chart(df_pd)
        .mark_circle(
            size=config.point_size,
            color=config.point_color,
            opacity=0.85,
            stroke="black",
            strokeWidth=0.75,
        )
        .encode(
            x=alt.X(
                "milo_mean_logfc:Q",
                title=x_title,
                axis=alt.Axis(
                    grid=True,
                    gridColor="#EEEEEE",
                    tickCount=3,
                    labelFontSize=MIN_FONT_PX,
                    titleFontSize=MIN_FONT_PX,
                    titlePadding=4,
                    labelPadding=2,
                ),
            ),
            y=alt.Y(
                "deconv_beta:Q",
                title=y_title,
                axis=alt.Axis(
                    grid=True,
                    gridColor="#EEEEEE",
                    tickCount=3,
                    labelFontSize=MIN_FONT_PX,
                    titleFontSize=MIN_FONT_PX,
                    titlePadding=4,
                    labelPadding=2,
                ),
            ),
            tooltip=[
                alt.Tooltip("cell_state:N", title="Cell State"),
                alt.Tooltip("sweep_idx:O", title="Sweep Replicate"),
                alt.Tooltip("milo_mean_logfc:Q", title="Milopy logFC", format="+.3f"),
                alt.Tooltip("deconv_beta:Q", title="Deconvolution β", format="+.3f"),
            ],
        )
    )

    # Dashed linear regression trend line
    trend_line = (
        c_base.transform_regression("milo_mean_logfc", "deconv_beta")
        .mark_line(color="#444444", strokeDash=[3, 3], strokeWidth=1.5)
    )

    # Reference crosshairs at x=0 and y=0
    h_rule = alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(y=alt.datum(0.0))
    v_rule = alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(x=alt.datum(0.0))

    # Multi-line compact title and subtitle so >=10pt text fits cleanly within a 43mm column
    if panel_data.panel_letter:
        title_lines = [
            f"({panel_data.panel_letter})  {panel_data.cancer_type} ({panel_data.year})",
            f"{panel_data.author}",
        ]
    else:
        title_lines = [
            f"{panel_data.cancer_type} ({panel_data.year})",
            f"{panel_data.author}",
        ]
    title_params = alt.TitleParams(
        text=title_lines,
        fontSize=14.0,
        anchor="start",
        offset=4,
    )

    panel_chart = (
        (h_rule + v_rule + trend_line + c_base)
        .properties(
            width=config.panel_width,
            height=config.panel_height,
            title=title_params,
        )
    )

    return panel_chart


def assemble_3x3_scatter_grid(
    panels: list[alt.Chart],
) -> alt.VConcatChart:
    """Assemble 9 panels into a 3x3 composite grid sized for 130mm width and >=10pt minimum text."""
    if len(panels) != 9:
        raise ValueError(f"Expected exactly 9 panels for a 3x3 grid, received {len(panels)}")

    row1 = alt.hconcat(panels[0], panels[1], panels[2], spacing=12)
    row2 = alt.hconcat(panels[3], panels[4], panels[5], spacing=12)
    row3 = alt.hconcat(panels[6], panels[7], panels[8], spacing=12)

    grid = (
        alt.vconcat(row1, row2, row3, spacing=8)
        .properties(
            padding={"left": 4, "right": 4, "top": 4, "bottom": 4},
        )
        .configure_view(strokeWidth=0.75, stroke="#CBD5E1")
        .configure_axis(labelFontSize=MIN_FONT_PX, titleFontSize=MIN_FONT_PX)
        .configure_title(fontSize=14.0, subtitleFontSize=MIN_FONT_PX)
    )

    return grid


def enforce_svg_physical_dimensions_mm(svg_text: str, width_mm: float = TARGET_WIDTH_MM) -> str:
    """Set exact physical width in mm (130mm) on the root <svg> element while preserving aspect ratio."""
    import re

    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    if not m:
        return svg_text
    vb_w = float(m.group(1))
    vb_h = float(m.group(2))
    height_mm = width_mm * (vb_h / max(vb_w, 1.0))

    svg_updated = re.sub(
        r'(<svg[^>]*?)\s+width="[^"]+"\s+height="[^"]+"',
        rf'\1 width="{width_mm:.1f}mm" height="{height_mm:.2f}mm"',
        svg_text,
        count=1,
    )
    return svg_updated


# ------------------------------------------------------------------------------
# Imperative Shell: Orchestration & File Exports
# ------------------------------------------------------------------------------

def generate_nine_cohorts_scatter_grid(config: PlotConfig) -> Result[Path, str]:
    """Execute complete ingestion, rendering, and vector/raster export pipeline."""
    config.out_svg.parent.mkdir(parents=True, exist_ok=True)
    config.article_png.parent.mkdir(parents=True, exist_ok=True)

    print("=== Constructing 3x3 Baseline Deconvolution Concordance Grid across 9 Cohorts ===")
    panels: list[alt.Chart] = []

    for meta in COHORT_ORDER_AND_METADATA:
        match load_cohort_panel_data(config.benchmark_dir, meta):
            case Failure(err):
                return Failure(f"Failed to load panel {meta['letter']} ({meta['cohort']}): {err}")
            case Success(panel_data):
                print(
                    f"  Panel {panel_data.panel_letter}: {panel_data.cohort_id:18s} | "
                    f"Points={panel_data.points_df.height:2d} ({panel_data.n_sweeps} sweeps x {panel_data.n_states} states) | "
                    f"rho={panel_data.rho_mean:+.3f} ± {panel_data.rho_std:.3f} | "
                    f"sign={panel_data.sign_mean:.0f}% ± {panel_data.sign_std:.0f}%"
                )
                chart = build_cohort_scatter_panel(panel_data, config)
                panels.append(chart)

    full_grid = assemble_3x3_scatter_grid(panels)
    grid_dict = full_grid.to_dict()

    # 1. Export SVG with exact 130mm physical width
    print(f"\nExporting publication vector SVG (width={TARGET_WIDTH_MM:.0f}mm, min font >= {MIN_FONT_PT:.0f}pt) to {config.out_svg}...")
    raw_svg = vlc.vegalite_to_svg(grid_dict)
    svg_str = enforce_svg_physical_dimensions_mm(raw_svg, width_mm=TARGET_WIDTH_MM)
    with open(config.out_svg, "w", encoding="utf-8") as f:
        f.write(svg_str)

    # 2. Export 300 DPI PNG at 130mm width (5.118 in * 300 DPI = 1535 px)
    print(f"Exporting 300 DPI raster PNG to {config.out_png}...")
    png_bytes = vlc.vegalite_to_png(grid_dict, scale=3.125)
    with open(config.out_png, "wb") as f:
        f.write(png_bytes)

    # 3. Copy to manuscript figures directory
    print(f"Copying figure to manuscript directory {config.article_png}...")
    with open(config.article_png, "wb") as f:
        f.write(png_bytes)

    print("\nSuccessfully exported all vector and raster assets!")
    return Success(config.out_svg)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 3x3 composite scatter grid of Milopy DA vs Pseudobulk Deconvolution across 9 cohorts."
    )
    parser.add_argument(
        "--benchmark-dir",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark"
        ),
        help="Root benchmark output directory containing cohort subfolders",
    )
    parser.add_argument(
        "--out-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
        ),
        help="Target output SVG path",
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.png"
        ),
        help="Target output PNG path",
    )
    parser.add_argument(
        "--article-png",
        type=Path,
        default=Path("article/figures/milo-deconv/nine_cohorts_deconv_scatter_grid.png"),
        help="Target manuscript figures PNG path",
    )
    args = parser.parse_args()

    config = PlotConfig(
        benchmark_dir=args.benchmark_dir,
        out_svg=args.out_svg,
        out_png=args.out_png,
        article_png=args.article_png,
    )

    match generate_nine_cohorts_scatter_grid(config):
        case Failure(err):
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"Complete: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
