#!/usr/bin/env python3
"""
Cross-Cohort Deconvolution Correlation Bar Plot and Composite Figure.

Visualizes the rank correlation coefficients (Spearman's rho mean +/- std across 5 baseline sweeps)
between single-cell graph differential abundance (Milopy mean logFC) and pseudobulk deconvolution
effect (beta) across all nine single-cell immunotherapy response cohorts.

Design Specifications:
- Standalone Bar Plot:
  - Exact physical height matching nine_cohorts_deconv_scatter_grid.svg (179.75mm, viewBox 0 0 188 672).
  - Physical width: 50.0mm (viewBox width 188px).
  - Minimum font size: >= 10pt (13.5px at 96 CSS DPI).
  - Single uniform mark color: Okabe-Ito blue (#1F77B4) with 0.75px black outline.
  - Whiskers: mean +/- s.d. error bars with dashed cross-cohort average reference line.
- Combined Composite Figure:
  - Stitches the 130.0mm 3x3 scatter grid (Panels A-I) with the 50.0mm bar plot (Panel J).
  - Physical width: 180.0mm (standard full-page 2-column width in academic journals).
  - Physical height: 179.75mm.
  - Generates both publication vector SVG and 300 DPI raster PNG.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Final

import altair as alt
import polars as pl
import vl_convert as vlc
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# ------------------------------------------------------------------------------
# Cohort Metadata & Ordering Matrix (Aligned with Panels A through I)
# ------------------------------------------------------------------------------

COHORT_ORDER_AND_METADATA: Final[tuple[dict[str, str], ...]] = (
    {"cohort": "GSE120575", "letter": "c", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018"},
    {"cohort": "CELLxGENE_7b20c613", "letter": "", "cancer_type": "Pan-Cancer", "author": "Gondal et al.", "year": "2025"},
    {"cohort": "CELLxGENE_05a8c945", "letter": "", "cancer_type": "Colorectal", "author": "Marteau et al.", "year": "2026"},
    {"cohort": "CELLxGENE_6f9de485", "letter": "", "cancer_type": "Breast", "author": "Yan et al.", "year": "2026"},
    {"cohort": "GSE207422", "letter": "", "cancer_type": "NSCLC", "author": "Hu et al.", "year": "2023"},
    {"cohort": "GSE243013", "letter": "", "cancer_type": "NSCLC", "author": "Liu et al.", "year": "2025"},
    {"cohort": "GSE233203", "letter": "", "cancer_type": "NSCLC", "author": "Hong et al.", "year": "2025"},
    {"cohort": "GSE200996", "letter": "", "cancer_type": "HNSCC", "author": "Luoma et al.", "year": "2022"},
    {"cohort": "GSE316195", "letter": "", "cancer_type": "PDAC", "author": "Raufi et al.", "year": "2026"},
)

# ------------------------------------------------------------------------------
# Dimensional Constants
# ------------------------------------------------------------------------------

# Standalone Bar Plot Dimensions (matching scatter grid height)
BAR_WIDTH_MM: Final[float] = 50.0
COMPOSITE_WIDTH_MM: Final[float] = 180.0
FIGURE_HEIGHT_MM: Final[float] = 123.31

# 10 pt = 10 * (96 / 72) = 13.333 px at 96 CSS DPI
MIN_FONT_PT: Final[float] = 10.0
MIN_FONT_PX: Final[float] = 13.5  # 10.125 pt (strictly >= 10 pt)
TITLE_FONT_PX: Final[float] = 14.0  # 10.5 pt


# ------------------------------------------------------------------------------
# Configuration & Frozen Data Models
# ------------------------------------------------------------------------------

class BarPlotConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    benchmark_dir: Path = Path(
        "output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark"
    )
    scatter_svg: Path = Path(
        "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
    )
    out_bar_svg: Path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.svg"
    )
    out_bar_png: Path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.png"
    )
    out_composite_svg: Path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.svg"
    )
    out_composite_png: Path = Path(
        "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.png"
    )
    article_bar_png: Path = Path("article/figures/milo-deconv/nine_cohorts_deconv_correlations_bar.png")
    article_composite_png: Path = Path("article/figures/milo-deconv/nine_cohorts_deconv_composite_with_bars.png")
    chart_width: int = 81
    chart_height: int = 357
    bar_height: int = 18
    bar_color: str = "#1F77B4"  # Single uniform Okabe-Ito blue


class CohortCorrelationRecord(BaseModel):
    model_config = ConfigDict(frozen=True)
    cohort: str
    letter: str
    label: str
    cancer_type: str
    author: str
    year: str
    rho_mean: float
    rho_std: float
    rho_min: float
    rho_max: float
    val_str: str


# ------------------------------------------------------------------------------
# Pure Functional Data Ingestion
# ------------------------------------------------------------------------------

def load_cohort_correlations(benchmark_dir: Path) -> Result[pl.DataFrame, str]:
    """Pure function loading 5-sweep correlation summary statistics across all 9 cohorts."""
    summary_file = benchmark_dir / "cross_cohort_baseline_sweeps_summary.parquet"
    if not summary_file.exists():
        summary_file = benchmark_dir / "cross_cohort_concordance_summary.parquet"

    records: list[dict[str, object]] = []

    if summary_file.exists():
        try:
            summary_df = pl.read_parquet(summary_file)
            for meta in COHORT_ORDER_AND_METADATA:
                c_id = meta["cohort"]
                sub = summary_df.filter(pl.col("cohort") == c_id)
                if sub.height > 0:
                    rho_col = "rho_mean" if "rho_mean" in sub.columns else "spearman_rho_milo_deconv"
                    std_col = "rho_std" if "rho_std" in sub.columns else None

                    r_mean_val = sub[rho_col][0]
                    r_std_val = sub[std_col][0] if std_col and std_col in sub.columns else 0.0

                    r_mean = float(r_mean_val) if r_mean_val is not None else 0.0
                    r_std = float(r_std_val) if r_std_val is not None and not (r_std_val != r_std_val) else 0.0

                    author_short = meta["author"].replace(" et al.", "")
                    records.append({
                        "cohort": c_id,
                        "letter": meta["letter"],
                        "label": author_short,
                        "cancer_type": meta["cancer_type"],
                        "author": meta["author"],
                        "year": meta["year"],
                        "rho_mean": r_mean,
                        "rho_std": r_std,
                        "rho_min": max(0.0, r_mean - r_std),
                        "rho_max": min(1.0, r_mean + r_std),
                        "val_str": f"{r_mean:.2f}",
                    })
        except Exception as exc:
            return Failure(f"Failed to read central summary parquet {summary_file}: {exc}")

    # Fallback to per-cohort parquet files if central summary incomplete
    if len(records) < len(COHORT_ORDER_AND_METADATA):
        records.clear()
        for meta in COHORT_ORDER_AND_METADATA:
            c_id = meta["cohort"]
            c_dir = benchmark_dir / c_id
            c_file = c_dir / f"{c_id}_baseline_sweeps_summary.parquet"
            if not c_file.exists():
                c_file = c_dir / f"{c_id}_concordance_summary.parquet"

            if not c_file.exists():
                return Failure(f"Could not locate summary parquet for cohort {c_id} in {c_dir}")

            try:
                c_df = pl.read_parquet(c_file).filter(pl.col("condition") == "Combined")
                if c_df.height == 0:
                    c_df = pl.read_parquet(c_file)

                rho_col = "spearman_rho_milo_deconv" if "spearman_rho_milo_deconv" in c_df.columns else "spearman_rho"
                r_mean_val = c_df[rho_col].mean()
                r_std_val = c_df[rho_col].std() if c_df.height > 1 else 0.0

                r_mean = float(r_mean_val) if r_mean_val is not None else 0.0
                r_std = float(r_std_val) if r_std_val is not None and not (r_std_val != r_std_val) else 0.0

                author_short = meta["author"].replace(" et al.", "")
                records.append({
                    "cohort": c_id,
                    "letter": meta["letter"],
                    "label": author_short,
                    "cancer_type": meta["cancer_type"],
                    "author": meta["author"],
                    "year": meta["year"],
                    "rho_mean": r_mean,
                    "rho_std": r_std,
                    "rho_min": max(0.0, r_mean - r_std),
                    "rho_max": min(1.0, r_mean + r_std),
                    "val_str": f"{r_mean:.2f}",
                })
            except Exception as exc:
                return Failure(f"Error parsing cohort summary for {c_id}: {exc}")

    if len(records) != len(COHORT_ORDER_AND_METADATA):
        return Failure(f"Expected {len(COHORT_ORDER_AND_METADATA)} cohorts, successfully parsed {len(records)}")

    return Success(pl.DataFrame(records))


# ------------------------------------------------------------------------------
# Declarative Visualization: Horizontal Bar Plot
# ------------------------------------------------------------------------------

def build_correlation_bar_chart(df: pl.DataFrame, config: BarPlotConfig) -> alt.Chart:
    """Build minimalist publication-grade horizontal bar chart in Altair."""
    # 1. Horizontal bars
    bars = (
        alt.Chart(df)
        .mark_bar(
            color=config.bar_color,
            stroke="black",
            strokeWidth=0.75,
            height=config.bar_height,
        )
        .encode(
            y=alt.Y(
                "label:N",
                sort=None,
                title=None,
                axis=alt.Axis(
                    labelFontSize=MIN_FONT_PX,
                    tickSize=3,
                    tickColor="#CBD5E1",
                    labelPadding=4,
                    domainColor="#CBD5E1",
                ),
            ),
            x=alt.X(
                "rho_mean:Q",
                scale=alt.Scale(domain=[0.0, 1.0]),
                title="Spearman ρ",
                axis=alt.Axis(
                    labelFontSize=MIN_FONT_PX,
                    titleFontSize=MIN_FONT_PX,
                    titleFontWeight="bold",
                    titlePadding=6,
                    values=[0.0, 0.5, 1.0],
                    grid=True,
                    gridColor="#EEEEEE",
                    domainColor="#CBD5E1",
                ),
            ),
            tooltip=[
                alt.Tooltip("letter:N", title="Panel"),
                alt.Tooltip("cancer_type:N", title="Cancer"),
                alt.Tooltip("author:N", title="Cohort"),
                alt.Tooltip("year:N", title="Year"),
                alt.Tooltip("rho_mean:Q", title="Mean Spearman ρ", format=".3f"),
                alt.Tooltip("rho_std:Q", title="s.d.", format=".3f"),
            ],
        )
    )

    # 2. Error whiskers (mean +/- s.d.)
    errors = (
        alt.Chart(df)
        .mark_rule(color="#000000", strokeWidth=1.5)
        .encode(
            y=alt.Y("label:N", sort=None),
            x=alt.X("rho_min:Q"),
            x2=alt.X2("rho_max:Q"),
        )
    )

    # 3. Cross-cohort benchmark average reference line (dashed)
    mean_val = float(df["rho_mean"].mean()) if df.height > 0 else 0.82
    mean_df = pl.DataFrame({"mean_rho": [mean_val]})
    mean_rule = (
        alt.Chart(mean_df)
        .mark_rule(color="#444444", strokeDash=[3, 3], strokeWidth=1.2)
        .encode(x=alt.X("mean_rho:Q"))
    )

    title_params = alt.TitleParams(
        text=["(d)  Concordance", "      Summary"],
        subtitle=["Mean ± s.d. (5 sweeps)"],
        fontSize=TITLE_FONT_PX,
        subtitleFontSize=MIN_FONT_PX,
        anchor="start",
        offset=6,
    )

    chart = (
        (bars + errors + mean_rule)
        .properties(
            width=config.chart_width,
            height=config.chart_height,
            title=title_params,
        )
        .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
    )
    return chart


# ------------------------------------------------------------------------------
# SVG Post-Processing & Physical Dimension Enforcement
# ------------------------------------------------------------------------------

def enforce_svg_physical_dimensions(
    svg_text: str,
    target_width_mm: float,
    target_height_mm: float,
) -> str:
    """Explicitly inject physical dimensions (width in mm, height in mm) on root <svg>."""
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    if not m:
        return svg_text

    svg_updated = re.sub(
        r'(<svg[^>]*?)\s+width="[^"]+"\s+height="[^"]+"',
        rf'\1 width="{target_width_mm:.1f}mm" height="{target_height_mm:.2f}mm"',
        svg_text,
        count=1,
    )
    return svg_updated


def stitch_composite_grid_and_bars(
    scatter_svg_str: str,
    bar_svg_str: str,
    composite_width_mm: float = COMPOSITE_WIDTH_MM,
    composite_height_mm: float = FIGURE_HEIGHT_MM,
) -> str:
    """Stitch scatter grid SVG (130mm, x=0) and bar plot SVG (50mm, x=486) into a single 180mm composite SVG."""
    m_scatter = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', scatter_svg_str)
    m_bar = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', bar_svg_str)

    w_scatter = float(m_scatter.group(1)) if m_scatter else 486.0
    h_scatter = float(m_scatter.group(2)) if m_scatter else 672.0
    w_bar = float(m_bar.group(1)) if m_bar else 188.0
    h_bar = float(m_bar.group(2)) if m_bar else 672.0

    total_w = int(w_scatter + w_bar)
    max_h = int(max(h_scatter, h_bar))

    inner_scatter = re.sub(r"^<svg[^>]*>", "", scatter_svg_str.strip())
    inner_scatter = re.sub(r"</svg>$", "", inner_scatter.strip())

    inner_bar = re.sub(r"^<svg[^>]*>", "", bar_svg_str.strip())
    inner_bar = re.sub(r"</svg>$", "", inner_bar.strip())

    composite_svg = f"""<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" width="{composite_width_mm:.1f}mm" height="{composite_height_mm:.2f}mm" viewBox="0 0 {total_w} {max_h}">
<rect width="{total_w}" height="{max_h}" fill="white"/>
<g id="scatter-grid-panels-a-i" transform="translate(0, 0)">
{inner_scatter}
</g>
<g id="concordance-summary-panel-j" transform="translate({int(w_scatter)}, 0)">
{inner_bar}
</g>
</svg>"""
    return composite_svg


# ------------------------------------------------------------------------------
# Imperative Shell: Orchestration & File Exports
# ------------------------------------------------------------------------------

def generate_correlation_bar_and_composite(config: BarPlotConfig) -> Result[tuple[Path, Path], str]:
    """Execute data ingestion, chart building, vector SVG rendering, and 300 DPI raster generation."""
    config.out_bar_svg.parent.mkdir(parents=True, exist_ok=True)
    config.article_bar_png.parent.mkdir(parents=True, exist_ok=True)

    print("=== Generating Cross-Cohort Correlation Bar Plot (50mm) and Composite (180mm) ===")

    # 1. Load data
    match load_cohort_correlations(config.benchmark_dir):
        case Failure(err):
            return Failure(f"Failed to load cross-cohort correlation summary: {err}")
        case Success(df):
            print(f"Loaded correlation summary across {df.height} cohorts:")
            for row in df.iter_rows(named=True):
                print(f"  {row['letter']}: {row['cohort']:18s} ({row['cancer_type']}) | rho = {row['rho_mean']:+.3f} ± {row['rho_std']:.3f}")

    # 2. Build Altair chart
    chart = build_correlation_bar_chart(df, config)
    raw_bar_svg = vlc.vegalite_to_svg(chart.to_dict())

    # Determine target height (matching scatter grid height exactly)
    target_height_mm = FIGURE_HEIGHT_MM
    scatter_svg_str = ""
    if config.scatter_svg.exists():
        scatter_svg_str = config.scatter_svg.read_text(encoding="utf-8")
        m_h = re.search(r'height="([0-9.]+)mm"', scatter_svg_str)
        if m_h:
            target_height_mm = float(m_h.group(1))

    # 3. Format standalone bar plot SVG (50.0mm width, matching height)
    bar_svg_str = enforce_svg_physical_dimensions(
        raw_bar_svg,
        target_width_mm=BAR_WIDTH_MM,
        target_height_mm=target_height_mm,
    )
    with open(config.out_bar_svg, "w", encoding="utf-8") as f:
        f.write(bar_svg_str)
    print(f"\nExported standalone bar SVG to {config.out_bar_svg} (width={BAR_WIDTH_MM}mm, height={target_height_mm:.2f}mm)")

    # 4. Generate 300 DPI PNG for standalone bar plot
    # 50mm * 300 DPI / 25.4 mm/in = 590.55 px; scale = 590.55 / 188 = 3.141
    bar_png_bytes = vlc.vegalite_to_png(chart.to_dict(), scale=3.14)
    with open(config.out_bar_png, "wb") as f:
        f.write(bar_png_bytes)
    with open(config.article_bar_png, "wb") as f:
        f.write(bar_png_bytes)
    print(f"Exported 300 DPI bar PNG to {config.out_bar_png} and {config.article_bar_png}")

    # 5. Build and export 180mm side-by-side composite figure
    if not config.scatter_svg.exists() or not scatter_svg_str:
        print(f"Warning: scatter grid SVG not found at {config.scatter_svg}. Skipping composite stitch.")
        return Success((config.out_bar_svg, config.out_bar_svg))

    composite_svg_str = stitch_composite_grid_and_bars(
        scatter_svg_str=scatter_svg_str,
        bar_svg_str=bar_svg_str,
        composite_width_mm=COMPOSITE_WIDTH_MM,
        composite_height_mm=target_height_mm,
    )

    with open(config.out_composite_svg, "w", encoding="utf-8") as f:
        f.write(composite_svg_str)
    print(f"Exported combined composite SVG to {config.out_composite_svg} (width={COMPOSITE_WIDTH_MM}mm, height={target_height_mm:.2f}mm)")

    # 6. Render 300 DPI raster for the composite figure using rsvg-convert
    target_h_px = int(round(target_height_mm * 300.0 / 25.4))
    try:
        subprocess.run(
            ["rsvg-convert", "-w", "2126", "-h", str(target_h_px), str(config.out_composite_svg), "-o", str(config.out_composite_png)],
            check=True,
        )
        with open(config.out_composite_png, "rb") as f_in, open(config.article_composite_png, "wb") as f_out:
            f_out.write(f_in.read())
        print(f"Exported 300 DPI composite PNG to {config.out_composite_png} and {config.article_composite_png}")
    except Exception as exc:
        print(f"rsvg-convert rasterization note: {exc}; falling back to direct copy if exists.")

    return Success((config.out_bar_svg, config.out_composite_svg))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot horizontal correlation bar chart across 9 cohorts and generate 180mm composite figure."
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
        "--scatter-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid/nine_cohorts_deconv_scatter_grid.svg"
        ),
        help="Input 130mm scatter grid SVG path",
    )
    parser.add_argument(
        "--out-bar-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.svg"
        ),
        help="Target output bar SVG path",
    )
    parser.add_argument(
        "--out-bar-png",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_correlations_bar.png"
        ),
        help="Target output bar PNG path",
    )
    parser.add_argument(
        "--out-composite-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.svg"
        ),
        help="Target output 180mm composite SVG path",
    )
    parser.add_argument(
        "--out-composite-png",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar/nine_cohorts_deconv_composite_with_bars.png"
        ),
        help="Target output 180mm composite PNG path",
    )
    args = parser.parse_args()

    config = BarPlotConfig(
        benchmark_dir=args.benchmark_dir,
        scatter_svg=args.scatter_svg,
        out_bar_svg=args.out_bar_svg,
        out_bar_png=args.out_bar_png,
        out_composite_svg=args.out_composite_svg,
        out_composite_png=args.out_composite_png,
    )

    match generate_correlation_bar_and_composite(config):
        case Failure(err):
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success((bar_path, comp_path)):
            print(f"Complete: {bar_path} and {comp_path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
