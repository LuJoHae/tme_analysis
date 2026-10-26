#!/usr/bin/env python3
"""
Publication-Grade Cross-Cohort Differential Abundance Diverging Bar Plot (Panel b).
Generates Panel b of the master 2x2 immunotherapy figure.
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

TARGET_WIDTH_MM: Final[float] = 50.0
TARGET_HEIGHT_MM: Final[float] = 123.31

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_VERMILION: Final[str] = "#D55E00"
COLOR_HAIRLINE: Final[str] = "#CBD5E1"

# ------------------------------------------------------------------------------
# Cohort Metadata & Ordering Matrix (Exact match with deconvolution Panel d)
# ------------------------------------------------------------------------------
COHORT_ORDER_AND_METADATA: Final[tuple[dict[str, str], ...]] = (
    {"cohort": "GSE120575", "letter": "b", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018"},
    {"cohort": "CELLxGENE_7b20c613", "letter": "", "cancer_type": "Pan-Cancer", "author": "Gondal et al.", "year": "2025"},
    {"cohort": "CELLxGENE_05a8c945", "letter": "", "cancer_type": "Colorectal", "author": "Marteau et al.", "year": "2026"},
    {"cohort": "CELLxGENE_6f9de485", "letter": "", "cancer_type": "Breast", "author": "Yan et al.", "year": "2026"},
    {"cohort": "GSE207422", "letter": "", "cancer_type": "NSCLC", "author": "Hu et al.", "year": "2023"},
    {"cohort": "GSE243013", "letter": "", "cancer_type": "NSCLC", "author": "Liu et al.", "year": "2025"},
    {"cohort": "GSE233203", "letter": "", "cancer_type": "NSCLC", "author": "Hong et al.", "year": "2025"},
    {"cohort": "GSE200996", "letter": "", "cancer_type": "HNSCC", "author": "Luoma et al.", "year": "2022"},
    {"cohort": "GSE316195", "letter": "", "cancer_type": "PDAC", "author": "Raufi et al.", "year": "2026"},
)


class DaBarConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    results_dir: Path = Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts")
    out_svg: Path = Path(
        "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.svg"
    )
    out_png: Path = Path(
        "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.png"
    )
    article_png: Path = Path("article/figures/single-cell/nine_cohorts_da_bars.png")
    chart_width: int = 81
    chart_height: int = 321
    bar_height: int = 18
    fdr_threshold: float = 0.10


def load_all_cohort_da_stats(config: DaBarConfig) -> Result[pd.DataFrame, str]:
    """Loads differential abundance statistics across all 9 cohorts."""
    rows: list[dict[str, object]] = []

    for meta in COHORT_ORDER_AND_METADATA:
        cohort_id = meta["cohort"]
        author_short = meta["author"].replace(" et al.", "")
        da_file = config.results_dir / cohort_id / "da_results.parquet"
        if not da_file.exists():
            return Failure(f"Missing DA results file at {da_file}")

        try:
            df = pl.read_parquet(da_file)
            n_total = len(df)
            if n_total == 0:
                return Failure(f"Empty DA results table for {cohort_id}")

            sig_mask = pl.col("FDR") < config.fdr_threshold
            up_count = len(df.filter(sig_mask & (pl.col("logFC") > 0)))
            down_count = len(df.filter(sig_mask & (pl.col("logFC") < 0)))

            pct_up = (up_count / n_total) * 100.0
            pct_down = (down_count / n_total) * 100.0

            # Two entries per cohort: one for Responders, one for Non-Responders
            rows.append({
                "cohort": cohort_id,
                "label": author_short,
                "direction": "R",
                "pct": pct_up,
                "author": meta["author"],
                "cancer_type": meta["cancer_type"],
            })
            rows.append({
                "cohort": cohort_id,
                "label": author_short,
                "direction": "NR",
                "pct": -pct_down,
                "author": meta["author"],
                "cancer_type": meta["cancer_type"],
            })
        except Exception as exc:
            return Failure(f"Error loading DA results for {cohort_id}: {exc}")

    return Success(pd.DataFrame(rows))


def build_da_diverging_bar_chart(mdf: pd.DataFrame, config: DaBarConfig) -> alt.Chart:
    """Builds diverging horizontal bar chart for Panel b."""
    cohort_order = [m["author"].replace(" et al.", "") for m in COHORT_ORDER_AND_METADATA]

    # Zero vertical hairline reference
    zero_rule = (
        alt.Chart(pd.DataFrame([{"x": 0.0}]))
        .mark_rule(color="#444444", strokeDash=[3, 3], strokeWidth=1.0)
        .encode(x=alt.X("x:Q"))
    )

    bars = (
        alt.Chart(mdf)
        .mark_bar(
            height=config.bar_height,
            stroke="black",
            strokeWidth=0.75,
        )
        .encode(
            y=alt.Y(
                "label:N",
                sort=cohort_order,
                title=None,
                axis=alt.Axis(
                    labelFontSize=MIN_FONT_PX,
                    tickSize=3,
                    tickColor=COLOR_HAIRLINE,
                    labelPadding=4,
                    domainColor=COLOR_HAIRLINE,
                ),
            ),
            x=alt.X(
                "pct:Q",
                scale=alt.Scale(domain=[-25.0, 75.0]),
                title="% Significant",
                axis=alt.Axis(
                    labelFontSize=MIN_FONT_PX,
                    titleFontSize=MIN_FONT_PX,
                    titleFontWeight="bold",
                    titlePadding=6,
                    values=[-20, 0, 20, 40, 60],
                    grid=True,
                    gridColor="#EEEEEE",
                    domainColor=COLOR_HAIRLINE,
                ),
            ),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(
                    domain=["R", "NR"],
                    range=[OKABE_VERMILION, OKABE_BLUE],
                ),
                legend=alt.Legend(
                    orient="bottom",
                    direction="horizontal",
                    title=None,
                    labelFontSize=MIN_FONT_PX,
                    symbolSize=90,
                    padding=2,
                    labelExpr="datum.label == 'R' ? 'R' : 'NR'",
                ),
            ),
            tooltip=[
                alt.Tooltip("label:N", title="Cohort"),
                alt.Tooltip("direction:N", title="Direction"),
                alt.Tooltip("pct:Q", title="% Neighborhoods", format="+.1f"),
            ],
        )
    )

    title_params = alt.TitleParams(
        text=["(b)  Differential", "      Abundance"],
        subtitle=["Milo Nhoods (FDR < 0.10)"],
        fontSize=TITLE_FONT_PX,
        subtitleFontSize=MIN_FONT_PX,
        anchor="start",
        offset=6,
    )

    chart = (
        (bars + zero_rule)
        .properties(
            width=config.chart_width,
            height=config.chart_height,
            title=title_params,
        )
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )
    return chart


def enforce_svg_physical_dimensions(
    svg_text: str,
    target_width_mm: float = TARGET_WIDTH_MM,
    target_height_mm: float = TARGET_HEIGHT_MM,
) -> str:
    """Enforces exact mm width and height on root <svg>."""
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


def generate_nine_cohorts_da_bars(config: DaBarConfig) -> Result[Path, str]:
    """Generates Panel b diverging bar chart vector SVG and 300 DPI raster PNG."""
    config.out_svg.parent.mkdir(parents=True, exist_ok=True)
    config.article_png.parent.mkdir(parents=True, exist_ok=True)

    print("=== Constructing Cross-Cohort Differential Abundance Diverging Bars (Panel b) ===")
    match load_all_cohort_da_stats(config):
        case Failure(err):
            return Failure(f"Failed to load DA statistics: {err}")
        case Success(mdf):
            for name, group in mdf.groupby("label", sort=False):
                r_val = group[group["direction"] == "R"]["pct"].values[0]
                nr_val = group[group["direction"] == "NR"]["pct"].values[0]
                print(f"  {name:15s}: Enriched R={r_val:5.1f}% | Enriched NR={abs(nr_val):5.1f}%")

            chart = build_da_diverging_bar_chart(mdf, config)
            chart_dict = chart.to_dict()

            print(f"\nExporting publication vector SVG (width={TARGET_WIDTH_MM:.0f}mm, height={TARGET_HEIGHT_MM:.2f}mm) to {config.out_svg}...")
            raw_svg = vlc.vegalite_to_svg(chart_dict)
            svg_str = enforce_svg_physical_dimensions(raw_svg, TARGET_WIDTH_MM, TARGET_HEIGHT_MM)

            with open(config.out_svg, "w", encoding="utf-8") as f:
                f.write(svg_str)

            print(f"Exporting 300 DPI raster PNG to {config.out_png}...")
            png_bytes = vlc.vegalite_to_png(chart_dict, scale=3.125)
            with open(config.out_png, "wb") as f:
                f.write(png_bytes)

            with open(config.article_png, "wb") as f:
                f.write(png_bytes)

            print("\nSuccessfully exported Panel b DA bar plot assets!")
            return Success(config.out_svg)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot cross-cohort differential abundance diverging bar chart.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts"),
        help="Root milopy analysis directory",
    )
    parser.add_argument(
        "--out-svg",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.svg"
        ),
        help="Output SVG path",
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/4a4_plot_cross_cohort_da_bars/nine_cohorts_da_bars.png"
        ),
        help="Output PNG path",
    )
    args = parser.parse_args()

    config = DaBarConfig(
        results_dir=args.results_dir,
        out_svg=args.out_svg,
        out_png=args.out_png,
    )

    match generate_nine_cohorts_da_bars(config):
        case Failure(err):
            print(f"Error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"Complete: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
