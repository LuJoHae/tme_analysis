#!/usr/bin/env python3
"""Cross-Cohort Immunotherapy Differential Abundance Portrait Figure Generator.

Generates a publication-grade Nature Methods portrait-mode diverging horizontal bar
chart summarizing single-cell differential abundance (Milo QL-GLM) across all
evaluated clinical immunotherapy cohorts.

Adheres strictly to:
- Pure functional programming (returns Result, immutable Pydantic models).
- Polars for declarative data manipulation.
- Altair + vl-convert for vector SVG and 300 DPI PNG rasterization.
- Nature Methods minimalist wireframe standard with Okabe-Ito colorblind-safe palette.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Final, Literal, Sequence

import altair as alt  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.result import Failure, Result, Success
import vl_convert as vlc  # type: ignore

from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    COHORT_METADATA,
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
    get_single_cell_cohort_display_label,
)

# Canonical ordering by cancer type and year for coherent biological comparison
CANCER_TYPE_ORDER: Final[tuple[str, ...]] = (
    "Melanoma",
    "Colorectal",
    "Breast",
    "NSCLC",
    "HNSCC",
    "PDAC",
)


class CohortSummaryRecord(BaseModel):
    """Immutable data record representing one cohort's DA statistics."""

    model_config = ConfigDict(frozen=True)

    accession: str
    cohort_label: str
    cancer_type: str
    author: str
    year: int
    n_total_nhoods: int
    n_sig_up: int
    n_sig_down: int
    pct_up: float
    pct_down: float


class PlotConfig(BaseModel):
    """Immutable configuration for cross-cohort summary plotting."""

    model_config = ConfigDict(frozen=True)

    results_dir: Path = Path("output/results/milopy_response_analysis")
    output_svg: Path = Path("output/reports/milopy_cohort_comparison_portrait.svg")
    output_png: Path | None = None
    sort_by: Literal["cancer_type", "response", "accession"] = "cancer_type"
    label_style: Literal["cancer_author_year_3line", "cancer_author_year", "author_year_cancer", "comprehensive"] = "cancer_author_year_3line"
    multiline_labels: bool = True
    chart_width: int = 180
    chart_height: int = 560
    bar_size: int = 18
    fdr_threshold: float = 0.10
    title: tuple[str, ...] = ("Cross-Cohort Immunotherapy", "Differential Abundance")
    subtitle: tuple[str, ...] = (
        "Percentage of neighborhoods significantly",
        "shifted by response outcome (Milo QL-GLM)",
    )


def load_cohort_da_record(
    accession: str,
    results_dir: Path,
    label_style: str = "cancer_author_year",
    fdr_threshold: float = 0.10,
) -> Result[CohortSummaryRecord, str]:
    """Loads and computes differential abundance statistics for a single cohort."""
    cohort_dir = results_dir / accession
    parquet_path = cohort_dir / "da_results.parquet"
    csv_path = cohort_dir / "da_results.csv"

    if parquet_path.exists():
        df = pl.read_parquet(parquet_path)
    elif csv_path.exists():
        df = pl.read_csv(csv_path)
    else:
        return Failure(f"No DA results found for {accession} in {cohort_dir}")

    meta = COHORT_METADATA.get(accession)
    cancer_type = meta["indication"] if meta else "Other"
    author = meta["author"] if meta else accession
    year = int(meta["year"]) if meta and meta.get("year") else 2020
    display_label = get_single_cell_cohort_display_label(accession, format_style=label_style)

    n_total = len(df)
    if n_total == 0:
        return Failure(f"DA results table for {accession} is empty")

    sig_mask = pl.col("FDR") < fdr_threshold
    up_count = len(df.filter(sig_mask & (pl.col("logFC") > 0)))
    down_count = len(df.filter(sig_mask & (pl.col("logFC") < 0)))

    pct_up = (up_count / n_total) * 100.0
    pct_down = (down_count / n_total) * 100.0

    return Success(
        CohortSummaryRecord(
            accession=accession,
            cohort_label=display_label,
            cancer_type=cancer_type,
            author=author,
            year=year,
            n_total_nhoods=n_total,
            n_sig_up=up_count,
            n_sig_down=down_count,
            pct_up=pct_up,
            pct_down=pct_down,
        )
    )


def collect_all_cohort_records(
    config: PlotConfig,
) -> Result[tuple[CohortSummaryRecord, ...], str]:
    """Collects DA summary records across all standard 9 single-cell cohorts."""
    cohort_accessions = tuple(p["acc"] for p in COHORT_PANELS)
    records: list[CohortSummaryRecord] = []
    errors: list[str] = []

    for acc in cohort_accessions:
        match load_cohort_da_record(acc, config.results_dir, config.label_style, config.fdr_threshold):
            case Success(record):
                records.append(record)
            case Failure(err):
                errors.append(err)

    if not records:
        return Failure(f"Failed to load any cohort DA records: {'; '.join(errors)}")

    # Sort cohorts
    match config.sort_by:
        case "cancer_type":
            cancer_rank = {c: i for i, c in enumerate(CANCER_TYPE_ORDER)}
            # Primary sort by cancer rank, secondary by year, then author
            records.sort(key=lambda r: (cancer_rank.get(r.cancer_type, 99), r.year, r.author))
        case "response":
            # Sort by net responder percentage (pct_up - pct_down) descending
            records.sort(key=lambda r: (r.pct_up - r.pct_down), reverse=True)
        case "accession":
            records.sort(key=lambda r: r.accession)

    return Success(tuple(records))


def build_plotting_dataframe(
    records: tuple[CohortSummaryRecord, ...],
) -> pl.DataFrame:
    """Transforms cohort records into long-format tidy Polars DataFrame for Altair."""
    rows: list[dict[str, object]] = []

    for r in records:
        # Responder-enriched (positive percentage)
        rows.append({
            "accession": r.accession,
            "cohort_label": r.cohort_label,
            "cancer_type": r.cancer_type,
            "author": r.author,
            "year": r.year,
            "direction": "Enriched in Responders",
            "percentage": r.pct_up,
            "percentage_plot": r.pct_up,
            "count": r.n_sig_up,
            "total_nhoods": r.n_total_nhoods,
            "label_text": f"+{r.pct_up:.1f}%" if r.pct_up > 0 else "",
            "label_dx": 5,
            "label_align": "left",
        })

        # Non-responder-enriched (negative percentage)
        rows.append({
            "accession": r.accession,
            "cohort_label": r.cohort_label,
            "cancer_type": r.cancer_type,
            "author": r.author,
            "year": r.year,
            "direction": "Enriched in Non-Responders",
            "percentage": r.pct_down,
            "percentage_plot": -r.pct_down,
            "count": r.n_sig_down,
            "total_nhoods": r.n_total_nhoods,
            "label_text": f"-{r.pct_down:.1f}%" if r.pct_down > 0 else "",
            "label_dx": -5,
            "label_align": "right",
        })

    return pl.DataFrame(rows)


def render_portrait_comparison_chart(
    data_df: pl.DataFrame,
    ordered_labels: Sequence[str],
    config: PlotConfig,
) -> alt.Chart:
    """Builds a Nature Methods compliant portrait diverging bar chart."""
    pdf = data_df.to_pandas()

    # Dynamic x-axis domain matching data extent with generous padding for bar-end labels on compressed scale
    min_val = float(pdf["percentage_plot"].min())
    max_val = float(pdf["percentage_plot"].max())
    min_pct = min(-45.0, min_val - 12.0)
    max_pct = max(80.0, max_val + 12.0)
    x_domain = [float(np.floor(min_pct / 5.0) * 5.0), float(np.ceil(max_pct / 5.0) * 5.0)]

    base = alt.Chart(pdf)

    # 1. Zero vertical reference line
    zero_rule = (
        alt.Chart(pd.DataFrame([{"x": 0.0}]))
        .mark_rule(color=COLOR_TEXT_MUTED, strokeWidth=0.85, strokeDash=[3, 3])
        .encode(x=alt.X("x:Q"))
    )

    # Label formatting for multi-line cohort labels:
    # Splits "Cancer Type · First Author et al. · (Year)" into 3 lines
    y_axis_kwargs: dict[str, object] = {
        "domain": True,
        "domainColor": COLOR_HAIRLINE,
        "ticks": True,
        "tickColor": COLOR_HAIRLINE,
        "tickSize": 5,
        "labelLimit": 0,
        "labelFont": "sans-serif",
        "labelFontSize": 9.0,
        "labelFontWeight": "normal",
        "labelColor": COLOR_TEXT_PRIMARY,
        "titleFont": "sans-serif",
        "titleFontSize": 10.5,
        "titleFontWeight": "bold",
        "titlePadding": 12,
    }
    if config.multiline_labels:
        y_axis_kwargs["labelExpr"] = "split(datum.label, ' · ')"
        y_axis_kwargs["labelLineHeight"] = 11.5

    # 2. Diverging Bars
    bars = (
        base.mark_bar(size=config.bar_size, cornerRadius=0)
        .encode(
            x=alt.X(
                "percentage_plot:Q",
                title="Significant Neighborhoods (%)",
                scale=alt.Scale(domain=x_domain),
                axis=alt.Axis(
                    values=[-40, -20, 0, 20, 40, 60, 80],
                    grid=True,
                    gridColor="#E2E8F0",
                    gridWidth=0.5,
                    domain=True,
                    domainColor=COLOR_HAIRLINE,
                    domainWidth=0.75,
                    tickColor=COLOR_HAIRLINE,
                    labelFont="sans-serif",
                    labelFontSize=8.5,
                    titleFont="sans-serif",
                    titleFontSize=10.5,
                    titleFontWeight="bold",
                ),
            ),
            y=alt.Y(
                "cohort_label:N",
                title="Cohort & Primary Reference",
                sort=list(ordered_labels),
                axis=alt.Axis(**y_axis_kwargs),
            ),
            color=alt.Color(
                "direction:N",
                title="Direction (FDR < 0.10)",
                scale=alt.Scale(
                    domain=["Enriched in Responders", "Enriched in Non-Responders"],
                    range=[OKABE_VERMILION, OKABE_BLUE],
                ),
                legend=alt.Legend(
                    orient="top",
                    titleFont="sans-serif",
                    titleFontSize=10.5,
                    titleFontWeight="bold",
                    labelFont="sans-serif",
                    labelFontSize=9,
                    symbolType="square",
                    symbolSize=90,
                    direction="vertical",
                    columns=1,
                    padding=8,
                ),
            ),
            tooltip=[
                alt.Tooltip("cohort_label:N", title="Cohort"),
                alt.Tooltip("accession:N", title="Accession"),
                alt.Tooltip("cancer_type:N", title="Cancer Type"),
                alt.Tooltip("author:N", title="Primary Author"),
                alt.Tooltip("year:Q", title="Publication Year"),
                alt.Tooltip("direction:N", title="Enrichment"),
                alt.Tooltip("percentage:Q", format=".1f", title="% Neighborhoods"),
                alt.Tooltip("count:Q", format=",d", title="Significant Count"),
                alt.Tooltip("total_nhoods:Q", format=",d", title="Total Evaluated Nhoods"),
            ],
        )
    )

    # 3. Direct numerical labels at bar ends
    text_up = (
        alt.Chart(pdf[pdf["percentage_plot"] > 0])
        .mark_text(
            font="sans-serif",
            fontSize=8.5,
            fontWeight="bold",
            color=COLOR_TEXT_PRIMARY,
            dx=4,
            align="left",
            baseline="middle",
        )
        .encode(
            x=alt.X("percentage_plot:Q"),
            y=alt.Y("cohort_label:N", sort=list(ordered_labels)),
            text=alt.Text("label_text:N"),
        )
    )

    text_down = (
        alt.Chart(pdf[pdf["percentage_plot"] < 0])
        .mark_text(
            font="sans-serif",
            fontSize=8.5,
            fontWeight="bold",
            color=COLOR_TEXT_PRIMARY,
            dx=-4,
            align="right",
            baseline="middle",
        )
        .encode(
            x=alt.X("percentage_plot:Q"),
            y=alt.Y("cohort_label:N", sort=list(ordered_labels)),
            text=alt.Text("label_text:N"),
        )
    )

    # Composite chart in portrait aspect ratio
    chart = (
        alt.layer(zero_rule, bars, text_up, text_down)
        .properties(
            width=config.chart_width,
            height=config.chart_height,
            title=alt.TitleParams(
                text=list(config.title) if isinstance(config.title, tuple) else config.title,
                subtitle=list(config.subtitle) if isinstance(config.subtitle, tuple) else config.subtitle,
                font="sans-serif",
                fontSize=12,
                fontWeight="bold",
                subtitleFont="sans-serif",
                subtitleFontSize=9.5,
                subtitleColor=COLOR_TEXT_SECONDARY,
                anchor="start",
                offset=14,
            ),
        )
    )

    return apply_nature_methods_theme(chart)


def export_figures(
    chart: alt.Chart,
    svg_path: Path,
    png_path: Path | None = None,
) -> Result[tuple[Path, Path], str]:
    """Exports the Altair chart to SVG and 300 DPI PNG via vl-convert."""
    try:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        chart_json = chart.to_json()

        # Render SVG
        svg_str = vlc.vegalite_to_svg(chart_json)
        svg_path.write_text(svg_str, encoding="utf-8")

        # Render 300 DPI PNG
        target_png = png_path if png_path is not None else svg_path.with_suffix(".png")
        png_bytes = vlc.vegalite_to_png(chart_json, scale=3.0)
        target_png.write_bytes(png_bytes)

        return Success((svg_path, target_png))
    except Exception as exc:
        return Failure(f"Failed to export figures: {exc}")


def parse_args() -> PlotConfig:
    """Parses command-line arguments into an immutable PlotConfig."""
    parser = argparse.ArgumentParser(
        description="Generate publication-grade portrait cross-cohort differential abundance summary chart."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/results/milopy_response_analysis"),
        help="Path to milopy response analysis results directory",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("output/reports/milopy_cohort_comparison_portrait.svg"),
        help="Target output SVG file path",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=None,
        help="Target output PNG file path (defaults to same stem as SVG)",
    )
    parser.add_argument(
        "--sort-by",
        type=str,
        choices=["cancer_type", "response", "accession"],
        default="cancer_type",
        help="Cohort sorting strategy: cancer_type (default), response, or accession",
    )
    parser.add_argument(
        "--label-style",
        type=str,
        choices=["cancer_author_year_3line", "cancer_author_year", "author_year_cancer", "comprehensive"],
        default="cancer_author_year_3line",
        help="Label format: cancer_author_year_3line (3 lines), cancer_author_year (2 lines), author_year_cancer, comprehensive",
    )
    parser.add_argument(
        "--multiline-labels",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render cohort labels across multiple lines (Cancer Type / Author et al. / Year)",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=180,
        help="Inner plot width in points (default: 180 for compressed portrait x-scale)",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=560,
        help="Inner plot height in points (default: 560 for ~1:2 portrait ratio)",
    )
    parser.add_argument(
        "--bar-size",
        type=int,
        default=18,
        help="Thickness of horizontal bars in points (default: 18)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Nominal significance FDR threshold (default: 0.10)",
    )

    args = parser.parse_args()
    return PlotConfig(
        results_dir=args.results_dir,
        output_svg=args.output_svg,
        output_png=args.output_png,
        sort_by=args.sort_by,
        label_style=args.label_style,
        multiline_labels=args.multiline_labels,
        chart_width=args.width,
        chart_height=args.height,
        bar_size=args.bar_size,
        fdr_threshold=args.fdr_threshold,
    )


def main() -> None:
    config = parse_args()
    print("=" * 80)
    print("CROSS-COHORT IMMUNOTHERAPY DIFFERENTIAL ABUNDANCE (PORTRAIT FIGURE)")
    print("=" * 80)
    print(f"Results Directory : {config.results_dir}")
    print(f"Output SVG        : {config.output_svg}")
    print(f"Sorting Strategy  : {config.sort_by}")
    print(f"Label Format      : {config.label_style}")
    print(f"Dimensions        : {config.chart_width} pt x {config.chart_height} pt (portrait)")

    records_res = collect_all_cohort_records(config)
    match records_res:
        case Failure(err):
            print(f"Error loading cohort records: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(records):
            print(f"\nSuccessfully collected statistics for {len(records)} cohorts:")
            for r in records:
                print(
                    f"  - {r.cohort_label:<40s} | Nhoods: {r.n_total_nhoods:5d} | "
                    f"Resp Up: {r.pct_up:5.1f}% ({r.n_sig_up:4d}) | "
                    f"NR Down: {r.pct_down:5.1f}% ({r.n_sig_down:4d})"
                )

    plot_df = build_plotting_dataframe(records)
    ordered_labels = [r.cohort_label for r in records]

    print("\nRendering portrait Altair chart...")
    chart = render_portrait_comparison_chart(plot_df, ordered_labels, config)

    print("Exporting vector SVG and 300 DPI PNG via vl-convert...")
    export_res = export_figures(chart, config.output_svg, config.output_png)
    match export_res:
        case Failure(err):
            print(f"Export error: {err}", file=sys.stderr)
            sys.exit(1)
        case Success((svg_p, png_p)):
            print(f"\n[EXPORTED SVG] {svg_p}")
            print(f"[EXPORTED PNG] {png_p}")
            print("=" * 80)
            print("Completed successfully.")


if __name__ == "__main__":
    main()
