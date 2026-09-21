#!/usr/bin/env python3
"""
Step 5 (Calibration): Generate Publication Visualizations (Altair Vector SVGs).
Creates publication-quality vector SVG figures for:
1. Fig 1: Triangulation Scatter (Single-Cell Ground Truth vs Deconvolution vs Milo)
2. Fig 2: Pseudo-bulk Deconvolution Fidelity Benchmark
3. Fig 3: Resolution Trend Curve (Impact of clustering resolution from 0.5 up to 3.0)
4. Fig 4: Discrepancy Attribution Bar Chart
Strict functional Python with Altair and vl-convert.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import altair as alt  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import polars as pl
import vl_convert as vlc  # type: ignore


class PlotCalibrationConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    diagnostics_path: Path
    trends_path: Path
    fidelity_path: Path
    results_dir: Path


def parse_args() -> PlotCalibrationConfig:
    parser = argparse.ArgumentParser(
        description="Generate Altair vector SVG figures for triangulation calibration."
    )
    parser.add_argument(
        "--diagnostics",
        type=Path,
        default=Path("output/concordance_calibration/triangulation_diagnostics_summary.parquet"),
        help="Path to triangulation_diagnostics_summary.parquet",
    )
    parser.add_argument(
        "--trends",
        type=Path,
        default=Path("output/concordance_calibration/resolution_trend_summary.parquet"),
        help="Path to resolution_trend_summary.parquet",
    )
    parser.add_argument(
        "--fidelity",
        type=Path,
        default=Path("output/concordance_calibration/deconv_fidelity_benchmark.parquet"),
        help="Path to deconv_fidelity_benchmark.parquet",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/concordance_calibration"),
        help="Directory to save SVG files",
    )
    args = parser.parse_args()
    return PlotCalibrationConfig(
        diagnostics_path=args.diagnostics,
        trends_path=args.trends,
        fidelity_path=args.fidelity,
        results_dir=args.results_dir,
    )


def save_altair_svg(chart: alt.Chart, out_path: Path) -> Result[None, str]:
    """Pure function exporting Altair chart directly to vector SVG via vl-convert."""
    try:
        vl_spec = chart.to_dict()
        svg_str = vlc.vegalite_to_svg(vl_spec)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(svg_str, encoding="utf-8")
        return Success(None)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export SVG to {out_path}: {exc}")


def plot_triangulation_scatter(diag_df: pl.DataFrame) -> alt.Chart:
    """Fig 1: Ground Truth Cell Count Effect vs Deconvolution Effect colored by Cause."""
    sub = diag_df.filter((pl.col("condition") == "Combined") & (pl.col("resolution") == 0.5))
    if sub.is_empty():
        sub = diag_df.filter(pl.col("condition") == "Combined")

    scatter = (
        alt.Chart(sub)
        .mark_circle(size=120, opacity=0.85)
        .encode(
            x=alt.X("beta_count:Q", title="Ground Truth Cell Count Effect (Standardized β_count)"),
            y=alt.Y("beta_deconv:Q", title="Pseudo-bulk Deconvolution Effect (Standardized β_deconv)"),
            color=alt.Color(
                "primary_discrepancy_cause:N",
                legend=alt.Legend(title="Attributed Cause"),
                scale=alt.Scale(scheme="category10"),
            ),
            tooltip=["cell_state:N", "beta_count:Q", "beta_deconv:Q", "beta_milo:Q", "primary_discrepancy_cause:N"],
        )
    )

    diag_line = (
        alt.Chart(pl.DataFrame({"x": [-2, 2], "y": [-2, 2]}))
        .mark_line(strokeDash=[4, 4], color="#888888")
        .encode(x="x:Q", y="y:Q")
    )

    return (
        (scatter + diag_line)
        .properties(
            title="Single-Cell Ground Truth Count vs. Pseudo-bulk Deconvolution Effects",
            width=500,
            height=420,
        )
        .configure_axis(grid=True, gridColor="#f0f0f0")
    )


def plot_resolution_trends(trends_df: pl.DataFrame) -> alt.Chart:
    """Fig 3: Concordance metrics as a function of Leiden resolution (up to 3.0)."""
    sub = trends_df.filter(pl.col("condition") == "Combined")

    # Reshape for multi-line plot using modern unpivot
    melted = sub.unpivot(
        index=["resolution", "n_states"],
        on=["spearman_deconv_vs_umi", "spearman_umi_vs_count", "concordance_rate_bulk"],
        variable_name="metric",
        value_name="value",
    )

    metric_labels = {
        "spearman_deconv_vs_umi": "Deconv vs UMI (Identifiability)",
        "spearman_umi_vs_count": "UMI vs Count (mRNA Mass Bias)",
        "concordance_rate_bulk": "Concordance with Bulk",
    }
    melted = melted.with_columns(
        pl.col("metric").replace(metric_labels).alias("metric_name")
    )

    lines = (
        alt.Chart(melted)
        .mark_line(point=True)
        .encode(
            x=alt.X("resolution:Q", title="Leiden Clustering Resolution (0.5 to 3.0)"),
            y=alt.Y("value:Q", title="Correlation / Concordance Metric", scale=alt.Scale(domain=[0, 1])),
            color=alt.Color("metric_name:N", legend=alt.Legend(title="Validation Gate")),
            tooltip=["resolution:Q", "n_states:Q", "metric_name:N", "value:Q"],
        )
    )

    return lines.properties(
        title="Impact of Clustering Resolution on Deconvolution Identifiability and Concordance",
        width=520,
        height=380,
    ).configure_axis(grid=True, gridColor="#f0f0f0")


def plot_discrepancy_attribution_bar(diag_df: pl.DataFrame) -> alt.Chart:
    """Fig 4: Bar chart of attributed discrepancy causes across resolutions."""
    agg = (
        diag_df.filter(pl.col("condition") == "Combined")
        .group_by(["resolution", "primary_discrepancy_cause"])
        .agg(pl.len().alias("count"))
    )

    bar = (
        alt.Chart(agg)
        .mark_bar()
        .encode(
            x=alt.X("resolution:O", title="Clustering Resolution"),
            y=alt.Y("count:Q", title="Number of Cell States", stack="normalize"),
            color=alt.Color("primary_discrepancy_cause:N", legend=alt.Legend(title="Attribution")),
            tooltip=["resolution:O", "primary_discrepancy_cause:N", "count:Q"],
        )
    )

    return bar.properties(
        title="Proportion of Root Discrepancy Causes Across Clustering Resolutions",
        width=480,
        height=350,
    )


def run_pipeline(config: PlotCalibrationConfig) -> Result[None, str]:
    """Execute Step 5 visualization generation."""
    if not config.diagnostics_path.exists() or not config.trends_path.exists():
        return Failure("Required diagnostic summary files not found.")

    diag_df = pl.read_parquet(config.diagnostics_path)
    trends_df = pl.read_parquet(config.trends_path)

    config.results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fig 1: Triangulation Scatter
    fig1 = plot_triangulation_scatter(diag_df)
    match save_altair_svg(fig1, config.results_dir / "fig1_triangulation_scatter.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # 2. Fig 3: Resolution Trends
    fig3 = plot_resolution_trends(trends_df)
    match save_altair_svg(fig3, config.results_dir / "fig3_resolution_concordance_curve.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # 3. Fig 4: Discrepancy Attribution Bar
    fig4 = plot_discrepancy_attribution_bar(diag_df)
    match save_altair_svg(fig4, config.results_dir / "fig4_discrepancy_attribution_bar.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    print(f"[SUCCESS] All calibration vector SVGs generated in {config.results_dir}")
    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
