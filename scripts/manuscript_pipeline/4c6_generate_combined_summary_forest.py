#!/usr/bin/env python3
"""Generate a publication-grade combined two-panel meta-analytic forest summary plot.

Combines ROC-AUC discrimination and PR-AUC precision into a single figure with a shared
y-axis (predictor labels rendered only once on the left panel).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.manuscript_pipeline._helpers.paths import resolve_upstream_artifact  # noqa: E402
from tme_response.visualization import create_combined_summary_forest_plot, export_chart_svg  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate two-panel combined summary forest plot (ROC-AUC and PR-AUC)."
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path(
            "output/manuscript_pipeline/3e5_integrate_deconvolution_our/iatlas_benchmark_metrics.parquet"
        ),
        help="Path to iAtlas 32-biomarker benchmark metrics parquet file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/4c6_generate_combined_summary_forest"),
        help="Directory to save the generated SVG and 300 DPI PNG.",
    )
    parser.add_argument(
        "--sort",
        type=str,
        default="roc_auc",
        choices=["roc_auc", "pr_auc", "alphabetical"],
        help="Ordering strategy for predictors along the shared y-axis.",
    )
    parser.add_argument(
        "--panel-width",
        type=int,
        default=180,
        help="Width in pixels for each individual panel.",
    )
    parser.add_argument(
        "--filename",
        type=str,
        default="fig_table_s2_combined_forest_summary.svg",
        help="Base filename for the output SVG.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    parquet_path = args.parquet
    if not parquet_path.exists():
        match resolve_upstream_artifact(
            "iatlas_benchmark_metrics.parquet",
            ("3e5_integrate_deconvolution_our", "3e4_update_texh_benchmark", "3e1_benchmark_iatlas_icb_predictors"),
            fallback_dirs=(Path("output/benchmarks/iatlas_roc_auc"),),
        ):
            case Some(found):
                parquet_path = found
            case Nothing:
                print(f"Error: Parquet file not found at {args.parquet}", file=sys.stderr)
                sys.exit(1)

    print(f"Loading benchmark metrics from: {parquet_path}")
    raw_df = pl.read_parquet(parquet_path)

    # Filter to pan-cancer standardized pre-treatment benchmark stratum
    pancan_df = raw_df.filter(
        (pl.col("cohort_id") == "PanCancer-Combined")
        & (pl.col("time_stratum") == "Pre")
        & (pl.col("response_stratum") == "Standard")
        & (pl.col("pooling_strategy") == "standardized")
    )

    if pancan_df.is_empty():
        print("Error: No PanCancer-Combined Pre_Standard standardized records found.", file=sys.stderr)
        sys.exit(1)

    print(f"Creating combined two-panel forest plot (sort_by={args.sort}, panel_width={args.panel_width})...")
    plot_result = create_combined_summary_forest_plot(
        pancan_df,
        title="Pan-Cancer Meta-Analytic ICB Response Performance Across 32 Biomarkers",
        panel_width=args.panel_width,
        roc_domain=[0.45, 0.80],
        pr_domain=[0.25, 0.60],
        sort_by=args.sort,
        spacing=15,
    )

    match plot_result:
        case Failure(err):
            print(f"Plot generation failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(chart):
            args.output_dir.mkdir(parents=True, exist_ok=True)
            svg_path = args.output_dir / args.filename

            match export_chart_svg(chart, svg_path):
                case Failure(err):
                    print(f"Failed to export SVG: {err}", file=sys.stderr)
                    sys.exit(1)
                case Success(saved_path):
                    print(f"Successfully saved SVG: {saved_path}")

            png_path = svg_path.with_suffix(".png")
            from plotting_utils import rasterize_svg_to_png

            match rasterize_svg_to_png(svg_path, png_path, scale=3.125, background="white"):
                case Success(_):
                    print(f"Rendered 300 DPI raster: {png_path}")
                case Failure(err):
                    print(f"Warning: Failed to render PNG: {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
