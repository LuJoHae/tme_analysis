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
from returns.result import Failure, Success

from tme_response.visualization import create_combined_summary_forest_plot, export_chart_svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate two-panel combined summary forest plot (ROC-AUC and PR-AUC)."
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc/iatlas_benchmark_metrics.parquet"),
        help="Path to iAtlas benchmark metrics parquet file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc"),
        help="Directory to save the generated SVG.",
    )
    parser.add_argument(
        "--brain-dir",
        type=Path,
        default=Path("/Users/halu/.gemini/antigravity/brain/8d054c32-dd8d-492f-aeba-2e31752acbd0"),
        help="Artifact directory to synchronize SVG and PNG.",
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

    if not args.parquet.exists():
        print(f"Error: Parquet file not found at {args.parquet}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading benchmark metrics from: {args.parquet}")
    raw_df = pl.read_parquet(args.parquet)

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

            # Copy to brain artifact directory
            if args.brain_dir.exists():
                dst_svg = args.brain_dir / args.filename
                shutil.copy2(svg_path, dst_svg)
                print(f"Copied to artifact directory: {dst_svg}")

                # Render 300 DPI PNG raster for validation
                dst_png = args.brain_dir / (args.filename.replace(".svg", ".png"))
                try:
                    subprocess.run(
                        ["rsvg-convert", "-d", "300", "-p", "300", str(dst_svg), "-o", str(dst_png)],
                        check=True,
                    )
                    print(f"Rendered 300 DPI raster: {dst_png}")
                except Exception as exc:
                    print(f"Warning: Failed to render PNG via rsvg-convert: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
