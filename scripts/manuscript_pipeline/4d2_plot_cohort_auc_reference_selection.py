#!/usr/bin/env python3
"""
Selected Predictive Capacity Subplots from Sade-Feldman Multi-Resolution Deconvolution Benchmark.

Derived from scripts/manuscript_pipeline/4d1_plot_sf_multires_validation_figures.py.
Plots exclusively in a single 1-row layout:
  - Subplot A: Mean Cohort AUC by Resolution
  - Subplot B: Mean Cohort AUC vs. Number of Cell States
  - Subplot C: Mean Cohort AUC vs. Number of Integrated Datasets (filtered to n_datasets in {1, 2, 3, 16})

Strict functional Python: Polars, Altair, returns.result, Pydantic (frozen=True).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Final

import altair as alt  # type: ignore[import-untyped]
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.manuscript_pipeline._helpers.paths import (  # noqa: E402
    MULTIRES_ALL_UPSTREAM_STAGES,
    resolve_upstream_artifact,
)

# Disable row limit for single-cell data rendering
alt.data_transformers.disable_max_rows()


class PlotSelectionConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    data_dir: Path
    results_dir: Path
    output_filename: str = "fig1_selected_cohort_auc_and_marginal_gain"
    marginal_metric: str = "melanoma"
    panel_width: int = 300
    panel_height: int = 240
    main_color: str = "#D55E00"  # Okabe-Ito Vermillion for distinct, high-contrast trend line


BOX_CONFIG: Final[dict[str, Any]] = {
    "box": {"fill": "#E2E8F0", "stroke": "#475569", "strokeWidth": 1.2},
    "median": {"stroke": "#0F172A", "strokeWidth": 2.0},
    "rule": {"stroke": "#64748B", "strokeWidth": 1.2},
    "ticks": False,
    "outliers": False,
}


REF_METADATA: Final[dict[str, dict[str, Any]]] = {
    "Sade-Feldman": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 16288},
    "Jerby-Arnon": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 7186},
    "Maynard-NSCLC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3000},
    "Ma-Liver": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 5115},
    "Yost-BCC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3500},
    "Combined-Atlas": {"category": "Criteria-Combined", "n_datasets": 16, "n_cells": 41284},
    "Melanoma-Duo": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10686},
    "SS2-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10186},
    "10x-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8615},
    "Cross-Tissue-Pair": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8115},
    "Triple-ICI": {"category": "Criteria-Combined", "n_datasets": 3, "n_cells": 13686},
    "Random-Pair-1": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 12301},
    "Random-Pair-2": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 6500},
    "Random-Pair-3": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 10686},
    "Random-Triplet-1": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15301},
    "Random-Triplet-2": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 11615},
    "Random-Triplet-3": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15801},
    "All-Datasets-Combined": {"category": "Criteria-Combined", "n_datasets": 4, "n_cells": 18801},
    "Random-Quadruplet-1": {"category": "Random-Combined", "n_datasets": 4, "n_cells": 18801},
}


def load_and_enrich_benchmark_data(data_dir: Path) -> Result[pl.DataFrame, str]:
    """Load benchmark summary parquet and enrich metadata attributes."""
    fallback_dirs = (
        Path("output/output/sade_feldman_deconv_validation"),
        Path("output/sade_feldman_deconv_validation"),
    )
    match resolve_upstream_artifact(
        "multi_resolution_benchmark_summary.parquet",
        MULTIRES_ALL_UPSTREAM_STAGES,
        explicit_dirs=(data_dir,),
        fallback_dirs=fallback_dirs,
    ):
        case Some(found):
            bench_file = found
        case Nothing:
            return Failure(
                f"Benchmark summary not found in {data_dir} or upstream stages {MULTIRES_ALL_UPSTREAM_STAGES}"
            )

    try:
        df = pl.read_parquet(bench_file)

        # Standardize criteria combination names to canonical publication terminology
        canonical_map = {
            "Cross-Tissue-Contrast": "Cross-Tissue-Pair",
            "Platform-SS2-Duo": "SS2-Cross-Cancer",
            "Platform-10x-Duo": "10x-Cross-Cancer",
            "ICI-Trio": "Triple-ICI",
            "Random-Quad-1": "Random-Quadruplet-1",
        }
        df = df.with_columns(
            pl.col("reference_type").replace(canonical_map)
        )

        if "subsample_fraction" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("subsample_fraction"))

        # Enrich metadata only if missing or zero
        df_enriched = df.with_columns(
            pl.when(pl.col("n_datasets").is_not_null() & (pl.col("n_datasets") > 0))
            .then(pl.col("n_datasets"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: REF_METADATA.get(str(r), {}).get("n_datasets", 1),
                    return_dtype=pl.Int64,
                )
            )
            .alias("n_datasets"),
            pl.when(pl.col("n_cells").is_not_null() & (pl.col("n_cells") > 0))
            .then(pl.col("n_cells"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: REF_METADATA.get(str(r), {}).get("n_cells", 0),
                    return_dtype=pl.Int64,
                )
            )
            .alias("n_cells"),
            pl.when(pl.col("reference_category").is_not_null() & (pl.col("reference_category") != "Unknown"))
            .then(pl.col("reference_category"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: str(REF_METADATA.get(str(r), {}).get("category", "Single Dataset")),
                    return_dtype=pl.String,
                )
            )
            .alias("reference_category"),
        )
        return Success(df_enriched)
    except Exception as exc:
        return Failure(f"Failed to load and enrich benchmark data: {exc}")


def build_selected_predictive_figure(
    df: pl.DataFrame,
    config: PlotSelectionConfig,
) -> Result[alt.VConcatChart, str]:
    """Construct the 4-panel figure containing subplots C, F, I (excluding n=4), and R."""
    try:
        # Full-cell configurations (exclude subsampling ladders)
        df_full = df.filter(
            (pl.col("subsample_fraction") >= 0.999)
            & (~pl.col("reference_type").str.contains(r"\(100% Cells\)"))
        )
        df_full_pd = df_full.to_pandas()

        color_scale = alt.Scale(scheme="tableau20")
        dash_scale = alt.Scale(
            domain=["Single Dataset", "Criteria-Combined", "Random-Combined"],
            range=[[1, 0], [6, 2], [2, 4]],
        )

        # Ensure sorted order for clean line trajectories
        df_full = df_full.sort(["reference_type", "resolution"])
        df_full_binned = df_full.with_columns(
            ((pl.col("n_clusters") // 5) * 5 + 2.5).alias("cluster_bin_center")
        )
        df_full_pd = df_full.to_pandas()
        df_binned_pd = df_full_binned.to_pandas()

        color_scale = alt.Scale(scheme="tableau20")

        base_a = alt.Chart(df_full_pd)

        panel_w = config.panel_width
        panel_h = config.panel_height

        # -----------------------------------------------------------------------
        # Subplot A: Mean Cohort AUC by Resolution (formerly C)
        # -----------------------------------------------------------------------
        box_a = base_a.mark_boxplot(size=26, **BOX_CONFIG).encode(
            x=alt.X(
                "resolution:O",
                title="Leiden Resolution",
                axis=alt.Axis(labelAngle=0),
            ),
            y=alt.Y(
                "mean_cohort_multivariate_auc:Q",
                title="Mean Cohort AUC",
                scale=alt.Scale(domain=[0.58, 0.78]),
            ),
        )
        lines_a = base_a.mark_line(strokeWidth=1.0, opacity=0.4, color="#94A3B8").encode(
            x=alt.X("resolution:O"),
            y=alt.Y("mean_cohort_multivariate_auc:Q"),
            detail="reference_type:N",
        )
        pts_a = base_a.mark_circle(size=16, color="#94A3B8", opacity=0.35).encode(
            x=alt.X("resolution:O"),
            y=alt.Y("mean_cohort_multivariate_auc:Q"),
            tooltip=[
                "reference_type",
                "resolution",
                "mean_cohort_multivariate_auc",
            ],
        )
        df_mean_res = (
            df_full.group_by("resolution")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_auc"))
            .sort("resolution")
            .to_pandas()
        )
        line_a = alt.Chart(df_mean_res).mark_line(
            color=config.main_color,
            strokeWidth=2.6,
            strokeDash=[4, 3],
            point=alt.OverlayMarkDef(color=config.main_color, size=65, filled=True),
        ).encode(
            x=alt.X("resolution:O"),
            y=alt.Y("mean_auc:Q"),
        )
        panel_a = alt.layer(box_a, lines_a, pts_a, line_a).properties(
            title="A. Mean Cohort AUC by Resolution",
            width=panel_w,
            height=panel_h,
        )

        # -----------------------------------------------------------------------
        # Subplot B: Mean Cohort AUC vs. Number of Cell States (formerly F)
        # -----------------------------------------------------------------------
        base_b = alt.Chart(df_full_pd)

        box_b = alt.Chart(df_binned_pd).mark_boxplot(size=20, **BOX_CONFIG).encode(
            x=alt.X(
                "cluster_bin_center:Q",
                title="Reference Cell States (Clusters)",
                scale=alt.Scale(domain=[4, 55], zero=False),
                axis=alt.Axis(tickMinStep=5),
            ),
            y=alt.Y(
                "mean_cohort_multivariate_auc:Q",
                title="Mean Cohort AUC",
                scale=alt.Scale(domain=[0.58, 0.78]),
            ),
        )
        lines_b = base_b.mark_line(strokeWidth=1.0, opacity=0.4, color="#94A3B8").encode(
            x=alt.X("n_clusters:Q", scale=alt.Scale(domain=[4, 55], zero=False)),
            y=alt.Y("mean_cohort_multivariate_auc:Q"),
            detail="reference_type:N",
        )
        pts_b = base_b.mark_circle(size=16, color="#94A3B8", opacity=0.35).encode(
            x=alt.X("n_clusters:Q"),
            y=alt.Y("mean_cohort_multivariate_auc:Q"),
            tooltip=[
                "reference_type",
                "n_clusters",
                "resolution",
                "mean_cohort_multivariate_auc",
            ],
        )
        df_b_binned = (
            df_full_binned.group_by("cluster_bin_center")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_auc"))
            .sort("cluster_bin_center")
            .to_pandas()
        )
        line_b = alt.Chart(df_b_binned).mark_line(
            color=config.main_color,
            strokeWidth=2.6,
            strokeDash=[4, 3],
            point=alt.OverlayMarkDef(color=config.main_color, size=65, filled=True),
        ).encode(
            x=alt.X("cluster_bin_center:Q"),
            y=alt.Y("mean_auc:Q"),
        )
        panel_b = alt.layer(box_b, lines_b, pts_b, line_b).properties(
            title="B. Mean Cohort AUC vs. Number of Cell States",
            width=panel_w,
            height=panel_h,
        )

        # -----------------------------------------------------------------------
        # Subplot C: Mean Cohort AUC vs. Number of Integrated Datasets (formerly I)
        # -----------------------------------------------------------------------
        # User requirement: "leave out dataset number equals 4 (so only 1 2 3 and 16)"
        allowed_datasets = [1, 2, 3, 16]
        df_filtered_ds = df_full.filter(pl.col("n_datasets").is_in(allowed_datasets))
        df_filtered_ds_pd = df_filtered_ds.to_pandas()

        df_ds_summary = (
            df_filtered_ds.group_by("n_datasets")
            .agg(
                pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"),
            )
            .sort("n_datasets")
            .to_pandas()
        )

        base_c = alt.Chart(df_filtered_ds_pd)
        base_ds = alt.Chart(df_ds_summary)

        panel_c = (
            alt.layer(
                base_c.mark_boxplot(size=30, **BOX_CONFIG).encode(
                    x=alt.X(
                        "n_datasets:O",
                        title="Number of Integrated Datasets",
                        sort=allowed_datasets,
                        axis=alt.Axis(labelAngle=0),
                    ),
                    y=alt.Y(
                        "mean_cohort_multivariate_auc:Q",
                        title="Mean Cohort AUC",
                        scale=alt.Scale(domain=[0.58, 0.78]),
                    ),
                ),
                base_c.transform_calculate(
                    jitter="0.5 + (random() - 0.5) * 0.4"
                ).mark_circle(size=36, opacity=0.65).encode(
                    x=alt.X(
                        "n_datasets:O",
                        title="Number of Integrated Datasets",
                        sort=allowed_datasets,
                    ),
                    xOffset=alt.XOffset("jitter:Q", scale=alt.Scale(domain=[0, 1])),
                    y=alt.Y(
                        "mean_cohort_multivariate_auc:Q",
                        scale=alt.Scale(domain=[0.58, 0.78]),
                    ),
                    color=alt.Color(
                        "reference_type:N",
                        title="Reference",
                        scale=color_scale,
                        legend=alt.Legend(columns=2, orient="right", titleFontWeight="bold"),
                    ),
                    tooltip=[
                        "reference_type",
                        "n_datasets",
                        "resolution",
                        "mean_cohort_multivariate_auc",
                    ],
                ),
                base_ds.mark_line(
                    color=config.main_color,
                    strokeWidth=2.6,
                    strokeDash=[4, 3],
                    point=alt.OverlayMarkDef(color=config.main_color, size=65, filled=True),
                ).encode(
                    x=alt.X("n_datasets:O", sort=allowed_datasets),
                    y=alt.Y("mean_coh_auc:Q"),
                ),
            )
            .properties(
                title="C. Mean Cohort AUC vs. Number of Integrated Datasets",
                width=panel_w,
                height=panel_h,
            )
        )

        # -----------------------------------------------------------------------
        # 1-Row Horizontal Composition (Subplots A, B, C)
        # -----------------------------------------------------------------------
        final_chart = (
            (panel_a | panel_b | panel_c)
            .resolve_scale(color="independent")
            .properties(
                title="Predictive Capacity Benchmark: Multi-Resolution Reference Selection (Subplots A–C)"
            )
            .configure_title(fontSize=15, anchor="start", font="Helvetica")
            .configure_axis(labelFontSize=10, titleFontSize=11, titleFontWeight="bold")
        )

        return Success(final_chart)
    except Exception as exc:
        return Failure(f"Failed to build selected predictive figure: {exc}")


def export_selected_figure(
    chart: alt.HConcatChart,
    config: PlotSelectionConfig,
) -> Result[tuple[Path, Path], str]:
    """Export chart to SVG and 300 DPI PNG."""
    try:
        config.results_dir.mkdir(parents=True, exist_ok=True)
        svg_path = config.results_dir / f"{config.output_filename}.svg"
        png_path = config.results_dir / f"{config.output_filename}.png"

        chart.save(str(svg_path))

        from plotting_utils import rasterize_svg_to_png

        match rasterize_svg_to_png(svg_path, png_path, scale=3.125, background="white"):
            case Success(_):
                pass
            case Failure(_):
                try:
                    chart.save(str(png_path), scale_factor=2.0)
                except Exception:
                    pass

        return Success((svg_path, png_path))
    except Exception as exc:
        return Failure(f"Failed to export figures: {exc}")


def run_pipeline(config: PlotSelectionConfig) -> Result[tuple[Path, Path], str]:
    """Execute complete functional pipeline: Load -> Build -> Export."""
    match load_and_enrich_benchmark_data(config.data_dir):
        case Failure(err):
            return Failure(err)
        case Success(df):
            match build_selected_predictive_figure(df, config):
                case Failure(err):
                    return Failure(err)
                case Success(chart):
                    return export_selected_figure(chart, config)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate 4-subplot predictive capacity benchmark (Subplots C, F, I without n=4, and R)."
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="output/manuscript_pipeline/3f8_compare_concordance_multires",
        help="Directory containing benchmark summary parquet file",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="output/manuscript_pipeline/4d2_plot_cohort_auc_reference_selection",
        help="Directory to save generated SVG and PNG figures",
    )
    parser.add_argument(
        "--output-filename",
        type=str,
        default="fig1_selected_cohort_auc_and_marginal_gain",
        help="Base name for exported SVG/PNG figures (without extension)",
    )
    parser.add_argument(
        "--marginal-metric",
        type=str,
        choices=["melanoma", "mean_cohort"],
        default="melanoma",
        help="Metric to compute marginal resolution gain: 'melanoma' (original line 1639) or 'mean_cohort'",
    )
    parser.add_argument(
        "--main-color",
        type=str,
        default="#D55E00",
        help="Hex color for main empirical mean trend line and overlay points (default: Okabe-Ito Vermillion #D55E00)",
    )
    args = parser.parse_args()

    config = PlotSelectionConfig(
        data_dir=Path(args.data_dir),
        results_dir=Path(args.results_dir),
        output_filename=args.output_filename,
        marginal_metric=args.marginal_metric,
        main_color=args.main_color,
    )

    match run_pipeline(config):
        case Success((svg_path, png_path)):
            print(f"[OK] Successfully generated:")
            print(f"     SVG: {svg_path}")
            print(f"     PNG: {png_path}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] Pipeline failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
