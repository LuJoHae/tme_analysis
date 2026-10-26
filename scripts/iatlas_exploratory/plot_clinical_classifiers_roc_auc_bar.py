#!/usr/bin/env python3
"""Plot a horizontal bar chart of ROC-AUC (±95% CI) for a curated clinical selection of ICB classifiers.

Reads `output/benchmarks/iatlas_roc_auc/iatlas_benchmark_metrics.parquet`, aggregates the
cross-cohort mean ROC-AUC (or standardized Pan-Cancer pool ROC-AUC) for a clinical selection
including TMB, checkpoint/infiltration IHC/RNA markers (CD274/PD-L1, PDCD1/PD-1, CD8A, CTLA4, CXCL9),
clinical trial signatures (Ayers GEP, Merck GEP, Fehrenbacher Teff, IFNG-6, CYT, TIDE), and
our deconvolution model, and exports a publication-grade vector SVG figure.

Adheres strictly to Functional Python Coding Style Rules:
- Pure functions and immutable frozen Pydantic models
- Result[T, str] and Maybe[T] monads for error and optionality handling
- Polars for declarative dataframe transformations
- Altair for declarative vector visualization exported to SVG
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
from typing import Protocol, Sequence, assert_never
import xml.etree.ElementTree as ET

import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from tme_response.visualization import (
    CLINICAL_DEFAULT_PREDICTORS,
    AggregationMode,
    ClinicalBarPlotConfig,
    aggregate_clinical_predictors,
    create_clinical_selection_bar_chart,
    export_chart_svg,
)


class ScriptConfig(BaseModel):
    """Immutable CLI configuration for the clinical classifier bar chart pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    parquet_path: Path
    output_svg_path: Path
    mode: AggregationMode
    predictors: tuple[str, ...]
    title: str
    subtitle: Maybe[str]
    width: int


class ParquetReader(Protocol):
    """Typeclass (Protocol) for reading a Parquet dataset into a Polars DataFrame."""

    def read(self, path: Path) -> Result[pl.DataFrame, str]: ...


class PolarsParquetReader(BaseModel):
    """Default immutable Parquet reader implementation."""

    model_config = ConfigDict(frozen=True)

    def read(self, path: Path) -> Result[pl.DataFrame, str]:
        return (
            Success(pl.read_parquet(path))
            if path.is_file()
            else Failure(f"Parquet file not found: {path}")
        )


def resolve_subtitle(mode: AggregationMode, custom_subtitle: Maybe[str]) -> Maybe[str]:
    """Pure function determining the subtitle based on aggregation mode and user override."""
    if isinstance(custom_subtitle, Some):
        return custom_subtitle
    match mode:
        case AggregationMode.COHORT_MEAN:
            return Some(
                "Unweighted mean ROC-AUC ± 95% CI across 9 pre-treatment clinical trials (dashed line = 0.50 random chance)"
            )
        case AggregationMode.PANCANCER_POOL:
            return Some(
                "Pan-cancer standardized pool ROC-AUC ± 95% CI across n = 1,015 patients (dashed line = 0.50 random chance)"
            )
        case _ as unreachable:
            assert_never(unreachable)


def validate_svg_file(svg_path: Path) -> Result[Path, str]:
    """Validate that the exported SVG file exists and parses cleanly as well-formed XML."""
    if not svg_path.is_file():
        return Failure(f"Exported SVG file does not exist: {svg_path}")
    try:
        ET.parse(svg_path)
        return Success(svg_path)
    except ET.ParseError as exc:
        return Failure(f"Invalid XML in exported SVG {svg_path}: {exc}")


def format_summary_table(summary_df: pl.DataFrame) -> str:
    """Pure formatter rendering the aggregated clinical classifier table for console inspection."""
    header = (
        f"{'Rank':<5} {'Classifier Display Name':<38} {'Category':<26} "
        f"{'ROC-AUC':>8} {'95% CI':>18} {'N':>4}"
    )
    rule = "-" * len(header)
    rows = tuple(
        f"{idx:<5} {str(row['predictor_display']):<38} {str(row['clinical_category']):<26} "
        f"{float(row['mean_roc_auc']):>8.4f} "
        f"[{float(row['ci_lower']):.3f}, {float(row['ci_upper']):.3f}] "
        f"{int(row['n_cohorts']):>4}"
        for idx, row in enumerate(summary_df.iter_rows(named=True), start=1)
    )
    return "\n".join((rule, header, rule, *rows, rule))


def run_bar_chart_pipeline(
    config: ScriptConfig,
    reader: ParquetReader,
) -> Result[tuple[Path, pl.DataFrame], str]:
    """Pure functional pipeline chaining data loading, aggregation, chart building, and SVG export."""
    plot_config = ClinicalBarPlotConfig(
        title=config.title,
        subtitle=resolve_subtitle(config.mode, config.subtitle),
        width=config.width,
    )
    return (
        reader.read(config.parquet_path)
        .bind(
            lambda raw_df: aggregate_clinical_predictors(
                raw_df,
                predictors=config.predictors,
                mode=config.mode,
            )
        )
        .bind(
            lambda summary_df: create_clinical_selection_bar_chart(
                summary_df,
                config=plot_config,
            )
            .bind(lambda chart: export_chart_svg(chart, config.output_svg_path))
            .bind(validate_svg_file)
            .map(lambda validated_path: (validated_path, summary_df))
        )
    )


def parse_predictors_arg(raw_predictors: str) -> tuple[str, ...]:
    """Parse a comma-separated string of predictor names or return the default clinical selection."""
    cleaned = tuple(p.strip() for p in raw_predictors.split(",") if p.strip())
    return cleaned if cleaned else CLINICAL_DEFAULT_PREDICTORS


def parse_cli_args(argv: Sequence[str]) -> ScriptConfig:
    """Parse command-line arguments into an immutable ScriptConfig."""
    parser = argparse.ArgumentParser(
        description="Generate a publication-grade horizontal ROC-AUC bar chart for clinical ICB classifiers."
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc/iatlas_benchmark_metrics.parquet"),
        help="Path to iAtlas benchmark metrics Parquet file.",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc/fig_clinical_classifiers_roc_auc_bar.svg"),
        help="Destination path for the exported vector SVG figure.",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=[m.value for m in AggregationMode],
        default=AggregationMode.COHORT_MEAN.value,
        help="Aggregation mode: 'cohort_mean' (default, 9-cohort unweighted mean) or 'pancancer_pool'.",
    )
    parser.add_argument(
        "--predictors",
        type=str,
        default="",
        help="Optional comma-separated list of predictor_name values (defaults to 13 clinical markers including TMB).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Clinical ICB Response Classifiers: Cross-Cohort ROC-AUC Comparison",
        help="Main title for the bar chart.",
    )
    parser.add_argument(
        "--subtitle",
        type=str,
        default="",
        help="Optional custom subtitle for the bar chart.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=380,
        help="Width of the bar chart canvas in pixels (default: 380).",
    )
    args = parser.parse_args(list(argv))

    subtitle_maybe: Maybe[str] = (
        Some(str(args.subtitle).strip())
        if str(args.subtitle).strip()
        else Nothing
    )

    return ScriptConfig(
        parquet_path=Path(args.parquet),
        output_svg_path=Path(args.output_svg),
        mode=AggregationMode(str(args.mode)),
        predictors=parse_predictors_arg(str(args.predictors)),
        title=str(args.title),
        subtitle=subtitle_maybe,
        width=int(args.width),
    )


def rasterize_svg_preview(svg_path: Path) -> Result[Path, str]:
    """Render a 300 DPI PNG raster alongside the SVG using rsvg-convert if available."""
    png_path = svg_path.with_suffix(".png")
    try:
        subprocess.run(
            ["rsvg-convert", "-d", "300", "-p", "300", str(svg_path), "-o", str(png_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return Success(png_path)
    except Exception as exc:
        return Failure(f"rsvg-convert rasterization skipped: {exc}")


def main(argv: Sequence[str] = ()) -> int:
    """Imperative shell entry point."""
    effective_argv = argv if argv else sys.argv[1:]
    config = parse_cli_args(effective_argv)
    reader = PolarsParquetReader()

    match run_bar_chart_pipeline(config, reader):
        case Failure(err):
            sys.stderr.write(f"[ERROR] Pipeline failed: {err}\n")
            return 1
        case Success((svg_path, summary_df)):
            sys.stdout.write(format_summary_table(summary_df) + "\n")
            sys.stdout.write(f"[OK] Exported vector SVG: {svg_path}\n")
            match rasterize_svg_preview(svg_path):
                case Success(png_path):
                    sys.stdout.write(f"[OK] Rendered 300 DPI PNG preview: {png_path}\n")
                case Failure(_):
                    pass
            return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
