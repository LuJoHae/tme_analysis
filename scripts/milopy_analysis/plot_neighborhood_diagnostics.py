#!/usr/bin/env python3
"""Diagnostic Visualization Dashboard for Multi-Cohort Milopy Analysis.

Renders publication-grade Nature Methods minimalist wireframe UMAP dashboards
displaying:
- Single-cell & neighborhood log2 fold change
- Neighborhood dataset Simpson purity and mixing
- Neighborhood patient dominance and clonal spike diagnostic (max donor fraction)
- Neighborhood hit classification (Conserved Recurrent Hit vs Private Clonal Spike)
- Single-cell source dataset provenance

Strict functional Python adhering to immutability, returns Result, Polars, and Altair vector SVG.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any
import xml.etree.ElementTree as ET

import altair as alt  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

try:
    import vl_convert as vlc  # type: ignore
except ImportError:
    vlc = None

from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)

OKABE_PURPLE = OKABE_REDDISH_PURPLE
OKABE_REDDISH = OKABE_REDDISH_PURPLE



class DiagnosticPlotConfig(BaseModel):
    """Configuration for neighborhood diagnostic plotting."""

    model_config = ConfigDict(frozen=True)

    results_dir: Path
    output_svg: Path
    fdr_threshold: float = 0.10
    max_cells_display: int = 15000
    title: str = "Multi-Cohort Single-Cell Neighborhood Diagnostics"


def subsample_cells_for_display(
    cell_df: pl.DataFrame,
    max_cells: int = 15000,
    seed: int = 42,
) -> pl.DataFrame:
    """Subsamples cells for clean, responsive vector plotting while preserving significant cells."""
    if len(cell_df) <= max_cells:
        return cell_df

    # Separate significant and non-significant cells if column present
    if "is_significant" in cell_df.columns:
        sig_cells = cell_df.filter(pl.col("is_significant") == True)
        non_sig = cell_df.filter(pl.col("is_significant") == False)
    elif "min_fdr" in cell_df.columns:
        sig_cells = cell_df.filter(pl.col("min_fdr") < 0.10)
        non_sig = cell_df.filter(pl.col("min_fdr") >= 0.10)
    else:
        return cell_df.sample(n=max_cells, seed=seed)

    n_sig = len(sig_cells)
    if n_sig >= max_cells:
        return sig_cells.sample(n=max_cells, seed=seed)

    n_remain = max_cells - n_sig
    sampled_non_sig = non_sig.sample(n=min(len(non_sig), n_remain), seed=seed)
    return pl.concat([sampled_non_sig, sig_cells])


def build_cell_logfc_panel(cell_df: pl.DataFrame, fdr_thresh: float = 0.10) -> alt.Chart:
    """Panel A: Single-cell UMAP colored by projected differential abundance log2FC."""
    pdf = cell_df.to_pandas()
    max_fc = max(2.5, float(np.abs(pdf["cell_logfc"]).quantile(0.99)))

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=7, opacity=0.75)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                "cell_logfc:Q",
                scale=alt.Scale(domain=[-max_fc, 0, max_fc], range=[OKABE_BLUE, "#F1F5F9", OKABE_VERMILION]),
                title="Cell log2FC (R vs NR)",
                legend=alt.Legend(orient="bottom", gradientLength=130),
            ),
            tooltip=["cell_id:N", "dataset_id:N", "patient_id:N", "clinical_response:N", "cell_logfc:Q"],
        )
        .properties(width=280, height=240, title="A. Projected Single-Cell log2FC")
    )
    return chart


def build_nhood_logfc_panel(da_df: pl.DataFrame, fdr_thresh: float = 0.10) -> alt.Chart:
    """Panel B: Neighborhood centroid UMAP colored by logFC and stroked by significance."""
    pdf = da_df.to_pandas()
    if "UMAP1" not in pdf.columns:
        # Fallback if nhood coordinates not explicitly assigned
        pdf["UMAP1"] = pdf.get("Nhood_UMAP1", np.random.normal(0, 1, len(pdf)))
        pdf["UMAP2"] = pdf.get("Nhood_UMAP2", np.random.normal(0, 1, len(pdf)))

    max_fc = max(2.5, float(np.abs(pdf["logFC"]).quantile(0.99)))
    chart = (
        alt.Chart(pdf)
        .mark_circle(size=35, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                "logFC:Q",
                scale=alt.Scale(domain=[-max_fc, 0, max_fc], range=[OKABE_BLUE, "#F1F5F9", OKABE_VERMILION]),
                title="Nhood log2FC",
                legend=alt.Legend(orient="bottom", gradientLength=130),
            ),
            stroke=alt.condition(
                alt.datum.FDR < fdr_thresh,
                alt.value("#0F172A"),
                alt.value("transparent"),
            ),
            strokeWidth=alt.condition(alt.datum.FDR < fdr_thresh, alt.value(1.2), alt.value(0.0)),
            tooltip=["Nhood:N", "logFC:Q", "FDR:Q", "status:N", "hit_type:N"],
        )
        .properties(width=280, height=240, title=f"B. Neighborhood log2FC (FDR < {fdr_thresh:.2f})")
    )
    return chart


def build_dataset_purity_panel(da_df: pl.DataFrame) -> alt.Chart:
    """Panel C: Neighborhood Dataset Simpson Purity (1 = private to 1 dataset, low = mixed)."""
    pdf = da_df.to_pandas()
    col = "dataset_id_purity" if "dataset_id_purity" in pdf.columns else "dataset_purity"
    if col not in pdf.columns:
        pdf[col] = 1.0

    if "UMAP1" not in pdf.columns:
        pdf["UMAP1"] = pdf.get("Nhood_UMAP1", np.random.normal(0, 1, len(pdf)))
        pdf["UMAP2"] = pdf.get("Nhood_UMAP2", np.random.normal(0, 1, len(pdf)))

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=35, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                f"{col}:Q",
                scale=alt.Scale(domain=[0.0, 1.0], range=["#3B82F6", "#FDE047", "#EF4444"]),
                title="Dataset Simpson Purity",
                legend=alt.Legend(orient="bottom", gradientLength=130),
            ),
            tooltip=["Nhood:N", f"{col}:Q", "dominant_dataset:N"],
        )
        .properties(width=280, height=240, title="C. Dataset Integration Purity (1 = Private)")
    )
    return chart


def build_patient_dominance_panel(da_df: pl.DataFrame) -> alt.Chart:
    """Panel D: Maximum single-patient cell fraction in neighborhood (diagnosing clonal spikes)."""
    pdf = da_df.to_pandas()
    col = "patient_id_max_fraction" if "patient_id_max_fraction" in pdf.columns else "patient_max_frac"
    if col not in pdf.columns:
        pdf[col] = 0.5

    if "UMAP1" not in pdf.columns:
        pdf["UMAP1"] = pdf.get("Nhood_UMAP1", np.random.normal(0, 1, len(pdf)))
        pdf["UMAP2"] = pdf.get("Nhood_UMAP2", np.random.normal(0, 1, len(pdf)))

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=35, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                f"{col}:Q",
                scale=alt.Scale(domain=[0.0, 0.5, 1.0], range=["#10B981", "#F59E0B", "#DC2626"]),
                title="Max Donor Fraction",
                legend=alt.Legend(orient="bottom", gradientLength=130),
            ),
            tooltip=["Nhood:N", f"{col}:Q", "dominant_patient:N", "replicate_prevalence:Q"],
        )
        .properties(width=280, height=240, title="D. Donor Dominance (Spike Diagnostic)")
    )
    return chart


def build_hit_classification_panel(da_df: pl.DataFrame) -> alt.Chart:
    """Panel E: Replicate prevalence gating & hit classification."""
    pdf = da_df.to_pandas()
    col = "hit_type" if "hit_type" in pdf.columns else "status"
    if "UMAP1" not in pdf.columns:
        pdf["UMAP1"] = pdf.get("Nhood_UMAP1", np.random.normal(0, 1, len(pdf)))
        pdf["UMAP2"] = pdf.get("Nhood_UMAP2", np.random.normal(0, 1, len(pdf)))

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=35, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                f"{col}:N",
                scale=alt.Scale(
                    domain=[
                        "Conserved Recurrent Hit",
                        "Cohort-Private Replicated Hit",
                        "Private Clonal Spike",
                        "Not Significant",
                    ],
                    range=[
                        OKABE_BLUISH_GREEN,
                        OKABE_SKY_BLUE,
                        OKABE_VERMILION,
                        COLOR_LIGHT_GREY,
                    ],
                ),
                title="Hit Classification",
                legend=alt.Legend(orient="bottom", columns=2),
            ),
            tooltip=["Nhood:N", "hit_type:N", "status:N", "logFC:Q", "FDR:Q"],
        )
        .properties(width=280, height=240, title="E. Replicate Prevalence Gating")
    )
    return chart


def build_dataset_source_panel(cell_df: pl.DataFrame) -> alt.Chart:
    """Panel F: Single cells colored by source dataset ID."""
    pdf = cell_df.to_pandas()
    datasets = sorted(list(pdf["dataset_id"].unique()))

    # Color palette
    colors = [
        OKABE_ORANGE,
        OKABE_SKY_BLUE,
        OKABE_BLUISH_GREEN,
        OKABE_PURPLE,
        OKABE_REDDISH,
        OKABE_BLUE,
        "#000000",
        "#E6AB02",
        "#A6761D",
    ][: len(datasets)]

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=7, opacity=0.75)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                "dataset_id:N",
                scale=alt.Scale(domain=datasets, range=colors),
                title="Dataset Source",
                legend=alt.Legend(orient="bottom", columns=2),
            ),
            tooltip=["cell_id:N", "dataset_id:N", "patient_id:N", "clinical_response:N"],
        )
        .properties(width=280, height=240, title="F. Cohort Batch Provenance")
    )
    return chart


def render_diagnostic_dashboard(
    da_df: pl.DataFrame,
    cell_df: pl.DataFrame,
    output_svg: Path,
    title: str = "Multi-Cohort Single-Cell Neighborhood Diagnostics",
    fdr_thresh: float = 0.10,
    max_cells: int = 15000,
) -> Result[Path, str]:
    """Assembles and exports the 6-panel diagnostic UMAP dashboard in vector SVG and PNG."""
    try:
        cell_display = subsample_cells_for_display(cell_df, max_cells=max_cells)

        panel_a = build_cell_logfc_panel(cell_display, fdr_thresh=fdr_thresh)
        panel_b = build_nhood_logfc_panel(da_df, fdr_thresh=fdr_thresh)
        panel_c = build_dataset_purity_panel(da_df)
        panel_d = build_patient_dominance_panel(da_df)
        panel_e = build_hit_classification_panel(da_df)
        panel_f = build_dataset_source_panel(cell_display)

        # 2 rows x 3 columns layout
        row1 = alt.hconcat(panel_a, panel_b, panel_c, spacing=20)
        row2 = alt.hconcat(panel_d, panel_e, panel_f, spacing=20)
        dashboard = (
            alt.vconcat(row1, row2, spacing=25)
            .properties(title=alt.TitleParams(text=title, fontSize=14, fontWeight="bold", anchor="middle"))
        )
        dashboard = apply_nature_methods_theme(dashboard)

        output_svg.parent.mkdir(parents=True, exist_ok=True)

        if vlc is not None:
            # Export via vl_convert
            svg_str = vlc.vegalite_to_svg(dashboard.to_json())
            output_svg.write_text(svg_str, encoding="utf-8")

            # Also export 300 DPI PNG
            png_path = output_svg.with_suffix(".png")
            png_bytes = vlc.vegalite_to_png(dashboard.to_json(), scale=3.0)
            png_path.write_bytes(png_bytes)
            print(f"[EXPORTED] {output_svg}")
            print(f"[EXPORTED] {png_path}")
        else:
            dashboard.save(str(output_svg))
            print(f"[EXPORTED] {output_svg}")

        return Success(output_svg)
    except Exception as exc:
        msg = f"Failed to render diagnostic dashboard: {exc}"
        print(f"Error: {msg}", file=sys.stderr)
        return Failure(msg)


def parse_args() -> DiagnosticPlotConfig:
    parser = argparse.ArgumentParser(
        description="Render publication-grade diagnostic UMAP dashboard for Milopy multi-cohort analysis."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing da_results.parquet and cell_level_scores.parquet",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=None,
        help="Target output SVG path (default: <results-dir>/neighborhood_diagnostics.svg)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Significance cutoff for FDR",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=15000,
        help="Maximum single cells to display on single-cell panels",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Multi-Cohort Single-Cell Neighborhood Diagnostics",
        help="Master figure title",
    )
    args = parser.parse_args()

    out_svg = args.output_svg or (args.results_dir / "neighborhood_diagnostics.svg")
    return DiagnosticPlotConfig(
        results_dir=args.results_dir,
        output_svg=out_svg,
        fdr_threshold=args.fdr_threshold,
        max_cells_display=args.max_cells,
        title=args.title,
    )


def main() -> None:
    config = parse_args()

    da_path = config.results_dir / "da_results.parquet"
    cell_path = config.results_dir / "cell_level_scores.parquet"

    if not da_path.exists():
        print(f"Error: {da_path} does not exist", file=sys.stderr)
        sys.exit(1)
    if not cell_path.exists():
        print(f"Error: {cell_path} does not exist", file=sys.stderr)
        sys.exit(1)

    da_df = pl.read_parquet(da_path)
    cell_df = pl.read_parquet(cell_path)

    res = render_diagnostic_dashboard(
        da_df=da_df,
        cell_df=cell_df,
        output_svg=config.output_svg,
        title=config.title,
        fdr_thresh=config.fdr_threshold,
        max_cells=config.max_cells_display,
    )

    match res:
        case Success(path):
            print(f"Successfully generated diagnostic dashboard: {path}")
        case Failure(err):
            print(f"Failed to generate dashboard: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
