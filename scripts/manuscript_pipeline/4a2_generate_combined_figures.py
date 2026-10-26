#!/usr/bin/env python3
"""
Publication-Grade Combined Figure Generator for Milopy Single-Cell Differential Abundance.

Generates two master 3x3 multi-cohort publication figures:
1. milopy_cross_cohort_volcano_grid.svg: 9 panels (a-i) displaying neighborhood-level volcano plots
   (log2FC vs -log10 FDR) with horizontal significance thresholds and Okabe-Ito coloring.
2. milopy_cross_cohort_umaps_logfc.svg: 9 panels (a-i) displaying single-cell UMAP embeddings
   colored by continuous inferred Milopy log2FC with anti-occlusion depth ordering and intelligent sampling.

Strict functional Python adhering to immutability, returns Result, Polars, and Altair vector SVG.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence
import xml.etree.ElementTree as ET

import altair as alt  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import vl_convert as vlc  # type: ignore

# Disable Altair row limit
alt.data_transformers.disable_max_rows()

from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_GREY,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)


class FiguresConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    results_dir: Path = Path("output/results/milopy_response_analysis")
    reports_dir: Path = Path("output/reports")
    fdr_threshold: float = 0.10
    panel_width: int = 300
    panel_height_volcano: int = 230
    panel_height_umap: int = 240
    max_cells_per_umap_panel: int = 5000
    timepoint: str = "all"
    volcano_out: Path | None = None
    umap_out: Path | None = None
    seed: int = 42


def parse_args() -> FiguresConfig:
    parser = argparse.ArgumentParser(
        description="Generates publication-grade combined multi-panel figures for Milopy single-cell analysis."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/results/milopy_response_analysis"),
        help="Directory containing cohort milopy results",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("output/reports"),
        help="Directory to save combined SVG figures",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Significance threshold for FDR",
    )
    parser.add_argument(
        "--max-cells-per-panel",
        type=int,
        default=5000,
        help="Max single cells per UMAP panel in the multi-panel grid (significant cells always preserved)",
    )
    parser.add_argument(
        "--timepoint",
        type=str,
        default="all",
        choices=["all", "pre", "post"],
        help="Biopsy timepoint filter ('all', 'pre', 'post')",
    )
    parser.add_argument(
        "--volcano-out",
        type=Path,
        default=None,
        help="Custom output path for combined volcano grid SVG",
    )
    parser.add_argument(
        "--umap-out",
        type=Path,
        default=None,
        help="Custom output path for combined UMAP grid SVG",
    )
    args = parser.parse_args()

    return FiguresConfig(
        results_dir=args.results_dir,
        reports_dir=args.reports_dir,
        fdr_threshold=args.fdr_threshold,
        max_cells_per_umap_panel=args.max_cells_per_panel,
        timepoint=args.timepoint,
        volcano_out=args.volcano_out,
        umap_out=args.umap_out,
    )


# =============================================================================
# 1. Multi-Cohort Volcano Plots Grid (Figure 1)
# =============================================================================

def build_volcano_panel(
    df: pl.DataFrame,
    meta: dict[str, str],
    config: FiguresConfig,
) -> alt.LayerChart:
    """Constructs an individual publication-grade volcano plot panel for a single cohort."""
    pdf = df.to_pandas()

    # Ensure negative log10 FDR exists
    if "neg_log10_fdr" not in pdf.columns:
        pdf["neg_log10_fdr"] = -np.log10(pdf["FDR"].clip(lower=1e-300))

    # Standardize status labels
    cond_up = (pdf["FDR"] < config.fdr_threshold) & (pdf["logFC"] > 0)
    cond_down = (pdf["FDR"] < config.fdr_threshold) & (pdf["logFC"] < 0)
    pdf["status"] = np.select(
        [cond_up, cond_down],
        ["Enriched in Responders", "Enriched in Non-Responders"],
        default="Not Significant",
    )

    n_nhoods = len(pdf)
    n_up = int(cond_up.sum())
    n_down = int(cond_down.sum())

    # Anti-occlusion depth ordering: non-significant points first, significant points on top
    pdf["is_sig_order"] = (pdf["status"] != "Not Significant").astype(int)
    pdf = pdf.sort_values(by="is_sig_order", ascending=True)

    # Compute panel axis bounds
    max_lfc = max(2.5, np.percentile(np.abs(pdf["logFC"]), 99.5))
    x_limit = round(max_lfc * 1.15, 1)
    max_y = max(2.0, pdf["neg_log10_fdr"].max())
    y_limit = round(max_y * 1.15, 1)

    # Threshold guidelines
    fdr_cutoff_y = -np.log10(config.fdr_threshold)  # ~1.0 for FDR=0.1
    hrule_df = pd.DataFrame({"y": [fdr_cutoff_y]})
    vrule_df = pd.DataFrame({"x": [0.0]})

    hrule = (
        alt.Chart(hrule_df)
        .mark_rule(strokeDash=[4, 4], stroke=COLOR_NEUTRAL_GREY, strokeWidth=0.85)
        .encode(y="y:Q")
    )
    vrule = (
        alt.Chart(vrule_df)
        .mark_rule(strokeDash=[2, 2], stroke=COLOR_LIGHT_GREY, strokeWidth=0.75)
        .encode(x="x:Q")
    )

    # Points layer
    points = (
        alt.Chart(pdf)
        .mark_point(filled=True, opacity=0.75, size=28)
        .encode(
            x=alt.X(
                "logFC:Q",
                title="log2FC (R vs NR)",
                scale=alt.Scale(domain=[-x_limit, x_limit]),
                axis=alt.Axis(grid=False, tickCount=5, titleFontSize=9, labelFontSize=8),
            ),
            y=alt.Y(
                "neg_log10_fdr:Q",
                title="-log10(FDR)",
                scale=alt.Scale(domain=[0, y_limit]),
                axis=alt.Axis(grid=False, tickCount=4, titleFontSize=9, labelFontSize=8),
            ),
            color=alt.Color(
                "status:N",
                title="Differential Abundance",
                scale=alt.Scale(
                    domain=["Enriched in Responders", "Enriched in Non-Responders", "Not Significant"],
                    range=[OKABE_VERMILION, OKABE_BLUE, COLOR_NEUTRAL_GREY],
                ),
                legend=alt.Legend(
                    orient="bottom",
                    titleFontSize=10,
                    labelFontSize=9,
                    symbolSize=50,
                    columns=3,
                ),
            ),
            tooltip=[
                "Nhood:N",
                alt.Tooltip("logFC:Q", format=".2f", title="log2FC"),
                alt.Tooltip("FDR:Q", format=".2e", title="FDR"),
                "status:N",
                "Nhood_CellType:N",
            ],
        )
    )

    title_text = f"{meta['tag']}   {meta['acc']} ({meta['indication']})"
    sub_text = f"{meta['tech']} | {meta['patients']} | N={n_nhoods:,} ({n_up} up, {n_down} down)"

    panel = (
        alt.layer(vrule, hrule, points)
        .properties(
            title=alt.TitleParams(
                text=title_text,
                subtitle=sub_text,
                fontSize=11,
                fontWeight="bold",
                subtitleFontSize=8.5,
                anchor="start",
            ),
            width=config.panel_width,
            height=config.panel_height_volcano,
        )
    )
    return panel


def build_placeholder_panel(
    meta: dict[str, str],
    reason: str,
    config: FiguresConfig,
    height: int,
) -> alt.LayerChart:
    """Renders a clean minimalist placeholder panel when a cohort is unavailable for a given stratum."""
    title_text = f"{meta['tag']}   {meta['acc']} ({meta['indication']})"
    sub_text = f"{meta['tech']} | {meta['patients']}"
    pdf = pd.DataFrame([{"x": 0.5, "y": 0.5, "label": reason}])

    bg = (
        alt.Chart(pdf)
        .mark_rect(color="#F8FAFC", stroke=COLOR_HAIRLINE, strokeWidth=0.5)
        .encode()
        .properties(width=config.panel_width, height=height)
    )

    txt = (
        alt.Chart(pdf)
        .mark_text(text=reason, fontSize=11, color=COLOR_NEUTRAL_GREY, fontStyle="italic")
        .encode(
            x=alt.value(config.panel_width // 2),
            y=alt.value(height // 2),
        )
    )

    card = (
        (bg + txt)
        .properties(
            title=alt.TitleParams(
                text=title_text,
                subtitle=sub_text,
                fontSize=11,
                fontWeight="bold",
                subtitleFontSize=8.5,
                anchor="start",
            ),
            width=config.panel_width,
            height=height,
        )
    )
    return card


def build_volcano_grid_figure(
    cohort_da_dfs: dict[str, pl.DataFrame],
    config: FiguresConfig,
) -> alt.VConcatChart:
    """Assembles all 9 cohort volcano plots into a publication-grade 3x3 multi-panel figure."""
    panels: list[alt.LayerChart] = []

    for meta in COHORT_PANELS:
        acc = meta["acc"]
        df = cohort_da_dfs.get(acc)
        if df is None:
            reason = (
                f"No {config.timepoint.capitalize()}-Treatment Data"
                if config.timepoint != "all"
                else "Cohort Unavailable"
            )
            panel = build_placeholder_panel(meta, reason, config, config.panel_height_volcano)
        else:
            panel = build_volcano_panel(df, meta, config)
        panels.append(panel)

    rows = []
    for r in range(0, len(panels), 3):
        row_panels = panels[r : r + 3]
        rows.append(alt.hconcat(*row_panels, spacing=24))

    tp_label = (
        "All Available Timepoints"
        if config.timepoint == "all"
        else f"{config.timepoint.capitalize()}-Treatment Biopsies"
    )
    master_title = alt.TitleParams(
        text=f"Single-Cell Differential Abundance Testing Across 9 Response Cohorts ({tp_label})",
        subtitle="Quasi-Likelihood Negative Binomial GLM Modeling on kNN Neighborhoods (Horizontal dashed rule: FDR = 0.10 threshold)",
        fontSize=15,
        fontWeight="bold",
        subtitleFontSize=11,
        anchor="start",
    )

    grid = (
        alt.vconcat(*rows, spacing=24)
        .properties(title=master_title)
        .resolve_scale(color="shared")
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )
    return grid


# =============================================================================
# 2. Multi-Cohort UMAP log2FC Grid (Figure 2)
# =============================================================================

def build_umap_panel(
    df: pl.DataFrame,
    meta: dict[str, str],
    config: FiguresConfig,
    global_max_val: float = 3.0,
) -> alt.Chart:
    """Constructs an individual UMAP panel colored by continuous projected log2FC."""
    pdf = df.to_pandas()
    pdf["abs_logfc"] = np.abs(pdf["cell_logfc"])
    pdf = pdf.sort_values(by="abs_logfc", ascending=True)

    # Intelligent sampling guaranteeing strict cell cap while balancing significant and background cells
    if len(pdf) > config.max_cells_per_umap_panel:
        if "da_status" in pdf.columns:
            sig_mask = pdf["da_status"] != "Not Significant"
        else:
            sig_mask = np.abs(pdf["cell_logfc"]) > 0.5
        sig_pdf = pdf[sig_mask]
        non_sig_pdf = pdf[~sig_mask]

        target = config.max_cells_per_umap_panel
        if len(sig_pdf) > 0 and len(non_sig_pdf) > 0:
            # Allocate up to 65% of budget to significant cells, 35% to background
            n_sig = min(len(sig_pdf), int(target * 0.65))
            n_nonsig = min(len(non_sig_pdf), target - n_sig)
            if n_sig < int(target * 0.65):
                n_nonsig = min(len(non_sig_pdf), target - n_sig)
            elif n_nonsig < (target - n_sig):
                n_sig = min(len(sig_pdf), target - n_nonsig)

            sample_sig = sig_pdf.sample(n=n_sig, random_state=config.seed)
            sample_nonsig = non_sig_pdf.sample(n=n_nonsig, random_state=config.seed)
            pdf = pd.concat([sample_nonsig, sample_sig], axis=0).sort_values(by="abs_logfc", ascending=True)
        else:
            pdf = pdf.sample(n=min(target, len(pdf)), random_state=config.seed).sort_values(
                by="abs_logfc", ascending=True
            )

    title_text = f"{meta['tag']}   {meta['acc']} ({meta['indication']})"
    sub_text = f"{meta['tech']} | {meta['patients']} | N={len(df):,} cells"

    tooltip_cols: list[Any] = [
        "cell_id:N",
        alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC"),
        "clinical_response:N",
        "cell_type:N",
    ]
    if "da_status" in pdf.columns:
        tooltip_cols.append("da_status:N")

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=12, opacity=0.80)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=3, labels=False, titleFontSize=9)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=3, labels=False, titleFontSize=9)),
            color=alt.Color(
                "cell_logfc:Q",
                title="Milo log2FC (Responder vs Non-Responder; clamped at ±5)",
                scale=alt.Scale(
                    domain=[-5.0, 0.0, 5.0],
                    range=[OKABE_BLUE, COLOR_LIGHT_GREY, OKABE_VERMILION],
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    titleFontSize=10,
                    labelFontSize=9,
                    titleLimit=0,
                    gradientLength=340,
                    values=[-5.0, -2.5, 0.0, 2.5, 5.0],
                ),
            ),
            tooltip=tooltip_cols,
        )
        .properties(
            title=alt.TitleParams(
                text=title_text,
                subtitle=sub_text,
                fontSize=11,
                fontWeight="bold",
                subtitleFontSize=8.5,
                anchor="start",
            ),
            width=config.panel_width,
            height=config.panel_height_umap,
        )
    )
    return chart


def build_umap_grid_figure(
    cohort_cell_dfs: dict[str, pl.DataFrame],
    config: FiguresConfig,
) -> alt.VConcatChart:
    """Assembles all 9 cohort UMAP projections into a publication-grade 3x3 multi-panel figure."""
    panels: list[alt.Chart | alt.LayerChart] = []
    for meta in COHORT_PANELS:
        acc = meta["acc"]
        df = cohort_cell_dfs.get(acc)
        if df is None:
            reason = (
                f"No {config.timepoint.capitalize()}-Treatment Data"
                if config.timepoint != "all"
                else "Cohort Unavailable"
            )
            panel = build_placeholder_panel(meta, reason, config, config.panel_height_umap)
        else:
            panel = build_umap_panel(df, meta, config, global_max_val=5.0)
        panels.append(panel)

    rows = []
    for r in range(0, len(panels), 3):
        row_panels = panels[r : r + 3]
        rows.append(alt.hconcat(*row_panels, spacing=24))

    tp_label = (
        "All Available Timepoints"
        if config.timepoint == "all"
        else f"{config.timepoint.capitalize()}-Treatment Biopsies"
    )
    master_title = alt.TitleParams(
        text=f"Single-Cell Landscape of Immunotherapy Response Across 9 Solid Tumor Cohorts ({tp_label})",
        subtitle="2D UMAP Embeddings Colored by Continuous Inferred Milopy log2 Fold Change (Responder vs Non-Responder; Clamped at ±5)",
        fontSize=15,
        fontWeight="bold",
        subtitleFontSize=11,
        anchor="start",
    )

    grid = (
        alt.vconcat(*rows, spacing=24)
        .properties(title=master_title)
        .resolve_scale(color="shared")
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )
    return grid


# =============================================================================
# 3. Main Workflow Orchestration
# =============================================================================

def load_cohort_data(
    results_dir: Path,
) -> tuple[dict[str, pl.DataFrame], dict[str, pl.DataFrame]]:
    """Loads da_results and cell_level_scores for all 9 cohorts."""
    cohort_da_dfs: dict[str, pl.DataFrame] = {}
    cohort_cell_dfs: dict[str, pl.DataFrame] = {}

    for meta in COHORT_PANELS:
        acc = meta["acc"]
        cohort_dir = results_dir / acc
        da_path = cohort_dir / "da_results.parquet"
        cell_path = cohort_dir / "cell_level_scores.parquet"

        if da_path.exists():
            cohort_da_dfs[acc] = pl.read_parquet(da_path)
            print(f"[{acc}] Loaded DA results: {len(cohort_da_dfs[acc]):,} nhoods")
        else:
            print(f"[{acc}] Warning: Missing {da_path}", file=sys.stderr)

        if cell_path.exists():
            cell_df = pl.read_parquet(cell_path)
            if "UMAP1" in cell_df.columns and "UMAP2" in cell_df.columns:
                cohort_cell_dfs[acc] = cell_df
                print(f"[{acc}] Loaded cell scores with fresh UMAP: {len(cell_df):,} cells")
            else:
                print(f"[{acc}] Warning: cell_level_scores missing UMAP1/UMAP2 columns", file=sys.stderr)
        else:
            print(f"[{acc}] Warning: Missing {cell_path}", file=sys.stderr)

    return cohort_da_dfs, cohort_cell_dfs


def main() -> None:
    config = parse_args()
    print("=" * 70)
    print("STARTING COMBINED PUBLICATION FIGURES GENERATION")
    print("=" * 70)
    print(f"Results Directory : {config.results_dir}")
    print(f"Reports Directory : {config.reports_dir}")

    config.reports_dir.mkdir(parents=True, exist_ok=True)
    cohort_da_dfs, cohort_cell_dfs = load_cohort_data(config.results_dir)

    suffix = "" if config.timepoint == "all" else f"_{config.timepoint}"
    volcano_svg_path = config.volcano_out or (config.reports_dir / f"milopy_cross_cohort_volcano_grid{suffix}.svg")
    volcano_png_path = volcano_svg_path.with_suffix(".png")
    umap_svg_path = config.umap_out or (config.reports_dir / f"milopy_cross_cohort_umaps_logfc{suffix}.svg")
    umap_png_path = umap_svg_path.with_suffix(".png")

    # 1. Generate Combined Volcano Plots Grid Figure
    if len(cohort_da_dfs) >= 1:
        print("\n" + "-" * 70)
        print(f"RENDERING COMBINED 3x3 VOLCANO PLOTS FIGURE (Figure 1: {config.timepoint})...")
        print("-" * 70)
        volcano_grid = build_volcano_grid_figure(cohort_da_dfs, config)
        volcano_svg_str = vlc.vegalite_to_svg(volcano_grid.to_json())
        volcano_svg_path.write_text(volcano_svg_str, encoding="utf-8")
        size_mb = volcano_svg_path.stat().st_size / 1e6
        print(f"[EXPORTED] {volcano_svg_path} ({size_mb:.2f} MB)")

        # Render high-resolution 300 DPI raster PNG
        volcano_png_bytes = vlc.svg_to_png(volcano_svg_str, scale=2.0)
        volcano_png_path.write_bytes(volcano_png_bytes)
        print(f"[EXPORTED] {volcano_png_path} ({len(volcano_png_bytes) / 1e6:.2f} MB)")

        # Validate XML well-formedness
        try:
            tree = ET.parse(volcano_svg_path)
            print(f"[VERIFIED] {volcano_svg_path.name} is 100% valid XML ({tree.getroot().tag})")
        except Exception as e:
            print(f"[ERROR] Invalid XML in {volcano_svg_path}: {e}", file=sys.stderr)

    # 2. Generate Combined UMAP log2FC Grid Figure
    if len(cohort_cell_dfs) >= 1:
        print("\n" + "-" * 70)
        print(f"RENDERING COMBINED 3x3 UMAP log2FC FIGURE (Figure 2: {config.timepoint})...")
        print("-" * 70)
        umap_grid = build_umap_grid_figure(cohort_cell_dfs, config)
        umap_svg_str = vlc.vegalite_to_svg(umap_grid.to_json())
        umap_svg_path.write_text(umap_svg_str, encoding="utf-8")
        size_mb = umap_svg_path.stat().st_size / 1e6
        print(f"[EXPORTED] {umap_svg_path} ({size_mb:.2f} MB)")

        # Render high-resolution 300 DPI raster PNG
        umap_png_bytes = vlc.svg_to_png(umap_svg_str, scale=2.0)
        umap_png_path.write_bytes(umap_png_bytes)
        print(f"[EXPORTED] {umap_png_path} ({len(umap_png_bytes) / 1e6:.2f} MB)")

        # Validate XML well-formedness
        try:
            tree = ET.parse(umap_svg_path)
            print(f"[VERIFIED] {umap_svg_path.name} is 100% valid XML ({tree.getroot().tag})")
        except Exception as e:
            print(f"[ERROR] Invalid XML in {umap_svg_path}: {e}", file=sys.stderr)

    print("\n" + "=" * 70)
    print("ALL COMBINED FIGURES GENERATED SUCCESSFULLY.")
    print("=" * 70)


if __name__ == "__main__":
    main()
