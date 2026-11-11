#!/usr/bin/env python3
"""Step 6: Publication Figure Generation for Multi-Resolution Reference Benchmarking.

Generates:
  - Figure 1: Comprehensive 18-panel predictive capacity benchmark across references,
              resolutions, cluster regimes, dataset count (up to n=9+), continuous single-cell scaling, and titrations.
  - Figure 2: Signature collinearity benchmark (condition number vs resolution and cluster count).
  - Figure 3: Classifier selection and optimization decision framework (strategy grid, Pareto frontier,
              cohort generalizability matrix, diversity scaling law).
  - Figure 4: Winning classifier biomarker drivers (cell state odds ratios and Milo DA concordance).
  - Figure 5: Cohort-averaged predictive capacity and marginal gain dynamics.

Enforces Nature Methods minimalist wireframe standard:
  - Fixed 180.00 mm width, unconstrained height
  - Typography scaled +30% with >= 13.0 pt floor guarantee
  - Subpanel letter badges (a) through (r) in top-left
  - Upright horizontal axis numbers (labelAngle=0)
  - Zero baseline gutter artifacts
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any, Final

import altair as alt
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats

# Local package imports
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "plotting_utils" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "plotting_utils" / "src"))
if str(REPO_ROOT / "packages" / "tme_response" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "tme_response" / "src"))

from plotting_utils import rasterize_svg_to_png
from tme_response.visualization.display_labels import with_display_labels

PANEL_W_FIG1: Final[int] = 280
PANEL_H_FIG1: Final[int] = 195
PANEL_W_2COL: Final[int] = 320
PANEL_H_2COL: Final[int] = 230


class PlotConfig(BaseModel):
    """Immutable configuration for benchmark figure plotting."""

    model_config = ConfigDict(frozen=True)

    data_dir: Path
    results_dir: Path


def format_corr_stat(val: float, p: float) -> str:
    """Format correlation coefficient and p-value with significance asterisks."""
    if p < 0.001:
        star = "***"
        p_str = f"{p:.1e}"
    elif p < 0.01:
        star = "**"
        p_str = f"{p:.3f}"
    elif p < 0.05:
        star = "*"
        p_str = f"{p:.3f}"
    else:
        star = "ns"
        p_str = f"{p:.3f}"
    return f"{val:+.2f}{star} (p = {p_str})"


def get_panel_corr_lines(
    x: np.ndarray[Any, Any],
    y: np.ndarray[Any, Any],
    x_log: bool = False,
    y_log: bool = False,
) -> tuple[str, str]:
    """Compute Pearson and Spearman correlations and return formatted 2-line subtitle."""
    x_arr = np.array(x, dtype=float)
    y_arr = np.array(y, dtype=float)
    x_lin = np.log10(x_arr) if x_log else x_arr
    y_lin = np.log10(y_arr) if y_log else y_arr
    valid = ~np.isnan(x_lin) & ~np.isnan(y_lin) & ~np.isinf(x_lin) & ~np.isinf(y_lin)
    if np.sum(valid) < 3:
        return ("Pearson r = N/A", "Spearman ρ = N/A")
    r_res = stats.pearsonr(x_lin[valid], y_lin[valid])
    rho_res = stats.spearmanr(x_arr[valid], y_arr[valid])
    return (
        f"Pearson r = {format_corr_stat(float(r_res.statistic), float(r_res.pvalue))}",
        f"Spearman ρ = {format_corr_stat(float(rho_res.statistic), float(rho_res.pvalue))}",
    )


def post_process_svg(
    svg_path: Path,
    n_panels: int,
    width_mm: float = 180.0,
    min_font_pt: float = 13.0,
) -> None:
    """Apply strict Nature Methods vector formatting to generated SVG files."""
    svg_content = svg_path.read_text(encoding="utf-8")
    letters = [f"({chr(ord('a') + i)})" for i in range(n_panels)]
    letter_idx = 0

    def inject_panel_letter(m: re.Match[str]) -> str:
        nonlocal letter_idx
        if letter_idx >= len(letters):
            return m.group(0)
        letter = letters[letter_idx]
        letter_idx += 1
        content = m.group(0)
        prefix = '<g class="mark-group role-title">'
        label_xml = (
            '<g class="mark-text role-panel-letter" pointer-events="none">'
            f'<text text-anchor="start" transform="translate(-60,-38)" '
            'font-family="sans-serif" font-size="16px" font-weight="bold" '
            f'fill="#000" opacity="1">{letter}</text></g>'
        )
        return prefix + label_xml + content[len(prefix):]

    svg_content = re.sub(
        r'<g class="mark-group role-title">.*?</g></g></g>',
        inject_panel_letter,
        svg_content,
    )

    # Enforce minimum font size
    svg_content = re.sub(
        r'font-size="([0-9\.]+)px"',
        lambda m: f'font-size="{max(min_font_pt, float(m.group(1))):.1f}px"',
        svg_content,
    )
    svg_content = re.sub(
        r'font-size:\s*([0-9\.]+)px',
        lambda m: f'font-size: {max(min_font_pt, float(m.group(1))):.1f}px',
        svg_content,
    )
    svg_content = re.sub(
        r'(<text[^>]*)(>Deconvolution Inference \(Our\)</text>)',
        lambda m: (
            m.group(1) + ' font-weight="bold"'
            if 'font-weight=' not in m.group(1)
            else re.sub(r'font-weight="[^"]*"', 'font-weight="bold"', m.group(1))
        ) + m.group(2),
        svg_content,
    )

    # Set tight cropped viewBox and exact physical width & height so zero outer whitespace remains
    match_dims = re.search(
        r'<svg([^>]+)width="([0-9\.]+)"\s+height="([0-9\.]+)"([^>]*)viewBox="0\s+0\s+([0-9\.]+)\s+([0-9\.]+)"',
        svg_content,
    )
    if match_dims:
        w = float(match_dims.group(2))
        h = float(match_dims.group(3))
        svg_path.write_text(svg_content, encoding="utf-8")
        vx, vy, vw, vh = 0.0, 0.0, w, h
        try:
            import io
            import subprocess
            from PIL import Image

            png_bytes = subprocess.check_output(
                ["rsvg-convert", "-w", str(int(round(w))), "-h", str(int(round(h))), "-b", "white", str(svg_path)]
            )
            im = Image.open(io.BytesIO(png_bytes)).convert("RGB")
            arr = np.asarray(im)
            non_white = np.any(arr < 252, axis=2)
            if np.any(non_white):
                ys, xs = np.where(non_white)
                pad = 1.0
                vx = max(0.0, float(xs.min()) - pad)
                vy = max(0.0, float(ys.min()) - pad)
                vw = min(w - vx, float(xs.max() - xs.min() + 1) + 2.0 * pad)
                vh = min(h - vy, float(ys.max() - ys.min() + 1) + 2.0 * pad)
        except Exception:  # noqa: BLE001
            pass
        height_mm = width_mm * (vh / max(vw, 1.0))
        prefix = (
            f'<svg{match_dims.group(1)}width="{width_mm:.2f}mm" height="{height_mm:.2f}mm"'
            f'{match_dims.group(4)}viewBox="{vx:.1f} {vy:.1f} {vw:.1f} {vh:.1f}"'
        )
        svg_content = svg_content[:match_dims.start()] + prefix + svg_content[match_dims.end():]
        svg_content = re.sub(
            r'<rect width="[0-9\.]+" height="[0-9\.]+" fill="#FFFFFF"/>',
            f'<rect x="{vx:.1f}" y="{vy:.1f}" width="{vw:.1f}" height="{vh:.1f}" fill="#FFFFFF"/>',
            svg_content,
            count=1,
        )

    svg_path.write_text(svg_content, encoding="utf-8")


def assign_cluster_regime(e: pl.Expr) -> pl.Expr:
    """Bin cluster counts into standardized granularity regimes."""
    return (
        pl.when(e < 15).then(pl.lit("< 15"))
        .when((e >= 15) & (e <= 19)).then(pl.lit("15–19"))
        .when((e >= 20) & (e <= 24)).then(pl.lit("20–24"))
        .when((e >= 25) & (e <= 29)).then(pl.lit("25–29"))
        .when((e >= 30) & (e <= 34)).then(pl.lit("30–34"))
        .when((e >= 35) & (e <= 39)).then(pl.lit("35–39"))
        .otherwise(pl.lit("≥ 40"))
    )


def plot_fig1_predictive_capacity(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 1: Comprehensive 18-panel predictive capacity benchmark."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        fallback_dirs = (
            Path("output/sc9_milo_deconv_validation/05_concordance"),
            Path("output/manuscript_pipeline/3f8_compare_concordance_multires"),
            Path("output/sade_feldman_deconv_validation"),
        )
        for fb in fallback_dirs:
            if (fb / "multi_resolution_benchmark_summary.parquet").exists():
                bench_file = fb / "multi_resolution_benchmark_summary.parquet"
                break
        if not bench_file.exists():
            return Failure(f"multi_resolution_benchmark_summary.parquet not found in {data_dir}")

    try:
        df = pl.read_parquet(bench_file)

        if "subsample_fraction" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("subsample_fraction"))

        df_full = df.filter(pl.col("reference_category") != "Cell-Subsample")

        bin_order = ["< 15", "15–19", "20–24", "25–29", "30–34", "35–39", "≥ 40"]
        df_full = df_full.with_columns(assign_cluster_regime(pl.col("n_clusters")).alias("cluster_bin"))
        df_full_pd = df_full.to_pandas()
        df_all_pd = df.to_pandas()

        np.random.seed(42)
        df_full_pd["res_jitter"] = np.random.uniform(-8.0, 8.0, len(df_full_pd))
        df_full_pd["cluster_jitter"] = np.random.uniform(-8.0, 8.0, len(df_full_pd))
        df_full_pd["dataset_jitter"] = np.random.uniform(-8.0, 8.0, len(df_full_pd))

        # Resolution summary
        df_res_summary = (
            df_full.group_by("resolution")
            .agg(
                pl.col("melanoma_multivariate_auc").mean().alias("mean_mel_auc"),
                pl.col("pancancer_multivariate_auc").mean().alias("mean_pan_auc"),
                pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"),
            )
            .sort("resolution")
            .to_pandas()
        )

        # Cluster regimes summary
        df_cluster_summary = (
            df_full.group_by("cluster_bin")
            .agg(
                pl.col("melanoma_multivariate_auc").mean().alias("mean_mel_auc"),
                pl.col("pancancer_multivariate_auc").mean().alias("mean_pan_auc"),
                pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"),
            )
            .to_pandas()
        )
        df_cluster_summary["cluster_bin"] = pd.Categorical(
            df_cluster_summary["cluster_bin"], categories=bin_order, ordered=True
        )
        df_cluster_summary = df_cluster_summary.sort_values("cluster_bin").reset_index(drop=True)

        # Dataset count summary
        df_ds_summary = (
            df_full.group_by("n_datasets")
            .agg(
                pl.col("melanoma_multivariate_auc").mean().alias("mean_mel_auc"),
                pl.col("pancancer_multivariate_auc").mean().alias("mean_pan_auc"),
                pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"),
            )
            .sort("n_datasets")
            .to_pandas()
        )

        # Smooth log-saturation fit curves
        valid_cells = df.filter(pl.col("n_cells") > 0)
        x_cells_arr = valid_cells["n_cells"].to_numpy().astype(np.float64)
        x_log = np.log(x_cells_arr)

        m_fit = np.polyfit(x_log, valid_cells["melanoma_multivariate_auc"].to_numpy(), 1)
        p_fit = np.polyfit(x_log, valid_cells["pancancer_multivariate_auc"].to_numpy(), 1)
        c_fit = np.polyfit(x_log, valid_cells["mean_cohort_multivariate_auc"].to_numpy(), 1)

        x_grid = np.geomspace(max(1000.0, float(np.min(x_cells_arr))), float(np.max(x_cells_arr)), 100)
        df_reg_mel = pd.DataFrame({"n_cells": x_grid, "melanoma_multivariate_auc": m_fit[0] * np.log(x_grid) + m_fit[1]})
        df_reg_pan = pd.DataFrame({"n_cells": x_grid, "pancancer_multivariate_auc": p_fit[0] * np.log(x_grid) + p_fit[1]})
        df_reg_coh = pd.DataFrame({"n_cells": x_grid, "mean_cohort_multivariate_auc": c_fit[0] * np.log(x_grid) + c_fit[1]})

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined", "Cell-Subsample"]
        ref_order = sorted(df_full["reference_type"].unique().to_list())
        color_scale = alt.Scale(scheme="tableau20", domain=ref_order)

        ref_legend = alt.Legend(
            orient="bottom",
            direction="horizontal",
            columns=5,
            titleAnchor="middle",
            labelLimit=250,
            titleLimit=400,
            title="Reference Dataset / Strategy",
            labelFontSize=13,
            titleFontSize=14,
            titleFontWeight="bold",
        )

        # Precompute df_delta for panel r
        df_delta = (
            df_full.sort("resolution")
            .with_columns(
                pl.col("melanoma_multivariate_auc")
                .shift(1)
                .over("reference_type")
                .alias("prev_auc")
            )
            .with_columns(
                (pl.col("melanoma_multivariate_auc") - pl.col("prev_auc")).alias("delta_auc")
            )
            .filter(pl.col("delta_auc").is_not_null())
            .group_by("resolution")
            .agg(pl.col("delta_auc").mean().alias("mean_delta_auc"))
            .sort("resolution")
            .to_pandas()
        )

        # Compute correlation statistics for Figure 1 panels
        r1_a_p, r1_a_s = get_panel_corr_lines(df_full["resolution"].to_numpy(), df_full["melanoma_multivariate_auc"].to_numpy())
        r1_b_p, r1_b_s = get_panel_corr_lines(df_full["resolution"].to_numpy(), df_full["pancancer_multivariate_auc"].to_numpy())
        r1_c_p, r1_c_s = get_panel_corr_lines(df_full["resolution"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())

        r2_d_p, r2_d_s = get_panel_corr_lines(df_full["n_clusters"].to_numpy(), df_full["melanoma_multivariate_auc"].to_numpy())
        r2_e_p, r2_e_s = get_panel_corr_lines(df_full["n_clusters"].to_numpy(), df_full["pancancer_multivariate_auc"].to_numpy())
        r2_f_p, r2_f_s = get_panel_corr_lines(df_full["n_clusters"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())

        r3_g_p, r3_g_s = get_panel_corr_lines(df_full["n_datasets"].to_numpy(), df_full["melanoma_multivariate_auc"].to_numpy())
        r3_h_p, r3_h_s = get_panel_corr_lines(df_full["n_datasets"].to_numpy(), df_full["pancancer_multivariate_auc"].to_numpy())
        r3_i_p, r3_i_s = get_panel_corr_lines(df_full["n_datasets"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())

        vc = df.filter(pl.col("n_cells") > 0)
        r4_j_p, r4_j_s = get_panel_corr_lines(vc["n_cells"].to_numpy(), vc["melanoma_multivariate_auc"].to_numpy(), x_log=True)
        r4_k_p, r4_k_s = get_panel_corr_lines(vc["n_cells"].to_numpy(), vc["pancancer_multivariate_auc"].to_numpy(), x_log=True)
        r4_l_p, r4_l_s = get_panel_corr_lines(vc["n_cells"].to_numpy(), vc["mean_cohort_multivariate_auc"].to_numpy(), x_log=True)

        df_sub_raw = df.filter(
            (pl.col("reference_category") == "Cell-Subsample")
            | (pl.col("reference_type") == "All-9-Milo-Combined")
        )
        r5_m_p, r5_m_s = get_panel_corr_lines(df_sub_raw["subsample_fraction"].to_numpy(), df_sub_raw["melanoma_multivariate_auc"].to_numpy())
        r5_n_p, r5_n_s = get_panel_corr_lines(df_sub_raw["subsample_fraction"].to_numpy(), df_sub_raw["pancancer_multivariate_auc"].to_numpy())
        r5_o_p, r5_o_s = get_panel_corr_lines(df_sub_raw["subsample_fraction"].to_numpy(), df_sub_raw["mean_cohort_multivariate_auc"].to_numpy())

        strat_map = {"Single Dataset": 1, "Criteria-Combined": 2, "Random-Combined": 3}
        strat_ranks = np.array([strat_map.get(c, 1) for c in df_full["reference_category"].to_list()])
        r6_p_p, r6_p_s = get_panel_corr_lines(strat_ranks, df_full["melanoma_multivariate_auc"].to_numpy())
        r6_q_p, r6_q_s = get_panel_corr_lines(strat_ranks, df_full["pancancer_multivariate_auc"].to_numpy())
        r6_r_p, r6_r_s = get_panel_corr_lines(df_delta["resolution"].to_numpy(), df_delta["mean_delta_auc"].to_numpy())

        base = alt.Chart(df_full_pd)
        base_all = alt.Chart(df_all_pd)

        # Row 1: AUC by Leiden Resolution (Boxplot + Jitter + Mean Trendline)
        panel_a = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    xOffset="res_jitter:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=ref_legend),
                    tooltip=["reference_type", "reference_category", "resolution", "melanoma_multivariate_auc", "n_clusters"],
                ),
                alt.Chart(df_res_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("resolution:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_mel_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
            )
            .properties(title=["Melanoma Response AUC", r1_a_p, r1_a_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_b = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer Response AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    xOffset="res_jitter:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "reference_category", "resolution", "pancancer_multivariate_auc", "n_clusters"],
                ),
                alt.Chart(df_res_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("resolution:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_pan_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
            )
            .properties(title=["Pan-Cancer Response AUC", r1_b_p, r1_b_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_c = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    xOffset="res_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "reference_category", "resolution", "mean_cohort_multivariate_auc", "n_clusters"],
                ),
                alt.Chart(df_res_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("resolution:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
            )
            .properties(title=["Mean Cohort Response AUC", r1_c_p, r1_c_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Row 2: AUC vs Number of Cell States (Boxplot + Jitter + Mean Trendline)
        panel_d = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    xOffset="cluster_jitter:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_clusters", "resolution", "melanoma_multivariate_auc"],
                ),
                alt.Chart(df_cluster_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("cluster_bin:O", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_mel_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
            )
            .properties(title=["Melanoma AUC vs Cell States", r2_d_p, r2_d_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_e = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer Response AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    xOffset="cluster_jitter:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_clusters", "resolution", "pancancer_multivariate_auc"],
                ),
                alt.Chart(df_cluster_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("cluster_bin:O", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_pan_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
            )
            .properties(title=["Pan-Cancer AUC vs Cell States", r2_e_p, r2_e_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_f = (
            alt.layer(
                base.mark_boxplot(size=20, outliers=False, opacity=0.8).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    xOffset="cluster_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_clusters", "resolution", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_cluster_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("cluster_bin:O", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
            )
            .properties(title=["Mean Cohort AUC vs Cell States", r2_f_p, r2_f_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Row 3: AUC vs Number of Integrated Datasets
        panel_g = (
            alt.layer(
                base.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    xOffset="dataset_jitter:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "melanoma_multivariate_auc"],
                ),
                alt.Chart(df_ds_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("n_datasets:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_mel_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
            )
            .properties(title=["Melanoma AUC vs Datasets", r3_g_p, r3_g_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_h = (
            alt.layer(
                base.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer Response AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    xOffset="dataset_jitter:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "pancancer_multivariate_auc"],
                ),
                alt.Chart(df_ds_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("n_datasets:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_pan_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
            )
            .properties(title=["Pan-Cancer AUC vs Datasets", r3_h_p, r3_h_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_i = (
            alt.layer(
                base.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.6).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    xOffset="dataset_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_ds_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("n_datasets:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
            )
            .properties(title=["Mean Cohort AUC vs Datasets", r3_i_p, r3_i_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Row 4: Continuous Single-Cell Scaling (Log-Fit)
        panel_j = (
            alt.layer(
                base_all.mark_circle(size=26, opacity=0.45).encode(
                    x=alt.X("n_cells:Q", title="Total Single Cells", scale=alt.Scale(type="log")),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "melanoma_multivariate_auc"],
                ),
                alt.Chart(df_reg_mel).mark_line(color="#1f77b4", strokeWidth=2.5).encode(
                    x="n_cells:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
            )
            .properties(title=["Melanoma AUC vs Cells (log)", r4_j_p, r4_j_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_k = (
            alt.layer(
                base_all.mark_circle(size=26, opacity=0.45).encode(
                    x=alt.X("n_cells:Q", title="Total Single Cells", scale=alt.Scale(type="log")),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer Response AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "pancancer_multivariate_auc"],
                ),
                alt.Chart(df_reg_pan).mark_line(color="#1f77b4", strokeWidth=2.5).encode(
                    x="n_cells:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
            )
            .properties(title=["Pan-Cancer AUC vs Cells (log)", r4_k_p, r4_k_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_l = (
            alt.layer(
                base_all.mark_circle(size=26, opacity=0.45).encode(
                    x=alt.X("n_cells:Q", title="Total Single Cells", scale=alt.Scale(type="log")),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_reg_coh).mark_line(color="#1f77b4", strokeWidth=2.5).encode(
                    x="n_cells:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
            )
            .properties(title=["Mean Cohort AUC vs Cells (log)", r4_l_p, r4_l_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Row 5: Titration Subsampling Curves
        df_sub = (
            df.filter(
                (pl.col("reference_category") == "Cell-Subsample")
                | (pl.col("reference_type") == "All-9-Milo-Combined")
            )
            .sort("subsample_fraction")
        )
        df_sub_mean = (
            df_sub.group_by("subsample_fraction")
            .agg(
                pl.col("melanoma_multivariate_auc").mean().alias("melanoma_multivariate_auc"),
                pl.col("pancancer_multivariate_auc").mean().alias("pancancer_multivariate_auc"),
                pl.col("mean_cohort_multivariate_auc").mean().alias("mean_cohort_multivariate_auc"),
            )
            .sort("subsample_fraction")
            .to_pandas()
        )
        base_sub_pts = alt.Chart(df_sub.to_pandas())
        base_sub_line = alt.Chart(df_sub_mean)

        panel_m = (
            alt.layer(
                base_sub_line.mark_line(strokeWidth=2.5, color="#1f77b4").encode(
                    x=alt.X("subsample_fraction:Q", title="Fraction of Sampled Cells", axis=alt.Axis(format="%")),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                ),
                base_sub_line.mark_point(size=70, filled=True, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    tooltip=["subsample_fraction", "melanoma_multivariate_auc"],
                ),
                base_sub_pts.mark_circle(size=25, opacity=0.35, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    tooltip=["reference_type", "resolution", "melanoma_multivariate_auc"],
                ),
            )
            .properties(title=["Melanoma AUC Titration", r5_m_p, r5_m_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_n = (
            alt.layer(
                base_sub_line.mark_line(strokeWidth=2.5, color="#1f77b4").encode(
                    x=alt.X("subsample_fraction:Q", title="Fraction of Sampled Cells", axis=alt.Axis(format="%")),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer Response AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                ),
                base_sub_line.mark_point(size=70, filled=True, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    tooltip=["subsample_fraction", "pancancer_multivariate_auc"],
                ),
                base_sub_pts.mark_circle(size=25, opacity=0.35, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    tooltip=["reference_type", "resolution", "pancancer_multivariate_auc"],
                ),
            )
            .properties(title=["Pan-Cancer AUC Titration", r5_n_p, r5_n_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_o = (
            alt.layer(
                base_sub_line.mark_line(strokeWidth=2.5, color="#1f77b4").encode(
                    x=alt.X("subsample_fraction:Q", title="Fraction of Sampled Cells", axis=alt.Axis(format="%")),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                base_sub_line.mark_point(size=70, filled=True, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    tooltip=["subsample_fraction", "mean_cohort_multivariate_auc"],
                ),
                base_sub_pts.mark_circle(size=25, opacity=0.35, color="#1f77b4").encode(
                    x="subsample_fraction:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    tooltip=["reference_type", "resolution", "mean_cohort_multivariate_auc"],
                ),
            )
            .properties(title=["Mean Cohort AUC Titration", r5_o_p, r5_o_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Row 6: Strategy Distributions & Marginal Gain
        panel_p = (
            alt.layer(
                base_all.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_category:N", title="Strategy", legend=None),
                ),
                base_all.mark_circle(size=25, opacity=0.5, xOffset=8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "resolution", "melanoma_multivariate_auc"],
                ),
            )
            .properties(title=["Melanoma AUC by Strategy", r6_p_p, r6_p_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_q = (
            alt.layer(
                base_all.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_category:N", title="Strategy", legend=None),
                ),
                base_all.mark_circle(size=25, opacity=0.5, xOffset=8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.50, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "resolution", "pancancer_multivariate_auc"],
                ),
            )
            .properties(title=["Pan-Cancer AUC by Strategy", r6_q_p, r6_q_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        panel_r = (
            alt.Chart(df_delta)
            .mark_bar()
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("mean_delta_auc:Q", title="Mean ΔAUC vs Prev. Res."),
                color=alt.condition(
                    alt.datum.mean_delta_auc > 0,
                    alt.value("#4daf4a"),
                    alt.value("#e41a1c"),
                ),
                tooltip=["resolution", "mean_delta_auc"],
            )
            .properties(title=["Marginal AUC Gain per Step", r6_r_p, r6_r_s], width=PANEL_W_FIG1, height=PANEL_H_FIG1)
        )

        # Assemble figure
        row1 = (panel_a | panel_b | panel_c).resolve_scale(color="shared")
        row2 = (panel_d | panel_e | panel_f).resolve_scale(color="shared")
        row3 = (panel_g | panel_h | panel_i).resolve_scale(color="shared")
        row4 = (panel_j | panel_k | panel_l).resolve_scale(color="shared")

        main_rows = alt.vconcat(row1, row2, row3, row4, spacing=28).resolve_scale(color="shared")
        row5 = (panel_m | panel_n | panel_o).resolve_scale(color="shared")
        row6 = (panel_p | panel_q | panel_r).resolve_scale(color="independent")

        final_chart = (
            alt.vconcat(main_rows, row5, row6, spacing=35)
            .resolve_scale(color="independent")
            .configure_axisY(minExtent=65, maxExtent=65, titlePadding=10)
            .configure_axis(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_legend(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_title(fontSize=14, fontWeight="bold")
        )

        out_path = results_dir / "fig1_predictive_capacity_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_path))

        post_process_svg(out_path, n_panels=18, width_mm=180.0, min_font_pt=13.0)

        png_path = results_dir / "fig1_predictive_capacity_benchmark.png"
        rasterize_svg_to_png(out_path, png_path, scale=3.125, background="white")

        return Success(out_path)

    except Exception as exc:
        return Failure(f"Failed to generate Figure 1: {exc}")


def plot_fig_sc9_summary_panels_cfi(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Summary Figure: Panels c, f, i from Figure 1 re-lettered as a, b, c without headings."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        fallback_dirs = (
            Path("output/sc9_milo_deconv_validation/05_concordance"),
            results_dir,
        )
        for fb in fallback_dirs:
            if (fb / "multi_resolution_benchmark_summary.parquet").exists():
                bench_file = fb / "multi_resolution_benchmark_summary.parquet"
                break
        else:
            return Failure(f"multi_resolution_benchmark_summary.parquet not found in {data_dir}")

    try:
        df = pl.read_parquet(bench_file)

        df_full = df.filter(
            pl.col("reference_category").is_in(["Single Dataset", "Criteria-Combined", "Random-Combined"])
        ).with_columns(
            assign_cluster_regime(pl.col("n_clusters")).alias("cluster_bin"),
        )

        np.random.seed(42)
        n_rows = df_full.height
        df_full = df_full.with_columns(
            pl.Series("res_jitter", np.random.uniform(-0.15, 0.15, n_rows)),
            pl.Series("cluster_jitter", np.random.uniform(-0.15, 0.15, n_rows)),
            pl.Series("dataset_jitter", np.random.uniform(-0.15, 0.15, n_rows)),
        )

        df_full_pd = df_full.to_pandas()

        bin_order = ["< 15", "15–19", "20–24", "25–29", "30–34", "35–39", "≥ 40"]

        df_res_summary = (
            df_full.group_by("resolution")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"))
            .sort("resolution")
            .to_pandas()
        )

        df_cluster_summary = (
            df_full.group_by("cluster_bin")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"))
            .to_pandas()
        )
        df_cluster_summary["bin_rank"] = df_cluster_summary["cluster_bin"].apply(
            lambda b: bin_order.index(b) if b in bin_order else 99
        )
        df_cluster_summary = df_cluster_summary.sort_values("bin_rank").drop(columns=["bin_rank"])

        df_ds_summary = (
            df_full.group_by("n_datasets")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"))
            .sort("n_datasets")
            .to_pandas()
        )

        r_a_p, r_a_s = get_panel_corr_lines(df_full["resolution"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())
        r_b_p, r_b_s = get_panel_corr_lines(df_full["n_clusters"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())
        r_c_p, r_c_s = get_panel_corr_lines(df_full["n_datasets"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())

        ref_order = sorted(df_full["reference_type"].unique().to_list())
        color_scale = alt.Scale(scheme="tableau20", domain=ref_order)

        ref_legend = alt.Legend(
            orient="bottom",
            direction="horizontal",
            columns=5,
            columnPadding=50,
            titleAnchor="middle",
            symbolLimit=100,
            labelLimit=250,
            titleLimit=400,
            title="Reference Dataset / Strategy",
            labelFontSize=11.5,
            titleFontSize=12.5,
            titleFontWeight="bold",
        )

        PANEL_W = 285
        PANEL_H = 220

        base = alt.Chart(df_full_pd)

        def make_panel_letter(letter: str) -> Any:
            df_lbl = pd.DataFrame({"letter": [letter]})
            return (
                alt.Chart(df_lbl)
                .mark_text(
                    align="left",
                    baseline="bottom",
                    fontSize=17,
                    fontWeight="bold",
                    color="#0F172A",
                    clip=False,
                )
                .encode(
                    x=alt.value(0),
                    y=alt.value(-14),
                    text="letter:N",
                )
            )

        def make_panel_stats(p_stat: str, s_stat: str, x_px: float = 12.0, y_px: float = 10.0) -> Any:
            p_short = p_stat.replace("Pearson ", "").replace("Spearman ", "")
            s_short = s_stat.replace("Pearson ", "").replace("Spearman ", "")
            df_stats = pd.DataFrame({
                "stat": [p_short, s_short],
                "y_offset": [y_px, y_px + 14.0],
            })
            return (
                alt.Chart(df_stats)
                .mark_text(
                    align="left",
                    baseline="top",
                    fontSize=11.0,
                    fontWeight="bold",
                    color="#1E293B",
                    clip=False,
                )
                .encode(
                    x=alt.value(x_px),
                    y=alt.Y("y_offset:Q", scale=None),
                    text="stat:N",
                )
            )

        panel_a = (
            alt.layer(
                base.mark_boxplot(size=24, outliers=False, opacity=0.8).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=34, opacity=0.6).encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    xOffset="res_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                    tooltip=["reference_type", "reference_category", "resolution", "mean_cohort_multivariate_auc", "n_clusters"],
                ),
                alt.Chart(df_res_summary).mark_line(color="#e41a1c", strokeWidth=2.8, point=alt.OverlayMarkDef(size=55, color="#e41a1c")).encode(
                    x=alt.X("resolution:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                make_panel_letter("(a)"),
                make_panel_stats(r_a_p, r_a_s, x_px=14, y_px=14),
            )
            .properties(width=PANEL_W, height=PANEL_H)
        )

        panel_b = (
            alt.layer(
                base.mark_boxplot(size=24, outliers=False, opacity=0.8).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title=None, scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=34, opacity=0.6).encode(
                    x=alt.X("cluster_bin:O", title="Number of Cell States (Regimes)", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    xOffset="cluster_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_clusters", "resolution", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_cluster_summary).mark_line(color="#e41a1c", strokeWidth=2.8, point=alt.OverlayMarkDef(size=55, color="#e41a1c")).encode(
                    x=alt.X("cluster_bin:O", sort=bin_order, axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                make_panel_letter("(b)"),
                make_panel_stats(r_b_p, r_b_s, x_px=14, y_px=14),
            )
            .properties(width=PANEL_W, height=PANEL_H)
        )

        panel_c = (
            alt.layer(
                base.mark_boxplot(size=24, outliers=False, opacity=0.8).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title=None, scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=34, opacity=0.6).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets", axis=alt.Axis(labelAngle=0)),
                    xOffset="dataset_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_ds_summary).mark_line(color="#e41a1c", strokeWidth=2.8, point=alt.OverlayMarkDef(size=55, color="#e41a1c")).encode(
                    x=alt.X("n_datasets:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                make_panel_letter("(c)"),
                make_panel_stats(r_c_p, r_c_s, x_px=26, y_px=14),
            )
            .properties(width=PANEL_W, height=PANEL_H)
        )

        row1 = alt.hconcat(panel_a, panel_b, panel_c, spacing=24).resolve_scale(color="shared")

        # Row 2 (Panel d: PR-AUC) & Row 3 (Panel e: ROC-AUC): 90°-rotated 32-biomarker forest plots
        iatlas_paths = (
            REPO_ROOT / "output/manuscript_pipeline/3e5_integrate_deconvolution_our/iatlas_benchmark_metrics.parquet",
            REPO_ROOT / "output/benchmarks/iatlas_roc_auc/iatlas_benchmark_metrics.parquet",
        )
        iatlas_parquet = next((p for p in iatlas_paths if p.exists()), iatlas_paths[0])
        if not iatlas_parquet.exists():
            return Failure(f"iatlas_benchmark_metrics.parquet not found at {iatlas_parquet}")

        raw_bench_df = pl.read_parquet(iatlas_parquet)
        pancan_bench = raw_bench_df.filter(
            (pl.col("cohort_id") == "PanCancer-Combined")
            & (pl.col("time_stratum") == "Pre")
            & (pl.col("response_stratum") == "Standard")
            & (pl.col("pooling_strategy") == "standardized")
        )
        enr_bench = with_display_labels(pancan_bench)
        pred_col = "predictor_display" if "predictor_display" in enr_bench.columns else "predictor_name"
        pr_summary_df = (
            enr_bench.group_by([pred_col, "category"])
            .agg(
                pl.col("pr_auc").mean().alias("pr_mean"),
                pl.col("pr_auc_ci_lower").mean().alias("pr_lower"),
                pl.col("pr_auc_ci_upper").mean().alias("pr_upper"),
                pl.col("roc_auc").mean().alias("roc_mean"),
                pl.col("roc_auc_ci_lower").mean().alias("roc_lower"),
                pl.col("roc_auc_ci_upper").mean().alias("roc_upper"),
                pl.col("baseline_prevalence").mean().alias("baseline_prevalence"),
            )
            .sort("pr_mean", descending=True)
        )
        predictor_order = pr_summary_df[pred_col].to_list()
        ref_val_pr = float(pr_summary_df["baseline_prevalence"].mean())
        ref_val_roc = 0.50
        pr_pd = pr_summary_df.to_pandas()

        cat_domain = ["Cell Deconvolution", "Composite Synergy", "Genomic", "Signature"]
        cat_colors = ["#4c78a8", "#f58518", "#e45756", "#72b7b2"]
        cat_scale = alt.Scale(domain=cat_domain, range=cat_colors)
        cat_legend = alt.Legend(
            orient="none",
            legendX=297,
            legendY=366,
            direction="horizontal",
            title="Biomarker Category",
            titleAnchor="middle",
            labelFontSize=11.5,
            titleFontSize=12.5,
            titleFontWeight="bold",
            symbolSize=95,
            columnPadding=18,
        )

        ROW2_W = 971
        FOREST_H = 180

        base_pr = alt.Chart(pr_pd).encode(
            x=alt.X(
                f"{pred_col}:N",
                title=None,
                sort=predictor_order,
                axis=alt.Axis(labels=False, ticks=True),
            ),
            color=alt.Color("category:N", title="Biomarker Category", scale=cat_scale, legend=None),
        )

        is_our_pred = alt.datum[pred_col] == "Deconvolution Inference (Our)"
        is_our_tick = alt.datum.value == "Deconvolution Inference (Our)"

        error_bars_pr = base_pr.mark_rule().encode(
            y=alt.Y(
                "pr_lower:Q",
                title="Mean PR-AUC (95% CI)",
                scale=alt.Scale(domain=[0.25, 0.60]),
                axis=alt.Axis(
                    values=[0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60],
                ),
            ),
            y2="pr_upper:Q",
            strokeWidth=alt.condition(is_our_pred, alt.value(3.0), alt.value(2.2)),
        )

        points_pr = base_pr.mark_point(filled=True, opacity=1.0).encode(
            y=alt.Y("pr_mean:Q", scale=alt.Scale(domain=[0.25, 0.60])),
            size=alt.condition(is_our_pred, alt.value(125), alt.value(80)),
            stroke=alt.condition(is_our_pred, alt.value("#0F172A"), alt.value("transparent")),
            strokeWidth=alt.condition(is_our_pred, alt.value(1.2), alt.value(0.0)),
            tooltip=[pred_col, "category", "pr_mean", "pr_lower", "pr_upper"],
        )

        ref_line_pr = (
            alt.Chart(pd.DataFrame({"y": [ref_val_pr]}))
            .mark_rule(strokeDash=[4, 4], color="#64748B", strokeWidth=1.5)
            .encode(y=alt.Y("y:Q", scale=alt.Scale(domain=[0.25, 0.60])))
        )

        panel_d = (
            alt.layer(ref_line_pr, error_bars_pr, points_pr, make_panel_letter("(d)"))
            .properties(width=ROW2_W, height=FOREST_H)
        )

        base_roc = alt.Chart(pr_pd).encode(
            x=alt.X(
                f"{pred_col}:N",
                title="Predictor / Biomarker",
                sort=predictor_order,
                axis=alt.Axis(
                    labelAngle=-40,
                    labelLimit=240,
                    labelFontSize=11.0,
                    labelFontWeight=alt.condition(is_our_tick, alt.value("bold"), alt.value("normal")),
                    labelColor=alt.condition(is_our_tick, alt.value("#000000"), alt.value("#334155")),
                    titleFontSize=13.0,
                    titleFontWeight="bold",
                    titlePadding=10,
                ),
            ),
            color=alt.Color("category:N", title="Biomarker Category", scale=cat_scale, legend=cat_legend),
        )

        error_bars_roc = base_roc.mark_rule().encode(
            y=alt.Y(
                "roc_lower:Q",
                title="Mean ROC-AUC (95% CI)",
                scale=alt.Scale(domain=[0.45, 0.80]),
                axis=alt.Axis(
                    values=[0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80],
                ),
            ),
            y2="roc_upper:Q",
            strokeWidth=alt.condition(is_our_pred, alt.value(3.0), alt.value(2.2)),
        )

        points_roc = base_roc.mark_point(filled=True, opacity=1.0).encode(
            y=alt.Y("roc_mean:Q", scale=alt.Scale(domain=[0.45, 0.80])),
            size=alt.condition(is_our_pred, alt.value(125), alt.value(80)),
            stroke=alt.condition(is_our_pred, alt.value("#0F172A"), alt.value("transparent")),
            strokeWidth=alt.condition(is_our_pred, alt.value(1.2), alt.value(0.0)),
            tooltip=[pred_col, "category", "roc_mean", "roc_lower", "roc_upper"],
        )

        ref_line_roc = (
            alt.Chart(pd.DataFrame({"y": [ref_val_roc]}))
            .mark_rule(strokeDash=[4, 4], color="#64748B", strokeWidth=1.5)
            .encode(y=alt.Y("y:Q", scale=alt.Scale(domain=[0.45, 0.80])))
        )

        panel_e = (
            alt.layer(ref_line_roc, error_bars_roc, points_roc, make_panel_letter("(e)"))
            .properties(width=ROW2_W, height=FOREST_H)
        )

        final_chart = (
            alt.vconcat(row1, panel_d, panel_e, spacing=26)
            .resolve_scale(color="independent")
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 24, "right": 15, "top": 28, "bottom": 18},
            )
            .configure_axisY(titlePadding=10)
            .configure_axis(labelFontSize=12.0, titleFontSize=13.0, titleFontWeight="bold")
            .configure_legend(labelFontSize=11.5, titleFontSize=12.5, titleFontWeight="bold")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
        )

        out_path = results_dir / "fig_sc9_resolution_states_datasets_summary.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_path))

        post_process_svg(out_path, n_panels=0, width_mm=180.0, min_font_pt=12.0)

        png_path = results_dir / "fig_sc9_resolution_states_datasets_summary.png"
        rasterize_svg_to_png(out_path, png_path, scale=3.125, background="white")

        return Success(out_path)

    except Exception as exc:
        return Failure(f"Failed to generate Summary Figure: {exc}")


def plot_fig2_signature_collinearity(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 2: Signature collinearity benchmark (condition number vs resolution and cell states)."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure("multi_resolution_benchmark_summary.parquet not found")
    try:
        df = pl.read_parquet(bench_file).filter(pl.col("condition_number").is_not_null())
        df_full = df.filter(pl.col("reference_category") != "Cell-Subsample")
        df_pd = df_full.to_pandas()

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]
        f2_a_p, f2_a_s = get_panel_corr_lines(df_full["resolution"].to_numpy(), df_full["condition_number"].to_numpy(), y_log=True)
        f2_b_p, f2_b_s = get_panel_corr_lines(df_full["n_clusters"].to_numpy(), df_full["condition_number"].to_numpy(), y_log=True)

        panel_a = (
            alt.Chart(df_pd)
            .mark_line(point=alt.OverlayMarkDef(size=40, filled=True), strokeWidth=1.8)
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("condition_number:Q", title="Condition Number (κ, log scale)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                detail="reference_type:N",
                tooltip=["reference_type", "resolution", "condition_number", "n_clusters"],
            )
            .properties(title=["Matrix Condition Number vs Resolution", f2_a_p, f2_a_s], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        panel_b = (
            alt.Chart(df_pd)
            .mark_circle(size=45, opacity=0.7)
            .encode(
                x=alt.X("n_clusters:Q", title="Number of Cell States (Clusters)"),
                y=alt.Y("condition_number:Q", title="Condition Number (κ, log scale)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                tooltip=["reference_type", "n_clusters", "condition_number"],
            )
            .properties(title=["Matrix Collinearity vs Cell States", f2_b_p, f2_b_s], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        final_chart = (
            alt.hconcat(panel_a, panel_b, spacing=25)
            .resolve_scale(color="shared")
            .configure_axis(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_legend(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_title(fontSize=15, fontWeight="bold")
        )

        out_svg = results_dir / "fig2_signature_collinearity_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))
        post_process_svg(out_svg, n_panels=2, width_mm=180.0, min_font_pt=13.0)

        png_path = results_dir / "fig2_signature_collinearity_benchmark.png"
        rasterize_svg_to_png(out_svg, png_path, scale=3.125, background="white")

        return Success(out_svg)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 2: {exc}")


def plot_fig3_reference_selection_decision_framework(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 3: Classifier selection and optimization decision framework."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure("multi_resolution_benchmark_summary.parquet not found")
    try:
        df = pl.read_parquet(bench_file)
        if "subsample_fraction" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("subsample_fraction"))

        df_full = df.filter(pl.col("reference_category") != "Cell-Subsample")

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]

        # Panel A: Strategy x Resolution Grid
        grid_df = (
            df_full.group_by(["reference_category", "resolution"])
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_cohort_auc"))
            .sort(["reference_category", "resolution"])
            .to_pandas()
        )
        grid_base = alt.Chart(grid_df)
        panel_a = (
            alt.layer(
                grid_base.mark_rect().encode(
                    x=alt.X("resolution:O", title="Leiden Resolution", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("reference_category:N", title="Strategy", sort=cat_order),
                    color=alt.Color("mean_cohort_auc:Q", title="Cohort AUC", scale=alt.Scale(scheme="viridis", domain=[0.60, 0.74])),
                ),
                grid_base.mark_text(baseline="middle", fontSize=13, fontWeight="bold").encode(
                    x="resolution:O",
                    y=alt.Y("reference_category:N", sort=cat_order),
                    text=alt.Text("mean_cohort_auc:Q", format=".3f"),
                    color=alt.condition(alt.datum.mean_cohort_auc > 0.69, alt.value("white"), alt.value("black")),
                ),
            )
            .properties(title=["Strategy × Resolution", "Cohort-Averaged AUC Grid"], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        # Panel B: Pareto Trade-off Frontier
        df_full_pd = df_full.to_pandas()
        pareto_df = df_full.filter(pl.col("is_pareto_optimal")).sort("condition_number").to_pandas() if "is_pareto_optimal" in df_full.columns else df_full_pd.head(5)
        pareto_chart = alt.Chart(pareto_df)

        panel_b = (
            alt.layer(
                alt.Chart(df_full_pd).mark_circle(size=40, opacity=0.45).encode(
                    x=alt.X("condition_number:Q", title="Signature Condition Number (κ, log scale)", scale=alt.Scale(type="log", domain=[1.5, 120])),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                    tooltip=["reference_type", "resolution", "condition_number", "mean_cohort_multivariate_auc"],
                ),
                pareto_chart.mark_line(color="#e41a1c", strokeWidth=2.2, strokeDash=[4, 2]).encode(
                    x="condition_number:Q",
                    y="mean_cohort_multivariate_auc:Q",
                ),
                pareto_chart.mark_circle(size=85, color="#e41a1c", filled=True).encode(
                    x="condition_number:Q",
                    y="mean_cohort_multivariate_auc:Q",
                    tooltip=["reference_type", "resolution", "condition_number", "mean_cohort_multivariate_auc"],
                ),
            )
            .properties(title=["Pareto Trade-off Frontier", "AUC vs Matrix Collinearity"], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        # Panel C: Top 10 Configurations Generalizability
        cohort_cols = [c for c in df.columns if c.startswith("auc_")]
        top10 = (
            df_full.sort("mean_cohort_multivariate_auc", descending=True)
            .head(10)
            .with_columns(
                (pl.col("reference_type") + " (res " + pl.col("resolution").cast(pl.String) + ")").alias("config_label")
            )
        )
        config_order = top10["config_label"].to_list()
        top10_melt = (
            top10.select(["config_label", "mean_cohort_multivariate_auc"] + cohort_cols)
            .unpivot(
                index=["config_label", "mean_cohort_multivariate_auc"],
                on=cohort_cols,
                variable_name="cohort_raw",
                value_name="auc",
            )
            .with_columns(
                pl.col("cohort_raw").str.replace("auc_", "").str.to_titlecase().alias("cohort")
            )
            .to_pandas()
        )
        c_base = alt.Chart(top10_melt)
        panel_c = (
            alt.layer(
                c_base.mark_rect().encode(
                    x=alt.X("cohort:N", title="Bulk ICB Cohort", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("config_label:N", title="Reference Configuration", sort=config_order),
                    color=alt.Color("auc:Q", title="Cohort AUC", scale=alt.Scale(scheme="plasma", domain=[0.45, 0.90])),
                ),
                c_base.mark_text(baseline="middle", fontSize=11, fontWeight="bold").encode(
                    x="cohort:N",
                    y=alt.Y("config_label:N", sort=config_order),
                    text=alt.Text("auc:Q", format=".2f"),
                    color=alt.condition(alt.datum.auc > 0.70, alt.value("white"), alt.value("black")),
                ),
            )
            .properties(title=["Multi-Cohort Generalizability", "Top 10 Reference Models across 9 Cohorts"], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        # Panel D: Cohort Diversity Scaling Curve
        f3_d_p, f3_d_s = get_panel_corr_lines(df_full["n_datasets"].to_numpy(), df_full["mean_cohort_multivariate_auc"].to_numpy())
        df_ds_summary = (
            df_full.group_by("n_datasets")
            .agg(pl.col("mean_cohort_multivariate_auc").mean().alias("mean_coh_auc"))
            .sort("n_datasets")
            .to_pandas()
        )
        df_full_pd["ds_jitter"] = np.random.uniform(-8.0, 8.0, len(df_full_pd))

        panel_d = (
            alt.layer(
                alt.Chart(df_full_pd).mark_boxplot(size=22, outliers=False, opacity=0.8).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets (K)", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                alt.Chart(df_full_pd).mark_circle(size=32, opacity=0.55).encode(
                    x=alt.X("n_datasets:O", title="Integrated Single-Cell Datasets (K)", axis=alt.Axis(labelAngle=0)),
                    xOffset="ds_jitter:Q",
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_category:N", title="Strategy", sort=cat_order, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "mean_cohort_multivariate_auc"],
                ),
                alt.Chart(df_ds_summary).mark_line(color="#e41a1c", strokeWidth=2.5, point=alt.OverlayMarkDef(size=50, color="#e41a1c")).encode(
                    x=alt.X("n_datasets:O", axis=alt.Axis(labelAngle=0)),
                    y=alt.Y("mean_coh_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
            )
            .properties(title=["Cohort Diversity Scaling", f3_d_p, f3_d_s], width=PANEL_W_2COL, height=PANEL_H_2COL)
        )

        row1 = (panel_a | panel_b).resolve_scale(color="independent")
        row2 = (panel_c | panel_d).resolve_scale(color="independent")
        final_chart = (
            alt.vconcat(row1, row2, spacing=25)
            .resolve_scale(color="independent")
            .configure_axis(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_legend(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_title(fontSize=15, fontWeight="bold")
        )

        out_svg = results_dir / "fig3_reference_selection_decision_framework.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))
        post_process_svg(out_svg, n_panels=4, width_mm=180.0, min_font_pt=13.0)

        png_path = results_dir / "fig3_reference_selection_decision_framework.png"
        rasterize_svg_to_png(out_svg, png_path, scale=3.125, background="white")

        return Success(out_svg)

    except Exception as exc:
        return Failure(f"Failed to generate Figure 3: {exc}")


def plot_fig4_winning_classifier_biomarker_drivers(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 4: Winning classifier biomarker drivers (odds ratios and Milo DA concordance)."""
    # Select winning model: all_9_combined_res0.75
    lr_file = data_dir.parent / "03_logistic_regression" / "logistic_regression_results_all_9_combined_res0.75.parquet"
    milo_file = data_dir.parent / "04_milopy" / "milopy_cell_state_da_all_resolutions.parquet"
    mg_file = data_dir.parent / "01_references" / "reference_marker_genes_all_9_combined_res0.75.parquet"
    deconv_file = data_dir.parent / "02_deconvolution" / "deconv_fractions_all_9_combined_res0.75.parquet"

    if not lr_file.exists():
        fallback_lr = Path("output/sc9_milo_deconv_validation/03_logistic_regression/logistic_regression_results_all_9_combined_res0.75.parquet")
        if fallback_lr.exists():
            lr_file = fallback_lr
            milo_file = Path("output/sc9_milo_deconv_validation/04_milopy/milopy_cell_state_da_all_resolutions.parquet")
            mg_file = Path("output/sc9_milo_deconv_validation/01_references/reference_marker_genes_all_9_combined_res0.75.parquet")
            deconv_file = Path("output/sc9_milo_deconv_validation/02_deconvolution/deconv_fractions_all_9_combined_res0.75.parquet")
        else:
            return Failure(f"Required logistic regression file not found: {lr_file}")

    try:
        lr = pl.read_parquet(lr_file).filter(pl.col("stratum") == "Melanoma")
        milo = pl.read_parquet(milo_file).filter(pl.col("resolution") == 0.75)
        mg = (
            pl.read_parquet(mg_file)
            .filter(pl.col("rank") <= 2)
            .group_by("cluster")
            .agg(pl.col("gene").str.join("+"))
            .rename({"cluster": "cell_state", "gene": "markers"})
        )

        stds_dict: dict[str, float] = {}
        if deconv_file.exists():
            deconv = pl.read_parquet(deconv_file).filter(pl.col("cancer_type") == "Melanoma")
            for cs in lr["cell_state"].to_list():
                if cs in deconv.columns:
                    stds_dict[cs] = max(1e-4, float(deconv[cs].std()))

        joined = lr.join(milo, on="cell_state").join(mg, on="cell_state")
        records: list[dict[str, object]] = []
        for row in joined.iter_rows(named=True):
            cs = str(row["cell_state"])
            markers = str(row["markers"])
            b = float(row["beta"])
            se = float(row["se"]) if not np.isnan(row["se"]) else 1.0
            p_val = float(row["p_value"])
            m_fc = float(row["milo_mean_logfc"])
            s = stds_dict.get(cs, 0.05)

            # Standardized Odds Ratio per 1 SD increase
            or_sd = float(np.exp(np.clip(b * s, -6, 6)))
            or_ci_l = float(np.exp(np.clip((b - 1.96 * se) * s, -6, 6)))
            or_ci_u = float(np.exp(np.clip((b + 1.96 * se) * s, -6, 6)))

            if p_val < 0.05 and or_sd > 1.0:
                sig_cat = "Responder Enriched (p < 0.05)"
            elif p_val < 0.05 and or_sd < 1.0:
                sig_cat = "Non-Responder Enriched (p < 0.05)"
            else:
                sig_cat = "Not Significant"

            label = f"CS{cs}: {markers}"
            records.append({
                "cell_state": cs,
                "markers": markers,
                "label": label,
                "beta": b,
                "se": se,
                "p_value": p_val,
                "or_sd": or_sd,
                "or_ci_lower": or_ci_l,
                "or_ci_upper": or_ci_u,
                "milo_mean_logfc": m_fc,
                "sig_category": sig_cat,
            })

        df_p = pd.DataFrame(records).sort_values("or_sd", ascending=False).reset_index(drop=True)
        label_order = df_p["label"].tolist()

        color_scale = alt.Scale(
            domain=["Responder Enriched (p < 0.05)", "Non-Responder Enriched (p < 0.05)", "Not Significant"],
            range=["#2b83ba", "#d7191c", "#94a3b8"],
        )

        ref_line = alt.Chart(pd.DataFrame({"x": [1.0]})).mark_rule(strokeDash=[4, 4], color="#475569", strokeWidth=1.5).encode(x="x:Q")

        chart_base = alt.Chart(df_p)
        panel_a = (
            alt.layer(
                ref_line,
                chart_base.mark_rule(strokeWidth=2.0, color="#64748b").encode(
                    y=alt.Y("label:N", title="Cell State (Cluster & Markers)", sort=label_order),
                    x=alt.X("or_ci_lower:Q", title="Standardized Odds Ratio (per 1 SD, log scale)", scale=alt.Scale(type="log", domain=[0.05, 20.0])),
                    x2="or_ci_upper:Q",
                ),
                chart_base.mark_circle(size=80, opacity=1.0).encode(
                    y=alt.Y("label:N", sort=label_order),
                    x=alt.X("or_sd:Q", scale=alt.Scale(type="log", domain=[0.05, 20.0])),
                    color=alt.Color("sig_category:N", title="Association", scale=color_scale, legend=alt.Legend(orient="bottom", direction="vertical")),
                    tooltip=["label", "or_sd", "or_ci_lower", "or_ci_upper", "p_value"],
                ),
            )
            .properties(title=["Cell State Odds Ratios (95% CI)", "Melanoma ICB Response (All-9 Reference)"], width=PANEL_W_2COL, height=340)
        )

        # Panel B: Milo DA vs Deconvolution Effect Size Concordance
        b_vals = df_p["beta"].to_numpy().astype(np.float64)
        m_vals = df_p["milo_mean_logfc"].to_numpy().astype(np.float64)
        rho_res = stats.spearmanr(b_vals, m_vals)
        rho = float(rho_res.statistic)
        rho_p = float(rho_res.pvalue)

        lin_res = stats.linregress(m_vals, b_vals)
        slope = float(lin_res.slope)
        intercept = float(lin_res.intercept)
        x_grid = np.linspace(float(np.min(m_vals)), float(np.max(m_vals)), 50)
        df_trend = pd.DataFrame({"milo_mean_logfc": x_grid, "beta": slope * x_grid + intercept})

        quad_x = alt.Chart(pd.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[3, 3], color="#cbd5e1", strokeWidth=1.2).encode(x="x:Q")
        quad_y = alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(strokeDash=[3, 3], color="#cbd5e1", strokeWidth=1.2).encode(y="y:Q")

        df_p["concordant"] = ((df_p["beta"] > 0) & (df_p["milo_mean_logfc"] > 0)) | ((df_p["beta"] < 0) & (df_p["milo_mean_logfc"] < 0))
        conc_color = alt.Scale(domain=[True, False], range=["#2b83ba", "#d7191c"])

        f4_b_p, f4_b_s = get_panel_corr_lines(m_vals, b_vals)

        panel_b = (
            alt.layer(
                quad_x,
                quad_y,
                alt.Chart(df_trend).mark_line(color="#1e293b", strokeWidth=2.0).encode(
                    x="milo_mean_logfc:Q",
                    y="beta:Q",
                ),
                chart_base.mark_circle(size=85, opacity=0.85).encode(
                    x=alt.X("milo_mean_logfc:Q", title="Single-Cell Milo DA log₂(Fold Change)"),
                    y=alt.Y("beta:Q", title="Bulk Deconvolution Logistic Regression β"),
                    color=alt.Color("concordant:N", title="Concordant Direction", scale=conc_color, legend=None),
                    tooltip=["label", "milo_mean_logfc", "beta", "p_value"],
                ),
                chart_base.mark_text(dx=8, dy=-4, fontSize=11, fontStyle="italic").encode(
                    x="milo_mean_logfc:Q",
                    y="beta:Q",
                    text="cell_state:N",
                ),
            )
            .properties(
                title=[
                    "Single-Cell DA vs Bulk β Concordance",
                    f4_b_p,
                    f4_b_s,
                ],
                width=PANEL_W_2COL,
                height=340,
            )
        )

        final_chart = (
            alt.hconcat(panel_a, panel_b, spacing=30)
            .resolve_scale(color="independent")
            .configure_axis(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_legend(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_title(fontSize=15, fontWeight="bold")
        )

        out_svg = results_dir / "fig4_winning_classifier_biomarker_drivers.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))
        post_process_svg(out_svg, n_panels=2, width_mm=180.0, min_font_pt=13.0)

        png_path = results_dir / "fig4_winning_classifier_biomarker_drivers.png"
        rasterize_svg_to_png(out_svg, png_path, scale=3.125, background="white")

        return Success(out_svg)

    except Exception as exc:
        return Failure(f"Failed to generate Figure 4: {exc}")


def plot_fig5_selected_cohort_auc(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 5: Selected cohort AUC and marginal gain per resolution step with statistical testing."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure("multi_resolution_benchmark_summary.parquet not found")
    try:
        df = pl.read_parquet(bench_file)
        df_full = df.filter(pl.col("reference_category") != "Cell-Subsample")
        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]
        colors = ["#4c78a8", "#f58518", "#54a24b"]
        color_scale = alt.Scale(domain=cat_order, range=colors)

        # 1. Compute linear regression slopes of AUC vs resolution for each strategy
        stats_records = []
        trend_records = []
        x_grid = np.linspace(0.25, 2.0, 50)

        for cat in cat_order:
            sub = df_full.filter(pl.col("reference_category") == cat)
            x = sub["resolution"].to_numpy()
            y = sub["mean_cohort_multivariate_auc"].to_numpy()
            res = stats.linregress(x, y)

            p_formatted = f"p = {res.pvalue:.1e}" if res.pvalue < 0.001 else f"p = {res.pvalue:.3f}"
            stat_label = f"{cat} (β = {res.slope:+.3f}, {p_formatted})"
            stats_records.append({
                "reference_category": cat,
                "slope": res.slope,
                "intercept": res.intercept,
                "rvalue": res.rvalue,
                "pvalue": res.pvalue,
                "stat_label": stat_label,
                "end_x": 2.0,
                "end_y": float(res.slope * 2.0 + res.intercept),
            })

            for xg in x_grid:
                trend_records.append({
                    "reference_category": cat,
                    "resolution": float(xg),
                    "mean_cohort_multivariate_auc": float(res.slope * xg + res.intercept),
                })

        df_stats = pd.DataFrame(stats_records)
        df_trend = pd.DataFrame(trend_records)

        # Calculate overall regression slope across all models
        x_all = df_full["resolution"].to_numpy()
        y_all = df_full["mean_cohort_multivariate_auc"].to_numpy()
        res_all = stats.linregress(x_all, y_all)
        all_p_str = f"p = {res_all.pvalue:.1e}" if res_all.pvalue < 0.001 else f"p = {res_all.pvalue:.3f}"
        f5_a_p, f5_a_s = get_panel_corr_lines(x_all, y_all)

        # Panel A: Trajectories + Regression Slopes + Slope Text Annotations
        raw_chart = (
            alt.Chart(df_full.to_pandas())
            .mark_line(point=alt.OverlayMarkDef(size=25, opacity=0.35), strokeWidth=1.0, opacity=0.25)
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort Response AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                color=alt.Color("reference_category:N", title="Strategy", scale=color_scale, legend=None),
                detail="reference_type:N",
                tooltip=["reference_type", "reference_category", "resolution", "mean_cohort_multivariate_auc"],
            )
        )

        trend_chart = (
            alt.Chart(df_trend)
            .mark_line(strokeWidth=3.2)
            .encode(
                x="resolution:Q",
                y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                color=alt.Color("reference_category:N", scale=color_scale, legend=None),
            )
        )

        text_annotations = (
            alt.Chart(df_stats)
            .mark_text(align="right", dx=-8, dy=-6, fontSize=11, fontWeight="bold")
            .encode(
                x="end_x:Q",
                y=alt.Y("end_y:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                text="stat_label:N",
                color=alt.Color("reference_category:N", scale=color_scale, legend=None),
            )
        )

        panel_a = (
            alt.layer(raw_chart, trend_chart, text_annotations)
            .properties(
                title=[
                    "Cohort-Averaged Predictive AUC across Resolutions",
                    f"Overall: β = {res_all.slope:+.3f} ({all_p_str})",
                    f"{f5_a_p} | {f5_a_s}",
                ],
                width=PANEL_W_2COL,
                height=PANEL_H_2COL,
            )
        )

        # Panel B: Marginal AUC Gain by Strategy per Resolution Step with p-values
        df_delta = (
            df_full.sort("resolution")
            .with_columns(
                pl.col("mean_cohort_multivariate_auc")
                .shift(1)
                .over("reference_type")
                .alias("prev_auc")
            )
            .with_columns(
                (pl.col("mean_cohort_multivariate_auc") - pl.col("prev_auc")).alias("delta_auc")
            )
            .filter(pl.col("delta_auc").is_not_null())
        )

        df_delta_agg = (
            df_delta.group_by(["reference_category", "resolution"])
            .agg(pl.col("delta_auc").mean().alias("mean_delta_auc"))
            .sort(["reference_category", "resolution"])
            .to_pandas()
        )

        # Compute step-wise p-values across all strategies combined
        step_stats = []
        for res_step in sorted(df_delta["resolution"].unique().to_list()):
            sub_step = df_delta.filter(pl.col("resolution") == res_step)
            vals = sub_step["delta_auc"].to_numpy()
            _, p_val = stats.ttest_1samp(vals, 0.0)
            mean_v = float(vals.mean())
            if p_val < 0.001:
                star = "***"
            elif p_val < 0.01:
                star = "**"
            elif p_val < 0.05:
                star = "*"
            else:
                star = "ns"
            step_stats.append({
                "resolution": float(res_step),
                "step_label": f"{res_step:.2f}",
                "pvalue": p_val,
                "significance": star,
                "y_pos": max(0.026, mean_v + 0.005),
            })

        df_step_stats = pd.DataFrame(step_stats)

        bar_chart = (
            alt.Chart(df_delta_agg)
            .mark_bar()
            .encode(
                x=alt.X("resolution:O", title="Leiden Resolution Step", axis=alt.Axis(labelAngle=0)),
                y=alt.Y("mean_delta_auc:Q", title="Mean ΔAUC vs Prev. Resolution"),
                color=alt.Color(
                    "reference_category:N",
                    title="Strategy",
                    sort=cat_order,
                    scale=color_scale,
                    legend=alt.Legend(orient="bottom", direction="horizontal", titleFontSize=13, labelFontSize=12),
                ),
                xOffset=alt.XOffset("reference_category:N", sort=cat_order),
                tooltip=["reference_category", "resolution", "mean_delta_auc"],
            )
        )

        step_text = (
            alt.Chart(df_step_stats)
            .mark_text(baseline="bottom", dy=-3, fontSize=12, fontWeight="bold", color="#1e293b")
            .encode(
                x=alt.X("resolution:O", axis=alt.Axis(labelAngle=0)),
                y=alt.Y("y_pos:Q"),
                text="significance:N",
            )
        )

        zero_rule = alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(color="#94a3b8", strokeWidth=1.0).encode(y="y:Q")

        panel_b = (
            alt.layer(bar_chart, step_text, zero_rule)
            .properties(
                title=[
                    "Marginal AUC Gain by Strategy per Resolution Step",
                    "Step Significance: *** p < 0.001, ** p < 0.01, * p < 0.05",
                ],
                width=PANEL_W_2COL,
                height=PANEL_H_2COL,
            )
        )

        final_chart = (
            alt.hconcat(panel_a, panel_b, spacing=25)
            .resolve_scale(color="shared")
            .configure_axis(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_legend(labelFontSize=13, titleFontSize=14, titleFontWeight="bold")
            .configure_title(fontSize=14, fontWeight="bold")
        )

        out_svg = results_dir / "fig5_selected_cohort_auc_and_marginal_gain.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))
        post_process_svg(out_svg, n_panels=2, width_mm=180.0, min_font_pt=13.0)

        png_path = results_dir / "fig5_selected_cohort_auc_and_marginal_gain.png"
        rasterize_svg_to_png(out_svg, png_path, scale=3.125, background="white")

        return Success(out_svg)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 5: {exc}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6: Generate publication benchmark figures for reference validation."
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/05_concordance",
        help="Directory containing multi_resolution_benchmark_summary.parquet",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="results/sc9_milo_deconv_validation",
        help="Directory to save generated publication figures",
    )
    args = parser.parse_args()

    data_path = Path(args.data_dir)
    results_path = Path(args.results_dir)
    results_path.mkdir(parents=True, exist_ok=True)

    print(f"=== Generating Multi-Resolution Benchmark Figures into {results_path} ===")

    print("Generating Figure 1 (18-panel Predictive Capacity Benchmark with Boxplots)...")
    match plot_fig1_predictive_capacity(data_path, results_path):
        case Success(p):
            print(f"  Saved Figure 1: {p}")
        case Failure(err):
            print(f"  Warning Figure 1: {err}", file=sys.stderr)

    print("Generating Summary Figure (Panels c, f, i: Resolution, Cell States, Datasets)...")
    match plot_fig_sc9_summary_panels_cfi(data_path, results_path):
        case Success(p):
            print(f"  Saved Summary Figure: {p}")
        case Failure(err):
            print(f"  Warning Summary Figure: {err}", file=sys.stderr)

    print("Generating Figure 2 (Signature Collinearity Benchmark)...")
    match plot_fig2_signature_collinearity(data_path, results_path):
        case Success(p):
            print(f"  Saved Figure 2: {p}")
        case Failure(err):
            print(f"  Warning Figure 2: {err}", file=sys.stderr)

    print("Generating Figure 3 (Reference Selection & Decision Framework)...")
    match plot_fig3_reference_selection_decision_framework(data_path, results_path):
        case Success(p):
            print(f"  Saved Figure 3: {p}")
        case Failure(err):
            print(f"  Warning Figure 3: {err}", file=sys.stderr)

    print("Generating Figure 4 (Winning Classifier Biomarker Drivers)...")
    match plot_fig4_winning_classifier_biomarker_drivers(data_path, results_path):
        case Success(p):
            print(f"  Saved Figure 4: {p}")
        case Failure(err):
            print(f"  Warning Figure 4: {err}", file=sys.stderr)

    print("Generating Figure 5 (Selected Cohort AUC & Marginal Gain)...")
    match plot_fig5_selected_cohort_auc(data_path, results_path):
        case Success(p):
            print(f"  Saved Figure 5: {p}")
        case Failure(err):
            print(f"  Warning Figure 5: {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
