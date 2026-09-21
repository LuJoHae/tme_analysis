#!/usr/bin/env python3
"""
Step 7: Generate Publication Visualizations (Altair Vector SVGs).
Creates publication-quality vector SVG charts for:
1. Fig 1: Concordance Scatter (Single-cell effect vs Bulk effect by quadrant)
2. Fig 2: Microenvironment Normalization Impact (Raw vs Purity-normalized vs mRNA-scaled)
3. Fig 3: Forest Plot of Cross-Cohort Effect Sizes
4. Fig 4: In Silico Deconvolution Recovery Calibration
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


class PlotConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    metrics_summary_path: Path
    state_details_path: Path
    insilico_path: Path
    results_dir: Path


def parse_args() -> PlotConfig:
    parser = argparse.ArgumentParser(
        description="Generate publication vector SVG figures for concordance analysis."
    )
    parser.add_argument("--metrics-summary", type=Path, required=True, help="Concordance metrics summary parquet")
    parser.add_argument("--state-details", type=Path, required=True, help="Concordance state details parquet")
    parser.add_argument("--insilico", type=Path, required=True, help="In silico recovery parquet")
    parser.add_argument("--results-dir", type=Path, required=True, help="Directory to save SVG files")
    args = parser.parse_args()
    return PlotConfig(
        metrics_summary_path=args.metrics_summary,
        state_details_path=args.state_details,
        insilico_path=args.insilico,
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


def plot_concordance_scatter(details_df: pl.DataFrame, summary_df: pl.DataFrame) -> alt.Chart:
    """Fig 1: Cross-cohort concordance scatter plot."""
    # Use normalized fractions
    norm_details = details_df.filter(pl.col("fraction_type") == "normalized")
    norm_summary = summary_df.filter(pl.col("fraction_type") == "normalized")

    rho_val = norm_summary["spearman_rho"][0] if norm_summary.height > 0 else 0.0
    pval = norm_summary["spearman_pval"][0] if norm_summary.height > 0 else 1.0

    color_scale = alt.Scale(
        domain=[
            "Concordant Responder",
            "Concordant Non-Responder",
            "Discordant (SC+, Bulk-)",
            "Discordant (SC-, Bulk+)",
        ],
        range=["#1b9e77", "#386cb0", "#d95f02", "#e7298a"],
    )

    scatter = (
        alt.Chart(norm_details)
        .mark_circle(size=120, opacity=0.85)
        .encode(
            x=alt.X("beta_sc:Q", title="Single-Cell Response Effect (CLR log-odds β)"),
            y=alt.Y("beta_bulk:Q", title="Bulk Deconvolution Response Effect (β)"),
            color=alt.Color("quadrant:N", scale=color_scale, legend=alt.Legend(title="Concordance Status")),
            tooltip=["cell_state:N", "beta_sc:Q", "beta_bulk:Q", "quadrant:N"],
        )
    )

    # Reference lines (x=0, y=0)
    rule_x = alt.Chart(pl.DataFrame({"x": [0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(x="x:Q")
    rule_y = alt.Chart(pl.DataFrame({"y": [0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(y="y:Q")

    chart = (
        (scatter + rule_x + rule_y)
        .properties(
            title=f"Single-Cell vs. Bulk Deconvolution Concordance (Spearman ρ = {rho_val:.3f}, p = {pval:.3e})",
            width=500,
            height=450,
        )
        .configure_axis(grid=True, gridColor="#f0f0f0")
    )
    return chart


def plot_purity_adjustment_impact(summary_df: pl.DataFrame) -> alt.Chart:
    """Fig 2: Bar chart showing change in Spearman rho across fraction types."""
    bar = (
        alt.Chart(summary_df)
        .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
        .encode(
            x=alt.X("fraction_type:N", title="Deconvolution Fraction Representation", sort=["raw", "normalized", "mrna_scaled"]),
            y=alt.Y("spearman_rho:Q", title="Cross-Cohort Spearman Rank Correlation (ρ)"),
            color=alt.Color("fraction_type:N", legend=None, scale=alt.Scale(scheme="category10")),
            tooltip=["fraction_type:N", "spearman_rho:Q", "concordance_rate:Q", "binomial_pval:Q"],
        )
    )

    text = (
        alt.Chart(summary_df)
        .mark_text(dy=-10, fontSize=12, fontWeight="bold")
        .encode(
            x=alt.X("fraction_type:N", sort=["raw", "normalized", "mrna_scaled"]),
            y=alt.Y("spearman_rho:Q"),
            text=alt.Text("spearman_rho:Q", format=".3f"),
        )
    )

    return (bar + text).properties(
        title="Impact of Tumor Purity Normalization and mRNA Scaling on Concordance",
        width=380,
        height=350,
    )


def plot_forest_effect_sizes(details_df: pl.DataFrame) -> alt.Chart:
    """Fig 3: Forest plot comparing single-cell vs bulk response associations."""
    norm_details = details_df.filter(pl.col("fraction_type") == "normalized")

    # Long format for comparison
    sc_subset = norm_details.select([
        pl.col("cell_state"),
        pl.col("beta_sc").alias("beta"),
        pl.col("se_sc").alias("se"),
        pl.lit("Single-Cell (Patient-level)").alias("modality"),
    ])
    bulk_subset = norm_details.select([
        pl.col("cell_state"),
        pl.col("beta_bulk").alias("beta"),
        pl.col("se_bulk").alias("se"),
        pl.lit("Bulk Deconvolution (Purity-norm)").alias("modality"),
    ])
    combined_long = pl.concat([sc_subset, bulk_subset]).with_columns([
        (pl.col("beta") - 1.96 * pl.col("se")).clip(-100.0, 100.0).alias("ci_lower"),
        (pl.col("beta") + 1.96 * pl.col("se")).clip(-100.0, 100.0).alias("ci_upper"),
    ])

    x_scale = alt.Scale(domain=[-100, 100], clamp=True)

    points = (
        alt.Chart(combined_long)
        .mark_point(filled=True, size=60)
        .encode(
            x=alt.X("beta:Q", title="Response Effect Size (log-odds β ± 95% CI)", scale=x_scale),
            y=alt.Y("cell_state:N", title="Cell State", sort="-x"),
            color=alt.Color("modality:N", scale=alt.Scale(range=["#7570b3", "#d95f02"])),
            tooltip=["cell_state:N", "modality:N", "beta:Q", "ci_lower:Q", "ci_upper:Q"],
        )
    )

    error_bars = (
        alt.Chart(combined_long)
        .mark_errorbar()
        .encode(
            x=alt.X("ci_lower:Q", scale=x_scale),
            x2="ci_upper:Q",
            y=alt.Y("cell_state:N", sort="-x"),
            color=alt.Color("modality:N"),
        )
    )

    rule_zero = alt.Chart(pl.DataFrame({"x": [0]})).mark_rule(strokeDash=[3, 3], color="#999999").encode(x="x:Q")

    return (error_bars + points + rule_zero).properties(
        title="Forest Plot: Cross-Cohort Effect Sizes per Cell State (Capped at ±100)",
        width=550,
        height=max(300, norm_details.height * 22),
    )


def plot_insilico_recovery(insilico_df: pl.DataFrame) -> alt.Chart:
    """Fig 4: In silico deconvolution recovery fidelity per cell state."""
    bar = (
        alt.Chart(insilico_df)
        .mark_bar()
        .encode(
            x=alt.X("pearson_r_recovery:Q", title="Deconvolution Recovery Pearson r", scale=alt.Scale(domain=[0, 1])),
            y=alt.Y("cell_state:N", title="Cell State", sort="-x"),
            color=alt.condition(
                alt.datum.pearson_r_recovery >= 0.5,
                alt.value("#1b9e77"),  # Identifiable
                alt.value("#e7298a"),  # Leakage prone
            ),
            tooltip=["cell_state:N", "pearson_r_recovery:Q", "rmse_recovery:Q", "is_identifiable:N"],
        )
    )

    threshold_rule = (
        alt.Chart(pl.DataFrame({"threshold": [0.5]}))
        .mark_rule(color="#d95f02", strokeDash=[4, 4], size=2)
        .encode(x="threshold:Q")
    )

    return (bar + threshold_rule).properties(
        title="In Silico Pseudo-bulk Recovery Identifiability (Threshold r >= 0.5)",
        width=480,
        height=max(250, insilico_df.height * 20),
    )


def run_pipeline(config: PlotConfig) -> Result[None, str]:
    """Execute Step 7 visualization generation."""
    for p in [config.metrics_summary_path, config.state_details_path, config.insilico_path]:
        if not p.exists():
            return Failure(f"Data file does not exist: {p}")

    summary_df = pl.read_parquet(config.metrics_summary_path)
    details_df = pl.read_parquet(config.state_details_path)
    insilico_df = pl.read_parquet(config.insilico_path)

    config.results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fig 1: Concordance Scatter
    fig1 = plot_concordance_scatter(details_df, summary_df)
    match save_altair_svg(fig1, config.results_dir / "fig1_concordance_scatter.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # 2. Fig 2: Normalization Impact
    fig2 = plot_purity_adjustment_impact(summary_df)
    match save_altair_svg(fig2, config.results_dir / "fig2_purity_adjustment_shift.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # 3. Fig 3: Forest Plot
    fig3 = plot_forest_effect_sizes(details_df)
    match save_altair_svg(fig3, config.results_dir / "fig3_forest_effect_sizes.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # 4. Fig 4: In Silico Calibration
    fig4 = plot_insilico_recovery(insilico_df)
    match save_altair_svg(fig4, config.results_dir / "fig4_insilico_calibration.svg"):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] All publication vector SVGs generated in {config.results_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
