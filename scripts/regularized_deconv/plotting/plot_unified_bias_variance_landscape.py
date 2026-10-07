#!/usr/bin/env python3
"""
Publication Figures: Unified Bias-Variance Landscape and Algorithm Dominance Map.

Visualizes the complete 3-factorial synthetic deconvolution experiment:
- Figure 1: 2D Faceted Bias-Variance Decomposition Grid (MSE = Bias² + Variance)
            across collinearity levels r and count parameter n_total for all 7 tools.
- Figure 2: Algorithm Dominance Heatmap & Decision Matrix:
            maps the winning deconvolution method (lowest MSE) across count depth and collinearity.

Strict functional Python: Polars, Altair, plotting_utils, and vl-convert.
"""

from __future__ import annotations

from pathlib import Path
import sys
import altair as alt  # type: ignore
import polars as pl

from plotting_utils import (
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    DECONV_COLOR_MAP,
    DECONV_METHOD_ORDER,
    export_altair_figure,
)


def load_unified_summary(
    data_dir: Path = Path("output/concordance"),
) -> pl.DataFrame:
    """Load unified 3-factorial summary from Parquet."""
    summary_path = data_dir / "unified_synthetic_bias_variance_summary.parquet"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Summary file not found at {summary_path}. "
            "Please run scripts/regularized_deconv/benchmark_unified_synthetic_experiment.py first."
        )
    return pl.read_parquet(summary_path)


def build_bias_variance_landscape_figure(df: pl.DataFrame, architecture: str = "skewed") -> alt.Chart:
    """
    Build a multi-panel stacked bar chart decomposing total MSE into Bias² and Variance
    across Collinearity r and Count Parameter n_total for all 7 tools under specified architecture.
    """
    arch_df = df.filter(pl.col("architecture") == architecture)
    present_methods = [m for m in DECONV_METHOD_ORDER if m in arch_df["method"].unique()]

    # Reshape into long format for stacked bars
    long_rows = []
    for row in arch_df.iter_rows(named=True):
        m = str(row["method"])
        if m in present_methods:
            long_rows.append({
                "target_r": f"r = {float(row['target_r']):.2f}",
                "count_n": f"n = {int(row['count_parameter_n']):,}",
                "architecture": str(row["architecture"]).capitalize(),
                "method": m,
                "component": "Variance",
                "value": float(row["variance"]),
            })
            long_rows.append({
                "target_r": f"r = {float(row['target_r']):.2f}",
                "count_n": f"n = {int(row['count_parameter_n']):,}",
                "architecture": str(row["architecture"]).capitalize(),
                "method": m,
                "component": "Squared Bias",
                "value": float(row["squared_bias"]),
            })

    long_df = pl.DataFrame(long_rows).to_pandas()

    comp_palette = alt.Scale(
        domain=["Variance", "Squared Bias"],
        range=["#7E22CE", "#EAB308"],  # Purple for variance, gold for squared bias
    )

    chart = (
        alt.Chart(long_df)
        .mark_bar(opacity=0.88)
        .encode(
            x=alt.X("method:N", sort=present_methods, title="Deconvolution Method", axis=alt.Axis(labelAngle=-35)),
            y=alt.Y("value:Q", title="MSE Decomposition (Bias² + Var)", stack="zero"),
            color=alt.Color("component:N", scale=comp_palette, title="Error Component", legend=alt.Legend(orient="top")),
            column=alt.Column("target_r:N", title="Collinearity Level r", header=alt.Header(titleFont="Arial", titleFontSize=11, titleFontWeight="bold")),
            row=alt.Row("count_n:N", sort=["n = 2,500", "n = 10,000", "n = 80,000"], title="Multinomial Count Parameter n_total", header=alt.Header(titleFont="Arial", titleFontSize=11, titleFontWeight="bold")),
        )
        .properties(width=160, height=130)
        .configure_axis(
            labelFont="Arial",
            titleFont="Arial",
            labelFontSize=8.5,
            titleFontSize=9.5,
            gridColor=COLOR_LIGHT_GREY,
            domainColor=COLOR_HAIRLINE,
            tickColor=COLOR_HAIRLINE,
        )
        .configure_title(
            font="Arial",
            fontSize=12,
            fontWeight="bold",
            anchor="start",
        )
        .configure_view(
            stroke=None,
        )
    )

    return chart


def build_dominance_matrix_figure(df: pl.DataFrame) -> alt.Chart:
    """
    Build a decision heatmap showing which deconvolution method achieves the minimum MSE
    for every combination of Collinearity r and Count Parameter n_total.
    """
    # Group by (target_r, count_parameter_n, architecture) and find argmin(state_mse)
    min_df = (
        df.sort("state_mse")
        .group_by(["target_r", "count_parameter_n", "architecture"])
        .first()
        .to_pandas()
    )

    present_methods = [m for m in DECONV_METHOD_ORDER if m in min_df["method"].unique()]
    color_scale = alt.Scale(
        domain=present_methods,
        range=[DECONV_COLOR_MAP[m] for m in present_methods],
    )

    min_df["count_label"] = min_df["count_parameter_n"].apply(lambda n: f"{n:,} counts")
    min_df["r_label"] = min_df["target_r"].apply(lambda r: f"r = {r:.2f}")

    chart = (
        alt.Chart(min_df)
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.5)
        .encode(
            x=alt.X("r_label:N", title="Collinearity Spectrum r", axis=alt.Axis(labelAngle=0)),
            y=alt.Y("count_label:N", sort=["2,500 counts", "10,000 counts", "80,000 counts"], title="Count Parameter n_total"),
            color=alt.Color("method:N", scale=color_scale, title="Lowest MSE Method (Winner)"),
            facet=alt.Facet("architecture:N", title="Biological State Architecture", header=alt.Header(titleFont="Arial", titleFontSize=11.5, titleFontWeight="bold")),
        )
        .properties(width=180, height=140)
        .configure_axis(
            labelFont="Arial",
            titleFont="Arial",
            labelFontSize=9.5,
            titleFontSize=10.5,
            gridColor=COLOR_LIGHT_GREY,
            domainColor=COLOR_HAIRLINE,
            tickColor=COLOR_HAIRLINE,
        )
        .configure_title(
            font="Arial",
            fontSize=12,
            fontWeight="bold",
            anchor="start",
        )
        .configure_legend(
            titleFont="Arial",
            labelFont="Arial",
            titleFontSize=10,
            labelFontSize=9,
            orient="bottom",
        )
        .configure_view(
            stroke=None,
        )
    )

    return chart


def main() -> None:
    data_dir = Path("output/concordance")
    fig_dir = Path("article/figures/deconvolution")
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("Loading unified benchmark summary...")
    df = load_unified_summary(data_dir)

    print("Building Bias-Variance Landscape Figure...")
    landscape_chart = build_bias_variance_landscape_figure(df)
    export_altair_figure(landscape_chart, fig_dir / "unified_bias_variance_landscape", scale=2.5)

    print("Building Algorithm Dominance Matrix Figure...")
    dominance_chart = build_dominance_matrix_figure(df)
    export_altair_figure(dominance_chart, fig_dir / "algorithm_dominance_matrix", scale=2.5)

    print("  [OK] Exported unified benchmark landscape figures successfully.")


if __name__ == "__main__":
    main()
