#!/usr/bin/env python3
"""
Publication Figure: Collinearity-Aware Regularized Deconvolution Benchmark (All 7 Methods).

Generates a clean 2x2 multi-panel statistical data figure:
- Panel A: Reference Hessian condition number κ(H) vs. collinearity regimes (log10 scale).
- Panel B: State-level proportion error (MSE, Mean +/- SD) vs. lineage error invariance.
- Panel C: Sibling cross-talk anti-correlation Corr(θ₁, θ₂) across regimes.
- Panel D: Spurious state dropout rate (% states collapsed to 0.0).

Strict functional Python: Polars, Altair, plotting_utils, and vl-convert.
"""

from __future__ import annotations

from pathlib import Path
import sys
import altair as alt  # type: ignore
import polars as pl

from plotting_utils import (
    DECONV_COLOR_MAP,
    DECONV_METHOD_ORDER,
    export_altair_figure,
)


COLOR_LINEAGE = "#64748B"  # Slate grey for invariant / target bounds


def load_data(
    data_dir: Path = Path("output/concordance"),
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Load benchmark summary and sample results from Parquet."""
    summary_path = data_dir / "synthetic_collinearity_benchmark_summary.parquet"
    sample_path = data_dir / "synthetic_collinearity_benchmark_results.parquet"

    if not (summary_path.exists() and sample_path.exists()):
        raise FileNotFoundError(
            f"Benchmark data not found in {data_dir}. "
            "Please run scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py first."
        )

    summary_df = pl.read_parquet(summary_path)
    sample_df = pl.read_parquet(sample_path)
    return summary_df, sample_df


def build_figure(summary_df: pl.DataFrame) -> alt.VConcatChart:
    """Build a publication-grade 2x2 multi-panel statistical figure using Altair."""
    # Prepare data for plotting
    df = summary_df.with_columns([
        (pl.col("state_mse_mean") - pl.col("state_mse_std")).clip(lower_bound=0.0).alias("state_mse_min"),
        (pl.col("state_mse_mean") + pl.col("state_mse_std")).alias("state_mse_max"),
    ]).to_pandas()

    present_methods = [m for m in DECONV_METHOD_ORDER if m in df["method"].unique()]

    color_scale = alt.Scale(
        domain=present_methods,
        range=[DECONV_COLOR_MAP[m] for m in present_methods],
    )

    x_encoding = alt.X(
        "collinearity_r:Q",
        title="Actual Sibling State Correlation r (0.0 → 0.99)",
        scale=alt.Scale(domain=[-0.05, 1.02]),
        axis=alt.Axis(values=[0.0, 0.2, 0.4, 0.6, 0.8, 0.99]),
    )

    # --------------------------------------------------------------------------
    # Panel A: Condition Number κ(H) on log10 scale
    # --------------------------------------------------------------------------
    p1_line = (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=50, filled=True), strokeWidth=2.2)
        .encode(
            x=x_encoding,
            y=alt.Y(
                "condition_number:Q",
                scale=alt.Scale(type="log", domain=[20, 6000]),
                title="Condition Number κ(H) [log10]",
                axis=alt.Axis(values=[20, 50, 100, 200, 500, 1000, 2000, 5000]),
            ),
            color=alt.Color("method:N", scale=color_scale, title="Method"),
        )
    )

    p1_rule = (
        alt.Chart()
        .mark_rule(color=COLOR_LINEAGE, strokeDash=[4, 4], strokeWidth=1.5)
        .encode(y=alt.datum(2000.0))
    )

    p1 = (p1_line + p1_rule).properties(
        width=380,
        height=240,
        title=alt.TitleParams(
            text="A: Reference Hessian Condition Number κ(H)",
            subtitle="Dashed line: Target bound κ_target = 2,000 (Adaptive Deficit)",
            subtitleColor=COLOR_LINEAGE,
            subtitleFontSize=9.5,
            subtitleFontWeight="bold",
        ),
    )

    # --------------------------------------------------------------------------
    # Panel B: State-Level MSE Explosion vs. Lineage Invariance
    # --------------------------------------------------------------------------
    p2_area = (
        alt.Chart(df)
        .mark_area(opacity=0.12)
        .encode(
            x=x_encoding,
            y=alt.Y("state_mse_min:Q", scale=alt.Scale(domain=[0.0, 0.0095])),
            y2=alt.Y2("state_mse_max:Q"),
            color=alt.Color("method:N", scale=color_scale, title="Method"),
        )
    )

    p2_line = (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=50, filled=True), strokeWidth=2.2)
        .encode(
            x=x_encoding,
            y=alt.Y(
                "state_mse_mean:Q",
                title="State Proportion MSE (Mean ± SD)",
                scale=alt.Scale(domain=[0.0, 0.0095]),
                axis=alt.Axis(values=[0.000, 0.002, 0.004, 0.006, 0.008]),
            ),
            color=alt.Color("method:N", scale=color_scale, title="Method"),
        )
    )

    mean_lineage_mse = float(summary_df["lineage_mse_mean"].mean())
    p2_lineage = (
        alt.Chart()
        .mark_rule(color=COLOR_LINEAGE, strokeDash=[3, 3], strokeWidth=1.5)
        .encode(y=alt.datum(mean_lineage_mse))
    )

    p2 = (p2_area + p2_line + p2_lineage).properties(
        width=380,
        height=240,
        title=alt.TitleParams(
            text="B: State-Level MSE Surge vs. Lineage Invariance",
            subtitle=f"Dashed line: Broad Lineage MSE ≈ {mean_lineage_mse:.1e} (Null-space invariant)",
            subtitleColor=COLOR_LINEAGE,
            subtitleFontSize=9.5,
            subtitleFontWeight="bold",
        ),
    )

    # --------------------------------------------------------------------------
    # Panel C: Sibling Cross-Talk Correlation
    # --------------------------------------------------------------------------
    p3_line = (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=50, filled=True), strokeWidth=2.2)
        .encode(
            x=x_encoding,
            y=alt.Y(
                "sibling_corr:Q",
                scale=alt.Scale(domain=[-0.25, 1.05]),
                title="Sibling States Correlation Corr(θ₁, θ₂)",
                axis=alt.Axis(values=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
            ),
            color=alt.Color("method:N", scale=color_scale, title="Method"),
        )
    )

    p3_zero = (
        alt.Chart()
        .mark_rule(color="#94A3B8", strokeDash=[4, 4], strokeWidth=1.0)
        .encode(y=alt.datum(0.0))
    )

    p3 = (p3_line + p3_zero).properties(
        width=380,
        height=240,
        title="C: Sibling Cross-Talk Correlation (Anti-Correlation Suppression)",
    )

    # --------------------------------------------------------------------------
    # Panel D: Spurious State Dropout Rate
    # --------------------------------------------------------------------------
    p4_line = (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(size=50, filled=True), strokeWidth=2.2)
        .encode(
            x=x_encoding,
            y=alt.Y(
                "dropout_rate_pct:Q",
                scale=alt.Scale(domain=[-0.5, 11.5]),
                title="Spurious State Dropouts (% States = 0.0)",
                axis=alt.Axis(values=[0, 2, 4, 6, 8, 10]),
            ),
            color=alt.Color("method:N", scale=color_scale, title="Method"),
        )
    )

    p4 = p4_line.properties(
        width=380,
        height=240,
        title="D: Spurious State Dropout Rate (% States Collapsed to 0.0)",
    )

    # Multi-panel assembly
    chart = ((p1 | p2) & (p3 | p4)).resolve_scale(
        color="shared",
    ).configure_axis(
        labelFont="Arial",
        titleFont="Arial",
        labelFontSize=10,
        titleFontSize=11,
        gridColor="#E2E8F0",
        domainColor="#334155",
        tickColor="#334155",
    ).configure_legend(
        titleFont="Arial",
        labelFont="Arial",
        titleFontSize=10.5,
        labelFontSize=9.5,
        orient="top",
        direction="horizontal",
    ).configure_title(
        font="Arial",
        fontSize=12,
        fontWeight="bold",
        anchor="start",
    ).configure_view(
        stroke=None,
    )

    return chart


def main() -> None:
    data_dir = Path("output/concordance")
    fig_dir = Path("article/figures/deconvolution")
    fig_dir.mkdir(parents=True, exist_ok=True)
    base_out = fig_dir / "collinearity_regularized_synthetic_benchmark"

    print("Loading benchmark data...")
    summary_df, _ = load_data(data_dir)

    print("Building Altair statistical figure (all 7 methods)...")
    chart = build_figure(summary_df)

    print("Exporting SVG & 300 DPI PNG via plotting_utils...")
    svg_out, png_out = export_altair_figure(chart, base_out, scale=2.5)
    print(f"  [OK] Exported SVG: {svg_out}")
    print(f"  [OK] Exported PNG: {png_out}")


if __name__ == "__main__":
    main()
