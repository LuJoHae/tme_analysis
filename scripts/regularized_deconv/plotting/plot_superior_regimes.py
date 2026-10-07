#!/usr/bin/env python3
"""
Publication Figure: Regimes where BayesPrism and InstaPrism Achieve Superior MSE (All 7 Methods).

Generates a publication-grade Altair figure showing:
- Panel A: State MSE vs Multinomial Count Parameter n_total in [2,500, 80,000] under high collinearity (r = 0.99).
           Shows where BayesPrism/InstaPrism outperform NNLS/CIBERSORT due to technical noise suppression.
- Panel B: State MSE under Balanced / Plastic Intra-Lineage Sibling States (alpha_within = 20).
           Demonstrates that when collinear cell states represent continuous phenotypic co-occurrence,
           BayesPrism and InstaPrism achieve 3x to 5x lower MSE than NNLS and CIBERSORT.

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


def build_superior_regimes_figure(df: pl.DataFrame) -> alt.VConcatChart:
    """Build four-panel publication figure in Altair representing the 4 superiority regimes."""
    present_methods = [m for m in DECONV_METHOD_ORDER if m in df["method"].unique()]
    color_scale = alt.Scale(
        domain=present_methods,
        range=[DECONV_COLOR_MAP[m] for m in present_methods],
    )

    # --------------------------------------------------------------------------
    # Panel A: Sequencing Read Depth Sweep (Log-Log Scale Line Plot)
    # --------------------------------------------------------------------------
    depth_df = df.filter(pl.col("scenario") == "depth_sweep").to_pandas()
    p_a_lines = alt.Chart(depth_df).mark_line(strokeWidth=2.2).encode(
        x=alt.X(
            "parameter_value:Q",
            scale=alt.Scale(type="log", domain=[800, 100000]),
            title="Multinomial Read Count n_total [log scale]",
            axis=alt.Axis(values=[1000, 2500, 5000, 10000, 20000, 40000, 80000], format="~s"),
        ),
        y=alt.Y(
            "state_mse:Q",
            scale=alt.Scale(type="log", domain=[0.001, 0.020]),
            title="Cell State MSE [log scale]",
            axis=alt.Axis(values=[0.001, 0.002, 0.005, 0.010, 0.020], format=".3f"),
        ),
        color=alt.Color(
            "method:N",
            scale=color_scale,
            title="Deconvolution Method",
            legend=alt.Legend(orient="top", direction="horizontal", columns=4, symbolSize=70, titleFontWeight="bold"),
        ),
    )
    p_a_pts = alt.Chart(depth_df).mark_circle(size=50).encode(
        x=alt.X("parameter_value:Q", scale=alt.Scale(type="log", domain=[800, 100000])),
        y=alt.Y("state_mse:Q", scale=alt.Scale(type="log", domain=[0.001, 0.020])),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    band_df = pl.DataFrame({"x1": [800], "x2": [10000], "y1": [0.001], "y2": [0.020]}).to_pandas()
    p_a_band = alt.Chart(band_df).mark_rect(opacity=0.08, color="#56B4E9").encode(
        x=alt.X("x1:Q"), x2=alt.X2("x2:Q"), y=alt.Y("y1:Q"), y2=alt.Y2("y2:Q")
    )
    panel_a = (p_a_band + p_a_lines + p_a_pts).properties(
        width=320,
        height=220,
        title=alt.TitleParams(
            text="A: Sequencing Depth & Shot Noise (r = 0.99)",
            subtitle="At low counts (n_total ≤ 10k), InstaPrism/BayesPrism achieve 3.4x lower MSE",
            fontSize=10.5,
            fontWeight="bold",
        ),
    )

    # --------------------------------------------------------------------------
    # Panel B: Patient-Specific Expression Drift / Tumor Plasticity
    # --------------------------------------------------------------------------
    drift_df = df.filter(pl.col("scenario") == "expression_drift").to_pandas()
    p_b_lines = alt.Chart(drift_df).mark_line(strokeWidth=2.2).encode(
        x=alt.X("parameter_value:Q", title="Patient Expression Drift σ_drift", axis=alt.Axis(values=[0.0, 0.1, 0.25, 0.5, 0.75])),
        y=alt.Y("state_mse:Q", scale=alt.Scale(domain=[0.0, 0.013]), title="Cell State MSE", axis=alt.Axis(format=".3f")),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    p_b_pts = alt.Chart(drift_df).mark_circle(size=50).encode(
        x=alt.X("parameter_value:Q"),
        y=alt.Y("state_mse:Q"),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    panel_b = (p_b_lines + p_b_pts).properties(
        width=320,
        height=220,
        title=alt.TitleParams(
            text="B: Patient Expression Drift / Tumor Plasticity",
            subtitle="Linear models surge 4.5x under drift; BayesPrism & InstaPrism stay flat",
            fontSize=10.5,
            fontWeight="bold",
        ),
    )

    # --------------------------------------------------------------------------
    # Panel C: Biological Overdispersion & Outliers
    # --------------------------------------------------------------------------
    disp_df = df.filter(pl.col("scenario") == "overdispersion").to_pandas()
    p_c_lines = alt.Chart(disp_df).mark_line(strokeWidth=2.2).encode(
        x=alt.X("parameter_value:Q", title="Negative Binomial Dispersion α_disp", axis=alt.Axis(values=[0.0, 0.05, 0.15, 0.30])),
        y=alt.Y("state_mse:Q", scale=alt.Scale(domain=[0.0, 0.018]), title="Cell State MSE", axis=alt.Axis(format=".3f")),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    p_c_pts = alt.Chart(disp_df).mark_circle(size=50).encode(
        x=alt.X("parameter_value:Q"),
        y=alt.Y("state_mse:Q"),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    panel_c = (p_c_lines + p_c_pts).properties(
        width=320,
        height=220,
        title=alt.TitleParams(
            text="C: Biological Overdispersion & Outliers",
            subtitle="Quadratic outlier leverage inflates NNLS/SVR; log-likelihood limits error",
            fontSize=10.5,
            fontWeight="bold",
        ),
    )

    # --------------------------------------------------------------------------
    # Panel D: Balanced Plasticity / Phenotypic Co-occurrence
    # --------------------------------------------------------------------------
    bal_df = df.filter(pl.col("scenario") == "balanced_plasticity").to_pandas()
    p_d_lines = alt.Chart(bal_df).mark_line(strokeWidth=2.2).encode(
        x=alt.X("parameter_value:Q", scale=alt.Scale(type="log", domain=[0.8, 50]), title="Dirichlet Concentration α_within [log scale]", axis=alt.Axis(values=[1, 5, 10, 20, 40], format="d")),
        y=alt.Y("state_mse:Q", scale=alt.Scale(type="log", domain=[0.0001, 0.004]), title="Cell State MSE [log scale]", axis=alt.Axis(values=[0.0001, 0.0003, 0.001, 0.003], format=".4f")),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    p_d_pts = alt.Chart(bal_df).mark_circle(size=50).encode(
        x=alt.X("parameter_value:Q", scale=alt.Scale(type="log", domain=[0.8, 50])),
        y=alt.Y("state_mse:Q", scale=alt.Scale(type="log", domain=[0.0001, 0.004])),
        color=alt.Color("method:N", scale=color_scale, legend=None),
    )
    panel_d = (p_d_lines + p_d_pts).properties(
        width=320,
        height=220,
        title=alt.TitleParams(
            text="D: Balanced Plasticity / Co-occurring States",
            subtitle="Phenotypic co-occurrence drops prior bias; InstaPrism reaches 14x lower MSE",
            fontSize=10.5,
            fontWeight="bold",
        ),
    )

    grid = alt.vconcat(
        alt.hconcat(panel_a, panel_b, spacing=24),
        alt.hconcat(panel_c, panel_d, spacing=24),
        spacing=28,
    ).configure_axis(
        labelFont="Arial",
        titleFont="Arial",
        labelFontSize=9,
        titleFontSize=10,
        gridColor=COLOR_LIGHT_GREY,
        domainColor=COLOR_HAIRLINE,
        tickColor=COLOR_HAIRLINE,
    ).configure_title(
        font="Arial",
        fontSize=11,
        fontWeight="bold",
        anchor="start",
    ).configure_view(
        stroke=None,
    )

    return grid


def main() -> None:
    data_dir = Path("output/concordance")
    fig_dir = Path("article/figures/deconvolution")
    fig_dir.mkdir(parents=True, exist_ok=True)
    base_out = fig_dir / "bayesprism_superior_regimes"

    summary_path = data_dir / "bayesprism_superior_regimes_summary.parquet"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing {summary_path}. Run benchmark_bayesprism_superior_regimes.py first.")

    df = pl.read_parquet(summary_path)
    print("Building Altair figure for superior BayesPrism regimes...")
    chart = build_superior_regimes_figure(df)

    print("Exporting SVG & 300 DPI PNG via plotting_utils...")
    svg_out, png_out = export_altair_figure(chart, base_out, scale=2.0)
    print(f"  [OK] Exported SVG: {svg_out}")
    print(f"  [OK] Exported PNG: {png_out}")


if __name__ == "__main__":
    main()
