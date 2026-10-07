#!/usr/bin/env python3
"""
Publication Figure: Estimator Bias-Variance Decomposition & True Dropout Diagnostics (All 7 Methods).

Generates a publication-grade two-panel figure:
- Panel A1: Estimator Variance vs Squared Bias across B=30 Replicate Sequencing Draws (r = 0.99).
- Panel A2: Estimator Variance across Replicates on log10 scale (Bar Chart).
- Panel B1: True State & Lineage Dropout Stress Test (False Positive Phantom Detection Rate % on true zeros).
- Panel B2: Active States RMSE.

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


def load_benchmark_data(
    data_dir: Path = Path("output/concordance"),
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Load variance summary and dropout summary from Parquet."""
    var_path = data_dir / "synthetic_estimator_variance_summary.parquet"
    drop_path = data_dir / "synthetic_dropout_stress_test_summary.parquet"

    if not (var_path.exists() and drop_path.exists()):
        raise FileNotFoundError(
            f"Benchmark files not found in {data_dir}. Run benchmark_variance_and_dropouts.py first."
        )

    var_df = pl.read_parquet(var_path)
    drop_df = pl.read_parquet(drop_path)
    return var_df, drop_df


def build_figure(var_df: pl.DataFrame, drop_df: pl.DataFrame) -> alt.VConcatChart:
    """Construct multi-panel diagnostic figure using Altair."""
    present_methods = [m for m in DECONV_METHOD_ORDER if m in var_df["method"].unique()]
    color_scale = alt.Scale(
        domain=present_methods,
        range=[DECONV_COLOR_MAP[m] for m in present_methods],
    )

    # --------------------------------------------------------------------------
    # Panel A1: Stacked Bar of Bias^2 and Variance (% Composition)
    # --------------------------------------------------------------------------
    decomp_rows = []
    for row in var_df.to_dicts():
        m = str(row["method"])
        if m in present_methods:
            decomp_rows.append({"method": m, "component": "Variance", "value": float(row["mean_variance"])})
            decomp_rows.append({"method": m, "component": "Squared Bias", "value": float(row["mean_bias_squared"])})
    decomp_df = pl.DataFrame(decomp_rows).to_pandas()

    comp_palette = alt.Scale(
        domain=["Variance", "Squared Bias"],
        range=["#7E22CE", "#EAB308"],  # Purple for variance, gold for bias
    )

    p_a1 = (
        alt.Chart(decomp_df)
        .mark_bar(opacity=0.88)
        .encode(
            x=alt.X("method:N", sort=present_methods, title="Deconvolution Method", axis=alt.Axis(labelAngle=-25)),
            y=alt.Y("value:Q", title="MSE Error Decomposition", stack="zero"),
            color=alt.Color("component:N", scale=comp_palette, title="MSE Component", legend=alt.Legend(orient="top")),
        )
        .properties(
            width=340,
            height=240,
            title=alt.TitleParams(
                text="A1: Estimator Error Decomposition (r = 0.99, B = 30 Replicates)",
                subtitle="NNLS/CIBERSORT are >89% variance; BayesPrism/InstaPrism are >99% centroid bias",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
                subtitleFontWeight="bold",
            ),
        )
    )

    # --------------------------------------------------------------------------
    # Panel A2: Estimator Variance Bar Chart on log10 Scale
    # --------------------------------------------------------------------------
    var_plot_df = var_df.with_columns(base=pl.lit(1e-6)).to_pandas()
    p_a2 = (
        alt.Chart(var_plot_df)
        .mark_bar(opacity=0.88)
        .encode(
            x=alt.X("method:N", sort=present_methods, title="Deconvolution Method", axis=alt.Axis(labelAngle=-25)),
            y=alt.Y(
                "mean_variance:Q",
                scale=alt.Scale(type="log", domain=[1e-6, 3e-3]),
                title="Estimator Variance Var(θ̂) [log10]",
                axis=alt.Axis(
                    values=[1e-6, 1e-5, 1e-4, 1e-3],
                    labelExpr="datum.value == 1e-6 ? '10⁻⁶' : datum.value == 1e-5 ? '10⁻⁵' : datum.value == 1e-4 ? '10⁻⁴' : datum.value == 1e-3 ? '10⁻³' : datum.label",
                ),
            ),
            y2=alt.Y2("base:Q"),
            color=alt.Color("method:N", scale=color_scale, legend=None),
        )
        .properties(
            width=340,
            height=240,
            title=alt.TitleParams(
                text="A2: Estimator Variance across Replicates (B = 30 Technical Draws)",
                subtitle="InstaPrism (3.0e-6) & BayesPrism (2.6e-5) suppress variance up to 500x",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
                subtitleFontWeight="bold",
            ),
        )
    )

    # --------------------------------------------------------------------------
    # Panel B1: False Positive Phantom Detection on True Zeros
    # --------------------------------------------------------------------------
    phantom_rows = []
    for row in drop_df.to_dicts():
        m = str(row["method"])
        if m in present_methods:
            phantom_rows.append({"method": m, "zero_type": "Absent Sibling State", "pct": float(row["phantom_state_pct"])})
            phantom_rows.append({"method": m, "zero_type": "Absent Lineage (3 States)", "pct": float(row["phantom_lineage_pct"])})
    phantom_df = pl.DataFrame(phantom_rows).to_pandas()

    phantom_palette = alt.Scale(
        domain=["Absent Sibling State", "Absent Lineage (3 States)"],
        range=["#0072B2", "#56B4E9"],
    )

    p_b1 = (
        alt.Chart(phantom_df)
        .mark_bar(opacity=0.88)
        .encode(
            x=alt.X("method:N", sort=present_methods, title="Deconvolution Method", axis=alt.Axis(labelAngle=-25)),
            y=alt.Y("pct:Q", title="False Detection Mass (% on True Zeros)", stack="zero"),
            color=alt.Color("zero_type:N", scale=phantom_palette, title="Dropout Context", legend=alt.Legend(orient="top")),
        )
        .properties(
            width=340,
            height=240,
            title=alt.TitleParams(
                text="B1: True Biological Absence: False Positive Phantom Mass",
                subtitle="Dirichlet priors leak 8.3% into missing sibling states; NNLS preserves zeros (<2.1%)",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
                subtitleFontWeight="bold",
            ),
        )
    )

    # --------------------------------------------------------------------------
    # Panel B2: Active States RMSE
    # --------------------------------------------------------------------------
    p_b2 = (
        alt.Chart(drop_df.to_pandas())
        .mark_bar(opacity=0.88)
        .encode(
            x=alt.X("method:N", sort=present_methods, title="Deconvolution Method", axis=alt.Axis(labelAngle=-25)),
            y=alt.Y("active_states_rmse:Q", title="Active States RMSE (Excluding Zeros)"),
            color=alt.Color("method:N", scale=color_scale, legend=None),
        )
        .properties(
            width=340,
            height=240,
            title=alt.TitleParams(
                text="B2: Reconstruction Fidelity on Non-Zero Active Subsets",
                subtitle="Sparsity-compatible models retain lower RMSE when true subsets drop out",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
                subtitleFontWeight="bold",
            ),
        )
    )

    # Assemble 2x2 grid with independent color scales
    row_a = alt.hconcat(p_a1, p_a2, spacing=24).resolve_scale(color="independent")
    row_b = alt.hconcat(p_b1, p_b2, spacing=24).resolve_scale(color="independent")
    chart = alt.vconcat(row_a, row_b, spacing=28).resolve_scale(color="independent").configure_axis(
        labelFont="Arial",
        titleFont="Arial",
        labelFontSize=9.5,
        titleFontSize=10.5,
        gridColor=COLOR_LIGHT_GREY,
        domainColor=COLOR_HAIRLINE,
        tickColor=COLOR_HAIRLINE,
    ).configure_title(
        font="Arial",
        fontSize=11.5,
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
    base_out = fig_dir / "variance_and_dropout_diagnostics"

    print("Loading benchmark data...")
    var_df, drop_df = load_benchmark_data(data_dir)

    print("Building diagnostic figure...")
    chart = build_figure(var_df, drop_df)

    print("Exporting SVG & 300 DPI PNG via plotting_utils...")
    svg_out, png_out = export_altair_figure(chart, base_out, scale=2.5)
    print(f"  [OK] Exported SVG: {svg_out}")
    print(f"  [OK] Exported PNG: {png_out}")


if __name__ == "__main__":
    main()
