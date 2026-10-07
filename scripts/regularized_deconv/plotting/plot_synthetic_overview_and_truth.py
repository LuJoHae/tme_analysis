#!/usr/bin/env python3
"""
Diagnostic Publication Figures: Synthetic Test Datasets and Deconvolution Ground-Truth Recovery.

Generates two publication-grade figures:
1. synthetic_dataset_overview.svg & .png:
   - Panel A: Reference signature correlation heatmaps across correlation levels (r = 0.0, 0.6, 0.9, 0.99)
              demonstrating intra-lineage block collinearity vs orthogonal cross-lineage profiles.
   - Panel B: Ground truth Dirichlet proportion distribution across the 12 synthetic states (N=50 samples).
   - Panel C: Bulk mixture read profile (80,000 counts) across representative marker genes.

2. deconvolution_results_vs_ground_truth.svg & .png:
   - Panel A: Scatter grid of Inferred Proportion (θ̂) vs Ground Truth (θ*) across all deconvolution tools
              at severe collinearity (r = 0.99) with y=x identity lines, R², and RMSE annotations.
   - Panel B: Paired stacked bar plots comparing Ground Truth vs Method Predictions across
              representative synthetic samples (N=6 samples) at severe collinearity.

Strict functional Python: Polars, Altair, plotting_utils, vl-convert.
"""

from __future__ import annotations

from pathlib import Path
import sys
import altair as alt  # type: ignore
import numpy as np
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


STATE_NAMES = [
    "L1_State_1", "L1_State_2", "L1_State_3",
    "L2_State_1", "L2_State_2", "L2_State_3",
    "L3_State_1", "L3_State_2", "L3_State_3",
    "L4_State_1", "L4_State_2", "L4_State_3",
]

STATE_COLORS = [
    "#0072B2", "#56B4E9", "#93C5FD",  # Lineage 1
    "#009E73", "#6EE7B7", "#A7F3D0",  # Lineage 2
    "#D55E00", "#E69F00", "#FDE047",  # Lineage 3
    "#7E22CE", "#CC79A7", "#F472B6",  # Lineage 4
]


def load_dataset_data(
    data_dir: Path = Path("output/concordance"),
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Load synthetic true proportions, estimates, and reference matrix."""
    true_path = data_dir / "synthetic_collinearity_true_proportions.parquet"
    est_path = data_dir / "synthetic_collinearity_estimated_proportions.parquet"
    ref_path = data_dir / "synthetic_collinearity_reference_matrices.parquet"

    if not (true_path.exists() and est_path.exists() and ref_path.exists()):
        raise FileNotFoundError(
            f"Benchmark data files not found in {data_dir}. "
            "Please run scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py first."
        )

    return pl.read_parquet(true_path), pl.read_parquet(est_path), pl.read_parquet(ref_path)


# ==============================================================================
# Figure 1: Synthetic Dataset Architecture
# ==============================================================================
def build_dataset_overview_figure(
    true_df: pl.DataFrame,
    ref_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Build multi-panel overview of synthetic reference and mixtures."""
    representative_r = [0.0, 0.6, 0.9, 0.99]
    n_states = len(STATE_NAMES)
    corr_rows: list[dict[str, object]] = []

    for r_val in representative_r:
        sub_ref = ref_df.filter(pl.col("target_r") == r_val)
        gene_cols = [c for c in sub_ref.columns if c.startswith("gene_")]
        if not gene_cols or sub_ref.height != n_states:
            continue
        sig_mat = sub_ref.select(gene_cols).to_numpy()
        corr_mat = np.corrcoef(sig_mat)

        for i in range(n_states):
            for j in range(n_states):
                corr_rows.append({
                    "target_r": f"r = {r_val:.2f}",
                    "state_i": STATE_NAMES[i],
                    "state_j": STATE_NAMES[j],
                    "correlation": float(corr_mat[i, j]),
                })

    corr_df = pl.DataFrame(corr_rows).to_pandas()

    heatmap = (
        alt.Chart(corr_df)
        .mark_rect()
        .encode(
            x=alt.X("state_i:N", sort=STATE_NAMES, title=None, axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("state_j:N", sort=list(reversed(STATE_NAMES)), title=None, axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                "correlation:Q",
                scale=alt.Scale(domain=[-0.3, 0.0, 1.0], range=["#0072B2", "#FFFFFF", "#D55E00"]),
                title="Signature Pearson r",
                legend=alt.Legend(orient="bottom", direction="horizontal", gradientLength=180),
            ),
            column=alt.Column(
                "target_r:N",
                title="A: Reference Signature Correlation Matrices Across Collinearity Levels",
                header=alt.Header(titleFont="Arial", titleFontSize=11.5, titleFontWeight="bold", labelFont="Arial", labelFontSize=10),
            ),
        )
        .properties(width=135, height=135)
    )

    # Panel B: Ground Truth Proportion Boxplot
    prop_long_rows: list[dict[str, object]] = []
    sub_prop = true_df.filter(pl.col("target_r") == 0.99)
    for row in sub_prop.iter_rows(named=True):
        for s_idx, s_name in enumerate(STATE_NAMES):
            prop_long_rows.append({
                "sample_id": row["sample_id"],
                "cell_state": s_name,
                "lineage": s_name.split("_")[0],
                "proportion": float(row[s_name]),
                "color_idx": s_idx,
            })
    prop_df = pl.DataFrame(prop_long_rows).to_pandas()

    state_palette = alt.Scale(domain=STATE_NAMES, range=STATE_COLORS)

    boxplot = (
        alt.Chart(prop_df)
        .mark_boxplot(extent="min-max", size=14)
        .encode(
            x=alt.X("cell_state:N", sort=STATE_NAMES, title="Fine Cell State (12 Synthetic States)", axis=alt.Axis(labelAngle=-45)),
            y=alt.Y("proportion:Q", title="Ground Truth Fraction (θ*)", scale=alt.Scale(domain=[0.0, 0.35])),
            color=alt.Color("cell_state:N", scale=state_palette, legend=None),
        )
        .properties(
            width=360,
            height=200,
            title=alt.TitleParams(
                text="B: Ground Truth Proportion Distribution (N = 50 Samples)",
                subtitle="Non-uniform Dirichlet prior across 12 states (Pure In Silico Generation)",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
            ),
        )
    )

    # Panel C: Representative Bulk Read Count Barplot
    rng = np.random.default_rng(42)
    G = 400
    T = 4
    S = 12
    base_profiles = [rng.gamma(2.0, 1.0, G) for _ in range(T)]
    phi_raw = np.zeros((S, G), dtype=np.float64)
    for t_idx in range(T):
        b_prof = base_profiles[t_idx]
        for s_idx in range(3):
            u_prof = rng.gamma(2.0, 1.0, G)
            phi_raw[t_idx * 3 + s_idx] = np.sqrt(0.8) * b_prof + np.sqrt(0.2) * u_prof
    phi = phi_raw / np.sum(phi_raw, axis=1, keepdims=True)
    sample_true_theta = rng.dirichlet(np.array([1.2, 1.0, 0.8, 1.0, 0.8, 0.6, 1.2, 1.0, 0.8, 0.8, 0.6, 0.4]))
    mixture_counts = rng.multinomial(80_000, sample_true_theta @ phi)

    top_genes_idx = np.argsort(mixture_counts)[-35:]
    mix_rows = [
        {"gene": f"G_{g:03d}", "count": int(mixture_counts[g]), "rank": rank}
        for rank, g in enumerate(top_genes_idx)
    ]
    mix_df = pl.DataFrame(mix_rows).to_pandas()

    barplot = (
        alt.Chart(mix_df)
        .mark_bar(color="#009E73", opacity=0.85)
        .encode(
            x=alt.X("gene:N", sort=alt.EncodingSortField(field="count", order="descending"), title="Representative Marker Genes", axis=alt.Axis(labelAngle=-45)),
            y=alt.Y("count:Q", title="Sequencing Reads (Poisson Counts)"),
        )
        .properties(
            width=360,
            height=200,
            title=alt.TitleParams(
                text="C: Synthetic Bulk Read Distribution (Sample 01, n_total = 80,000 Counts)",
                subtitle="Multinomial sequencing sampling: y_n ~ Multinomial(80,000, θ* · Φ)",
                subtitleColor=COLOR_TEXT_SECONDARY,
                subtitleFontSize=9.5,
            ),
        )
    )

    lower_row = alt.hconcat(boxplot, barplot).resolve_scale(color="independent")
    full_chart = alt.vconcat(heatmap, lower_row).configure_axis(
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

    return full_chart


# ==============================================================================
# Figure 2: Deconvolution Results vs Ground Truth (All Methods)
# ==============================================================================
def build_ground_truth_comparison_figure(
    true_df: pl.DataFrame,
    est_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Build multi-panel comparison of inferred cell fractions vs ground truth across all methods."""
    severe_target_r = 0.99
    sub_true = true_df.filter(pl.col("target_r") == severe_target_r)
    sub_est = est_df.filter(pl.col("target_r") == severe_target_r)

    present_methods = [m for m in DECONV_METHOD_ORDER if m in sub_est["method"].unique()]

    # 1. Build paired scatter dataframe
    scatter_rows: list[dict[str, object]] = []
    true_by_sample: dict[str, dict[str, float]] = {}
    for row in sub_true.iter_rows(named=True):
        sid = str(row["sample_id"])
        true_by_sample[sid] = {s: float(row[s]) for s in STATE_NAMES}

    for row in sub_est.iter_rows(named=True):
        sid = str(row["sample_id"])
        method = str(row["method"])
        if sid in true_by_sample and method in present_methods:
            t_dict = true_by_sample[sid]
            for s_name in STATE_NAMES:
                scatter_rows.append({
                    "sample_id": sid,
                    "method": method,
                    "cell_state": s_name,
                    "lineage": s_name.split("_")[0],
                    "true_proportion": t_dict[s_name],
                    "estimated_proportion": float(row[s_name]),
                })

    sc_df = pl.DataFrame(scatter_rows).to_pandas()

    # Compute R² and RMSE per method
    max_val = 0.45
    r2_rows: list[dict[str, object]] = []
    for m in present_methods:
        m_data = sc_df[sc_df["method"] == m]
        y_true = m_data["true_proportion"].values
        y_est = m_data["estimated_proportion"].values
        ss_res = np.sum((y_true - y_est) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        r2 = max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0
        corr = float(np.corrcoef(y_true, y_est)[0, 1]) if np.std(y_est) > 1e-8 else 0.0
        rmse = float(np.sqrt(np.mean((y_true - y_est) ** 2)))
        r2_rows.append({
            "method": m,
            "label": f"R² = {r2:.3f}\nRMSE = {rmse:.4f}\nr = {corr:.3f}",
            "x": 0.01,
            "y": 0.43,
        })
    annot_df = pl.DataFrame(r2_rows).to_pandas()

    identity_df = pl.DataFrame({"x": [0.0, max_val], "y": [0.0, max_val]}).to_pandas()
    identity_line = (
        alt.Chart(identity_df)
        .mark_line(color="#94A3B8", strokeDash=[4, 4], strokeWidth=1.5, clip=True)
        .encode(x="x:Q", y="y:Q")
    )

    method_charts: list[alt.Chart] = []
    for m_idx, m in enumerate(present_methods):
        m_sc = sc_df[sc_df["method"] == m]
        m_annot = annot_df[annot_df["method"] == m]
        is_first_in_row = (m_idx == 0 or m_idx == 4)
        is_bottom_row = (m_idx >= 4)

        sc_base = alt.Chart(m_sc).encode(
            x=alt.X(
                "true_proportion:Q",
                title="Ground Truth Fraction (θ*)" if is_bottom_row else "",
                scale=alt.Scale(domain=[-0.02, max_val]),
                axis=alt.Axis(labels=True, ticks=True, values=[0.0, 0.1, 0.2, 0.3, 0.4]),
            ),
            y=alt.Y(
                "estimated_proportion:Q",
                title="Inferred Fraction (θ̂)" if is_first_in_row else "",
                scale=alt.Scale(domain=[-0.02, max_val]),
                axis=alt.Axis(labels=True, ticks=True, values=[0.0, 0.1, 0.2, 0.3, 0.4]),
            ),
        )
        sc_pts = sc_base.mark_circle(size=26, opacity=0.65, color=DECONV_COLOR_MAP[m], clip=True)

        annot_lbl = (
            alt.Chart(m_annot)
            .mark_text(align="left", baseline="top", fontSize=8.5, fontWeight="bold", color="#1E293B", lineHeight=11)
            .encode(x="x:Q", y="y:Q", text="label:N")
        )

        m_chart = (sc_pts + identity_line + annot_lbl).properties(
            width=150,
            height=140,
            title=alt.TitleParams(text=m, fontSize=9.5, fontWeight="bold", anchor="middle"),
        )
        method_charts.append(m_chart)

    row1 = alt.hconcat(*method_charts[:4], spacing=16)
    row2 = alt.hconcat(*method_charts[4:], spacing=16)
    scatter_grid = alt.vconcat(row1, row2, spacing=16).properties(
        title=alt.TitleParams(
            text="A: Deconvolution Ground-Truth Recovery Scatter Grid (Severe Collinearity, r = 0.99)",
            subtitle="Inferred Cell Fraction θ̂ vs. Ground Truth θ* across all 7 deconvolution tools (2-Row Grid)",
            fontSize=11.5,
            fontWeight="bold",
            anchor="start",
        )
    )

    # Panel B: Representative Sample Cellular Composition
    rep_samples = [f"Synthetic_Sample_{i:02d}" for i in range(6)]
    comp_rows: list[dict[str, object]] = []

    for sid in rep_samples:
        if sid in true_by_sample:
            t_dict = true_by_sample[sid]
            for s_name in STATE_NAMES:
                comp_rows.append({
                    "sample_id": sid.replace("Synthetic_", ""),
                    "source": "Ground Truth",
                    "cell_state": s_name,
                    "fraction": t_dict[s_name],
                })

    for row in sub_est.iter_rows(named=True):
        sid = str(row["sample_id"])
        method = str(row["method"])
        if sid in rep_samples and method in present_methods:
            for s_name in STATE_NAMES:
                comp_rows.append({
                    "sample_id": sid.replace("Synthetic_", ""),
                    "source": method,
                    "cell_state": s_name,
                    "fraction": float(row[s_name]),
                })

    comp_df = pl.DataFrame(comp_rows).to_pandas()
    source_order = ["Ground Truth"] + present_methods
    state_palette = alt.Scale(domain=STATE_NAMES, range=STATE_COLORS)

    composition_bars = (
        alt.Chart(comp_df)
        .mark_bar(stroke="#FFFFFF", strokeWidth=0.3)
        .encode(
            x=alt.X("source:N", sort=source_order, title="Method / Ground Truth", axis=alt.Axis(labelAngle=-30)),
            y=alt.Y("fraction:Q", stack="zero", title="Cellular Fraction (Sum = 1.0)", scale=alt.Scale(domain=[0.0, 1.0])),
            color=alt.Color("cell_state:N", scale=state_palette, title="Fine Cell State", legend=alt.Legend(orient="bottom", columns=4)),
            column=alt.Column("sample_id:N", title="Representative Synthetic Mixtures (Samples 00 to 05)", header=alt.Header(titleFont="Arial", titleFontSize=10.5, titleFontWeight="bold")),
        )
        .properties(width=110, height=220)
    )

    full_figure = alt.vconcat(scatter_grid, composition_bars, spacing=24).resolve_scale(color="independent").configure_axis(
        labelFont="Arial",
        titleFont="Arial",
        labelFontSize=9,
        titleFontSize=10,
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

    return full_figure


def main() -> None:
    data_dir = Path("output/concordance")
    fig_dir = Path("article/figures/deconvolution")
    fig_dir.mkdir(parents=True, exist_ok=True)

    true_df, est_df, ref_df = load_dataset_data(data_dir)

    print("Building Synthetic Dataset Overview figure...")
    overview_chart = build_dataset_overview_figure(true_df, ref_df)
    export_altair_figure(overview_chart, fig_dir / "synthetic_dataset_overview", scale=2.5)

    print("Building Deconvolution Results vs Ground Truth figure (all methods)...")
    truth_chart = build_ground_truth_comparison_figure(true_df, est_df)
    export_altair_figure(truth_chart, fig_dir / "deconvolution_results_vs_ground_truth", scale=2.5)

    print("  [OK] Exported overview and ground truth figures successfully.")


if __name__ == "__main__":
    main()
