#!/usr/bin/env python3
"""
Step 6: Pure Altair Plotting Engine for Sade-Feldman Deconvolution & Validation Pipeline.
Generates publication-quality resolution-independent vector SVGs (and PNGs) for all analytical steps,
including multi-cohort concordance comparisons, Pre/Post stratification, and Harmony single-cell integration.
Outputs to results/sade_feldman_deconv_validation/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import altair as alt  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# Disable row limit for single-cell UMAP coordinates
alt.data_transformers.disable_max_rows()


class PlotConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    data_dir: Path
    results_dir: Path
    stratum: str = "Melanoma"


def plot_step1_reference_heatmap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot marker gene expression heatmap across reference cell states."""
    tidy_path = data_dir / "reference_phi_tidy.parquet"
    markers_path = data_dir / "reference_marker_genes.parquet"

    if not tidy_path.exists() or not markers_path.exists():
        return Failure(f"Step 1 files missing in {data_dir}")

    try:
        df_markers = pl.read_parquet(markers_path)
        top_genes = df_markers.filter(pl.col("rank") <= 3)["gene"].unique().to_list()

        df_tidy = pl.read_parquet(tidy_path)
        df_sub = df_tidy.filter(pl.col("gene").is_in(top_genes))

        df_plot = (
            df_sub.with_columns(
                (
                    (pl.col("linear_mean") - pl.col("linear_mean").mean().over("gene"))
                    / (pl.col("linear_mean").std().over("gene") + 1e-6)
                ).alias("z_score")
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_plot)
            .mark_rect()
            .encode(
                x=alt.X("gene:N", title="Marker Gene", sort=top_genes),
                y=alt.Y("cluster:N", title="Cell State", sort="ascending"),
                color=alt.Color(
                    "z_score:Q",
                    title="Relative Expression (Z-score)",
                    scale=alt.Scale(scheme="viridis"),
                ),
                tooltip=["cluster", "gene", "linear_mean", "z_score"],
            )
            .properties(
                title="Sade-Feldman Deconvolution Reference: Top Marker Genes",
                width=650,
                height=350,
            )
        )

        out_file = results_dir / "step01_reference_marker_heatmap.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 1 heatmap: {exc}")


def plot_step1b_integration_umap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot Harmony batch correction UMAP comparing platform integration."""
    meta_path = data_dir / "integrated_cell_metadata.parquet"
    if not meta_path.exists():
        return Failure(f"Integrated cell metadata missing: {meta_path}")

    try:
        df_meta = pl.read_parquet(meta_path)
        # Subsample for lightweight rendering if large
        if df_meta.height > 15000:
            df_plot = df_meta.sample(n=15000, seed=42).to_pandas()
        else:
            df_plot = df_meta.to_pandas()

        # Panel 1: Colored by Sequencing Technology / Dataset (Batch)
        panel_batch = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.7)
            .encode(
                x=alt.X("umap_1:Q", title="Integrated UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="Integrated UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color(
                    "sequencing_tech:N",
                    title="Platform",
                    scale=alt.Scale(domain=["Smart-seq2", "10x_Chromium"], range=["#e41a1c", "#377eb8"]),
                ),
                tooltip=["dataset", "sequencing_tech", "integrated_cluster"],
            )
            .properties(title="A. Platform Integration (Harmony Corrected)", width=340, height=320)
        )

        # Panel 2: Colored by Integrated Clusters
        panel_clusters = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.7)
            .encode(
                x=alt.X("umap_1:Q", title="Integrated UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="Integrated UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("integrated_cluster:N", title="Integrated Cluster", scale=alt.Scale(scheme="tableau20")),
                tooltip=["integrated_cluster", "dataset"],
            )
            .properties(title="B. Integrated Biological Cell Clusters", width=340, height=320)
        )

        composite = alt.hconcat(panel_batch, panel_clusters).properties(
            title="Harmony Integration: Sade-Feldman (Smart-seq2) & 10x Single-Cell Atlas"
        )

        out_file = results_dir / "step01b_integrated_umap_batch_correction.svg"
        composite.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 1b integrated UMAP: {exc}")


def plot_step2_fractions_distribution(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot cell state fraction distributions across cohorts."""
    fracs_path = data_dir / "deconv_fractions.parquet"
    if not fracs_path.exists():
        return Failure(f"Fractions file missing: {fracs_path}")

    try:
        df_fracs = pl.read_parquet(fracs_path)
        meta_cols = ["sample_id", "cohort", "cancer_type"]
        cell_states = [c for c in df_fracs.columns if c not in meta_cols]

        df_long = (
            df_fracs.unpivot(
                index=meta_cols,
                on=cell_states,
                variable_name="cell_state",
                value_name="fraction",
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_long)
            .mark_boxplot(extent="min-max", size=15)
            .encode(
                x=alt.X("cell_state:N", title="Cell State", axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("fraction:Q", title="Inferred Cell Fraction", scale=alt.Scale(zero=True)),
                color=alt.Color("cohort:N", title="Cohort"),
            )
            .properties(
                title="Deconvoluted Cell State Fractions across iAtlas Cohorts",
                width=750,
                height=350,
            )
        )

        out_file = results_dir / "step02_deconv_fractions_distribution.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 2 fractions: {exc}")


def plot_step2b_cohort_fractions_grid(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot faceted individual cohort fraction distributions stratified by responder status."""
    sample_fracs_path = data_dir / "sample_fractions_with_response.parquet"
    if not sample_fracs_path.exists():
        return Failure(f"Sample fractions with response file missing: {sample_fracs_path}")

    try:
        df_samples = pl.read_parquet(sample_fracs_path)
        meta_cols = ["sample_id", "cohort", "cancer_type", "response"]
        cell_states = [c for c in df_samples.columns if c not in meta_cols]

        df_long = (
            df_samples.unpivot(
                index=meta_cols,
                on=cell_states,
                variable_name="cell_state",
                value_name="fraction",
            )
            .with_columns(
                pl.when(pl.col("response") == 1)
                .then(pl.lit("Responder"))
                .otherwise(pl.lit("Non-Responder"))
                .alias("response_status")
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_long)
            .mark_boxplot(extent="min-max", size=10)
            .encode(
                x=alt.X("cell_state:N", title="Cell State", axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("fraction:Q", title="Inferred Fraction", scale=alt.Scale(zero=True)),
                color=alt.Color(
                    "response_status:N",
                    title="Response",
                    scale=alt.Scale(domain=["Responder", "Non-Responder"], range=["#e41a1c", "#377eb8"]),
                ),
            )
            .properties(width=280, height=200)
            .facet(facet=alt.Facet("cohort:N", title="iAtlas Cohort"), columns=3)
            .properties(
                title="Deconvoluted Immune Cell State Fractions by Individual iAtlas Cohort and Response Status"
            )
        )

        out_file = results_dir / "step02b_cohort_deconv_fractions_grid.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 2b cohort fractions grid: {exc}")


def plot_step2c_cohort_stacked_compositions(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot horizontal 100% stacked bar chart of average immune cell state composition across cohorts."""
    fracs_path = data_dir / "deconv_fractions.parquet"
    if not fracs_path.exists():
        return Failure(f"Fractions file missing: {fracs_path}")

    try:
        df_fracs = pl.read_parquet(fracs_path)
        meta_cols = ["sample_id", "cohort", "cancer_type"]
        cell_states = [c for c in df_fracs.columns if c not in meta_cols]

        df_mean = (
            df_fracs.group_by(["cohort", "cancer_type"])
            .agg([pl.col(c).mean() for c in cell_states])
            .unpivot(
                index=["cohort", "cancer_type"],
                on=cell_states,
                variable_name="cell_state",
                value_name="mean_fraction",
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_mean)
            .mark_bar()
            .encode(
                y=alt.Y("cohort:N", title="iAtlas Cohort", sort=alt.EncodingSortField(field="cancer_type", order="ascending")),
                x=alt.X("mean_fraction:Q", stack="normalize", title="Relative Immune Proportion", axis=alt.Axis(format="%")),
                color=alt.Color(
                    "cell_state:N",
                    title="Cell State",
                    scale=alt.Scale(scheme="tableau20"),
                ),
                tooltip=["cohort", "cancer_type", "cell_state", "mean_fraction"],
            )
            .properties(
                title="Average Deconvoluted Immune Cell Composition Across 9 iAtlas Cohorts",
                width=650,
                height=320,
            )
        )

        out_file = results_dir / "step02c_cohort_cell_state_stacked_bars.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 2c stacked compositions: {exc}")


def plot_step3_logistic_regression(
    data_dir: Path,
    results_dir: Path,
    stratum: str,
) -> Result[tuple[Path, Path], str]:
    """Plot volcano and forest plot for logistic regression response association."""
    res_path = data_dir / "logistic_regression_results.parquet"
    if not res_path.exists():
        return Failure(f"Logistic results missing: {res_path}")

    try:
        df_res = pl.read_parquet(res_path).filter(pl.col("stratum") == stratum).to_pandas()
        if len(df_res) == 0:
            return Failure(f"No results found for stratum '{stratum}'")

        volcano = (
            alt.Chart(df_res)
            .mark_circle(size=90, opacity=0.85)
            .encode(
                x=alt.X("beta:Q", title="Log Odds Ratio (Effect Size β)"),
                y=alt.Y("log10_pval:Q", title="-log10(P-value)"),
                color=alt.Color(
                    "significant_fdr01:N",
                    title="FDR < 0.1",
                    scale=alt.Scale(domain=[True, False], range=["#e41a1c", "#377eb8"]),
                ),
                tooltip=["cell_state", "beta", "or", "p_value", "fdr", "auc"],
            )
            .properties(
                title=f"Immunotherapy Response Association: {stratum} Cohorts",
                width=450,
                height=380,
            )
        )
        hline = (
            alt.Chart(pd.DataFrame({"y": [-np.log10(0.05)]}))
            .mark_rule(strokeDash=[4, 4], color="gray")
            .encode(y="y:Q")
        )
        volcano_final = volcano + hline
        out_volcano = results_dir / "step03_logistic_regression_volcano.svg"
        volcano_final.save(str(out_volcano))

        points = (
            alt.Chart(df_res)
            .mark_point(filled=True, size=60)
            .encode(
                x=alt.X("or:Q", title="Odds Ratio (95% CI)", scale=alt.Scale(type="log")),
                y=alt.Y("cell_state:N", title="Cell State", sort=alt.EncodingSortField(field="or", order="descending")),
                color=alt.Color("beta:Q", title="Beta", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-3, 3])),
                tooltip=["cell_state", "or", "or_ci_lower", "or_ci_upper", "p_value"],
            )
        )
        error_bars = (
            alt.Chart(df_res)
            .mark_rule()
            .encode(
                x="or_ci_lower:Q",
                x2="or_ci_upper:Q",
                y=alt.Y("cell_state:N", sort=alt.EncodingSortField(field="or", order="descending")),
                color=alt.value("#555555"),
            )
        )
        vline = (
            alt.Chart(pd.DataFrame({"x": [1.0]}))
            .mark_rule(strokeDash=[3, 3], color="black")
            .encode(x="x:Q")
        )
        forest_final = (error_bars + points + vline).properties(
            title=f"Forest Plot: Odds of Immunotherapy Response ({stratum})",
            width=500,
            height=380,
        )
        out_forest = results_dir / "step03_logistic_regression_forest.svg"
        forest_final.save(str(out_forest))

        return Success((out_volcano, out_forest))
    except Exception as exc:
        return Failure(f"Failed to plot step 3 logistic regression: {exc}")


def plot_step4_milopy(data_dir: Path, results_dir: Path) -> Result[tuple[Path, Path], str]:
    """Plot milopy neighborhood volcano and cell-state differential abundance."""
    nhoods_path = data_dir / "milopy_nhoods_results.parquet"
    states_path = data_dir / "milopy_cell_state_da.parquet"

    if not nhoods_path.exists() or not states_path.exists():
        return Failure(f"Milo files missing in {data_dir}")

    try:
        df_nhoods = pl.read_parquet(nhoods_path).with_columns(
            (-pl.col("p_value").log10()).alias("log10_pval"),
            (pl.col("fdr") < 0.1).alias("significant"),
        ).to_pandas()

        volcano = (
            alt.Chart(df_nhoods)
            .mark_circle(size=40, opacity=0.7)
            .encode(
                x=alt.X("logfc:Q", title="Log2 Fold Change (Responder vs Non-Responder)"),
                y=alt.Y("log10_pval:Q", title="-log10(P-value)"),
                color=alt.Color(
                    "significant:N",
                    title="FDR < 0.1",
                    scale=alt.Scale(domain=[True, False], range=["#e41a1c", "#9ecae1"]),
                ),
                tooltip=["nhood_id", "logfc", "p_value", "fdr"],
            )
            .properties(
                title="Milopy Single-Cell Differential Abundance (GSE120575)",
                width=450,
                height=380,
            )
        )
        out_volcano = results_dir / "step04_milopy_nhood_volcano.svg"
        volcano.save(str(out_volcano))

        df_states = pl.read_parquet(states_path).to_pandas()
        bar = (
            alt.Chart(df_states)
            .mark_bar()
            .encode(
                x=alt.X("milo_mean_logfc:Q", title="Mean Milo Log2FC"),
                y=alt.Y("cell_state:N", title="Cell State", sort=alt.EncodingSortField(field="milo_mean_logfc", order="descending")),
                color=alt.Color(
                    "milo_mean_logfc:Q",
                    title="Milo Log2FC",
                    scale=alt.Scale(scheme="redblue", domain=[-1.5, 1.5], reverse=True),
                ),
                tooltip=["cell_state", "n_cells", "milo_mean_logfc", "milo_median_logfc", "milo_wilcoxon_pval"],
            )
            .properties(
                title="Single-Cell Differential Abundance by Reference Cell State",
                width=500,
                height=380,
            )
        )
        out_bar = results_dir / "step04_milopy_cell_state_da.svg"
        bar.save(str(out_bar))

        return Success((out_volcano, out_bar))
    except Exception as exc:
        return Failure(f"Failed to plot step 4 milopy: {exc}")


def plot_step4b_pre_vs_post_milo(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot comparison of Pre-treatment vs Post-treatment vs Combined Milo differential abundance."""
    all_states_path = data_dir / "milopy_cell_state_da_all.parquet"
    if not all_states_path.exists():
        return Failure(f"All conditions Milo states table missing: {all_states_path}")

    try:
        df_all = pl.read_parquet(all_states_path).to_pandas()

        chart = (
            alt.Chart(df_all)
            .mark_bar()
            .encode(
                x=alt.X("condition:N", title="Treatment Status", axis=alt.Axis(labels=True)),
                y=alt.Y("milo_mean_logfc:Q", title="Mean Milo Log2FC (~Response)"),
                color=alt.Color(
                    "condition:N",
                    title="Cohort",
                    scale=alt.Scale(domain=["Pre", "Post", "Combined"], range=["#377eb8", "#e41a1c", "#4daf4a"]),
                ),
                column=alt.Column("cell_state:N", title="Cell State", header=alt.Header(labelAngle=-45)),
                tooltip=["cell_state", "condition", "milo_mean_logfc", "n_cells", "milo_wilcoxon_pval"],
            )
            .properties(width=55, height=260)
        )

        out_file = results_dir / "step04b_milopy_pre_vs_post.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 4b Pre vs Post Milo: {exc}")


def plot_step5_concordance_scatter(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot primary concordance scatter plot between bulk deconv Beta and Milo logFC."""
    conc_path = data_dir / "concordance_metrics.parquet"
    summary_path = data_dir / "concordance_summary.parquet"

    if not conc_path.exists():
        return Failure(f"Concordance file missing: {conc_path}")

    try:
        df_conc = pl.read_parquet(conc_path).to_pandas()
        rho_val, p_val = 0.0, 1.0
        if summary_path.exists():
            df_sum = pl.read_parquet(summary_path)
            rho_val = float(df_sum["spearman_rho"][0])
            p_val = float(df_sum["spearman_pvalue"][0])

        hline = alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(strokeDash=[3, 3], color="gray").encode(y="y:Q")
        vline = alt.Chart(pd.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[3, 3], color="gray").encode(x="x:Q")

        scatter = (
            alt.Chart(df_conc)
            .mark_circle(size=120, opacity=0.9)
            .encode(
                x=alt.X("beta:Q", title="Bulk Deconvolution Logistic Regression Beta (Log Odds Ratio)"),
                y=alt.Y("milo_mean_logfc:Q", title="Single-Cell Milo Mean Log2FC (~Response)"),
                color=alt.Color(
                    "quadrant:N",
                    title="Concordance Category",
                    scale=alt.Scale(
                        domain=[
                            "Concordant Responder",
                            "Concordant Non-Responder",
                            "Discordant (Bulk+, SC-)",
                            "Discordant (Bulk-, SC+)",
                        ],
                        range=["#e41a1c", "#377eb8", "#ff7f00", "#984ea3"],
                    ),
                ),
                tooltip=["cell_state", "beta", "or", "milo_mean_logfc", "quadrant"],
            )
        )

        labels = (
            alt.Chart(df_conc)
            .mark_text(align="left", baseline="middle", dx=7, fontSize=10)
            .encode(
                x="beta:Q",
                y="milo_mean_logfc:Q",
                text="cell_state:N",
                color=alt.value("#333333"),
            )
        )

        trend = (
            alt.Chart(df_conc)
            .transform_regression("beta", "milo_mean_logfc")
            .mark_line(color="black", strokeDash=[4, 4])
            .encode(x="beta:Q", y="milo_mean_logfc:Q")
        )

        chart = (hline + vline + trend + scatter + labels).properties(
            title=f"Cross-Modality Concordance: Bulk Deconv Beta vs. Single-Cell Milo DA (Spearman rho = {rho_val:.2f}, p = {p_val:.3e})",
            width=580,
            height=460,
        )

        out_file = results_dir / "step05_concordance_scatter.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 5 concordance scatter: {exc}")


def plot_step5b_cohort_concordance(data_dir: Path, results_dir: Path) -> Result[tuple[Path, Path], str]:
    """Plot multi-cohort concordance comparison forest/bar chart and individual cohort scatter matrix."""
    summary_path = data_dir / "cohort_level_concordance_summary.parquet"
    metrics_path = data_dir / "concordance_metrics_full.parquet"

    if not summary_path.exists() or not metrics_path.exists():
        return Failure(f"Cohort concordance tables missing in {data_dir}")

    try:
        df_sum = pl.read_parquet(summary_path).to_pandas()

        # Chart 1: Cohort Concordance Summary Forest / Bar Plot
        bar_rho = (
            alt.Chart(df_sum)
            .mark_bar()
            .encode(
                x=alt.X("spearman_rho:Q", title="Spearman Rank Correlation (rho)"),
                y=alt.Y("cohort:N", title="iAtlas Cohort", sort=alt.EncodingSortField(field="spearman_rho", order="descending")),
                color=alt.Color(
                    "cancer_type:N",
                    title="Cancer Type",
                    scale=alt.Scale(domain=["Melanoma", "Bladder", "Pancreatic", "Breast", "Renal Cell"], scheme="category10"),
                ),
                tooltip=["cohort", "cancer_type", "spearman_rho", "spearman_pvalue", "concordance_percentage", "concordant_states_count"],
            )
        )
        vline_0 = alt.Chart(pd.DataFrame({"x": [0.0]})).mark_rule(color="black").encode(x="x:Q")

        summary_chart = (bar_rho + vline_0).properties(
            title="Individual Cohort Concordance: Bulk Deconvolution vs. Sade-Feldman Milo DA",
            width=550,
            height=320,
        )
        out_summary = results_dir / "step05b_cohort_concordance_comparison.svg"
        summary_chart.save(str(out_summary))

        # Chart 2: Multi-panel scatter plot grid across cohorts
        df_metrics = (
            pl.read_parquet(metrics_path)
            .filter(pl.col("stratum").str.starts_with("Cohort_") & (pl.col("condition") == "Combined"))
            .to_pandas()
        )

        hline = alt.Chart(df_metrics).mark_rule(strokeDash=[3, 3], color="gray").encode(y=alt.datum(0.0))
        vline = alt.Chart(df_metrics).mark_rule(strokeDash=[3, 3], color="gray").encode(x=alt.datum(0.0))

        scatters = (
            alt.Chart(df_metrics)
            .mark_circle(size=70, opacity=0.85)
            .encode(
                x=alt.X("beta:Q", title="Bulk Deconv Beta"),
                y=alt.Y("milo_mean_logfc:Q", title="Milo Log2FC"),
                color=alt.Color("is_concordant:N", title="Concordant", scale=alt.Scale(domain=[True, False], range=["#e41a1c", "#377eb8"])),
                tooltip=["cohort", "cell_state", "beta", "milo_mean_logfc", "quadrant"],
            )
        )
        trends = (
            alt.Chart(df_metrics)
            .transform_regression("beta", "milo_mean_logfc", groupby=["cohort"])
            .mark_line(color="black", strokeDash=[3, 3])
            .encode(x="beta:Q", y="milo_mean_logfc:Q")
        )

        grid = (
            alt.layer(hline, vline, trends, scatters, data=df_metrics)
            .properties(width=170, height=170)
            .facet(facet=alt.Facet("cohort:N", title="iAtlas Cohort"), columns=3)
            .properties(title="Individual Cohort Scatter Grid: Bulk Deconvolution vs. Milo Differential Abundance")
        )

        out_grid = results_dir / "step05b_individual_cohort_scatters.svg"
        grid.save(str(out_grid))

        return Success((out_summary, out_grid))
    except Exception as exc:
        return Failure(f"Failed to plot step 5b cohort concordance: {exc}")


def plot_step6_dual_umap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot multi-panel UMAP: Reference Cell States, Single-Cell Milo DA, and Bulk Deconvolution Beta."""
    cells_path = data_dir / "milopy_cell_level_scores.parquet"
    conc_path = data_dir / "concordance_metrics.parquet"

    if not cells_path.exists():
        return Failure(f"Cell level scores missing: {cells_path}")

    try:
        df_cells = pl.read_parquet(cells_path)
        if conc_path.exists():
            df_conc = pl.read_parquet(conc_path).select(["cell_state", "beta"])
            df_cells = df_cells.join(df_conc, on="cell_state", how="left")
        else:
            df_cells = df_cells.with_columns(pl.lit(0.0).alias("beta"))

        df_plot = df_cells.fill_null(0.0).to_pandas()

        panel_a = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.8)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("cell_state:N", title="Reference Cell State", scale=alt.Scale(scheme="tableau20")),
                tooltip=["cell_state", "sample_id", "treatment_status"],
            )
            .properties(title="A. Reference Cell States (Sade-Feldman)", width=320, height=320)
        )

        panel_b = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color(
                    "milo_logfc:Q",
                    title="Milo Log2FC (~Response)",
                    scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0]),
                ),
                tooltip=["cell_state", "milo_logfc"],
            )
            .properties(title="B. Single-Cell Milo DA (~Response)", width=320, height=320)
        )

        panel_c = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color(
                    "beta:Q",
                    title="Bulk Deconv Beta",
                    scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0]),
                ),
                tooltip=["cell_state", "beta"],
            )
            .properties(title="C. Bulk Deconv Beta (Mapped to Manifold)", width=320, height=320)
        )

        composite = alt.hconcat(panel_a, panel_b, panel_c).properties(
            title="Comparison of Single-Cell DA vs. Bulk Deconvolution Response Predictors on Sade-Feldman Manifold"
        )

        out_svg = results_dir / "step06_dual_umap_validation.svg"
        composite.save(str(out_svg))

        try:
            out_png = results_dir / "step06_dual_umap_validation.png"
            composite.save(str(out_png))
        except Exception:
            pass

        return Success(out_svg)
    except Exception as exc:
        return Failure(f"Failed to plot step 6 dual UMAP: {exc}")


def plot_step6b_stratified_umap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot 4-panel UMAP: Reference, Pre-treatment Milo, Post-treatment Milo, and Bulk Deconv."""
    cells_path = data_dir / "milopy_cell_level_scores.parquet"
    conc_path = data_dir / "concordance_metrics.parquet"

    if not cells_path.exists():
        return Failure(f"Cell level scores missing: {cells_path}")

    try:
        df_cells = pl.read_parquet(cells_path)
        if conc_path.exists():
            df_conc = pl.read_parquet(conc_path).select(["cell_state", "beta"])
            df_cells = df_cells.join(df_conc, on="cell_state", how="left")
        else:
            df_cells = df_cells.with_columns(pl.lit(0.0).alias("beta"))

        df_plot = df_cells.fill_null(0.0).to_pandas()

        # Panel A: Reference
        p_ref = (
            alt.Chart(df_plot)
            .mark_circle(size=10, opacity=0.8)
            .encode(
                x=alt.X("umap_1:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("cell_state:N", title="Cell State", scale=alt.Scale(scheme="tableau20")),
            )
            .properties(title="A. Reference Cell States", width=250, height=250)
        )

        # Panel B: Pre-treatment Milo DA
        p_pre = (
            alt.Chart(df_plot.dropna(subset=["milo_logfc_pre"]))
            .mark_circle(size=10, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("milo_logfc_pre:Q", title="Pre Milo Log2FC", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0])),
            )
            .properties(title="B. Baseline Pre-Treatment DA", width=250, height=250)
        )

        # Panel C: Post-treatment Milo DA
        p_post = (
            alt.Chart(df_plot.dropna(subset=["milo_logfc_post"]))
            .mark_circle(size=10, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("milo_logfc_post:Q", title="Post Milo Log2FC", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0])),
            )
            .properties(title="C. On-Treatment Post DA", width=250, height=250)
        )

        # Panel D: Bulk Deconvolution Beta
        p_deconv = (
            alt.Chart(df_plot)
            .mark_circle(size=10, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title=None, axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("beta:Q", title="Bulk Deconv Beta", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0])),
            )
            .properties(title="D. Bulk Deconvolution Beta", width=250, height=250)
        )

        row1 = alt.hconcat(p_ref, p_pre)
        row2 = alt.hconcat(p_post, p_deconv)
        composite = alt.vconcat(row1, row2).properties(
            title="Sade-Feldman Single-Cell Manifold: Pre vs. Post DA and Bulk Deconvolution Response Predictors"
        )

        out_svg = results_dir / "step06b_stratified_umap_validation.svg"
        composite.save(str(out_svg))

        try:
            out_png = results_dir / "step06b_stratified_umap_validation.png"
            composite.save(str(out_png))
        except Exception:
            pass

        return Success(out_svg)
    except Exception as exc:
        return Failure(f"Failed to plot step 6b stratified UMAP: {exc}")


def run_all_plots(config: PlotConfig) -> Result[list[Path], str]:
    """Generate all figures across all steps."""
    config.results_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    print("Generating Step 1 Reference Marker Heatmap...")
    match plot_step1_reference_heatmap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 1b Harmony Integration UMAP...")
    match plot_step1b_integration_umap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 2 Deconvolution Fractions Distribution...")
    match plot_step2_fractions_distribution(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 2b Individual Cohort Deconvolution Fractions Grid...")
    match plot_step2b_cohort_fractions_grid(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 2c Cohort Stacked Cell Composition Bars...")
    match plot_step2c_cohort_stacked_compositions(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 3 Logistic Regression Volcano & Forest Plots...")
    match plot_step3_logistic_regression(config.data_dir, config.results_dir, config.stratum):
        case Success((v_path, f_path)):
            generated.extend([v_path, f_path])
            print(f"  Saved: {v_path.name}, {f_path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 4 Milopy Differential Abundance Plots...")
    match plot_step4_milopy(config.data_dir, config.results_dir):
        case Success((v_path, b_path)):
            generated.extend([v_path, b_path])
            print(f"  Saved: {v_path.name}, {b_path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 4b Pre vs. Post Milo Comparison Plot...")
    match plot_step4b_pre_vs_post_milo(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 5 Primary Concordance Scatter Plot...")
    match plot_step5_concordance_scatter(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 5b Multi-Cohort Concordance Plots...")
    match plot_step5b_cohort_concordance(config.data_dir, config.results_dir):
        case Success((s_path, g_path)):
            generated.extend([s_path, g_path])
            print(f"  Saved: {s_path.name}, {g_path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 6 Multi-Panel UMAP Manifold Visualization...")
    match plot_step6_dual_umap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 6b Stratified Pre vs. Post UMAP Manifold Visualization...")
    match plot_step6b_stratified_umap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    if not generated:
        return Failure("No figures could be generated.")

    return Success(generated)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6: Pure Altair SVG Plotting Engine for Sade-Feldman Pipeline."
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory containing output parquet files",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="results/sade_feldman_deconv_validation",
        help="Directory to save vector SVG figures",
    )
    parser.add_argument(
        "--stratum",
        type=str,
        default="Melanoma",
        help="Stratum to highlight in logistic regression plots",
    )
    args = parser.parse_args()

    config = PlotConfig(
        data_dir=Path(args.data_dir),
        results_dir=Path(args.results_dir),
        stratum=args.stratum,
    )

    match run_all_plots(config):
        case Success(paths):
            print(f"\nStep 6 completed successfully. Generated {len(paths)} SVG figures in: {config.results_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 6 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
