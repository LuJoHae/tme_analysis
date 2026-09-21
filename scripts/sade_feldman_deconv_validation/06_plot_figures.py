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
    """Plot cell state fraction distributions across cohorts, split into subplots per cell state."""
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
            .mark_boxplot(extent="min-max", size=12)
            .encode(
                x=alt.X("cohort:N", title="Cohort", axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("fraction:Q", title="Inferred Fraction", scale=alt.Scale(zero=True)),
                color=alt.Color("cancer_type:N", title="Cancer Type"),
                tooltip=["cohort", "cancer_type", "cell_state", "fraction"],
            )
            .properties(width=220, height=180)
            .facet(facet=alt.Facet("cell_state:N", title="Cell State"), columns=4)
            .resolve_scale(y="independent")
            .properties(
                title="Deconvoluted Cell State Fractions Across iAtlas Cohorts (Stratified by Cell State)"
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
        p_perm = 1.0
        kappa_val = 0.0
        cos_sim = 0.0
        if summary_path.exists():
            df_sum = pl.read_parquet(summary_path)
            rho_val = float(df_sum["spearman_rho"][0])
            p_val = float(df_sum["spearman_pvalue"][0])
            if "spearman_perm_pvalue" in df_sum.columns:
                p_perm = float(df_sum["spearman_perm_pvalue"][0])
            if "cohen_kappa" in df_sum.columns:
                kappa_val = float(df_sum["cohen_kappa"][0])
            if "cosine_similarity" in df_sum.columns:
                cos_sim = float(df_sum["cosine_similarity"][0])

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
            title=f"Cross-Modality Concordance: Bulk Deconv Beta vs. Single-Cell Milo DA (Spearman rho = {rho_val:.2f}, p_perm = {p_perm:.3f}, Cosine Sim = {cos_sim:.2f})",
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

    print("Generating Step 2 Multi-Resolution Fractions (Sade-Feldman Standalone)...")
    match plot_multi_resolution_fractions(config.data_dir, config.results_dir, "sf"):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 2 Multi-Resolution Fractions (Combined Multi-Atlas)...")
    match plot_multi_resolution_fractions(config.data_dir, config.results_dir, "comb"):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Multi-Resolution Predictive Benchmark (AUC & Collinearity)...")
    match plot_resolution_predictive_benchmark(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Multi-Resolution Marker Specificity Heatmaps...")
    match plot_multi_resolution_marker_heatmaps(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Combined Reference UMAP across Resolutions...")
    match plot_combined_reference_umap_resolutions(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Figure 1 — Predictive Capacity Benchmark...")
    match plot_fig1_predictive_capacity(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Figure 2 — Signature Collinearity Benchmark...")
    match plot_fig2_signature_collinearity(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Figure 3 — Cellular Granularity Benchmark...")
    match plot_fig3_cellular_granularity(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Figure 4 — Pareto Tradeoff Optimization...")
    match plot_fig4_pareto_tradeoff(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    if not generated:
        return Failure("No figures could be generated.")

    return Success(generated)


def plot_multi_resolution_fractions(data_dir: Path, results_dir: Path, ref_type: str) -> Result[Path, str]:
    """Plot cell state fraction distributions across 9 iAtlas cohorts for each clustering resolution."""
    try:
        dfs: list[pl.DataFrame] = []
        for res in (0.5, 1.0, 1.5, 2.0):
            p = data_dir / f"deconv_fractions_{ref_type}_res{res}.parquet"
            if not p.exists() and abs(res - 0.5) < 1e-4 and ref_type == "sf":
                p = data_dir / "deconv_fractions.parquet"
            if p.exists():
                df = pl.read_parquet(p)
                meta_cols = {"sample_id", "cohort", "cancer_type"}
                states = [c for c in df.columns if c not in meta_cols]
                df_long = df.melt(
                    id_vars=["cohort", "cancer_type", "sample_id"],
                    value_vars=states,
                    variable_name="cell_state",
                    value_name="fraction",
                ).with_columns(pl.lit(f"res={res}").alias("resolution"))
                dfs.append(df_long)

        if not dfs:
            return Failure(f"No deconv fractions parquets found for ref_type={ref_type}")

        df_all = pl.concat(dfs, how="vertical")
        df_sample = df_all.sample(n=min(5000, df_all.height), seed=42).to_pandas()

        box = (
            alt.Chart(df_sample)
            .mark_boxplot(size=14, opacity=0.85)
            .encode(
                x=alt.X("resolution:N", title="Leiden Resolution"),
                y=alt.Y("fraction:Q", title="Inferred State Fraction", scale=alt.Scale(zero=True)),
                color=alt.Color("resolution:N", title="Resolution", scale=alt.Scale(scheme="viridis")),
            )
        )

        grid = (
            box.properties(width=160, height=140)
            .facet(facet=alt.Facet("cohort:N", title="iAtlas Cohort"), columns=3)
            .properties(
                title=f"Cell State Fraction Distributions across Resolutions ({'Sade-Feldman Standalone' if ref_type == 'sf' else 'Combined Multi-Atlas'})"
            )
        )

        out_name = f"step02_deconv_fractions_by_resolution_{ref_type}.svg"
        out_path = results_dir / out_name
        results_dir.mkdir(parents=True, exist_ok=True)
        grid.save(str(out_path))
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to plot multi-resolution fractions for {ref_type}: {exc}")


def plot_resolution_predictive_benchmark(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot multi-resolution comparative benchmark: Response AUC, Condition Number (Collinearity), and Cluster Count."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure(f"Benchmark summary file not found: {bench_file}")

    try:
        df_bench = pl.read_parquet(bench_file).to_pandas()
        ref_order = sorted(df_bench["reference_type"].unique().tolist())
        
        # Priority mapping for consistent, intuitive reference colors
        preferred_colors = {
            "Sade-Feldman": "#e41a1c",      # Red
            "Combined-Atlas": "#377eb8",    # Blue
            "Jerby-Arnon": "#4daf4a",       # Green
            "Maynard-NSCLC": "#984ea3",     # Purple
            "Ma-Liver": "#ff7f00",          # Orange
            "Yost-BCC": "#a65628",          # Brown
            "PanCancer-Atlas": "#f781bf",   # Pink
        }
        fallback_colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"]
        color_range = [
            preferred_colors.get(ref, fallback_colors[i % len(fallback_colors)])
            for i, ref in enumerate(ref_order)
        ]
        color_scale = alt.Scale(domain=ref_order, range=color_range)

        # Panel A: Melanoma Response Multivariate AUC vs Resolution
        panel_a = (
            alt.Chart(df_bench)
            .mark_line(point=alt.OverlayMarkDef(size=80, filled=True), strokeWidth=2.5)
            .encode(
                x=alt.X("resolution:O", title="Leiden Clustering Resolution"),
                y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.45, 0.85])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale),
                tooltip=["reference_type", "resolution", "n_clusters", "melanoma_multivariate_auc", "pancancer_multivariate_auc", "condition_number"],
            )
            .properties(title="A. Predictive Capacity (Melanoma AUC)", width=240, height=200)
        )
        rule_05 = alt.Chart(pd.DataFrame({"y": [0.5]})).mark_rule(strokeDash=[3, 3], color="gray").encode(y="y:Q")
        panel_a = panel_a + rule_05

        # Panel B: Condition Number (Collinearity) vs Resolution
        panel_b = (
            alt.Chart(df_bench)
            .mark_line(point=alt.OverlayMarkDef(size=80, filled=True), strokeWidth=2.5)
            .encode(
                x=alt.X("resolution:O", title="Leiden Clustering Resolution"),
                y=alt.Y("condition_number:Q", title="Signature Condition Number (kappa)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale),
                tooltip=["reference_type", "resolution", "condition_number"],
            )
            .properties(title="B. Signature Collinearity (Condition Number)", width=240, height=200)
        )

        # Panel C: Number of Resolved Clusters vs Resolution
        panel_c = (
            alt.Chart(df_bench)
            .mark_line(point=alt.OverlayMarkDef(size=80, filled=True), strokeWidth=2.5)
            .encode(
                x=alt.X("resolution:O", title="Leiden Clustering Resolution"),
                y=alt.Y("n_clusters:Q", title="Number of Resolved Clusters"),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale),
                tooltip=["reference_type", "resolution", "n_clusters"],
            )
            .properties(title="C. Cellular Granularity (Clusters)", width=240, height=200)
        )

        chart = (panel_a | panel_b | panel_c).resolve_scale(y="independent").properties(
            title="Multi-Resolution Reference Benchmarking: Predictive Power vs. Matrix Collinearity across Granularities"
        )

        out_path = results_dir / "step05c_resolution_predictive_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))

        # Also attempt static PNG export if vl-convert is available
        png_path = results_dir / "step05c_resolution_predictive_benchmark.png"
        try:
            chart.save(str(png_path), scale_factor=2.0)
            print(f"Exported PNG visualization to: {png_path}")
        except Exception as png_err:
            print(f"Notice: PNG export skipped ({png_err}), SVG saved successfully.")

        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to plot resolution predictive benchmark: {exc}")


def plot_multi_resolution_marker_heatmaps(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot marker gene fold changes across multi-resolution references."""
    try:
        records: list[dict[str, object]] = []
        for res in (0.5, 1.0, 1.5, 2.0):
            p = data_dir / f"reference_marker_genes_res{res}.parquet"
            if not p.exists() and abs(res - 0.5) < 1e-4:
                p = data_dir / "reference_marker_genes.parquet"
            if p.exists():
                df_m = pl.read_parquet(p)
                sub = df_m.filter(pl.col("rank") <= 2)
                for r in sub.iter_rows(named=True):
                    records.append(
                        {
                            "resolution": f"res={res}",
                            "cluster": str(r["cluster"]),
                            "gene": str(r["gene"]),
                            "log2fc": float(r["log2fc"]),
                            "rank": int(r["rank"]),
                        }
                    )

        if not records:
            return Failure(f"No multi-resolution marker gene parquets found in {data_dir}")

        df_plot = pl.DataFrame(records).to_pandas()

        chart = (
            alt.Chart(df_plot)
            .mark_circle(size=60)
            .encode(
                x=alt.X("gene:N", title="Canonical Marker Gene", sort=alt.EncodingSortField(field="log2fc", order="descending")),
                y=alt.Y("cluster:N", title="Cell State Cluster"),
                color=alt.Color("log2fc:Q", title="Log2 FC", scale=alt.Scale(scheme="redblue", domainMid=0)),
                size=alt.Size("log2fc:Q", title="Fold Change", scale=alt.Scale(range=[20, 100])),
                tooltip=["resolution", "cluster", "gene", "log2fc", "rank"],
            )
            .properties(width=500, height=160)
            .facet(facet=alt.Facet("resolution:N", title="Clustering Resolution"), columns=1)
            .properties(title="Multi-Resolution Marker Gene Specificity (Sade-Feldman Reference)")
        )

        out_path = results_dir / "step01c_multi_resolution_marker_heatmaps.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to plot multi-resolution marker heatmaps: {exc}")


def plot_combined_reference_umap_resolutions(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot multi-panel UMAP of Harmony-integrated reference colored by clusters at resolutions 0.5, 1.0, 1.5, 2.0."""
    meta_path = data_dir / "integrated_cell_metadata.parquet"
    if not meta_path.exists():
        return Failure(f"Integrated cell metadata missing: {meta_path}")

    try:
        df_meta = pl.read_parquet(meta_path)
        n_sub = min(3500, df_meta.height)
        df_sub = df_meta.sample(n=n_sub, seed=42)

        long_records: list[dict[str, object]] = []
        for res in (0.5, 1.0, 1.5, 2.0):
            r_col = f"integrated_leiden_{res}"
            if r_col not in df_sub.columns:
                r_col = "integrated_cluster"
            clusters = [str(x) for x in df_sub[r_col].to_list()]
            u1 = df_sub["umap_1"].to_list()
            u2 = df_sub["umap_2"].to_list()
            tech = df_sub["sequencing_tech"].to_list()

            for i in range(n_sub):
                long_records.append(
                    {
                        "umap_1": float(u1[i]),
                        "umap_2": float(u2[i]),
                        "cluster": str(clusters[i]),
                        "resolution": f"Resolution = {res}",
                        "tech": str(tech[i]),
                    }
                )

        df_long = pl.DataFrame(long_records).to_pandas()

        chart = (
            alt.Chart(df_long)
            .mark_circle(size=10, opacity=0.75)
            .encode(
                x=alt.X("umap_1:Q", title="Harmony UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="Harmony UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("cluster:N", title="Integrated Cluster", legend=None, scale=alt.Scale(scheme="tableau20")),
                tooltip=["resolution", "cluster", "tech"],
            )
            .properties(width=240, height=220)
            .facet(facet=alt.Facet("resolution:N", title="Harmony-Integrated Subclustering Resolution"), columns=2)
            .properties(
                title="Harmony Pan-Cancer Reference Manifold across Finer Leiden Clustering Resolutions"
            )
        )

        out_path = results_dir / "step01d_combined_reference_umap_resolutions.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to plot combined reference UMAP across resolutions: {exc}")



# ---------------------------------------------------------------------------
# Benchmark Figure 1 — Predictive Capacity
# ---------------------------------------------------------------------------

def plot_fig1_predictive_capacity(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 1: Comprehensive 18-panel predictive capacity benchmark across references, resolutions,
    dataset integration count (up to n=16), continuous single-cell scaling, and cell subsampling titration ladders.
    """
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        fallback = Path("output/output/sade_feldman_deconv_validation/multi_resolution_benchmark_summary.parquet")
        if fallback.exists():
            bench_file = fallback
        else:
            return Failure(f"Benchmark summary not found: {bench_file}")

    try:
        ref_metadata: dict[str, dict[str, Any]] = {
            "Sade-Feldman": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 16288},
            "Jerby-Arnon": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 7186},
            "Maynard-NSCLC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3000},
            "Ma-Liver": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 5115},
            "Yost-BCC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3500},
            "Combined-Atlas": {"category": "Criteria-Combined", "n_datasets": 16, "n_cells": 41284},
            "Melanoma-Duo": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10686},
            "SS2-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10186},
            "10x-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8615},
            "Cross-Tissue-Pair": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8115},
            "Triple-ICI": {"category": "Criteria-Combined", "n_datasets": 3, "n_cells": 13686},
            "Random-Pair-1": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 12301},
            "Random-Pair-2": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 6500},
            "Random-Pair-3": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 10686},
            "Random-Triplet-1": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15301},
            "Random-Triplet-2": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 11615},
            "Random-Triplet-3": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15801},
            "All-Datasets-Combined": {"category": "Criteria-Combined", "n_datasets": 4, "n_cells": 18801},
            "Random-Quadruplet-1": {"category": "Random-Combined", "n_datasets": 4, "n_cells": 18801},
        }

        df = pl.read_parquet(bench_file)

        if "subsample_fraction" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("subsample_fraction"))

        # Enrich metadata only if missing or zero
        df = df.with_columns(
            pl.when(pl.col("n_datasets").is_not_null() & (pl.col("n_datasets") > 0))
            .then(pl.col("n_datasets"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: ref_metadata.get(r, {}).get("n_datasets", 1),
                    return_dtype=pl.Int64,
                )
            )
            .alias("n_datasets"),
            pl.when(pl.col("n_cells").is_not_null() & (pl.col("n_cells") > 0))
            .then(pl.col("n_cells"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: ref_metadata.get(r, {}).get("n_cells", 0),
                    return_dtype=pl.Int64,
                )
            )
            .alias("n_cells"),
            pl.when(pl.col("reference_category").is_not_null() & (pl.col("reference_category") != "Unknown"))
            .then(pl.col("reference_category"))
            .otherwise(
                pl.col("reference_type").map_elements(
                    lambda r: str(ref_metadata.get(r, {}).get("category", "Single Dataset")),
                    return_dtype=pl.String,
                )
            )
            .alias("reference_category"),
        )

        # Full-cell configurations for resolution, cluster, and dataset-level sweeps
        df_full = df.filter(pl.col("subsample_fraction") >= 0.999)
        df_full_pd = df_full.to_pandas()
        df_all_pd = df.to_pandas()

        # Dataset-level mean summary for Row 3
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

        # Smooth logarithmic saturation fit curves across ALL single cell points (including subsampled)
        valid_cells = df.filter(pl.col("n_cells") > 0)
        x_cells_arr = valid_cells["n_cells"].to_numpy().astype(np.float64)
        x_log = np.log(x_cells_arr)

        m_fit = np.polyfit(x_log, valid_cells["melanoma_multivariate_auc"].to_numpy(), 1)
        p_fit = np.polyfit(x_log, valid_cells["pancancer_multivariate_auc"].to_numpy(), 1)
        c_fit = np.polyfit(x_log, valid_cells["mean_cohort_multivariate_auc"].to_numpy(), 1)

        x_grid = np.geomspace(max(1000.0, float(valid_cells["n_cells"].min())), float(valid_cells["n_cells"].max()), 100)
        df_reg_mel = pd.DataFrame({"n_cells": x_grid, "melanoma_multivariate_auc": m_fit[0] * np.log(x_grid) + m_fit[1]})
        df_reg_pan = pd.DataFrame({"n_cells": x_grid, "pancancer_multivariate_auc": p_fit[0] * np.log(x_grid) + p_fit[1]})
        df_reg_coh = pd.DataFrame({"n_cells": x_grid, "mean_cohort_multivariate_auc": c_fit[0] * np.log(x_grid) + c_fit[1]})

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined", "Cell-Subsample"]
        color_scale = alt.Scale(scheme="tableau20")
        dash_scale = alt.Scale(
            domain=["Single Dataset", "Criteria-Combined", "Random-Combined"],
            range=[[1, 0], [6, 2], [2, 4]],
        )

        base = alt.Chart(df_full_pd)
        base_all = alt.Chart(df_all_pd)
        rule_05 = alt.Chart(pd.DataFrame({"y": [0.5]})).mark_rule(strokeDash=[3, 3], color="gray", opacity=0.7).encode(y="y:Q")

        PANEL_WIDTH = 280
        PANEL_HEIGHT = 200

        # ---------------------------------------------------------------------------
        # ROW 1: AUC by Leiden Resolution
        # ---------------------------------------------------------------------------
        panel_a = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=alt.Legend(columns=2)),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale),
                tooltip=["reference_type", "reference_category", "resolution", "melanoma_multivariate_auc", "n_clusters"],
            )
            .properties(title="A. Melanoma Response AUC by Resolution", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_b = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "resolution", "pancancer_multivariate_auc"],
            )
            .properties(title="B. Pan-Cancer AUC by Resolution", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_c = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "resolution", "mean_cohort_multivariate_auc"],
            )
            .properties(title="C. Mean Cohort AUC by Resolution", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        # ---------------------------------------------------------------------------
        # ROW 2: AUC against Number of Reference Cell States (n_clusters)
        # ---------------------------------------------------------------------------
        panel_d = (
            base.mark_line(strokeWidth=1.8, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("n_clusters:Q", title="Reference Cell States (Clusters)", axis=alt.Axis(tickMinStep=5)),
                y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "n_clusters", "resolution", "melanoma_multivariate_auc"],
            )
            .properties(title="D. Melanoma AUC vs. Number of Cell States", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_e = (
            base.mark_line(strokeWidth=1.8, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("n_clusters:Q", title="Reference Cell States (Clusters)", axis=alt.Axis(tickMinStep=5)),
                y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "n_clusters", "resolution", "pancancer_multivariate_auc"],
            )
            .properties(title="E. Pan-Cancer AUC vs. Number of Cell States", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_f = (
            base.mark_line(strokeWidth=1.8, point=alt.OverlayMarkDef(size=45, filled=True))
            .encode(
                x=alt.X("n_clusters:Q", title="Reference Cell States (Clusters)", axis=alt.Axis(tickMinStep=5)),
                y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "n_clusters", "resolution", "mean_cohort_multivariate_auc"],
            )
            .properties(title="F. Mean Cohort AUC vs. Number of Cell States", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        # ---------------------------------------------------------------------------
        # ROW 3: AUC across Number of Datasets Integrated (n up to 16)
        # ---------------------------------------------------------------------------
        base_ds = alt.Chart(df_ds_summary)

        panel_g = (
            alt.layer(
                base.mark_boxplot(size=26, opacity=0.75, color="#e0e0e0").encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                ),
                base.mark_circle(size=32, opacity=0.65).encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    xOffset=alt.XOffset("resolution:Q", scale=alt.Scale(range=[-10, 10])),
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.55, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "melanoma_multivariate_auc"],
                ),
                base_ds.mark_line(color="#2b5c8f", strokeWidth=2.2, strokeDash=[3, 2], point=alt.OverlayMarkDef(color="#2b5c8f", size=50)).encode(
                    x=alt.X("n_datasets:O"),
                    y=alt.Y("mean_mel_auc:Q"),
                ),
            )
            .properties(title="G. Melanoma AUC vs. Number of Integrated Datasets", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        panel_h = (
            alt.layer(
                base.mark_boxplot(size=26, opacity=0.75, color="#e0e0e0").encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                ),
                base.mark_circle(size=32, opacity=0.65).encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    xOffset=alt.XOffset("resolution:Q", scale=alt.Scale(range=[-10, 10])),
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "pancancer_multivariate_auc"],
                ),
                base_ds.mark_line(color="#2b5c8f", strokeWidth=2.2, strokeDash=[3, 2], point=alt.OverlayMarkDef(color="#2b5c8f", size=50)).encode(
                    x=alt.X("n_datasets:O"),
                    y=alt.Y("mean_pan_auc:Q"),
                ),
            )
            .properties(title="H. Pan-Cancer AUC vs. Number of Integrated Datasets", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        panel_i = (
            alt.layer(
                base.mark_boxplot(size=26, opacity=0.75, color="#e0e0e0").encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                ),
                base.mark_circle(size=32, opacity=0.65).encode(
                    x=alt.X("n_datasets:O", title="Number of Integrated Datasets"),
                    xOffset=alt.XOffset("resolution:Q", scale=alt.Scale(range=[-10, 10])),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_datasets", "resolution", "mean_cohort_multivariate_auc"],
                ),
                base_ds.mark_line(color="#2b5c8f", strokeWidth=2.2, strokeDash=[3, 2], point=alt.OverlayMarkDef(color="#2b5c8f", size=50)).encode(
                    x=alt.X("n_datasets:O"),
                    y=alt.Y("mean_coh_auc:Q"),
                ),
            )
            .properties(title="I. Mean Cohort AUC vs. Number of Integrated Datasets", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        # ---------------------------------------------------------------------------
        # ROW 4: AUC across Continuous Single Cell Scaling (Logarithmic Fit)
        # ---------------------------------------------------------------------------
        chart_reg_mel = alt.Chart(df_reg_mel).mark_line(color="#1f77b4", strokeWidth=2.2, strokeDash=[4, 3]).encode(
            x="n_cells:Q",
            y="melanoma_multivariate_auc:Q",
        )
        panel_j = (
            alt.layer(
                base_all.mark_circle(size=40, opacity=0.7).encode(
                    x=alt.X("n_cells:Q", title="Number of Single Cells Integrated", axis=alt.Axis(format="~s")),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "resolution", "subsample_fraction", "melanoma_multivariate_auc"],
                ),
                chart_reg_mel,
            )
            .properties(title="J. Melanoma AUC vs. Single Cells (Saturation Curve)", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        chart_reg_pan = alt.Chart(df_reg_pan).mark_line(color="#ff7f0e", strokeWidth=2.2, strokeDash=[4, 3]).encode(
            x="n_cells:Q",
            y="pancancer_multivariate_auc:Q",
        )
        panel_k = (
            alt.layer(
                base_all.mark_circle(size=40, opacity=0.7).encode(
                    x=alt.X("n_cells:Q", title="Number of Single Cells Integrated", axis=alt.Axis(format="~s")),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "resolution", "subsample_fraction", "pancancer_multivariate_auc"],
                ),
                chart_reg_pan,
            )
            .properties(title="K. Pan-Cancer AUC vs. Single Cells (Saturation Curve)", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        chart_reg_coh = alt.Chart(df_reg_coh).mark_line(color="#2ca02c", strokeWidth=2.2, strokeDash=[4, 3]).encode(
            x="n_cells:Q",
            y="mean_cohort_multivariate_auc:Q",
        )
        panel_l = (
            alt.layer(
                base_all.mark_circle(size=40, opacity=0.7).encode(
                    x=alt.X("n_cells:Q", title="Number of Single Cells Integrated", axis=alt.Axis(format="~s")),
                    y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "n_cells", "resolution", "subsample_fraction", "mean_cohort_multivariate_auc"],
                ),
                chart_reg_coh,
            )
            .properties(title="L. Mean Cohort AUC vs. Single Cells (Saturation Curve)", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        # ---------------------------------------------------------------------------
        # ROW 5: AUC across Cell Subsampling Fraction (10% to 100% Titration Curves)
        # ---------------------------------------------------------------------------
        df_sub = (
            df.filter(
                pl.col("reference_type").str.contains("Cells")
                | pl.col("reference_type").is_in(["Combined-Atlas", "All-Datasets-Combined", "Random-Triplet-1"])
            )
            .with_columns(
                pl.col("reference_type")
                .map_elements(lambda r: r.split(" (")[0] if " (" in r else r, return_dtype=pl.String)
                .alias("subsample_series")
            )
            .filter(pl.col("resolution") == 0.5)
            .sort(["subsample_series", "subsample_fraction"])
        )
        df_sub_pd = df_sub.to_pandas()
        base_sub = alt.Chart(df_sub_pd)
        series_scale = alt.Scale(
            domain=["Combined-Atlas", "All-Datasets-Combined", "Random-Triplet-1"],
            range=["#d95f02", "#7570b3", "#1b9e77"],
        )

        panel_m = (
            base_sub.mark_line(strokeWidth=2.2, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("subsample_fraction:Q", title="Cell Subsampling Fraction", axis=alt.Axis(format="%", tickMinStep=0.1)),
                y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                color=alt.Color("subsample_series:N", title="Atlas Series", scale=series_scale),
                tooltip=["subsample_series", "subsample_fraction", "n_cells", "melanoma_multivariate_auc"],
            )
            .properties(title="M. Melanoma AUC vs. Cell Subsampling Fraction", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_n = (
            base_sub.mark_line(strokeWidth=2.2, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("subsample_fraction:Q", title="Cell Subsampling Fraction", axis=alt.Axis(format="%", tickMinStep=0.1)),
                y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                color=alt.Color("subsample_series:N", title="Atlas Series", scale=series_scale),
                tooltip=["subsample_series", "subsample_fraction", "n_cells", "pancancer_multivariate_auc"],
            )
            .properties(title="N. Pan-Cancer AUC vs. Cell Subsampling Fraction", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        ) + rule_05

        panel_o = (
            base_sub.mark_line(strokeWidth=2.2, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("subsample_fraction:Q", title="Cell Subsampling Fraction", axis=alt.Axis(format="%", tickMinStep=0.1)),
                y=alt.Y("mean_cohort_multivariate_auc:Q", title="Mean Cohort AUC", scale=alt.Scale(domain=[0.58, 0.78])),
                color=alt.Color("subsample_series:N", title="Atlas Series", scale=series_scale),
                tooltip=["subsample_series", "subsample_fraction", "n_cells", "mean_cohort_multivariate_auc"],
            )
            .properties(title="O. Mean Cohort AUC vs. Cell Subsampling Fraction", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        # ---------------------------------------------------------------------------
        # ROW 6: Strategy Distribution & Marginal Gain
        # ---------------------------------------------------------------------------
        panel_p = (
            alt.layer(
                base.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma AUC", scale=alt.Scale(domain=[0.55, 0.74])),
                    color=alt.Color("reference_category:N", title="Strategy", legend=None),
                ),
                base.mark_circle(size=25, opacity=0.5, xOffset=8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("melanoma_multivariate_auc:Q", scale=alt.Scale(domain=[0.55, 0.74])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "resolution", "melanoma_multivariate_auc"],
                ),
            )
            .properties(title="P. Melanoma AUC Distribution by Strategy", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        panel_q = (
            alt.layer(
                base.mark_boxplot(size=28, outliers=False, opacity=0.8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("pancancer_multivariate_auc:Q", title="Pan-Cancer AUC", scale=alt.Scale(domain=[0.54, 0.68])),
                    color=alt.Color("reference_category:N", title="Strategy", legend=None),
                ),
                base.mark_circle(size=25, opacity=0.5, xOffset=8).encode(
                    x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                    y=alt.Y("pancancer_multivariate_auc:Q", scale=alt.Scale(domain=[0.54, 0.68])),
                    color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                    tooltip=["reference_type", "resolution", "pancancer_multivariate_auc"],
                ),
            )
            .properties(title="Q. Pan-Cancer AUC Distribution by Strategy", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

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
        panel_r = (
            alt.Chart(df_delta)
            .mark_bar()
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("mean_delta_auc:Q", title="Mean ΔAUC vs Previous Resolution"),
                color=alt.condition(
                    alt.datum.mean_delta_auc > 0,
                    alt.value("#4daf4a"),
                    alt.value("#e41a1c"),
                ),
                tooltip=["resolution", "mean_delta_auc"],
            )
            .properties(title="R. Marginal AUC Gain per Resolution Step", width=PANEL_WIDTH, height=PANEL_HEIGHT)
        )

        row1 = (panel_a | panel_b | panel_c).resolve_scale(color="shared", strokeDash="shared")
        row2 = (panel_d | panel_e | panel_f).resolve_scale(color="shared", strokeDash="shared")
        row3 = (panel_g | panel_h | panel_i).resolve_scale(color="shared")
        row4 = (panel_j | panel_k | panel_l).resolve_scale(color="shared")
        row5 = (panel_m | panel_n | panel_o).resolve_scale(color="shared")
        row6 = (panel_p | panel_q | panel_r).resolve_scale(color="independent")

        final_chart = (
            alt.vconcat(row1, row2, row3, row4, row5, row6, spacing=25)
            .properties(
                title="Figure 1 — Predictive Capacity Benchmark: Multi-Resolution Reference & Subsampling Comparison"
            )
            .configure_title(fontSize=16, anchor="start", font="Helvetica")
            .configure_axis(labelFontSize=10, titleFontSize=11, titleFontWeight="bold")
        )

        out_path = results_dir / "fig1_predictive_capacity_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_path))
        png_path = results_dir / "fig1_predictive_capacity_benchmark.png"
        try:
            final_chart.save(str(png_path), scale_factor=2.0)
        except Exception:
            pass
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 1: {exc}")


# ---------------------------------------------------------------------------
# Benchmark Figure 2 — Signature Collinearity
# ---------------------------------------------------------------------------

def plot_fig2_signature_collinearity(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 2: 4-panel signature collinearity analysis across all references and resolutions."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure(f"Benchmark summary not found: {bench_file}")

    try:
        df = pl.read_parquet(bench_file).filter(pl.col("condition_number").is_not_null() & pl.col("condition_number").is_finite())
        df_pd = df.to_pandas()

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]
        color_scale = alt.Scale(scheme="tableau20")

        base = alt.Chart(df_pd)

        # Panel 2A: Log-scale κ trajectory with threshold guidelines
        kappa_thresholds = pd.DataFrame({"y": [30.0, 100.0], "label": ["κ=30 (Moderate)", "κ=100 (High)"]})
        rules_2a = (
            alt.Chart(kappa_thresholds)
            .mark_rule(strokeDash=[4, 3], opacity=0.6)
            .encode(
                y="y:Q",
                color=alt.Color("label:N", scale=alt.Scale(range=["#ff7f00", "#e41a1c"]), title="Threshold"),
            )
        )
        panel_a = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("condition_number:Q", title="Condition Number (κ)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=alt.Legend(columns=2)),
                tooltip=["reference_type", "reference_category", "resolution", "condition_number", "n_clusters"],
            )
            .properties(title="A. Signature Collinearity (κ) Trajectory", width=340, height=240)
        )
        panel_a = (panel_a + rules_2a).resolve_scale(color="independent")

        # Panel 2B: Cluster count vs κ phase-space scatter
        panel_b = (
            base.mark_circle(opacity=0.75)
            .encode(
                x=alt.X("n_clusters:Q", title="Number of Clusters"),
                y=alt.Y("condition_number:Q", title="Condition Number (κ)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                size=alt.Size("resolution:Q", title="Resolution", scale=alt.Scale(range=[30, 180])),
                tooltip=["reference_type", "reference_category", "resolution", "n_clusters", "condition_number"],
            )
            .properties(title="B. Cluster Count vs κ Phase Space", width=280, height=240)
        )
        kappa_line = alt.Chart(pd.DataFrame({"x": [0, 60]})).mark_rule(strokeDash=[4, 4], color="#e41a1c", opacity=0.5).encode(
            y=alt.datum(100)
        )
        panel_b = panel_b + kappa_line

        # Panel 2C: Condition number κ across references at coarse (0.25), medium (1.0), fine (2.0) resolution
        df_bars = (
            df.filter(pl.col("resolution").is_in([0.25, 1.0, 2.0]))
            .with_columns(pl.col("resolution").cast(pl.Utf8).alias("res_label"))
            .sort(["reference_category", "reference_type"])
            .to_pandas()
        )
        panel_c = (
            alt.Chart(df_bars)
            .mark_circle(size=80, opacity=0.9)
            .encode(
                x=alt.X("reference_type:N", title=None, sort=alt.EncodingSortField(field="reference_category"), axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("condition_number:Q", title="Condition Number (κ)", scale=alt.Scale(type="log", domain=[3, 400])),
                color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                column=alt.Column("res_label:N", title="Resolution"),
                tooltip=["reference_type", "reference_category", "res_label", "condition_number"],
            )
            .properties(title="C. κ at Coarse / Medium / Fine Resolution", width=180, height=200)
        )

        # Panel 2D: Distribution violin/strip of κ per strategy
        panel_d = alt.layer(
            base.mark_boxplot(size=35, outliers=False, opacity=0.85).encode(
                x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                y=alt.Y("condition_number:Q", title="κ (log scale)", scale=alt.Scale(type="log")),
                color=alt.Color("reference_category:N", legend=None, sort=cat_order),
            ),
            base.mark_circle(size=25, opacity=0.55, xOffset=8).encode(
                x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                y=alt.Y("condition_number:Q", scale=alt.Scale(type="log")),
                color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                tooltip=["reference_type", "resolution", "condition_number"],
            ),
        ).properties(title="D. κ Distribution by Strategy", width=240, height=220)

        top_row = (panel_a | panel_b).resolve_scale(color="independent")
        bot_row = (panel_d).resolve_scale(color="independent")
        chart = alt.vconcat(
            top_row,
            alt.hconcat(panel_c, panel_d).resolve_scale(color="independent"),
        ).properties(
            title="Figure 2 — Signature Collinearity Analysis: Condition Number across References and Resolutions"
        )

        out_path = results_dir / "fig2_signature_collinearity_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))
        png_path = results_dir / "fig2_signature_collinearity_benchmark.png"
        try:
            chart.save(str(png_path), scale_factor=2.0)
        except Exception:
            pass
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 2: {exc}")


# ---------------------------------------------------------------------------
# Benchmark Figure 3 — Cellular Granularity
# ---------------------------------------------------------------------------

def plot_fig3_cellular_granularity(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 3: 4-panel cellular granularity analysis: cluster counts, signature genes, per-cluster AUC."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure(f"Benchmark summary not found: {bench_file}")

    try:
        df = pl.read_parquet(bench_file)
        df_pd = df.to_pandas()

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]
        color_scale = alt.Scale(scheme="tableau20")
        dash_scale = alt.Scale(domain=cat_order, range=[[1, 0], [6, 2], [2, 4]])

        base = alt.Chart(df_pd)

        # Panel 3A: Cluster count scaling curves
        panel_a = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("n_clusters:Q", title="Number of Clusters"),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=alt.Legend(columns=2)),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale),
                tooltip=["reference_type", "reference_category", "resolution", "n_clusters"],
            )
            .properties(title="A. Cluster Count Scaling", width=300, height=210)
        )

        # Panel 3B: Signature gene burden across resolutions
        df_sig = df.filter(pl.col("n_signature_genes") > 0).to_pandas()
        if df_sig.empty:
            df_sig = df_pd.copy()
        base_sig = alt.Chart(df_sig)
        panel_b = (
            base_sig.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("n_signature_genes:Q", title="Signature Gene Count"),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "resolution", "n_signature_genes"],
            )
            .properties(title="B. Signature Gene Burden", width=300, height=210)
        )

        # Panel 3C: AUC per cluster — efficiency metric (melanoma AUC / n_clusters)
        df_eff = df.with_columns(
            (pl.col("melanoma_multivariate_auc") / pl.col("n_clusters").cast(pl.Float64)).alias("auc_per_cluster")
        ).to_pandas()
        base_eff = alt.Chart(df_eff)
        panel_c = (
            base_eff.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("auc_per_cluster:Q", title="AUC per Cluster"),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=None),
                strokeDash=alt.StrokeDash("reference_category:N", title="Strategy", scale=dash_scale, legend=None),
                tooltip=["reference_type", "reference_category", "resolution", "auc_per_cluster", "n_clusters"],
            )
            .properties(title="C. AUC Efficiency per Cluster", width=300, height=210)
        )

        # Panel 3D: Violin/strip of n_clusters by strategy
        panel_d = alt.layer(
            base.mark_boxplot(size=35, outliers=False, opacity=0.85).encode(
                x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                y=alt.Y("n_clusters:Q", title="Number of Clusters"),
                color=alt.Color("reference_category:N", legend=None, sort=cat_order),
            ),
            base.mark_circle(size=28, opacity=0.55, xOffset=8).encode(
                x=alt.X("reference_category:N", title="Strategy", sort=cat_order),
                y=alt.Y("n_clusters:Q"),
                color=alt.Color("reference_type:N", scale=color_scale, legend=None),
                tooltip=["reference_type", "resolution", "n_clusters"],
            ),
        ).properties(title="D. Cluster Count Distribution by Strategy", width=240, height=210)

        top_row = (panel_a | panel_b | panel_c).resolve_scale(color="shared", strokeDash="shared")
        chart = (top_row & panel_d).resolve_scale(color="independent").properties(
            title="Figure 3 — Cellular Granularity: Cluster Count, Signature Burden, AUC Efficiency"
        )

        out_path = results_dir / "fig3_cellular_granularity_benchmark.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))
        png_path = results_dir / "fig3_cellular_granularity_benchmark.png"
        try:
            chart.save(str(png_path), scale_factor=2.0)
        except Exception:
            pass
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 3: {exc}")


# ---------------------------------------------------------------------------
# Benchmark Figure 4 — Pareto Tradeoff Optimization
# ---------------------------------------------------------------------------

def plot_fig4_pareto_tradeoff(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Figure 4: 3-panel Pareto tradeoff analysis: Goldilocks scatter, heatmap, efficiency trajectory."""
    bench_file = data_dir / "multi_resolution_benchmark_summary.parquet"
    if not bench_file.exists():
        return Failure(f"Benchmark summary not found: {bench_file}")

    try:
        df = pl.read_parquet(bench_file).filter(
            pl.col("condition_number").is_not_null() & pl.col("condition_number").is_finite() & (pl.col("condition_number") > 0)
        )
        df_pareto = df.with_columns(
            pl.col("condition_number").log(base=10.0).alias("log10_kappa"),
            (pl.col("melanoma_multivariate_auc") / pl.col("condition_number").log(base=10.0)).alias("efficiency"),
        )
        df_pd = df_pareto.to_pandas()

        cat_order = ["Single Dataset", "Criteria-Combined", "Random-Combined"]
        color_scale = alt.Scale(scheme="tableau20")

        base = alt.Chart(df_pd)

        # Panel 4A: Goldilocks scatter — log10(κ) vs Melanoma AUC
        # Optimal zone shading: AUC>0.65 & κ<30 → log10(κ)<1.477
        opt_zone = alt.Chart(pd.DataFrame({"x1": [0.4], "x2": [1.477], "y1": [0.65], "y2": [0.735]}))
        opt_rect = opt_zone.mark_rect(opacity=0.10, color="#4daf4a").encode(
            x="x1:Q", x2="x2:Q", y="y1:Q", y2="y2:Q"
        )
        opt_label = alt.Chart(pd.DataFrame({"x": [0.9], "y": [0.725], "label": ["Optimal Zone (AUC>0.65, κ<30)"]})).mark_text(
            color="#2c8e2c", fontSize=11, fontStyle="italic", fontWeight="bold"
        ).encode(x="x:Q", y="y:Q", text="label:N")
        scatter_4a = (
            base.mark_circle(opacity=0.8)
            .encode(
                x=alt.X("log10_kappa:Q", title="log₁₀(κ) — Signature Collinearity", scale=alt.Scale(domain=[0.4, 2.6])),
                y=alt.Y("melanoma_multivariate_auc:Q", title="Melanoma Response AUC", scale=alt.Scale(domain=[0.54, 0.74])),
                color=alt.Color("reference_category:N", title="Strategy", sort=cat_order),
                size=alt.Size("resolution:Q", title="Resolution", scale=alt.Scale(range=[30, 220])),
                shape=alt.Shape("reference_category:N", title="Strategy", sort=cat_order),
                tooltip=["reference_type", "reference_category", "resolution", "melanoma_multivariate_auc", "condition_number", "log10_kappa"],
            )
            .properties(title="A. Goldilocks Scatter: AUC vs log₁₀(κ)", width=350, height=270)
        )
        kappa30_line = alt.Chart(pd.DataFrame({"x": [np.log10(30)]})).mark_rule(
            strokeDash=[4, 3], color="#ff7f00", opacity=0.7
        ).encode(x="x:Q")
        kappa100_line = alt.Chart(pd.DataFrame({"x": [np.log10(100)]})).mark_rule(
            strokeDash=[4, 3], color="#e41a1c", opacity=0.7
        ).encode(x="x:Q")
        panel_a = alt.layer(opt_rect, opt_label, scatter_4a, kappa30_line, kappa100_line).resolve_scale(color="independent", size="independent", shape="independent")

        # Panel 4B: AUC × Condition-number heatmap — reference_type × resolution
        df_heat = df.with_columns(
            pl.col("resolution").cast(pl.Utf8).alias("res_str")
        ).sort(["reference_category", "reference_type", "resolution"]).to_pandas()
        panel_b = (
            alt.Chart(df_heat)
            .mark_rect()
            .encode(
                x=alt.X("res_str:N", title="Resolution", sort=sorted(df_heat["res_str"].unique().tolist())),
                y=alt.Y("reference_type:N", title="Reference",
                        sort=df_heat.drop_duplicates("reference_type").sort_values("reference_category")["reference_type"].tolist()),
                color=alt.Color("melanoma_multivariate_auc:Q", title="Melanoma AUC",
                                scale=alt.Scale(scheme="viridis", domain=[0.5, 0.85])),
                tooltip=["reference_type", "reference_category", "res_str", "melanoma_multivariate_auc", "condition_number"],
            )
            .properties(title="B. AUC Heatmap: Reference × Resolution", width=280, height=320)
        )

        # Panel 4C: Efficiency metric trajectory — AUC / log10(κ) by resolution
        panel_c = (
            base.mark_line(strokeWidth=2.0, point=alt.OverlayMarkDef(size=50, filled=True))
            .encode(
                x=alt.X("resolution:Q", title="Leiden Resolution", axis=alt.Axis(tickMinStep=0.25)),
                y=alt.Y("efficiency:Q", title="AUC / log₁₀(κ)  — Efficiency"),
                color=alt.Color("reference_type:N", title="Reference", scale=color_scale, legend=alt.Legend(columns=2)),
                tooltip=["reference_type", "reference_category", "resolution", "efficiency", "melanoma_multivariate_auc", "condition_number"],
            )
            .properties(title="C. Efficiency Metric: AUC / log₁₀(κ) by Resolution", width=330, height=270)
        )

        chart = (
            alt.hconcat(panel_a, panel_b, panel_c)
            .resolve_scale(color="independent", size="independent", shape="independent")
            .properties(title="Figure 4 — Pareto Tradeoff Optimization: Predictive Power vs. Signature Collinearity")
        )

        out_path = results_dir / "fig4_pareto_tradeoff_optimization.svg"
        results_dir.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_path))
        png_path = results_dir / "fig4_pareto_tradeoff_optimization.png"
        try:
            chart.save(str(png_path), scale_factor=2.0)
        except Exception:
            pass
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to generate Figure 4: {exc}")


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
