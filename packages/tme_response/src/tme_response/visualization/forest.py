"""Altair cross-cohort summary forest and error-bar visualizations with SVG export."""

from __future__ import annotations

import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


def create_summary_forest_plot(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair forest plot showing mean metric and 95% CI across cohorts."""
    try:
        is_pr = metric == "pr_auc"
        metric_col = "pr_auc" if is_pr else "roc_auc"
        lower_col = "pr_auc_ci_lower" if is_pr else "roc_auc_ci_lower"
        upper_col = "pr_auc_ci_upper" if is_pr else "roc_auc_ci_upper"

        default_title = (
            "Cross-Cohort Summary Predictor Performance (Mean PR-AUC ± 95% CI)"
            if is_pr
            else "Cross-Cohort Summary Predictor Performance (Mean ROC-AUC ± 95% CI)"
        )
        plot_title = title or default_title

        # Aggregate mean and CI bounds per predictor
        summary_df = (
            benchmark_df.group_by(["predictor_name", "category"])
            .agg([
                pl.col(metric_col).mean().alias("mean_metric"),
                pl.col(lower_col).mean().alias("mean_ci_lower"),
                pl.col(upper_col).mean().alias("mean_ci_upper"),
                pl.len().alias("n_cohorts"),
            ])
            .sort("mean_metric", descending=True)
        )

        domain = [0.0, 0.85] if is_pr else [0.35, 0.85]
        ref_val = float(benchmark_df["baseline_prevalence"].mean()) if is_pr and "baseline_prevalence" in benchmark_df.columns else 0.50

        base = alt.Chart(summary_df).encode(
            y=alt.Y("predictor_name:N", title="Predictor / Biomarker", sort="-x"),
            color=alt.Color("category:N", title="Category"),
        )

        error_bars = base.mark_rule(strokeWidth=2).encode(
            x=alt.X(
                "mean_ci_lower:Q",
                title="Mean PR-AUC (95% CI)" if is_pr else "Mean ROC-AUC (95% CI)",
                scale=alt.Scale(domain=domain),
            ),
            x2="mean_ci_upper:Q",
        )

        points = base.mark_point(filled=True, size=70).encode(
            x="mean_metric:Q",
            tooltip=[
                "predictor_name",
                "category",
                alt.Tooltip("mean_metric:Q", format=".3f", title="Mean Metric"),
                alt.Tooltip("mean_ci_lower:Q", format=".3f", title="CI Lower"),
                alt.Tooltip("mean_ci_upper:Q", format=".3f", title="CI Upper"),
                "n_cohorts",
            ],
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"x": [ref_val]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        chart = (error_bars + points + ref_line).properties(
            title=plot_title,
            width=380,
            height=alt.Step(26),
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create summary forest plot: {exc}")


def create_cohort_predictability_forest_plot(
    predictability_df: pl.DataFrame,
    title: str = "Cohort Predictability Index (Mean RNA Biomarker ROC-AUC ± SD)",
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair horizontal forest plot ranking cohorts by predictability."""
    try:
        if predictability_df.is_empty():
            return Failure("Predictability DataFrame is empty.")

        # Compute lower and upper error bounds
        annotated = predictability_df.with_columns([
            pl.max_horizontal([pl.col("mean_roc_auc_rna") - pl.col("std_roc_auc_rna"), pl.lit(0.20)]).alias("ci_lower"),
            pl.min_horizontal([pl.col("mean_roc_auc_rna") + pl.col("std_roc_auc_rna"), pl.lit(1.00)]).alias("ci_upper"),
            pl.format("{} (AUC: {}, prev: {}%)", pl.col("cohort_id"), pl.col("mean_roc_auc_rna").round(2), (pl.col("baseline_prevalence") * 100).round(0)).alias("label"),
        ])

        base = alt.Chart(annotated).encode(
            y=alt.Y(
                "cohort_id:N",
                title="Cohort / Clinical Pool",
                sort=alt.EncodingSortField(field="mean_roc_auc_rna", order="descending"),
            ),
            color=alt.Color("cancer_type:N", title="Cancer Type"),
        )

        error_bars = base.mark_rule(strokeWidth=2.2).encode(
            x=alt.X("ci_lower:Q", title="Mean RNA-Seq ROC-AUC (± 1 SD)", scale=alt.Scale(domain=[0.30, 0.90])),
            x2="ci_upper:Q",
        )

        points = base.mark_point(filled=True, size=85).encode(
            x="mean_roc_auc_rna:Q",
            tooltip=[
                "cohort_id",
                "cancer_type",
                "pooling_strategy",
                "n_samples",
                "n_responders",
                alt.Tooltip("baseline_prevalence:Q", format=".1%", title="Prevalence"),
                alt.Tooltip("mean_roc_auc_rna:Q", format=".3f", title="Mean RNA AUC"),
                alt.Tooltip("std_roc_auc_rna:Q", format=".3f", title="SD RNA AUC"),
                alt.Tooltip("mean_pr_auc_rna:Q", format=".3f", title="Mean RNA PR-AUC"),
                alt.Tooltip("mean_delta_pr_auc_rna:Q", format=".3f", title="Mean ΔPR-AUC"),
                "best_predictor_rna",
                alt.Tooltip("max_roc_auc_rna:Q", format=".3f", title="Max RNA AUC"),
            ],
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"x": [0.50]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        chart = (error_bars + points + ref_line).properties(
            title=title,
            width=400,
            height=alt.Step(26),
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create cohort predictability forest plot: {exc}")


def create_pooling_comparison_chart(
    benchmark_df: pl.DataFrame,
    title: str = "Impact of Score Standardization on Pooled ROC-AUC",
) -> Result[alt.Chart, str]:
    """Create a publication-grade grouped comparison chart evaluating raw vs standardized pooling."""
    try:
        pool_df = benchmark_df.filter(
            pl.col("pooling_strategy").is_in(["raw", "standardized"])
        )
        if pool_df.is_empty():
            return Failure("No pooled cohort records with pooling_strategy in ['raw', 'standardized'].")

        # Aggregate mean ROC-AUC and average confidence bounds across predictors per pool and strategy
        summary_df = (
            pool_df.group_by(["cohort_id", "pooling_strategy"])
            .agg([
                pl.col("roc_auc").mean().alias("mean_auc"),
                pl.col("roc_auc_ci_lower").mean().alias("ci_lower"),
                pl.col("roc_auc_ci_upper").mean().alias("ci_upper"),
                pl.col("pr_auc").mean().alias("mean_pr"),
                pl.col("delta_pr_auc").mean().alias("mean_delta_pr"),
                pl.len().alias("n_predictors"),
            ])
            .with_columns(
                pl.format("{}", pl.col("mean_auc").round(3)).alias("display_text")
            )
            .sort("cohort_id")
        )

        base = alt.Chart(summary_df).encode(
            x=alt.X(
                "cohort_id:N",
                title="Combined Cohort Pool",
                axis=alt.Axis(labelAngle=0, titleFontSize=12, labelFontSize=11),
            ),
            xOffset=alt.XOffset("pooling_strategy:N", sort=["raw", "standardized"]),
            color=alt.Color(
                "pooling_strategy:N",
                title="Pooling Strategy",
                scale=alt.Scale(
                    domain=["raw", "standardized"],
                    range=["#8c96c6", "#2b8cbe"],
                ),
                legend=alt.Legend(titleFontSize=11, labelFontSize=10),
            ),
        )

        bars = base.mark_bar(width=28, cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
            y=alt.Y(
                "mean_auc:Q",
                title="Mean ROC-AUC across Predictors",
                scale=alt.Scale(domain=[0.50, 0.65], zero=False, clamp=True),
                axis=alt.Axis(titleFontSize=12, labelFontSize=11),
            ),
            tooltip=[
                "cohort_id",
                "pooling_strategy",
                alt.Tooltip("mean_auc:Q", format=".3f", title="Mean ROC-AUC"),
                alt.Tooltip("ci_lower:Q", format=".3f", title="CI Lower"),
                alt.Tooltip("ci_upper:Q", format=".3f", title="CI Upper"),
                alt.Tooltip("mean_pr:Q", format=".3f", title="Mean PR-AUC"),
                alt.Tooltip("mean_delta_pr:Q", format=".3f", title="Mean ΔPR-AUC"),
                "n_predictors",
            ],
        )

        error_bars = base.mark_rule(strokeWidth=1.8, color="#252525").encode(
            y="ci_lower:Q",
            y2="ci_upper:Q",
        )

        text = base.mark_text(
            baseline="bottom",
            dy=-5,
            fontSize=11,
            fontWeight="bold",
        ).encode(
            y="mean_auc:Q",
            text="display_text:N",
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"y": [0.50]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(y="y:Q")
        )

        chart = (
            (bars + error_bars + text + ref_line)
            .properties(
                title=title,
                width=360,
                height=260,
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create pooling comparison chart: {exc}")
