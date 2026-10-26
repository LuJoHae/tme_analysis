"""Altair cross-cohort summary forest and error-bar visualizations with SVG export."""

from __future__ import annotations

import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


from .display_labels import with_display_labels


def create_summary_forest_plot(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str | None = None,
    width: int | None = None,
    domain: list[float] | None = None,
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

        # Enrich with human-readable display labels
        df_enriched = with_display_labels(benchmark_df)
        pred_label_col = "predictor_display" if "predictor_display" in df_enriched.columns else "predictor_name"

        # Aggregate mean and CI bounds per predictor
        summary_df = (
            df_enriched.group_by([pred_label_col, "category"])
            .agg([
                pl.col(metric_col).mean().alias("mean_metric"),
                pl.col(lower_col).mean().alias("mean_ci_lower"),
                pl.col(upper_col).mean().alias("mean_ci_upper"),
                pl.len().alias("n_cohorts"),
            ])
            .sort("mean_metric", descending=True)
        )

        default_domain = [0.25, 0.60] if is_pr else [0.44, 0.80]
        plot_domain = domain if domain is not None else default_domain
        ref_val = float(benchmark_df["baseline_prevalence"].mean()) if is_pr and "baseline_prevalence" in benchmark_df.columns else 0.50

        base = alt.Chart(summary_df).encode(
            y=alt.Y(
                f"{pred_label_col}:N",
                title="Predictor / Biomarker",
                sort="-x",
                axis=alt.Axis(labelLimit=380, labelFontSize=10, titlePadding=15),
            ),
            color=alt.Color("category:N", title="Category"),
        )

        error_bars = base.mark_rule(strokeWidth=2).encode(
            x=alt.X(
                "mean_ci_lower:Q",
                title="Mean PR-AUC (95% CI)" if is_pr else "Mean ROC-AUC (95% CI)",
                scale=alt.Scale(domain=plot_domain),
            ),
            x2="mean_ci_upper:Q",
        )

        points = base.mark_point(filled=True, size=70).encode(
            x="mean_metric:Q",
            tooltip=[
                pred_label_col,
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
            title=alt.TitleParams(text=plot_title, fontSize=13, anchor="start"),
            width=width if width is not None else 260,
            height=alt.Step(26),
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create summary forest plot: {exc}")


def create_cohort_predictability_forest_plot(
    predictability_df: pl.DataFrame,
    title: str = "Cohort Predictability Index (Mean RNA Biomarker ROC-AUC ± SD)",
    width: int | None = None,
    domain: list[float] | None = None,
    show_y_title: bool = False,
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair horizontal forest plot ranking cohorts by predictability."""
    try:
        if predictability_df.is_empty():
            return Failure("Predictability DataFrame is empty.")

        # Enrich with human-readable display labels
        df_enriched = with_display_labels(predictability_df)
        cohort_label_col = "cohort_display" if "cohort_display" in df_enriched.columns else "cohort_id"

        # Compute lower and upper error bounds
        annotated = df_enriched.with_columns([
            pl.max_horizontal([pl.col("mean_roc_auc_rna") - pl.col("std_roc_auc_rna"), pl.lit(0.20)]).alias("ci_lower"),
            pl.min_horizontal([pl.col("mean_roc_auc_rna") + pl.col("std_roc_auc_rna"), pl.lit(1.00)]).alias("ci_upper"),
            pl.format("{} (AUC: {}, prev: {}%)", pl.col(cohort_label_col), pl.col("mean_roc_auc_rna").round(2), (pl.col("baseline_prevalence") * 100).round(0)).alias("label"),
        ])

        base = alt.Chart(annotated).encode(
            y=alt.Y(
                f"{cohort_label_col}:N",
                title="Cohort / Clinical Pool" if show_y_title else None,
                sort=alt.EncodingSortField(field="mean_roc_auc_rna", order="descending"),
                axis=alt.Axis(labelLimit=380, labelFontSize=10, titlePadding=110 if show_y_title else 0),
            ),
            color=alt.Color("cancer_type:N", title="Cancer Type"),
        )

        plot_domain = domain if domain is not None else [0.38, 0.88]
        plot_width = width if width is not None else 280

        error_bars = base.mark_rule(strokeWidth=2.2).encode(
            x=alt.X("ci_lower:Q", title="Mean RNA-Seq ROC-AUC (± 1 SD)", scale=alt.Scale(domain=plot_domain)),
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
            width=plot_width,
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

        # Enrich with human-readable display labels
        pool_enriched = with_display_labels(pool_df)
        c_label = "cohort_display" if "cohort_display" in pool_enriched.columns else "cohort_id"

        # Aggregate mean ROC-AUC and average confidence bounds across predictors per pool and strategy
        summary_df = (
            pool_enriched.group_by([c_label, "cohort_id", "pooling_strategy"])
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
            .sort(c_label)
        )

        base = alt.Chart(summary_df).encode(
            x=alt.X(
                f"{c_label}:N",
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


def create_combined_summary_forest_plot(
    benchmark_df: pl.DataFrame,
    title: str | None = "Pan-Cancer Meta-Analytic ICB Response Performance Across 32 Biomarkers",
    panel_width: int = 180,
    roc_domain: list[float] | None = None,
    pr_domain: list[float] | None = None,
    sort_by: str = "roc_auc",
    spacing: int = 15,
) -> Result[alt.HConcatChart, str]:
    """Create a publication-grade two-panel forest plot (ROC-AUC and PR-AUC) with a single shared y-axis.

    Panel A (Left): Discrimination (ROC-AUC) with complete predictor y-axis labels.
    Panel B (Right): Precision (PR-AUC) with y-axis labels and ticks suppressed.
    """
    try:
        if benchmark_df.is_empty():
            return Failure("Benchmark dataframe is empty.")

        # Enrich with human-readable display labels
        df_enriched = with_display_labels(benchmark_df)
        pred_label_col = "predictor_display" if "predictor_display" in df_enriched.columns else "predictor_name"

        # Required columns check
        req_cols = {"roc_auc", "roc_auc_ci_lower", "roc_auc_ci_upper", "pr_auc", "pr_auc_ci_lower", "pr_auc_ci_upper"}
        missing = req_cols - set(df_enriched.columns)
        if missing:
            return Failure(f"Missing required metric columns: {missing}")

        # Aggregate mean and CI bounds per predictor
        agg_exprs = [
            pl.col("roc_auc").mean().alias("roc_mean"),
            pl.col("roc_auc_ci_lower").mean().alias("roc_lower"),
            pl.col("roc_auc_ci_upper").mean().alias("roc_upper"),
            pl.col("pr_auc").mean().alias("pr_mean"),
            pl.col("pr_auc_ci_lower").mean().alias("pr_lower"),
            pl.col("pr_auc_ci_upper").mean().alias("pr_upper"),
            pl.len().alias("n_cohorts"),
        ]
        if "baseline_prevalence" in df_enriched.columns:
            agg_exprs.append(pl.col("baseline_prevalence").mean().alias("baseline_prevalence"))

        summary_df = df_enriched.group_by([pred_label_col, "category"]).agg(agg_exprs)

        # Determine predictor ordering for y-axis (top-to-bottom in Altair domain)
        match sort_by:
            case "roc_auc":
                sorted_df = summary_df.sort("roc_mean", descending=True)
            case "pr_auc":
                sorted_df = summary_df.sort("pr_mean", descending=True)
            case "alphabetical":
                sorted_df = summary_df.sort(pred_label_col, descending=False)
            case _:
                sorted_df = summary_df.sort("roc_mean", descending=True)

        predictor_order = sorted_df[pred_label_col].to_list()

        # Scale domains and reference lines
        roc_plot_domain = roc_domain if roc_domain is not None else [0.45, 0.80]
        pr_plot_domain = pr_domain if pr_domain is not None else [0.25, 0.60]

        ref_val_roc = 0.50
        ref_val_pr = (
            float(sorted_df["baseline_prevalence"].mean())
            if "baseline_prevalence" in sorted_df.columns
            else 0.314
        )

        legend = alt.Legend(
            orient="bottom",
            direction="horizontal",
            title="Category",
            titleOrient="left",
            titleFontSize=12.5,
            titleFontWeight="bold",
            labelFontSize=11.5,
            symbolSize=100,
            padding=12,
            columnPadding=16,
        )

        # Panel A: ROC-AUC (Left Panel, with Y-axis labels)
        base_roc = alt.Chart(sorted_df).encode(
            y=alt.Y(
                f"{pred_label_col}:N",
                title="Predictor / Biomarker",
                scale=alt.Scale(domain=predictor_order),
                axis=alt.Axis(
                    labelLimit=380,
                    labelFontSize=11.5,
                    titleFontSize=13,
                    titleFontWeight="bold",
                    titlePadding=15,
                ),
            ),
            color=alt.Color("category:N", title="Category", legend=legend),
        )

        error_bars_roc = base_roc.mark_rule(strokeWidth=2).encode(
            x=alt.X(
                "roc_lower:Q",
                title="Mean ROC-AUC (95% CI)",
                scale=alt.Scale(domain=roc_plot_domain),
                axis=alt.Axis(
                    values=[0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80],
                    labelFontSize=11,
                    titleFontSize=12.5,
                    titleFontWeight="bold",
                    titlePadding=10,
                ),
            ),
            x2="roc_upper:Q",
        )

        points_roc = base_roc.mark_point(filled=True, size=75).encode(
            x="roc_mean:Q",
            tooltip=[
                pred_label_col,
                "category",
                alt.Tooltip("roc_mean:Q", format=".3f", title="Mean ROC-AUC"),
                alt.Tooltip("roc_lower:Q", format=".3f", title="ROC-AUC CI Lower"),
                alt.Tooltip("roc_upper:Q", format=".3f", title="ROC-AUC CI Upper"),
                "n_cohorts",
            ],
        )

        ref_line_roc = (
            alt.Chart(pl.DataFrame({"x": [ref_val_roc]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        panel_roc = (error_bars_roc + points_roc + ref_line_roc).properties(
            width=panel_width,
            height=alt.Step(26),
        )

        # Panel B: PR-AUC (Right Panel, Y-axis labels and ticks suppressed)
        base_pr = alt.Chart(sorted_df).encode(
            y=alt.Y(
                f"{pred_label_col}:N",
                title=None,
                scale=alt.Scale(domain=predictor_order),
                axis=alt.Axis(labels=False, ticks=False, domain=False),
            ),
            color=alt.Color("category:N", title="Category", legend=legend),
        )

        error_bars_pr = base_pr.mark_rule(strokeWidth=2).encode(
            x=alt.X(
                "pr_lower:Q",
                title="Mean PR-AUC (95% CI)",
                scale=alt.Scale(domain=pr_plot_domain),
                axis=alt.Axis(
                    values=[0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60],
                    labelFontSize=11,
                    titleFontSize=12.5,
                    titleFontWeight="bold",
                    titlePadding=10,
                ),
            ),
            x2="pr_upper:Q",
        )

        points_pr = base_pr.mark_point(filled=True, size=75).encode(
            x="pr_mean:Q",
            tooltip=[
                pred_label_col,
                "category",
                alt.Tooltip("pr_mean:Q", format=".3f", title="Mean PR-AUC"),
                alt.Tooltip("pr_lower:Q", format=".3f", title="PR-AUC CI Lower"),
                alt.Tooltip("pr_upper:Q", format=".3f", title="PR-AUC CI Upper"),
                "n_cohorts",
            ],
        )

        ref_line_pr = (
            alt.Chart(pl.DataFrame({"x": [ref_val_pr]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        panel_pr = (error_bars_pr + points_pr + ref_line_pr).properties(
            width=panel_width,
            height=alt.Step(26),
        )

        combined = (
            alt.hconcat(panel_roc, panel_pr, spacing=spacing)
            .resolve_scale(y="shared")
            .resolve_legend(color="shared")
        )

        if title:
            combined = combined.properties(
                title=alt.TitleParams(text=title, fontSize=15, fontWeight="bold", anchor="start")
            )

        return Success(combined)
    except Exception as exc:
        return Failure(f"Failed to create combined summary forest plot: {exc}")
