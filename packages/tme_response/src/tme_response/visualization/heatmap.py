"""Altair publication-grade ROC-AUC heatmap visualizations with SVG export."""

from __future__ import annotations

import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


def create_auc_heatmap(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str = "Immunotherapy Response Predictor Performance",
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair heatmap displaying discrimination metrics across predictors and cohorts.

    Supports metric="roc_auc", metric="pr_auc", or metric="delta_pr_auc".
    Includes text overlays with formatted metric values and significance asterisks (* p < 0.05, ** p < 0.01).
    """
    try:
        if metric not in benchmark_df.columns:
            return Failure(f"Metric '{metric}' not found in benchmark data.")

        p_col = (
            "p_value_prauc"
            if metric in ("pr_auc", "delta_pr_auc") and "p_value_prauc" in benchmark_df.columns
            else "p_value_vs_half"
        )
        has_p = p_col in benchmark_df.columns

        # Add label text: e.g. "0.72*"
        if has_p:
            df_annotated = benchmark_df.with_columns(
                pl.when(pl.col(p_col) < 0.01)
                .then(pl.format("{}**", pl.col(metric).round(2)))
                .when(pl.col(p_col) < 0.05)
                .then(pl.format("{}*", pl.col(metric).round(2)))
                .otherwise(pl.format("{}", pl.col(metric).round(2)))
                .alias("display_text")
            )
        else:
            df_annotated = benchmark_df.with_columns(
                pl.format("{}", pl.col(metric).round(2)).alias("display_text")
            )

        match metric:
            case "roc_auc":
                color_scale = alt.Scale(scheme="redblue", domain=[0.35, 0.85], domainMid=0.50)
                color_title = "ROC-AUC"
                text_color_cond = (alt.datum[metric] > 0.70) | (alt.datum[metric] < 0.40)
            case "pr_auc":
                color_scale = alt.Scale(scheme="viridis", domain=[0.0, 1.0])
                color_title = "PR-AUC"
                text_color_cond = (alt.datum[metric] > 0.65) | (alt.datum[metric] < 0.25)
            case "delta_pr_auc":
                color_scale = alt.Scale(scheme="redblue", domain=[-0.25, 0.35], domainMid=0.0)
                color_title = "ΔPR-AUC"
                text_color_cond = (alt.datum[metric] > 0.18) | (alt.datum[metric] < -0.10)
            case _:
                color_scale = alt.Scale(scheme="blues")
                color_title = metric
                text_color_cond = alt.datum[metric] > 0.5

        base = alt.Chart(df_annotated).encode(
            x=alt.X("cohort_id:N", title="Cohort", axis=alt.Axis(labelAngle=-45)),
            y=alt.Y(
                "predictor_name:N",
                title="Predictor",
                sort=alt.EncodingSortField(field=metric, op="mean", order="descending"),
            ),
        )

        tooltip_candidates = [
            "cohort_id",
            "cancer_type",
            "predictor_name",
            "category",
            "time_stratum",
            "response_stratum",
            "pooling_strategy",
            "roc_auc",
            "roc_auc_ci_lower",
            "roc_auc_ci_upper",
            "pr_auc",
            "pr_auc_ci_lower",
            "pr_auc_ci_upper",
            "delta_pr_auc",
            "baseline_prevalence",
            "n_samples",
            "n_responders",
        ]
        active_tooltips = [c for c in tooltip_candidates if c in benchmark_df.columns]

        heatmap = base.mark_rect().encode(
            color=alt.Color(f"{metric}:Q", title=color_title, scale=color_scale),
            tooltip=active_tooltips,
        )

        text = base.mark_text(baseline="middle", fontSize=9.5, fontWeight="bold").encode(
            text="display_text:N",
            color=alt.condition(
                text_color_cond,
                alt.value("white"),
                alt.value("black"),
            ),
        )

        chart = (heatmap + text).properties(
            title=title,
            width=alt.Step(55),
            height=alt.Step(26),
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create AUC heatmap: {exc}")
