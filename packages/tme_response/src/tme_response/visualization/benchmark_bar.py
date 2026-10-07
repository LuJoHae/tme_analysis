"""Altair benchmark comparison bar charts with SVG export."""

from __future__ import annotations

import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


def create_benchmark_bar_chart(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str = "Cross-Cohort Predictor Performance Comparison",
) -> Result[alt.Chart, str]:
    """Pure functional constructor of grouped Altair bar chart comparing predictor performance."""
    try:
        if metric not in benchmark_df.columns:
            return Failure(f"Metric '{metric}' not found in benchmark DataFrame.")

        chart = (
            alt.Chart(benchmark_df)
            .mark_bar()
            .encode(
                x=alt.X("predictor_name:N", title="Predictor / Biomarker", sort="-y"),
                y=alt.Y(f"{metric}:Q", title=metric.upper(), scale=alt.Scale(domain=[0, 1])),
                color=alt.Color("category:N", title="Category"),
                column=alt.Column("cohort_id:N", title="Validation Cohort"),
                tooltip=["cohort_id", "predictor_name", "category", metric, "n_samples"],
            )
            .properties(
                title=title,
                width=160,
                height=260,
            )
        )
        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create benchmark bar chart: {exc}")
