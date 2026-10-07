"""Altair declarative survival (C-index, Cox HR) and Decision Curve Analysis (DCA) visualizations."""

from __future__ import annotations

from typing import Sequence
import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


def create_c_index_forest_plot(
    survival_df: pl.DataFrame,
    endpoint: str = "OS",
    title: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair forest plot showing Harrell's C-index (± 95% bootstrap CI).

    Parameters
    ----------
    survival_df : pl.DataFrame
        DataFrame with columns: predictor_name, endpoint, c_index, ci_lower, ci_upper, p_value, category (optional).
    endpoint : str
        Endpoint to filter on, typically 'OS' or 'PFS'.
    title : str | None
        Chart title. Defaults to standard endpoint title.
    """
    try:
        sub_df = survival_df.filter(pl.col("endpoint") == endpoint)
        if sub_df.is_empty():
            return Failure(f"No survival records found for endpoint '{endpoint}'.")

        # If data is cohort-specific, average across cohorts per predictor; if already summarized, use as-is
        group_cols = ["predictor_name"]
        if "category" in sub_df.columns:
            group_cols.append("category")

        if "cohort_id" in sub_df.columns and sub_df["cohort_id"].n_unique() > 1:
            summary_df = (
                sub_df.group_by(group_cols)
                .agg([
                    pl.col("c_index").mean().alias("mean_c_index"),
                    pl.col("ci_lower").mean().alias("mean_ci_lower"),
                    pl.col("ci_upper").mean().alias("mean_ci_upper"),
                    pl.len().alias("n_cohorts"),
                ])
                .sort("mean_c_index", descending=True)
            )
            val_col = "mean_c_index"
            lower_col = "mean_ci_lower"
            upper_col = "mean_ci_upper"
        else:
            summary_df = sub_df.sort("c_index", descending=True)
            val_col = "c_index"
            lower_col = "ci_lower"
            upper_col = "ci_upper"

        plot_title = title or f"Predictor Discrimination for {endpoint} (Harrell's C-index ± 95% CI)"

        base = alt.Chart(summary_df).encode(
            y=alt.Y("predictor_name:N", title="Predictor / Biomarker", sort="-x"),
        )
        if "category" in summary_df.columns:
            base = base.encode(color=alt.Color("category:N", title="Category"))

        error_bars = base.mark_rule(strokeWidth=2.2).encode(
            x=alt.X(
                f"{lower_col}:Q",
                title=f"Harrell's C-index ({endpoint}) [95% Bootstrap CI]",
                scale=alt.Scale(domain=[0.35, 0.80]),
                axis=alt.Axis(titleFontSize=12, labelFontSize=11),
            ),
            x2=f"{upper_col}:Q",
        )

        tooltip_cols = [
            "predictor_name",
            alt.Tooltip(f"{val_col}:Q", format=".3f", title="C-index"),
            alt.Tooltip(f"{lower_col}:Q", format=".3f", title="CI Lower"),
            alt.Tooltip(f"{upper_col}:Q", format=".3f", title="CI Upper"),
        ]
        if "category" in summary_df.columns:
            tooltip_cols.append("category")
        if "n_cohorts" in summary_df.columns:
            tooltip_cols.append("n_cohorts")

        points = base.mark_point(filled=True, size=75).encode(
            x=f"{val_col}:Q",
            tooltip=tooltip_cols,
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"x": [0.50]}))
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
        return Failure(f"Failed to create C-index forest plot: {exc}")


def create_cox_hr_forest_plot(
    cox_df: pl.DataFrame,
    endpoint: str = "OS",
    title: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair forest plot of Cox Proportional Hazard Ratios on log scale.

    Hazard Ratio < 1.0 signifies higher score correlates with reduced hazard / prolonged survival.
    Hazard Ratio > 1.0 signifies higher score correlates with elevated hazard / shortened survival.
    """
    try:
        sub_df = cox_df.filter(pl.col("endpoint") == endpoint)
        if sub_df.is_empty():
            return Failure(f"No Cox regression records found for endpoint '{endpoint}'.")

        group_cols = ["predictor_name"]
        if "category" in sub_df.columns:
            group_cols.append("category")

        if "cohort_id" in sub_df.columns and sub_df["cohort_id"].n_unique() > 1:
            summary_df = (
                sub_df.group_by(group_cols)
                .agg([
                    # Geometric mean for hazard ratios
                    pl.col("hazard_ratio").log().mean().exp().alias("mean_hr"),
                    pl.col("hr_ci_lower").log().mean().exp().alias("mean_hr_lower"),
                    pl.col("hr_ci_upper").log().mean().exp().alias("mean_hr_upper"),
                    pl.len().alias("n_cohorts"),
                ])
                .sort("mean_hr")
            )
            hr_col = "mean_hr"
            lower_col = "mean_hr_lower"
            upper_col = "mean_hr_upper"
        else:
            summary_df = sub_df.sort("hazard_ratio")
            hr_col = "hazard_ratio"
            lower_col = "hr_ci_lower"
            upper_col = "hr_ci_upper"

        plot_title = title or f"Hazard Ratio per 1-SD Increase in Score for {endpoint} (Univariable Cox Model)"

        base = alt.Chart(summary_df).encode(
            y=alt.Y("predictor_name:N", title="Predictor / Biomarker", sort="x"),
        )
        if "category" in summary_df.columns:
            base = base.encode(color=alt.Color("category:N", title="Category"))

        error_bars = base.mark_rule(strokeWidth=2.2).encode(
            x=alt.X(
                f"{lower_col}:Q",
                title=f"Hazard Ratio (95% CI) per 1-SD [{endpoint}]",
                scale=alt.Scale(type="log", domain=[0.25, 3.5]),
                axis=alt.Axis(titleFontSize=12, labelFontSize=11),
            ),
            x2=f"{upper_col}:Q",
        )

        tooltip_cols = [
            "predictor_name",
            alt.Tooltip(f"{hr_col}:Q", format=".3f", title="Hazard Ratio"),
            alt.Tooltip(f"{lower_col}:Q", format=".3f", title="95% CI Lower"),
            alt.Tooltip(f"{upper_col}:Q", format=".3f", title="95% CI Upper"),
        ]
        if "category" in summary_df.columns:
            tooltip_cols.append("category")
        if "n_cohorts" in summary_df.columns:
            tooltip_cols.append("n_cohorts")

        points = base.mark_point(filled=True, size=75).encode(
            x=f"{hr_col}:Q",
            tooltip=tooltip_cols,
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"x": [1.0]}))
            .mark_rule(strokeDash=[4, 4], color="black", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        chart = (error_bars + points + ref_line).properties(
            title=plot_title,
            width=380,
            height=alt.Step(26),
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create Cox HR forest plot: {exc}")


def create_dca_net_benefit_chart(
    dca_df: pl.DataFrame,
    title: str = "Decision Curve Analysis (Net Clinical Benefit)",
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair Decision Curve Analysis chart comparing models against Treat All/None."""
    try:
        if dca_df.is_empty():
            return Failure("DCA DataFrame is empty.")

        # Identify models vs reference baselines
        strategies = dca_df["strategy"].unique().to_list()
        model_strategies = [s for s in strategies if s not in ("Treat All", "Treat None")]

        # Determine y scale range
        max_nb = float(dca_df["net_benefit"].max())
        y_max = max(0.40, round(max_nb + 0.05, 2))

        # Model lines
        model_df = dca_df.filter(pl.col("strategy").is_in(model_strategies))
        model_lines = (
            alt.Chart(model_df)
            .mark_line(strokeWidth=2.2)
            .encode(
                x=alt.X(
                    "threshold:Q",
                    title="Threshold Probability (pt)",
                    scale=alt.Scale(domain=[0.05, 0.75]),
                    axis=alt.Axis(titleFontSize=12, labelFontSize=11),
                ),
                y=alt.Y(
                    "net_benefit:Q",
                    title="Net Benefit (True Positives - Weight * False Positives)",
                    scale=alt.Scale(domain=[-0.05, y_max]),
                    axis=alt.Axis(titleFontSize=12, labelFontSize=11),
                ),
                color=alt.Color(
                    "strategy:N",
                    title="Strategy / Model",
                    legend=alt.Legend(symbolOpacity=1.0, titleFontSize=11, labelFontSize=10),
                ),
                tooltip=[
                    "strategy",
                    alt.Tooltip("threshold:Q", format=".2f", title="Threshold pt"),
                    alt.Tooltip("net_benefit:Q", format=".3f", title="Net Benefit"),
                    alt.Tooltip("interventions_avoided_per_100:Q", format=".1f", title="Avoided / 100 pts"),
                ],
            )
        )

        # Baseline: Treat All
        treat_all_df = dca_df.filter(pl.col("strategy") == "Treat All")
        treat_all_line = (
            alt.Chart(treat_all_df)
            .mark_line(strokeWidth=1.8, strokeDash=[5, 5], color="#6b7280")
            .encode(
                x="threshold:Q",
                y="net_benefit:Q",
                tooltip=[
                    "strategy",
                    alt.Tooltip("threshold:Q", format=".2f", title="Threshold pt"),
                    alt.Tooltip("net_benefit:Q", format=".3f", title="Net Benefit"),
                ],
            )
        )

        # Baseline: Treat None (y = 0)
        treat_none_line = (
            alt.Chart(pl.DataFrame({"y": [0.0]}))
            .mark_rule(strokeWidth=1.5, color="#111827")
            .encode(y="y:Q")
        )

        chart = (
            (model_lines + treat_all_line + treat_none_line)
            .properties(
                title=title,
                width=450,
                height=320,
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create DCA chart: {exc}")
