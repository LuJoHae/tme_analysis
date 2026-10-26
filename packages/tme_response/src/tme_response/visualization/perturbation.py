"""Altair publication-grade decay curves and resilience ranking visualizations for perturbations."""

from __future__ import annotations

from typing import Literal, assert_never
import altair as alt
import polars as pl
from returns.result import Failure, Result, Success

from .display_labels import get_cohort_display_label, with_display_labels

UncertaintyMode = Literal["shaded", "errorbars", "none"]


def create_perturbation_decay_chart(
    perturbation_df: pl.DataFrame,
    perturbation_type: str,
    metric: str = "roc_auc",
    uncertainty: UncertaintyMode = "shaded",
    title: str | None = None,
    subtitle: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade multi-line decay plot showing metric vs intensity.

    Supports 3 uncertainty modes:
    - 'shaded': Semi-transparent confidence area band.
    - 'errorbars': Explicit vertical error whiskers for each point.
    - 'none': Clean curves showing point estimates only.

    The legend strictly displays vibrant point colors matching marker points,
    unaffected by shaded area opacity.
    """
    try:
        df_enriched = with_display_labels(perturbation_df)
        sub_df = df_enriched.filter(pl.col("perturbation_type") == perturbation_type)
        if sub_df.is_empty():
            return Failure(f"No records found for perturbation_type='{perturbation_type}'.")

        pred_label_col = "predictor_display" if "predictor_display" in sub_df.columns else "predictor_name"
        is_pr = metric == "pr_auc"
        metric_col = "pr_auc" if is_pr else "roc_auc"
        lower_col = "pr_auc_ci_lower" if is_pr else "roc_auc_ci_lower"
        upper_col = "pr_auc_ci_upper" if is_pr else "roc_auc_ci_upper"

        x_title_map = {
            "jitter": "Multiplicative Expression Jitter (σ)",
            "dropout": "Gene Dropout Rate (Probability)",
            "dilution": "Immune Effector Dilution Factor (α)",
            "label_noise": "Response Label Noise Rate (η)",
        }
        x_title = x_title_map.get(perturbation_type, "Perturbation Intensity")

        cohort_ids = sub_df["cohort_id"].unique().to_list()
        cohort_str = get_cohort_display_label(cohort_ids[0]) if len(cohort_ids) == 1 else f"{len(cohort_ids)} Cohorts"

        p_name = perturbation_type.replace("_", " ").title()
        main_title = title or f"Predictor Decay under {p_name} ({cohort_str})"

        metric_name = "PR-AUC" if is_pr else "ROC-AUC"
        if subtitle is not None:
            sub_title = subtitle
        else:
            match uncertainty:
                case "shaded":
                    sub_title = f"Cohort: {cohort_str} | Y-Axis: Empirical {metric_name} (Shaded Band = 95% Bootstrap CI; Dashed = 0.50 Chance)"
                case "errorbars":
                    sub_title = f"Cohort: {cohort_str} | Y-Axis: Empirical {metric_name} (Errorbars = 95% Bootstrap CI; Dashed = 0.50 Chance)"
                case "none":
                    sub_title = f"Cohort: {cohort_str} | Y-Axis: Empirical {metric_name} (Point Estimates; Dashed = 0.50 Chance)"
                case _ as unreachable:
                    assert_never(unreachable)

        plot_title = alt.TitleParams(
            text=main_title,
            subtitle=sub_title,
            fontSize=13,
            subtitleFontSize=10.5,
            subtitleColor="#555555",
        )

        color_scale = alt.Scale(scheme="category10")
        pred_legend = alt.Legend(
            title="Predictor",
            symbolType="circle",
            symbolOpacity=1.0,
            symbolSize=70,
            symbolStrokeWidth=0,
            labelLimit=300,
        )

        base = alt.Chart(sub_df).encode(
            x=alt.X("intensity:Q", title=x_title),
            color=alt.Color(
                f"{pred_label_col}:N",
                title="Predictor",
                scale=color_scale,
                legend=pred_legend,
            ),
        )

        lines = base.mark_line(strokeWidth=2).encode(
            y=alt.Y(
                f"{metric_col}:Q",
                title="Empirical ROC-AUC" if not is_pr else "Empirical PR-AUC",
                scale=alt.Scale(domain=[0.35, 0.92], zero=False),
            ),
        )

        points = base.mark_point(size=55, filled=True, opacity=1.0).encode(
            y=f"{metric_col}:Q",
            tooltip=[
                "cohort_id",
                pred_label_col,
                "perturbation_type",
                alt.Tooltip("intensity:Q", format=".2f", title="Intensity"),
                alt.Tooltip(f"{metric_col}:Q", format=".3f", title="Score"),
                alt.Tooltip(f"{lower_col}:Q", format=".3f", title="CI Lower"),
                alt.Tooltip(f"{upper_col}:Q", format=".3f", title="CI Upper"),
            ],
        )

        areas = base.mark_area(opacity=0.18).encode(
            y=f"{lower_col}:Q",
            y2=f"{upper_col}:Q",
        )

        error_bars = base.mark_rule(strokeWidth=1.3, opacity=0.85).encode(
            y=f"{lower_col}:Q",
            y2=f"{upper_col}:Q",
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"y": [0.50]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(y="y:Q")
        )

        match uncertainty:
            case "shaded":
                core = areas + lines + points
            case "errorbars":
                core = lines + error_bars + points
            case "none":
                core = lines + points
            case _ as unreachable:
                assert_never(unreachable)

        chart = (
            (core + ref_line)
            .properties(
                title=plot_title,
                width=440,
                height=280,
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create perturbation decay chart: {exc}")


def create_resilience_ranking_chart(
    resilience_df: pl.DataFrame,
    title: str = "Predictor Perturbation Resilience Index (PRI)",
    subtitle: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade horizontal bar chart ranking predictors by their resilience score."""
    try:
        if resilience_df.is_empty():
            return Failure("Resilience DataFrame is empty.")

        df_enriched = with_display_labels(resilience_df)
        pred_label_col = "predictor_display" if "predictor_display" in df_enriched.columns else "predictor_name"
        annotated = df_enriched.with_columns([
            pl.format("{}", pl.col("pri_score").round(2)).alias("display_text"),
            pl.col("perturbation_type").str.replace_all("_", " ").str.to_titlecase().alias("perturbation_type"),
        ])

        cohort_ids = resilience_df["cohort_id"].unique().to_list()
        cohort_str = get_cohort_display_label(cohort_ids[0]) if len(cohort_ids) == 1 else f"{len(cohort_ids)} Cohorts"
        sub_title = subtitle or f"Cohort: {cohort_str} | Metric: Normalized Trapezoidal AUC Retention (1.0 = Complete Resilience)"

        plot_title = alt.TitleParams(
            text=title,
            subtitle=sub_title,
            fontSize=13,
            subtitleFontSize=10.5,
            subtitleColor="#555555",
        )

        base = alt.Chart(annotated).encode(
            y=alt.Y(
                f"{pred_label_col}:N",
                title="Predictor / Biomarker",
                sort=alt.EncodingSortField(field="pri_score", op="mean", order="descending"),
                axis=alt.Axis(labelLimit=380, labelFontSize=10),
            ),
            color=alt.Color("perturbation_type:N", title="Perturbation Modality", scale=alt.Scale(scheme="tableau10")),
            xOffset=alt.XOffset("perturbation_type:N"),
        )

        bars = base.mark_bar(height=14, cornerRadiusEnd=3).encode(
            x=alt.X(
                "pri_score:Q",
                title="Perturbation Resilience Index (PRI)",
                scale=alt.Scale(domain=[0.0, 1.05]),
            ),
            tooltip=[
                "cohort_id",
                pred_label_col,
                "perturbation_type",
                alt.Tooltip("baseline_roc_auc:Q", format=".3f", title="Baseline AUC"),
                alt.Tooltip("min_roc_auc:Q", format=".3f", title="Min AUC"),
                alt.Tooltip("pri_score:Q", format=".3f", title="PRI Score"),
                alt.Tooltip("relative_auc_retained:Q", format=".1%", title="Retention"),
            ],
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"x": [1.0]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(x="x:Q")
        )

        chart = (
            (bars + ref_line)
            .properties(
                title=plot_title,
                width=450,
                height=alt.Step(26),
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create resilience ranking chart: {exc}")


def create_pan_cohort_faceted_decay_chart(
    perturbation_df: pl.DataFrame,
    perturbation_type: str,
    metric: str = "roc_auc",
    uncertainty: UncertaintyMode = "shaded",
    title: str | None = None,
    subtitle: str | None = None,
    columns: int = 3,
) -> Result[alt.Chart, str]:
    """Create a publication-grade 3x3 faceted grid showing decay curves across all cohorts.

    Supports 3 uncertainty modes:
    - 'shaded': Semi-transparent confidence area band.
    - 'errorbars': Explicit vertical error whiskers for each point.
    - 'none': Clean curves showing point estimates only.

    Legend symbols strictly reflect the solid, full-opacity colors of the points.
    """
    try:
        df_enriched = with_display_labels(perturbation_df)
        sub_df = df_enriched.filter(pl.col("perturbation_type") == perturbation_type)
        if sub_df.is_empty():
            return Failure(f"No records found for perturbation_type='{perturbation_type}'.")

        pred_label_col = "predictor_display" if "predictor_display" in sub_df.columns else "predictor_name"
        cohort_label_col = "cohort_display" if "cohort_display" in sub_df.columns else "cohort_id"

        is_pr = metric == "pr_auc"
        metric_col = "pr_auc" if is_pr else "roc_auc"
        lower_col = "pr_auc_ci_lower" if is_pr else "roc_auc_ci_lower"
        upper_col = "pr_auc_ci_upper" if is_pr else "roc_auc_ci_upper"

        x_title_map = {
            "jitter": "Multiplicative Expression Jitter (σ)",
            "dropout": "Gene Dropout Rate (Probability)",
            "dilution": "Immune Effector Dilution Factor (α)",
            "label_noise": "Response Label Noise Rate (η)",
        }
        x_title = x_title_map.get(perturbation_type, "Perturbation Intensity")
        p_name = perturbation_type.replace("_", " ").title()

        n_cohorts = sub_df["cohort_id"].n_unique()
        main_title = title or f"Pan-Cohort Performance Decay: {p_name} Perturbation ({n_cohorts} Cohorts)"

        metric_name = "PR-AUC" if is_pr else "ROC-AUC"
        if subtitle is not None:
            sub_title = subtitle
        else:
            match uncertainty:
                case "shaded":
                    sub_title = (
                        f"Panels: Clinical Cohorts | Y-Axis: Empirical {metric_name} "
                        "(Shaded Band = 95% Bootstrap CI; Dashed Line = 0.50 Random Chance)"
                    )
                case "errorbars":
                    sub_title = (
                        f"Panels: Clinical Cohorts | Y-Axis: Empirical {metric_name} "
                        "(Errorbars = 95% Bootstrap CI; Dashed Line = 0.50 Random Chance)"
                    )
                case "none":
                    sub_title = (
                        f"Panels: Clinical Cohorts | Y-Axis: Empirical {metric_name} "
                        "(Point Estimates; Dashed Line = 0.50 Random Chance)"
                    )
                case _ as unreachable:
                    assert_never(unreachable)

        plot_title = alt.TitleParams(
            text=main_title,
            subtitle=sub_title,
            fontSize=14,
            subtitleFontSize=11,
            subtitleColor="#444444",
        )

        color_scale = alt.Scale(scheme="category10")
        pred_legend = alt.Legend(
            title="Predictor",
            symbolType="circle",
            symbolOpacity=1.0,
            symbolSize=65,
            symbolStrokeWidth=0,
            labelLimit=300,
        )

        base = alt.Chart(sub_df).encode(
            x=alt.X("intensity:Q", title=x_title),
            color=alt.Color(
                f"{pred_label_col}:N",
                title="Predictor",
                scale=color_scale,
                legend=pred_legend,
            ),
        )

        lines = base.mark_line(strokeWidth=1.8).encode(
            y=alt.Y(
                f"{metric_col}:Q",
                title="Empirical ROC-AUC" if not is_pr else "Empirical PR-AUC",
                scale=alt.Scale(domain=[0.25, 0.95], zero=False),
            ),
        )

        points = base.mark_point(size=35, filled=True, opacity=1.0).encode(
            y=f"{metric_col}:Q",
            tooltip=[
                "cohort_id",
                cohort_label_col,
                pred_label_col,
                "perturbation_type",
                alt.Tooltip("intensity:Q", format=".2f", title="Intensity"),
                alt.Tooltip(f"{metric_col}:Q", format=".3f", title="Score"),
                alt.Tooltip(f"{lower_col}:Q", format=".3f", title="CI Lower"),
                alt.Tooltip(f"{upper_col}:Q", format=".3f", title="CI Upper"),
            ],
        )

        areas = base.mark_area(opacity=0.18).encode(
            y=f"{lower_col}:Q",
            y2=f"{upper_col}:Q",
        )

        error_bars = base.mark_rule(strokeWidth=1.2, opacity=0.80).encode(
            y=f"{lower_col}:Q",
            y2=f"{upper_col}:Q",
        )

        ref_line = base.mark_rule(
            strokeDash=[4, 4], color="darkgray", strokeWidth=1.2
        ).encode(y=alt.datum(0.50))

        match uncertainty:
            case "shaded":
                core = areas + lines + points
            case "errorbars":
                core = lines + error_bars + points
            case "none":
                core = lines + points
            case _ as unreachable:
                assert_never(unreachable)

        faceted = (
            (core + ref_line)
            .properties(width=220, height=160)
            .facet(
                facet=alt.Facet(f"{cohort_label_col}:N", title=None, header=alt.Header(labelFontSize=10.5, labelFontWeight="bold")),
                columns=columns,
            )
            .properties(title=plot_title)
        )

        return Success(faceted)
    except Exception as exc:
        return Failure(f"Failed to create pan-cohort faceted decay chart: {exc}")


def create_cross_cohort_resilience_heatmap(
    resilience_df: pl.DataFrame,
    title: str | None = None,
    subtitle: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a publication-grade cross-cohort heatmap of Perturbation Resilience Index (PRI)."""
    try:
        if resilience_df.is_empty():
            return Failure("Resilience DataFrame is empty.")

        # Enrich with human-readable display labels
        df_enriched = with_display_labels(resilience_df)
        c_label = "cohort_display" if "cohort_display" in df_enriched.columns else "cohort_id"
        p_label = "predictor_display" if "predictor_display" in df_enriched.columns else "predictor_name"

        # Aggregate across perturbation types per cohort and predictor
        agg_df = (
            df_enriched.group_by([c_label, p_label, "cohort_id", "predictor_name", "category"])
            .agg(
                pl.col("pri_score").mean().alias("mean_pri"),
                pl.col("relative_auc_retained").mean().alias("mean_retention"),
            )
            .with_columns(
                pl.format("{}", pl.col("mean_pri").round(2)).alias("display_text")
            )
        )

        main_title = title or "Cross-Cohort Perturbation Resilience Index (PRI)"
        sub_title = subtitle or "Mean Normalized Area Under Decay Curve (1.0 = Complete Retention across 4 Perturbation Types)"

        plot_title = alt.TitleParams(
            text=main_title,
            subtitle=sub_title,
            fontSize=14,
            subtitleFontSize=11,
            subtitleColor="#444444",
        )

        base = alt.Chart(agg_df).encode(
            x=alt.X(
                f"{c_label}:N",
                title="Clinical Cohort",
                axis=alt.Axis(labelAngle=-40, labelLimit=320, labelFontSize=10),
            ),
            y=alt.Y(
                f"{p_label}:N",
                title="Predictor / Biomarker",
                sort=alt.EncodingSortField(field="mean_pri", op="mean", order="descending"),
                axis=alt.Axis(labelLimit=380, labelFontSize=10),
            ),
        )

        heatmap = base.mark_rect().encode(
            color=alt.Color(
                "mean_pri:Q",
                title="Mean PRI",
                scale=alt.Scale(scheme="viridis", domain=[0.80, 1.02]),
            ),
            tooltip=[
                c_label,
                p_label,
                "cohort_id",
                "predictor_name",
                "category",
                alt.Tooltip("mean_pri:Q", format=".3f", title="Mean PRI"),
                alt.Tooltip("mean_retention:Q", format=".1%", title="Mean Retention"),
            ],
        )

        text = base.mark_text(baseline="middle", fontSize=10, fontWeight="bold").encode(
            text="display_text:N",
            color=alt.condition(
                alt.datum.mean_pri < 0.88,
                alt.value("white"),
                alt.value("black"),
            ),
        )

        chart = (
            (heatmap + text)
            .properties(
                title=plot_title,
                width=alt.Step(48),
                height=alt.Step(28),
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create cross-cohort resilience heatmap: {exc}")


def create_meta_analytic_decay_chart(
    perturbation_df: pl.DataFrame,
    perturbation_type: str,
    metric: str = "roc_auc",
    uncertainty: UncertaintyMode = "shaded",
    title: str | None = None,
    subtitle: str | None = None,
) -> Result[alt.Chart, str]:
    """Create a meta-analytic decay curve averaging empirical ROC-AUC across all cohorts +/- SEM.

    Supports 3 uncertainty modes:
    - 'shaded': Semi-transparent SEM area band.
    - 'errorbars': Explicit vertical error whiskers for each point (+/- 1 SEM).
    - 'none': Clean curves showing mean point estimates only.

    Legend symbols strictly reflect the solid, full-opacity colors of the points.
    """
    try:
        df_enriched = with_display_labels(perturbation_df)
        sub_df = df_enriched.filter(pl.col("perturbation_type") == perturbation_type)
        if sub_df.is_empty():
            return Failure(f"No records found for perturbation_type='{perturbation_type}'.")

        pred_label_col = "predictor_display" if "predictor_display" in sub_df.columns else "predictor_name"
        is_pr = metric == "pr_auc"
        metric_col = "pr_auc" if is_pr else "roc_auc"

        # Compute mean +/- SEM across cohorts per (intensity, predictor_display)
        agg_df = (
            sub_df.group_by(["intensity", pred_label_col])
            .agg(
                pl.col(metric_col).mean().alias("mean_score"),
                pl.col(metric_col).std().alias("std_score"),
                pl.len().alias("n_cohorts"),
            )
            .with_columns(
                (pl.col("std_score") / (pl.col("n_cohorts").sqrt().clip(lower_bound=1))).alias("sem_score"),
            )
            .with_columns(
                (pl.col("mean_score") - pl.col("sem_score")).alias("ci_lo"),
                (pl.col("mean_score") + pl.col("sem_score")).alias("ci_hi"),
            )
        )

        x_title_map = {
            "jitter": "Multiplicative Expression Jitter (σ)",
            "dropout": "Gene Dropout Rate (Probability)",
            "dilution": "Immune Effector Dilution Factor (α)",
            "label_noise": "Response Label Noise Rate (η)",
        }
        x_title = x_title_map.get(perturbation_type, "Perturbation Intensity")
        p_name = perturbation_type.replace("_", " ").title()

        n_c = sub_df["cohort_id"].n_unique()
        main_title = title or f"Meta-Analytic Performance Decay: {p_name} ({n_c} Cohorts)"

        metric_name = "PR-AUC" if is_pr else "ROC-AUC"
        if subtitle is not None:
            sub_title = subtitle
        else:
            match uncertainty:
                case "shaded":
                    sub_title = (
                        f"Mean Empirical {metric_name} across {n_c} Clinical Cohorts "
                        "(Shaded Band = ±1 SEM; Dashed Line = 0.50 Random Chance)"
                    )
                case "errorbars":
                    sub_title = (
                        f"Mean Empirical {metric_name} ± 1 SEM across {n_c} Clinical Cohorts "
                        "(Dashed Line = 0.50 Random Chance)"
                    )
                case "none":
                    sub_title = (
                        f"Mean Empirical {metric_name} across {n_c} Clinical Cohorts "
                        "(Mean Trajectories; Dashed Line = 0.50 Random Chance)"
                    )
                case _ as unreachable:
                    assert_never(unreachable)

        plot_title = alt.TitleParams(
            text=main_title,
            subtitle=sub_title,
            fontSize=13,
            subtitleFontSize=10.5,
            subtitleColor="#444444",
        )

        color_scale = alt.Scale(scheme="category10")
        pred_legend = alt.Legend(
            title="Predictor",
            symbolType="circle",
            symbolOpacity=1.0,
            symbolSize=70,
            symbolStrokeWidth=0,
            labelLimit=300,
        )

        base = alt.Chart(agg_df).encode(
            x=alt.X("intensity:Q", title=x_title),
            color=alt.Color(
                f"{pred_label_col}:N",
                title="Predictor",
                scale=color_scale,
                legend=pred_legend,
            ),
        )

        lines = base.mark_line(strokeWidth=2.2).encode(
            y=alt.Y(
                "mean_score:Q",
                title="Mean Empirical ROC-AUC" if not is_pr else "Mean Empirical PR-AUC",
                scale=alt.Scale(domain=[0.40, 0.75], zero=False),
            ),
        )

        points = base.mark_point(size=50, filled=True, opacity=1.0).encode(
            y="mean_score:Q",
            tooltip=[
                pred_label_col,
                alt.Tooltip("intensity:Q", format=".2f", title="Intensity"),
                alt.Tooltip("mean_score:Q", format=".3f", title="Mean Score"),
                alt.Tooltip("ci_lo:Q", format=".3f", title="-1 SEM"),
                alt.Tooltip("ci_hi:Q", format=".3f", title="+1 SEM"),
                alt.Tooltip("n_cohorts:Q", format="d", title="N Cohorts"),
            ],
        )

        areas = base.mark_area(opacity=0.18).encode(
            y="ci_lo:Q",
            y2="ci_hi:Q",
        )

        error_bars = base.mark_rule(strokeWidth=1.3, opacity=0.85).encode(
            y="ci_lo:Q",
            y2="ci_hi:Q",
        )

        ref_line = (
            alt.Chart(pl.DataFrame({"y": [0.50]}))
            .mark_rule(strokeDash=[4, 4], color="darkgray", strokeWidth=1.5)
            .encode(y="y:Q")
        )

        match uncertainty:
            case "shaded":
                core = areas + lines + points
            case "errorbars":
                core = lines + error_bars + points
            case "none":
                core = lines + points
            case _ as unreachable:
                assert_never(unreachable)

        chart = (
            (core + ref_line)
            .properties(
                title=plot_title,
                width=450,
                height=280,
            )
        )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create meta-analytic decay chart: {exc}")
