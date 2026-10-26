"""Altair publication-grade ROC-AUC heatmap visualizations with SVG export."""

from __future__ import annotations

from typing import Literal, Sequence
import altair as alt
import polars as pl
from returns.result import Failure, Result, Success


from .display_labels import (
    CANONICAL_CANCER_ORDER,
    COHORT_METADATA_MAP,
    get_canonical_cohort_order,
    normalize_cancer_type,
    with_display_labels,
)


def create_cancer_group_track(
    cohort_order: Sequence[str],
    df: pl.DataFrame,
    x_col: str = "cohort_display",
    col_width: float = 68.0,
    track_height: float = 40.0,
) -> alt.Chart:
    """Create a publication-grade spanning cancer type grouping track to place above an AUC heatmap.

    Groups consecutive cohorts sharing the same cancer type (or pooled category) into
    labeled rectangular wireframe boxes aligned with the heatmap columns.
    """
    cohort_to_group: dict[str, str] = {}
    cols_to_select = [
        c for c in [x_col, "cohort_id", "cancer_type", "pooling_strategy"] if c in df.columns
    ]
    sub_df = df.select(cols_to_select).unique()

    for row in sub_df.iter_rows(named=True):
        label = str(row.get(x_col) or row.get("cohort_id") or "")
        cid = str(row.get("cohort_id") or "")
        strat = str(row.get("pooling_strategy") or "")
        ctype = str(row.get("cancer_type") or "")
        meta = COHORT_METADATA_MAP.get(cid, {})

        if strat == "standardized" or "Combined" in label or "Combined" in cid or "prefix" in meta:
            group = "Combined Pools"
        elif ctype:
            group = normalize_cancer_type(ctype)
        elif "cancer" in meta:
            group = normalize_cancer_type(str(meta["cancer"]))
        else:
            group = "Other"

        cohort_to_group[label] = group
        if cid:
            cohort_to_group[cid] = group

    def _lookup_group(cohort_str: str) -> str:
        if cohort_str in cohort_to_group:
            return cohort_to_group[cohort_str]
        if "Combined" in cohort_str:
            return "Combined Pools"
        for cid, meta in COHORT_METADATA_MAP.items():
            if cid in cohort_str or str(meta.get("author", "___")) in cohort_str:
                return normalize_cancer_type(str(meta.get("cancer", "Other")))
        for canonical in CANONICAL_CANCER_ORDER:
            if canonical.lower() in cohort_str.lower():
                return canonical
        return "Other"

    # Find contiguous runs of the same cancer group in cohort_order
    runs: list[tuple[str, int, int]] = []
    for idx, c in enumerate(cohort_order):
        grp = _lookup_group(c)
        if not runs or runs[-1][0] != grp:
            runs.append((grp, idx, idx + 1))
        else:
            prev_grp, start_idx, _ = runs[-1]
            runs[-1] = (prev_grp, start_idx, idx + 1)

    total_width = float(len(cohort_order) * col_width)

    line1_list: list[str] = []
    line2_list: list[str] = []
    x_start_list: list[float] = []
    x_end_list: list[float] = []
    x_mid_list: list[float] = []
    fill_list: list[str] = []

    for i, (grp, start_idx, end_idx) in enumerate(runs):
        xs = float(start_idx * col_width)
        xe = float(end_idx * col_width)
        xm = (xs + xe) / 2.0

        match grp:
            case "Clear Cell Renal Cell Carcinoma" | "Clear Cell Renal Cell":
                l1, l2 = "Clear Cell Renal", "Cell Carcinoma"
            case "Melanoma":
                l1, l2 = "Melanoma", ""
            case "Pancreatic Ductal":
                l1, l2 = "Pancreatic", "Ductal"
            case "Triple Negative Breast":
                l1, l2 = "Triple-Neg.", "Breast"
            case "Urothelial Bladder":
                l1, l2 = "Urothelial", "Bladder"
            case "Combined Pools":
                l1, l2 = "Combined Standardized", "Cohort Pools"
            case _ if " " in grp:
                parts = grp.split(" ", 1)
                l1, l2 = parts[0], parts[1]
            case _:
                l1, l2 = grp, ""

        fill = (
            "#E2E8F0"
            if grp == "Combined Pools"
            else ("#F1F5F9" if i % 2 == 0 else "#FFFFFF")
        )

        line1_list.append(l1)
        line2_list.append(l2)
        x_start_list.append(xs)
        x_end_list.append(xe)
        x_mid_list.append(xm)
        fill_list.append(fill)

    groups_data = pl.DataFrame({
        "line1": line1_list,
        "line2": line2_list,
        "x_start": x_start_list,
        "x_end": x_end_list,
        "x_mid": x_mid_list,
        "fill_color": fill_list,
    })

    base_track = alt.Chart(groups_data).encode(
        x=alt.X("x_start:Q", scale=alt.Scale(domain=[0.0, total_width]), axis=None),
        x2="x_end:Q",
    )
    rect = base_track.mark_rect(stroke="#94A3B8", strokeWidth=0.75).encode(
        color=alt.Color("fill_color:N", scale=None),
    )

    layers: list[alt.Chart] = [rect]

    single_line_data = groups_data.filter(pl.col("line2") == "")
    if not single_line_data.is_empty():
        t_single = alt.Chart(single_line_data).mark_text(
            baseline="middle",
            fontSize=12.0,
            fontWeight="bold",
            color="#0F172A",
            dy=0,
        ).encode(
            x=alt.X("x_mid:Q", scale=alt.Scale(domain=[0.0, total_width]), axis=None),
            text="line1:N",
        )
        layers.append(t_single)

    two_lines_data = groups_data.filter(pl.col("line2") != "")
    if not two_lines_data.is_empty():
        t1 = alt.Chart(two_lines_data).mark_text(
            baseline="middle",
            fontSize=11.0,
            fontWeight="bold",
            color="#0F172A",
            dy=-7.5,
        ).encode(
            x=alt.X("x_mid:Q", scale=alt.Scale(domain=[0.0, total_width]), axis=None),
            text="line1:N",
        )
        t2 = alt.Chart(two_lines_data).mark_text(
            baseline="middle",
            fontSize=11.0,
            fontWeight="bold",
            color="#0F172A",
            dy=7.5,
        ).encode(
            x=alt.X("x_mid:Q", scale=alt.Scale(domain=[0.0, total_width]), axis=None),
            text="line2:N",
        )
        layers.extend([t1, t2])

    track_chart = alt.layer(*layers).properties(
        width=total_width,
        height=track_height,
        title=alt.TitleParams(
            text="Cancer Type Grouping",
            fontSize=12.5,
            fontWeight="bold",
            color="#475569",
            anchor="start",
        ),
    )
    return track_chart


def create_auc_heatmap(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str = "Immunotherapy Response Predictor Performance",
    subtitle: str | None = None,
    cohort_sort: Sequence[str] | Literal["auto"] | None = "auto",
    show_cancer_groups: bool = False,
) -> Result[alt.Chart, str]:
    """Create a publication-grade Altair heatmap displaying discrimination metrics across predictors and cohorts.

    Supports metric="roc_auc", metric="pr_auc", or metric="delta_pr_auc".
    Includes text overlays with formatted metric values and significance asterisks (* p < 0.05, ** p < 0.01).
    Uses human-readable labels without underscores for predictors, and without hyphens for cohorts.
    """
    try:
        if metric not in benchmark_df.columns:
            return Failure(f"Metric '{metric}' not found in benchmark data.")

        # Enrich with human-readable display labels
        df_enriched = with_display_labels(benchmark_df)

        p_col = (
            "p_value_prauc"
            if metric in ("pr_auc", "delta_pr_auc") and "p_value_prauc" in df_enriched.columns
            else "p_value_vs_half"
        )
        has_p = p_col in df_enriched.columns

        # Add label text: e.g. "0.72*"
        if has_p:
            df_annotated = df_enriched.with_columns(
                pl.when(pl.col(p_col) < 0.01)
                .then(pl.format("{}**", pl.col(metric).round(2)))
                .when(pl.col(p_col) < 0.05)
                .then(pl.format("{}*", pl.col(metric).round(2)))
                .otherwise(pl.format("{}", pl.col(metric).round(2)))
                .alias("display_text")
            )
        else:
            df_annotated = df_enriched.with_columns(
                pl.format("{}", pl.col(metric).round(2)).alias("display_text")
            )

        match metric:
            case "roc_auc":
                color_scale = alt.Scale(scheme="redblue", domain=[0.35, 0.85], domainMid=0.50)
                color_title = "ROC-AUC"
                text_color_cond = (alt.datum[metric] > 0.70) | (alt.datum[metric] < 0.40)
            case "pr_auc":
                color_scale = alt.Scale(scheme="redblue", domain=[0.10, 0.85], domainMid=0.35)
                color_title = "PR-AUC"
                text_color_cond = (alt.datum[metric] > 0.65) | (alt.datum[metric] < 0.22)
            case "delta_pr_auc":
                color_scale = alt.Scale(scheme="redblue", domain=[-0.25, 0.35], domainMid=0.0)
                color_title = "ΔPR-AUC"
                text_color_cond = (alt.datum[metric] > 0.18) | (alt.datum[metric] < -0.10)
            case _:
                color_scale = alt.Scale(scheme="blues")
                color_title = metric
                text_color_cond = alt.datum[metric] > 0.5

        x_col = "cohort_display" if "cohort_display" in df_annotated.columns else "cohort_id"
        y_col = "predictor_display" if "predictor_display" in df_annotated.columns else "predictor_name"

        x_sort: list[str] | None
        match cohort_sort:
            case "auto":
                x_sort = get_canonical_cohort_order(df_annotated, x_col=x_col)
            case seq if isinstance(seq, (list, tuple)):
                x_sort = list(seq)
            case _:
                x_sort = None

        base = alt.Chart(df_annotated).encode(
            x=alt.X(
                f"{x_col}:N",
                title="Cohort",
                sort=x_sort,
                axis=alt.Axis(labelAngle=-45, labelLimit=420, labelFontSize=12.0, titleFontSize=13.5),
            ),
            y=alt.Y(
                f"{y_col}:N",
                title="Predictor / Biomarker",
                sort=alt.EncodingSortField(field=metric, op="mean", order="descending"),
                axis=alt.Axis(labelLimit=380, labelFontSize=12.0, titleFontSize=13.5),
            ),
        )

        tooltip_candidates = [
            "predictor_display",
            "cohort_display",
            "predictor_name",
            "cohort_id",
            "cancer_type",
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
        active_tooltips = [c for c in tooltip_candidates if c in df_annotated.columns]

        heatmap = base.mark_rect().encode(
            color=alt.Color(
                f"{metric}:Q",
                title=color_title,
                scale=color_scale,
                legend=alt.Legend(titleFontSize=13.5, labelFontSize=12.0),
            ),
            tooltip=active_tooltips,
        )

        text = base.mark_text(baseline="middle", fontSize=11.5, fontWeight="bold").encode(
            text="display_text:N",
            color=alt.condition(
                text_color_cond,
                alt.value("white"),
                alt.value("black"),
            ),
        )

        title_params = (
            alt.TitleParams(
                text=title,
                subtitle=subtitle,
                fontSize=15.0,
                subtitleFontSize=12.5,
                subtitleColor="#555555",
                anchor="start",
            )
            if subtitle
            else alt.TitleParams(text=title, fontSize=15.0, anchor="start")
        )

        if show_cancer_groups and x_sort:
            main_heat = (heatmap + text).properties(
                width=alt.Step(68),
                height=alt.Step(27),
            )
            track = create_cancer_group_track(
                cohort_order=x_sort,
                df=df_annotated,
                x_col=x_col,
                col_width=68.0,
                track_height=40.0,
            )
            chart = alt.vconcat(track, main_heat, spacing=4).properties(
                title=title_params,
            )
        else:
            chart = (heatmap + text).properties(
                title=title_params,
                width=alt.Step(68),
                height=alt.Step(27),
            )

        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create AUC heatmap: {exc}")

