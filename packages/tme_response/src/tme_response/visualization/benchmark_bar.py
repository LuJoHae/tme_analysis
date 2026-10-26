"""Altair benchmark comparison bar charts with SVG export."""

from __future__ import annotations

from enum import Enum
from typing import Final, Mapping, Sequence, assert_never

import altair as alt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing
from returns.result import Failure, Result, Success

from tme_response.visualization.display_labels import with_display_labels

CLINICAL_DEFAULT_PREDICTORS: Final[tuple[str, ...]] = (
    "TMB_Genomic",
    "Gene_CD274",
    "Gene_PDCD1",
    "Gene_CD8A",
    "Gene_CTLA4",
    "Gene_CXCL9",
    "Ayers_GEP",
    "Cristescu_GEP",
    "Fehrenbacher_Teff",
    "Ayers_IFNG_6",
    "CYT",
    "TIDE_Inverted",
)

CLINICAL_DISPLAY_OVERRIDES: Final[Mapping[str, str]] = {
    "Deconvolution_Inference_Our": "Deconvolution Inference (Our)",
    "TMB_Genomic": "TMB (Genomic Mutational Burden)",
    "Gene_CD274": "CD274 / PD L1 (Single Gene)",
    "Gene_PDCD1": "PDCD1 / PD 1 (Single Gene)",
    "Gene_CD8A": "CD8A (Single Gene)",
    "Gene_CTLA4": "CTLA4 (Single Gene)",
    "Gene_CXCL9": "CXCL9 (Single Gene)",
    "Ayers_GEP": "Ayers GEP (Ayers et al., 2017)",
    "Cristescu_GEP": "Merck GEP (Cristescu et al., 2018)",
    "Fehrenbacher_Teff": "Teff (Fehrenbacher et al., 2016)",
    "Ayers_IFNG_6": "IFNG 6 (Ayers et al., 2017)",
    "CYT": "CYT (Rooney et al., 2015)",
    "TIDE_Inverted": "Inverted TIDE (Jiang et al., 2018)",
}

CLINICAL_CATEGORY_MAP: Final[Mapping[str, str]] = {
    "Deconvolution_Inference_Our": "Cell Deconvolution (Our)",
    "TMB_Genomic": "Genomic Marker (TMB)",
    "Gene_CD274": "Single Gene / IHC Marker",
    "Gene_PDCD1": "Single Gene / IHC Marker",
    "Gene_CD8A": "Single Gene / IHC Marker",
    "Gene_CTLA4": "Single Gene / IHC Marker",
    "Gene_CXCL9": "Single Gene / IHC Marker",
    "Ayers_GEP": "Clinical Trial Signature",
    "Cristescu_GEP": "Clinical Trial Signature",
    "Fehrenbacher_Teff": "Clinical Trial Signature",
    "Ayers_IFNG_6": "Clinical Trial Signature",
    "CYT": "Clinical Trial Signature",
    "TIDE_Inverted": "Clinical Trial Signature",
    "Composite_Ayers_GEP_TMB": "Multi Omic Composite",
    "Composite_Gene_CXCL9_TMB": "Multi Omic Composite",
    "Gene_CXCL9_AP_Gated": "Multi Omic Composite",
}

CLINICAL_CATEGORY_ORDER: Final[tuple[str, ...]] = (
    "Cell Deconvolution (Our)",
    "Clinical Trial Signature",
    "Single Gene / IHC Marker",
    "Genomic Marker (TMB)",
    "Multi Omic Composite",
)

# Okabe-Ito colorblind-safe palette mapped to clinical categories
CLINICAL_CATEGORY_COLORS: Final[tuple[str, ...]] = (
    "#0072B2",  # Okabe-Ito Blue: Cell Deconvolution (Our)
    "#009E73",  # Okabe-Ito Bluish Green: Clinical Trial Signature
    "#E69F00",  # Okabe-Ito Orange: Single Gene / IHC Marker
    "#D55E00",  # Okabe-Ito Vermillion: Genomic Marker (TMB)
    "#CC79A7",  # Okabe-Ito Reddish Purple: Multi Omic Composite
)


class AggregationMode(str, Enum):
    """Strategy for aggregating classifier ROC-AUC performance across cohorts."""

    COHORT_MEAN = "cohort_mean"
    PANCANCER_POOL = "pancancer_pool"


class ClinicalBarPlotConfig(BaseModel):
    """Immutable configuration for the clinical classifier ROC-AUC horizontal bar chart."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    title: str = "Clinical ICB Response Classifiers: Cross-Cohort ROC-AUC Comparison"
    subtitle: Maybe[str] = Nothing
    x_domain: tuple[float, float] = (0.35, 0.90)
    width: int = 380
    row_step: int = 28
    reference_auc: float = 0.50


def classify_clinical_category(predictor_name: str, fallback_category: str) -> str:
    """Pure mapping from predictor identifier to clinical biomarker category."""
    return CLINICAL_CATEGORY_MAP.get(
        predictor_name,
        "Multi Omic Composite"
        if "composite" in fallback_category.lower()
        else (
            "Genomic Marker (TMB)"
            if "genomic" in fallback_category.lower()
            else (
                "Cell Deconvolution (Our)"
                if "deconv" in fallback_category.lower()
                else "Clinical Trial Signature"
            )
        ),
    )


def filter_stratum_by_mode(df: pl.DataFrame, mode: AggregationMode) -> pl.DataFrame:
    """Filter benchmark rows according to the chosen aggregation mode."""
    base_pre_std = df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Standard")
    )
    match mode:
        case AggregationMode.COHORT_MEAN:
            return base_pre_std.filter(pl.col("pooling_strategy") == "cohort")
        case AggregationMode.PANCANCER_POOL:
            return base_pre_std.filter(
                (pl.col("cohort_id") == "PanCancer-Combined")
                & (pl.col("pooling_strategy") == "standardized")
            )
        case _ as unreachable:
            assert_never(unreachable)


def aggregate_clinical_predictors(
    benchmark_df: pl.DataFrame,
    predictors: Sequence[str] = CLINICAL_DEFAULT_PREDICTORS,
    mode: AggregationMode = AggregationMode.COHORT_MEAN,
) -> Result[pl.DataFrame, str]:
    """Pure Polars aggregation of ROC-AUC and 95% CIs for selected clinical classifiers."""
    required_cols = {
        "predictor_name",
        "category",
        "time_stratum",
        "response_stratum",
        "pooling_strategy",
        "roc_auc",
        "roc_auc_ci_lower",
        "roc_auc_ci_upper",
    }
    missing_cols = required_cols - set(benchmark_df.columns)
    if missing_cols:
        return Failure(f"Missing required columns in benchmark DataFrame: {sorted(missing_cols)}")

    predictor_list = list(predictors)
    if not predictor_list:
        return Failure("Predictor selection list cannot be empty.")

    filtered = filter_stratum_by_mode(benchmark_df, mode).filter(
        pl.col("predictor_name").is_in(predictor_list)
    )
    if filtered.is_empty():
        return Failure(f"No matching records found for mode '{mode.value}' and requested predictors.")

    enriched = with_display_labels(filtered).with_columns([
        pl.when(pl.col("predictor_name").is_in(list(CLINICAL_DISPLAY_OVERRIDES.keys())))
        .then(pl.col("predictor_name").replace(CLINICAL_DISPLAY_OVERRIDES))
        .otherwise(pl.col("predictor_display"))
        .alias("predictor_display"),
        pl.struct(["predictor_name", "category"])
        .map_elements(
            lambda s: classify_clinical_category(
                str(s["predictor_name"]), str(s["category"])
            ),
            return_dtype=pl.String,
        )
        .alias("clinical_category"),
    ])

    summary_df = (
        enriched.group_by(["predictor_name", "predictor_display", "clinical_category"])
        .agg([
            pl.col("roc_auc").mean().alias("mean_roc_auc"),
            pl.col("roc_auc_ci_lower").mean().alias("ci_lower"),
            pl.col("roc_auc_ci_upper").mean().alias("ci_upper"),
            pl.col("roc_auc").std().fill_null(0.0).alias("std_roc_auc"),
            pl.len().alias("n_cohorts"),
        ])
        .with_columns(
            pl.format("{}", pl.col("mean_roc_auc").round(3)).alias("auc_label")
        )
        .sort("mean_roc_auc", descending=True)
    )
    return Success(summary_df)


def create_clinical_selection_bar_chart(
    summary_df: pl.DataFrame,
    config: ClinicalBarPlotConfig = ClinicalBarPlotConfig(),
) -> Result[alt.LayerChart, str]:
    """Construct a Nature Methods horizontal bar chart of ROC-AUC (±95% CI) for clinical classifiers."""
    required = {
        "predictor_display",
        "clinical_category",
        "mean_roc_auc",
        "ci_lower",
        "ci_upper",
        "n_cohorts",
        "auc_label",
    }
    missing = required - set(summary_df.columns)
    if missing:
        return Failure(f"Missing required summary columns: {sorted(missing)}")
    if summary_df.is_empty():
        return Failure("Summary DataFrame is empty.")

    sorted_df = summary_df.sort("mean_roc_auc", descending=True)
    predictor_order = sorted_df["predictor_display"].to_list()
    present_categories = [
        cat
        for cat in CLINICAL_CATEGORY_ORDER
        if cat in set(sorted_df["clinical_category"].to_list())
    ]
    present_colors = [
        CLINICAL_CATEGORY_COLORS[CLINICAL_CATEGORY_ORDER.index(cat)]
        for cat in present_categories
    ]

    x_min, x_max = config.x_domain
    x_scale = alt.Scale(domain=[x_min, x_max], zero=False, clamp=True)
    y_encoding = alt.Y(
        "predictor_display:N",
        title="Clinical Classifier / Biomarker",
        scale=alt.Scale(domain=predictor_order),
        axis=alt.Axis(
            labelLimit=360,
            labelFontSize=11,
            titleFontSize=12,
            titleFontWeight="bold",
            titlePadding=12,
            tickColor="#CBD5E1",
            domainColor="#CBD5E1",
        ),
    )

    base = alt.Chart(sorted_df).encode(y=y_encoding)

    bars = base.mark_bar(
        cornerRadius=0,
        height=18,
        stroke="#1E293B",
        strokeWidth=0.5,
        opacity=0.92,
    ).encode(
        x=alt.X(
            "mean_roc_auc:Q",
            title="Cross-Cohort ROC-AUC (95% CI)",
            scale=x_scale,
            axis=alt.Axis(
                values=[0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90],
                labelFontSize=10.5,
                titleFontSize=12,
                titleFontWeight="bold",
                gridColor="#F1F5F9",
                tickColor="#1E293B",
                domainColor="#CBD5E1",
            ),
        ),
        color=alt.Color(
            "clinical_category:N",
            title="Biomarker Class",
            scale=alt.Scale(domain=present_categories, range=present_colors),
            legend=alt.Legend(
                orient="bottom",
                direction="horizontal",
                columns=2,
                titleFontSize=11,
                titleFontWeight="bold",
                labelFontSize=10.5,
                symbolType="square",
                symbolSize=110,
                labelLimit=260,
            ),
        ),
        tooltip=[
            alt.Tooltip("predictor_display:N", title="Classifier"),
            alt.Tooltip("clinical_category:N", title="Category"),
            alt.Tooltip("mean_roc_auc:Q", format=".4f", title="Mean ROC-AUC"),
            alt.Tooltip("ci_lower:Q", format=".3f", title="95% CI Lower"),
            alt.Tooltip("ci_upper:Q", format=".3f", title="95% CI Upper"),
            alt.Tooltip("n_cohorts:Q", title="Cohorts Evaluated"),
        ],
    )

    error_bars = base.mark_rule(color="#0F172A", strokeWidth=1.5).encode(
        x=alt.X("ci_lower:Q", scale=x_scale),
        x2="ci_upper:Q",
    )

    caps_lower = base.mark_tick(
        color="#0F172A", thickness=1.5, size=8, orient="vertical"
    ).encode(
        x=alt.X("ci_lower:Q", scale=x_scale),
    )

    caps_upper = base.mark_tick(
        color="#0F172A", thickness=1.5, size=8, orient="vertical"
    ).encode(
        x=alt.X("ci_upper:Q", scale=x_scale),
    )

    value_labels = base.mark_text(
        align="left",
        baseline="middle",
        dx=5,
        fontSize=10,
        fontWeight="bold",
        color="#0F172A",
    ).encode(
        x=alt.X("ci_upper:Q", scale=x_scale),
        text="auc_label:N",
    )

    chance_rule = (
        alt.Chart(pl.DataFrame({"x": [config.reference_auc]}))
        .mark_rule(strokeDash=[4, 4], color="#475569", strokeWidth=1.5)
        .encode(x=alt.X("x:Q", scale=x_scale))
    )

    subtitle_text = config.subtitle.value_or(
        "Unweighted mean ROC-AUC ± 95% CI across 9 pre-treatment clinical trials (dashed line = 0.50 random chance)"
    )
    title_params = alt.TitleParams(
        text=config.title,
        subtitle=subtitle_text,
        fontSize=13.5,
        subtitleFontSize=10.5,
        subtitleColor="#475569",
        fontWeight="bold",
        color="#0F172A",
        anchor="start",
        offset=10,
    )

    layered = (
        (bars + error_bars + caps_lower + caps_upper + chance_rule + value_labels)
        .properties(
            title=title_params,
            width=config.width,
            height=alt.Step(config.row_step),
        )
        .configure_view(stroke="#CBD5E1", strokeWidth=0.75, fill="#FFFFFF")
        .configure_axis(labelFont="sans-serif", titleFont="sans-serif")
        .configure_legend(labelFont="sans-serif", titleFont="sans-serif")
    )
    return Success(layered)



def create_benchmark_bar_chart(
    benchmark_df: pl.DataFrame,
    metric: str = "roc_auc",
    title: str = "Cross-Cohort Predictor Performance Comparison",
) -> Result[alt.Chart, str]:
    """Pure functional constructor of grouped Altair bar chart comparing predictor performance."""
    try:
        if metric not in benchmark_df.columns:
            return Failure(f"Metric '{metric}' not found in benchmark DataFrame.")

        df_enriched = with_display_labels(benchmark_df)
        pred_label_col = "predictor_display" if "predictor_display" in df_enriched.columns else "predictor_name"
        cohort_label_col = "cohort_display" if "cohort_display" in df_enriched.columns else "cohort_id"
        metric_title = "ROC-AUC" if metric.lower() == "roc_auc" else "PR-AUC"

        chart = (
            alt.Chart(df_enriched)
            .mark_bar()
            .encode(
                x=alt.X(f"{pred_label_col}:N", title="Predictor / Biomarker", sort="-y", axis=alt.Axis(labelLimit=300, labelAngle=-45)),
                y=alt.Y(f"{metric}:Q", title=metric_title, scale=alt.Scale(domain=[0, 1])),
                color=alt.Color("category:N", title="Category"),
                column=alt.Column(f"{cohort_label_col}:N", title="Validation Cohort"),
                tooltip=[cohort_label_col, pred_label_col, "category", metric],
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
