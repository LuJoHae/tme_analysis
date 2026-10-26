"""Property-based and unit tests for clinical ICB classifier ROC-AUC bar chart (tme_response.visualization.benchmark_bar)."""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st
import polars as pl
from returns.maybe import Some
from returns.result import Failure, Success

from tme_response.visualization.benchmark_bar import (
    CLINICAL_DEFAULT_PREDICTORS,
    AggregationMode,
    ClinicalBarPlotConfig,
    aggregate_clinical_predictors,
    classify_clinical_category,
    create_clinical_selection_bar_chart,
)


def _make_sample_clinical_df() -> pl.DataFrame:
    """Create a synthetic benchmark DataFrame containing TMB and clinical predictors across cohorts and pools."""
    rows: list[dict[str, object]] = []
    cohorts = ("Hugo-iAtlas", "Riaz-iAtlas", "Gide-iAtlas")
    predictors = (
        ("Deconvolution_Inference_Our", "Cell Deconvolution", 0.72, 0.60, 0.84),
        ("TMB_Genomic", "genomic", 0.56, 0.42, 0.70),
        ("Gene_CD274", "signature", 0.58, 0.44, 0.72),
        ("Gene_CD8A", "signature", 0.65, 0.51, 0.79),
        ("Ayers_GEP", "signature", 0.63, 0.49, 0.77),
    )

    for cid in cohorts:
        for p_name, cat, auc_val, ci_l, ci_u in predictors:
            rows.append({
                "cohort_id": cid,
                "cancer_type": "Melanoma",
                "predictor_name": p_name,
                "category": cat,
                "time_stratum": "Pre",
                "response_stratum": "Standard",
                "pooling_strategy": "cohort",
                "roc_auc": auc_val,
                "roc_auc_ci_lower": ci_l,
                "roc_auc_ci_upper": ci_u,
            })

    for p_name, cat, auc_val, ci_l, ci_u in predictors:
        rows.append({
            "cohort_id": "PanCancer-Combined",
            "cancer_type": "Pan-Cancer",
            "predictor_name": p_name,
            "category": cat,
            "time_stratum": "Pre",
            "response_stratum": "Standard",
            "pooling_strategy": "standardized",
            "roc_auc": auc_val + 0.01,
            "roc_auc_ci_lower": ci_l + 0.01,
            "roc_auc_ci_upper": ci_u + 0.01,
        })

    return pl.DataFrame(rows)


def test_default_predictors_contain_tmb_and_clinical_markers() -> None:
    """Ensure CLINICAL_DEFAULT_PREDICTORS includes TMB and clinical IHC/signature markers and excludes our classifier."""
    assert "TMB_Genomic" in CLINICAL_DEFAULT_PREDICTORS
    assert "Deconvolution_Inference_Our" not in CLINICAL_DEFAULT_PREDICTORS
    assert "Gene_CD274" in CLINICAL_DEFAULT_PREDICTORS
    assert "Gene_PDCD1" in CLINICAL_DEFAULT_PREDICTORS
    assert "Gene_CD8A" in CLINICAL_DEFAULT_PREDICTORS
    assert "Ayers_GEP" in CLINICAL_DEFAULT_PREDICTORS
    assert "Fehrenbacher_Teff" in CLINICAL_DEFAULT_PREDICTORS


def test_aggregate_clinical_predictors_cohort_mean_and_pancancer() -> None:
    """Verify both COHORT_MEAN and PANCANCER_POOL aggregation modes produce sorted summaries with zero underscores."""
    df = _make_sample_clinical_df()

    res_cohort = aggregate_clinical_predictors(df, mode=AggregationMode.COHORT_MEAN)
    assert isinstance(res_cohort, Success)
    summary_cohort = res_cohort.unwrap()
    assert summary_cohort.height == 4
    assert "TMB_Genomic" in summary_cohort["predictor_name"].to_list()
    assert "Deconvolution_Inference_Our" not in summary_cohort["predictor_name"].to_list()
    assert summary_cohort["n_cohorts"].to_list() == [3, 3, 3, 3]

    # Verify descending sort by mean_roc_auc
    aucs = summary_cohort["mean_roc_auc"].to_list()
    assert aucs == sorted(aucs, reverse=True)

    # Verify zero underscores in display names
    for label in summary_cohort["predictor_display"].to_list():
        assert "_" not in str(label)

    # Verify PANCANCER_POOL mode
    res_pool = aggregate_clinical_predictors(df, mode=AggregationMode.PANCANCER_POOL)
    assert isinstance(res_pool, Success)
    summary_pool = res_pool.unwrap()
    assert summary_pool.height == 4
    assert summary_pool["n_cohorts"].to_list() == [1, 1, 1, 1]


def test_create_clinical_selection_bar_chart_spec() -> None:
    """Verify Altair chart construction produces a 6-layer horizontal bar chart with 0.50 reference line."""
    df = _make_sample_clinical_df()
    summary_df = aggregate_clinical_predictors(df).unwrap()

    config = ClinicalBarPlotConfig(
        title="Custom Clinical Bar Title",
        subtitle=Some("Custom Subtitle"),
        x_domain=(0.35, 0.90),
        width=400,
    )
    res_chart = create_clinical_selection_bar_chart(summary_df, config=config)
    assert isinstance(res_chart, Success)
    chart_dict = res_chart.unwrap().to_dict()

    assert chart_dict["width"] == 400
    assert chart_dict["title"]["text"] == "Custom Clinical Bar Title"
    assert chart_dict["title"]["subtitle"] == "Custom Subtitle"
    assert len(chart_dict["layer"]) == 6

    # First layer is the bar mark with cornerRadius=0
    bar_layer = chart_dict["layer"][0]
    assert bar_layer["mark"]["type"] == "bar"
    assert bar_layer["mark"]["cornerRadius"] == 0
    assert bar_layer["encoding"]["x"]["scale"]["domain"] == [0.35, 0.90]


def test_aggregate_clinical_predictors_failures() -> None:
    """Verify Failure monad is returned on missing columns or empty predictor lists."""
    df = _make_sample_clinical_df()
    assert isinstance(aggregate_clinical_predictors(df, predictors=()), Failure)
    assert isinstance(aggregate_clinical_predictors(pl.DataFrame({"a": [1]})), Failure)
    assert isinstance(create_clinical_selection_bar_chart(pl.DataFrame()), Failure)


@given(
    auc_tmb=st.floats(min_value=0.35, max_value=0.85, allow_nan=False, allow_infinity=False),
    auc_gep=st.floats(min_value=0.35, max_value=0.85, allow_nan=False, allow_infinity=False),
    delta_ci=st.floats(min_value=0.01, max_value=0.10, allow_nan=False, allow_infinity=False),
)
@settings(max_examples=25)
def test_property_aggregation_invariants(
    auc_tmb: float,
    auc_gep: float,
    delta_ci: float,
) -> None:
    """Property: Aggregated mean_roc_auc stays within [ci_lower, ci_upper] and is sorted descending."""
    df = pl.DataFrame([
        {
            "cohort_id": "Hugo-iAtlas",
            "cancer_type": "Melanoma",
            "predictor_name": "TMB_Genomic",
            "category": "genomic",
            "time_stratum": "Pre",
            "response_stratum": "Standard",
            "pooling_strategy": "cohort",
            "roc_auc": auc_tmb,
            "roc_auc_ci_lower": max(0.0, auc_tmb - delta_ci),
            "roc_auc_ci_upper": min(1.0, auc_tmb + delta_ci),
        },
        {
            "cohort_id": "Hugo-iAtlas",
            "cancer_type": "Melanoma",
            "predictor_name": "Ayers_GEP",
            "category": "signature",
            "time_stratum": "Pre",
            "response_stratum": "Standard",
            "pooling_strategy": "cohort",
            "roc_auc": auc_gep,
            "roc_auc_ci_lower": max(0.0, auc_gep - delta_ci),
            "roc_auc_ci_upper": min(1.0, auc_gep + delta_ci),
        },
    ])

    res = aggregate_clinical_predictors(df)
    assert isinstance(res, Success)
    out = res.unwrap()
    assert out.height == 2

    means = out["mean_roc_auc"].to_list()
    lowers = out["ci_lower"].to_list()
    uppers = out["ci_upper"].to_list()

    assert means[0] >= means[1] - 1e-12
    for m, lo, hi in zip(means, lowers, uppers):
        assert lo <= m + 1e-12
        assert m <= hi + 1e-12
    assert classify_clinical_category("TMB_Genomic", "genomic") == "Genomic Marker (TMB)"
