"""Unit tests for forest plot visualizations (tme_response.visualization.forest)."""

from __future__ import annotations

import polars as pl
from returns.result import Success

from tme_response.visualization.forest import (
    create_cohort_predictability_forest_plot,
    create_combined_summary_forest_plot,
    create_summary_forest_plot,
)


def _sample_benchmark_df() -> pl.DataFrame:
    return pl.DataFrame({
        "cohort_id": ["PanCancer-Combined", "PanCancer-Combined"],
        "cancer_type": ["Pan-Cancer", "Pan-Cancer"],
        "predictor_name": ["Deconvolution_Inference_Our", "Gene_CD8A"],
        "category": ["Cell Deconvolution", "Single Gene"],
        "roc_auc": [0.7201, 0.6200],
        "roc_auc_ci_lower": [0.6850, 0.5810],
        "roc_auc_ci_upper": [0.7560, 0.6560],
        "pr_auc": [0.5250, 0.4080],
        "pr_auc_ci_lower": [0.4450, 0.3280],
        "pr_auc_ci_upper": [0.6050, 0.4880],
        "baseline_prevalence": [0.314, 0.314],
        "delta_pr_auc": [0.2110, 0.0940],
    })


def _sample_predictability_df() -> pl.DataFrame:
    return pl.DataFrame({
        "cohort_id": ["Gide-iAtlas", "Hugo-iAtlas"],
        "cancer_type": ["Melanoma", "Melanoma"],
        "pooling_strategy": ["cohort", "cohort"],
        "n_samples": [91, 27],
        "n_responders": [49, 14],
        "baseline_prevalence": [0.538, 0.519],
        "mean_roc_auc_rna": [0.7488, 0.5166],
        "std_roc_auc_rna": [0.1146, 0.1041],
        "mean_pr_auc_rna": [0.6200, 0.4500],
        "mean_delta_pr_auc_rna": [0.0820, -0.0690],
        "best_predictor_rna": ["Gene_CD8A", "CYT"],
        "max_roc_auc_rna": [0.8500, 0.6800],
    })


def test_create_summary_forest_plot_roc_auc_defaults() -> None:
    """Verify that create_summary_forest_plot defaults to compact width 260 and trimmed domain [0.44, 0.80]."""
    df = _sample_benchmark_df()
    res = create_summary_forest_plot(df, metric="roc_auc")
    assert isinstance(res, Success)
    chart = res.unwrap()
    chart_dict = chart.to_dict()

    assert chart_dict["width"] == 260
    # Inspect x scale domain in error_bars layer
    layer_0 = chart_dict["layer"][0]
    x_scale = layer_0["encoding"]["x"]["scale"]
    assert x_scale["domain"] == [0.44, 0.80]


def test_create_summary_forest_plot_overrides() -> None:
    """Verify that custom width and domain arguments are respected."""
    df = _sample_benchmark_df()
    res = create_summary_forest_plot(df, metric="roc_auc", width=320, domain=[0.40, 0.85])
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()

    assert chart_dict["width"] == 320
    assert chart_dict["layer"][0]["encoding"]["x"]["scale"]["domain"] == [0.40, 0.85]


def test_create_summary_forest_plot_pr_auc_defaults() -> None:
    """Verify PR-AUC defaults to compact width 260 and trimmed domain [0.25, 0.60]."""
    df = _sample_benchmark_df()
    res = create_summary_forest_plot(df, metric="pr_auc")
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()

    assert chart_dict["width"] == 260
    assert chart_dict["layer"][0]["encoding"]["x"]["scale"]["domain"] == [0.25, 0.60]


def test_create_cohort_predictability_forest_plot_defaults() -> None:
    """Verify predictability forest plot defaults to compact width 280, trimmed domain [0.38, 0.88], and title=None."""
    df = _sample_predictability_df()
    res = create_cohort_predictability_forest_plot(df)
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()

    assert chart_dict["width"] == 280
    layer_0 = chart_dict["layer"][0]
    assert layer_0["encoding"]["x"]["scale"]["domain"] == [0.38, 0.88]
    # Check y-axis title is None (clean, collision-free)
    y_enc = layer_0["encoding"]["y"]
    assert y_enc.get("title") is None


def test_create_cohort_predictability_forest_plot_overrides() -> None:
    """Verify custom width, domain, and show_y_title on predictability plot."""
    df = _sample_predictability_df()
    res = create_cohort_predictability_forest_plot(
        df,
        width=350,
        domain=[0.35, 0.95],
        show_y_title=True,
    )
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()

    assert chart_dict["width"] == 350
    assert chart_dict["layer"][0]["encoding"]["x"]["scale"]["domain"] == [0.35, 0.95]
    assert chart_dict["layer"][0]["encoding"]["y"]["title"] == "Cohort / Clinical Pool"


def test_create_combined_summary_forest_plot_defaults() -> None:
    """Verify combined forest plot defaults: 2 panels, performance ranking (ROC-AUC desc), shared y-axis."""
    df = _sample_benchmark_df()
    res = create_combined_summary_forest_plot(df)
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()

    # Must be an hconcat chart
    assert "hconcat" in chart_dict
    panels = chart_dict["hconcat"]
    assert len(panels) == 2

    # Panel A (ROC-AUC)
    panel_a = panels[0]
    assert panel_a["width"] == 180
    assert "title" not in panel_a
    scale_a = panel_a["layer"][0]["encoding"]["x"]["scale"]
    assert scale_a["domain"] == [0.45, 0.80]
    y_enc_a = panel_a["layer"][0]["encoding"]["y"]
    assert y_enc_a["title"] == "Predictor / Biomarker"
    # Deconvolution Inference (Our) has roc_auc 0.7201 > Gene_CD8A 0.6200, so it must be first
    assert y_enc_a["scale"]["domain"][0] == "Deconvolution Inference (Our)"

    # Panel B (PR-AUC)
    panel_b = panels[1]
    assert panel_b["width"] == 180
    assert "title" not in panel_b
    scale_b = panel_b["layer"][0]["encoding"]["x"]["scale"]
    assert scale_b["domain"] == [0.25, 0.60]
    y_enc_b = panel_b["layer"][0]["encoding"]["y"]
    assert y_enc_b.get("title") is None
    assert y_enc_b["axis"]["labels"] is False
    assert y_enc_b["axis"]["ticks"] is False

    # Resolves and bottom legend
    assert chart_dict["resolve"]["scale"]["y"] == "shared"
    assert chart_dict["resolve"]["legend"]["color"] == "shared"
    legend = panel_a["layer"][0]["encoding"]["color"]["legend"]
    assert legend["orient"] == "bottom"
    assert legend["direction"] == "horizontal"


def test_create_combined_summary_forest_plot_sorting() -> None:
    """Verify sorting parameter works for alphabetical and pr_auc orders."""
    df = _sample_benchmark_df()

    # Alphabetical order: CD8A comes before Deconvolution
    res_alpha = create_combined_summary_forest_plot(df, sort_by="alphabetical")
    assert isinstance(res_alpha, Success)
    domain_alpha = res_alpha.unwrap().to_dict()["hconcat"][0]["layer"][0]["encoding"]["y"]["scale"]["domain"]
    assert domain_alpha[0] == "CD8A (Single Gene)"

    # PR-AUC order: Deconvolution (0.5250) > CD8A (0.4080)
    res_pr = create_combined_summary_forest_plot(df, sort_by="pr_auc")
    assert isinstance(res_pr, Success)
    domain_pr = res_pr.unwrap().to_dict()["hconcat"][0]["layer"][0]["encoding"]["y"]["scale"]["domain"]
    assert domain_pr[0] == "Deconvolution Inference (Our)"


def test_create_combined_summary_forest_plot_custom_domains_and_width() -> None:
    """Verify custom domains and panel width overrides."""
    df = _sample_benchmark_df()
    res = create_combined_summary_forest_plot(
        df,
        panel_width=200,
        roc_domain=[0.40, 0.85],
        pr_domain=[0.20, 0.65],
    )
    assert isinstance(res, Success)
    chart_dict = res.unwrap().to_dict()
    panels = chart_dict["hconcat"]
    assert panels[0]["width"] == 200
    assert panels[0]["layer"][0]["encoding"]["x"]["scale"]["domain"] == [0.40, 0.85]
    assert panels[1]["width"] == 200
    assert panels[1]["layer"][0]["encoding"]["x"]["scale"]["domain"] == [0.20, 0.65]


def test_create_combined_summary_forest_plot_empty_failure() -> None:
    """Verify empty DataFrame returns Failure."""
    empty_df = pl.DataFrame()
    res = create_combined_summary_forest_plot(empty_df)
    assert not isinstance(res, Success)

