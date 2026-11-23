#!/usr/bin/env python3
"""3e2_plot_iatlas_icb_predictors.py: Generate publication SVGs from iAtlas predictor benchmarks.

Reads tabular benchmark outputs (.parquet) produced by
`3e1_benchmark_iatlas_icb_predictors.py` and exports all baseline,
predictability, pooling-comparison, and perturbation-resilience vector figures
(SVG) via `tme_response`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import polars as pl
from returns.result import Failure, Success

from tme_response import (
    create_auc_heatmap,
    create_benchmark_bar_chart,
    create_cohort_predictability_forest_plot,
    create_cross_cohort_resilience_heatmap,
    create_meta_analytic_decay_chart,
    create_pan_cohort_faceted_decay_chart,
    create_perturbation_decay_chart,
    create_pooling_comparison_chart,
    create_resilience_ranking_chart,
    create_summary_forest_plot,
    export_chart_svg,
    get_cohort_display_label,
)

GENOMIC_MUTATION_PREDICTORS = {
    "TMB_Genomic",
    "Composite_Gene_CXCL9_TMB",
    "Composite_Ayers_GEP_TMB",
    "Gene_CXCL9_AP_Gated",
}

PERTURBATION_TYPES = ("jitter", "dropout", "dilution", "label_noise")
UNCERTAINTY_MODES = (
    ("shaded", "_shaded"),
    ("none", "_no_uncertainty"),
    ("errorbars", "_errorbars"),
    ("shaded", ""),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot publication SVGs for iAtlas ICB predictor benchmarks."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3e1_benchmark_iatlas_icb_predictors"),
        help="Directory containing benchmark .parquet files from 3e1.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for publication SVGs (defaults to --input-dir).",
    )
    return parser.parse_args()


def _export_svg(
    chart_res,
    svg_path: Path,
    label: str,
    *,
    target_width_mm: float | None = None,
    min_font_size: float | None = None,
) -> None:
    match chart_res:
        case Success(chart):
            kwargs = {}
            if target_width_mm is not None:
                kwargs["target_width_mm"] = target_width_mm
            if min_font_size is not None:
                kwargs["min_font_size"] = min_font_size
            match export_chart_svg(chart, svg_path, **kwargs):
                case Success(p):
                    print(f"Saved: {p}")
                case Failure(err):
                    print(f"Failed to export {label} SVG: {err}")
        case Failure(err):
            print(f"{label} generation failed: {err}")


def plot_baseline_figures(
    results_df: pl.DataFrame,
    predictability_df: pl.DataFrame,
    output_dir: Path,
) -> None:
    """Generate and export baseline discrimination heatmaps, forest plots, and bar charts."""
    print("\n--- Generating Publication Vector Figures (SVG) ---")

    primary_display_df = results_df.filter(
        pl.col("pooling_strategy").is_in(["cohort", "standardized"])
    )
    pre_std_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Standard")
    )
    extreme_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Extreme")
    )
    pre_std_transcriptomic = pre_std_primary.filter(
        ~pl.col("predictor_name").is_in(GENOMIC_MUTATION_PREDICTORS)
    )
    extreme_transcriptomic = extreme_primary.filter(
        ~pl.col("predictor_name").is_in(GENOMIC_MUTATION_PREDICTORS)
    )

    heatmap_specs = [
        (
            pre_std_primary,
            "roc_auc",
            "Pre-Treatment Standard ICB Response (Multi-Omic, ROC-AUC)",
            "All 31 predictors across 9 cohorts + 3 pools (gaps in Gide et al., 2019 indicate unprofiled WES/CNA data)",
            "fig_roc_auc_heatmap_pre_standard.svg",
            "ROC-AUC heatmap",
        ),
        (
            pre_std_transcriptomic,
            "roc_auc",
            "Pre-Treatment Standard ICB Response (Transcriptomic Complete Matrix, ROC-AUC)",
            "27 RNA-seq gene expression predictors evaluated across all 9 cohorts (100% complete data, zero missing values)",
            "fig_roc_auc_heatmap_transcriptomic_pre_standard.svg",
            "Transcriptomic ROC-AUC heatmap",
        ),
        (
            extreme_primary,
            "roc_auc",
            "Pre-Treatment Extreme Responders: CR/PR vs PD (Multi-Omic, ROC-AUC)",
            "All 31 predictors across 9 cohorts + 3 pools (gaps denote unprofiled WES/CNA)",
            "fig_roc_auc_heatmap_pre_extreme.svg",
            "Extreme ROC-AUC heatmap",
        ),
        (
            extreme_transcriptomic,
            "roc_auc",
            "Pre-Treatment Extreme Responders: CR/PR vs PD (Transcriptomic Complete Matrix, ROC-AUC)",
            "27 RNA-seq gene expression predictors (100% complete data, zero missing values)",
            "fig_roc_auc_heatmap_transcriptomic_pre_extreme.svg",
            "Transcriptomic extreme ROC-AUC heatmap",
        ),
        (
            pre_std_primary,
            "pr_auc",
            "Pre-Treatment Standard ICB Response (Multi-Omic, PR-AUC)",
            "All 31 predictors across 9 cohorts + 3 pools (gaps in Gide et al., 2019 indicate unprofiled WES/CNA data)",
            "fig_pr_auc_heatmap_pre_standard.svg",
            "PR-AUC heatmap",
        ),
        (
            pre_std_transcriptomic,
            "pr_auc",
            "Pre-Treatment Standard ICB Response (Transcriptomic Complete Matrix, PR-AUC)",
            "27 RNA-seq gene expression predictors evaluated across all 9 cohorts (100% complete data, zero missing values)",
            "fig_pr_auc_heatmap_transcriptomic_pre_standard.svg",
            "Transcriptomic PR-AUC heatmap",
        ),
        (
            extreme_primary,
            "pr_auc",
            "Pre-Treatment Extreme Responders: CR/PR vs PD (Multi-Omic, PR-AUC)",
            "All 31 predictors across 9 cohorts + 3 pools (gaps denote unprofiled WES/CNA)",
            "fig_pr_auc_heatmap_pre_extreme.svg",
            "Extreme PR-AUC heatmap",
        ),
        (
            extreme_transcriptomic,
            "pr_auc",
            "Pre-Treatment Extreme Responders: CR/PR vs PD (Transcriptomic Complete Matrix, PR-AUC)",
            "27 RNA-seq gene expression predictors (100% complete data, zero missing values)",
            "fig_pr_auc_heatmap_transcriptomic_pre_extreme.svg",
            "Transcriptomic extreme PR-AUC heatmap",
        ),
    ]

    for df_sub, metric, title, subtitle, filename, label in heatmap_specs:
        if not df_sub.is_empty():
            _export_svg(
                create_auc_heatmap(df_sub, metric=metric, title=title, subtitle=subtitle),
                output_dir / filename,
                label,
                target_width_mm=180.0,
                min_font_size=10.0,
            )

    if not pre_std_primary.is_empty():
        _export_svg(
            create_summary_forest_plot(
                pre_std_primary,
                metric="roc_auc",
                title="Cross-Cohort Predictor Performance (Mean ROC-AUC ± 95% CI)",
            ),
            output_dir / "fig_forest_summary_pre_standard.svg",
            "Predictor ROC-AUC forest plot",
        )
        _export_svg(
            create_summary_forest_plot(
                pre_std_primary,
                metric="pr_auc",
                title="Cross-Cohort Predictor Performance (Mean PR-AUC ± 95% CI)",
            ),
            output_dir / "fig_forest_summary_prauc_pre_standard.svg",
            "Predictor PR-AUC forest plot",
        )

    if not predictability_df.is_empty():
        pred_pre_std = predictability_df.filter(
            (pl.col("time_stratum") == "Pre")
            & (pl.col("response_stratum") == "Standard")
            & (pl.col("pooling_strategy").is_in(["cohort", "standardized"]))
        )
        if not pred_pre_std.is_empty():
            _export_svg(
                create_cohort_predictability_forest_plot(
                    pred_pre_std,
                    title="Cohort Predictability Index (Mean RNA Biomarker ROC-AUC ± SD)",
                ),
                output_dir / "fig_cohort_predictability_forest.svg",
                "Cohort predictability forest plot",
            )

    pre_std_all_pools = results_df.filter(
        (pl.col("time_stratum") == "Pre")
        & (pl.col("response_stratum") == "Standard")
        & (pl.col("pooling_strategy").is_in(["raw", "standardized"]))
    )
    if not pre_std_all_pools.is_empty():
        _export_svg(
            create_pooling_comparison_chart(
                pre_std_all_pools,
                title="Impact of Score Standardization on Pooled ROC-AUC",
            ),
            output_dir / "fig_pooling_comparison_pre_standard.svg",
            "Pooling comparison chart",
        )

    if not pre_std_primary.is_empty():
        _export_svg(
            create_benchmark_bar_chart(pre_std_primary, metric="roc_auc"),
            output_dir / "fig_benchmark_bars_pre_standard.svg",
            "Bar chart",
        )


def plot_perturbation_figures(
    pert_benchmark_df: pl.DataFrame,
    resilience_df: pl.DataFrame,
    output_dir: Path,
) -> None:
    """Generate and export perturbation resilience heatmaps, rankings, and decay curves."""
    if not resilience_df.is_empty():
        _export_svg(
            create_cross_cohort_resilience_heatmap(resilience_df),
            output_dir / "fig_perturbation_resilience_cross_cohort_heatmap.svg",
            "Cross-cohort resilience heatmap",
        )

        for cid in pert_benchmark_df["cohort_id"].unique().to_list():
            cid_resil = resilience_df.filter(pl.col("cohort_id") == cid)
            if not cid_resil.is_empty():
                cohort_dir = output_dir / "cohorts" / cid
                cohort_dir.mkdir(parents=True, exist_ok=True)
                cid_disp = get_cohort_display_label(cid)
                chart_res = create_resilience_ranking_chart(
                    cid_resil,
                    title=f"Predictor Perturbation Resilience Index (PRI) - {cid_disp}",
                    subtitle=f"Cohort: {cid_disp} | Metric: Normalized Area Retention (1.0 = Complete Resilience)",
                )
                _export_svg(
                    chart_res,
                    cohort_dir / "fig_perturbation_resilience_ranking.svg",
                    f"Resilience ranking ({cid})",
                )
                if cid == "Gide-iAtlas":
                    _export_svg(
                        chart_res,
                        output_dir / "fig_perturbation_resilience_ranking.svg",
                        f"Resilience ranking ({cid})",
                    )

    for p_type in PERTURBATION_TYPES:
        for u_mode, suffix in UNCERTAINTY_MODES:
            _export_svg(
                create_pan_cohort_faceted_decay_chart(
                    pert_benchmark_df, perturbation_type=p_type, columns=3, uncertainty=u_mode
                ),
                output_dir / f"fig_perturbation_decay_{p_type}_all_cohorts{suffix}.svg",
                f"Pan-cohort faceted decay ({p_type}, {u_mode})",
            )
            _export_svg(
                create_meta_analytic_decay_chart(
                    pert_benchmark_df, perturbation_type=p_type, uncertainty=u_mode
                ),
                output_dir / f"fig_perturbation_decay_pan_cohort_mean_{p_type}{suffix}.svg",
                f"Meta-analytic decay ({p_type}, {u_mode})",
            )

    for cid in pert_benchmark_df["cohort_id"].unique().to_list():
        cid_pert = pert_benchmark_df.filter(pl.col("cohort_id") == cid)
        cohort_dir = output_dir / "cohorts" / cid
        cohort_dir.mkdir(parents=True, exist_ok=True)
        for p_type in PERTURBATION_TYPES:
            for u_mode, suffix in UNCERTAINTY_MODES:
                chart_res = create_perturbation_decay_chart(
                    cid_pert, perturbation_type=p_type, uncertainty=u_mode
                )
                _export_svg(
                    chart_res,
                    cohort_dir / f"fig_perturbation_decay_{p_type}{suffix}.svg",
                    f"Cohort decay ({cid}, {p_type}, {u_mode})",
                )
                if cid == "Gide-iAtlas":
                    _export_svg(
                        chart_res,
                        output_dir / f"fig_perturbation_decay_{p_type}{suffix}.svg",
                        f"Cohort decay ({cid}, {p_type}, {u_mode})",
                    )


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics_pq = input_dir / "iatlas_benchmark_metrics.parquet"
    if not metrics_pq.exists():
        raise FileNotFoundError(f"Missing baseline benchmark metrics: {metrics_pq}")

    results_df = pl.read_parquet(metrics_pq)
    pred_pq = input_dir / "cohort_predictability_metrics.parquet"
    predictability_df = pl.read_parquet(pred_pq) if pred_pq.exists() else pl.DataFrame()

    plot_baseline_figures(results_df, predictability_df, output_dir)

    pert_pq = input_dir / "perturbation_benchmark_metrics.parquet"
    resil_pq = input_dir / "perturbation_resilience_metrics.parquet"
    if pert_pq.exists():
        pert_df = pl.read_parquet(pert_pq)
        resil_df = pl.read_parquet(resil_pq) if resil_pq.exists() else pl.DataFrame()
        plot_perturbation_figures(pert_df, resil_df, output_dir)

    print("\n=== Plotting Finished Successfully ===")


if __name__ == "__main__":
    main()
