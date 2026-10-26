#!/usr/bin/env python3
"""update_texh_benchmark.py: Targeted benchmark update for Jiang T-cell Exhaustion.

Re-evaluates Jiang_Texh across all individual and pooled cohorts with positive
orientation (invert_for_response=False), updating:
  1. iatlas_benchmark_metrics.parquet / .csv
  2. cohort_predictability_metrics.parquet / .csv
  3. survival_benchmark_metrics.parquet / .csv
  4. dca_benchmark_metrics.parquet / .csv
  5. All associated publication vector figures (SVG) in output/benchmarks/iatlas_roc_auc/
     and copying them directly to the conversation artifact directory.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Sequence
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from tme_datasets.paths import find_repo_root
from tme_response import (
    PredictionResult,
    PredictorCategory,
    compute_cohort_predictability,
    compute_jiang_texh_score,
    create_auc_heatmap,
    create_benchmark_bar_chart,
    create_c_index_forest_plot,
    create_cohort_predictability_forest_plot,
    create_dca_net_benefit_chart,
    create_summary_forest_plot,
    evaluate_prediction,
    export_chart_svg,
    load_multi_omic_cohort,
    pool_cohort_stratum_data,
    standardize_prediction_scores,
)
from scripts.benchmark_iatlas_icb_predictors import (
    COMBINED_COHORT_SPECS,
    DEFAULT_IATLAS_COHORTS,
    generate_combinatorial_strata,
)
from scripts.benchmark_survival_and_dca import (
    evaluate_cohort_dca,
    evaluate_cohort_survival,
)


def main() -> None:
    repo_root = find_repo_root()
    output_dir = repo_root / "output/benchmarks/iatlas_roc_auc"
    output_dir.mkdir(parents=True, exist_ok=True)
    brain_dir = Path("/Users/halu/.gemini/antigravity/brain/8d054c32-dd8d-492f-aeba-2e31752acbd0")

    print("=== Starting Fast Targeted Update for Jiang_Texh ===")
    cohort_list = list(DEFAULT_IATLAS_COHORTS)

    cohort_cache: dict[str, dict[str, object]] = {}
    new_response_records: list[dict[str, object]] = []
    new_survival_records: list[dict[str, object]] = []
    new_dca_records: list[dict[str, object]] = []

    # 1. Process Individual Cohorts
    for cohort_id in cohort_list:
        print(f"\n--- Processing Cohort: {cohort_id} ---")
        match load_multi_omic_cohort(cohort_id, repo_root=repo_root, filter_baseline=False):
            case Failure(err):
                print(f"  [ERROR] Failed to load {cohort_id}: {err}")
                continue
            case Success(cohort):
                match compute_jiang_texh_score(cohort, invert_for_response=False):
                    case Failure(err):
                        print(f"  [WARN] Failed Texh score for {cohort_id}: {err}")
                        continue
                    case Success(texh_pred):
                        predictors = [texh_pred]

                strata = generate_combinatorial_strata(cohort.clinical_annotations)
                cohort_cache[cohort_id] = {
                    "cancer_type": cohort.cancer_type,
                    "clinical": cohort.clinical_annotations,
                    "strata": strata,
                    "predictors": predictors,
                }

                # Evaluate response discrimination
                for stratum_key, stratum_clinical in strata.items():
                    t_stratum, r_stratum = stratum_key.split("_", 1)
                    match evaluate_prediction(
                        prediction=texh_pred,
                        clinical_df=stratum_clinical,
                        cohort_id=cohort.cohort_id,
                        cancer_type=cohort.cancer_type,
                        time_stratum=t_stratum,
                        response_stratum=r_stratum,
                        pooling_strategy="cohort",
                        n_bootstrap=1000,
                    ):
                        case Success(res):
                            new_response_records.append({
                                "cohort_id": res.cohort_id,
                                "cancer_type": res.cancer_type,
                                "predictor_name": res.predictor_name,
                                "category": res.category.value,
                                "time_stratum": res.time_stratum,
                                "response_stratum": res.response_stratum,
                                "pooling_strategy": res.pooling_strategy,
                                "n_samples": res.n_samples,
                                "n_responders": res.n_responders,
                                "n_non_responders": res.n_non_responders,
                                "roc_auc": res.roc_auc,
                                "roc_auc_ci_lower": res.roc_auc_ci_lower,
                                "roc_auc_ci_upper": res.roc_auc_ci_upper,
                                "p_value_vs_half": res.p_value_vs_half,
                                "pr_auc": res.pr_auc,
                                "pr_auc_ci_lower": res.pr_auc_ci_lower,
                                "pr_auc_ci_upper": res.pr_auc_ci_upper,
                                "baseline_prevalence": res.baseline_prevalence,
                                "delta_pr_auc": res.delta_pr_auc,
                                "p_value_prauc": res.p_value_prauc,
                                "brier_score": res.brier_score,
                            })
                        case Failure(err):
                            pass

                # Evaluate Survival & DCA
                surv_records = evaluate_cohort_survival(
                    predictors=predictors,
                    clinical_df=cohort.clinical_annotations,
                    cohort_id=cohort.cohort_id,
                    cancer_type=cohort.cancer_type,
                    n_bootstrap=500,
                )
                new_survival_records.extend(surv_records)

                dca_records = evaluate_cohort_dca(
                    predictors=predictors,
                    clinical_df=cohort.clinical_annotations,
                    cohort_id=cohort.cohort_id,
                )
                new_dca_records.extend(dca_records)

    # 2. Process Combined Cohort Pools
    print("\n=== Processing Combined Pools ===")
    for pool_id, pool_spec in COMBINED_COHORT_SPECS.items():
        pool_cancer_type = pool_spec["cancer_type"]
        member_cohorts = [cid for cid in pool_spec["cohort_ids"] if cid in cohort_cache]
        if len(member_cohorts) < 2:
            continue

        all_stratum_keys: set[str] = set()
        for cid in member_cohorts:
            all_stratum_keys.update(cohort_cache[cid]["strata"].keys())

        for stratum_key in sorted(all_stratum_keys):
            t_stratum, r_stratum = stratum_key.split("_", 1)
            strata_tuples: list[tuple[str, pl.DataFrame, Sequence[PredictionResult]]] = []
            for cid in member_cohorts:
                c_strata = cohort_cache[cid]["strata"]
                if stratum_key in c_strata:
                    strata_tuples.append((
                        cid,
                        c_strata[stratum_key],
                        cohort_cache[cid]["predictors"],
                    ))

            if len(strata_tuples) < 2:
                continue

            for strategy in ["standardized", "raw"]:
                do_std = strategy == "standardized"
                match pool_cohort_stratum_data(
                    strata_tuples,
                    group_id=pool_id,
                    cancer_type=pool_cancer_type,
                    standardize=do_std,
                ):
                    case Success((pooled_clin, pooled_preds)):
                        for p in pooled_preds:
                            match evaluate_prediction(
                                prediction=p,
                                clinical_df=pooled_clin,
                                cohort_id=pool_id,
                                cancer_type=pool_cancer_type,
                                time_stratum=t_stratum,
                                response_stratum=r_stratum,
                                pooling_strategy=strategy,
                                n_bootstrap=1000,
                            ):
                                case Success(res):
                                    new_response_records.append({
                                        "cohort_id": res.cohort_id,
                                        "cancer_type": res.cancer_type,
                                        "predictor_name": res.predictor_name,
                                        "category": res.category.value,
                                        "time_stratum": res.time_stratum,
                                        "response_stratum": res.response_stratum,
                                        "pooling_strategy": res.pooling_strategy,
                                        "n_samples": res.n_samples,
                                        "n_responders": res.n_responders,
                                        "n_non_responders": res.n_non_responders,
                                        "roc_auc": res.roc_auc,
                                        "roc_auc_ci_lower": res.roc_auc_ci_lower,
                                        "roc_auc_ci_upper": res.roc_auc_ci_upper,
                                        "p_value_vs_half": res.p_value_vs_half,
                                        "pr_auc": res.pr_auc,
                                        "pr_auc_ci_lower": res.pr_auc_ci_lower,
                                        "pr_auc_ci_upper": res.pr_auc_ci_upper,
                                        "baseline_prevalence": res.baseline_prevalence,
                                        "delta_pr_auc": res.delta_pr_auc,
                                        "p_value_prauc": res.p_value_prauc,
                                        "brier_score": res.brier_score,
                                    })
                                case Failure(_):
                                    pass
                    case Failure(_):
                        pass

        # Pool survival data
        pooled_clin_list = []
        for cid in member_cohorts:
            c_df = cohort_cache[cid]["clinical"]
            c_df_clean = c_df.with_columns([
                pl.col(c).cast(pl.String)
                for c, dt in zip(c_df.columns, c_df.dtypes)
                if dt in (pl.Categorical, pl.Object)
            ])
            pooled_clin_list.append(c_df_clean)
        pooled_clin = pl.concat(pooled_clin_list, how="diagonal_relaxed")

        pooled_preds_pieces = []
        cat = PredictorCategory.SIGNATURE
        for cid in member_cohorts:
            for p in cohort_cache[cid]["predictors"]:
                if p.predictor_name == "Jiang_Texh":
                    cat = p.category
                    std_p = standardize_prediction_scores(p)
                    pooled_preds_pieces.append(std_p.predictions)

        if len(pooled_preds_pieces) >= 2:
            pooled_texh_pred = PredictionResult(
                predictor_name="Jiang_Texh",
                category=cat,
                predictions=pl.concat(pooled_preds_pieces, how="vertical"),
            )
            surv_records = evaluate_cohort_survival(
                predictors=[pooled_texh_pred],
                clinical_df=pooled_clin,
                cohort_id=pool_id,
                cancer_type=pool_cancer_type,
                n_bootstrap=500,
            )
            new_survival_records.extend(surv_records)

            dca_records = evaluate_cohort_dca(
                predictors=[pooled_texh_pred],
                clinical_df=pooled_clin,
                cohort_id=pool_id,
            )
            new_dca_records.extend(dca_records)

    print(f"\nGenerated {len(new_response_records)} updated response records for Jiang_Texh.")
    print(f"Generated {len(new_survival_records)} updated survival records for Jiang_Texh.")
    print(f"Generated {len(new_dca_records)} updated DCA records for Jiang_Texh.")

    # 3. Update Response Parquet & CSV
    resp_pq = output_dir / "iatlas_benchmark_metrics.parquet"
    resp_csv = output_dir / "iatlas_benchmark_metrics.csv"
    existing_resp_df = pl.read_parquet(resp_pq)
    filtered_resp_df = existing_resp_df.filter(
        ~pl.col("predictor_name").is_in(["Jiang_Texh_Inverted", "Jiang_Texh"])
    )
    updated_resp_df = pl.concat([filtered_resp_df, pl.DataFrame(new_response_records)])
    updated_resp_df.write_parquet(resp_pq)
    updated_resp_df.write_csv(resp_csv)
    print(f"Saved updated response benchmarks ({len(updated_resp_df)} records) to {resp_pq}")

    # Recompute Predictability
    match compute_cohort_predictability(updated_resp_df):
        case Success(pred_df):
            pred_pq = output_dir / "cohort_predictability_metrics.parquet"
            pred_csv = output_dir / "cohort_predictability_metrics.csv"
            pred_df.write_parquet(pred_pq)
            pred_df.write_csv(pred_csv)
            print(f"Saved updated cohort predictability ({len(pred_df)} rows) to {pred_pq}")
        case Failure(err):
            print(f"Failed predictability: {err}")

    # 4. Update Survival & DCA Parquets
    surv_pq = output_dir / "survival_benchmark_metrics.parquet"
    surv_csv = output_dir / "survival_benchmark_metrics.csv"
    if surv_pq.exists():
        existing_surv_df = pl.read_parquet(surv_pq)
        filtered_surv_df = existing_surv_df.filter(
            ~pl.col("predictor_name").is_in(["Jiang_Texh_Inverted", "Jiang_Texh"])
        )
        updated_surv_df = pl.concat([filtered_surv_df, pl.DataFrame(new_survival_records)])
        updated_surv_df.write_parquet(surv_pq)
        updated_surv_df.write_csv(surv_csv)
        print(f"Saved updated survival benchmarks to {surv_pq}")

    dca_pq = output_dir / "dca_benchmark_metrics.parquet"
    dca_csv = output_dir / "dca_benchmark_metrics.csv"
    if dca_pq.exists():
        existing_dca_df = pl.read_parquet(dca_pq)
        filtered_dca_df = existing_dca_df.filter(
            ~pl.col("strategy").is_in(["Jiang_Texh_Inverted", "Jiang_Texh"])
        )
        updated_dca_df = pl.concat([filtered_dca_df, pl.DataFrame(new_dca_records)])
        updated_dca_df.write_parquet(dca_pq)
        updated_dca_df.write_csv(dca_csv)
        print(f"Saved updated DCA benchmarks to {dca_pq}")

    # 5. Regenerate Heatmaps and Forest Figures
    print("\n=== Regenerating Publication Altair Figures ===")
    primary_display_df = updated_resp_df.filter(
        pl.col("pooling_strategy").is_in(["cohort", "standardized"])
    )
    pre_std_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Standard")
    )
    extreme_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Extreme")
    )

    GENOMIC_MUTATION_PREDICTORS = {
        "TMB_Genomic",
        "Composite_Gene_CXCL9_TMB",
        "Composite_Ayers_GEP_TMB",
        "Gene_CXCL9_AP_Gated",
    }
    pre_std_transcriptomic = pre_std_primary.filter(
        ~pl.col("predictor_name").is_in(GENOMIC_MUTATION_PREDICTORS)
    )
    extreme_transcriptomic = extreme_primary.filter(
        ~pl.col("predictor_name").is_in(GENOMIC_MUTATION_PREDICTORS)
    )

    # 1A. Full Multi-Omic ROC-AUC Heatmap
    match create_auc_heatmap(
        pre_std_primary,
        metric="roc_auc",
        title="Pre-Treatment Standard ICB Response (Multi-Omic, ROC-AUC)",
        subtitle="31 predictors across 9 clinical cohorts and 3 standardized pools (gray denotes missing data)",
    ):
        case Success(chart):
            svg_p = output_dir / "fig_roc_auc_heatmap_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # 1B. Pure Transcriptomic Complete-Matrix ROC-AUC Heatmap
    match create_auc_heatmap(
        pre_std_transcriptomic,
        metric="roc_auc",
        title="Pre-Treatment Standard ICB Response (Transcriptomic Complete Matrix, ROC-AUC)",
        subtitle="27 RNA-seq gene expression predictors evaluated across all 9 cohorts (100% complete data, zero missing values)",
    ):
        case Success(chart):
            svg_p = output_dir / "fig_roc_auc_heatmap_transcriptomic_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # 2. Pure Transcriptomic Extreme ROC-AUC Heatmap
    match create_auc_heatmap(
        extreme_transcriptomic,
        metric="roc_auc",
        title="Pre-Treatment Extreme ICB Response (Transcriptomic Complete Matrix, ROC-AUC)",
        subtitle="27 RNA-seq gene expression predictors evaluated across extreme responders (CR/PR vs PD)",
    ):
        case Success(chart):
            svg_p = output_dir / "fig_roc_auc_heatmap_transcriptomic_pre_extreme.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # 3A. Full Multi-Omic PR-AUC Heatmap
    match create_auc_heatmap(
        pre_std_primary,
        metric="pr_auc",
        title="Pre-Treatment Standard ICB Response (Multi-Omic, PR-AUC)",
        subtitle="31 predictors across 9 clinical cohorts and 3 standardized pools",
    ):
        case Success(chart):
            svg_p = output_dir / "fig_pr_auc_heatmap_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # 3B. Pure Transcriptomic Complete-Matrix PR-AUC Heatmap
    match create_auc_heatmap(
        pre_std_transcriptomic,
        metric="pr_auc",
        title="Pre-Treatment Standard ICB Response (Transcriptomic Complete Matrix, PR-AUC)",
        subtitle="27 RNA-seq gene expression predictors evaluated across all 9 cohorts (100% complete data, zero missing values)",
    ):
        case Success(chart):
            svg_p = output_dir / "fig_pr_auc_heatmap_transcriptomic_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # Summary Forest Plots
    pancan_pre_std = updated_resp_df.filter(
        (pl.col("cohort_id") == "PanCancer-Combined")
        & (pl.col("time_stratum") == "Pre")
        & (pl.col("response_stratum") == "Standard")
        & (pl.col("pooling_strategy") == "standardized")
    )
    match create_summary_forest_plot(pancan_pre_std, metric="roc_auc", title="Pan-Cancer Meta-Analytic ICB Response Discrimination (ROC-AUC)"):
        case Success(chart):
            svg_p = output_dir / "fig_forest_summary_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    match create_summary_forest_plot(pancan_pre_std, metric="pr_auc", title="Pan-Cancer Meta-Analytic ICB Response Precision (PR-AUC)"):
        case Success(chart):
            svg_p = output_dir / "fig_forest_summary_prauc_pre_standard.svg"
            export_chart_svg(chart, svg_p)
            print(f"Saved: {svg_p}")
        case Failure(err):
            print(f"Error: {err}")

    # Cohort Predictability Forest
    pred_pq = output_dir / "cohort_predictability_metrics.parquet"
    if pred_pq.exists():
        predictability_df = pl.read_parquet(pred_pq)
        pre_pred_df = predictability_df.filter(
            (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Standard")
        )
        match create_cohort_predictability_forest_plot(pre_pred_df):
            case Success(chart):
                svg_p = output_dir / "fig_cohort_predictability_forest.svg"
                export_chart_svg(chart, svg_p)
                print(f"Saved: {svg_p}")
            case Failure(err):
                print(f"Error: {err}")

    # Survival Forest & DCA
    if surv_pq.exists():
        surv_df = pl.read_parquet(surv_pq)
        pancan_os = surv_df.filter(
            (pl.col("cohort_id") == "PanCancer-Combined") & (pl.col("endpoint") == "OS")
        )
        if not pancan_os.is_empty():
            match create_c_index_forest_plot(pancan_os, title="Pan-Cancer Overall Survival Discrimination (Harrell's C-index)"):
                case Success(chart):
                    svg_p = output_dir / "fig_survival_c_index_os_forest.svg"
                    export_chart_svg(chart, svg_p)
                    print(f"Saved: {svg_p}")
                case Failure(err):
                    print(f"Error: {err}")

    if dca_pq.exists():
        dca_df = pl.read_parquet(dca_pq)
        pancan_dca = dca_df.filter(pl.col("cohort_id") == "PanCancer-Combined")
        if not pancan_dca.is_empty():
            match create_dca_net_benefit_chart(pancan_dca, title="Pan-Cancer Decision Curve Analysis: Net Clinical Benefit"):
                case Success(chart):
                    svg_p = output_dir / "fig_dca_net_benefit_pan_cohort.svg"
                    export_chart_svg(chart, svg_p)
                    print(f"Saved: {svg_p}")
                case Failure(err):
                    print(f"Error: {err}")

    # 6. Copy to Brain Artifact Directory
    copy_map = {
        "fig_roc_auc_heatmap_pre_standard.svg": "fig_table_s2_roc_auc_heatmap.svg",
        "fig_roc_auc_heatmap_transcriptomic_pre_standard.svg": "fig_table_s2_roc_auc_heatmap_transcriptomic.svg",
        "fig_pr_auc_heatmap_pre_standard.svg": "fig_table_s2_pr_auc_heatmap.svg",
        "fig_pr_auc_heatmap_transcriptomic_pre_standard.svg": "fig_table_s2_pr_auc_heatmap_transcriptomic.svg",
        "fig_forest_summary_pre_standard.svg": "fig_table_s2_forest_summary_roc_auc.svg",
        "fig_forest_summary_prauc_pre_standard.svg": "fig_table_s2_forest_summary_pr_auc.svg",
        "fig_survival_c_index_os_forest.svg": "fig_table_s2_survival_c_index_os.svg",
        "fig_dca_net_benefit_pan_cohort.svg": "fig_table_s2_dca_net_benefit_pan_cohort.svg",
        "fig_cohort_predictability_forest.svg": "fig_cohort_predictability_forest.svg",
    }

    print("\n=== Copying Updated SVGs to Conversation Artifacts ===")
    for src_name, dst_name in copy_map.items():
        src_path = output_dir / src_name
        if src_path.exists():
            dst_path = brain_dir / dst_name
            shutil.copy2(src_path, dst_path)
            print(f"Copied {src_name} -> {dst_path.name}")

    print("\n=== Fast Targeted Update Completed Successfully ===")


if __name__ == "__main__":
    main()
