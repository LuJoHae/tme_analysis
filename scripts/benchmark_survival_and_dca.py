#!/usr/bin/env python3
"""benchmark_survival_and_dca.py: Time-to-Event Survival and Decision Curve Analysis Benchmark.

Evaluates Harrell's C-index (with bootstrap 95% CIs) and Univariable Cox Proportional Hazards
models (Hazard Ratios per 1-SD increase) for Overall Survival (OS) and Progression-Free Survival (PFS)
across all cBioPortal iAtlas cohorts and pooled cancer types.
Computes Decision Curve Analysis (DCA) Net Benefit and Avoidable Interventions across clinical threshold
probabilities for individual cohorts and pan-cohort averages.
Exports tabular summaries (.parquet, .csv) and publication-grade Altair vector SVGs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence
import altair as alt
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from tme_datasets.paths import find_repo_root
from tme_response import (
    PredictionResult,
    PredictorCategory,
    calculate_c_index_bootstrap,
    calculate_decision_curve,
    compute_ayers_ifng6_score,
    compute_cd8_duo_score,
    compute_cd8_infiltrate_predictor,
    compute_cristescu_gep_score,
    compute_cyt_score,
    compute_davoli_cis_score,
    compute_dna_rna_composite,
    compute_fehrenbacher_teff_score,
    compute_freeman_pgm_score,
    compute_genebio_target_score,
    compute_gep_score,
    compute_huang_nrs_score,
    compute_impres_score,
    compute_ipres_score,
    compute_jiang_ctls_score,
    compute_jiang_tams_score,
    compute_jiang_texh_score,
    compute_kong_netbio_score,
    compute_messina_cks_score,
    compute_nurmik_cafs_score,
    compute_roh_is_score,
    compute_single_gene_score,
    compute_tide_score,
    compute_wu_mias_score,
    create_c_index_forest_plot,
    create_cox_hr_forest_plot,
    create_dca_net_benefit_chart,
    export_chart_svg,
    fit_univariable_cox,
    load_multi_omic_cohort,
    standardize_prediction_scores,
)

DEFAULT_IATLAS_COHORTS = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "McDermott-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "Choueiri-iAtlas",
)

COMBINED_COHORT_SPECS = {
    "Melanoma-Combined": {
        "cancer_type": "Melanoma",
        "cohort_ids": ["Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas"],
    },
    "RenalCell-Combined": {
        "cancer_type": "Renal Cell Car RCC",
        "cohort_ids": ["McDermott-iAtlas", "Choueiri-iAtlas"],
    },
    "PanCancer-Combined": {
        "cancer_type": "Pan-Cancer",
        "cohort_ids": list(DEFAULT_IATLAS_COHORTS),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark Survival (C-index, Cox HR) and DCA on iAtlas cohorts."
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default=",".join(DEFAULT_IATLAS_COHORTS),
        help="Comma-separated list of iAtlas cohort IDs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc"),
        help="Output directory for tabular results and publication SVGs.",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=500,
        help="Number of bootstrap resamples for C-index 95%% CI (default: 500).",
    )
    parser.add_argument(
        "--enable-combined-pools",
        action="store_true",
        default=True,
        help="Enable combined cancer-type and pan-cancer cohort pools (default: True).",
    )
    return parser.parse_args()


def compute_standard_predictors(cohort) -> list[PredictionResult]:
    """Compute standard ICB predictors for a given cohort."""
    predictors: list[PredictionResult] = []

    # 1. Single-Gene Predictors
    for gene in ["CXCL9", "CD8A", "PDCD1", "CD274", "CTLA4"]:
        match compute_single_gene_score(cohort, gene):
            case Success(pred):
                predictors.append(pred)
            case Failure(_):
                pass

    # 2. Curated Table S2 Baseline Signatures
    for fn in [
        compute_cyt_score,
        compute_impres_score,
        compute_gep_score,
        lambda c: compute_ipres_score(c, invert_for_response=True),
        compute_genebio_target_score,
        compute_cd8_duo_score,
        compute_davoli_cis_score,
        compute_fehrenbacher_teff_score,
        compute_freeman_pgm_score,
        compute_huang_nrs_score,
        compute_ayers_ifng6_score,
        compute_jiang_ctls_score,
        lambda c: compute_jiang_tams_score(c, invert_for_response=True),
        lambda c: compute_jiang_texh_score(c, invert_for_response=True),
        compute_messina_cks_score,
        lambda c: compute_nurmik_cafs_score(c, invert_for_response=True),
        compute_roh_is_score,
        compute_wu_mias_score,
        compute_cristescu_gep_score,
        compute_kong_netbio_score,
    ]:
        match fn(cohort):
            case Success(pred):
                predictors.append(pred)
            case Failure(_):
                pass

    match compute_tide_score(cohort, invert_for_response=True):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_cd8_infiltrate_predictor(cohort):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    if isinstance(cohort.tmb_scores, Some):
        tmb_df = cohort.tmb_scores.unwrap()
        predictors.append(
            PredictionResult(
                predictor_name="TMB_Genomic",
                category=PredictorCategory.GENOMIC,
                predictions=tmb_df.select([
                    pl.col("sample_id"),
                    pl.col("tmb_per_mb").alias("score"),
                ]),
            )
        )
        match compute_single_gene_score(cohort, "CXCL9"):
            case Success(cxcl9_pred):
                match compute_dna_rna_composite(cohort, cxcl9_pred):
                    case Success(comp_pred):
                        predictors.append(comp_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    return predictors


def evaluate_cohort_survival(
    predictors: Sequence[PredictionResult],
    clinical_df: pl.DataFrame,
    cohort_id: str,
    cancer_type: str,
    n_bootstrap: int = 500,
) -> list[dict[str, object]]:
    """Evaluate Harrell's C-index and Univariable Cox Models for OS and PFS."""
    records: list[dict[str, object]] = []

    endpoints = [
        ("OS", "os_months", "os_event"),
        ("PFS", "pfs_months", "pfs_event"),
    ]

    for ep_name, time_col, event_col in endpoints:
        if time_col not in clinical_df.columns or event_col not in clinical_df.columns:
            continue

        clin_sub = clinical_df.filter(
            pl.col(time_col).is_not_null()
            & pl.col(time_col).is_finite()
            & (pl.col(time_col) > 0)
            & pl.col(event_col).is_not_null()
            & pl.col(event_col).is_finite()
        )

        n_samples = len(clin_sub)
        if n_samples < 5:
            continue

        n_events = int(clin_sub[event_col].cast(pl.Int64).sum())
        if n_events < 2:
            continue


        for pred in predictors:
            merged = clin_sub.join(pred.predictions, on="sample_id", how="inner").filter(
                pl.col("score").is_not_null() & pl.col("score").is_finite()
            )
            if len(merged) < 5:
                continue

            times = merged[time_col].to_numpy()
            events = merged[event_col].to_numpy()
            scores = merged["score"].to_numpy()

            # 1. C-index bootstrap
            c_res = calculate_c_index_bootstrap(
                times, events, scores, n_bootstrap=n_bootstrap, seed=42
            )
            # 2. Cox Proportional Hazards
            cox_res = fit_univariable_cox(times, events, scores, standardize_score=True)

            c_val, ci_l, ci_u, p_c = np.nan, np.nan, np.nan, np.nan
            if isinstance(c_res, Success):
                ci_model = c_res.unwrap()
                c_val = ci_model.c_index
                ci_l = ci_model.ci_lower
                ci_u = ci_model.ci_upper
                p_c = ci_model.p_value_vs_half

            hr_val, hr_l, hr_u, coef, bse, p_cox = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            if isinstance(cox_res, Success):
                cox_model = cox_res.unwrap()
                hr_val = cox_model.hazard_ratio
                hr_l = cox_model.hr_ci_lower
                hr_u = cox_model.hr_ci_upper
                coef = cox_model.coefficient
                bse = cox_model.se_coefficient
                p_cox = cox_model.p_value

            if np.isfinite(c_val) or np.isfinite(hr_val):
                records.append({
                    "cohort_id": cohort_id,
                    "cancer_type": cancer_type,
                    "predictor_name": pred.predictor_name,
                    "category": pred.category.value,
                    "endpoint": ep_name,
                    "n_samples": len(merged),
                    "n_events": int(np.sum(events)),
                    "c_index": c_val,
                    "ci_lower": ci_l,
                    "ci_upper": ci_u,
                    "p_value_c_index": p_c,
                    "hazard_ratio": hr_val,
                    "hr_ci_lower": hr_l,
                    "hr_ci_upper": hr_u,
                    "coefficient": coef,
                    "se_coefficient": bse,
                    "p_value_cox": p_cox,
                })

    return records


def evaluate_cohort_dca(
    predictors: Sequence[PredictionResult],
    clinical_df: pl.DataFrame,
    cohort_id: str,
    thresholds: Sequence[float] | None = None,
) -> list[dict[str, object]]:
    """Evaluate Decision Curve Analysis (DCA) for standard Pre-treatment binary response."""
    records: list[dict[str, object]] = []

    if "response_binary" not in clinical_df.columns:
        return records

    # Keep baseline / pre-treatment samples
    has_tp = "biopsy_timepoint" in clinical_df.columns
    resp_sub = clinical_df.filter(
        pl.col("response_binary").is_not_null()
        & pl.col("response_binary").is_finite()
        & (
            pl.col("biopsy_timepoint").is_in(["Pre", "Unknown"])
            if has_tp
            else pl.lit(True)
        )
    )



    n_samples = len(resp_sub)
    if n_samples < 10:
        return records

    n_resp = int(resp_sub["response_binary"].cast(pl.Int64).sum())
    n_non = n_samples - n_resp
    if n_resp < 3 or n_non < 3:
        return records

    # Evaluate each predictor
    seen_baselines = False
    for pred in predictors:
        merged = resp_sub.join(pred.predictions, on="sample_id", how="inner").filter(
            pl.col("score").is_not_null() & pl.col("score").is_finite()
        )
        if len(merged) < 10:
            continue

        match calculate_decision_curve(
            merged["response_binary"].to_numpy(),
            merged["score"].to_numpy(),
            predictor_name=pred.predictor_name,
            thresholds=thresholds,
            calibration="logistic",
        ):
            case Success(dca_df):
                for row in dca_df.to_dicts():
                    strat = str(row["strategy"])
                    # Avoid duplicating Treat All / Treat None multiple times per cohort
                    if strat in ("Treat All", "Treat None"):
                        if not seen_baselines:
                            row["cohort_id"] = cohort_id
                            records.append(row)
                    else:
                        row["cohort_id"] = cohort_id
                        records.append(row)
                seen_baselines = True
            case Failure(_):
                pass

    return records


def run_survival_and_dca_benchmark(
    cohorts: Sequence[str] = DEFAULT_IATLAS_COHORTS,
    output_dir: Path = Path("output/benchmarks/iatlas_roc_auc"),
    n_bootstrap: int = 500,
    enable_combined_pools: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Execute end-to-end survival (OS/PFS) and Decision Curve Analysis across cohorts."""
    repo_root = find_repo_root()
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_list = list(cohorts)
    print(f"=== Starting Survival & Decision Curve Analysis Benchmark ===")
    print(f"Target Cohorts ({len(cohort_list)}): {', '.join(cohort_list)}")
    print(f"Output Directory: {output_dir}")
    print(f"Bootstrap Resamples: {n_bootstrap}")

    all_survival_records: list[dict[str, object]] = []
    all_dca_records: list[dict[str, object]] = []


    # Cache for cohort pooling
    cohort_cache: dict[str, dict[str, object]] = {}

    # 1. Process Individual Cohorts
    for cohort_id in cohort_list:
        print(f"\n--- Loading Cohort: {cohort_id} ---")
        match load_multi_omic_cohort(cohort_id, repo_root=repo_root, filter_baseline=False):
            case Failure(err):
                print(f"  [SKIPPED] Failed to load {cohort_id}: {err}")
                continue
            case Success(cohort):
                predictors = compute_standard_predictors(cohort)
                print(f"  Calculated {len(predictors)} predictors.")

                # Save into cache
                cohort_cache[cohort_id] = {
                    "cancer_type": cohort.cancer_type,
                    "clinical": cohort.clinical_annotations,
                    "predictors": predictors,
                }

                # Evaluate Survival
                surv_records = evaluate_cohort_survival(
                    predictors=predictors,
                    clinical_df=cohort.clinical_annotations,
                    cohort_id=cohort.cohort_id,
                    cancer_type=cohort.cancer_type,
                    n_bootstrap=n_bootstrap,
                )
                print(f"  Survival evaluations: {len(surv_records)} records.")
                all_survival_records.extend(surv_records)

                # Evaluate DCA
                dca_records = evaluate_cohort_dca(
                    predictors=predictors,
                    clinical_df=cohort.clinical_annotations,
                    cohort_id=cohort.cohort_id,
                )
                print(f"  DCA evaluations: {len(dca_records)} points.")
                all_dca_records.extend(dca_records)

    # 2. Process Combined Cohort Pools with Standardized Scores
    if enable_combined_pools:
        print("\n=== Constructing Combined Cohort Pools (Standardized Scores) ===")
        for pool_name, pool_spec in COMBINED_COHORT_SPECS.items():

            pool_cohort_ids = [cid for cid in pool_spec["cohort_ids"] if cid in cohort_cache]
            if len(pool_cohort_ids) < 2:
                continue

            cancer_type = pool_spec["cancer_type"]
            print(f"\n--- Pool: {pool_name} ({len(pool_cohort_ids)} cohorts: {', '.join(pool_cohort_ids)}) ---")

            # Combine clinical annotations
            pooled_clinical_list = []
            for cid in pool_cohort_ids:
                c_df = cohort_cache[cid]["clinical"]
                c_df_clean = c_df.with_columns([
                    pl.col(c).cast(pl.String)
                    for c, dt in zip(c_df.columns, c_df.dtypes)
                    if dt in (pl.Categorical, pl.Object)
                ])
                pooled_clinical_list.append(c_df_clean)
            pooled_clinical = pl.concat(pooled_clinical_list, how="diagonal_relaxed")


            # Standardize and pool predictors
            pooled_predictors: list[PredictionResult] = []
            all_pred_names = {
                p.predictor_name for cid in pool_cohort_ids for p in cohort_cache[cid]["predictors"]
            }

            for p_name in sorted(all_pred_names):
                pred_pieces = []
                cat = PredictorCategory.SIGNATURE
                for cid in pool_cohort_ids:
                    matching = [p for p in cohort_cache[cid]["predictors"] if p.predictor_name == p_name]
                    if matching:
                        cat = matching[0].category
                        std_p = standardize_prediction_scores(matching[0])
                        pred_pieces.append(std_p.predictions)

                if len(pred_pieces) >= 2:
                    pooled_preds_df = pl.concat(pred_pieces, how="vertical")
                    pooled_predictors.append(
                        PredictionResult(
                            predictor_name=p_name,
                            category=cat,
                            predictions=pooled_preds_df,
                        )
                    )

            # Evaluate pooled survival
            surv_records = evaluate_cohort_survival(
                predictors=pooled_predictors,
                clinical_df=pooled_clinical,
                cohort_id=pool_name,
                cancer_type=cancer_type,
                n_bootstrap=n_bootstrap,
            )

            print(f"  Pooled Survival evaluations: {len(surv_records)} records.")
            all_survival_records.extend(surv_records)

            # Evaluate pooled DCA
            dca_records = evaluate_cohort_dca(
                predictors=pooled_predictors,
                clinical_df=pooled_clinical,
                cohort_id=pool_name,
            )
            print(f"  Pooled DCA evaluations: {len(dca_records)} points.")
            all_dca_records.extend(dca_records)

    # 3. Export Tabular Datasets
    survival_df = pl.DataFrame(all_survival_records)
    dca_df = pl.DataFrame(all_dca_records)

    if not survival_df.is_empty():
        surv_parquet = output_dir / "survival_benchmark_metrics.parquet"
        surv_csv = output_dir / "survival_benchmark_metrics.csv"
        survival_df.write_parquet(surv_parquet)
        survival_df.write_csv(surv_csv)
        print(f"\nSaved Survival benchmark dataset: {surv_parquet} ({len(survival_df)} records)")

    if not dca_df.is_empty():
        dca_parquet = output_dir / "dca_benchmark_metrics.parquet"
        dca_csv = output_dir / "dca_benchmark_metrics.csv"
        dca_df.write_parquet(dca_parquet)
        dca_df.write_csv(dca_csv)
        print(f"Saved DCA benchmark dataset: {dca_parquet} ({len(dca_df)} records)")

    # 4. Generate Publication Vector Visualizations (Altair -> SVG)
    print("\n=== Generating Publication Vector Graphics (SVGs) ===")

    # A. Survival Forest Plots (OS and PFS)
    if not survival_df.is_empty():
        primary_cohorts_surv = survival_df.filter(
            ~pl.col("cohort_id").is_in(["Melanoma-Combined", "RenalCell-Combined", "PanCancer-Combined"])
        )

        for ep in ["OS", "PFS"]:
            # C-index forest
            match create_c_index_forest_plot(
                primary_cohorts_surv,
                endpoint=ep,
                title=f"Cross-Cohort Predictor Discrimination for {ep} (Harrell's C-index ± 95% CI)",
            ):
                case Success(chart):
                    svg_p = output_dir / f"fig_survival_c_index_{ep.lower()}_forest.svg"
                    export_chart_svg(chart, svg_p)
                    print(f"Saved: {svg_p}")
                case Failure(err):
                    print(f"Failed C-index forest for {ep}: {err}")

            # Cox HR forest
            match create_cox_hr_forest_plot(
                primary_cohorts_surv,
                endpoint=ep,
                title=f"Hazard Ratio per 1-SD Increase for {ep} (Univariable Cox Model)",
            ):
                case Success(chart):
                    svg_p = output_dir / f"fig_survival_cox_hr_{ep.lower()}_forest.svg"
                    export_chart_svg(chart, svg_p)
                    print(f"Saved: {svg_p}")
                case Failure(err):
                    print(f"Failed Cox HR forest for {ep}: {err}")

    # B. Decision Curve Analysis (DCA) Charts
    if not dca_df.is_empty():
        # Pan-cohort mean DCA: average net benefit and avoided interventions per threshold & strategy
        pan_dca = (
            dca_df.filter(
                ~pl.col("cohort_id").is_in(["Melanoma-Combined", "RenalCell-Combined", "PanCancer-Combined"])
            )
            .group_by(["threshold", "strategy"])
            .agg([
                pl.col("net_benefit").mean().alias("net_benefit"),
                pl.col("interventions_avoided_per_100").mean().alias("interventions_avoided_per_100"),
                pl.len().alias("n_cohorts"),
            ])
            .sort("threshold")
        )

        match create_dca_net_benefit_chart(
            pan_dca,
            title="Pan-Cohort Decision Curve Analysis (Mean Net Clinical Benefit)",
        ):
            case Success(chart):
                svg_p = output_dir / "fig_dca_net_benefit_pan_cohort.svg"
                export_chart_svg(chart, svg_p)
                print(f"Saved: {svg_p}")
            case Failure(err):
                print(f"Failed Pan-Cohort DCA chart: {err}")

        # Individual Cohort DCA Charts
        for cid in dca_df["cohort_id"].unique().to_list():
            cid_dca = dca_df.filter(pl.col("cohort_id") == cid)
            if not cid_dca.is_empty():
                cohort_dir = output_dir / "cohorts" / cid
                cohort_dir.mkdir(parents=True, exist_ok=True)
                match create_dca_net_benefit_chart(
                    cid_dca,
                    title=f"Decision Curve Analysis (Net Clinical Benefit) - {cid}",
                ):
                    case Success(chart):
                        cohort_svg = cohort_dir / "fig_dca_net_benefit.svg"
                        export_chart_svg(chart, cohort_svg)
                        if cid == "Gide-iAtlas":
                            export_chart_svg(chart, output_dir / "fig_dca_net_benefit_gide.svg")
                    case Failure(err):
                        print(f"Failed DCA chart for {cid}: {err}")

    print("\n=== Survival & DCA Benchmark Finished Successfully ===")
    return survival_df, dca_df


def main() -> None:
    args = parse_args()
    cohort_list = [c.strip() for c in args.cohorts.split(",") if c.strip()]
    run_survival_and_dca_benchmark(
        cohorts=cohort_list,
        output_dir=args.output_dir,
        n_bootstrap=args.n_bootstrap,
        enable_combined_pools=args.enable_combined_pools,
    )


if __name__ == "__main__":
    main()

