#!/usr/bin/env python3
"""benchmark_iatlas_icb_predictors.py: Systematic benchmark of ICB predictors on iAtlas cohorts.

Applies all transcriptomic, systems, and multi-omic prediction methods from tme_response
to the cBioPortal iAtlas cohorts in tme_datasets. Evaluates ROC-AUC, 95% bootstrap CIs,
PR-AUC, excess precision (ΔPR-AUC), and cohort predictability.
Constructs pooled cancer cohorts (Melanoma-Combined, RenalCell-Combined, PanCancer-Combined)
under both raw and z-score standardized pooling strategies.
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
from tme_datasets.query import load_dataset
from tme_response import (
    CohortBenchmarkResult,
    CohortPredictabilityResult,
    MultiOmicCohort,
    PredictionResult,
    PredictorCategory,
    UNIVERSAL_RNA_PREDICTORS,
    apply_antigen_presentation_gating,
    calculate_discrimination_bootstrap,
    calculate_roc_auc,
    compute_ayers_ifng6_score,
    compute_cd8_duo_score,
    compute_cd8_infiltrate_predictor,
    compute_cohort_predictability,
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
    compute_meta_analytic_auc,
    compute_nurmik_cafs_score,
    compute_perturbation_resilience,
    compute_roh_is_score,
    compute_single_gene_score,
    compute_tide_score,
    compute_wu_mias_score,
    create_auc_heatmap,
    create_benchmark_bar_chart,
    create_cohort_predictability_forest_plot,
    create_cross_cohort_resilience_heatmap,
    create_meta_analytic_decay_chart,
    create_pan_cohort_faceted_decay_chart,
    create_perturbation_decay_chart,
    create_pooling_comparison_chart,
    create_resilience_ranking_chart,
    create_roc_chart,
    create_summary_forest_plot,
    evaluate_prediction,
    export_chart_svg,
    load_multi_omic_cohort,
    pool_cohort_stratum_data,
    run_cohort_perturbation_sweep,
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
        "cancer_type": "Renal Cell Carcinoma",
        "cohort_ids": ["McDermott-iAtlas", "Choueiri-iAtlas"],
    },
    "PanCancer-Combined": {
        "cancer_type": "Pan-Cancer",
        "cohort_ids": list(DEFAULT_IATLAS_COHORTS),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark ICB response predictors on cBioPortal iAtlas cohorts."
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default=",".join(DEFAULT_IATLAS_COHORTS),
        help="Comma-separated list of iAtlas cohort IDs to benchmark.",
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
        default=1000,
        help="Number of bootstrap resamples for 95%% confidence intervals (default: 1000).",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=10,
        help="Minimum required annotated samples per stratum (default: 10).",
    )
    parser.add_argument(
        "--min-class-samples",
        type=int,
        default=3,
        help="Minimum required samples per binary class (default: 3).",
    )
    parser.add_argument(
        "--enable-combined-pools",
        action="store_true",
        default=True,
        help="Enable combined cancer-type and pan-cancer cohort pools (default: True).",
    )
    parser.add_argument(
        "--run-perturbations",
        action="store_true",
        default=True,
        help="Run systematic data perturbation sweeps (default: True).",
    )
    parser.add_argument(
        "--perturbation-cohorts",
        type=str,
        default=",".join(DEFAULT_IATLAS_COHORTS),
        help="Comma-separated cohorts to stress-test under perturbations (default: all cohorts).",
    )
    parser.add_argument(
        "--perturbation-bootstrap",
        type=int,
        default=200,
        help="Number of bootstrap iterations for perturbation points (default: 200).",
    )
    parser.add_argument(
        "--only-perturbations",
        action="store_true",
        default=False,
        help="Skip baseline recomputation and only run perturbation sweeps.",
    )
    parser.add_argument(
        "--force-perturbations",
        action="store_true",
        default=False,
        help="Force recomputation of perturbation sweeps even if cached evaluations exist.",
    )
    parser.add_argument(
        "--run-survival",
        action="store_true",
        default=False,
        help="Run Time-to-Event Survival (OS/PFS) and Decision Curve Analysis (DCA).",
    )
    parser.add_argument(
        "--survival-bootstrap",
        type=int,
        default=500,
        help="Number of bootstrap resamples for survival C-index (default: 500).",
    )
    return parser.parse_args()




def generate_combinatorial_strata(
    clinical_df: pl.DataFrame,
    min_samples: int = 10,
    min_class_samples: int = 3,
) -> dict[str, pl.DataFrame]:
    """Generate all valid combinatorial strata of timepoints x response definitions.

    Timepoints:
      - 'Pre': Baseline pre-treatment / unknown
      - 'On': On-treatment biopsies
      - 'All': All available biopsies
    Response:
      - 'Standard': Binary response (CR/PR/DCB vs SD/PD/NDB)
      - 'Extreme': Extreme responders (CR/PR vs PD, excluding SD)
    """
    strata: dict[str, pl.DataFrame] = {}

    has_tp = "biopsy_timepoint" in clinical_df.columns
    has_recist = "response_recist" in clinical_df.columns

    time_filters = {
        "Pre": (
            (pl.col("biopsy_timepoint").is_in(["Pre", "Unknown"]) | pl.col("biopsy_timepoint").is_null())
            if has_tp
            else pl.lit(True)
        ),
        "On": (
            pl.col("biopsy_timepoint") == "On-Treatment"
            if has_tp
            else pl.lit(False)
        ),
        "All": pl.lit(True),
    }

    resp_filters = {
        "Standard": pl.col("response_binary").is_finite(),
        "Extreme": (
            pl.col("response_recist").is_in(["CR", "PR", "PD"])
            if has_recist
            else pl.col("response_binary").is_finite()
        ),
    }

    for t_name, t_cond in time_filters.items():
        for r_name, r_cond in resp_filters.items():
            sub = clinical_df.filter(t_cond & r_cond)
            if sub.is_empty():
                continue

            # Standardize binary response for extreme responders if needed
            if r_name == "Extreme" and has_recist:
                sub = sub.with_columns(
                    pl.when(pl.col("response_recist").is_in(["CR", "PR"]))
                    .then(1.0)
                    .otherwise(0.0)
                    .alias("response_binary")
                )

            # Check sample and class count constraints
            n_total = len(sub)
            n_resp = int(sub["response_binary"].sum() or 0)
            n_non = n_total - n_resp

            if (
                n_total >= min_samples
                and n_resp >= min_class_samples
                and n_non >= min_class_samples
            ):
                stratum_key = f"{t_name}_{r_name}"
                strata[stratum_key] = sub

    return strata


def compute_all_cohort_predictors(cohort: MultiOmicCohort) -> list[PredictionResult]:
    """Compute all available transcriptomic, systems, and multi-omic predictors for a cohort."""
    predictors: list[PredictionResult] = []

    # 1. Single-Gene Transcriptomic Predictors
    for g in ["CXCL9", "CD8A", "PDCD1", "CD274", "CTLA4"]:
        match compute_single_gene_score(cohort, g):
            case Success(pred):
                predictors.append(pred)
            case Failure(_):
                pass

    # 2. Curated Transcriptomic & COMPASS Table S2 Baselines
    signature_callers = [
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
    ]

    for caller in signature_callers:
        match caller(cohort):
            case Success(pred):
                predictors.append(pred)
            case Failure(_):
                pass

    # 2. Systems Evasion Model (TIDE)
    match compute_tide_score(cohort, invert_for_response=True):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    # 3. Cellular Infiltrate (MCP-counter CD8 T cells)
    match compute_cd8_infiltrate_predictor(cohort):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    # 4. Standalone Genomic Baseline: TMB
    if isinstance(cohort.tmb_scores, Some):
        tmb_df = cohort.tmb_scores.unwrap()
        tmb_pred = PredictionResult(
            predictor_name="TMB_Genomic",
            category=PredictorCategory.GENOMIC,
            predictions=tmb_df.select([
                pl.col("sample_id"),
                pl.col("tmb_per_mb").alias("score"),
            ]),
        )
        predictors.append(tmb_pred)

        # Multi-omic composite: CXCL9 + TMB
        match compute_single_gene_score(cohort, "CXCL9"):
            case Success(cxcl9_pred):
                match compute_dna_rna_composite(cohort, cxcl9_pred):
                    case Success(comp_pred):
                        predictors.append(comp_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

        # Multi-omic composite: GEP + TMB
        match compute_gep_score(cohort):
            case Success(gep_pred):
                match compute_dna_rna_composite(cohort, gep_pred):
                    case Success(comp_pred):
                        predictors.append(comp_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    # 5. Antigen Presentation Gating (B2M / JAK1)
    if isinstance(cohort.driver_mutations, Some):
        match compute_single_gene_score(cohort, "CXCL9"):
            case Success(cxcl9_pred):
                match apply_antigen_presentation_gating(cohort, cxcl9_pred):
                    case Success(gated_pred):
                        predictors.append(gated_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    return predictors


def main() -> None:
    args = parse_args()
    repo_root = find_repo_root()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_list = [c.strip() for c in args.cohorts.split(",") if c.strip()]
    print(f"=== Starting Systematic ICB Response Benchmark ===")
    print(f"Target Individual Cohorts ({len(cohort_list)}): {', '.join(cohort_list)}")
    print(f"Output Directory: {output_dir}")
    print(f"Bootstrap Resamples: {args.n_bootstrap}")

    all_benchmark_records: list[dict[str, object]] = []

    if args.only_perturbations:
        parquet_path = output_dir / "iatlas_benchmark_metrics.parquet"
        if not parquet_path.exists():
            print(f"[ERROR] Cannot run only perturbations: {parquet_path} does not exist.")
            return
        results_df = pl.read_parquet(parquet_path)
        print(f"Loaded existing baseline benchmark dataset ({len(results_df)} evaluations).")
        pred_parquet = output_dir / "cohort_predictability_metrics.parquet"
        predictability_df = pl.read_parquet(pred_parquet) if pred_parquet.exists() else pl.DataFrame()
    else:
        # Cache for cohort pooling: cohort_id -> dict(cancer_type, strata, predictors)
        cohort_cache: dict[str, dict[str, object]] = {}

        # ---------------------------------------------------------
        # 1. Process Individual Cohorts
        # ---------------------------------------------------------
        for cohort_id in cohort_list:
            print(f"\n--- Loading Cohort: {cohort_id} ---")
            match load_multi_omic_cohort(cohort_id, repo_root=repo_root, filter_baseline=False):
                case Failure(err):
                    print(f"  [SKIPPED] Failed to load {cohort_id}: {err}")
                    continue
                case Success(cohort):
                    n_samples = len(cohort.sample_ids)
                    has_tmb = isinstance(cohort.tmb_scores, Some)
                    has_mut = isinstance(cohort.driver_mutations, Some)
                    print(
                        f"  Loaded {cohort_id} ({cohort.cancer_type}): "
                        f"{n_samples} total samples (TMB: {has_tmb}, Mutations: {has_mut})"
                    )

                    predictors = compute_all_cohort_predictors(cohort)
                    print(f"  Calculated {len(predictors)} predictors.")

                    strata = generate_combinatorial_strata(
                        cohort.clinical_annotations,
                        min_samples=args.min_samples,
                        min_class_samples=args.min_class_samples,
                    )
                    print(f"  Valid Strata ({len(strata)}): {', '.join(strata.keys())}")

                    # Save into cache for combined pool processing
                    cohort_cache[cohort_id] = {
                        "cancer_type": cohort.cancer_type,
                        "strata": strata,
                        "predictors": predictors,
                    }

                    # Evaluate individual cohort
                    for stratum_key, stratum_clinical in strata.items():
                        t_stratum, r_stratum = stratum_key.split("_", 1)

                        for pred in predictors:
                            match evaluate_prediction(
                                prediction=pred,
                                clinical_df=stratum_clinical,
                                cohort_id=cohort.cohort_id,
                                cancer_type=cohort.cancer_type,
                                time_stratum=t_stratum,
                                response_stratum=r_stratum,
                                pooling_strategy="cohort",
                                n_bootstrap=args.n_bootstrap,
                            ):
                                case Success(res):
                                    all_benchmark_records.append({
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
                                    continue

        # ---------------------------------------------------------
        # 2. Process Combined Cohort Pools (Raw vs. Standardized)
        # ---------------------------------------------------------
        if args.enable_combined_pools and cohort_cache:
            print(f"\n=== Processing Combined Cohort Pools ===")
            for pool_id, pool_spec in COMBINED_COHORT_SPECS.items():
                pool_cancer_type = pool_spec["cancer_type"]
                member_cohorts = [cid for cid in pool_spec["cohort_ids"] if cid in cohort_cache]

                if len(member_cohorts) < 2:
                    print(f"  [SKIPPED] {pool_id}: Fewer than 2 available member cohorts.")
                    continue

                print(f"\n--- Constructing Pool: {pool_id} ({pool_cancer_type}) from {len(member_cohorts)} cohorts ---")

                # Discover all distinct stratum keys across member cohorts
                all_stratum_keys: set[str] = set()
                for cid in member_cohorts:
                    all_stratum_keys.update(cohort_cache[cid]["strata"].keys())

                for stratum_key in sorted(all_stratum_keys):
                    t_stratum, r_stratum = stratum_key.split("_", 1)

                    # Assemble cohort strata
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

                    # Run both pooling strategies: standardized (z-score) and raw
                    for strategy in ["standardized", "raw"]:
                        do_std = strategy == "standardized"
                        match pool_cohort_stratum_data(
                            strata_tuples,
                            group_id=pool_id,
                            cancer_type=pool_cancer_type,
                            standardize=do_std,
                        ):
                            case Failure(err):
                                print(f"    Failed pooling {pool_id} ({stratum_key}, {strategy}): {err}")
                                continue
                            case Success((pooled_clin, pooled_preds)):
                                n_total = len(pooled_clin)
                                n_resp = int(pooled_clin["response_binary"].sum() or 0)
                                n_non = n_total - n_resp

                                if (
                                    n_total < args.min_samples
                                    or n_resp < args.min_class_samples
                                    or n_non < args.min_class_samples
                                ):
                                    continue

                                print(
                                    f"    Evaluated {pool_id} ({stratum_key}, strategy={strategy}): "
                                    f"{n_total} samples ({n_resp} responders, {len(pooled_preds)} predictors)"
                                )

                                for p in pooled_preds:
                                    match evaluate_prediction(
                                        prediction=p,
                                        clinical_df=pooled_clin,
                                        cohort_id=pool_id,
                                        cancer_type=pool_cancer_type,
                                        time_stratum=t_stratum,
                                        response_stratum=r_stratum,
                                        pooling_strategy=strategy,
                                        n_bootstrap=args.n_bootstrap,
                                    ):
                                        case Success(res):
                                            all_benchmark_records.append({
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
                                            continue

        if not all_benchmark_records:
            print("\n[ERROR] No valid benchmark results generated across cohorts.")
            return

        # Consolidate results into a Polars DataFrame
        results_df = pl.DataFrame(all_benchmark_records)
        print(f"\n=== Benchmark Complete: {len(results_df)} evaluations generated ===")

        # ---------------------------------------------------------
        # 3. Export Tabular Datasets
        # ---------------------------------------------------------
        parquet_path = output_dir / "iatlas_benchmark_metrics.parquet"
        csv_path = output_dir / "iatlas_benchmark_metrics.csv"
        results_df.write_parquet(parquet_path)
        results_df.write_csv(csv_path)
        print(f"Saved: {parquet_path}")
        print(f"Saved: {csv_path}")

        # Compute and Export Cohort Predictability Index
        match compute_cohort_predictability(results_df):
            case Success(predictability_df):
                pred_parquet = output_dir / "cohort_predictability_metrics.parquet"
                pred_csv = output_dir / "cohort_predictability_metrics.csv"
                predictability_df.write_parquet(pred_parquet)
                predictability_df.write_csv(pred_csv)
                print(f"Saved: {pred_parquet}")
                print(f"Saved: {pred_csv}")
            case Failure(err):
                print(f"Failed to compute cohort predictability: {err}")
                predictability_df = pl.DataFrame()

    # ---------------------------------------------------------
    # 4. Generate and Export Publication Vector Figures (SVG)
    # ---------------------------------------------------------
    print("\n--- Generating Publication Vector Figures (SVG) ---")

    # Primary display cohort subset: Individual cohorts + Standardized pools (exclude unnormalized raw pools from primary heatmaps)
    primary_display_df = results_df.filter(
        pl.col("pooling_strategy").is_in(["cohort", "standardized"])
    )

    pre_std_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Standard")
    )
    extreme_primary = primary_display_df.filter(
        (pl.col("time_stratum") == "Pre") & (pl.col("response_stratum") == "Extreme")
    )

    # 1. Heatmap: ROC-AUC Pre-Treatment Standard Response
    if not pre_std_primary.is_empty():
        match create_auc_heatmap(pre_std_primary, metric="roc_auc", title="Pre-Treatment Standard ICB Response (ROC-AUC)"):
            case Success(chart):
                svg_path = output_dir / "fig_roc_auc_heatmap_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export ROC-AUC heatmap SVG: {err}")
            case Failure(err):
                print(f"ROC-AUC heatmap generation failed: {err}")

    # 2. Heatmap: ROC-AUC Extreme Responders (CR/PR vs PD)
    if not extreme_primary.is_empty():
        match create_auc_heatmap(extreme_primary, metric="roc_auc", title="Pre-Treatment Extreme Responders: CR/PR vs PD (ROC-AUC)"):
            case Success(chart):
                svg_path = output_dir / "fig_roc_auc_heatmap_pre_extreme.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export extreme ROC-AUC heatmap SVG: {err}")
            case Failure(err):
                print(f"Extreme ROC-AUC heatmap generation failed: {err}")

    # 3. Heatmap: PR-AUC Pre-Treatment Standard Response
    if not pre_std_primary.is_empty():
        match create_auc_heatmap(pre_std_primary, metric="pr_auc", title="Pre-Treatment Standard ICB Response (PR-AUC)"):
            case Success(chart):
                svg_path = output_dir / "fig_pr_auc_heatmap_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export PR-AUC heatmap SVG: {err}")
            case Failure(err):
                print(f"PR-AUC heatmap generation failed: {err}")

    # 4. Heatmap: PR-AUC Extreme Responders
    if not extreme_primary.is_empty():
        match create_auc_heatmap(extreme_primary, metric="pr_auc", title="Pre-Treatment Extreme Responders: CR/PR vs PD (PR-AUC)"):
            case Success(chart):
                svg_path = output_dir / "fig_pr_auc_heatmap_pre_extreme.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export extreme PR-AUC heatmap SVG: {err}")
            case Failure(err):
                print(f"Extreme PR-AUC heatmap generation failed: {err}")

    # 5. Forest Plot: Predictor Summary ROC-AUC ± 95% CI
    if not pre_std_primary.is_empty():
        match create_summary_forest_plot(pre_std_primary, metric="roc_auc", title="Cross-Cohort Predictor Performance (Mean ROC-AUC ± 95% CI)"):
            case Success(chart):
                svg_path = output_dir / "fig_forest_summary_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export predictor ROC-AUC forest SVG: {err}")
            case Failure(err):
                print(f"Predictor ROC-AUC forest plot failed: {err}")

    # 6. Forest Plot: Predictor Summary PR-AUC ± 95% CI
    if not pre_std_primary.is_empty():
        match create_summary_forest_plot(pre_std_primary, metric="pr_auc", title="Cross-Cohort Predictor Performance (Mean PR-AUC ± 95% CI)"):
            case Success(chart):
                svg_path = output_dir / "fig_forest_summary_prauc_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export predictor PR-AUC forest SVG: {err}")
            case Failure(err):
                print(f"Predictor PR-AUC forest plot failed: {err}")

    # 7. Forest Plot: Cohort Predictability Index (Pre_Standard)
    if not predictability_df.is_empty():
        pred_pre_std = predictability_df.filter(
            (pl.col("time_stratum") == "Pre")
            & (pl.col("response_stratum") == "Standard")
            & (pl.col("pooling_strategy").is_in(["cohort", "standardized"]))
        )
        if not pred_pre_std.is_empty():
            match create_cohort_predictability_forest_plot(
                pred_pre_std, title="Cohort Predictability Index (Mean RNA Biomarker ROC-AUC ± SD)"
            ):
                case Success(chart):
                    svg_path = output_dir / "fig_cohort_predictability_forest.svg"
                    match export_chart_svg(chart, svg_path):
                        case Success(p):
                            print(f"Saved: {p}")
                        case Failure(err):
                            print(f"Failed to export cohort predictability forest SVG: {err}")
                case Failure(err):
                    print(f"Cohort predictability forest plot failed: {err}")

    # 8. Comparison Plot: Raw vs Standardized Pooling
    pre_std_all_pools = results_df.filter(
        (pl.col("time_stratum") == "Pre")
        & (pl.col("response_stratum") == "Standard")
        & (pl.col("pooling_strategy").is_in(["raw", "standardized"]))
    )
    if not pre_std_all_pools.is_empty():
        match create_pooling_comparison_chart(
            pre_std_all_pools, title="Impact of Score Standardization on Pooled ROC-AUC"
        ):
            case Success(chart):
                svg_path = output_dir / "fig_pooling_comparison_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export pooling comparison SVG: {err}")
            case Failure(err):
                print(f"Pooling comparison chart failed: {err}")

    # 9. Grouped Bar Chart
    if not pre_std_primary.is_empty():
        match create_benchmark_bar_chart(pre_std_primary, metric="roc_auc"):
            case Success(chart):
                svg_path = output_dir / "fig_benchmark_bars_pre_standard.svg"
                match export_chart_svg(chart, svg_path):
                    case Success(p):
                        print(f"Saved: {p}")
                    case Failure(err):
                        print(f"Failed to export bar chart SVG: {err}")
            case Failure(err):
                print(f"Bar chart generation failed: {err}")

    # ---------------------------------------------------------
    # 5. Perturbation & Robustness Benchmarking
    # ---------------------------------------------------------
    if args.run_perturbations:
        print("\n=== Running Systematic Data Perturbation & Robustness Benchmarks ===")
        pert_cohort_ids = [c.strip() for c in args.perturbation_cohorts.split(",") if c.strip()]

        pert_pq = output_dir / "perturbation_benchmark_metrics.parquet"
        pert_csv = output_dir / "perturbation_benchmark_metrics.csv"

        cached_records: list[dict[str, object]] = []
        already_computed: set[str] = set()
        if pert_pq.exists() and not args.force_perturbations:
            cached_df = pl.read_parquet(pert_pq)
            cached_records = cached_df.to_dicts()
            already_computed = set(cached_df["cohort_id"].unique().to_list())
            print(f"Loaded {len(cached_df)} existing perturbation evaluations across {len(already_computed)} cohorts: {', '.join(sorted(already_computed))}")

        new_pert_records: list[dict[str, object]] = []
        for p_cid in pert_cohort_ids:
            if p_cid in already_computed:
                print(f"\n--- Cohort {p_cid}: Already evaluated in cache; skipping recomputation ---")
                continue

            print(f"\n--- Perturbation Sweep on: {p_cid} ---")
            match load_multi_omic_cohort(p_cid, repo_root=repo_root, filter_baseline=True):
                case Failure(err):
                    print(f"  [SKIPPED] Failed loading {p_cid}: {err}")
                    continue
                case Success(p_cohort):
                    match run_cohort_perturbation_sweep(
                        p_cohort, n_bootstrap=args.perturbation_bootstrap
                    ):
                        case Success(pert_df):
                            print(f"  Generated {len(pert_df)} perturbation evaluations on {p_cid}.")
                            new_pert_records.extend(pert_df.to_dicts())
                        case Failure(err):
                            print(f"  Perturbation sweep failed on {p_cid}: {err}")

        all_pert_records = cached_records + new_pert_records

        if all_pert_records:
            pert_benchmark_df = pl.DataFrame(all_pert_records)
            pert_benchmark_df.write_parquet(pert_pq)
            pert_benchmark_df.write_csv(pert_csv)
            print(f"Saved: {pert_pq}")
            print(f"Saved: {pert_csv}")

            # Compute resilience index across all cohorts
            match compute_perturbation_resilience(pert_benchmark_df):
                case Success(resilience_df):
                    resil_pq = output_dir / "perturbation_resilience_metrics.parquet"
                    resil_csv = output_dir / "perturbation_resilience_metrics.csv"
                    resilience_df.write_parquet(resil_pq)
                    resilience_df.write_csv(resil_csv)
                    print(f"Saved: {resil_pq}")
                    print(f"Saved: {resil_csv}")

                    # 1. Cross-Cohort PRI Heatmap
                    match create_cross_cohort_resilience_heatmap(resilience_df):
                        case Success(chart):
                            svg_p = output_dir / "fig_perturbation_resilience_cross_cohort_heatmap.svg"
                            export_chart_svg(chart, svg_p)
                            print(f"Saved: {svg_p}")
                        case Failure(err):
                            print(f"Failed to export cross-cohort resilience heatmap: {err}")

                    # 2. Cohort-Specific Resilience Ranking SVGs
                    for cid in pert_benchmark_df["cohort_id"].unique().to_list():
                        cid_resil = resilience_df.filter(pl.col("cohort_id") == cid)
                        if not cid_resil.is_empty():
                            cohort_dir = output_dir / "cohorts" / cid
                            cohort_dir.mkdir(parents=True, exist_ok=True)
                            match create_resilience_ranking_chart(
                                cid_resil,
                                title=f"Predictor Perturbation Resilience Index (PRI) - {cid}",
                                subtitle=f"Cohort: {cid} | Metric: Normalized Area Retention (1.0 = Complete Resilience)",
                            ):
                                case Success(chart):
                                    export_chart_svg(chart, cohort_dir / "fig_perturbation_resilience_ranking.svg")
                                    if cid == "Gide-iAtlas":
                                        export_chart_svg(chart, output_dir / "fig_perturbation_resilience_ranking.svg")
                                case Failure(err):
                                    print(f"Failed to export resilience ranking chart for {cid}: {err}")
                case Failure(err):
                    print(f"Failed to compute resilience index: {err}")

            # 3. Pan-Cohort 3x3 Faceted Decay Grids (Export shaded, no_uncertainty, and errorbars)
            for p_type in ["jitter", "dropout", "dilution", "label_noise"]:
                for u_mode, suffix in [("shaded", "_shaded"), ("none", "_no_uncertainty"), ("errorbars", "_errorbars"), ("shaded", "")]:
                    match create_pan_cohort_faceted_decay_chart(
                        pert_benchmark_df, perturbation_type=p_type, columns=3, uncertainty=u_mode
                    ):
                        case Success(chart):
                            svg_p = output_dir / f"fig_perturbation_decay_{p_type}_all_cohorts{suffix}.svg"
                            export_chart_svg(chart, svg_p)
                            print(f"Saved: {svg_p}")
                        case Failure(err):
                            print(f"Failed to export pan-cohort faceted decay chart for {p_type} ({u_mode}): {err}")

            # 4. Pan-Cohort Meta-Analytic Mean Decay Curves (Export shaded, no_uncertainty, and errorbars)
            for p_type in ["jitter", "dropout", "dilution", "label_noise"]:
                for u_mode, suffix in [("shaded", "_shaded"), ("none", "_no_uncertainty"), ("errorbars", "_errorbars"), ("shaded", "")]:
                    match create_meta_analytic_decay_chart(
                        pert_benchmark_df, perturbation_type=p_type, uncertainty=u_mode
                    ):
                        case Success(chart):
                            svg_p = output_dir / f"fig_perturbation_decay_pan_cohort_mean_{p_type}{suffix}.svg"
                            export_chart_svg(chart, svg_p)
                            print(f"Saved: {svg_p}")
                        case Failure(err):
                            print(f"Failed to export meta-analytic decay chart for {p_type} ({u_mode}): {err}")

            # 5. Individual Cohort Decay SVGs for ALL Cohorts
            for cid in pert_benchmark_df["cohort_id"].unique().to_list():
                cid_pert = pert_benchmark_df.filter(pl.col("cohort_id") == cid)
                cohort_dir = output_dir / "cohorts" / cid
                cohort_dir.mkdir(parents=True, exist_ok=True)
                for p_type in ["jitter", "dropout", "dilution", "label_noise"]:
                    for u_mode, suffix in [("shaded", "_shaded"), ("none", "_no_uncertainty"), ("errorbars", "_errorbars"), ("shaded", "")]:
                        match create_perturbation_decay_chart(
                            cid_pert,
                            perturbation_type=p_type,
                            uncertainty=u_mode,
                        ):
                            case Success(chart):
                                export_chart_svg(chart, cohort_dir / f"fig_perturbation_decay_{p_type}{suffix}.svg")
                                if cid == "Gide-iAtlas":
                                    export_chart_svg(chart, output_dir / f"fig_perturbation_decay_{p_type}{suffix}.svg")
                            case Failure(err):
                                print(f"Failed to export decay chart for {p_type} on {cid} ({u_mode}): {err}")

    # ---------------------------------------------------------
    # 5. Time-to-Event Survival (OS/PFS) and Decision Curve Analysis (DCA)
    # ---------------------------------------------------------
    if args.run_survival:
        from benchmark_survival_and_dca import run_survival_and_dca_benchmark

        print("\n=== Running Survival & Decision Curve Analysis Benchmark ===")
        run_survival_and_dca_benchmark(
            cohorts=cohort_list,
            output_dir=output_dir,
            n_bootstrap=args.survival_bootstrap,
            enable_combined_pools=args.enable_combined_pools,
        )

    print("\n=== Benchmark Pipeline Finished Successfully ===")


if __name__ == "__main__":
    main()

