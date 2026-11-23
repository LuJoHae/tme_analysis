#!/usr/bin/env python3
"""benchmark_iatlas_icb_predictors.py: Systematic benchmark of ICB predictors on iAtlas cohorts.

Applies all transcriptomic, systems, and multi-omic prediction methods from tme_response
to the cBioPortal iAtlas cohorts in tme_datasets. Evaluates ROC-AUC, 95% bootstrap CIs,
PR-AUC, excess precision (ΔPR-AUC), and cohort predictability.
Constructs pooled cancer cohorts (Melanoma-Combined, RenalCell-Combined, PanCancer-Combined)
under both raw and z-score standardized pooling strategies.
Exports tabular summaries (.parquet, .csv). Publication SVGs are generated separately
by `3e2_plot_iatlas_icb_predictors.py`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence
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
    compute_roh_is_score,
    compute_single_gene_score,
    compute_tide_score,
    compute_wu_mias_score,
    evaluate_prediction,
    load_multi_omic_cohort,
    pool_cohort_stratum_data,
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
        default=Path("output/manuscript_pipeline/3e1_benchmark_iatlas_icb_predictors"),
        help="Output directory for tabular benchmark results (.parquet, .csv).",
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
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable combined cancer-type and pan-cancer cohort pools (default: True).",
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
        compute_jiang_texh_score,
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


    print("\n=== Benchmark Pipeline Finished Successfully ===")


if __name__ == "__main__":
    main()
