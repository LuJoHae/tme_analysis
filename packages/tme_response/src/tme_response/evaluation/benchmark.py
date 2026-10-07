"""Orchestrator for multi-cohort benchmark evaluations."""

from __future__ import annotations

from typing import Sequence
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from ..data.adapter import load_multi_omic_cohort
from ..deconvolution.wrapper import compute_cd8_infiltrate_predictor
from ..schemas import CohortBenchmarkResult, MultiOmicCohort, PredictionResult
from ..models.tide import compute_tide_score
from ..signatures.cyt import compute_cyt_score
from ..signatures.gep import compute_gep_score
from ..signatures.impres import compute_impres_score
from ..signatures.ipres import compute_ipres_score
from ..signatures.single_gene import compute_single_gene_score
from ..synergy.dna_rna import compute_dna_rna_composite
from ..synergy.gating import apply_antigen_presentation_gating
from .metrics import evaluate_prediction


def run_cohort_benchmark(cohort_id: str) -> Result[tuple[CohortBenchmarkResult, ...], str]:
    """Execute all signature, systems, and multi-omic predictors on a single cohort."""
    match load_multi_omic_cohort(cohort_id):
        case Failure(err):
            return Failure(err)
        case Success(cohort):
            pass

    clin_df = cohort.clinical_annotations
    predictions: list[PredictionResult] = []

    # 1. Transcriptomic Signatures
    match compute_cyt_score(cohort):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    match compute_impres_score(cohort):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    match compute_gep_score(cohort):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    match compute_single_gene_score(cohort, "CXCL9"):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    match compute_single_gene_score(cohort, "CD8A"):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    match compute_ipres_score(cohort, invert_for_response=True):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    # 2. TIDE Evasion / Response
    match compute_tide_score(cohort, invert_for_response=True):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    # 3. Cellular Infiltration
    match compute_cd8_infiltrate_predictor(cohort):
        case Success(pred):
            predictions.append(pred)
        case Failure(_):
            pass

    # 4. Multi-Omic DNA-RNA Synergy (if TMB is available)
    if isinstance(cohort.tmb_scores, Some):
        match compute_single_gene_score(cohort, "CXCL9"):
            case Success(cxcl9_pred):
                match compute_dna_rna_composite(cohort, cxcl9_pred):
                    case Success(comp_pred):
                        predictions.append(comp_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    # 5. Antigen Presentation Gating (if mutations available)
    if isinstance(cohort.driver_mutations, Some):
        match compute_single_gene_score(cohort, "CXCL9"):
            case Success(cxcl9_pred):
                match apply_antigen_presentation_gating(cohort, cxcl9_pred):
                    case Success(gated_pred):
                        predictions.append(gated_pred)
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    # Evaluate each prediction
    results: list[CohortBenchmarkResult] = []
    for pred in predictions:
        match evaluate_prediction(pred, clin_df, cohort.cohort_id, cohort.cancer_type):
            case Success(res):
                results.append(res)
            case Failure(_):
                pass

    return Success(tuple(results))


def run_multi_cohort_benchmark(
    cohort_ids: Sequence[str],
) -> Result[pl.DataFrame, str]:
    """Run systematic benchmark across a sequence of registered tme_datasets cohorts."""
    all_results: list[CohortBenchmarkResult] = []

    for cid in cohort_ids:
        match run_cohort_benchmark(cid):
            case Success(res_tuple):
                all_results.extend(res_tuple)
            case Failure(_):
                # Skip cohorts that fail loading / filtering gracefully
                continue

    if not all_results:
        return Failure("No benchmark results could be calculated across the provided cohorts.")

    df = pl.DataFrame([
        {
            "cohort_id": r.cohort_id,
            "cancer_type": r.cancer_type,
            "predictor_name": r.predictor_name,
            "category": r.category.value,
            "time_stratum": r.time_stratum,
            "response_stratum": r.response_stratum,
            "n_samples": r.n_samples,
            "n_responders": r.n_responders,
            "n_non_responders": r.n_non_responders,
            "roc_auc": r.roc_auc,
            "roc_auc_ci_lower": r.roc_auc_ci_lower,
            "roc_auc_ci_upper": r.roc_auc_ci_upper,
            "p_value_vs_half": r.p_value_vs_half,
            "pr_auc": r.pr_auc,
            "brier_score": r.brier_score,
        }
        for r in all_results
    ])

    return Success(df)
