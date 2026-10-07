"""Orchestration engine for systematic perturbation sweeps across MultiOmicCohorts."""

from __future__ import annotations

from typing import Sequence
import polars as pl
from returns.result import Failure, Result, Success

from ..evaluation.metrics import evaluate_prediction
from ..schemas import MultiOmicCohort, PredictionResult
from ..signatures.cyt import compute_cyt_score
from ..signatures.gep import compute_gep_score
from ..signatures.impres import compute_impres_score
from ..signatures.ipres import compute_ipres_score
from ..signatures.single_gene import compute_single_gene_score
from ..models.tide import compute_tide_score
from ..deconvolution.wrapper import compute_cd8_infiltrate_predictor
from .expression import (
    apply_expression_jitter,
    apply_gene_dropout,
    apply_immune_dilution,
)
from .labels import apply_label_noise
from .models import (
    DilutionConfig,
    DropoutConfig,
    JitterConfig,
    LabelNoiseConfig,
    PerturbationSweepConfig,
)


def compute_standard_cohort_predictors(cohort: MultiOmicCohort) -> list[PredictionResult]:
    """Compute the 8 universal transcriptomic predictors for a cohort."""
    predictors: list[PredictionResult] = []

    match compute_cyt_score(cohort):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_impres_score(cohort):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_gep_score(cohort):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_single_gene_score(cohort, "CXCL9"):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_single_gene_score(cohort, "CD8A"):
        case Success(pred):
            predictors.append(pred)
        case Failure(_):
            pass

    match compute_ipres_score(cohort, invert_for_response=True):
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

    return predictors


def run_cohort_perturbation_sweep(
    cohort: MultiOmicCohort,
    sweep_config: PerturbationSweepConfig | None = None,
    n_bootstrap: int = 200,
) -> Result[pl.DataFrame, str]:
    """Execute a systematic parametric sweep across 4 perturbation modalities for a cohort."""
    cfg = sweep_config or PerturbationSweepConfig()
    clin_df = cohort.clinical_annotations

    if "response_binary" not in clin_df.columns:
        return Failure(f"Cohort {cohort.cohort_id} missing binary response annotations.")

    records: list[dict[str, object]] = []

    def evaluate_and_record(
        preds: Sequence[PredictionResult],
        clin: pl.DataFrame,
        p_type: str,
        intensity_val: float,
    ) -> None:
        for p in preds:
            match evaluate_prediction(
                prediction=p,
                clinical_df=clin,
                cohort_id=cohort.cohort_id,
                cancer_type=cohort.cancer_type,
                time_stratum="Pre",
                response_stratum="Standard",
                n_bootstrap=n_bootstrap,
            ):
                case Success(res):
                    records.append({
                        "cohort_id": res.cohort_id,
                        "cancer_type": res.cancer_type,
                        "predictor_name": res.predictor_name,
                        "category": res.category.value,
                        "perturbation_type": p_type,
                        "intensity": float(intensity_val),
                        "roc_auc": res.roc_auc,
                        "roc_auc_ci_lower": res.roc_auc_ci_lower,
                        "roc_auc_ci_upper": res.roc_auc_ci_upper,
                        "pr_auc": res.pr_auc,
                        "pr_auc_ci_lower": res.pr_auc_ci_lower,
                        "pr_auc_ci_upper": res.pr_auc_ci_upper,
                        "delta_pr_auc": res.delta_pr_auc,
                        "baseline_prevalence": res.baseline_prevalence,
                        "n_samples": res.n_samples,
                    })
                case Failure(_):
                    pass

    # 1. Expression Jitter Sweep
    for sigma in cfg.jitter_sigmas:
        p_cohort = apply_expression_jitter(cohort, JitterConfig(sigma=sigma))
        p_preds = compute_standard_cohort_predictors(p_cohort)
        evaluate_and_record(p_preds, clin_df, "jitter", sigma)

    # 2. Gene Dropout Sweep
    for rate in cfg.dropout_rates:
        if rate == 0.0:
            continue  # Already captured at baseline jitter=0
        p_cohort = apply_gene_dropout(cohort, DropoutConfig(dropout_rate=rate))
        p_preds = compute_standard_cohort_predictors(p_cohort)
        evaluate_and_record(p_preds, clin_df, "dropout", rate)

    # 3. Immune Dilution Sweep
    for alpha in cfg.dilution_factors:
        if alpha == 1.0:
            continue  # Already captured at baseline
        p_cohort = apply_immune_dilution(cohort, DilutionConfig(dilution_factor=alpha))
        p_preds = compute_standard_cohort_predictors(p_cohort)
        evaluate_and_record(p_preds, clin_df, "dilution", alpha)

    # 4. Label Noise Sweep
    base_preds = compute_standard_cohort_predictors(cohort)
    for eta in cfg.label_noise_rates:
        if eta == 0.0:
            continue  # Baseline
        p_clin = apply_label_noise(clin_df, LabelNoiseConfig(noise_rate=eta))
        evaluate_and_record(base_preds, p_clin, "label_noise", eta)

    if not records:
        return Failure(f"No valid perturbation evaluations generated for {cohort.cohort_id}.")

    return Success(pl.DataFrame(records))
