"""Pure functional single-cohort hyperparameter optimization for deconvolution reference construction and classifiers."""

from __future__ import annotations

import gc
import time
from typing import Any, Mapping, Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..deconvolution.models import DeconvolutionReferenceConfig
from ..deconvolution.builder import build_deconvolution_reference
from ..logging import get_logger
from .classifier_hpo import (
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    InnerOptimizationResult,
    TransformType,
    evaluate_single_cohort_stratified_cv,
    run_single_cohort_classifier_hpo,
)
from .cohort_sampler import sample_single_cell_cohorts
from .hpo_models import (
    FidelityRung,
    HPORunResult,
    HPOSearchSpace,
    HPOTrialConfig,
    TrialEvaluationResult,
)
from .hpo_optimizer import (
    batch_richardson_lucy_deconvolution,
    generate_trial_config,
)

logger = get_logger("sampling.single_cohort_hpo")


class TrialReferenceContext(BaseModel):
    """Immutable cached reference and bulk matrices for a single HPO trial."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    norm_ref: np.ndarray
    sub_bulk: np.ndarray
    state_labels: tuple[str, ...]
    ref_genes: tuple[str, ...]
    collinearity_max: float
    condition_number: float
    total_cells_sampled: int


def prepare_trial_reference_context(
    config: HPOTrialConfig,
    bulk_mat: np.ndarray,
    bulk_genes: Sequence[str],
    target_cohort_name: str = "target",
    bulk_ensembl_genes: Sequence[str] | None = None,
    bulk_symbol_genes: Sequence[str] | None = None,
) -> Result[TrialReferenceContext, str]:
    """Sample single-cell cohorts, cluster, build reference, and align gene spaces with bulk cohort."""
    # Step 1: Sample and harmonize single cells
    sample_res = sample_single_cell_cohorts(
        sampling_spec=config.sampling_spec,
        cluster_spec=config.cluster_spec,
        harmonize_config=config.harmonize_config,
    )
    match sample_res:
        case Failure(err):
            return Failure(f"Sampling pipeline failed: {err}")
        case Success(sampled_data):
            pass

    # Step 2: Build deconvolution reference profile
    ref_res = build_deconvolution_reference(
        data_source=sampled_data.adata,
        config=config.reference_config,
    )
    match ref_res:
        case Failure(err):
            return Failure(f"Reference builder failed: {err}")
        case Success(ref_result):
            pass

    ref_mat, state_labels = ref_result.to_instaprism()
    ref_genes = ref_result.gene_names

    # Check collinearity and condition number
    if ref_result.collinearity.height > 0 and "correlation" in ref_result.collinearity.columns:
        collinearity_max = float(ref_result.collinearity["correlation"].max() or 0.0)
    else:
        collinearity_max = 0.0

    try:
        s_vals = np.linalg.svd(ref_mat, compute_uv=False)
        condition_number = float(s_vals[0] / (s_vals[-1] + 1e-12))
    except Exception:
        condition_number = 1e6

    # Step 3: Align gene spaces
    ref_symbols = getattr(ref_result, "gene_symbols", ())
    ref_ens_to_idx = {g.upper(): idx for idx, g in enumerate(ref_genes)}
    ref_sym_to_idx = {s.upper(): idx for idx, s in enumerate(ref_symbols) if s}

    bulk_ens = bulk_ensembl_genes if bulk_ensembl_genes is not None else bulk_genes
    bulk_sym = bulk_symbol_genes if bulk_symbol_genes is not None else bulk_genes

    # Evaluate overlap across all 4 potential space pairings
    c_ens_ens = sum(1 for g in bulk_ens if str(g).upper() in ref_ens_to_idx)
    c_sym_sym = sum(1 for g in bulk_sym if str(g).upper() in ref_sym_to_idx) if ref_sym_to_idx else 0
    c_sym_ens = sum(1 for g in bulk_sym if str(g).upper() in ref_ens_to_idx)
    c_ens_sym = sum(1 for g in bulk_ens if str(g).upper() in ref_sym_to_idx) if ref_sym_to_idx else 0

    pairs = [
        (c_ens_ens, bulk_ens, ref_ens_to_idx),
        (c_sym_sym, bulk_sym, ref_sym_to_idx),
        (c_sym_ens, bulk_sym, ref_ens_to_idx),
        (c_ens_sym, bulk_ens, ref_sym_to_idx),
    ]
    best_count, chosen_bulk_genes, chosen_ref_map = max(pairs, key=lambda x: x[0])

    shared_genes: list[str] = []
    ref_sub_idx: list[int] = []
    bulk_sub_idx: list[int] = []
    seen_ref_idx: set[int] = set()

    bulk_gene_to_idx = {g: idx for idx, g in enumerate(chosen_bulk_genes)}
    for g in chosen_bulk_genes:
        g_up = str(g).upper()
        if g_up in chosen_ref_map:
            r_idx = chosen_ref_map[g_up]
            if r_idx not in seen_ref_idx:
                seen_ref_idx.add(r_idx)
                shared_genes.append(g)
                ref_sub_idx.append(r_idx)
                bulk_sub_idx.append(bulk_gene_to_idx[g])

    if len(shared_genes) < 50:
        return Failure(
            f"Insufficient gene overlap ({len(shared_genes)}) between reference and cohort '{target_cohort_name}'. "
            f"Overlap counts: ens-ens={c_ens_ens}, sym-sym={c_sym_sym}, sym-ens={c_sym_ens}, ens-sym={c_ens_sym}"
        )

    sub_ref = ref_mat[:, ref_sub_idx]
    sub_bulk = bulk_mat[:, bulk_sub_idx]

    r_sums = sub_ref.sum(axis=1, keepdims=True)
    r_sums[r_sums == 0] = 1.0
    norm_ref = sub_ref / r_sums

    return Success(
        TrialReferenceContext(
            norm_ref=norm_ref,
            sub_bulk=sub_bulk,
            state_labels=tuple(state_labels),
            ref_genes=tuple(ref_genes),
            collinearity_max=collinearity_max,
            condition_number=condition_number,
            total_cells_sampled=sampled_data.total_cells,
        )
    )


def evaluate_single_cohort_rung_with_context(
    ref_ctx: TrialReferenceContext,
    trial_id: int,
    y: np.ndarray,
    target_cohort_name: str = "target",
    rung: FidelityRung = FidelityRung.RUNG_0_SCREENING,
    initial_fractions: np.ndarray | None = None,
    early_prune_threshold: float = 0.52,
    n_jobs: int = 32,
) -> Result[tuple[TrialEvaluationResult, np.ndarray], str]:
    """Evaluate a fidelity rung using precomputed reference profile and warm-started deconvolution."""
    t0 = time.time()

    match rung:
        case FidelityRung.RUNG_0_SCREENING:
            iter_count = 20
            n_repeats = 1
        case FidelityRung.RUNG_1_REFINEMENT:
            iter_count = 20 if initial_fractions is not None else 40
            n_repeats = 2
        case FidelityRung.RUNG_2_FULL:
            iter_count = 20 if initial_fractions is not None else 60
            n_repeats = 3

    # Warm-started or cold-start batch Richardson-Lucy deconvolution
    frac_mat = batch_richardson_lucy_deconvolution(
        bulk=ref_ctx.sub_bulk,
        reference=ref_ctx.norm_ref,
        n_iter=iter_count,
        initial_theta=initial_fractions,
    )

    inner_opt_maybe: Maybe[InnerOptimizationResult] = Nothing

    match rung:
        case FidelityRung.RUNG_0_SCREENING:
            # Fast screening with default CLR Logistic Regression
            default_cfg = InnerClassifierConfig(
                transform=TransformType.CLR,
                selector=FeatureSelectorType.NONE,
                classifier=ClassifierType.LOGISTIC_REGRESSION,
                c_regularization=1.0,
                n_jobs=n_jobs,
            )
            cv_res = evaluate_single_cohort_stratified_cv(
                frac_mat=frac_mat,
                y=y,
                feature_names=ref_ctx.state_labels,
                config=default_cfg,
                n_splits=5,
                n_repeats=n_repeats,
                seed=trial_id * 100 + 42,
            )
            match cv_res:
                case Failure(err):
                    return Failure(f"Rung 0 CV failed: {err}")
                case Success((mean_auc, mean_pr, sel_feats)):
                    pass

        case FidelityRung.RUNG_1_REFINEMENT:
            # Screening transforms and feature selectors with ElasticNet
            refine_cfgs = (
                InnerClassifierConfig(
                    transform=TransformType.CLR,
                    selector=FeatureSelectorType.ELASTICNET,
                    classifier=ClassifierType.LOGISTIC_REGRESSION,
                    c_regularization=1.0,
                    l1_ratio=0.5,
                    n_jobs=n_jobs,
                ),
                InnerClassifierConfig(
                    transform=TransformType.LOGIT,
                    selector=FeatureSelectorType.ELASTICNET,
                    classifier=ClassifierType.LOGISTIC_REGRESSION,
                    c_regularization=1.0,
                    l1_ratio=0.5,
                    n_jobs=n_jobs,
                ),
                InnerClassifierConfig(
                    transform=TransformType.CLR_RATIOS,
                    selector=FeatureSelectorType.ELASTICNET,
                    classifier=ClassifierType.LOGISTIC_REGRESSION,
                    c_regularization=1.0,
                    l1_ratio=0.5,
                    n_jobs=n_jobs,
                ),
            )
            opt_res = run_single_cohort_classifier_hpo(
                frac_mat=frac_mat,
                y=y,
                feature_names=ref_ctx.state_labels,
                candidate_configs=refine_cfgs,
                n_splits=5,
                n_repeats=n_repeats,
                seed=trial_id * 100 + 42,
                n_jobs=n_jobs,
            )
            match opt_res:
                case Failure(err):
                    return Failure(f"Rung 1 inner HPO failed: {err}")
                case Success(opt_val):
                    inner_opt_maybe = Some(opt_val)
                    mean_auc = opt_val.best_loco_auc
                    mean_pr = opt_val.best_loco_pr_auc

        case FidelityRung.RUNG_2_FULL:
            # Full classifier HPO across all candidate transforms, feature selectors, and classifiers
            opt_res = run_single_cohort_classifier_hpo(
                frac_mat=frac_mat,
                y=y,
                feature_names=ref_ctx.state_labels,
                candidate_configs=None,
                n_splits=5,
                n_repeats=n_repeats,
                seed=trial_id * 100 + 42,
                n_jobs=n_jobs,
            )
            match opt_res:
                case Failure(err):
                    return Failure(f"Rung 2 inner HPO failed: {err}")
                case Success(opt_val):
                    inner_opt_maybe = Some(opt_val)
                    mean_auc = opt_val.best_loco_auc
                    mean_pr = opt_val.best_loco_pr_auc

    elapsed = round(time.time() - t0, 3)

    is_pruned = False
    prune_reason: Maybe[str] = Nothing

    if rung == FidelityRung.RUNG_0_SCREENING and mean_auc < early_prune_threshold:
        is_pruned = True
        prune_reason = Some(f"Screening AUC ({mean_auc:.3f}) below threshold {early_prune_threshold:.3f}")
    elif rung == FidelityRung.RUNG_1_REFINEMENT and mean_auc < (early_prune_threshold + 0.01):
        is_pruned = True
        prune_reason = Some(f"Refinement AUC ({mean_auc:.3f}) below threshold {early_prune_threshold + 0.01:.3f}")
    elif ref_ctx.condition_number > 5e5:
        is_pruned = True
        prune_reason = Some(f"Severe reference condition number: {ref_ctx.condition_number:.1e}")

    result = TrialEvaluationResult(
        trial_id=trial_id,
        rung=rung,
        mean_loco_auc=mean_auc,
        mean_loco_pr_auc=mean_pr,
        cohort_aucs={target_cohort_name: mean_auc},
        collinearity_max=ref_ctx.collinearity_max,
        condition_number=ref_ctx.condition_number,
        n_cell_states=len(ref_ctx.state_labels),
        n_shared_genes=len(ref_ctx.ref_genes),
        total_cells_sampled=ref_ctx.total_cells_sampled,
        elapsed_seconds=elapsed,
        inner_optimization=inner_opt_maybe,
        is_pruned=is_pruned,
        prune_reason=prune_reason,
    )

    return Success((result, frac_mat))


def evaluate_single_cohort_trial_with_fidelity(
    config: HPOTrialConfig,
    bulk_mat: np.ndarray,
    bulk_genes: Sequence[str],
    y: np.ndarray,
    target_cohort_name: str = "target",
    bulk_ensembl_genes: Sequence[str] | None = None,
    bulk_symbol_genes: Sequence[str] | None = None,
    rung: FidelityRung = FidelityRung.RUNG_0_SCREENING,
    early_prune_threshold: float = 0.52,
    n_jobs: int = 32,
) -> Result[TrialEvaluationResult, str]:
    """Evaluate a single trial configuration on a single target clinical cohort under specified fidelity rung."""
    prep_res = prepare_trial_reference_context(
        config=config,
        bulk_mat=bulk_mat,
        bulk_genes=bulk_genes,
        target_cohort_name=target_cohort_name,
        bulk_ensembl_genes=bulk_ensembl_genes,
        bulk_symbol_genes=bulk_symbol_genes,
    )
    match prep_res:
        case Failure(err):
            return Failure(err)
        case Success(ref_ctx):
            pass

    eval_res = evaluate_single_cohort_rung_with_context(
        ref_ctx=ref_ctx,
        trial_id=config.trial_id,
        y=y,
        target_cohort_name=target_cohort_name,
        rung=rung,
        initial_fractions=None,
        early_prune_threshold=early_prune_threshold,
        n_jobs=n_jobs,
    )
    match eval_res:
        case Failure(err):
            return Failure(err)
        case Success((eval_val, _)):
            return Success(eval_val)


def run_single_cohort_sampling_hpo(
    target_cohort_name: str,
    target_bulk_mat: np.ndarray,
    target_bulk_genes: Sequence[str],
    target_y: np.ndarray,
    search_space: HPOSearchSpace,
    target_bulk_ensembl_genes: Sequence[str] | None = None,
    target_bulk_symbol_genes: Sequence[str] | None = None,
    n_trials: int = 15,
    seed: int = 42,
    early_prune_threshold: float = 0.52,
    n_jobs: int = 32,
) -> Result[HPORunResult, str]:
    """Execute complete multi-fidelity ASHA hyperparameter optimization campaign for a single target cohort."""
    t_start = time.time()
    rng = np.random.default_rng(seed)

    evaluations: list[TrialEvaluationResult] = []
    configs: dict[int, HPOTrialConfig] = {}
    r0_scores: list[float] = []

    logger.info(
        "Initiating cohort-specific reference sampling HPO for '%s' (%d trials planned, n_jobs=%d)...",
        target_cohort_name,
        n_trials,
        n_jobs,
    )

    for tid in range(1, n_trials + 1):
        cfg = generate_trial_config(tid, search_space, rng)
        configs[tid] = cfg

        # Step 1: Build reference profile ONCE per trial
        prep_res = prepare_trial_reference_context(
            config=cfg,
            bulk_mat=target_bulk_mat,
            bulk_genes=target_bulk_genes,
            target_cohort_name=target_cohort_name,
            bulk_ensembl_genes=target_bulk_ensembl_genes,
            bulk_symbol_genes=target_bulk_symbol_genes,
        )
        match prep_res:
            case Failure(err):
                logger.warning("[Trial %03d] Reference preparation failed: %s", tid, err)
                continue
            case Success(ref_ctx):
                pass

        # Compute dynamic adaptive pruning threshold based on running median
        if len(r0_scores) >= 3:
            dynamic_threshold = max(early_prune_threshold, float(np.median(r0_scores)))
        else:
            dynamic_threshold = early_prune_threshold

        # Step 2: Evaluate at Rung 0 (Screening)
        r0_res = evaluate_single_cohort_rung_with_context(
            ref_ctx=ref_ctx,
            trial_id=tid,
            y=target_y,
            target_cohort_name=target_cohort_name,
            rung=FidelityRung.RUNG_0_SCREENING,
            initial_fractions=None,
            early_prune_threshold=dynamic_threshold,
            n_jobs=n_jobs,
        )

        match r0_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 0 evaluation failed: %s", tid, err)
                continue
            case Success((r0_eval, frac_20)):
                r0_scores.append(r0_eval.mean_loco_auc)
                if r0_eval.is_pruned:
                    logger.info(
                        "[Trial %03d] Pruned at Rung 0: %s (AUC=%.3f vs Cutoff=%.3f)",
                        tid,
                        r0_eval.prune_reason.value_or("below cutoff"),
                        r0_eval.mean_loco_auc,
                        dynamic_threshold,
                    )
                    evaluations.append(r0_eval)
                    continue

        # Step 3: Promote to Rung 1 (Refinement) with warm-started deconvolution
        r1_res = evaluate_single_cohort_rung_with_context(
            ref_ctx=ref_ctx,
            trial_id=tid,
            y=target_y,
            target_cohort_name=target_cohort_name,
            rung=FidelityRung.RUNG_1_REFINEMENT,
            initial_fractions=frac_20,
            early_prune_threshold=dynamic_threshold,
            n_jobs=n_jobs,
        )

        match r1_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 1 evaluation failed: %s", tid, err)
                continue
            case Success((r1_eval, frac_40)):
                if r1_eval.is_pruned:
                    logger.info("[Trial %03d] Pruned at Rung 1: %s", tid, r1_eval.prune_reason.value_or("below cutoff"))
                    evaluations.append(r1_eval)
                    continue

        # Step 4: Promote to Rung 2 (Full Evaluation) with warm-started deconvolution
        r2_res = evaluate_single_cohort_rung_with_context(
            ref_ctx=ref_ctx,
            trial_id=tid,
            y=target_y,
            target_cohort_name=target_cohort_name,
            rung=FidelityRung.RUNG_2_FULL,
            initial_fractions=frac_40,
            early_prune_threshold=dynamic_threshold,
            n_jobs=n_jobs,
        )

        match r2_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 2 evaluation failed: %s", tid, err)
                evaluations.append(r1_eval)
            case Success((r2_eval, _)):
                logger.info(
                    "[Trial %03d] Rung 2 Completed! %s CV-AUC=%.3f, CondNo=%.1e, States=%d",
                    tid,
                    target_cohort_name,
                    r2_eval.mean_loco_auc,
                    r2_eval.condition_number,
                    r2_eval.n_cell_states,
                )
                evaluations.append(r2_eval)

        gc.collect()

    if not evaluations:
        return Failure(f"All HPO trials failed for cohort '{target_cohort_name}'.")

    rows = [
        {
            "trial_id": ev.trial_id,
            "rung": ev.rung.value,
            "mean_loco_auc": ev.mean_loco_auc,  # reuse column name for consistent schema
            "mean_loco_pr_auc": ev.mean_loco_pr_auc,
            "collinearity_max": ev.collinearity_max,
            "condition_number": ev.condition_number,
            "n_cell_states": ev.n_cell_states,
            "n_shared_genes": ev.n_shared_genes,
            "total_cells_sampled": ev.total_cells_sampled,
            "elapsed_seconds": ev.elapsed_seconds,
            "is_pruned": ev.is_pruned,
            "prune_reason": ev.prune_reason.value_or(""),
        }
        for ev in evaluations
    ]
    eval_df = pl.DataFrame(rows)

    non_pruned = [ev for ev in evaluations if not ev.is_pruned]
    best_eval = max(non_pruned, key=lambda x: x.mean_loco_auc) if non_pruned else max(evaluations, key=lambda x: x.mean_loco_auc)
    best_cfg = configs[best_eval.trial_id]

    pareto_trials: list[int] = []
    candidates = non_pruned if non_pruned else evaluations
    for cand in candidates:
        is_dominated = False
        for other in candidates:
            if other.trial_id == cand.trial_id:
                continue
            if (other.mean_loco_auc >= cand.mean_loco_auc and other.collinearity_max <= cand.collinearity_max) and (
                other.mean_loco_auc > cand.mean_loco_auc or other.collinearity_max < cand.collinearity_max
            ):
                is_dominated = True
                break
        if not is_dominated:
            pareto_trials.append(cand.trial_id)

    best_classifier_cfg: Maybe[InnerClassifierConfig] = Nothing
    if isinstance(best_eval.inner_optimization, Some):
        best_classifier_cfg = Some(best_eval.inner_optimization.unwrap().best_config)

    total_time = round(time.time() - t_start, 2)
    n_pruned = sum(1 for ev in evaluations if ev.is_pruned)

    result = HPORunResult(
        best_trial_id=best_eval.trial_id,
        best_config=best_cfg,
        best_loco_auc=best_eval.mean_loco_auc,
        evaluations=tuple(evaluations),
        evaluations_df=eval_df,
        pareto_trials=tuple(sorted(set(pareto_trials))),
        best_classifier_config=best_classifier_cfg,
        total_trials=len(evaluations),
        pruned_trials=n_pruned,
        total_elapsed_seconds=total_time,
    )

    return Success(result)
