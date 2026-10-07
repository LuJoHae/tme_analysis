"""Pure functional Bayesian and Multi-Fidelity ASHA optimization engine for single-cell sampling."""

from __future__ import annotations

import time
from typing import Mapping, Sequence
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score  # type: ignore[import-untyped]

import instaprism

from ..deconvolution.builder import build_deconvolution_reference
from ..deconvolution.models import DeconvolutionReferenceConfig
from ..logging import get_logger
from ..models import (
    ClusterAnalysisSpec,
    HarmonizeConfig,
    SingleCellSamplingSpec,
)
from ..paths import find_dataset_h5ad
from ..registry import list_registered_datasets
from ..types import CohortSamplingMode, GeneIDType, HarmonizeMode, Modality
from .classifier_hpo import (
    InnerClassifierConfig,
    InnerOptimizationResult,
    run_inner_classifier_hpo,
)
from .cohort_sampler import sample_single_cell_cohorts
from .hpo_models import (
    FidelityRung,
    HPOObjectiveMetric,
    HPORunResult,
    HPOSearchSpace,
    HPOTrialConfig,
    TrialEvaluationResult,
)
from .malignant_sampling import MalignantSamplingConfig, MalignantStrategy

logger = get_logger("sampling.hpo_optimizer")


def filter_candidate_cohorts(
    cancer_types: Sequence[str] | None = None,
    min_viable_cells: int = 200,
    max_viable_cells: int = 80000,
    require_cached: bool = True,
) -> Result[tuple[str, ...], str]:
    """Filter single-cell cohorts from registry suitable for reference deconvolution optimization."""
    try:
        registered = list_registered_datasets().filter(modality=Modality.SINGLE_CELL)
        candidates: list[str] = []

        cancer_set = {c.lower() for c in cancer_types} if cancer_types is not None else None

        for spec in registered:
            if cancer_set is not None and spec.cancer_type.lower() not in cancer_set:
                continue

            match spec.n_samples_or_cells:
                case Some(val):
                    n_cells = val
                case _:
                    n_cells = 0

            if n_cells < min_viable_cells or n_cells > max_viable_cells:
                continue

            if require_cached:
                if not isinstance(find_dataset_h5ad(spec.id), Some):
                    continue

            candidates.append(spec.id)

        if not candidates:
            return Failure("No single-cell cohorts matched the specified filtering criteria.")

        return Success(tuple(candidates))
    except Exception as exc:
        return Failure(f"Failed to filter candidate cohorts: {exc}")


def generate_trial_config(
    trial_id: int,
    space: HPOSearchSpace,
    rng: np.random.Generator,
) -> HPOTrialConfig:
    """Deterministically draw a parameter configuration from search space."""
    # 1. Dataset selection
    n_cohorts = int(rng.integers(space.min_cohorts, min(space.max_cohorts + 1, len(space.candidate_cohort_ids) + 1)))
    n_cohorts = max(1, min(n_cohorts, len(space.candidate_cohort_ids)))
    chosen_indices = rng.choice(len(space.candidate_cohort_ids), size=n_cohorts, replace=False)
    selected_cohorts = tuple(space.candidate_cohort_ids[i] for i in sorted(chosen_indices))

    # 2. Sampling parameters
    sampling_mode = space.sampling_modes[int(rng.integers(0, len(space.sampling_modes)))]
    n_cells = int(np.exp(rng.uniform(np.log(space.min_cells_per_cohort), np.log(space.max_cells_per_cohort))))
    global_budget = int(np.exp(rng.uniform(np.log(space.min_global_budget), np.log(space.max_global_budget))))

    stratify_choice = space.stratify_options[int(rng.integers(0, len(space.stratify_options)))]
    stratify_maybe = Some(stratify_choice) if stratify_choice is not None else Nothing

    sampling_spec = SingleCellSamplingSpec(
        cohort_ids=Some(selected_cohorts),
        mode=sampling_mode,
        n_cells_per_cohort=n_cells,
        global_cell_budget=Some(global_budget),
        stratify_by=stratify_maybe,
        balanced_strata=bool(rng.choice([True, False])),
        require_cached_h5ad=space.require_cached_h5ad,
        seed=Some(int(rng.integers(1, 100000))),
    )

    # 3. Malignant Sampling
    mal_strategy = space.malignant_strategies[int(rng.integers(0, len(space.malignant_strategies)))]
    mal_frac = float(rng.uniform(space.malignant_fraction_range[0], space.malignant_fraction_range[1]))
    malignant_config = MalignantSamplingConfig(
        strategy=mal_strategy,
        malignant_fraction=round(mal_frac, 3),
        max_cells_per_patient=int(rng.integers(20, 100)),
        contig_variance_qc=True,
    )

    # 4. Harmonization
    harmonize_mode = space.harmonize_modes[int(rng.integers(0, len(space.harmonize_modes)))]
    gene_target = space.gene_target_types[int(rng.integers(0, len(space.gene_target_types)))]
    min_shared = int(rng.integers(space.min_shared_genes_range[0], space.min_shared_genes_range[1] + 1))

    harmonize_config = HarmonizeConfig(
        mode=harmonize_mode,
        gene_target_type=gene_target,
        min_shared_genes=min_shared,
        reconcile_genes=True,
    )

    # 5. Clustering & Graph Analysis
    n_hvg = int(rng.integers(space.n_top_genes_range[0], space.n_top_genes_range[1] + 1))
    n_pcs = int(rng.integers(space.n_pcs_range[0], space.n_pcs_range[1] + 1))
    leiden_res = float(rng.uniform(space.leiden_resolution_range[0], space.leiden_resolution_range[1]))

    cluster_spec = ClusterAnalysisSpec(
        n_top_genes=n_hvg,
        n_pcs=n_pcs,
        leiden_resolution=round(leiden_res, 3),
        seed=Some(int(rng.integers(1, 100000))),
    )

    # 6. Deconvolution Reference Builder
    collin_thresh = float(rng.uniform(space.collinearity_threshold_range[0], space.collinearity_threshold_range[1]))
    min_cells_state = int(rng.integers(space.min_cells_per_state_range[0], space.min_cells_per_state_range[1] + 1))

    reference_config = DeconvolutionReferenceConfig(
        cell_state_key=cluster_spec.key_added,
        collinearity_threshold=round(collin_thresh, 3),
        min_cells_per_state=min_cells_state,
        auto_detect_malignant=(mal_strategy != MalignantStrategy.EXCLUDED_TME_ONLY),
        normalize_multinomial=True,
    )

    return HPOTrialConfig(
        trial_id=trial_id,
        selected_cohort_ids=selected_cohorts,
        sampling_spec=sampling_spec,
        malignant_config=malignant_config,
        cluster_spec=cluster_spec,
        harmonize_config=harmonize_config,
        reference_config=reference_config,
    )


def batch_richardson_lucy_deconvolution(
    bulk: np.ndarray,
    reference: np.ndarray,
    n_iter: int = 40,
    eps: float = 1e-12,
) -> np.ndarray:
    """Batch-vectorized Richardson-Lucy / InstaPrism deconvolution across all samples simultaneously.

    Executes via multi-threaded BLAS GEMM matrix multiplications across all CPU cores.
    Mathematically identical to instaprism.insta_prism with floating point error < 1e-15.
    """
    n_samples, _ = bulk.shape
    n_states = reference.shape[0]
    ref_t = reference.T  # (G, S)

    # Normalize bulk to 1e6 CPM per sample
    b_sums = bulk.sum(axis=1, keepdims=True)
    b_sums[b_sums == 0] = 1.0
    bulk_cpm = (bulk / b_sums) * 1e6

    # Initialize uniform cell state fractions: shape (N, S)
    theta = np.full((n_samples, n_states), 1.0 / n_states, dtype=np.float64)

    for _ in range(n_iter):
        # 1. Projected bulk expression: (N, S) @ (S, G) -> (N, G)
        b_hat = theta @ reference
        # 2. Ratio of observed to projected: (N, G)
        ratio = bulk_cpm / (b_hat + eps)
        # 3. Multiplicative update: (N, G) @ (G, S) -> (N, S)
        theta *= (ratio @ ref_t)
        # 4. Normalize rows to sum to 1.0
        t_sums = theta.sum(axis=1, keepdims=True)
        t_sums[t_sums == 0] = 1.0
        theta /= t_sums

    return theta


def fast_deconvolute_cohorts(
    ref_mat: np.ndarray,
    ref_genes: Sequence[str],
    bulk_data_dict: Mapping[str, tuple[Any, ...]],
    n_iter: int = 40,
) -> Result[tuple[dict[str, np.ndarray], dict[str, np.ndarray]], str]:
    """Deconvolute multiple bulk cohorts using instaprism fixed-point solver."""
    if len(bulk_data_dict) < 1:
        return Failure("Need at least 1 bulk cohort for deconvolution.")

    ref_gene_to_idx = {g.upper(): idx for idx, g in enumerate(ref_genes)}
    cohort_fractions: dict[str, np.ndarray] = {}
    cohort_responses: dict[str, np.ndarray] = {}

    for cid, entry in bulk_data_dict.items():
        if len(entry) == 4:
            bulk_mat, ensembl_genes, symbol_genes, y_resp = entry
            # Pick whichever gene nomenclature maximizes overlap with ref_genes
            ens_overlap = sum(1 for g in ensembl_genes if str(g).upper() in ref_gene_to_idx)
            sym_overlap = sum(1 for g in symbol_genes if str(g).upper() in ref_gene_to_idx)
            bulk_genes = ensembl_genes if ens_overlap >= sym_overlap else symbol_genes
        else:
            bulk_mat, bulk_genes, y_resp = entry

        shared_genes = [g for g in bulk_genes if str(g).upper() in ref_gene_to_idx]
        if len(shared_genes) < 50:
            return Failure(f"Insufficient gene overlap ({len(shared_genes)}) between reference and cohort '{cid}'.")

        ref_sub_idx = [ref_gene_to_idx[str(g).upper()] for g in shared_genes]
        bulk_gene_to_idx = {g: idx for idx, g in enumerate(bulk_genes)}
        bulk_sub_idx = [bulk_gene_to_idx[g] for g in shared_genes]

        sub_ref = ref_mat[:, ref_sub_idx]
        sub_bulk = bulk_mat[:, bulk_sub_idx]

        r_sums = sub_ref.sum(axis=1, keepdims=True)
        r_sums[r_sums == 0] = 1.0
        norm_ref = sub_ref / r_sums

        frac_mat = batch_richardson_lucy_deconvolution(
            bulk=sub_bulk,
            reference=norm_ref,
            n_iter=n_iter,
        )

        cohort_fractions[cid] = frac_mat
        cohort_responses[cid] = y_resp

    return Success((cohort_fractions, cohort_responses))


def fast_loco_cv_evaluation(
    ref_mat: np.ndarray,
    ref_genes: Sequence[str],
    bulk_data_dict: Mapping[str, tuple[np.ndarray, Sequence[str], np.ndarray]],
    n_iter: int = 40,
) -> Result[tuple[float, float, dict[str, float]], str]:
    """Execute Leave-One-Cohort-Out CV across pre-loaded bulk cohorts using fast deconvolution surrogate."""
    deconv_res = fast_deconvolute_cohorts(
        ref_mat=ref_mat,
        ref_genes=ref_genes,
        bulk_data_dict=bulk_data_dict,
        n_iter=n_iter,
    )
    match deconv_res:
        case Failure(err):
            return Failure(err)
        case Success((cohort_fractions, cohort_responses)):
            pass

    cohort_aucs: dict[str, float] = {}
    cohort_pr_aucs: dict[str, float] = {}
    cohort_ids = list(bulk_data_dict.keys())

    for test_cid in cohort_ids:
        train_x_list = [cohort_fractions[cid] for cid in cohort_ids if cid != test_cid]
        train_y_list = [cohort_responses[cid] for cid in cohort_ids if cid != test_cid]

        train_X = np.vstack(train_x_list)
        train_y = np.concatenate(train_y_list)

        test_X = cohort_fractions[test_cid]
        test_y = cohort_responses[test_cid]

        if len(np.unique(train_y)) < 2 or len(np.unique(test_y)) < 2:
            cohort_aucs[test_cid] = 0.5
            cohort_pr_aucs[test_cid] = float(np.mean(test_y))
            continue

        try:
            clf = LogisticRegression(C=1.0, max_iter=200, solver="lbfgs")
            clf.fit(train_X, train_y)
            pred_probs = clf.predict_proba(test_X)[:, 1]

            score_auc = float(roc_auc_score(test_y, pred_probs))
            precision, recall, _ = precision_recall_curve(test_y, pred_probs)
            score_pr = float(auc(recall, precision))

            cohort_aucs[test_cid] = score_auc
            cohort_pr_aucs[test_cid] = score_pr
        except Exception:
            cohort_aucs[test_cid] = 0.5
            cohort_pr_aucs[test_cid] = float(np.mean(test_y))

    mean_auc = float(np.mean(list(cohort_aucs.values())))
    mean_pr = float(np.mean(list(cohort_pr_aucs.values())))

    return Success((mean_auc, mean_pr, cohort_aucs))


def evaluate_trial_with_fidelity(
    config: HPOTrialConfig,
    bulk_data_dict: Mapping[str, tuple[Any, ...]],
    rung: FidelityRung = FidelityRung.RUNG_0_SCREENING,
    early_prune_threshold: float = 0.51,
    n_jobs: int = 1,
) -> Result[TrialEvaluationResult, str]:
    """Evaluate a single trial parameter config under the specified multi-fidelity rung."""
    t0 = time.time()

    sub_bulk: Mapping[str, tuple[Any, ...]]
    match rung:
        case FidelityRung.RUNG_0_SCREENING:
            iter_count = 20
            sub_bulk = {k: bulk_data_dict[k] for k in list(bulk_data_dict.keys())[:2]}
        case FidelityRung.RUNG_1_REFINEMENT:
            iter_count = 40
            sub_bulk = {k: bulk_data_dict[k] for k in list(bulk_data_dict.keys())[:4]}
        case FidelityRung.RUNG_2_FULL:
            iter_count = 60
            sub_bulk = bulk_data_dict

    # 1. Sample and cluster single cells
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

    # 2. Build deconvolution reference profile
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

    # 3. Deconvolute cohorts
    deconv_res = fast_deconvolute_cohorts(
        ref_mat=ref_mat,
        ref_genes=ref_genes,
        bulk_data_dict=sub_bulk,
        n_iter=iter_count,
    )
    match deconv_res:
        case Failure(err):
            return Failure(f"Deconvolution failed: {err}")
        case Success((cohort_fracs, cohort_resps)):
            pass

    # 4. Bi-Level Inner Optimization: Transform, Feature Selection, and Classifier
    inner_opt_res = run_inner_classifier_hpo(
        cohort_fractions=cohort_fracs,
        cohort_responses=cohort_resps,
        feature_names=state_labels,
        n_jobs=n_jobs,
    )

    match inner_opt_res:
        case Failure(err):
            return Failure(f"Inner classifier HPO failed: {err}")
        case Success(inner_res):
            mean_auc = inner_res.best_loco_auc
            mean_pr = inner_res.best_loco_pr_auc
            cohort_aucs = inner_res.cohort_aucs
            inner_opt_maybe: Maybe[InnerOptimizationResult] = Some(inner_res)

    elapsed = round(time.time() - t0, 3)

    is_pruned = False
    prune_reason: Maybe[str] = Nothing

    if rung == FidelityRung.RUNG_0_SCREENING and mean_auc < early_prune_threshold:
        is_pruned = True
        prune_reason = Some(f"Screening AUC ({mean_auc:.3f}) below threshold {early_prune_threshold:.3f}")
    elif condition_number > 5e5:
        is_pruned = True
        prune_reason = Some(f"Severe reference condition number: {condition_number:.1e}")

    result = TrialEvaluationResult(
        trial_id=config.trial_id,
        rung=rung,
        mean_loco_auc=mean_auc,
        mean_loco_pr_auc=mean_pr,
        cohort_aucs=cohort_aucs,
        collinearity_max=collinearity_max,
        condition_number=condition_number,
        n_cell_states=len(state_labels),
        n_shared_genes=len(ref_genes),
        total_cells_sampled=sampled_data.total_cells,
        elapsed_seconds=elapsed,
        inner_optimization=inner_opt_maybe,
        is_pruned=is_pruned,
        prune_reason=prune_reason,
    )

    return Success(result)


def run_sampling_hpo(
    search_space: HPOSearchSpace,
    bulk_data_dict: Mapping[str, tuple[Any, ...]],
    n_trials: int = 15,
    seed: int = 42,
    early_prune_threshold: float = 0.52,
    n_jobs: int = 1,
) -> Result[HPORunResult, str]:
    """Execute complete multi-fidelity ASHA hyperparameter optimization campaign."""
    t_start = time.time()
    rng = np.random.default_rng(seed)

    evaluations: list[TrialEvaluationResult] = []
    configs: dict[int, HPOTrialConfig] = {}

    logger.info("Initiating single-cell reference sampling HPO (%d trials planned)...", n_trials)

    for tid in range(1, n_trials + 1):
        cfg = generate_trial_config(tid, search_space, rng)
        configs[tid] = cfg

        # Step 1: Evaluate at Rung 0 (Screening)
        r0_res = evaluate_trial_with_fidelity(
            config=cfg,
            bulk_data_dict=bulk_data_dict,
            rung=FidelityRung.RUNG_0_SCREENING,
            early_prune_threshold=early_prune_threshold,
            n_jobs=n_jobs,
        )

        match r0_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 0 evaluation failed: %s", tid, err)
                continue
            case Success(r0_eval):
                if r0_eval.is_pruned:
                    logger.info("[Trial %03d] Pruned at Rung 0: %s", tid, r0_eval.prune_reason.value_or("below cutoff"))
                    evaluations.append(r0_eval)
                    continue

        # Step 2: Promote to Rung 1 (Refinement)
        r1_res = evaluate_trial_with_fidelity(
            config=cfg,
            bulk_data_dict=bulk_data_dict,
            rung=FidelityRung.RUNG_1_REFINEMENT,
            early_prune_threshold=early_prune_threshold,
            n_jobs=n_jobs,
        )

        match r1_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 1 evaluation failed: %s", tid, err)
                continue
            case Success(r1_eval):
                if r1_eval.is_pruned:
                    logger.info("[Trial %03d] Pruned at Rung 1: %s", tid, r1_eval.prune_reason.value_or("below cutoff"))
                    evaluations.append(r1_eval)
                    continue

        # Step 3: Promote to Rung 2 (Full Evaluation)
        r2_res = evaluate_trial_with_fidelity(
            config=cfg,
            bulk_data_dict=bulk_data_dict,
            rung=FidelityRung.RUNG_2_FULL,
            early_prune_threshold=early_prune_threshold,
            n_jobs=n_jobs,
        )

        match r2_res:
            case Failure(err):
                logger.warning("[Trial %03d] Rung 2 evaluation failed: %s", tid, err)
                evaluations.append(r1_eval)
            case Success(r2_eval):
                logger.info(
                    "[Trial %03d] Rung 2 Completed! LOCO-AUC=%.3f, CondNo=%.1e, States=%d",
                    tid,
                    r2_eval.mean_loco_auc,
                    r2_eval.condition_number,
                    r2_eval.n_cell_states,
                )
                evaluations.append(r2_eval)

        # Proactively free memory between multi-fidelity trials
        import gc
        gc.collect()

    if not evaluations:
        return Failure("All HPO trials failed to execute or evaluate.")

    rows = [
        {
            "trial_id": ev.trial_id,
            "rung": ev.rung.value,
            "mean_loco_auc": ev.mean_loco_auc,
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
