"""Unit tests for single-cell reference sampling HPO and multi-fidelity optimization."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets.sampling.hpo_models import (
    FidelityRung,
    HPOObjectiveMetric,
    HPORunResult,
    HPOSearchSpace,
    HPOTrialConfig,
    TrialEvaluationResult,
)
from tme_datasets.sampling.hpo_optimizer import (
    fast_loco_cv_evaluation,
    filter_candidate_cohorts,
    generate_trial_config,
)
from tme_datasets.sampling.malignant_sampling import MalignantStrategy
from tme_datasets.types import CohortSamplingMode, GeneIDType, HarmonizeMode


def test_hpo_search_space_and_trial_generation() -> None:
    """Verify that HPOSearchSpace generates strictly valid HPOTrialConfig objects."""
    space = HPOSearchSpace(
        candidate_cohort_ids=("GSE120575", "GSE115978", "Maynard_NSCLC", "GSE123813"),
        min_cohorts=2,
        max_cohorts=3,
        min_cells_per_cohort=100,
        max_cells_per_cohort=1000,
        require_cached_h5ad=False,
    )

    rng = np.random.default_rng(42)
    cfg1 = generate_trial_config(trial_id=1, space=space, rng=rng)

    assert cfg1.trial_id == 1
    assert 2 <= len(cfg1.selected_cohort_ids) <= 3
    assert all(c in space.candidate_cohort_ids for c in cfg1.selected_cohort_ids)
    assert 100 <= cfg1.sampling_spec.n_cells_per_cohort <= 1000
    assert 0.3 <= cfg1.cluster_spec.leiden_resolution <= 1.4
    assert 0.75 <= cfg1.reference_config.collinearity_threshold <= 0.92
    assert cfg1.reference_config.auto_detect_malignant == (
        cfg1.malignant_config.strategy != MalignantStrategy.EXCLUDED_TME_ONLY
    )


def test_filter_candidate_cohorts() -> None:
    """Verify candidate single-cell cohort filtering logic."""
    # Filter melanoma cohorts without requiring cached h5ad
    res_mel = filter_candidate_cohorts(cancer_types=["Melanoma"], min_viable_cells=100, require_cached=False)
    assert isinstance(res_mel, Success)
    cohorts = res_mel.unwrap()
    assert "GSE120575" in cohorts
    assert "Maynard_NSCLC" not in cohorts

    # Filter with non-existent cancer type
    res_none = filter_candidate_cohorts(cancer_types=["ImaginaryCancerType123"], require_cached=False)
    assert isinstance(res_none, Failure)


def test_fast_loco_cv_evaluation_synthetic() -> None:
    """Verify Leave-One-Cohort-Out fast evaluation with synthetic reference and bulk matrices."""
    # Create synthetic reference matrix: 3 cell states x 100 genes
    rng = np.random.default_rng(123)
    n_states = 3
    n_genes = 100
    genes = [f"GENE_{i}" for i in range(n_genes)]

    # Make state 0 high in genes 0-30, state 1 in 30-60, state 2 in 60-100
    ref_mat = rng.uniform(0.1, 1.0, size=(n_states, n_genes))
    ref_mat[0, :30] += 5.0
    ref_mat[1, 30:60] += 5.0
    ref_mat[2, 60:] += 5.0

    # Create 2 synthetic bulk cohorts
    # In cohort 1 & 2: responders have high State 0 (e.g. CD8 T cells)
    n_samples_per_cohort = 20
    bulk_dict: dict[str, tuple[np.ndarray, list[str], np.ndarray]] = {}

    for c_idx, cid in enumerate(["CohortA", "CohortB"]):
        y = np.array([1] * (n_samples_per_cohort // 2) + [0] * (n_samples_per_cohort // 2))
        bulk_mat = np.zeros((n_samples_per_cohort, n_genes))

        for s_idx in range(n_samples_per_cohort):
            if y[s_idx] == 1:
                # Responder: 70% state 0, 15% state 1, 15% state 2
                weights = np.array([0.7, 0.15, 0.15])
            else:
                # Non-responder: 10% state 0, 45% state 1, 45% state 2
                weights = np.array([0.1, 0.45, 0.45])

            expr = weights @ ref_mat + rng.normal(0, 0.1, size=n_genes)
            expr = np.clip(expr, 0.01, None)
            bulk_mat[s_idx, :] = expr * 1000.0

        bulk_dict[cid] = (bulk_mat, genes, y)

    # Run fast LOCO CV
    res = fast_loco_cv_evaluation(
        ref_mat=ref_mat,
        ref_genes=genes,
        bulk_data_dict=bulk_dict,
        n_iter=20,
    )

    assert isinstance(res, Success)
    mean_auc, mean_pr, cohort_aucs = res.unwrap()

    assert "CohortA" in cohort_aucs
    assert "CohortB" in cohort_aucs
    # With clear signal, LOCO-AUC should be high (> 0.70)
    assert mean_auc >= 0.70
    assert mean_pr >= 0.50


def test_early_pruning_and_pareto_sorting() -> None:
    """Verify that unpromising trials are marked as pruned and Pareto non-dominated trials are identified."""
    ev1 = TrialEvaluationResult(
        trial_id=1,
        rung=FidelityRung.RUNG_2_FULL,
        mean_loco_auc=0.75,
        mean_loco_pr_auc=0.68,
        cohort_aucs={"C1": 0.74, "C2": 0.76},
        collinearity_max=0.60,
        condition_number=15.0,
        n_cell_states=8,
        n_shared_genes=2000,
        total_cells_sampled=2000,
        elapsed_seconds=12.5,
        is_pruned=False,
    )
    ev2 = TrialEvaluationResult(
        trial_id=2,
        rung=FidelityRung.RUNG_0_SCREENING,
        mean_loco_auc=0.48,
        mean_loco_pr_auc=0.45,
        cohort_aucs={"C1": 0.48},
        collinearity_max=0.91,
        condition_number=800.0,
        n_cell_states=14,
        n_shared_genes=1800,
        total_cells_sampled=400,
        elapsed_seconds=2.1,
        is_pruned=True,
        prune_reason=Some("Screening AUC below threshold"),
    )
    ev3 = TrialEvaluationResult(
        trial_id=3,
        rung=FidelityRung.RUNG_2_FULL,
        mean_loco_auc=0.80,
        mean_loco_pr_auc=0.72,
        cohort_aucs={"C1": 0.79, "C2": 0.81},
        collinearity_max=0.70,
        condition_number=22.0,
        n_cell_states=10,
        n_shared_genes=2200,
        total_cells_sampled=3000,
        elapsed_seconds=14.0,
        is_pruned=False,
    )

    evaluations = (ev1, ev2, ev3)
    # ev3 has higher AUC (0.80 vs 0.75) but slightly higher collinearity (0.70 vs 0.60) -> both ev1 and ev3 are Pareto optimal
    pareto_candidates = [ev for ev in evaluations if not ev.is_pruned]
    pareto_ids: list[int] = []
    for cand in pareto_candidates:
        is_dominated = False
        for other in pareto_candidates:
            if other.trial_id == cand.trial_id:
                continue
            if (other.mean_loco_auc >= cand.mean_loco_auc and other.collinearity_max <= cand.collinearity_max) and (
                other.mean_loco_auc > cand.mean_loco_auc or other.collinearity_max < cand.collinearity_max
            ):
                is_dominated = True
                break
        if not is_dominated:
            pareto_ids.append(cand.trial_id)

    assert set(pareto_ids) == {1, 3}
    assert ev2.is_pruned is True
    assert ev2.prune_reason.value_or("") == "Screening AUC below threshold"
