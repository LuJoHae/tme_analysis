"""Property-based testing and invariant verification for tme_response."""

import numpy as np
import polars as pl
import pytest
from returns.maybe import Nothing
from returns.result import Success

from tme_response.evaluation.metrics import calculate_roc_auc
from tme_response.schemas import MultiOmicCohort
from tme_response.signatures.impres import IMPRES_PAIRS, compute_impres_score

try:
    from hypothesis import given, strategies as st
    HAS_HYPOTHESIS = True
except ImportError:
    HAS_HYPOTHESIS = False


def test_impres_bounds_numpy_random() -> None:
    """Invariant: IMPRES score must always be integer in [0, 15] for any arbitrary expression."""
    rng = np.random.default_rng(42)
    all_genes = list({g for pair in IMPRES_PAIRS for g in pair})

    for _ in range(25):
        random_vals = rng.uniform(0.0, 1000.0, size=len(all_genes))
        expr_dict: dict[str, list[object]] = {"sample_id": ["S1"]}
        for idx, g in enumerate(all_genes):
            expr_dict[g] = [float(random_vals[idx])]

        cohort = MultiOmicCohort(
            cohort_id="PropTest",
            cancer_type="Melanoma",
            sample_ids=("S1",),
            expression_tpm=pl.DataFrame(expr_dict),
            clinical_annotations=pl.DataFrame({"sample_id": ["S1"], "response_binary": [1.0]}),
            tmb_scores=Nothing,
            driver_mutations=Nothing,
            cna_scores=Nothing,
        )

        res = compute_impres_score(cohort)
        assert isinstance(res, Success)
        score = res.unwrap().predictions["score"][0]
        assert 0.0 <= score <= 15.0


@pytest.mark.parametrize("scale_factor", [0.1, 0.5, 2.0, 10.0, 100.0])
def test_impres_scale_invariance_param(scale_factor: float) -> None:
    """Invariant: Multiplying expression by any positive scalar strictly preserves IMPRES score."""
    all_genes = list({g for pair in IMPRES_PAIRS for g in pair})
    base_vals = [float(i + 1) for i in range(len(all_genes))]
    scaled_vals = [v * scale_factor for v in base_vals]

    cohort_base = MultiOmicCohort(
        cohort_id="Base",
        cancer_type="Melanoma",
        sample_ids=("S1",),
        expression_tpm=pl.DataFrame({"sample_id": ["S1"], **{g: [base_vals[i]] for i, g in enumerate(all_genes)}}),
        clinical_annotations=pl.DataFrame({"sample_id": ["S1"], "response_binary": [1.0]}),
        tmb_scores=Nothing,
        driver_mutations=Nothing,
        cna_scores=Nothing,
    )

    cohort_scaled = MultiOmicCohort(
        cohort_id="Scaled",
        cancer_type="Melanoma",
        sample_ids=("S1",),
        expression_tpm=pl.DataFrame({"sample_id": ["S1"], **{g: [scaled_vals[i]] for i, g in enumerate(all_genes)}}),
        clinical_annotations=pl.DataFrame({"sample_id": ["S1"], "response_binary": [1.0]}),
        tmb_scores=Nothing,
        driver_mutations=Nothing,
        cna_scores=Nothing,
    )

    res_base = compute_impres_score(cohort_base).unwrap().predictions["score"][0]
    res_scaled = compute_impres_score(cohort_scaled).unwrap().predictions["score"][0]

    assert np.isclose(res_base, res_scaled)


def test_roc_auc_bounds_random() -> None:
    """Invariant: ROC-AUC is always bounded between 0.0 and 1.0."""
    rng = np.random.default_rng(123)
    for _ in range(25):
        y_true = rng.integers(0, 2, size=30)
        y_score = rng.normal(0, 1, size=30)
        auc_val = calculate_roc_auc(y_true, y_score)
        assert 0.0 <= auc_val <= 1.0
