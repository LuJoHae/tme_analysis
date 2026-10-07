"""Tests for vectorized Rao Score permutation testing and null calibration in Milopy.
"""

from __future__ import annotations

import numpy as np
import pytest
from returns.result import Success

from milopy import (
    PermutationConfig,
    compute_benjamini_hochberg_fdr,
    compute_score_permutation_null,
)


def test_benjamini_hochberg_fdr() -> None:
    """Verifies that Benjamini-Hochberg FDR maintains monotonicity and bounded range."""
    p_vals = np.array([0.001, 0.01, 0.04, 0.05, 0.50, 0.90])
    q_vals = compute_benjamini_hochberg_fdr(p_vals)

    assert len(q_vals) == len(p_vals)
    assert np.all(q_vals >= 0.0)
    assert np.all(q_vals <= 1.0)
    # Check monotonicity
    assert np.all(np.diff(q_vals) >= -1e-12)
    assert q_vals[0] < q_vals[-1]


def test_permutation_null_deflates_single_patient_spike() -> None:
    """Verifies that an extreme count spike in a single patient cannot achieve significance under permutation."""
    rng = np.random.default_rng(42)
    n_nhoods = 100
    n_patients = 17  # 15 R, 2 NR (matching PDAC GSE316195)
    resp_labels = ["responder"] * 15 + ["non-responder"] * 2

    # Baseline random counts
    counts = rng.poisson(lam=5.0, size=(n_nhoods, n_patients))

    # Nhood 0: Extreme single-patient spike in patient 15 (NR1)
    # 0 in all 15 responders, 50 in NR1, 0 in NR2
    counts[0, :15] = 0
    counts[0, 15] = 50
    counts[0, 16] = 0

    # Nhood 1: Recurrent signal across BOTH non-responders
    # 0 in all 15 responders, 25 in NR1, 25 in NR2
    counts[1, :15] = 0
    counts[1, 15] = 25
    counts[1, 16] = 25

    res = compute_score_permutation_null(
        count_matrix=counts,
        response_labels=resp_labels,
        config=PermutationConfig(n_permutations=500, seed=42, fdr_threshold=0.10),
    )

    assert isinstance(res, Success)
    p_out = res.unwrap()

    p_spike = p_out.permutation_pvalues[0]
    p_recurrent = p_out.permutation_pvalues[1]

    # For a single patient spike with 2 NR out of 17, whenever NR1 is randomly assigned
    # to NR, the score statistic is identically large. With sampling and library size variation,
    # p_spike is heavily deflated (>= 0.04 vs ~1e-5 in parametric tests) and FDR is deflated to > 0.80.
    assert p_spike >= 0.04, f"Single-patient spike p-value ({p_spike}) must be deflated!"
    assert p_recurrent < p_spike, f"Recurrent hit ({p_recurrent}) must have smaller p-value than spike ({p_spike})!"

    # Under 100 tests, a p-value >= 0.04 cannot pass FDR < 0.10
    fdr_spike = p_out.permutation_fdr[0]
    assert fdr_spike > 0.10, f"Single-patient spike FDR ({fdr_spike}) must not pass significance threshold!"


def test_permutation_null_identifies_true_recurrent_signal() -> None:
    """Verifies that a true signal shared across multiple donors achieves significant permutation p-value."""
    n_nhoods = 50
    n_patients = 20  # 10 R, 10 NR
    resp_labels = ["responder"] * 10 + ["non-responder"] * 10

    counts = np.random.default_rng(42).poisson(lam=5.0, size=(n_nhoods, n_patients))

    # Nhood 0: Clean recurrent signal present in all 10 responders and 0 non-responders
    counts[0, :10] = 30
    counts[0, 10:] = 0

    res = compute_score_permutation_null(
        count_matrix=counts,
        response_labels=resp_labels,
        config=PermutationConfig(n_permutations=500, seed=42, fdr_threshold=0.10),
    )

    assert isinstance(res, Success)
    p_out = res.unwrap()

    assert p_out.permutation_pvalues[0] < 0.01
    assert p_out.permutation_fdr[0] < 0.05
    assert p_out.is_perm_significant[0]
