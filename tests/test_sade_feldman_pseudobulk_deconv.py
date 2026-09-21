"""
Unit and property-based tests for Sade-Feldman benchmarking pipeline.
Tests pure functions for cell fraction computation, deconvolution, logistic regression, and concordance metrics.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import polars as pl
import pytest
from returns.result import Success

from scripts.sade_feldman_pseudobulk_deconv_benchmark import (
    classify_quadrant,
    compute_true_cell_fractions,
    evaluate_concordance,
    fit_logistic_stratum,
    run_pseudobulk_deconvolution,
)


def test_compute_true_cell_fractions_simplex() -> None:
    """Verify that compute_true_cell_fractions produces non-negative proportions summing to 1.0."""
    df_mock_cells = pl.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(10)],
        "sample_id": ["S1"] * 6 + ["S2"] * 4,
        "cell_state": ["StateA", "StateA", "StateB", "StateB", "StateB", "StateC", "StateA", "StateC", "StateC", "StateC"],
        "response": ["Responder"] * 6 + ["Non-responder"] * 4,
        "treatment_status": ["Pre"] * 6 + ["Pre"] * 4,
    })

    fractions_df = compute_true_cell_fractions(df_mock_cells)
    assert fractions_df.height == 6  # 2 samples * 3 states

    # Check simplex constraint for each sample
    for sid in ["S1", "S2"]:
        sub = fractions_df.filter(pl.col("sample_id") == sid)
        total_frac = sub["true_fraction"].sum()
        assert pytest.approx(total_frac, abs=1e-6) == 1.0
        assert (sub["true_fraction"] >= 0.0).all()

    # S1 has 2 StateA / 6 = 0.3333, 3 StateB / 6 = 0.5, 1 StateC / 6 = 0.1667
    s1_b = fractions_df.filter((pl.col("sample_id") == "S1") & (pl.col("cell_state") == "StateB"))["true_fraction"][0]
    assert pytest.approx(s1_b, abs=1e-4) == 0.5


def test_classify_quadrant() -> None:
    """Verify directional quadrant classification rules."""
    assert classify_quadrant(1.5, 2.0) == "Concordant Responder"
    assert classify_quadrant(-1.2, -0.8) == "Concordant Non-Responder"
    assert classify_quadrant(1.2, -0.5) == "Discordant (Milo+, Deconv-)"
    assert classify_quadrant(-0.9, 1.1) == "Discordant (Milo-, Deconv+)"
    assert classify_quadrant(0.02, 0.05) == "Concordant Neutral"


def test_fit_logistic_stratum() -> None:
    """Verify standardized effect size estimation in fit_logistic_stratum."""
    # State1 enriched in Responders, State2 enriched in Non-responders
    df_frac = pl.DataFrame({
        "sample_id": ["S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4"],
        "response": ["Responder", "Responder", "Non-responder", "Non-responder"] * 2,
        "cell_state": ["State1"] * 4 + ["State2"] * 4,
        "inferred_fraction": [0.4, 0.5, 0.1, 0.15, 0.05, 0.08, 0.35, 0.42],
    })

    res_df = fit_logistic_stratum(df_frac, fraction_col="inferred_fraction", prefix="deconv")
    assert res_df.height == 2

    s1_row = res_df.filter(pl.col("cell_state") == "State1").to_dicts()[0]
    s2_row = res_df.filter(pl.col("cell_state") == "State2").to_dicts()[0]

    # State1 should have positive beta (enriched in Responder)
    assert s1_row["deconv_beta"] > 0.0
    assert s1_row["deconv_delta_mean"] > 0.0

    # State2 should have negative beta (depleted in Responder)
    assert s2_row["deconv_beta"] < 0.0
    assert s2_row["deconv_delta_mean"] < 0.0


def test_run_pseudobulk_deconvolution_reconstruction() -> None:
    """Verify that run_pseudobulk_deconvolution recovers simplex proportions."""
    clusters = ["StateA", "StateB", "StateC"]
    genes = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7", "Gene8", "Gene9", "Gene10"]

    # Synthetic reference matrix: 3 states x 10 genes
    phi_data: dict[str, object] = {"cluster": clusters}
    np.random.seed(42)
    for g_idx, g in enumerate(genes):
        # Distinct gene signatures
        weights = [10.0 if i == (g_idx % 3) else 0.5 for i in range(3)]
        phi_data[g] = weights

    phi_df = pl.DataFrame(phi_data)

    # Synthetic pseudobulk: 2 samples with known mixture
    sample_ids = ["Sample1", "Sample2"]
    known_fracs = np.array([
        [0.6, 0.3, 0.1],
        [0.1, 0.2, 0.7],
    ])

    ref_sub = phi_df.select(genes).to_numpy()
    norm_phi = ref_sub / ref_sub.sum(axis=1, keepdims=True)
    bulk_mat = (norm_phi.T @ known_fracs.T).T * 1e4  # (2, 10)

    deconv_res = run_pseudobulk_deconvolution(
        sample_ids=sample_ids,
        bulk_genes=genes,
        bulk_mat=bulk_mat,
        phi_df=phi_df,
        marker_genes=genes,
        backend="instaprism",
        n_iter=50,
    )

    assert isinstance(deconv_res, Success)
    df_inf = deconv_res.unwrap()
    assert df_inf.height == 6  # 2 samples * 3 states

    for s_idx, sid in enumerate(sample_ids):
        sub = df_inf.filter(pl.col("sample_id") == sid)
        assert pytest.approx(sub["inferred_fraction"].sum(), abs=1e-4) == 1.0
        inf_arr = sub["inferred_fraction"].to_numpy()
        # High correlation with known fractions
        corr = np.corrcoef(known_fracs[s_idx, :], inf_arr)[0, 1]
        assert corr > 0.90


def test_evaluate_concordance() -> None:
    """Verify evaluation of concordance metrics and summary calculations."""
    clusters = [f"C{i}" for i in range(4)]
    logreg_df = pl.DataFrame({
        "condition": ["Pre"] * 4,
        "cell_state": clusters,
        "deconv_beta": [1.2, 0.8, -1.0, -0.6],
        "deconv_se": [0.3] * 4,
        "deconv_pval": [0.01] * 4,
        "true_beta": [1.1, 0.7, -0.9, -0.5],
        "true_se": [0.3] * 4,
        "true_pval": [0.01] * 4,
        "n_samples": [20] * 4,
    })

    milo_df = pl.DataFrame({
        "condition": ["Pre"] * 4,
        "cell_state": clusters,
        "milo_mean_logfc": [1.5, 0.9, -1.4, -0.7],
        "milo_std_logfc": [0.4] * 4,
        "milo_pval": [0.005] * 4,
        "pct_positive_cells": [0.9, 0.8, 0.1, 0.2],
        "pct_negative_cells": [0.1, 0.2, 0.9, 0.8],
    })

    eval_res = evaluate_concordance(logreg_df, milo_df)
    assert isinstance(eval_res, Success)
    detailed_df, summary_df = eval_res.unwrap()

    assert detailed_df.height == 4
    assert summary_df.height == 1

    sum_row = summary_df.to_dicts()[0]
    # Perfect directional concordance
    assert sum_row["sign_concordance_pct"] == 100.0
    assert sum_row["spearman_rho_milo_deconv"] == 1.0
    assert sum_row["fidelity_spearman_rho"] == 1.0
    assert sum_row["n_concordant"] == 4
