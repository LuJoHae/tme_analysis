"""
Tests for Synthetic Single-Cell Generation, Milopy DA, and Pseudobulk Deconvolution Benchmark.
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import polars as pl
import pytest
from returns.result import Success

from scripts.synthetic_milo_deconv_benchmark import (
    BenchmarkConfig,
    assign_cells_to_patients,
    classify_quadrant,
    cluster_and_assign_labels,
    evaluate_concordance,
    generate_synthetic_scrna,
    run_pseudobulk_deconvolution,
)


def test_generate_synthetic_scrna() -> None:
    config = BenchmarkConfig(n_cells=600, n_genes=50, n_clusters=3, n_patients=6)
    rng = np.random.default_rng(42)
    adata = generate_synthetic_scrna(config, rng)

    assert adata.n_obs == 600
    assert adata.n_vars == 50
    assert adata.raw is not None
    assert (adata.X.toarray() >= 0).all()


def test_cluster_and_assign_labels() -> None:
    config = BenchmarkConfig(n_cells=600, n_genes=60, n_clusters=3, n_patients=6)
    rng = np.random.default_rng(42)
    adata = generate_synthetic_scrna(config, rng)

    result = cluster_and_assign_labels(
        adata,
        resolution=0.4,
        variable_ratios=False,
        min_prob=0.15,
        max_prob=0.85,
        rng=rng,
    )
    assert isinstance(result, Success)

    adata_clustered, cluster_probs = result.unwrap()
    assert "leiden" in adata_clustered.obs.columns
    assert "cluster_label" in adata_clustered.obs.columns
    assert "cluster_true_prob" in adata_clustered.obs.columns
    assert len(cluster_probs) >= 2


def test_assign_cells_to_patients() -> None:
    config = BenchmarkConfig(n_cells=600, n_genes=50, n_clusters=3, n_patients=6, match_prob=0.9)
    rng = np.random.default_rng(42)
    adata = generate_synthetic_scrna(config, rng)
    adata_clustered, _ = cluster_and_assign_labels(
        adata,
        resolution=0.4,
        variable_ratios=False,
        min_prob=0.1,
        max_prob=0.9,
        rng=rng,
    ).unwrap()

    result = assign_cells_to_patients(adata_clustered, 6, rng)
    assert isinstance(result, Success)

    adata_assigned = result.unwrap()
    assert "patient" in adata_assigned.obs.columns
    assert "response" in adata_assigned.obs.columns
    assert set(adata_assigned.obs["response"].unique()) == {"R", "NR"}
    assert len(adata_assigned.obs["patient"].unique()) == 6


def test_classify_quadrant() -> None:
    assert classify_quadrant(1.5, 0.8, 0.8) == "Concordant Responder"
    assert classify_quadrant(-1.2, -0.4, 0.2) == "Concordant Non-Responder"
    assert classify_quadrant(0.5, -0.3, 0.8) == "Discordant (Milo+, Deconv-)"
    assert classify_quadrant(-0.5, 0.3, 0.2) == "Discordant (Milo-, Deconv+)"
    assert classify_quadrant(0.1, 0.1, 0.5) == "Concordant Neutral"


def test_evaluate_concordance() -> None:
    milo_df = pl.DataFrame({
        "cluster": ["0", "1", "2"],
        "milo_mean_logFC": [1.5, -1.2, 0.0],
        "milo_std_logFC": [0.2, 0.3, 0.1],
        "n_cells": [200, 200, 200],
    })
    deconv_df = pl.DataFrame({
        "cluster": ["0", "1", "2"],
        "deconv_beta": [2.0, -1.5, 0.5],
        "deconv_se": [0.3, 0.4, 0.2],
        "deconv_pval": [0.01, 0.02, 0.5],
        "delta_mean_fraction": [0.15, -0.12, 0.02],
        "mean_fraction_responder": [0.25, 0.05, 0.10],
        "mean_fraction_non_responder": [0.10, 0.17, 0.08],
    })
    probs = {"0": 0.85, "1": 0.15, "2": 0.50}

    result = evaluate_concordance(milo_df, deconv_df, probs)
    assert isinstance(result, Success)

    res_df = result.unwrap()
    assert res_df.height == 3
    assert res_df["quadrant"][0] == "Concordant Responder"
    assert res_df["quadrant"][1] == "Concordant Non-Responder"
    assert bool(res_df["is_concordant"][0]) is True
    assert bool(res_df["is_concordant"][1]) is True


def test_pseudobulk_simplex_constraint() -> None:
    config = BenchmarkConfig(n_cells=600, n_genes=50, n_clusters=3, n_patients=6)
    rng = np.random.default_rng(42)
    adata = generate_synthetic_scrna(config, rng)
    adata_clustered, _ = cluster_and_assign_labels(
        adata,
        resolution=0.4,
        variable_ratios=False,
        min_prob=0.15,
        max_prob=0.85,
        rng=rng,
    ).unwrap()
    adata_assigned = assign_cells_to_patients(adata_clustered, 6, rng).unwrap()

    result = run_pseudobulk_deconvolution(adata_assigned, n_iter=50)
    assert isinstance(result, Success)

    frac_df = result.unwrap()
    sums = frac_df.group_by("patient").agg(pl.col("inferred_fraction").sum().alias("sum_frac"))
    for s in sums["sum_frac"].to_list():
        assert pytest.approx(s, rel=1e-3) == 1.0
