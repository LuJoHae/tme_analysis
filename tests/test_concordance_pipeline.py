"""
Unit and Property-based Tests for the Concordance Pipeline.
Tests mathematical invariants, CLR simplex geometry, purity normalization,
quadrant classification, and SVG generation.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import polars as pl
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success, Failure

# Import functions from scripts
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "concordance"))

from importlib import import_module

step01 = import_module("01_sc_patient_response_association")
step02 = import_module("02_build_curated_reference")
step04 = import_module("04_deconvolve_bulk_bayesprism")
step06 = import_module("06_evaluate_concordance")
step07 = import_module("07_plot_concordance_figures")


def test_binarize_response() -> None:
    assert step01.binarize_response("Responder") == Some(1)
    assert step01.binarize_response("r") == Some(1)
    assert step01.binarize_response("CR") == Some(1)
    assert step01.binarize_response("Non-Responder") == Some(0)
    assert step01.binarize_response("NR") == Some(0)
    assert step01.binarize_response("PD") == Some(0)
    assert step01.binarize_response("Unknown") == Nothing


def test_clr_transformation_invariants() -> None:
    # CLR rows must sum to 0
    rng = np.random.default_rng(42)
    counts = rng.integers(0, 500, size=(10, 5)).astype(np.float64)
    clr = step01.compute_clr_matrix(counts, pseudocount=1e-5)

    assert clr.shape == (10, 5)
    row_sums = clr.sum(axis=1)
    np.testing.assert_allclose(row_sums, np.zeros(10), atol=1e-10)


def test_filter_uninformative_genes() -> None:
    genes = ["CD8A", "MT-CO1", "RPS6", "RPL13A", "HSPA1A", "CXCL13", "IFNG"]
    mask = step02.filter_uninformative_genes(genes)
    kept = [g for g, k in zip(genes, mask, strict=True) if k]
    assert kept == ["CD8A", "CXCL13", "IFNG"]


def test_collinearity_detection() -> None:
    # Create two identical profiles and one orthogonal profile
    profile_a = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    profile_b = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # identical to a
    profile_c = np.array([5.0, 4.0, 1.0, 2.0, 0.0])

    profiles = np.vstack([profile_a, profile_b, profile_c])
    states = ["StateA", "StateB", "StateC"]

    collin_df = step02.evaluate_collinearity(profiles, states, threshold=0.85)
    assert collin_df.height == 3

    # Pair A-B must have correlation ~ 1.0 and is_collinear == True
    pair_ab = collin_df.filter((pl.col("state_a") == "StateA") & (pl.col("state_b") == "StateB"))
    assert pair_ab["is_collinear"][0] is True
    assert pair_ab["correlation"][0] > 0.99


def test_tumor_purity_normalization() -> None:
    # 3 samples, 3 states: Malignant, CD8_T, Macrophage
    states = ["Malignant", "CD8_T", "Macrophage"]
    raw_fractions = pl.DataFrame({
        "sample_id": ["s1", "s2", "s3"],
        "Malignant": [0.8, 0.5, 0.2],
        "CD8_T": [0.1, 0.25, 0.4],
        "Macrophage": [0.1, 0.25, 0.4],
    })

    mrna_df = pl.DataFrame({
        "cell_state": states,
        "relative_rna_content": [2.0, 0.5, 1.5],
    })

    norm_df, scaled_df = step04.compute_normalized_fractions(
        raw_fractions,
        states,
        malignant_label="Malignant",
        mrna_df=mrna_df,
    )

    # In norm_df, Malignant must be 0, and CD8_T + Macrophage must sum to 1.0
    for row in norm_df.iter_rows(named=True):
        assert row["Malignant"] == 0.0
        assert np.isclose(row["CD8_T"] + row["Macrophage"], 1.0, atol=1e-6)

    # In scaled_df, all rows must sum to 1.0
    for row in scaled_df.iter_rows(named=True):
        row_sum = sum(row[s] for s in states)
        assert np.isclose(row_sum, 1.0, atol=1e-6)


def test_quadrant_classification() -> None:
    assert step06.classify_quadrant(1.2, 0.8) == "Concordant Responder"
    assert step06.classify_quadrant(-0.5, -1.1) == "Concordant Non-Responder"
    assert step06.classify_quadrant(0.9, -0.4) == "Discordant (SC+, Bulk-)"
    assert step06.classify_quadrant(-0.8, 0.3) == "Discordant (SC-, Bulk+)"
    assert step06.classify_quadrant(0.0, 0.5) == "Neutral"


def test_concordance_metrics_calculation() -> None:
    # Test perfect concordance
    sc_df = pl.DataFrame({
        "cell_state": ["State1", "State2", "State3", "State4"],
        "beta_sc": [1.5, 0.8, -0.7, -1.2],
    })
    bulk_df = pl.DataFrame({
        "cell_state": ["State1", "State2", "State3", "State4"],
        "fraction_type": ["normalized"] * 4,
        "beta_bulk": [2.1, 0.5, -0.4, -0.9],
    })

    summary, details = step06.compute_concordance_for_fraction_type(sc_df, bulk_df, "normalized")
    assert summary["concordance_rate"] == 1.0
    assert summary["n_concordant"] == 4
    assert summary["n_discordant"] == 0
    assert summary["spearman_rho"] == 1.0
    assert summary["binomial_pval"] < 0.1


def test_altair_svg_export(tmp_path: Path) -> None:
    df = pl.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
    import altair as alt

    chart = alt.Chart(df).mark_point().encode(x="x:Q", y="y:Q")
    out_svg = tmp_path / "test_chart.svg"

    res = step07.save_altair_svg(chart, out_svg)
    assert isinstance(res, Success)
    assert out_svg.exists()
    content = out_svg.read_text()
    assert "<svg" in content
