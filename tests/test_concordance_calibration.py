"""
Unit Tests for Single-Cell Ground-Truth Calibration and Triangulation Pipeline.
Tests simplex normalization, standardized logistic regression,
diagnostic attribution gating, and Altair SVG rendering.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from returns.result import Success

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "concordance_calibration"))

from importlib import import_module

step02 = import_module("02_pseudobulk_self_deconvolution")
step03 = import_module("03_fit_triangulation_models")
step04 = import_module("04_diagnose_discrepancies")
step05 = import_module("05_plot_triangulation_figures")


def test_pseudobulk_deconv_simplex() -> None:
    # Test that simulated deconvolution fractions sum to 1 per patient
    umi_df = pl.DataFrame({
        "patient": ["p1", "p1", "p2", "p2"],
        "response": ["Responder", "Responder", "Non-Responder", "Non-Responder"],
        "condition": ["Combined", "Combined", "Combined", "Combined"],
        "resolution": [0.5, 0.5, 0.5, 0.5],
        "cell_state": ["StateA", "StateB", "StateA", "StateB"],
        "total_umi": [1000.0, 3000.0, 2000.0, 2000.0],
        "proportion_umi": [0.25, 0.75, 0.5, 0.5],
    })

    deconv_df, fidelity_df = step02.simulate_deconvolution_fractions(umi_df, noise_level=0.01)

    assert deconv_df.height == 4
    for pat in ["p1", "p2"]:
        sub = deconv_df.filter(pl.col("patient") == pat)
        tot = sub["proportion_deconv"].sum()
        assert np.isclose(tot, 1.0, atol=1e-5)


def test_standardized_logistic_fit() -> None:
    # 20 samples, test that standardized beta_z is finite and SE is valid
    df = pl.DataFrame({
        "condition": ["Combined"] * 20,
        "resolution": [0.5] * 20,
        "cell_state": ["State1"] * 20,
        "response": ["Responder"] * 10 + ["Non-Responder"] * 10,
        "val": list(np.linspace(0.1, 0.9, 10)) + list(np.linspace(0.01, 0.2, 10)),
    })

    res_df = step03.fit_standardized_logistic(df, "val", "test")
    assert res_df.height == 1
    row = res_df.to_dicts()[0]
    assert row["beta_test"] > 0.0  # Positively associated with responder
    assert row["se_test"] > 0.0
    assert row["pval_test"] < 0.05


def test_diagnostic_attribution_gates() -> None:
    # Case 1: Concordant
    assert step04.attribute_primary_cause(0.8, 0.85, 0.75, 0.9, 0.7, 0.9) == "Robustly Concordant"

    # Case 2: Deconvolution leakage (low fidelity or deconv flips relative to UMI)
    assert step04.attribute_primary_cause(0.5, 0.6, -0.4, 0.5, -0.3, 0.3) == "Deconvolution Collinear Leakage"

    # Case 3: mRNA mass disparity (UMI flips relative to cell counts)
    assert step04.attribute_primary_cause(0.6, -0.5, -0.4, 0.5, -0.3, 0.8) == "mRNA Mass / Cell Size Disparity"

    # Case 4: Milo local graph vs cluster discrepancy (Milo flips relative to cell counts)
    assert step04.attribute_primary_cause(-0.5, -0.4, -0.4, 0.8, -0.3, 0.8) == "Milo Local Graph vs Cluster Discrepancy"

    # Case 5: Cross-cohort biological heterogeneity (all SC agree, but bulk flips)
    assert step04.attribute_primary_cause(0.6, 0.6, 0.55, 0.7, -0.5, 0.85) == "Cross-Cohort Biological Heterogeneity"


def test_altair_calibration_svg(tmp_path: Path) -> None:
    import altair as alt
    df = pl.DataFrame({"x": [1, 2, 3], "y": [2, 4, 6]})
    chart = alt.Chart(df).mark_line().encode(x="x:Q", y="y:Q")

    out_file = tmp_path / "test_calib.svg"
    res = step05.save_altair_svg(chart, out_file)
    assert isinstance(res, Success)
    assert out_file.exists()
    assert "<svg" in out_file.read_text()
