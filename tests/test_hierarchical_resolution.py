"""
Unit Tests for Step 9: Hierarchical BayesPrism and Collinearity Consolidation.
Tests coarse lineage mapping, collinear cluster merging, condition number reduction,
and Altair SVG export.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from returns.result import Success

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "concordance"))

from importlib import import_module

step09 = import_module("09_hierarchical_bayesprism_resolution")


def test_map_coarse_lineage() -> None:
    assert step09.map_coarse_lineage("00_Tem/Trm cytotoxic T cells") == "Cytotoxic_T"
    assert step09.map_coarse_lineage("06_Tem/Temra cytotoxic T cells") == "Cytotoxic_T"
    assert step09.map_coarse_lineage("02_Regulatory T cells") == "Treg"
    assert step09.map_coarse_lineage("04_Naive B cells") == "B_lineage"
    assert step09.map_coarse_lineage("09_Plasma cells") == "B_lineage"
    assert step09.map_coarse_lineage("05_Macrophages") == "Macrophage"
    assert step09.map_coarse_lineage("10_pDC") == "pDC"


def test_consolidate_collinear_reference() -> None:
    # 4 clusters, 100 genes:
    # cluster 0 and 1 belong to Cytotoxic_T with r = 0.99
    # cluster 2 is Macrophage (independent)
    # cluster 3 is B_lineage (independent)
    rng = np.random.default_rng(42)
    g1 = rng.gamma(2.0, 1.0, size=100)
    g2 = g1 + rng.normal(0.0, 0.05, size=100)  # r ~ 0.99
    g3 = rng.gamma(3.0, 1.0, size=100)
    g4 = rng.gamma(1.0, 2.0, size=100)

    ref_mat = np.vstack([g1, g2, g3, g4])
    clusters = [
        "00_Tem/Trm cytotoxic T cells",
        "01_Tem/Trm cytotoxic T cells",
        "05_Macrophages",
        "04_Naive B cells",
    ]

    kappa_before = np.linalg.cond(ref_mat.T)
    cons_mat, cons_names, cluster_map = step09.consolidate_collinear_reference(
        ref_mat, clusters, threshold=0.85
    )
    kappa_after = np.linalg.cond(cons_mat.T)

    # Verifies clusters 0 and 1 merged into 1 consolidated state
    assert len(cons_names) == 3
    assert "Cytotoxic_T_Consolidated" in cons_names[0]
    # Condition number should drop significantly
    assert kappa_after < kappa_before


def test_hierarchical_svg_export(tmp_path: Path) -> None:
    effects_df = pl.DataFrame({
        "evaluation_tier": ["Tier 1: Major Lineage (Hierarchical)"] * 3,
        "entity": ["Cytotoxic_T", "B_lineage", "Macrophage"],
        "beta_sc": [0.5, 1.2, -0.8],
        "beta_bulk": [0.4, 0.9, -0.6],
        "p_value_bulk": [0.05, 0.01, 0.02],
    })
    metrics_df = pl.DataFrame({
        "evaluation_tier": ["Tier 1: Major Lineage (Hierarchical)"],
        "n_entities": [3],
        "spearman_rho": [1.0],
        "spearman_pval": [0.0],
        "pearson_r": [0.99],
        "concordance_rate": [1.0],
    })

    out_file = tmp_path / "test_fig6.svg"
    res = step09.plot_hierarchical_comparison(effects_df, metrics_df, out_file)
    assert isinstance(res, Success)
    assert out_file.exists()
    assert "<svg" in out_file.read_text()
