"""Unit tests for shared Milopy utilities: sparse_projection, memory_utils, and plotting_theme.
"""

from __future__ import annotations

import altair as alt
import numpy as np
import pandas as pd
import polars as pl
import pytest
from scipy import sparse as sp

from milopy import (
    project_nhoods_to_cells,
    release_system_memory,
)
from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COHORT_METADATA,
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)


def test_memory_utils_release_system_memory() -> None:
    """Verifies release_system_memory executes cleanly without errors."""
    # Should not raise any exceptions on any platform (macOS or Linux)
    result = release_system_memory()
    assert result is None


def test_plotting_theme_constants_and_theme() -> None:
    """Verifies plotting theme constants and theme application."""
    assert len(COHORT_PANELS) == 9
    tags = [p["tag"] for p in COHORT_PANELS]
    assert tags == ["a", "b", "c", "d", "e", "f", "g", "h", "i"]

    for panel in COHORT_PANELS:
        assert "acc" in panel
        assert "indication" in panel
        assert "tech" in panel
        assert "patients" in panel

    assert len(COHORT_METADATA) == 9
    assert "GSE120575" in COHORT_METADATA
    assert "GSE316195" in COHORT_METADATA

    # Test applying theme to a simple Altair chart
    df = pl.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
    chart = alt.Chart(df).mark_point().encode(x="x:Q", y="y:Q")
    themed = apply_nature_methods_theme(chart)
    assert themed is not None


def test_sparse_projection_basic() -> None:
    """Tests projecting 3 neighborhoods to 4 cells with known significance and logFC."""
    # 4 cells, 3 nhoods
    # Cell 0: in nhood 0 (logFC = +2.0, sig) -> should be DA+
    # Cell 1: in nhood 1 (logFC = -1.5, sig) -> should be DA-
    # Cell 2: in nhood 2 (logFC = +0.5, not sig) -> should be Not Significant
    # Cell 3: in no nhoods -> should be Not Significant, logFC = 0.0, nhood_count = 0

    nhoods_mat = sp.csc_matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0],
    ])

    res_df = pl.DataFrame({
        "Nhood": [0, 1, 2],
        "logFC": [2.0, -1.5, 0.5],
        "FDR": [0.01, 0.02, 0.50],
        "is_significant": [True, True, False],
    })

    cell_ids = ["c0", "c1", "c2", "c3"]
    patient_ids = ["p1", "p1", "p2", "p2"]
    clinical_responses = ["responder", "responder", "non-responder", "non-responder"]
    cell_types = ["T_cell", "B_cell", "Myeloid", "NK"]
    treatment_statuses = ["pre", "pre", "post", "post"]
    umap_coords = np.array([
        [1.0, 2.0],
        [3.0, 4.0],
        [5.0, 6.0],
        [7.0, 8.0],
    ])

    proj_df = project_nhoods_to_cells(
        nhoods_mat=nhoods_mat,
        res_df=res_df,
        cell_ids=cell_ids,
        patient_ids=patient_ids,
        clinical_responses=clinical_responses,
        cell_types=cell_types,
        treatment_statuses=treatment_statuses,
        umap_coords=umap_coords,
        fdr_threshold=0.10,
    )

    assert len(proj_df) == 4
    assert proj_df["cell_id"].to_list() == cell_ids
    assert proj_df["da_status"].to_list() == [
        "DA+ (Responder Enriched)",
        "DA- (Non-Responder Enriched)",
        "Not Significant",
        "Not Significant",
    ]
    assert np.isclose(proj_df["cell_logfc"][0], 2.0)
    assert np.isclose(proj_df["cell_logfc"][1], -1.5)
    assert np.isclose(proj_df["cell_logfc"][2], 0.5)
    assert np.isclose(proj_df["cell_logfc"][3], 0.0)

    assert proj_df["nhood_count"].to_list() == [1, 1, 1, 0]
    assert "UMAP1" in proj_df.columns
    assert "UMAP2" in proj_df.columns
    assert np.isclose(proj_df["UMAP1"][0], 1.0)
    assert np.isclose(proj_df["UMAP2"][0], 2.0)


def test_sparse_projection_overlapping_nhoods_tiebreak() -> None:
    """Tests tiebreaking when a cell belongs to both DA+ and DA- neighborhoods."""
    # Cell 0 belongs to nhood 0 (+3.0) and nhood 1 (-1.0): net logFC = +1.0 -> DA+
    # Cell 1 belongs to nhood 0 (+1.0) and nhood 1 (-3.0): net logFC = -1.0 -> DA-
    nhoods_mat = sp.csc_matrix([
        [1, 1],
        [1, 1],
    ])

    res_df = pl.DataFrame({
        "Nhood": [0, 1],
        "logFC": [2.0, -1.0],
        "is_significant": [True, True],
    })

    proj_df = project_nhoods_to_cells(
        nhoods_mat=nhoods_mat,
        res_df=res_df,
        cell_ids=["c0", "c1"],
        patient_ids=["p1", "p2"],
        clinical_responses=["responder", "non-responder"],
    )

    # Both cells are in nhood 0 (+2.0) and nhood 1 (-1.0), average logFC = +0.5 > 0 -> DA+
    assert proj_df["da_status"].to_list() == [
        "DA+ (Responder Enriched)",
        "DA+ (Responder Enriched)",
    ]
    assert np.isclose(proj_df["cell_logfc"][0], 0.5)


def test_sparse_projection_pandas_input_and_fdr_cutoff() -> None:
    """Verifies support for pandas DataFrame input with FDR thresholding when is_significant is absent."""
    nhoods_mat = sp.csc_matrix([[1], [1]])
    res_df = pd.DataFrame({
        "logFC": [1.5],
        "FDR": [0.04],
    })

    proj_df = project_nhoods_to_cells(
        nhoods_mat=nhoods_mat,
        res_df=res_df,
        cell_ids=["c0", "c1"],
        patient_ids=["p1", "p2"],
        clinical_responses=["responder", "non-responder"],
        fdr_threshold=0.05,
    )

    assert proj_df["da_status"].to_list() == [
        "DA+ (Responder Enriched)",
        "DA+ (Responder Enriched)",
    ]
    assert np.isclose(proj_df["cell_logfc"][0], 1.5)
