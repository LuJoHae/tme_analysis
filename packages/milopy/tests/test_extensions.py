"""Unit tests for milopy extensions: permutation, prevalence, composition, projection, memory, and harmony.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
from returns.result import Success
from scipy import sparse as sp

from milopy import (
    NeighborhoodCompositionConfig,
    NeighborhoodCompositionResult,
    PermutationConfig,
    PermutationResult,
    ReplicatePrevalenceConfig,
    annotate_nhoods_with_metadata,
    compute_benjamini_hochberg_fdr,
    compute_neighborhood_composition,
    compute_score_permutation_null,
    ensure_pca_and_graph,
    evaluate_replicate_prevalence,
    project_nhoods_to_cells,
    release_system_memory,
)


def test_milopy_memory_release() -> None:
    # Execution should succeed without error
    release_system_memory()


def test_milopy_benjamini_hochberg_fdr() -> None:
    p_vals = np.array([0.001, 0.01, 0.04, 0.05, 0.50, 0.90])
    q_vals = compute_benjamini_hochberg_fdr(p_vals)
    assert len(q_vals) == len(p_vals)
    assert np.all(q_vals >= 0.0)
    assert np.all(q_vals <= 1.0)
    assert np.all(np.diff(q_vals) >= -1e-12)


def test_milopy_score_permutation_null() -> None:
    # 2 neighborhoods, 20 patients (10 R, 10 NR)
    rng = np.random.default_rng(42)
    counts = np.zeros((2, 20), dtype=int)
    # Nhood 0: Strong recurrent signal across all 10 responders
    counts[0, :10] = 30
    counts[0, 10:] = 0
    # Nhood 1: Balanced null counts
    counts[1, :] = 10
    labels = ["responder"] * 10 + ["non-responder"] * 10

    res = compute_score_permutation_null(
        count_matrix=counts,
        response_labels=labels,
        library_sizes=np.full(20, 1000.0),
        config=PermutationConfig(n_permutations=200, seed=42),
    )
    assert isinstance(res, Success)
    pres: PermutationResult = res.unwrap()
    assert pres.n_nhoods == 2
    assert pres.permutation_pvalues[0] < 0.05
    assert pres.permutation_pvalues[1] > 0.20


def test_milopy_replicate_prevalence() -> None:
    # 6 cells, 2 nhoods
    nhoods = sp.csc_matrix([
        [1, 0],
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
        [0, 1],
    ])
    obs = pd.DataFrame({
        "patient": ["P1", "P2", "P3", "P4", "P4", "P4"],  # Nhood 0 has 3 pts; Nhood 1 has only 1 pt (P4)
        "response": ["responder", "responder", "responder", "non-responder", "non-responder", "non-responder"],
    })
    adata = ad.AnnData(X=np.zeros((6, 10)), obs=obs)
    adata.obsm["nhoods"] = nhoods

    res_df = pd.DataFrame({
        "logFC": [2.0, -2.0],
        "FDR": [0.01, 0.01],
    })

    res = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col="patient",
        design_col="response",
        config=ReplicatePrevalenceConfig(min_replicates=2),
    )
    assert isinstance(res, Success)
    out_df = res.unwrap()
    assert out_df.loc[0, "status"] == "Enriched in Responders"
    assert out_df.loc[0, "is_significant"] == True
    # Nhood 1 has only 1 donor (P4), so it must be flagged as a private clonal spike
    assert out_df.loc[1, "status"] == "Private Clonal Spike"
    assert out_df.loc[1, "is_significant"] == False


def test_milopy_neighborhood_composition() -> None:
    nhoods = sp.csc_matrix([
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
    ])
    obs = pd.DataFrame({
        "dataset_id": ["DS1", "DS1", "DS1", "DS2"],
    })
    adata = ad.AnnData(X=np.zeros((4, 10)), obs=obs)
    adata.obsm["nhoods"] = nhoods

    res = compute_neighborhood_composition(adata, category_col="dataset_id")
    assert isinstance(res, Success)
    comp: NeighborhoodCompositionResult = res.unwrap()
    df = comp.metrics_df
    assert len(df) == 2
    assert df["dataset_id_purity"][0] == 1.0  # pure DS1
    assert df["dataset_id_purity"][1] == 0.5  # 50% DS1, 50% DS2


def test_milopy_sparse_projection() -> None:
    nhoods = sp.csc_matrix([
        [1, 0],
        [0, 1],
    ])
    res_df = pl.DataFrame({
        "Nhood": [0, 1],
        "logFC": [2.5, -1.8],
        "FDR": [0.01, 0.02],
        "is_significant": [True, True],
    })
    cell_ids = ["c1", "c2"]
    patient_ids = ["p1", "p2"]
    responses = ["responder", "non-responder"]

    proj_df = project_nhoods_to_cells(
        nhoods_mat=nhoods,
        res_df=res_df,
        cell_ids=cell_ids,
        patient_ids=patient_ids,
        clinical_responses=responses,
    )
    assert len(proj_df) == 2
    assert proj_df["da_status"][0] == "DA+ (Responder Enriched)"
    assert proj_df["da_status"][1] == "DA- (Non-Responder Enriched)"
