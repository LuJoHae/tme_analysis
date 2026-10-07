"""Tests for single-cell dataset harmonization (harmonypy) and replicate prevalence filtering in Milopy analysis.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse as sp

from milopy import (
    annotate_nhoods_with_metadata,
    ensure_pca_and_graph,
)


@pytest.fixture
def synthetic_multibatch_adata() -> ad.AnnData:
    """Creates a synthetic multi-patient single-cell AnnData object with batch differences."""
    rng = np.random.default_rng(42)
    n_cells = 300
    n_genes = 100

    # 3 patients: PtA (responder, 100 cells), PtB (responder, 100 cells), PtC (non-responder, 100 cells)
    pts = ["PtA"] * 100 + ["PtB"] * 100 + ["PtC"] * 100
    responses = ["responder"] * 200 + ["non-responder"] * 100
    cell_types = (["CD8_T"] * 50 + ["Macrophage"] * 50) * 3

    # Generate count matrix with batch shifts
    X = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
    # Add batch effect to PtC
    X[200:, :20] += 5.0

    obs = pd.DataFrame(
        {
            "patient_id": pd.Categorical(pts),
            "clinical_response": responses,
            "cell_type": cell_types,
            "treatment_status": ["pre-treatment"] * n_cells,
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_genes)])
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)


def test_ensure_pca_and_graph_harmony(synthetic_multibatch_adata: ad.AnnData) -> None:
    """Verifies that ensure_pca_and_graph integrates multi-patient data using harmonypy."""
    adata = synthetic_multibatch_adata.copy()
    basis = ensure_pca_and_graph(
        adata=adata,
        patient_col="patient_id",
        k=15,
        d=10,
        use_harmony=True,
        seed=42,
    )

    assert basis == "X_pca_harmony"
    assert "X_pca_raw" in adata.obsm
    assert "X_pca_harmony" in adata.obsm
    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape == (300, 10)
    assert adata.obsm["X_pca_harmony"].shape == (300, 10)
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (300, 2)


def test_ensure_pca_and_graph_single_batch(synthetic_multibatch_adata: ad.AnnData) -> None:
    """Verifies that single-patient datasets gracefully skip harmony and return X_pca."""
    adata = synthetic_multibatch_adata[:50].copy()
    adata.obs["patient_id"] = "SinglePt"

    basis = ensure_pca_and_graph(
        adata=adata,
        patient_col="patient_id",
        k=10,
        d=10,
        use_harmony=True,
        seed=42,
    )

    assert basis == "X_pca"
    assert "X_pca" in adata.obsm
    assert "X_pca_harmony" not in adata.obsm


def test_annotate_nhoods_patient_diversity_and_prevalence() -> None:
    """Verifies that annotate_nhoods_with_metadata accurately calculates patient counts and filters private spikes."""
    n_cells = 60
    # 4 patients: 2 responders (R1, R2), 2 non-responders (NR1, NR2)
    patients = ["R1"] * 15 + ["R2"] * 15 + ["NR1"] * 15 + ["NR2"] * 15
    responses = ["responder"] * 30 + ["non-responder"] * 30
    cell_types = ["CD8_T"] * 60

    obs = pd.DataFrame(
        {
            "patient_id": patients,
            "clinical_response": responses,
            "cell_type": cell_types,
        },
        index=[f"c_{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(X=sp.csr_matrix((n_cells, 10)), obs=obs)

    # Construct 3 synthetic neighborhoods:
    # Nhood 0: Cells from NR1 only (Private Spike: 1 non-responder patient)
    # Nhood 1: Cells from NR1 and NR2 (Recurrent non-responder: 2 non-responder patients)
    # Nhood 2: Cells from R1 and R2 (Recurrent responder: 2 responder patients)
    cols = []
    rows = []

    # Nhood 0: cells 30..44 (all NR1)
    for r in range(30, 45):
        rows.append(r)
        cols.append(0)

    # Nhood 1: cells 30..35 (NR1) and 45..50 (NR2)
    for r in list(range(30, 36)) + list(range(45, 51)):
        rows.append(r)
        cols.append(1)

    # Nhood 2: cells 0..10 (R1) and 15..25 (R2)
    for r in list(range(0, 11)) + list(range(15, 26)):
        rows.append(r)
        cols.append(2)

    nhoods_mat = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cells, 3)).tocsc()
    adata.obsm["nhoods"] = nhoods_mat

    # Raw results dataframe with FDR < 0.10 for all 3
    res_df = pd.DataFrame(
        {
            "logFC": [-3.5, -3.2, +2.8],
            "FDR": [0.01, 0.02, 0.03],
            "PValue": [0.001, 0.002, 0.003],
        },
        index=[0, 1, 2],
    )

    annotated = annotate_nhoods_with_metadata(
        adata=adata,
        res_df=res_df,
        patient_col="patient_id",
        design_col="clinical_response",
        min_replicates=2,
        fdr_threshold=0.10,
    )

    # Check patient counts
    assert annotated.loc[0, "n_patients_total"] == 1
    assert annotated.loc[0, "n_patients_non_responder"] == 1
    assert annotated.loc[0, "n_patients_responder"] == 0

    assert annotated.loc[1, "n_patients_total"] == 2
    assert annotated.loc[1, "n_patients_non_responder"] == 2
    assert annotated.loc[1, "n_patients_responder"] == 0

    assert annotated.loc[2, "n_patients_total"] == 2
    assert annotated.loc[2, "n_patients_non_responder"] == 0
    assert annotated.loc[2, "n_patients_responder"] == 2

    # Check prevalence filter results:
    # Nhood 0 is a private spike: only 1 non-responder patient -> should NOT be significant!
    assert annotated.loc[0, "status"] == "Private Clonal Spike"
    assert not annotated.loc[0, "is_significant"]

    # Nhood 1 has 2 non-responder patients -> Enriched in Non-Responders
    assert annotated.loc[1, "status"] == "Enriched in Non-Responders"
    assert annotated.loc[1, "is_significant"]

    # Nhood 2 has 2 responder patients -> Enriched in Responders
    assert annotated.loc[2, "status"] == "Enriched in Responders"
    assert annotated.loc[2, "is_significant"]
