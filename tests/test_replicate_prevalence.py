"""Tests for biological replicate prevalence and multi-cohort recurrence filtering in Milopy.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from returns.result import Failure, Success
from scipy import sparse as sp

from milopy import (
    ReplicatePrevalenceConfig,
    evaluate_replicate_prevalence,
)


def create_synthetic_cohort_adata() -> tuple[ad.AnnData, pd.DataFrame]:
    """Builds a minimal synthetic AnnData and Milo result table for prevalence testing."""
    # 20 cells: 10 responder (pts: R1, R2, R3, R4, R5), 10 non-responder (pts: NR1, NR2, NR3)
    patient_ids = [
        "R1", "R1", "R2", "R2", "R3", "R3", "R4", "R4", "R5", "R5",
        "NR1", "NR1", "NR1", "NR1", "NR2", "NR2", "NR3", "NR3", "NR3", "NR3",
    ]
    responses = (
        ["responder"] * 10 + ["non-responder"] * 10
    )
    dataset_ids = (
        ["DS1"] * 5 + ["DS2"] * 5 + ["DS1"] * 5 + ["DS2"] * 5
    )

    obs = pd.DataFrame(
        {
            "patient_id": patient_ids,
            "clinical_response": responses,
            "dataset_id": dataset_ids,
        },
        index=[f"cell_{i}" for i in range(20)],
    )

    # 4 Neighborhoods:
    # Nhood 0: Cells 10, 11, 12, 13 (all NR1) -> Single-donor spike
    # Nhood 1: Cells 10 (NR1), 14 (NR2), 16 (NR3) -> Multi-donor NR hit
    # Nhood 2: Cells 0 (R1), 2 (R2), 4 (R3), 6 (R4) -> Multi-donor R hit
    # Nhood 3: Cells 0, 10 -> Mixed / non-significant
    nhoods = np.zeros((20, 4), dtype=np.int32)
    nhoods[[10, 11, 12, 13], 0] = 1
    nhoods[[10, 14, 16], 1] = 1
    nhoods[[0, 2, 4, 6], 2] = 1
    nhoods[[0, 10], 3] = 1

    adata = ad.AnnData(
        X=np.zeros((20, 10), dtype=np.float32),
        obs=obs,
        obsm={"nhoods": sp.csc_matrix(nhoods)},
    )

    res_df = pd.DataFrame(
        {
            "logFC": [-2.5, -2.0, 2.0, 0.1],
            "FDR": [0.005, 0.01, 0.02, 0.65],
            "dataset_id_purity": [1.0, 0.50, 0.50, 0.50],
        },
        index=[0, 1, 2, 3],
    )

    return adata, res_df


def test_evaluate_replicate_prevalence_filters_single_donor_spike() -> None:
    """Verifies that a single-donor spike is caught and labeled as Private Clonal Spike."""
    adata, res_df = create_synthetic_cohort_adata()

    config = ReplicatePrevalenceConfig(
        min_replicates=2,
        min_prevalence_frac=0.10,
        fdr_threshold=0.10,
    )

    res = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col="patient_id",
        design_col="clinical_response",
        config=config,
    )

    assert isinstance(res, Success)
    out_df = res.unwrap()

    # Nhood 0 has only NR1 (1 patient) -> Private Clonal Spike
    assert out_df.loc[0, "n_patients_non_responder"] == 1
    assert out_df.loc[0, "status"] == "Private Clonal Spike"
    assert bool(out_df.loc[0, "is_significant"]) is False

    # Nhood 1 has NR1, NR2, NR3 (3 patients) -> Enriched in Non-Responders
    assert out_df.loc[1, "n_patients_non_responder"] == 3
    assert out_df.loc[1, "status"] == "Enriched in Non-Responders"
    assert bool(out_df.loc[1, "is_significant"]) is True

    # Nhood 2 has R1, R2, R3, R4 (4 patients) -> Enriched in Responders
    assert out_df.loc[2, "n_patients_responder"] == 4
    assert out_df.loc[2, "status"] == "Enriched in Responders"
    assert bool(out_df.loc[2, "is_significant"]) is True

    # Nhood 3 has FDR = 0.65 -> Not Significant
    assert out_df.loc[3, "status"] == "Not Significant"
    assert bool(out_df.loc[3, "is_significant"]) is False


def test_evaluate_replicate_prevalence_multi_cohort_recurrence() -> None:
    """Verifies that multi-cohort hits are partitioned into Conserved vs Cohort-Private."""
    adata, res_df = create_synthetic_cohort_adata()

    config = ReplicatePrevalenceConfig(
        min_replicates=2,
        min_prevalence_frac=0.10,
        min_datasets=2,
        dataset_purity_threshold=0.85,
        fdr_threshold=0.10,
    )

    res = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col="patient_id",
        design_col="clinical_response",
        dataset_col="dataset_id",
        config=config,
    )

    assert isinstance(res, Success)
    out_df = res.unwrap()

    # Nhood 1 spans DS1 and DS2 with purity 0.50 < 0.85 -> Conserved Recurrent Hit
    assert out_df.loc[1, "hit_type"] == "Conserved Recurrent Hit"

    # Nhood 0 failed replicate prevalence -> Private Clonal Spike
    assert out_df.loc[0, "hit_type"] == "Private Clonal Spike"


def test_evaluate_replicate_prevalence_missing_columns() -> None:
    """Verifies that missing columns or invalid inputs return Failure."""
    adata, res_df = create_synthetic_cohort_adata()

    # Missing patient column
    res_bad_pt = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col="non_existent_patient",
        design_col="clinical_response",
    )
    assert isinstance(res_bad_pt, Failure)
    assert "Patient column 'non_existent_patient' not found" in res_bad_pt.failure()

    # Missing design column
    res_bad_design = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col="patient_id",
        design_col="non_existent_design",
    )
    assert isinstance(res_bad_design, Failure)
    assert "Design column 'non_existent_design' not found" in res_bad_design.failure()
