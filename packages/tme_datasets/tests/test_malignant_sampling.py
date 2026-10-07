"""Unit tests for malignant vs microenvironment sampling."""

from __future__ import annotations

import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import pandas as pd
import pytest
from returns.result import Failure, Success

from tme_datasets.sampling.malignant_sampling import (
    MalignantSamplingConfig,
    MalignantStrategy,
    sample_malignant_and_tme_cells,
)


@pytest.fixture
def synthetic_tme_adata() -> ad.AnnData:
    """Fixture creating AnnData with known malignant and TME cell populations."""
    n_cells = 300
    n_genes = 50
    rng = np.random.default_rng(42)

    X = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
    # 100 malignant cells, 200 TME cells
    is_mal = np.array([True] * 100 + [False] * 200)

    # 4 distinct patients for malignant cells
    patient_ids = [f"Pat_{i % 4}" for i in range(100)] + ["TME_donor"] * 200

    obs = pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(n_cells)],
        "is_malignant": is_mal,
        "patient_id": patient_ids,
        "cell_type": ["Tumor"] * 100 + ["CD8_T"] * 100 + ["Macrophage"] * 100,
    })
    var = pd.DataFrame(index=[f"GENE_{i}" for i in range(n_genes)])

    return ad.AnnData(X=X, obs=obs, var=var)


def test_excluded_tme_only(synthetic_tme_adata: ad.AnnData) -> None:
    """Verify that EXCLUDED_TME_ONLY samples exclusively non-malignant cells."""
    cfg = MalignantSamplingConfig(strategy=MalignantStrategy.EXCLUDED_TME_ONLY)
    res = sample_malignant_and_tme_cells(synthetic_tme_adata, config=cfg, target_total_cells=80)

    assert isinstance(res, Success)
    sampled = res.unwrap()
    assert sampled.n_obs == 80
    assert np.all(~sampled.obs["is_malignant"].to_numpy())
    assert all(ct in ("CD8_T", "Macrophage") for ct in sampled.obs["cell_type"])


def test_pooled_generic_fraction(synthetic_tme_adata: ad.AnnData) -> None:
    """Verify that POOLED_GENERIC respects requested malignant cell fraction."""
    cfg = MalignantSamplingConfig(
        strategy=MalignantStrategy.POOLED_GENERIC,
        malignant_fraction=0.25,
    )
    res = sample_malignant_and_tme_cells(synthetic_tme_adata, config=cfg, target_total_cells=100)

    assert isinstance(res, Success)
    sampled = res.unwrap()
    assert sampled.n_obs == 100

    n_mal = int(np.sum(sampled.obs["is_malignant"].to_numpy()))
    assert n_mal == 25
    assert (100 - n_mal) == 75


def test_patient_stratified_cap(synthetic_tme_adata: ad.AnnData) -> None:
    """Verify that PATIENT_STRATIFIED distributes tumor cells evenly across available patients."""
    cfg = MalignantSamplingConfig(
        strategy=MalignantStrategy.PATIENT_STRATIFIED,
        malignant_fraction=0.40,  # 40 tumor cells out of 100
        max_cells_per_patient=15,
    )
    res = sample_malignant_and_tme_cells(synthetic_tme_adata, config=cfg, target_total_cells=100)

    assert isinstance(res, Success)
    sampled = res.unwrap()
    assert sampled.n_obs == 100

    mal_sub = sampled[sampled.obs["is_malignant"].to_numpy()]
    counts_by_pat = mal_sub.obs["patient_id"].value_counts().to_dict()

    # Should have drawn from all 4 patients
    assert len(counts_by_pat) == 4
    for pat, count in counts_by_pat.items():
        assert count <= 15
