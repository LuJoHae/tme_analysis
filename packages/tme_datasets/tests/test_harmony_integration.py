"""Tests for tme_datasets harmony batch integration and milopy cohort resolution."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Success
from tme_datasets import (
    MILOPY_RESPONSE_COHORTS,
    resolve_milopy_cohorts,
    run_harmony_batch_integration,
)


def test_resolve_milopy_cohorts_all() -> None:
    res = resolve_milopy_cohorts("all")
    assert isinstance(res, Success)
    cohorts = res.unwrap()
    assert len(cohorts) == 9
    assert "GSE120575" in cohorts
    assert "GSE316195" in cohorts


def test_resolve_milopy_cohorts_by_cancer_type() -> None:
    res_mel = resolve_milopy_cohorts("cancer_type", "melanoma")
    assert isinstance(res_mel, Success)
    assert set(res_mel.unwrap()) == {"GSE120575", "CELLxGENE_7b20c613"}

    res_lung = resolve_milopy_cohorts("cancer_type", "nsclc")
    assert isinstance(res_lung, Success)
    assert set(res_lung.unwrap()) == {"GSE207422", "GSE243013", "GSE233203"}


def test_resolve_milopy_cohorts_by_treatment() -> None:
    res_mono = resolve_milopy_cohorts("treatment", "anti-pd1_monotherapy")
    assert isinstance(res_mono, Success)
    assert set(res_mono.unwrap()) == {"CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485"}


def test_resolve_milopy_cohorts_by_platform() -> None:
    res_10x = resolve_milopy_cohorts("platform", "10x")
    assert isinstance(res_10x, Success)
    # Matches 10x_chromium_3p, 10x_chromium_5p, 10x_chromium
    assert "CELLxGENE_7b20c613" in res_10x.unwrap()
    assert "GSE243013" in res_10x.unwrap()


def test_resolve_milopy_cohorts_custom() -> None:
    res = resolve_milopy_cohorts("custom", "GSE120575, GSE200996")
    assert isinstance(res, Success)
    assert res.unwrap() == ("GSE120575", "GSE200996")


def test_run_harmony_batch_integration_synthetic() -> None:
    # 60 cells, 2 batches, 20 genes
    rng = np.random.default_rng(42)
    n_cells = 60
    n_genes = 20
    X = rng.poisson(lam=3.0, size=(n_cells, n_genes)).astype(np.float32)
    batches = ["batch_A"] * 30 + ["batch_B"] * 30

    obs = pd.DataFrame({
        "cell_id": [f"c_{i}" for i in range(n_cells)],
        "batch_id": batches,
    })
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    res = run_harmony_batch_integration(adata, batch_key="batch_id", n_comps=5, max_iter=5, seed=42)
    assert isinstance(res, Success)
    harm_adata = res.unwrap()
    assert "X_pca_raw" in harm_adata.obsm
    assert "X_pca_harmony" in harm_adata.obsm
    assert "X_pca" in harm_adata.obsm
    assert harm_adata.obsm["X_pca_harmony"].shape[0] == n_cells
    assert harm_adata.obsm["X_pca_harmony"].shape[1] == 5
