"""Tests for load_iatlas_cohort_or_combined and immunotherapy gene collections."""

import pytest
from returns.result import Success, Failure

from tme_datasets import (
    load_iatlas_cohort_or_combined,
    IMMUNE_CHECKPOINT_GENES,
    IMMUNOTHERAPY_GENE_PANEL,
)


def test_immunotherapy_gene_collections():
    """Verify standard gene panel collections are defined and non-empty."""
    assert len(IMMUNE_CHECKPOINT_GENES) == 13
    assert "CD274" in IMMUNE_CHECKPOINT_GENES
    assert "PDCD1" in IMMUNE_CHECKPOINT_GENES
    assert "CTLA4" in IMMUNE_CHECKPOINT_GENES

    assert len(IMMUNOTHERAPY_GENE_PANEL) == 101
    assert "CD274" in IMMUNOTHERAPY_GENE_PANEL
    assert "IFNG" in IMMUNOTHERAPY_GENE_PANEL
    assert "CXCL9" in IMMUNOTHERAPY_GENE_PANEL
    assert "TBX21" in IMMUNOTHERAPY_GENE_PANEL


def test_load_iatlas_cohort_or_combined_single():
    """Verify loading a single cohort returns a Success[AnnData]."""
    result = load_iatlas_cohort_or_combined("Hugo-iAtlas")
    assert isinstance(result, Success)
    adata = result.unwrap()
    assert adata.n_obs == 27
    assert adata.n_vars > 1000


def test_load_iatlas_cohort_or_combined_combined():
    """Verify loading combined melanoma group returns Success[AnnData] with cohort column."""
    result = load_iatlas_cohort_or_combined("melanoma")
    assert isinstance(result, Success)
    adata = result.unwrap()
    assert adata.n_obs >= 338
    assert "dataset_id" in adata.obs.columns
