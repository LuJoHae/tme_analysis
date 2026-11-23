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
    """Verify loading a single cohort returns a Success[AnnData] with no layers and normalization metadata in uns."""
    result = load_iatlas_cohort_or_combined("Hugo-iAtlas")
    assert isinstance(result, Success)
    adata = result.unwrap()
    assert adata.n_obs == 27
    assert adata.n_vars > 1000
    assert list(adata.layers.keys()) == []
    assert adata.uns.get("is_raw_counts") is False
    assert adata.uns.get("expression_type") == "log_normalized"
    assert adata.uns.get("normalization_method") == "log_TPM"
    assert adata.uns.get("normalization", {}).get("unit") == "TPM"
    assert adata.uns.get("normalization", {}).get("scale") == "log"


def test_load_iatlas_cohort_or_combined_combined():
    """Verify loading combined melanoma group returns Success[AnnData] with cohort column and no layers."""
    result = load_iatlas_cohort_or_combined("melanoma")
    assert isinstance(result, Success)
    adata = result.unwrap()
    assert adata.n_obs >= 338
    assert "dataset_id" in adata.obs.columns
    assert list(adata.layers.keys()) == []


def test_infer_iatlas_normalization_metadata_units_and_scales():
    """Verify _infer_iatlas_normalization_metadata classifies TPM, FPKM, RPKM, and upper-quartile counts."""
    from tme_datasets.preprocessing.matrix_inspection import ExpressionInspectionResult, ExpressionType
    from tme_datasets.providers.bulk_iatlas import _infer_iatlas_normalization_metadata

    log_insp = ExpressionInspectionResult(
        expression_type=ExpressionType.LOG_NORMALIZED,
        is_raw_counts=False,
        is_integer=False,
        min_value=0.0,
        max_value=22.0,
        mean_library_depth=150000.0,
        cv_library_depth=0.05,
        n_zeros=100,
        sparsity=0.2,
    )
    lin_insp = ExpressionInspectionResult(
        expression_type=ExpressionType.LINEAR_NORMALIZED,
        is_raw_counts=False,
        is_integer=False,
        min_value=0.0,
        max_value=45000.0,
        mean_library_depth=1_000_000.0,
        cv_library_depth=0.0001,
        n_zeros=100,
        sparsity=0.2,
    )

    fpkm_meta = _infer_iatlas_normalization_metadata("data_mrna_seq_fpkm.txt", {"profile_name": "mRNA (FPKM)"}, lin_insp)
    assert fpkm_meta["unit"] == "FPKM"
    assert fpkm_meta["scale"] == "linear"
    assert fpkm_meta["method"] == "FPKM"

    rpkm_meta = _infer_iatlas_normalization_metadata("data_mrna_seq_rpkm.txt", {"profile_name": "mRNA (log RPKM)"}, log_insp)
    assert rpkm_meta["unit"] == "RPKM"
    assert rpkm_meta["scale"] == "log"
    assert rpkm_meta["method"] == "log_RPKM"

    uq_meta = _infer_iatlas_normalization_metadata(
        "data_mrna_seq_expression.txt",
        {"profile_name": "mRNA expression (log2 upper quartile normalized counts)"},
        log_insp,
    )
    assert uq_meta["unit"] == "upper_quartile_counts"
    assert uq_meta["scale"] == "log2"
    assert uq_meta["method"] == "log2_upper_quartile_counts"
