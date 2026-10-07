"""Unit tests for the 17 Pan-Cancer Atlas dataset loaders, registry, and matrix inspection."""

from pathlib import Path
import gzip
import tempfile
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.io as sio
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Success

from tme_datasets.preprocessing.matrix_inspection import (
    ExpressionInspectionResult,
    ExpressionType,
    inspect_expression_type,
    tag_expression_metadata,
)
from tme_datasets.providers.atlas_single_cell import (
    _apply_subset_and_subsample,
    _standardize_obs,
    load_azizi_brca,
    load_becker_coad,
    load_biermann_brainmet,
    load_borcherding_ccrcc,
    load_cheng_pancancer,
    load_durante_uvm,
    load_khaliq_cc,
    load_kim_luad,
    load_leader_nsclc,
    load_lu_hcc,
    load_pelka_crc,
    load_pu_ptc,
    load_qian_pancancer,
    load_sharma_hcc,
    load_vazquez_ov,
    load_zhang_myeloid,
    load_zhang_tnbc,
)
from tme_datasets.query import load_dataset
from tme_datasets.registry import (
    DATASET_ALIASES,
    DATASET_REGISTRY,
    get_dataset_spec,
    list_registered_datasets,
    resolve_dataset_id,
)
from tme_datasets.transforms.perturbations import randomize_negative_binomial
from tme_datasets.models import NegativeBinomialConfig


# =========================================================================
# 1. Tests for Matrix Inspection and Expression Classification
# =========================================================================
def test_inspect_expression_type_raw_counts() -> None:
    """Verifies that integer UMI count matrices are accurately classified as RAW_COUNTS."""
    # Synthetic raw counts with integer values and varying cell depths
    rng = np.random.default_rng(42)
    raw_mat = sp.csr_matrix(rng.poisson(lam=5.0, size=(50, 100)).astype(np.float32))
    adata = ad.AnnData(X=raw_mat)

    res = inspect_expression_type(adata)
    assert res.expression_type == ExpressionType.RAW_COUNTS
    assert res.is_raw_counts is True
    assert res.is_integer is True
    assert res.cv_library_depth > 0.01

    tagged = tag_expression_metadata(adata)
    assert tagged.uns["is_raw_counts"] is True
    assert tagged.uns["expression_type"] == "raw_counts"
    assert "counts" in tagged.layers


def test_inspect_expression_type_log_normalized() -> None:
    """Verifies that log-normalized matrices (values <= 30, fractional) are classified correctly."""
    rng = np.random.default_rng(42)
    # log1p of normalized counts produces small non-integer values <= 10.0
    log_mat = sp.csr_matrix(np.log1p(rng.exponential(scale=2.0, size=(50, 100)).astype(np.float32)))
    adata = ad.AnnData(X=log_mat)

    res = inspect_expression_type(adata)
    assert res.expression_type == ExpressionType.LOG_NORMALIZED
    assert res.is_raw_counts is False
    assert res.is_integer is False
    assert res.max_value <= 30.0

    tagged = tag_expression_metadata(adata)
    assert tagged.uns["is_raw_counts"] is False
    assert tagged.uns["expression_type"] == "log_normalized"
    assert "log1p" in tagged.layers
    assert "linear" in tagged.layers  # verifies expm1 linear layer was generated


def test_negative_binomial_safeguard_on_normalized_data() -> None:
    """Verifies that randomize_negative_binomial safely fails on pre-normalized matrices."""
    rng = np.random.default_rng(42)
    log_mat = sp.csr_matrix(np.log1p(rng.exponential(scale=2.0, size=(20, 20)).astype(np.float32)))
    adata = ad.AnnData(X=log_mat)
    adata.uns["is_raw_counts"] = False

    cfg = NegativeBinomialConfig(dispersion=0.1)
    res = randomize_negative_binomial(adata, cfg)
    assert isinstance(res, Failure)
    assert "Cannot apply Negative Binomial count perturbation on pre-normalized expression" in res.failure()


# =========================================================================
# 2. Tests for Registry and Alias Resolution
# =========================================================================
def test_all_17_atlas_datasets_registered() -> None:
    """Verifies all 17 Pan-Cancer Atlas datasets are present in DATASET_REGISTRY."""
    expected_accessions = [
        "GSE178341",  # Pelka
        "GSE114727",  # Azizi
        "E-MTAB-8107", # Qian
        "GSE154763",  # Cheng
        "GSE154826",  # Leader
        "GSE131907",  # Kim
        "GSE201349",  # Becker
        "GSE200997",  # Khaliq
        "GSE121638",  # Borcherding
        "GSE156625",  # Sharma
        "GSE149614",  # Lu
        "GSE184362",  # Pu
        "GSE139829",  # Durante
        "GSE200218",  # Biermann
        "GSE180661",  # Vazquez
        "GSE169246",  # Zhang 2021
        "GSE215120",  # Zhang 2022
    ]

    for acc in expected_accessions:
        spec = get_dataset_spec(acc)
        assert isinstance(spec, Some), f"Accession {acc} not found in registry"
        unwrapped = spec.unwrap()
        assert unwrapped.id == acc
        assert unwrapped.cancer_type is not None

    reg_datasets = list_registered_datasets()
    assert len(reg_datasets) >= 35


def test_alias_resolution_for_atlas_cohorts() -> None:
    """Verifies study aliases map seamlessly to primary accessions."""
    test_cases = [
        ("Pelka", "GSE178341"),
        ("Pelka_CRC", "GSE178341"),
        ("Pelka2021", "GSE178341"),
        ("Azizi_BRCA", "GSE114727"),
        ("Qian_PanCancer", "E-MTAB-8107"),
        ("Cheng_PanCancer", "GSE154763"),
        ("Leader_NSCLC", "GSE154826"),
        ("Kim_LUAD", "GSE131907"),
        ("Becker_COAD", "GSE201349"),
        ("Khaliq_CC", "GSE200997"),
        ("Borcherding_ccRCC", "GSE121638"),
        ("Sharma_HCC", "GSE156625"),
        ("Lu_HCC", "GSE149614"),
        ("Pu_PTC", "GSE184362"),
        ("Durante_UVM", "GSE139829"),
        ("Biermann_BrainMet", "GSE200218"),
        ("Vazquez_OV", "GSE180661"),
        ("Zhang_TNBC", "GSE169246"),
        ("Zhang_Myeloid", "GSE215120"),
    ]

    for alias, expected_id in test_cases:
        resolved = resolve_dataset_id(alias)
        assert resolved == expected_id, f"Alias {alias} resolved to {resolved}, expected {expected_id}"
        spec = get_dataset_spec(alias)
        assert isinstance(spec, Some)
        assert spec.unwrap().id == expected_id


# =========================================================================
# 3. Tests for Subsetting, Subsampling, and Obs Standardization
# =========================================================================
def test_standardize_obs_preserves_columns() -> None:
    """Verifies that _standardize_obs adds canonical columns without overwriting custom ones."""
    obs_df = pd.DataFrame({
        "donor_id": ["P1", "P2", "P3"],
        "batch_id": ["B1", "B2", "B3"],
        "annotation": ["T_cell", "B_cell", "Myeloid"],
        "custom_metric": [0.1, 0.5, 0.9],
    }, index=["cell_A", "cell_B", "cell_C"])

    adata = ad.AnnData(X=sp.csr_matrix((3, 5)), obs=obs_df)
    std_adata = _standardize_obs(
        adata,
        dataset="TestDataset",
        organ="Lung",
        cancer_type="NSCLC",
        cancer_code="LUAD",
        patient_col="donor_id",
        sample_col="batch_id",
        cell_type_col="annotation",
    )

    # Check canonical columns
    assert std_adata.obs["cell_id"].tolist() == ["cell_A", "cell_B", "cell_C"]
    assert std_adata.obs["dataset"].unique().tolist() == ["TestDataset"]
    assert std_adata.obs["organ"].unique().tolist() == ["Lung"]
    assert std_adata.obs["cancer_code"].unique().tolist() == ["LUAD"]
    assert std_adata.obs["patient"].tolist() == ["P1", "P2", "P3"]
    assert std_adata.obs["cell_type_author"].tolist() == ["T_cell", "B_cell", "Myeloid"]

    # Verify custom metric preserved
    assert std_adata.obs["custom_metric"].tolist() == [0.1, 0.5, 0.9]


def test_apply_subset_and_subsample() -> None:
    """Verifies subset filtering and random downsampling."""
    obs_df = pd.DataFrame({
        "cancer_type": ["CRC"] * 60 + ["NSCLC"] * 40,
        "patient": ["P1"] * 50 + ["P2"] * 50,
    }, index=[f"cell_{i}" for i in range(100)])

    adata = ad.AnnData(X=sp.csr_matrix((100, 10)), obs=obs_df)

    # 1. Subset by cancer_type
    sub_adata = _apply_subset_and_subsample(adata, subset={"cancer_type": "CRC"})
    assert sub_adata.n_obs == 60
    assert np.all(sub_adata.obs["cancer_type"] == "CRC")

    # 2. Downsample
    sampled = _apply_subset_and_subsample(adata, subsample_n=25)
    assert sampled.n_obs == 25


# =========================================================================
# 4. Tests for Provider Parsers with Synthetic Mock Files
# =========================================================================
def test_load_sharma_hcc_parser(tmp_path: Path) -> None:
    """Tests load_sharma_hcc on mock matrix and feature files."""
    raw_dir = tmp_path / "sharma_mock"
    raw_dir.mkdir()

    # Generate mock MTX, barcodes, and genes
    rng = np.random.default_rng(42)
    mat = rng.poisson(lam=2.0, size=(10, 5))  # 10 genes x 5 cells
    mtx_p = raw_dir / "GSE156625_HCCmatrix.mtx.gz"
    with gzip.open(mtx_p, "wb") as f:
        sio.mmwrite(f, mat)

    bc_p = raw_dir / "GSE156625_HCCbarcodes.tsv.gz"
    with gzip.open(bc_p, "wt") as f:
        for i in range(5):
            f.write(f"BC_{i}\n")

    gene_p = raw_dir / "GSE156625_HCCgenes.tsv.gz"
    with gzip.open(gene_p, "wt") as f:
        for i in range(10):
            f.write(f"ENSG0000000{i:04d}\tGENE_{i}\n")

    res = load_sharma_hcc(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 5
    assert adata.n_vars == 10
    assert adata.obs["dataset"].iloc[0] == "Sharma2020"
    assert adata.obs["organ"].iloc[0] == "Liver"
    assert adata.obs["cancer_code"].iloc[0] == "LIHC"
    assert adata.var_names[0] == "ENSG00000000000"
    assert adata.var["gene_name"].iloc[0] == "GENE_0"
    assert adata.uns["is_raw_counts"] is True


def test_load_khaliq_cc_parser(tmp_path: Path) -> None:
    """Tests load_khaliq_cc on mock count CSV and annotation CSV."""
    raw_dir = tmp_path / "khaliq_mock"
    raw_dir.mkdir()

    # Mock count CSV: genes as column 0, cell barcodes as remaining columns
    df_counts = pd.DataFrame({
        "gene": [f"ENSG0000000{i:04d}" for i in range(8)],
        "cell_1": [0, 5, 2, 0, 10, 1, 0, 3],
        "cell_2": [1, 0, 3, 4, 0, 2, 0, 0],
        "cell_3": [0, 2, 0, 1, 5, 0, 8, 2],
    })
    counts_p = raw_dir / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
    df_counts.to_csv(counts_p, index=False, compression="gzip")

    df_meta = pd.DataFrame({
        "patient": ["P1", "P1", "P2"],
        "CellType": ["T", "B", "Myeloid"],
    }, index=["cell_1", "cell_2", "cell_3"])
    meta_p = raw_dir / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
    df_meta.to_csv(meta_p, compression="gzip")

    res = load_khaliq_cc(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 3
    assert adata.n_vars == 8
    assert adata.obs["dataset"].iloc[0] == "Khaliq2022"
    assert adata.obs["organ"].iloc[0] == "Colon"
    assert adata.obs["cell_type_author"].tolist() == ["T", "B", "Myeloid"]
    assert adata.uns["is_raw_counts"] is True
