"""Unit tests for Tier 1 single-cell and spatial transcriptomics data loaders."""

from __future__ import annotations

import gzip
from pathlib import Path
import tarfile

import anndata as ad
import numpy as np
import pandas as pd
from returns.maybe import Some
from returns.result import Failure, Success
import scipy.io as sio
import scipy.sparse as sp

from tme_datasets.providers.tier1_single_cell import (
    TIER1_COHORTS,
    _harmonize_single_response,
    is_tier1_dataset,
    load_tier1_cohort,
)
from tme_datasets.registry import get_dataset_spec, resolve_dataset_id


# =========================================================================
# 1. Configuration & Registry Integrity
# =========================================================================
def test_tier1_cohorts_catalog_completeness() -> None:
    """Verify that all 84 Tier 1 cohorts are configured and registered."""
    assert len(TIER1_COHORTS) == 84
    for acc, cfg in TIER1_COHORTS.items():
        assert cfg.accession == acc
        assert is_tier1_dataset(acc)
        # Spec registered
        spec_opt = get_dataset_spec(acc)
        assert isinstance(spec_opt, Some)
        spec = spec_opt.unwrap()
        assert spec.id == acc
        assert spec.cancer_type == cfg.indication


def test_tier1_aliases_resolution() -> None:
    """Verify alias resolution for representative Tier 1 cohorts."""
    assert resolve_dataset_id("GSE274141-Breast") == "GSE274141"
    assert resolve_dataset_id("GSE205506-CRC") == "GSE205506"
    assert resolve_dataset_id("CELLxGENE_2554a654-CRC") == "CELLxGENE_2554a654"


# =========================================================================
# 2. Clinical Response Harmonization
# =========================================================================
def test_response_harmonization() -> None:
    """Verify pure mapping of clinical response terminology."""
    assert _harmonize_single_response("CR")[0] == "responder"
    assert _harmonize_single_response("PR")[0] == "responder"
    assert _harmonize_single_response("pCR")[0] == "responder"
    assert _harmonize_single_response("Favorable")[0] == "responder"
    assert _harmonize_single_response("PD")[0] == "non-responder"
    assert _harmonize_single_response("Residual Disease")[0] == "non-responder"
    assert _harmonize_single_response("Unfavourable")[0] == "non-responder"
    assert _harmonize_single_response("SD")[0] == "stable"
    assert _harmonize_single_response(None)[0] == "not-evaluable"
    assert _harmonize_single_response("NA")[0] == "not-evaluable"


# =========================================================================
# 3. Archetype A: GEO 10x RAW Tarball Loader
# =========================================================================
def test_ingest_geo_10x_tar(tmp_path: Path) -> None:
    """Test pure ingestion of a synthetic 10x MTX tarball."""
    accession = "GSE274139"
    raw_dir = tmp_path / accession
    raw_dir.mkdir(parents=True)

    # Build mock sample files
    sample_dir = tmp_path / "mock_10x"
    sample_dir.mkdir()

    # 10 genes x 5 cells count matrix (cells x genes = 5 x 10)
    counts = sp.csr_matrix(np.random.randint(10, 50, size=(5, 10), dtype=np.int32))
    # sio.mmwrite writes rows x cols, standard 10x has rows = genes, cols = cells
    sio.mmwrite(sample_dir / "matrix.mtx", counts.T)
    with open(sample_dir / "matrix.mtx", "rb") as f_in, gzip.open(sample_dir / "matrix.mtx.gz", "wb") as f_out:
        f_out.writelines(f_in)

    with gzip.open(sample_dir / "barcodes.tsv.gz", "wt") as f:
        for i in range(5):
            f.write(f"AAACCC_{i}-1\n")

    with gzip.open(sample_dir / "features.tsv.gz", "wt") as f:
        for i in range(10):
            f.write(f"ENSG0000000{i}\tGene_{i}\tGene Expression\n")

    # Package into GSE274139_RAW.tar
    tar_path = raw_dir / f"{accession}_RAW.tar"
    with tarfile.open(tar_path, "w") as archive:
        archive.add(sample_dir / "matrix.mtx.gz", arcname="sample1_matrix.mtx.gz")
        archive.add(sample_dir / "barcodes.tsv.gz", arcname="sample1_barcodes.tsv.gz")
        archive.add(sample_dir / "features.tsv.gz", arcname="sample1_features.tsv.gz")

    res = load_tier1_cohort(accession, raw_dir)
    assert isinstance(res, Success), f"Expected Success, got: {res}"
    adata = res.unwrap()

    assert adata.n_obs == 5
    assert adata.n_vars == 10
    assert sp.issparse(adata.X)
    assert "counts" in adata.layers
    assert "dataset_id" in adata.obs.columns
    assert adata.obs["dataset_id"].iloc[0] == accession
    assert adata.obs["indication"].iloc[0] == "Breast"
    assert "cell_id" in adata.obs.columns
    assert "patient_id" in adata.obs.columns
    assert "treatment_status" in adata.obs.columns
    assert "clinical_response" in adata.obs.columns


# =========================================================================
# 4. Archetype B: CELLxGENE H5AD Loader
# =========================================================================
def test_ingest_cellxgene_h5ad(tmp_path: Path) -> None:
    """Test pure ingestion of a synthetic CELLxGENE H5AD restoring raw counts."""
    accession = "CELLxGENE_2554a654"
    raw_dir = tmp_path / accession
    raw_dir.mkdir(parents=True)

    n_cells, n_genes = 8, 15
    raw_mat = sp.csr_matrix(np.random.randint(1, 100, size=(n_cells, n_genes), dtype=np.int32))
    norm_mat = sp.csr_matrix(np.random.uniform(0.1, 5.0, size=(n_cells, n_genes)).astype(np.float32))

    # Raw counts in adata.raw, normalized in adata.X
    mock_adata = ad.AnnData(
        X=norm_mat,
        obs=pd.DataFrame(
            {
                "donor_id": [f"Donor_{i % 2}" for i in range(n_cells)],
                "cell_type": ["T cell"] * n_cells,
                "treatment_response": ["Responder"] * (n_cells // 2) + ["Non-responder"] * (n_cells // 2),
            },
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=[f"GENE_{j}" for j in range(n_genes)]),
    )
    raw_ad = ad.AnnData(X=raw_mat, var=mock_adata.var)
    mock_adata.raw = raw_ad

    h5ad_path = raw_dir / f"{accession}.h5ad"
    mock_adata.write_h5ad(h5ad_path)

    res = load_tier1_cohort(accession, raw_dir)
    assert isinstance(res, Success), f"Expected Success, got: {res}"
    adata = res.unwrap()

    assert adata.n_obs == n_cells
    assert adata.n_vars == n_genes
    # Raw counts restored to adata.X
    assert np.allclose(adata.X.data, np.round(adata.X.data))
    assert adata.uns["is_raw_counts"] is True
    assert "counts" in adata.layers
    assert adata.obs["dataset_id"].iloc[0] == accession
    assert adata.obs["indication"].iloc[0] == "CRC"
    assert "clinical_response" in adata.obs.columns


# =========================================================================
# 5. Archetype C: Flat Count Table (CSV / TSV)
# =========================================================================
def test_ingest_flat_count_table(tmp_path: Path) -> None:
    """Test pure ingestion of a gzipped CSV count matrix."""
    accession = "GSE274141"
    raw_dir = tmp_path / accession
    raw_dir.mkdir(parents=True)

    n_genes, n_cells = 20, 6
    # rows = genes, cols = cells
    df = pd.DataFrame(
        np.random.randint(0, 80, size=(n_genes, n_cells), dtype=np.int32),
        index=[f"ENSG0000_{i}" for i in range(n_genes)],
        columns=[f"Cell_Barcode_{j}" for j in range(n_cells)],
    )

    csv_path = raw_dir / "GSE274141_read-counts-n54-new.csv.gz"
    df.to_csv(csv_path, compression="gzip")

    res = load_tier1_cohort(accession, raw_dir)
    assert isinstance(res, Success), f"Expected Success, got: {res}"
    adata = res.unwrap()

    assert adata.n_obs == n_cells
    assert adata.n_vars == n_genes
    assert sp.issparse(adata.X)
    assert "counts" in adata.layers
    assert adata.obs["dataset_id"].iloc[0] == accession
    assert adata.obs["indication"].iloc[0] == "Breast"


# =========================================================================
# 6. Archetype D: Spatial Assay (Coordinates retention)
# =========================================================================
def test_ingest_spatial_assay(tmp_path: Path) -> None:
    """Test pure ingestion of a spatial assay retaining coordinates in obsm['spatial']."""
    accession = "GSE301720"
    raw_dir = tmp_path / accession
    raw_dir.mkdir(parents=True)

    # Flat table
    n_genes, n_spots = 12, 4
    df = pd.DataFrame(
        np.random.randint(5, 50, size=(n_genes, n_spots), dtype=np.int32),
        index=[f"Gene_{i}" for i in range(n_genes)],
        columns=[f"spot_{j}" for j in range(n_spots)],
    )
    df.to_csv(raw_dir / "counts.csv.gz", compression="gzip")

    # Spatial coordinates
    pos_df = pd.DataFrame(
        {
            0: [f"spot_{j}" for j in range(n_spots)],
            1: [1] * n_spots,
            2: [0, 0, 1, 1],
            3: [0, 1, 0, 1],
            4: [100.0, 150.0, 200.0, 250.0],
            5: [10.0, 20.0, 30.0, 40.0],
        }
    )
    pos_df.to_csv(raw_dir / "tissue_positions_list.csv", header=False, index=False)

    res = load_tier1_cohort(accession, raw_dir)
    assert isinstance(res, Success), f"Expected Success, got: {res}"
    adata = res.unwrap()

    assert adata.n_obs == n_spots
    assert adata.n_vars == n_genes
    assert "spatial" in adata.obsm
    assert adata.obsm["spatial"].shape == (n_spots, 2)


# =========================================================================
# 7. Railway Oriented Error Handling
# =========================================================================
def test_error_handling_nonexistent_cohort(tmp_path: Path) -> None:
    """Verify Failure result returned on unregistered cohort."""
    res = load_tier1_cohort("GSE999999", tmp_path)
    assert isinstance(res, Failure)
    assert "not a registered Tier 1 cohort" in res.failure()


def test_error_handling_missing_archive(tmp_path: Path) -> None:
    """Verify Failure result returned on missing raw files."""
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    res = load_tier1_cohort("GSE274139", empty_dir)
    assert isinstance(res, Failure)
    assert "Tar archive for GSE274139 not found" in res.failure()


# =========================================================================
# 8. Query Dispatcher Integration
# =========================================================================
def test_query_dispatcher_tier1(tmp_path: Path) -> None:
    """Verify that query.load_dataset dispatches correctly to Tier 1 loader."""
    from tme_datasets.query import load_dataset

    accession = "GSE274141"
    raw_dir = tmp_path / accession
    raw_dir.mkdir(parents=True)

    n_genes, n_cells = 10, 4
    df = pd.DataFrame(
        np.random.randint(1, 50, size=(n_genes, n_cells), dtype=np.int32),
        index=[f"ENSG000000000{i}" for i in range(n_genes)],
        columns=[f"Cell_{j}" for j in range(n_cells)],
    )
    df.to_csv(raw_dir / "GSE274141_read-counts-n54-new.csv.gz", compression="gzip")

    res = load_dataset(
        accession,
        base_dir=tmp_path,
        raw_dir=raw_dir,
        auto_download=False,
        normalize_ensembl=False,
        cache_h5ad=False,
        apply_qc=False,
    )
    assert isinstance(res, Success), f"Expected Success from query.load_dataset, got: {res}"
    adata = res.unwrap()
    assert adata.n_obs == n_cells
    assert adata.obs["dataset_id"].iloc[0] == accession
