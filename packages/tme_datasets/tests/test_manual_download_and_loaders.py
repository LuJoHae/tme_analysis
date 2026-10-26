"""Unit tests for manual download resolution, informative error handling, and new loaders."""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets import (
    find_dataset_h5ad,
    get_manual_download_dir,
    load_dataset,
)
from tme_datasets.providers.bulk_papers import (
    load_genentech_egad,
    load_paper_h5ad,
)
from tme_datasets.providers.single_cell import (
    load_gse179994,
    load_ma_liver,
    load_maynard,
    load_yost,
)


def test_get_manual_download_dir(tmp_path: Path) -> None:
    """Verify get_manual_download_dir resolves to data/manual_download."""
    manual_dir = get_manual_download_dir(repo_root=tmp_path)
    assert manual_dir == tmp_path / "data/manual_download"


def test_missing_paper_h5ad_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing paper H5AD returns structured Failure with actionable instructions."""
    missing_file = tmp_path / "data/manual_download/Auslander.h5ad"
    res = load_paper_h5ad(missing_file)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "Paper H5AD file not found" in err
    assert "data/manual_download/Auslander.h5ad" in err
    assert "setup_manual_downloads.py" in err


def test_missing_egad_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing EGAD alignment directory returns structured Failure."""
    missing_dir = tmp_path / "data/manual_download/EGAD00001006631-align"
    res = load_genentech_egad(missing_dir)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "EGAD alignment directory not found" in err
    assert "EGAD00001006631-align" in err


def test_missing_maynard_returns_informative_failure(tmp_path: Path) -> None:
    """Verify missing Maynard dataset returns structured Failure."""
    from returns.maybe import Nothing
    with patch("tme_datasets.paths.find_dataset_h5ad", return_value=Nothing):
        res = load_maynard(tmp_path, auto_download=False)
    assert isinstance(res, Failure)
    err = res.failure()
    assert "Maynard NSCLC dataset not found" in err
    assert "Maynard_NSCLC.h5ad" in err


def test_download_and_build_maynard_full_with_mock_data(tmp_path: Path) -> None:
    """Verify download_and_build_maynard_full parses expression CSV and metadata."""
    from tme_datasets.providers.single_cell import download_and_build_maynard_full

    raw_dir = tmp_path / "Maynard_NSCLC"
    raw_dir.mkdir(parents=True)

    # 1. Mock metacells CSV
    meta_csv = raw_dir / "S01_metacells.csv"
    meta_csv.write_text(
        ",cell_id,patient_id,sample_name,biopsy_time_status\n"
        "0,c1,P01,S01,PR\n"
        "1,c2,P02,S02,PD\n"
    )

    # 2. Mock expression CSV (2 genes x 2 cells)
    data_csv = raw_dir / "S01_datafinal.csv"
    data_csv.write_text(
        '"","c1","c2"\n'
        '"CD8A",5.0,0.0\n'
        '"PDCD1",2.5,8.0\n'
    )

    out_h5ad = tmp_path / "Maynard_NSCLC.h5ad"
    res = download_and_build_maynard_full(raw_dir=raw_dir, output_h5ad=out_h5ad, force_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 2
    assert adata.n_vars == 2
    assert list(adata.obs_names) == ["c1", "c2"]
    assert list(adata.var_names) == ["ENSG00000153563", "ENSG00000188389"]
    assert list(adata.var["gene_name"]) == ["CD8A", "PDCD1"]
    assert adata.obs.loc["c1", "response_binary"] == 1.0
    assert adata.obs.loc["c2", "response_binary"] == 0.0
    assert out_h5ad.exists()


def test_load_ma_liver_with_mock_files(tmp_path: Path) -> None:
    """Verify load_ma_liver correctly reads 10x matrix, genes, and barcodes."""
    raw_dir = tmp_path / "GSE125449"
    raw_dir.mkdir(parents=True)

    # 1. Create dummy MatrixMarket file (3 genes x 2 cells)
    mm_content = (
        "%%MatrixMarket matrix coordinate real general\n"
        "%test\n"
        "3 2 4\n"
        "1 1 5.0\n"
        "2 1 10.0\n"
        "1 2 2.0\n"
        "3 2 8.0\n"
    )
    with gzip.open(raw_dir / "GSE125449_Set1_matrix.mtx.gz", "wt") as f:
        f.write(mm_content)

    with gzip.open(raw_dir / "GSE125449_Set1_genes.tsv.gz", "wt") as f:
        f.write("CD8A\nCD4\nFOXP3\n")

    with gzip.open(raw_dir / "GSE125449_Set1_barcodes.tsv.gz", "wt") as f:
        f.write("cell1\ncell2\n")

    res = load_ma_liver(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert list(adata.obs_names) == ["cell1", "cell2"]
    assert list(adata.var_names) == ["CD8A", "CD4", "FOXP3"]


def test_load_yost_with_mock_counts(tmp_path: Path) -> None:
    """Verify load_yost parses tab-delimited counts and annotates patient and response."""
    raw_dir = tmp_path / "GSE123813"
    raw_dir.mkdir(parents=True)

    # Header: Gene, bcc.su001.pre, bcc.su005.post
    counts_content = (
        "Gene\tbcc.su001.pre\tbcc.su005.post\n"
        "CD8A\t10\t0\n"
        "PDCD1\t5\t2\n"
        "CTLA4\t8\t1\n"
    )
    with gzip.open(raw_dir / "GSE123813_bcc_scRNA_counts.txt.gz", "wt") as f:
        f.write(counts_content)

    res = load_yost(raw_dir, auto_download=False)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert "patient" in adata.obs.columns
    assert "response_binary" in adata.obs.columns
    # su001 is Responder (1), su005 is Non-responder (0)
    assert adata.obs.loc["bcc.su001.pre", "response_binary"] == 1
    assert adata.obs.loc["bcc.su005.post", "response_binary"] == 0


def test_load_gse179994_with_mock_rds(tmp_path: Path) -> None:
    """Verify load_gse179994 reads RDS dataframe and attaches condition metadata."""
    raw_dir = tmp_path / "GSE179994"
    raw_dir.mkdir(parents=True)

    rds_file = raw_dir / "GSE179994_all.Tcell.rawCounts.rds"
    # Create mock RDS object via pyreadr mock
    mock_df = pd.DataFrame(
        [[5, 0], [2, 8]],
        index=["CD8A", "CD4"],
        columns=["cell1", "cell2"],
    )

    with patch("pyreadr.read_r", return_value={"counts": mock_df}):
        rds_file.touch()
        res = load_gse179994(raw_dir, auto_download=False)
        assert isinstance(res, Success)
        adata = res.unwrap()
        assert adata.n_obs == 2
        assert adata.n_vars == 2


def test_find_dataset_h5ad_finds_manual_download(tmp_path: Path) -> None:
    """Verify find_dataset_h5ad strictly looks in preprocessed_dir and does not hunt through manual_download."""
    manual_dir = tmp_path / "data/manual_download"
    manual_dir.mkdir(parents=True)
    fake_h5ad = manual_dir / "Auslander.h5ad"
    fake_h5ad.write_bytes(b"mock_h5ad_data")

    # Under strict single-path lookup, searching in manual_download is disabled
    found = find_dataset_h5ad("Auslander", repo_root=tmp_path)
    assert found is Nothing

    # When placed in the exact preprocessed path, it is discovered
    preprocessed_dir = tmp_path / "data/preprocessed"
    preprocessed_dir.mkdir(parents=True)
    target_h5ad = preprocessed_dir / "Auslander.h5ad"
    target_h5ad.write_bytes(b"mock_h5ad_data")
    found_prep = find_dataset_h5ad("Auslander", repo_root=tmp_path)
    assert isinstance(found_prep, Some)
    assert found_prep.unwrap() == target_h5ad


def test_bulk_lung_registry_and_aliases() -> None:
    """Verify Ravi, Prat, Jung (GSE135222), and Cho (GSE126044) are registered and resolved."""
    from tme_datasets import get_dataset_spec, resolve_dataset_id

    assert resolve_dataset_id("GSE135222") == "Jung"
    assert resolve_dataset_id("Jung_NSCLC") == "Jung"
    assert resolve_dataset_id("GSE126044") == "Cho"
    assert resolve_dataset_id("Cho_NSCLC") == "Cho"
    assert resolve_dataset_id("GSE93157") == "Prat"
    assert resolve_dataset_id("Ravi_NSCLC") == "Ravi"

    for cohort in ("Ravi", "Prat", "Jung", "Cho"):
        spec_maybe = get_dataset_spec(cohort)
        assert isinstance(spec_maybe, Some)
        spec = spec_maybe.unwrap()
        assert spec.id == cohort
        assert spec.has_response_labels is True


def test_preprocess_ravi_with_mock_excel(tmp_path: Path) -> None:
    """Verify preprocess_ravi parses Table_S1_Clinical_Annotations and Table_S13_RNA_TPM."""
    from tme_datasets.providers.bulk_papers import preprocess_ravi

    raw_dir = tmp_path / "Ravi"
    raw_dir.mkdir(parents=True)
    xlsx_path = raw_dir / "41588_2023_1355_MOESM3_ESM.xlsx"

    clin_df = pd.DataFrame(
        {
            "Harmonized_SU2C_Participant_ID_v2": ["P1", "P2"],
            "Harmonized_SU2C_RNA_Tumor_Sample_ID_v2": ["P1-T1", "P2-T1"],
            "Patient_Age_at_Diagnosis": [65.0, 58.0],
            "Patient_Sex": ["M", "F"],
            "Histology_Harmonized": ["Adenocarcinoma", "Squamous"],
            "Initial_Stage": [4.0, 3.0],
            "Agent_PD1": ["Pembrolizumab", "Nivolumab"],
            "Agent_PD1_Category": ["PD-1", "PD-1"],
            "Line_of_Therapy": [1.0, 2.0],
            "Harmonized_Confirmed_BOR": ["PR", "PD"],
            "Harmonized_PFS_Days": [300.0, 45.0],
            "Harmonized_PFS_Event": [0.0, 1.0],
            "Harmonized_OS_Days": [600.0, 120.0],
            "Harmonized_OS_Event": [0.0, 1.0],
        }
    )
    tpm_df = pd.DataFrame(
        {
            "Name": ["ENSG00000153563.15", "ENSG00000188389.10"],
            "Description": ["CD8A", "PDCD1"],
            "P1-T1": [42.5, 12.0],
            "P2-T1": [3.1, 0.5],
        }
    )

    with pd.ExcelWriter(xlsx_path) as writer:
        # Both sheets have 2 header title/blank rows before column names at row index 2
        pd.DataFrame([["Title"], [""]]).to_excel(
            writer, sheet_name="Table_S1_Clinical_Annotations", index=False, header=False
        )
        clin_df.to_excel(
            writer, sheet_name="Table_S1_Clinical_Annotations", startrow=2, index=False
        )
        pd.DataFrame([["Title"], [""]]).to_excel(
            writer, sheet_name="Table_S13_RNA_TPM", index=False, header=False
        )
        tpm_df.to_excel(
            writer, sheet_name="Table_S13_RNA_TPM", startrow=2, index=False
        )

    res = preprocess_ravi(xlsx_path)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.shape == (2, 2)
    assert list(adata.obs_names) == ["P1-T1", "P2-T1"]
    assert list(adata.var_names) == ["ENSG00000153563", "ENSG00000188389"]
    assert list(adata.var["gene_name"]) == ["CD8A", "PDCD1"]
    assert adata.obs.loc["P1-T1", "response"] == "R"
    assert adata.obs.loc["P2-T1", "response"] == "NR"
    assert adata.obs.loc["P1-T1", "cancer_type"] == "NSCLC"
    assert adata.uns["expression_type"] == "tpm"


def test_preprocess_prat_with_mock_geo_files(tmp_path: Path) -> None:
    """Verify preprocess_prat parses GSE93157 expression and series matrix metadata."""
    from tme_datasets.providers.bulk_papers import preprocess_prat

    raw_dir = tmp_path / "Prat"
    raw_dir.mkdir(parents=True)

    expr_content = (
        "ID_REF\tSAMPLE 1\tSAMPLE 2\n"
        "CD8A\t8.5\t2.1\n"
        "PDCD1\t6.2\t1.0\n"
    )
    with gzip.open(raw_dir / "GSE93157_raw_data_values.txt.gz", "wt") as f:
        f.write(expr_content)

    matrix_content = (
        '!Sample_title\t"Patient1"\t"Patient2"\n'
        '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
        '!Sample_source_name_ch1\t"LUNG NON-SQUAMOUS CANCER"\t"MELANOMA"\n'
        '!Sample_characteristics_ch1\t"age: 64"\t"age: 52"\n'
        '!Sample_characteristics_ch1\t"Sex: M"\t"Sex: F"\n'
        '!Sample_characteristics_ch1\t"drug: NIVOLUMAB"\t"drug: PEMBROLIZUMAB"\n'
        '!Sample_characteristics_ch1\t"best.resp: PR"\t"best.resp: PD"\n'
        '!Sample_characteristics_ch1\t"response: 1"\t"response: 0"\n'
        '!Sample_characteristics_ch1\t"pfs: 14.5"\t"pfs: 2.1"\n'
        '!Sample_characteristics_ch1\t"pfse: 0"\t"pfse: 1"\n'
    )
    with gzip.open(raw_dir / "GSE93157_series_matrix.txt.gz", "wt") as f:
        f.write(matrix_content)

    res = preprocess_prat(raw_dir)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.shape == (2, 2)
    assert list(adata.obs_names) == ["Patient1", "Patient2"]
    assert adata.obs.loc["Patient1", "cancer_type"] == "NSCLC"
    assert adata.obs.loc["Patient2", "cancer_type"] == "Melanoma"
    assert adata.obs.loc["Patient1", "response"] == "R"
    assert adata.obs.loc["Patient2", "response"] == "NR"


def test_preprocess_jung_with_mock_geo_files(tmp_path: Path) -> None:
    """Verify preprocess_jung parses GSE135222 TPM matrix and attaches Table S3 response."""
    from tme_datasets.providers.bulk_papers import preprocess_jung

    raw_dir = tmp_path / "Jung"
    raw_dir.mkdir(parents=True)

    expr_content = (
        "gene_id\tNSCLC378\tNSCLC556\n"
        "ENSG00000153563.15\t55.0\t4.0\n"
        "ENSG00000188389.10\t20.0\t1.5\n"
    )
    with gzip.open(raw_dir / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", "wt") as f:
        f.write(expr_content)

    matrix_content = (
        '!Sample_title\t"Macrogen_NSCLC378_T"\t"Macrogen_NSCLC556_T"\n'
        '!Sample_geo_accession\t"GSM3995420"\t"GSM3995421"\n'
        '!Sample_characteristics_ch1\t"age: 60"\t"age: 65"\n'
        '!Sample_characteristics_ch1\t"Sex: M"\t"Sex: F"\n'
        '!Sample_characteristics_ch1\t"pfs time: 461"\t"pfs time: 42"\n'
        '!Sample_characteristics_ch1\t"progression: 0"\t"progression: 1"\n'
    )
    with gzip.open(raw_dir / "GSE135222_series_matrix.txt.gz", "wt") as f:
        f.write(matrix_content)

    res = preprocess_jung(raw_dir)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.shape == (2, 2)
    assert list(adata.var_names) == ["ENSG00000153563", "ENSG00000188389"]
    assert adata.obs.loc["NSCLC378", "response"] == "R"
    assert adata.obs.loc["NSCLC556", "response"] == "NR"
    assert adata.obs.loc["NSCLC378", "cancer_type"] == "NSCLC"
    assert adata.uns["expression_type"] == "tpm"


def test_preprocess_cho_with_mock_geo_files(tmp_path: Path) -> None:
    """Verify preprocess_cho parses GSE126044 counts and attaches Table S2A clinical metadata."""
    from tme_datasets.providers.bulk_papers import preprocess_cho

    raw_dir = tmp_path / "Cho"
    raw_dir.mkdir(parents=True)

    counts_content = (
        "\tDis_01\tDis_02\n"
        "CD8A\t12\t450\n"
        "PDCD1\t3\t180\n"
    )
    with gzip.open(raw_dir / "GSE126044_counts.txt.gz", "wt") as f:
        f.write(counts_content)

    matrix_content = (
        '!Sample_title\t"Dis_01"\t"Dis_02"\n'
        '!Sample_geo_accession\t"GSM3589673"\t"GSM3589674"\n'
        '!Sample_characteristics_ch1\t"disease state: non-small cell lung cancer"\t"disease state: non-small cell lung cancer"\n'
    )
    with gzip.open(raw_dir / "GSE126044_series_matrix.txt.gz", "wt") as f:
        f.write(matrix_content)

    res = preprocess_cho(raw_dir)
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.shape == (2, 2)
    assert list(adata.obs_names) == ["Dis_01", "Dis_02"]
    assert adata.obs.loc["Dis_01", "response"] == "NR"
    assert adata.obs.loc["Dis_02", "response"] == "R"
    assert adata.obs.loc["Dis_02", "RECIST"] == "PR"
    assert adata.obs.loc["Dis_01", "cancer_type"] == "NSCLC"
    assert adata.uns["expression_type"] == "counts"

