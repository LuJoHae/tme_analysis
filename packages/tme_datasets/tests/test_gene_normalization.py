"""Unit tests for Ensembl gene ID normalization, conflict resolution, attribute enrichment, and caching."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
from returns.result import Success

from gene_utils import normalize_genes_to_ensembl
from tme_datasets.paths import get_ensembl_dir
from tme_datasets.preprocessing.gene_normalization import normalize_dataset_to_ensembl


class DummyGene:
    """Mock pyensembl Gene object."""

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        contig: str,
        start: int = 1000,
        end: int = 5000,
        strand: str = "+",
        biotype: str = "protein_coding",
    ) -> None:
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.biotype = biotype


class DummySpecies:
    latin_name = "homo_sapiens"


class MockEnsemblRelease:
    """Mock pyensembl.EnsemblRelease for fast, deterministic unit testing."""

    def __init__(self, release: int = 111, species: str = "human") -> None:
        self.release = release
        self.species = DummySpecies()
        self.download_cache = MagicMock()
        self.db = MagicMock()
        self.db._database_file_exists.return_value = True

        # Pre-configured test database
        self._genes_by_id = {
            "ENSG00000153563": DummyGene("ENSG00000153563", "CD8A", "2", 86782800, 86807185, "+", "protein_coding"),
            "ENSG00000204531": DummyGene("ENSG00000204531", "POU5F1", "6", 31164335, 31170700, "+", "protein_coding"),
            "ENSG00000233911": DummyGene("ENSG00000233911", "POU5F1P1", "CHR_HSCHR6_MHC_APD_CTG1", 1000, 2000, "+", "pseudogene"),
            "ENSG00000012048": DummyGene("ENSG00000012048", "BRCA1", "17", 43044295, 43125483, "-", "protein_coding"),
            "ENSG00000999999": DummyGene("ENSG00000999999", "AMBIGUOUS", "CHR_PATCH_1", 500, 1500, "+", "pseudogene"),
            "ENSG00000111111": DummyGene("ENSG00000111111", "AMBIGUOUS", "1", 10000, 25000, "+", "protein_coding"),
        }

        self._names_to_ids = {
            "CD8A": ["ENSG00000153563"],
            "POU5F1": ["ENSG00000204531"],
            "BRCA1": ["ENSG00000012048"],
            "AMBIGUOUS": ["ENSG00000999999", "ENSG00000111111"],
        }

    def required_local_files_exist(self) -> bool:
        return True

    def gene_by_id(self, gene_id: str) -> DummyGene:
        if gene_id in self._genes_by_id:
            return self._genes_by_id[gene_id]
        raise ValueError(f"Gene ID {gene_id} not found")

    def gene_ids_of_gene_name(self, gene_name: str) -> list[str]:
        if gene_name in self._names_to_ids:
            return self._names_to_ids[gene_name]
        return []


def test_get_ensembl_dir() -> None:
    """Verify get_ensembl_dir returns configured data/ensembl directory."""
    dir_path = get_ensembl_dir()
    assert isinstance(dir_path, Path)
    assert dir_path.name == "ensembl"
    assert dir_path.parent.name == "data"


def test_normalize_genes_exact_and_conflict_resolution(tmp_path: Path) -> None:
    """Verify 1-to-1 matching and 1-to-many canonical contig prioritization."""
    mock_ens = MockEnsemblRelease(111)

    # 3 genes: CD8A (exact 1-to-1), AMBIGUOUS (has 2 IDs, one canonical chr1 and one patch), UNKNOWN_XYZ (unmapped)
    genes = ["CD8A", "AMBIGUOUS", "UNKNOWN_XYZ"]
    counts = np.array([[10.0, 5.0, 1.0], [20.0, 15.0, 2.0]], dtype=np.float32)
    obs = pd.DataFrame(index=["cell_1", "cell_2"])
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=counts, obs=obs, var=var)

    with patch("gene_utils._gene_utils.ensure_ensembl_release_installed", return_value=mock_ens):
        norm = normalize_genes_to_ensembl(
            adata,
            release=111,
            ensembl_dir=tmp_path,
            drop_unmapped=True,
            aggregation="sum",
        )

    # UNKNOWN_XYZ dropped, leaving 2 genes
    assert norm.n_vars == 2
    assert "ENSG00000153563" in norm.var_names  # CD8A
    assert "ENSG00000111111" in norm.var_names  # AMBIGUOUS (prioritized chr1 protein_coding over patch)

    # Verify conflict metadata
    ambig_row = norm.var.loc["ENSG00000111111"]
    assert ambig_row["contig"] == "1"
    assert ambig_row["biotype"] == "protein_coding"
    assert ambig_row["mapping_status"] == "contig_prioritized"
    assert "ENSG00000999999" in ambig_row["alternative_ensembl_ids"]

    # Verify unmapped genes recorded in uns
    assert "unmapped_genes" in norm.uns
    assert "UNKNOWN_XYZ" in norm.uns["unmapped_genes"]
    assert norm.uns["n_unmapped_genes"] == 1


def test_normalize_genes_alias_resolution(tmp_path: Path) -> None:
    """Verify alias resolution maps previous symbol OCT4 -> POU5F1 -> ENSG00000204531."""
    mock_ens = MockEnsemblRelease(111)

    genes = ["OCT4"]
    counts = np.array([[42.0]], dtype=np.float32)
    adata = ad.AnnData(X=counts, obs=pd.DataFrame(index=["cell_1"]), var=pd.DataFrame(index=genes))

    # Mock mygene returning alias mapping
    mock_mg = MagicMock()
    mock_mg.querymany.return_value = [{"query": "OCT4", "symbol": "POU5F1"}]

    with patch("gene_utils._gene_utils.ensure_ensembl_release_installed", return_value=mock_ens), \
         patch("mygene.MyGeneInfo", return_value=mock_mg):
        norm = normalize_genes_to_ensembl(
            adata,
            release=111,
            ensembl_dir=tmp_path,
            drop_unmapped=True,
        )

    assert norm.n_vars == 1
    assert norm.var_names[0] == "ENSG00000204531"
    row = norm.var.loc["ENSG00000204531"]
    assert row["gene_name"] == "POU5F1"
    assert row["original_id"] == "OCT4"
    assert row["mapping_status"] == "alias_resolved"


def test_normalize_genes_duplicate_aggregation(tmp_path: Path) -> None:
    """Verify duplicate Ensembl IDs are aggregated correctly (sum and mean)."""
    mock_ens = MockEnsemblRelease(111)

    # Suppose raw dataset has two probes: "CD8A_probe1" and "CD8A_probe2" that both map to CD8A
    # We configure mock_ens to map both probe names to CD8A's ENSG
    mock_ens._names_to_ids["CD8A_probe1"] = ["ENSG00000153563"]
    mock_ens._names_to_ids["CD8A_probe2"] = ["ENSG00000153563"]

    genes = ["CD8A_probe1", "CD8A_probe2"]
    counts = np.array([[10.0, 20.0], [5.0, 15.0]], dtype=np.float32)
    adata = ad.AnnData(X=counts, obs=pd.DataFrame(index=["c1", "c2"]), var=pd.DataFrame(index=genes))

    with patch("gene_utils._gene_utils.ensure_ensembl_release_installed", return_value=mock_ens):
        # 1. Sum aggregation
        norm_sum = normalize_genes_to_ensembl(
            adata,
            release=111,
            ensembl_dir=tmp_path,
            aggregation="sum",
        )
        assert norm_sum.n_vars == 1
        assert norm_sum.var_names[0] == "ENSG00000153563"
        # 10 + 20 = 30, 5 + 15 = 20
        np.testing.assert_allclose(norm_sum.X.flatten(), [30.0, 20.0])

        # 2. Mean aggregation
        norm_mean = normalize_genes_to_ensembl(
            adata,
            release=111,
            ensembl_dir=tmp_path,
            aggregation="mean",
        )
        # (10 + 20) / 2 = 15, (5 + 15) / 2 = 10
        np.testing.assert_allclose(norm_mean.X.flatten(), [15.0, 10.0])


def test_persistent_parquet_caching(tmp_path: Path) -> None:
    """Verify mapping cache is persisted as Parquet and reused on subsequent calls."""
    mock_ens = MockEnsemblRelease(111)
    cache_file = tmp_path / "gene_mapping_cache_release_111.parquet"
    assert not cache_file.exists()

    genes = ["CD8A", "BRCA1"]
    adata = ad.AnnData(
        X=np.ones((2, 2), dtype=np.float32),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=genes),
    )

    with patch("gene_utils._gene_utils.ensure_ensembl_release_installed", return_value=mock_ens):
        # First call populates cache
        _ = normalize_genes_to_ensembl(adata, release=111, ensembl_dir=tmp_path)
        assert cache_file.exists()

        # Inspect parquet cache via Polars
        df_cached = pl.read_parquet(cache_file)
        assert len(df_cached) == 2
        assert "CD8A" in df_cached["query"].to_list()
        assert "BRCA1" in df_cached["query"].to_list()

        # Second call reads from parquet cache
        norm2 = normalize_genes_to_ensembl(adata, release=111, ensembl_dir=tmp_path)
        assert norm2.n_vars == 2
        assert set(norm2.var_names) == {"ENSG00000153563", "ENSG00000012048"}


def test_tme_datasets_wrapper_success(tmp_path: Path) -> None:
    """Verify tme_datasets.preprocessing.normalize_dataset_to_ensembl returns Result monad."""
    mock_ens = MockEnsemblRelease(111)
    adata = ad.AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["CD8A"]),
    )

    with patch("gene_utils._gene_utils.ensure_ensembl_release_installed", return_value=mock_ens):
        res = normalize_dataset_to_ensembl(adata, release=111, ensembl_dir=tmp_path)
        assert isinstance(res, Success)
        norm_adata = res.unwrap()
        assert norm_adata.var_names[0] == "ENSG00000153563"
        assert "contig" in norm_adata.var.columns
        assert "biotype" in norm_adata.var.columns
        assert "gene_name" in norm_adata.var.columns
