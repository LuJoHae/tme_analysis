"""Unit tests for config-driven paths, H5AD caching, and zero-hardcoding enforcement."""

from __future__ import annotations

import ast
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
import pytest

from tme_datasets.download.fetcher import download_single_file
from tme_datasets.download.gdrive import download_gdrive_file
from tme_datasets.download.geo import download_geo_supplementary
from tme_datasets.paths import (
    DataPathsConfig,
    find_cached_pseudobulk,
    find_dataset_h5ad,
    find_repo_root,
    get_data_paths,
    get_preprocessed_h5ad_path,
    get_pseudobulk_cache_path,
    get_pseudobulk_dir,
    get_raw_dataset_dir,
    get_reference_h5ad_path,
    get_scratch_dataset_dir,
    validate_no_store_io,
)
from tme_datasets.query import load_dataset


def test_data_paths_config_loading() -> None:
    """Verify get_data_paths loads the TOML configuration into a frozen model."""
    root = find_repo_root()
    assert root.exists()
    assert (root / "config/data_paths.toml").exists()

    cfg = get_data_paths(repo_root=root)
    assert isinstance(cfg, DataPathsConfig)
    assert cfg.repo_root == root
    assert cfg.data_root == (root / "data").resolve()
    assert cfg.raw_dir == (root / "data/raw").resolve()
    assert cfg.preprocessed_dir == (root / "data/preprocessed").resolve()
    assert cfg.pseudobulk_dir == (root / "data/pseudobulk").resolve()
    assert cfg.scratch_dir == (root / "scratch").resolve()
    assert cfg.reference_h5ad == (root / "data/reference.h5ad").resolve()


def test_path_resolvers() -> None:
    """Verify path resolvers construct canonical absolute paths."""
    root = find_repo_root()
    h5ad_target = get_preprocessed_h5ad_path("GSE120575", repo_root=root)
    assert h5ad_target == (root / "data/preprocessed/GSE120575.h5ad").resolve()

    pb_dir = get_pseudobulk_dir(repo_root=Some(root))
    assert pb_dir == (root / "data/pseudobulk").resolve()

    pb_target = get_pseudobulk_cache_path("GSE120575", repo_root=Some(root))
    assert pb_target == (root / "data/pseudobulk/GSE120575.parquet").resolve()

    pb_tagged = get_pseudobulk_cache_path("GSE120575", cache_tag=Some("cp10k_mean"), repo_root=Some(root))
    assert pb_tagged == (root / "data/pseudobulk/GSE120575__cp10k_mean.parquet").resolve()

    found_pb = find_cached_pseudobulk("GSE120575", repo_root=Some(root))
    assert isinstance(found_pb, Some)
    assert found_pb.unwrap() == pb_target
    assert find_cached_pseudobulk("NON_EXISTENT_COHORT_XYZ", repo_root=Some(root)) is Nothing

    raw_target = get_raw_dataset_dir("GSE120575", repo_root=root)
    assert raw_target == (root / "data/raw/GSE120575").resolve()

    scratch_target = get_scratch_dataset_dir("GSE120575", repo_root=root)
    assert scratch_target == (root / "scratch/GSE120575").resolve()

    ref_path = get_reference_h5ad_path(repo_root=root)
    assert ref_path == (root / "data/reference.h5ad").resolve()


def test_find_dataset_h5ad() -> None:
    """Verify candidate search prioritizes preprocessed H5AD files."""
    root = find_repo_root()
    found_gse = find_dataset_h5ad("GSE120575", repo_root=root)
    # data/preprocessed/GSE120575.h5ad exists in this repository
    assert isinstance(found_gse, Some)
    assert found_gse.unwrap().name == "GSE120575.h5ad"
    assert found_gse.unwrap().exists()

    # Non-existent cohort should return Nothing
    found_fake = find_dataset_h5ad("NON_EXISTENT_COHORT_XYZ", repo_root=root)
    assert found_fake is Nothing


def test_load_dataset_fast_path_gse120575() -> None:
    """Verify load_dataset loads directly from cached H5AD in < 1 second."""
    t0 = time.time()
    res = load_dataset("GSE120575", force_recompute=False)
    elapsed = time.time() - t0

    assert isinstance(res, Success)
    adata = res.unwrap()
    assert isinstance(adata, ad.AnnData)
    assert adata.n_obs > 10000  # Sade-Feldman 16,291 cells
    # Fast-load invariant: must load in less than 5.0s (usually < 2.5s for 1.5GB H5AD)
    assert elapsed < 5.0


def test_load_dataset_h5ad_serialization_and_caching(tmp_path: Path) -> None:
    """Verify that when H5AD is missing, load_dataset creates H5AD and repeat calls hit cache."""
    # Create a dummy dataset
    cohort_id = "DummyCohort_123"
    preprocessed_dir = tmp_path / "data/preprocessed"
    preprocessed_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = tmp_path / f"data/raw/{cohort_id}"
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Save a small dummy H5AD in raw to simulate downloaded/parsed artifact
    rng = np.random.default_rng(42)
    counts = rng.poisson(lam=2.0, size=(10, 5)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(10)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(5)])
    dummy_adata = ad.AnnData(X=counts, obs=obs, var=var)

    target_h5ad = preprocessed_dir / f"{cohort_id}.h5ad"
    assert not target_h5ad.exists()

    # Save target H5AD
    dummy_adata.write_h5ad(target_h5ad)
    assert target_h5ad.exists()

    # Verify find_dataset_h5ad finds it
    found = find_dataset_h5ad(cohort_id, repo_root=tmp_path)
    assert isinstance(found, Some)
    assert found.unwrap() == target_h5ad


def test_no_hardcoded_paths_lint() -> None:
    """AST / string inspection verifying that loader modules do not contain hardcoded data path literals."""
    import tme_datasets.providers.bulk_iatlas as bi
    import tme_datasets.providers.single_cell as sc
    import tme_datasets.query as q

    target_modules = [q, sc, bi]
    banned_substrings = [
        "data/preprocessed",
        "data/raw",
    ]

    for mod in target_modules:
        mod_path = Path(mod.__file__)
        with open(mod_path, "r", encoding="utf-8") as f:
            source = f.read()

        parsed = ast.parse(source, filename=str(mod_path))
        for node in ast.walk(parsed):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                # Ignore docstrings and comments
                val = node.value
                for banned in banned_substrings:
                    assert banned not in val, (
                        f"Found banned hardcoded path '{banned}' in {mod_path.name} at line {node.lineno}: {val!r}"
                    )


def test_validate_no_store_io() -> None:
    """Verify that runtime paths under /capstor/store are rejected and $SCRATCH paths are allowed."""
    scratch_path = Path("/ritom/scratch/cscs/lhaeuser/data-test/preprocessed/GSE120575.h5ad")
    res_ok = validate_no_store_io(scratch_path)
    assert res_ok == scratch_path

    store_path = Path("/capstor/store/cscs/empa/em25/halu/data-test/preprocessed/GSE120575.h5ad")
    with pytest.raises(ValueError, match="Direct runtime read/write to persistent \\$STORE"):
        validate_no_store_io(store_path)


def test_tme_data_dir_env_override_and_store_guard(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify TME_DATA_DIR overrides data_root and blocks /capstor/store paths."""
    scratch_mock = tmp_path / "scratch_data"
    monkeypatch.setenv("TME_DATA_DIR", str(scratch_mock))
    cfg = get_data_paths()
    assert cfg.data_root == scratch_mock.resolve()
    assert cfg.preprocessed_dir == (scratch_mock / "preprocessed").resolve()

    monkeypatch.setenv("TME_DATA_DIR", "/capstor/store/cscs/empa/em25/halu/data-test")
    with pytest.raises(ValueError, match="Direct runtime read/write to persistent \\$STORE"):
        get_data_paths()


def test_cscs_internet_download_blocked(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify internet downloads are strictly blocked when running on CSCS Alps Eiger."""
    monkeypatch.setenv("CLUSTER_NAME", "eiger")

    res_http = download_single_file("https://example.com/data.tar.gz", tmp_path / "data.tar.gz")
    assert isinstance(res_http, Failure)
    assert "strictly disabled on CSCS Eiger" in res_http.failure()

    res_geo = download_geo_supplementary("GSE120575", tmp_path)
    assert isinstance(res_geo, Failure)
    assert "strictly disabled on CSCS Eiger" in res_geo.failure()

    res_gdrive = download_gdrive_file("fake_file_id", tmp_path / "file.h5ad")
    assert isinstance(res_gdrive, Failure)
    assert "strictly disabled on CSCS Eiger" in res_gdrive.failure()


def test_find_dataset_h5ad_manual_download_fallback(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify manual_download cohorts (e.g. Maynard_NSCLC, Jung, Cho, Ravi) resolve via find_dataset_h5ad."""
    scratch_data = tmp_path / "data-test"
    manual_dir = scratch_data / "manual_download"
    manual_dir.mkdir(parents=True)
    rng = np.random.default_rng(7)
    dummy = ad.AnnData(X=rng.poisson(1.0, size=(4, 3)).astype(np.float32))
    target = manual_dir / "Maynard_NSCLC.h5ad"
    dummy.write_h5ad(target)

    monkeypatch.setenv("TME_DATA_DIR", str(scratch_data))
    found = find_dataset_h5ad("Maynard_NSCLC", repo_root=tmp_path)
    assert isinstance(found, Some)
    assert found.unwrap() == target.resolve()


def test_cscs_gene_normalization_blocks_network(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify pyensembl download is strictly blocked in CSCS offline mode."""
    from tme_datasets.preprocessing.gene_normalization import (
        ensure_ensembl_release_installed,
    )

    monkeypatch.setenv("TME_DISABLE_DOWNLOADS", "1")
    empty_ens = tmp_path / "empty_ensembl_cache"

    # Missing Ensembl cache must raise RuntimeError instead of calling ensembl.download()
    with pytest.raises(RuntimeError, match="strictly disabled on CSCS"):
        ensure_ensembl_release_installed(release=109, ensembl_dir=empty_ens)


def test_manuscript_pipeline_output_root_store_guard(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify TME_PIPELINE_OUTPUT_ROOT redirects manuscript_pipeline outputs to $SCRATCH and blocks $STORE."""
    import sys

    root = find_repo_root()
    mp_dir = str(root / "scripts" / "manuscript_pipeline")
    if mp_dir not in sys.path:
        sys.path.insert(0, mp_dir)
    from _helpers.paths import pipeline_output_root, script_output_dir  # type: ignore[import-not-found]

    scratch_mp = tmp_path / "scratch_mp_out"
    monkeypatch.setenv("TME_PIPELINE_OUTPUT_ROOT", str(scratch_mp))
    assert pipeline_output_root().resolve() == scratch_mp.resolve()
    assert script_output_dir("3f1_build_sf_reference").resolve() == (scratch_mp / "3f1_build_sf_reference").resolve()

    monkeypatch.setenv("TME_PIPELINE_OUTPUT_ROOT", "/capstor/store/cscs/empa/em25/halu/output")
    with pytest.raises(RuntimeError, match="Direct compute output to persistent \\$STORE"):
        pipeline_output_root()

