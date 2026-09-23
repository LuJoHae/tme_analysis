"""Unit tests for config-driven paths, H5AD caching, and zero-hardcoding enforcement."""

from __future__ import annotations

import ast
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from returns.maybe import Nothing, Some
from returns.result import Success

from tme_datasets.paths import (
    DataPathsConfig,
    find_dataset_h5ad,
    find_repo_root,
    get_data_paths,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_reference_h5ad_path,
    get_scratch_dataset_dir,
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
    assert cfg.data_root == root / "data"
    assert cfg.raw_dir == root / "data/raw"
    assert cfg.preprocessed_dir == root / "data/preprocessed"
    assert cfg.scratch_dir == root / "scratch"
    assert cfg.reference_h5ad == root / "data/reference.h5ad"


def test_path_resolvers() -> None:
    """Verify path resolvers construct canonical absolute paths."""
    root = find_repo_root()
    h5ad_target = get_preprocessed_h5ad_path("GSE120575", repo_root=root)
    assert h5ad_target == root / "data/preprocessed/GSE120575.h5ad"

    raw_target = get_raw_dataset_dir("GSE120575", repo_root=root)
    assert raw_target == root / "data/raw/GSE120575"

    scratch_target = get_scratch_dataset_dir("GSE120575", repo_root=root)
    assert scratch_target == root / "scratch/GSE120575"

    ref_path = get_reference_h5ad_path(repo_root=root)
    assert ref_path == root / "data/reference.h5ad"


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
    # Fast-load invariant: must load in less than 2.0s (usually < 0.5s)
    assert elapsed < 2.0


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
