"""Declarative, config-driven path resolution for tme_analysis data assets."""

from __future__ import annotations

from functools import lru_cache
import os
from pathlib import Path
import tomllib
from typing import Mapping
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some

from .logging import get_logger

logger = get_logger("paths")


class DataPathsConfig(BaseModel):
    """Immutable data paths configuration model."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    repo_root: Path
    data_root: Path
    raw_dir: Path
    preprocessed_dir: Path
    dataset_papers_dir: Path
    scratch_dir: Path
    reference_h5ad: Path
    preprocessed_template: str = "{dataset_id}.h5ad"
    raw_dataset_template: str = "{dataset_id}"
    candidate_patterns: tuple[str, ...] = (
        "data/preprocessed/{dataset_id}.h5ad",
        "dataset_papers/{dataset_id}.h5ad",
        "data/{dataset_id}.h5ad",
        "scratch/lair/ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized/{dataset_id}.h5ad",
        "scratch/lair/CBioPortalDataset-{dataset_id}/{dataset_id}.h5ad",
        "scratch/{dataset_id}/{dataset_id}.h5ad",
        "data/raw/{dataset_id}/{dataset_id}.h5ad",
    )
    overrides: Mapping[str, str] = Field(default_factory=dict)


def find_repo_root(start_dir: Path | None = None) -> Path:
    """Find repository root by walking upward until config or marker file is found."""
    current = (start_dir or Path.cwd()).resolve()
    for parent in [current] + list(current.parents):
        if (parent / "config/data_paths.toml").exists():
            return parent
        if (parent / "pyproject.toml").exists() and (parent / "packages").exists():
            return parent
        if (parent / ".git").exists():
            return parent
    return current


@lru_cache(maxsize=4)
def get_data_paths(
    repo_root: Path | None = None,
    config_path: Path | None = None,
) -> DataPathsConfig:
    """Load configuration from config/data_paths.toml, with deterministic fallback.

    Args:
        repo_root: Optional root directory of the repository. If None, auto-detected.
        config_path: Optional explicit path to the TOML configuration file.

    Returns:
        Immutable DataPathsConfig containing resolved directory paths.
    """
    root = (repo_root or find_repo_root()).resolve()

    # Determine configuration file location
    cfg_file = config_path
    if cfg_file is None:
        env_cfg = os.environ.get("TME_DATA_PATHS_CONFIG")
        if env_cfg:
            cfg_file = Path(env_cfg)
        else:
            candidate = root / "config/data_paths.toml"
            if candidate.exists():
                cfg_file = candidate

    raw_cfg: dict = {}
    if cfg_file and Path(cfg_file).exists():
        try:
            with open(cfg_file, "rb") as f:
                raw_cfg = tomllib.load(f)
            logger.debug("Loaded data paths configuration from %s", cfg_file)
        except Exception as exc:
            logger.warning("Failed to parse %s: %s. Using default paths.", cfg_file, exc)

    paths_sec = raw_cfg.get("paths", {})
    patterns_sec = paths_sec.get("candidate_h5ad_search_patterns", {})
    candidate_patterns = tuple(patterns_sec.get("patterns", [
        "data/preprocessed/{dataset_id}.h5ad",
        "dataset_papers/{dataset_id}.h5ad",
        "data/{dataset_id}.h5ad",
        "scratch/lair/ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized/{dataset_id}.h5ad",
        "scratch/lair/CBioPortalDataset-{dataset_id}/{dataset_id}.h5ad",
        "scratch/{dataset_id}/{dataset_id}.h5ad",
        "data/raw/{dataset_id}/{dataset_id}.h5ad",
    ]))

    overrides = dict(paths_sec.get("overrides", {}))

    def _resolve(rel_or_abs: str) -> Path:
        p = Path(rel_or_abs)
        return p if p.is_absolute() else (root / p).resolve()

    return DataPathsConfig(
        repo_root=root,
        data_root=_resolve(paths_sec.get("data_root", "data")),
        raw_dir=_resolve(paths_sec.get("raw_dir", "data/raw")),
        preprocessed_dir=_resolve(paths_sec.get("preprocessed_dir", "data/preprocessed")),
        dataset_papers_dir=_resolve(paths_sec.get("dataset_papers_dir", "dataset_papers")),
        scratch_dir=_resolve(paths_sec.get("scratch_dir", "scratch")),
        reference_h5ad=_resolve(paths_sec.get("reference_h5ad", "data/reference.h5ad")),
        preprocessed_template=paths_sec.get("preprocessed_template", "{dataset_id}.h5ad"),
        raw_dataset_template=paths_sec.get("raw_dataset_template", "{dataset_id}"),
        candidate_patterns=candidate_patterns,
        overrides=overrides,
    )


def get_preprocessed_h5ad_path(dataset_id: str, repo_root: Path | None = None) -> Path:
    """Return the canonical target path for caching a dataset as H5AD."""
    cfg = get_data_paths(repo_root=repo_root)
    filename = cfg.preprocessed_template.format(dataset_id=dataset_id)
    return cfg.preprocessed_dir / filename


def get_raw_dataset_dir(dataset_id: str, repo_root: Path | None = None) -> Path:
    """Return the canonical directory for raw files belonging to a dataset."""
    cfg = get_data_paths(repo_root=repo_root)
    # Check if an explicit override exists
    if dataset_id in cfg.overrides:
        override_val = cfg.overrides[dataset_id]
        p = Path(override_val)
        return p if p.is_absolute() else (cfg.repo_root / p).resolve()
    dirname = cfg.raw_dataset_template.format(dataset_id=dataset_id)
    return cfg.raw_dir / dirname


def get_scratch_dataset_dir(dataset_id: str, repo_root: Path | None = None) -> Path:
    """Return the ephemeral scratch directory for a dataset."""
    cfg = get_data_paths(repo_root=repo_root)
    return cfg.scratch_dir / dataset_id


def get_reference_h5ad_path(repo_root: Path | None = None) -> Path:
    """Return the path to the merged cross-cohort reference H5AD."""
    cfg = get_data_paths(repo_root=repo_root)
    return cfg.reference_h5ad


def find_dataset_h5ad(dataset_id: str, repo_root: Path | None = None) -> Maybe[Path]:
    """Find an existing H5AD file across configured candidate directories.

    Returns:
        Some(path) if a valid non-empty H5AD file exists, else Nothing.
    """
    cfg = get_data_paths(repo_root=repo_root)

    # 1. Check override if defined
    if dataset_id in cfg.overrides:
        override_val = cfg.overrides[dataset_id]
        override_path = Path(override_val)
        resolved = override_path if override_path.is_absolute() else (cfg.repo_root / override_path).resolve()
        if resolved.exists() and resolved.is_file() and resolved.stat().st_size > 0:
            return Some(resolved)

    # 2. Check candidate search patterns in priority order
    for pattern in cfg.candidate_patterns:
        formatted_rel = pattern.format(dataset_id=dataset_id)
        candidate = (cfg.repo_root / formatted_rel).resolve()
        if candidate.exists() and candidate.is_file() and candidate.stat().st_size > 0:
            return Some(candidate)

    return Nothing
