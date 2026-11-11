"""Pure output and upstream artifact resolution helpers for scripts/manuscript_pipeline/."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Final

from returns.maybe import Maybe, Nothing, Some

PIPELINE_OUTPUT_ROOT: Final[Path] = Path("output/manuscript_pipeline")
FORBIDDEN_STORE_PREFIX: Final[str] = "/capstor/store"

MULTIRES_REF_STAGES: Final[tuple[str, ...]] = (
    "3f1_build_sf_reference",
    "3f2_build_integrated_reference",
    "3f3_build_dataset_reference",
    "3f4_build_combined_references",
)

MULTIRES_ALL_UPSTREAM_STAGES: Final[tuple[str, ...]] = (
    "3f1_build_sf_reference",
    "3f2_build_integrated_reference",
    "3f3_build_dataset_reference",
    "3f4_build_combined_references",
    "3f5_deconvolute_iatlas_multires",
    "3f6_logistic_regression_multires",
    "3f7_analyze_milopy_multires",
    "3f8_compare_concordance_multires",
)


def pipeline_output_root() -> Path:
    """Return the common output directory root for all manuscript_pipeline scripts.

    Honors `TME_PIPELINE_OUTPUT_ROOT` when set (e.g. pointing to `$SCRATCH` on CSCS Eiger)
    and rejects direct `/capstor/store` (`$STORE`) paths during compute jobs.
    """
    env_override = os.environ.get("TME_PIPELINE_OUTPUT_ROOT", "").strip()
    root = Path(env_override) if env_override else PIPELINE_OUTPUT_ROOT
    if str(root).startswith(FORBIDDEN_STORE_PREFIX):
        raise RuntimeError(
            f"Direct compute output to persistent $STORE ('{root}') is strictly prohibited. "
            "Write outputs to $SCRATCH and persist via 'make stage-out'."
        )
    return root


def script_output_dir(script_file: str | Path) -> Path:
    """Return `<pipeline_output_root>/<script_stem>` for the given script file."""
    return pipeline_output_root() / Path(script_file).stem


def stage_output_dir(stage_stem: str) -> Path:
    """Return `<pipeline_output_root>/<stage_stem>` for an upstream stage script."""
    return pipeline_output_root() / stage_stem


def stage_search_dirs(
    stage_stems: tuple[str, ...],
    explicit_dirs: tuple[Path, ...] = (),
    fallback_dirs: tuple[Path, ...] = (),
) -> tuple[Path, ...]:
    """Build a deduplicated tuple of existing/candidate directories in priority order."""
    ordered = (
        *explicit_dirs,
        *(stage_output_dir(stem) for stem in stage_stems),
        *fallback_dirs,
    )
    seen: set[Path] = set()
    unique: list[Path] = []
    for candidate in ordered:
        resolved = candidate.resolve() if candidate.exists() else candidate
        if resolved not in seen:
            seen.add(resolved)
            unique.append(candidate)
    return tuple(unique)


def resolve_upstream_artifact(
    filename: str,
    stage_stems: tuple[str, ...],
    explicit_dirs: tuple[Path, ...] = (),
    fallback_dirs: tuple[Path, ...] = (),
) -> Maybe[Path]:
    """Return Some(path) for the first upstream directory containing `filename`, else Nothing."""
    for directory in stage_search_dirs(stage_stems, explicit_dirs, fallback_dirs):
        candidate = directory / filename
        if candidate.exists():
            return Some(candidate)
    return Nothing


def discover_upstream_glob(
    pattern: str,
    stage_stems: tuple[str, ...],
    explicit_dirs: tuple[Path, ...] = (),
    fallback_dirs: tuple[Path, ...] = (),
) -> tuple[Path, ...]:
    """Collect deduplicated matching files by basename across upstream stage directories."""
    by_name: dict[str, Path] = {}
    for directory in stage_search_dirs(stage_stems, explicit_dirs, fallback_dirs):
        if directory.exists():
            for match_path in sorted(directory.glob(pattern)):
                if match_path.name not in by_name:
                    by_name[match_path.name] = match_path
    return tuple(by_name[k] for k in sorted(by_name))
