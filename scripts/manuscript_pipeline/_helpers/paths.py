"""Pure output directory resolution helpers for scripts/manuscript_pipeline/."""

from __future__ import annotations

from pathlib import Path
from typing import Final

PIPELINE_OUTPUT_ROOT: Final[Path] = Path("output/manuscript_pipeline")


def pipeline_output_root() -> Path:
    """Return the common output directory root for all manuscript_pipeline scripts."""
    return PIPELINE_OUTPUT_ROOT


def script_output_dir(script_file: str | Path) -> Path:
    """Return `output/manuscript_pipeline/<script_stem>` for the given script file."""
    return PIPELINE_OUTPUT_ROOT / Path(script_file).stem


def stage_output_dir(stage_stem: str) -> Path:
    """Return `output/manuscript_pipeline/<stage_stem>` for an upstream stage script."""
    return PIPELINE_OUTPUT_ROOT / stage_stem
