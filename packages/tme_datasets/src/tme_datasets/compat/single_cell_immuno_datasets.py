"""Compatibility shim mirroring legacy single_cell_immuno_datasets interface."""

from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some

from ..models import DatasetSpec as TmeDatasetSpec
from ..models import QualityControlSpec as TmeQCSpec


class DataDirectories(BaseModel):
    """Legacy DataDirectories model."""

    model_config = ConfigDict(frozen=True)

    base_dir: Path
    raw_dir: Path
    preprocessed_dir: Path
    reports_dir: Path

    @classmethod
    def from_base(cls, base_dir: Path | str) -> DataDirectories:
        p = Path(base_dir).resolve()
        return cls(
            base_dir=p,
            raw_dir=p / "raw",
            preprocessed_dir=p / "preprocessed",
            reports_dir=p / "reports",
        )


TIER_1_DATASETS = (
    "GSE120575",
    "GSE115978",
    "GSE123139",
    "GSE123813",
    "GSE125449",
    "GSE179994",
    "GSE159115",
    "GSE171306",
    "Gondal2025",
)
