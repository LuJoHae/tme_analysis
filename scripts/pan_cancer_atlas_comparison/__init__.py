"""Pan-Cancer T-Cell, CAF & Myeloid Triad Atlas Comparison Suite (`Consensus-Meta42` Validation).

Strictly adheres to functional Python standards (.agents/rules/code-style-guide.md):
- Functional Core (`core.py`, `plotting.py`) / Imperative Shell (`compare_meta42_external_atlases.py`, `render_figures.py`)
- Frozen Pydantic models (`models.py`), `typing.Protocol` typeclasses, exhaustive `match...case`
- Zero `None` / `Optional` / `del` / `global` / `nonlocal` (`returns.maybe.Maybe` & `returns.result.Result`)
- 100% `tme_datasets` public API exclusivity for all single-cell and bulk RNA-seq dataset loading
"""

from __future__ import annotations

from .models import (
    AtlasCentroidAccumulator,
    AtlasCentroidBundle,
    BulkIcbCohortBundle,
    ExternalAtlasId,
    IciDirection,
    PanAtlasConfig,
    PanAtlasOutputs,
    ProfileOrigin,
    SignificanceTier,
    StromalGatingResult,
    TriadCompartment,
)

__all__ = [
    "AtlasCentroidAccumulator",
    "AtlasCentroidBundle",
    "BulkIcbCohortBundle",
    "ExternalAtlasId",
    "IciDirection",
    "PanAtlasConfig",
    "PanAtlasOutputs",
    "ProfileOrigin",
    "SignificanceTier",
    "StromalGatingResult",
    "TriadCompartment",
]
