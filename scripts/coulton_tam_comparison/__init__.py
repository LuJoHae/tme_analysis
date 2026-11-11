"""Functional comparison of Consensus-Meta42 states against Coulton et al. (2024) TAM substates and ICI outcomes."""

from __future__ import annotations

from .models import (
    AssaySource,
    BiclusterResult,
    BulkIcbCohortBundle,
    ComparisonConfig,
    ComparisonOutputs,
    CoultonCentroidAccumulator,
    CoultonCentroidBundle,
    CoultonClusterRole,
    CoultonLiteratureSpec,
    DendrogramBracket,
    IciOutcomeDirection,
    LineageCompartment,
    Meta42MacrophageTier,
    Meta42ReferenceBundle,
    OutcomeConcordanceClass,
    SignificanceTier,
    TissueStratum,
)

__all__ = [
    "AssaySource",
    "BiclusterResult",
    "BulkIcbCohortBundle",
    "ComparisonConfig",
    "ComparisonOutputs",
    "CoultonCentroidAccumulator",
    "CoultonCentroidBundle",
    "CoultonClusterRole",
    "CoultonLiteratureSpec",
    "DendrogramBracket",
    "IciOutcomeDirection",
    "LineageCompartment",
    "Meta42MacrophageTier",
    "Meta42ReferenceBundle",
    "OutcomeConcordanceClass",
    "SignificanceTier",
    "TissueStratum",
]
