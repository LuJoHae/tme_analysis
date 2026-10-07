"""Deconvolution reference building and modeling framework."""

from __future__ import annotations

from .adapters import (
    export_to_bayesprism,
    export_to_instaprism,
)
from .builder import build_deconvolution_reference
from .malignant import (
    calculate_cnv_proxy_scores,
    detect_malignant_cells,
)
from .models import (
    DeconvolutionReferenceConfig,
    DeconvolutionReferenceResult,
)

__all__ = [
    "build_deconvolution_reference",
    "DeconvolutionReferenceConfig",
    "DeconvolutionReferenceResult",
    "detect_malignant_cells",
    "calculate_cnv_proxy_scores",
    "export_to_bayesprism",
    "export_to_instaprism",
]
