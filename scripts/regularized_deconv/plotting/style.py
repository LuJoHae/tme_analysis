#!/usr/bin/env python3
"""
Plotting style adapter for deconvolution benchmark figures.

Integrates the in-house 'plotting_utils' package:
- Re-exports Nature Methods themes and axes.
- Re-exports Okabe-Ito deconvolution palettes (DECONV_METHOD_ORDER, DECONV_COLOR_MAP).
- Provides standardized export_altair_figure.
"""

from __future__ import annotations

import altair as alt  # type: ignore
from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    DECONV_COLOR_MAP,
    DECONV_METHOD_ORDER,
    OKABE_PALETTE,
    apply_nature_methods_theme,
    create_nature_axis,
    export_altair_figure,
)

__all__ = [
    "COLOR_BG_WHITE",
    "COLOR_HAIRLINE",
    "COLOR_LIGHT_GREY",
    "COLOR_TEXT_PRIMARY",
    "COLOR_TEXT_SECONDARY",
    "DECONV_COLOR_MAP",
    "DECONV_METHOD_ORDER",
    "OKABE_PALETTE",
    "apply_nature_methods_theme",
    "create_nature_axis",
    "export_altair_figure",
]
