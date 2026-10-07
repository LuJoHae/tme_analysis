"""Unit tests for plotting_utils package.
"""

from __future__ import annotations

import altair as alt
import polars as pl
from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COHORT_METADATA,
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_PALETTE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)


def test_plotting_utils_palettes_and_cohorts() -> None:
    assert len(COHORT_PANELS) == 9
    assert len(COHORT_METADATA) == 9
    assert len(OKABE_PALETTE) >= 8
    assert OKABE_BLUE == "#0072B2"
    assert OKABE_VERMILION == "#D55E00"


def test_apply_nature_methods_theme() -> None:
    df = pl.DataFrame({"x": [1, 2, 3], "y": [10, 20, 30]})
    chart = alt.Chart(df).mark_circle().encode(x="x:Q", y="y:Q")
    themed = apply_nature_methods_theme(chart)
    assert themed is not None
    chart_dict = themed.to_dict()
    assert "config" in chart_dict
    assert chart_dict["config"]["view"]["fill"] == COLOR_BG_WHITE
