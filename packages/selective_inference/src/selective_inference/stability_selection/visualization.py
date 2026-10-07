"""Declarative Altair visualizations for Stability Selection.

Follows strict guidelines:
- Declarative Altair (based on Vega-Lite)
- Polars data source integration
- Okabe-Ito colorblind-safe palette
- Clean minimalist layout suitable for SVG publication export
"""

from typing import cast
import altair as alt
import numpy as np
import polars as pl
from returns.maybe import Some

from selective_inference.stability_selection.types import StabilityResult


# Okabe-Ito color palette
COLOR_SELECTED = "#D55E00"      # Vermilion / Red-Orange
COLOR_UNSELECTED = "#94A3B8"    # Slate hairline
COLOR_THRESHOLD = "#C00000"     # Dark crimson rule
COLOR_CUTOFF = "#0284C7"        # Sky blue for regularization budget cutoff


def plot_stability_paths(
    result: StabilityResult,
    title: str = "Stability Paths",
    width: int = 500,
    height: int = 350,
) -> alt.LayerChart:
    """Create an interactive Altair chart of the stability selection paths.

    Parameters
    ----------
    result : StabilityResult
        Fitted stability selection result.
    title : str
        Chart title.
    width : int
        Chart width in pixels.
    height : int
        Chart height in pixels.

    Returns
    -------
    alt.LayerChart
        Declarative Altair chart object.
    """
    df = result.to_path_polars()
    cutoff = result.parameters.cutoff

    # Base line chart of paths
    lines = (
        alt.Chart(df)
        .mark_line(strokeWidth=1.5, opacity=0.85)
        .encode(
            x=alt.X(
                "log_lambda:Q",
                title="log₁₀(λ) Regularization Penalty",
                scale=alt.Scale(reverse=True),
            ),
            y=alt.Y(
                "selection_probability:Q",
                title="Selection Probability Π(λ)",
                scale=alt.Scale(domain=[0, 1.05]),
            ),
            color=alt.Color(
                "selected:N",
                title="Status",
                scale=alt.Scale(
                    domain=[True, False],
                    range=[COLOR_SELECTED, COLOR_UNSELECTED],
                ),
                legend=alt.Legend(
                    labelExpr="datum.value ? 'Stable Feature' : 'Noise / Below Threshold'"
                ),
            ),
            detail="feature:N",
            tooltip=[
                alt.Tooltip("feature:N", title="Feature"),
                alt.Tooltip("selection_probability:Q", format=".3f", title="Probability"),
                alt.Tooltip("lambda:Q", format=".4e", title="λ"),
            ],
        )
    )

    # Threshold horizontal reference line
    threshold_df = pl.DataFrame({"cutoff": [cutoff], "label": [f"Threshold π_thr = {cutoff:.2f}"]})
    rule = (
        alt.Chart(threshold_df)
        .mark_rule(
            color=COLOR_THRESHOLD,
            strokeDash=[5, 5],
            strokeWidth=1.5,
        )
        .encode(y="cutoff:Q")
    )

    text = (
        alt.Chart(threshold_df)
        .mark_text(
            align="right",
            dx=-10,
            dy=-8,
            fontSize=11,
            color=COLOR_THRESHOLD,
            fontWeight="bold",
        )
        .encode(
            y="cutoff:Q",
            text="label:N",
        )
    )

    layers: list[alt.Chart] = [lines, rule, text]
    subtitle_elements = [
        f"Selected {len(result.selected_features)} / {result.parameters.p} features",
        f"PFER bound ≤ {result.parameters.pfer:.2f} ({result.parameters.assumption})",
    ]

    # Vertical regularization budget cutoff line
    match result.lambda_cutoff:
        case Some(lam_cut) if lam_cut > 0:
            log_cut = float(np.log10(lam_cut))
            label_text = f"λ cutoff = {lam_cut:.4g}"
            cutoff_df = pl.DataFrame({
                "log_lambda": [log_cut],
                "label": [label_text],
            })
            cut_rule = (
                alt.Chart(cutoff_df)
                .mark_rule(
                    color=COLOR_CUTOFF,
                    strokeDash=[4, 3],
                    strokeWidth=1.5,
                )
                .encode(x="log_lambda:Q")
            )
            cut_text = (
                alt.Chart(cutoff_df)
                .mark_text(
                    align="right",
                    dx=-6,
                    dy=10,
                    fontSize=10,
                    color=COLOR_CUTOFF,
                    fontWeight="bold",
                )
                .encode(
                    x="log_lambda:Q",
                    y=alt.value(10),
                    text="label:N",
                )
            )
            layers.extend([cut_rule, cut_text])
            subtitle_elements.insert(0, f"λ cutoff = {lam_cut:.4g} (q budget ≤ {result.parameters.q:.1f})")
        case _:
            pass

    return cast(
        alt.LayerChart,
        alt.layer(*layers)
        .properties(
            title=alt.TitleParams(
                text=title,
                subtitle=[" | ".join(subtitle_elements)],
                anchor="start",
            ),
            width=width,
            height=height,
        )
        .interactive(),
    )


def plot_stability_scores(
    result: StabilityResult,
    max_features: int = 25,
    title: str = "Feature Stability Scores",
    width: int = 400,
    height: int = 400,
) -> alt.LayerChart:
    """Create a horizontal bar chart of the highest stability scores.

    Parameters
    ----------
    result : StabilityResult
        Fitted stability selection result.
    max_features : int
        Maximum number of top features to display.
    title : str
        Chart title.

    Returns
    -------
    alt.Chart
    """
    df = result.to_polars().head(max_features)
    cutoff = result.parameters.cutoff

    bars = (
        alt.Chart(df)
        .mark_bar(cornerRadiusEnd=0)
        .encode(
            x=alt.X(
                "stability_score:Q",
                title="Max Selection Probability Π_max",
                scale=alt.Scale(domain=[0, 1.05]),
            ),
            y=alt.Y(
                "feature:N",
                title="Feature",
                sort="-x",
            ),
            color=alt.Color(
                "selected:N",
                title="Selected",
                scale=alt.Scale(
                    domain=[True, False],
                    range=[COLOR_SELECTED, COLOR_UNSELECTED],
                ),
            ),
            tooltip=[
                alt.Tooltip("feature:N", title="Feature"),
                alt.Tooltip("stability_score:Q", format=".3f", title="Stability Score"),
                alt.Tooltip("selected:N", title="Is Selected"),
            ],
        )
    )

    threshold_df = pl.DataFrame({"cutoff": [cutoff]})
    rule = (
        alt.Chart(threshold_df)
        .mark_rule(
            color=COLOR_THRESHOLD,
            strokeDash=[5, 5],
            strokeWidth=1.5,
        )
        .encode(x="cutoff:Q")
    )

    return cast(
        alt.LayerChart,
        (bars + rule)
        .properties(
            title=alt.TitleParams(text=title, anchor="start"),
            width=width,
            height=height,
        ),
    )
