"""Altair declarative ROC and PR curve visualizations with SVG export."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import altair as alt
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success
from sklearn.metrics import precision_recall_curve, roc_curve


def create_roc_chart(
    y_true: Sequence[float],
    y_score: Sequence[float],
    model_name: str,
    title: str = "Receiver Operating Characteristic (ROC)",
) -> Result[alt.Chart, str]:
    """Pure functional constructor of Altair ROC curve."""
    try:
        fpr, tpr, _ = roc_curve(np.asarray(y_true), np.asarray(y_score))
        df = pl.DataFrame({
            "fpr": fpr,
            "tpr": tpr,
            "model": model_name,
        })

        line = (
            alt.Chart(df)
            .mark_line(strokeWidth=2.5)
            .encode(
                x=alt.X("fpr:Q", title="False Positive Rate (1 - Specificity)", scale=alt.Scale(domain=[0, 1])),
                y=alt.Y("tpr:Q", title="True Positive Rate (Sensitivity)", scale=alt.Scale(domain=[0, 1])),
                color=alt.Color("model:N", title="Predictor"),
            )
        )

        diagonal = (
            alt.Chart(pl.DataFrame({"x": [0, 1], "y": [0, 1]}))
            .mark_line(strokeDash=[4, 4], color="gray")
            .encode(x="x:Q", y="y:Q")
        )

        chart = (line + diagonal).properties(
            title=title,
            width=360,
            height=320,
        )
        return Success(chart)
    except Exception as exc:
        return Failure(f"Failed to create ROC chart: {exc}")


def export_chart_svg(chart: alt.Chart, output_svg_path: Path) -> Result[Path, str]:
    """Export Altair chart directly to high-resolution SVG using vl-convert."""
    try:
        output_svg_path.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(output_svg_path), format="svg")
        return Success(output_svg_path)
    except Exception as exc:
        return Failure(f"Failed to export chart to SVG: {exc}")
