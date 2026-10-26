"""Altair declarative ROC and PR curve visualizations with SVG export."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Sequence
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


def export_chart_svg(
    chart: alt.Chart | alt.HConcatChart | alt.VConcatChart | alt.LayerChart | Any,
    output_svg_path: Path,
    target_width_mm: float | None = None,
    min_font_size: float | None = None,
) -> Result[Path, str]:
    """Export Altair chart directly to high-resolution SVG using vl-convert with optional mm scaling and min font size."""
    try:
        output_svg_path.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(output_svg_path), format="svg")

        if target_width_mm is not None or min_font_size is not None:
            svg_content = output_svg_path.read_text(encoding="utf-8")

            if min_font_size is not None:
                svg_content = re.sub(
                    r'font-size="([0-9\.]+)px"',
                    lambda m: f'font-size="{max(min_font_size, float(m.group(1)))}px"',
                    svg_content,
                )
                svg_content = re.sub(
                    r'font-size:\s*([0-9\.]+)px',
                    lambda m: f'font-size: {max(min_font_size, float(m.group(1)))}px',
                    svg_content,
                )

            if target_width_mm is not None:
                match_dims = re.search(
                    r'<svg([^>]+)width="([0-9\.]+)"\s+height="([0-9\.]+)"([^>]*)viewBox="0\s+0\s+([0-9\.]+)\s+([0-9\.]+)"',
                    svg_content,
                )
                if match_dims:
                    w = float(match_dims.group(2))
                    h = float(match_dims.group(3))
                    h_mm = target_width_mm * h / w
                    prefix = f'<svg{match_dims.group(1)}width="{target_width_mm:.2f}mm" height="{h_mm:.2f}mm"{match_dims.group(4)}viewBox="0 0 {w} {h}"'
                    svg_content = svg_content[:match_dims.start()] + prefix + svg_content[match_dims.end():]

            output_svg_path.write_text(svg_content, encoding="utf-8")

        return Success(output_svg_path)
    except Exception as exc:
        return Failure(f"Failed to export chart to SVG: {exc}")
