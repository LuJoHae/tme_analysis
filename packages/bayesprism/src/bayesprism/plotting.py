from typing import Sequence
import numpy as np
import polars as pl
import altair as alt
import vl_convert as vlc
from returns.result import Result, Success, Failure


def plot_cor_phi(phi: np.ndarray, cell_names: Sequence[str]) -> alt.Chart:
    """
    Generate correlation heatmap between reference cell state/type expression profiles.
    Returns Altair Chart object.
    """
    cor_matrix = np.corrcoef(phi)
    K = len(cell_names)

    records = []
    for i in range(K):
        for j in range(K):
            records.append({
                "cell_type_1": cell_names[i],
                "cell_type_2": cell_names[j],
                "correlation": float(cor_matrix[i, j]),
            })

    df = pl.DataFrame(records)

    chart = (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X("cell_type_1:N", title="Cell Type 1"),
            y=alt.Y("cell_type_2:N", title="Cell Type 2"),
            color=alt.Color("correlation:Q", scale=alt.Scale(scheme="viridis"), title="Pearson r"),
            tooltip=["cell_type_1", "cell_type_2", "correlation"],
        )
        .properties(
            title="Reference Expression Profile Correlation",
            width=400,
            height=400,
        )
    )

    return chart


def export_chart_svg(chart: alt.Chart, output_filepath: str) -> Result[str, str]:
    """
    Export Altair chart to scalable vector graphic (SVG) file using vl-convert.
    """
    try:
        vega_spec = chart.to_dict()
        svg_str = vlc.vegalite_to_svg(vega_spec)
        with open(output_filepath, "w", encoding="utf-8") as f:
            f.write(svg_str)
        return Success(output_filepath)
    except Exception as e:
        return Failure(f"Failed to export chart to SVG: {str(e)}")
