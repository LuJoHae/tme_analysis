"""Altair-based publication quality control diagnostic visualizations following Luecken & Theis (2019)."""

from __future__ import annotations

import json
from pathlib import Path
import altair as alt
import numpy as np
import polars as pl
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import QualityControlSpec

# Allow Altair to process single-cell datasets with >5,000 observations
alt.data_transformers.disable_max_rows()

logger = get_logger("preprocessing.qc_plots")


def build_qc_dashboard(
    metrics_df: pl.DataFrame,
    thresholds: QualityControlSpec,
    dataset_name: str = "Single-Cell Dataset",
) -> alt.VConcatChart:
    """Construct the 5-panel Luecken & Theis (2019) Figure 2 QC dashboard using Altair.

    Panels:
        Panel A: Histogram of count depth per cell (total_counts) with threshold line(s).
        Panel B: Histogram of detected genes per cell (n_genes_by_counts) with threshold lines.
        Panel C: Barcode rank vs. count depth (log-log elbow plot) showing droplet inflection.
        Panel D: Joint bivariate scatter of count depth vs. detected genes, colored by mito %.
        Panel E: Histogram of mitochondrial count percentage (pct_counts_mt) with max cutoff.

    Args:
        metrics_df: Polars DataFrame containing 'total_counts', 'n_genes_by_counts', 'pct_counts_mt'.
        thresholds: QualityControlSpec containing chosen thresholds.
        dataset_name: Name of the dataset for plot titles.

    Returns:
        Altair VConcatChart dashboard ready for display or SVG export.
    """
    total_cells = len(metrics_df)
    # Downsample points for scatter plot if >10,000 cells to prevent huge SVG file size
    scatter_df = metrics_df.sample(n=10000, seed=42) if total_cells > 10000 else metrics_df

    # Convert to pandas for Altair compatibility
    df_pd = metrics_df.to_pandas()
    scatter_pd = scatter_df.to_pandas()

    # 1. Panel A: Count depth histogram
    count_hist = (
        alt.Chart(df_pd)
        .mark_bar(opacity=0.7, color="#4c78a8")
        .encode(
            x=alt.X("total_counts:Q", bin=alt.Bin(maxbins=40), title="Count depth (total_counts)"),
            y=alt.Y("count():Q", title="Number of cells"),
            tooltip=[alt.Tooltip("count():Q", title="Cells")],
        )
        .properties(title="A. Count Depth Distribution", width=360, height=220)
    )

    rules_a = [{"threshold": float(thresholds.min_counts_per_cell), "label": f"Min: {thresholds.min_counts_per_cell}"}]
    if isinstance(thresholds.max_counts_per_cell, Some):
        rules_a.append({"threshold": float(thresholds.max_counts_per_cell.unwrap()), "label": f"Max: {thresholds.max_counts_per_cell.unwrap()}"})
    df_rules_a = pl.DataFrame(rules_a).to_pandas()

    rule_chart_a = (
        alt.Chart(df_rules_a)
        .mark_rule(color="#d62728", strokeDash=[4, 4], size=2)
        .encode(x="threshold:Q")
    )
    panel_a = count_hist + rule_chart_a

    # 2. Panel B: Detected genes histogram
    gene_hist = (
        alt.Chart(df_pd)
        .mark_bar(opacity=0.7, color="#2ca02c")
        .encode(
            x=alt.X("n_genes_by_counts:Q", bin=alt.Bin(maxbins=40), title="Genes detected (n_genes)"),
            y=alt.Y("count():Q", title="Number of cells"),
            tooltip=[alt.Tooltip("count():Q", title="Cells")],
        )
        .properties(title="B. Detected Genes Distribution", width=360, height=220)
    )

    df_rules_b = pl.DataFrame([
        {"threshold": float(thresholds.min_genes_per_cell), "label": f"Min: {thresholds.min_genes_per_cell}"},
        {"threshold": float(thresholds.max_genes_per_cell), "label": f"Max: {thresholds.max_genes_per_cell}"},
    ]).to_pandas()

    rule_chart_b = (
        alt.Chart(df_rules_b)
        .mark_rule(color="#d62728", strokeDash=[4, 4], size=2)
        .encode(x="threshold:Q")
    )
    panel_b = gene_hist + rule_chart_b

    # 3. Panel C: Barcode rank vs. Count depth (log-log elbow plot)
    sorted_counts = np.sort(df_pd["total_counts"].values)[::-1]
    rank_df = pl.DataFrame({
        "barcode_rank": np.arange(1, len(sorted_counts) + 1),
        "count_depth": sorted_counts,
    }).to_pandas()

    # Subsample elbow curve for clean vector rendering if very large
    if len(rank_df) > 5000:
        step = max(1, len(rank_df) // 2000)
        rank_df = rank_df.iloc[::step]

    elbow_line = (
        alt.Chart(rank_df)
        .mark_line(color="#1f77b4", size=2)
        .encode(
            x=alt.X("barcode_rank:Q", scale=alt.Scale(type="log"), title="Barcode rank (log scale)"),
            y=alt.Y("count_depth:Q", scale=alt.Scale(type="log"), title="Count depth (log scale)"),
            tooltip=["barcode_rank:Q", "count_depth:Q"],
        )
        .properties(title="C. Barcode Rank vs. Count Depth (Elbow)", width=360, height=220)
    )

    df_rule_c = pl.DataFrame([{"threshold": float(thresholds.min_counts_per_cell)}]).to_pandas()
    rule_chart_c = alt.Chart(df_rule_c).mark_rule(color="#d62728", strokeDash=[4, 4], size=2).encode(y="threshold:Q")
    panel_c = elbow_line + rule_chart_c

    # 4. Panel E: Mitochondrial percentage histogram
    mt_hist = (
        alt.Chart(df_pd)
        .mark_bar(opacity=0.7, color="#e377c2")
        .encode(
            x=alt.X("pct_counts_mt:Q", bin=alt.Bin(maxbins=40), title="Mitochondrial counts %"),
            y=alt.Y("count():Q", title="Number of cells"),
            tooltip=[alt.Tooltip("count():Q", title="Cells")],
        )
        .properties(title="E. Mitochondrial Fraction Distribution", width=360, height=220)
    )

    df_rule_e = pl.DataFrame([{"threshold": float(thresholds.max_pct_mitochondrial)}]).to_pandas()
    rule_chart_e = alt.Chart(df_rule_e).mark_rule(color="#d62728", strokeDash=[4, 4], size=2).encode(x="threshold:Q")
    panel_e = mt_hist + rule_chart_e

    # 5. Panel D: Bivariate joint scatter plot (Counts vs Genes colored by Mito %)
    scatter_chart = (
        alt.Chart(scatter_pd)
        .mark_circle(size=18, opacity=0.6)
        .encode(
            x=alt.X("total_counts:Q", title="Count depth per cell"),
            y=alt.Y("n_genes_by_counts:Q", title="Number of genes detected"),
            color=alt.Color("pct_counts_mt:Q", scale=alt.Scale(scheme="viridis"), title="Mito %"),
            tooltip=["total_counts:Q", "n_genes_by_counts:Q", "pct_counts_mt:Q"],
        )
        .properties(
            title=f"D. Bivariate Joint QC Distribution ({dataset_name}) - Colored by Mito %",
            width=760,
            height=300,
        )
    )

    # Add bounding threshold rules to bivariate plot
    x_rule = alt.Chart(df_rules_a).mark_rule(color="#d62728", strokeDash=[4, 4], size=1.5).encode(x="threshold:Q")
    y_rule = alt.Chart(df_rules_b).mark_rule(color="#d62728", strokeDash=[4, 4], size=1.5).encode(y="threshold:Q")
    panel_d = scatter_chart + x_rule + y_rule

    # Compose 5-panel dashboard
    top_row = alt.hconcat(panel_a, panel_b)
    mid_row = alt.hconcat(panel_c, panel_e)
    dashboard = alt.vconcat(top_row, mid_row, panel_d).properties(
        title=alt.TitleParams(
            text=f"Quality Control Diagnostics & Thresholds (Luecken & Theis 2019): {dataset_name}",
            fontSize=16,
            anchor="middle",
        )
    )
    return dashboard


def export_qc_plots(
    chart: alt.Chart | alt.VConcatChart | alt.HConcatChart,
    output_dir: Path,
    dataset_name: str = "dataset",
) -> Result[tuple[Path, Path], str]:
    """Export Altair QC dashboard chart to vector SVG and PNG via vl-convert.

    Args:
        chart: Altair chart or compound chart to export.
        output_dir: Directory where SVG and PNG files will be written.
        dataset_name: Dataset name identifier for file naming.

    Returns:
        Success((svg_path, png_path)) or Failure(err).
    """
    try:
        import vl_convert as vlc

        output_dir.mkdir(parents=True, exist_ok=True)
        svg_path = output_dir / f"{dataset_name}_qc_dashboard.svg"
        png_path = output_dir / f"{dataset_name}_qc_dashboard.png"

        vl_spec = chart.to_dict()
        svg_str = vlc.vegalite_to_svg(vl_spec)
        svg_path.write_text(svg_str, encoding="utf-8")

        png_bytes = vlc.vegalite_to_png(vl_spec, scale=2.0)
        png_path.write_bytes(png_bytes)

        logger.info("Successfully exported QC dashboard to %s and %s", svg_path.name, png_path.name)
        return Success((svg_path, png_path))
    except Exception as exc:
        msg = f"Failed to export QC plots for {dataset_name}: {exc}"
        logger.warning(msg)
        return Failure(msg)


def display_qc_plots_inline(chart: alt.Chart | alt.VConcatChart | alt.HConcatChart) -> None:
    """Attempt to render the Altair chart inline if running inside a Jupyter notebook or IPython."""
    try:
        from IPython import get_ipython

        ip = get_ipython()
        if ip is not None:
            from IPython.display import display

            display(chart)
    except Exception as exc:
        logger.debug("Inline plot display not available: %s", exc)


__all__ = [
    "build_qc_dashboard",
    "display_qc_plots_inline",
    "export_qc_plots",
]
