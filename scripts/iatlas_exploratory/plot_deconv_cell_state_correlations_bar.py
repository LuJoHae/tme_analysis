"""Generate publication-grade Altair vertical bar chart of deconvoluted cell state correlations with ICB response.

Reproduces the Gide et al. (2019) deconvoluted cell-state point-biserial correlation plot
from jupyter/deconv.ipynb (Cell 23) in the Nature Methods wireframe style used across
tme_response.visualization, with a two-line title and explicit legend for Cytotoxic T Cells
(top positive correlations) and Cancer-Associated Fibroblasts (most negative correlations).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
from typing import Final

import altair as alt
import h5py  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats
import vl_convert as vlc

GROUP_CYTOTOXIC_T: Final[str] = "Cytotoxic T Cells"
GROUP_CAF: Final[str] = "Cancer-Associated Fibroblasts (CAFs)"
GROUP_OTHER: Final[str] = "Other Cell States"

GROUP_ORDER: Final[tuple[str, ...]] = (
    GROUP_CYTOTOXIC_T,
    GROUP_CAF,
    GROUP_OTHER,
)

# Okabe-Ito colorblind-safe Vermillion (red), Blue, and neutral slate grey
GROUP_COLORS: Final[tuple[str, ...]] = (
    "#D55E00",  # Red / Vermillion: Cytotoxic T Cells
    "#0072B2",  # Blue: Cancer-Associated Fibroblasts (CAFs)
    "#94A3B8",  # Neutral Grey: Other Cell States
)


class CellStateCorrelationPlotConfig(BaseModel):
    """Immutable configuration for the deconvoluted cell-state correlation vertical bar chart."""

    model_config = ConfigDict(frozen=True)

    gide_h5ad: Path = Path("output/deconv_gide_correlations/gide.h5ad")
    log_fractions_npy: Path = Path("output/deconv_gide_correlations/log_fractions.npy")
    output_svg: Path = Path("output/deconv_gide_correlations/features.svg")
    output_parquet: Path = Path(
        "output/deconv_gide_correlations/gide_cell_state_correlations.parquet"
    )
    n_top_positive: int = 15
    n_bottom_negative: int = 15
    width: int = 580
    height: int = 190


def _read_gide_response_numeric(gide_h5ad: Path) -> Result[npt.NDArray[np.float64], str]:
    """Read binary response vector (NR=0.0, R=1.0, NaN otherwise) from Gide h5ad via h5py."""
    try:
        with h5py.File(gide_h5ad, "r") as f:
            resp_group = f["obs"]["response"]
            if isinstance(resp_group, h5py.Group) and "categories" in resp_group:
                cats = [
                    c.decode("utf-8") if isinstance(c, bytes) else str(c)
                    for c in resp_group["categories"][:]
                ]
                codes = resp_group["codes"][:]
                values = [cats[int(c)] if int(c) >= 0 else "" for c in codes]
            else:
                raw = f["obs"]["response"][:]
                values = [
                    v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in raw
                ]
        mapped = np.array(
            [1.0 if v == "R" else (0.0 if v == "NR" else np.nan) for v in values],
            dtype=float,
        )
        return Success(mapped)
    except Exception as exc:
        return Failure(f"Failed to read response from {gide_h5ad}: {exc}")


def compute_gide_cell_state_correlations(
    gide_h5ad: Path,
    log_fractions_npy: Path,
    n_top_positive: int = 15,
    n_bottom_negative: int = 15,
) -> Result[pl.DataFrame, str]:
    """Compute point-biserial correlations between deconvoluted cell states and ICB response."""
    if not gide_h5ad.exists():
        return Failure(f"Gide AnnData file not found: {gide_h5ad}")
    if not log_fractions_npy.exists():
        return Failure(f"Deconvoluted log-fractions file not found: {log_fractions_npy}")

    resp_res = _read_gide_response_numeric(gide_h5ad)
    if isinstance(resp_res, Failure):
        return resp_res
    y_numeric = resp_res.unwrap()

    try:
        log_fractions = np.load(log_fractions_npy)
    except Exception as exc:
        return Failure(f"Failed to load log-fractions array: {exc}")

    fractions = np.exp(log_fractions)
    x_matrix = np.log(fractions + 1e-6).astype(float)
    valid_sample_mask = ~np.isnan(y_numeric)
    y_clean = y_numeric[valid_sample_mask]

    rows: list[dict[str, object]] = []
    for idx in range(x_matrix.shape[1]):
        col_vals = x_matrix[valid_sample_mask, idx]
        if float(np.std(col_vals)) == 0.0:
            continue
        res = stats.pearsonr(y_clean, col_vals)
        r_val = float(res.statistic)
        p_val = float(res.pvalue)
        if np.isfinite(r_val) and np.isfinite(p_val):
            rows.append({
                "cell_state_index": idx,
                "cell_state_id": f"cell_type_{idx:03d}",
                "point_biserial_r": r_val,
                "p_value": p_val,
            })

    if not rows:
        return Failure("No finite cell-state correlations computed.")

    total_states = len(rows)
    df = (
        pl.DataFrame(rows)
        .sort("point_biserial_r", descending=True)
        .with_row_index(name="rank_index")
        .with_columns(
            pl.when(pl.col("rank_index") < n_top_positive)
            .then(pl.lit(GROUP_CYTOTOXIC_T))
            .when(pl.col("rank_index") >= (total_states - n_bottom_negative))
            .then(pl.lit(GROUP_CAF))
            .otherwise(pl.lit(GROUP_OTHER))
            .alias("cell_state_group")
        )
    )
    return Success(df)


def create_cell_state_correlation_bar_chart(
    corr_df: pl.DataFrame,
    config: CellStateCorrelationPlotConfig = CellStateCorrelationPlotConfig(),
) -> Result[alt.LayerChart, str]:
    """Build a Nature Methods vertical bar chart of sorted cell-state correlations."""
    required_cols = {"rank_index", "point_biserial_r", "cell_state_group"}
    missing = required_cols - set(corr_df.columns)
    if missing:
        return Failure(f"Missing required columns in correlation DataFrame: {sorted(missing)}")
    if corr_df.is_empty():
        return Failure("Correlation DataFrame is empty.")

    n_states = corr_df.height
    x_scale = alt.Scale(domain=[0, n_states - 1], nice=False)
    y_scale = alt.Scale(domain=[-0.45, 0.55], nice=False)

    bars = (
        alt.Chart(corr_df)
        .mark_bar(width=1.45, cornerRadius=0, opacity=0.95)
        .encode(
            x=alt.X(
                "rank_index:Q",
                title="Deconvoluted Cell State (Sorted by Correlation)",
                scale=x_scale,
                axis=alt.Axis(
                    values=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450],
                    labelFontSize=10.5,
                    titleFontSize=12,
                    titleFontWeight="bold",
                    grid=False,
                    tickColor="#1E293B",
                    domainColor="#CBD5E1",
                ),
            ),
            y=alt.Y(
                "point_biserial_r:Q",
                title="Point-Biserial Correlation (r)",
                scale=y_scale,
                axis=alt.Axis(
                    values=[-0.4, -0.2, 0.0, 0.2, 0.4],
                    labelFontSize=10.5,
                    titleFontSize=12,
                    titleFontWeight="bold",
                    gridColor="#F1F5F9",
                    tickColor="#1E293B",
                    domainColor="#CBD5E1",
                ),
            ),
            color=alt.Color(
                "cell_state_group:N",
                title="Highlighted Cell State Lineage",
                scale=alt.Scale(domain=list(GROUP_ORDER), range=list(GROUP_COLORS)),
                legend=alt.Legend(
                    orient="bottom",
                    direction="horizontal",
                    columns=3,
                    titleFontSize=11,
                    titleFontWeight="bold",
                    labelFontSize=10.5,
                    symbolType="square",
                    symbolSize=110,
                    labelLimit=260,
                ),
            ),
            tooltip=[
                alt.Tooltip("rank_index:Q", title="Sorted Rank"),
                alt.Tooltip("cell_state_id:N", title="Cell State"),
                alt.Tooltip("cell_state_group:N", title="Lineage Group"),
                alt.Tooltip("point_biserial_r:Q", format=".4f", title="Point-Biserial r"),
                alt.Tooltip("p_value:Q", format=".3e", title="p-value"),
            ],
        )
    )

    zero_rule = (
        alt.Chart(pl.DataFrame({"y": [0.0]}))
        .mark_rule(color="#0F172A", strokeWidth=1.0)
        .encode(y=alt.Y("y:Q", scale=y_scale))
    )

    annotations_df = pl.DataFrame([
        {
            "x": 20.0,
            "y": 0.46,
            "label": "Cytotoxic T Cells (Responder-Enriched)",
            "group": GROUP_CYTOTOXIC_T,
        },
        {
            "x": float(n_states - 20),
            "y": -0.34,
            "label": "Cancer-Associated Fibroblasts (Non-Responder)",
            "group": GROUP_CAF,
        },
    ])

    callout_positive = (
        alt.Chart(annotations_df.filter(pl.col("group") == GROUP_CYTOTOXIC_T))
        .mark_text(
            align="left",
            baseline="middle",
            fontSize=10.5,
            fontWeight="bold",
            color="#D55E00",
        )
        .encode(
            x=alt.X("x:Q", scale=x_scale),
            y=alt.Y("y:Q", scale=y_scale),
            text="label:N",
        )
    )

    callout_negative = (
        alt.Chart(annotations_df.filter(pl.col("group") == GROUP_CAF))
        .mark_text(
            align="right",
            baseline="middle",
            fontSize=10.5,
            fontWeight="bold",
            color="#0072B2",
        )
        .encode(
            x=alt.X("x:Q", scale=x_scale),
            y=alt.Y("y:Q", scale=y_scale),
            text="label:N",
        )
    )

    title_params = alt.TitleParams(
        text=[
            "Point-Biserial Correlations of Deconvoluted Cell States",
            "with ICB Response (Sorted Most Positive to Most Negative)",
        ],
        subtitle="Gide et al. (2019) melanoma cohort (n = 89 patients, 452 deconvoluted cell states)",
        fontSize=13,
        subtitleFontSize=10.5,
        subtitleColor="#475569",
        fontWeight="bold",
        color="#0F172A",
        anchor="start",
        offset=10,
    )

    layered = (
        (bars + zero_rule + callout_positive + callout_negative)
        .properties(
            title=title_params,
            width=config.width,
            height=config.height,
        )
        .configure_view(stroke="#CBD5E1", strokeWidth=0.75, fill="#FFFFFF")
        .configure_axis(labelFont="sans-serif", titleFont="sans-serif")
        .configure_legend(labelFont="sans-serif", titleFont="sans-serif")
    )
    return Success(layered)


def save_parquet(df: pl.DataFrame, parquet_path: Path) -> Result[pl.DataFrame, str]:
    """Save Polars DataFrame to Parquet and return the DataFrame wrapped in Result."""
    try:
        parquet_path.parent.mkdir(parents=True, exist_ok=True)
        df.write_parquet(parquet_path)
        return Success(df)
    except Exception as exc:
        return Failure(f"Failed to write Parquet to {parquet_path}: {exc}")


def export_chart_svg_and_png(chart: alt.LayerChart, output_svg: Path) -> Result[Path, str]:
    """Export Altair chart to vector SVG via vl-convert and render 300 DPI PNG via rsvg-convert or vl-convert."""
    try:
        from plotting_utils import rasterize_svg_to_png

        output_svg.parent.mkdir(parents=True, exist_ok=True)
        svg_str = vlc.vegalite_to_svg(chart.to_dict())
        output_svg.write_text(svg_str, encoding="utf-8")

        png_path = output_svg.with_suffix(".png")
        match rasterize_svg_to_png(output_svg, png_path, scale=3.125, background=None):
            case Failure(err):
                return Failure(err)
            case Success(_):
                return Success(output_svg)
    except Exception as exc:
        return Failure(f"Failed to export SVG/PNG to '{output_svg}': {exc}")


def run_pipeline(
    config: CellStateCorrelationPlotConfig,
) -> Result[tuple[pl.DataFrame, Path], str]:
    """Monadic pipeline computing correlations, saving Parquet, and exporting SVG/PNG."""
    return (
        compute_gide_cell_state_correlations(
            gide_h5ad=config.gide_h5ad,
            log_fractions_npy=config.log_fractions_npy,
            n_top_positive=config.n_top_positive,
            n_bottom_negative=config.n_bottom_negative,
        )
        .bind(lambda corr_df: save_parquet(corr_df, config.output_parquet))
        .bind(
            lambda corr_df: create_cell_state_correlation_bar_chart(
                corr_df, config=config
            )
            .bind(lambda chart: export_chart_svg_and_png(chart, config.output_svg))
            .map(lambda saved_svg: (corr_df, saved_svg))
        )
    )


def parse_args(argv: list[str]) -> CellStateCorrelationPlotConfig:
    """Parse CLI arguments into an immutable CellStateCorrelationPlotConfig."""
    parser = argparse.ArgumentParser(
        description="Plot deconvoluted cell-state point-biserial correlations with ICB response."
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=Path("output/deconv_gide_correlations/features.svg"),
        help="Destination path for exported SVG figure.",
    )
    args = parser.parse_args(argv)
    return CellStateCorrelationPlotConfig(output_svg=args.output_svg)


def main() -> int:
    """Imperative shell entrypoint."""
    config = parse_args(sys.argv[1:])
    result = run_pipeline(config)
    if isinstance(result, Failure):
        sys.stderr.write(f"[ERROR] {result.failure()}\n")
        return 1
    corr_df, svg_path = result.unwrap()
    sys.stdout.write(
        f"[OK] Computed correlations for {corr_df.height} cell states.\n"
        f"[OK] Saved Parquet: {config.output_parquet}\n"
        f"[OK] Exported vector SVG: {svg_path}\n"
        f"[OK] Rendered 300 DPI PNG: {svg_path.with_suffix('.png')}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
