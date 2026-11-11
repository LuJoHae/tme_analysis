#!/usr/bin/env python3
"""Native single-figure redraw combining Focused 14-state Consensus-Meta42 DEG markers + bulk t-stat (panel a),
exact-circle Observed Clinical RECIST Radar Immunogram (panel b), and 90-degree-rotated Ward-biclustered
TME Triad alignments + vertical Meta42 bulk ICB forest + row/column dendrograms in a single row across
Myeloid <-> Coulton TAM (panel c), T/NK <-> Chu T-cell (panel d), and Stromal/CAF <-> Shi CAF (panel e).

All panels are heading-free (displaying only bold panel letters a-e) and accompanied by an auto-generated
importable Typst caption module at ``article/captions/combined_markers_radar_triad_main_caption.typ``.
"""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Any, Final, cast
import xml.etree.ElementTree as ET

import altair as alt
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

_REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plotting_utils import (
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_VERMILION,
    export_altair_figure,
)
from scripts.coulton_tam_comparison.core import compute_ward_biclustering, strip_coulton_cluster_number
from scripts.coulton_tam_comparison.models import DendrogramBracket
from scripts.coulton_tam_comparison.plotting import (
    _brackets_to_polyline_df,
    _extract_simplified_bicluster_tables,
    format_significance_stars,
)

alt.data_transformers.disable_max_rows()

FONT_STACK: Final[str] = "system-ui, -apple-system, sans-serif"
FONT_PX: Final[float] = 12.6
LETTER_PX: Final[float] = 15.5


class CombinedMarkersConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    benchmark_dir: Path = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    coulton_dir: Path = Path("output/results/coulton_tam_comparison")
    pan_atlas_dir: Path = Path("output/results/pan_cancer_atlas_comparison")
    out_dir: Path = Path("output/manuscript_pipeline/4j_plot_combined_markers_radar_coulton")
    caption_path: Path = Path("article/captions/combined_markers_radar_triad_main_caption.typ")
    scripts_dir: Path = Path("scripts")
    row1_h: int = 248


def load_sibling(path: Path, name: str) -> Result[ModuleType, str]:
    if not path.exists():
        return Failure(f"Missing module file: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        return Failure(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return Success(module)


def panel_letter_layer(letter: str, x_px: float = -34.0, y_px: float = -24.0) -> Any:
    """Standalone bold lowercase panel letter placed at an explicit pixel offset."""
    return (
        alt.Chart(pl.DataFrame({"lbl": [letter]}).to_pandas())
        .mark_text(
            align="left",
            baseline="bottom",
            fontSize=LETTER_PX,
            fontWeight="bold",
            color=COLOR_TEXT_PRIMARY,
            clip=False,
        )
        .encode(x=alt.value(x_px), y=alt.value(y_px), text="lbl:N")
    )


def independent(chart: Any) -> Any:
    return chart.resolve_scale(
        color="independent", shape="independent", x="independent", y="independent"
    ).resolve_legend(color="independent", shape="independent")


# ---------------------------------------------------------------------------
# Row 1: (a) Combined 14-state Top-4 DEG Marker Matrix + Bulk Welch t-stat
#        (b) Exact-Circle Observed Clinical RECIST R vs. NR 12-Spoke Radar
# ---------------------------------------------------------------------------

def _shorten_state_label(lbl: str) -> str:
    """Shorten state label by stripping 'M13:' prefix and compacting 'CXCL13+PDCD1+ Tex' to 'CXCL13+ Tex'."""
    import re as _re

    cleaned = _re.sub(r"^M\d+:?\s*", "", str(lbl).strip())
    parts = cleaned.split(" ")
    if len(parts) >= 2 and "+" in parts[0]:
        first_gene = parts[0].split("+")[0] + "+"
        return f"{first_gene} {' '.join(parts[1:])}"
    return cleaned


def build_combined_deg_and_tstat_panel(
    deg_df: pl.DataFrame, letter: str, mat_w: int, bar_w: int, height: int
) -> Any:
    """(a) Combined 14 x 4 protein-coding HGNC marker specificity matrix and paired bulk Welch t-stat bars."""
    deg_short_df = deg_df.with_columns(
        pl.col("state_label").map_elements(_shorten_state_label, return_dtype=pl.Utf8).alias("state_label")
    )
    state_order_df = (
        deg_short_df.select(["cell_state", "state_label", "state_bulk_beta"])
        .unique()
        .sort("state_bulk_beta", descending=True)
    )
    ordered_state_labels = state_order_df["state_label"].to_list()
    deg_plot_df = deg_short_df.with_columns(
        (pl.lit("#") + pl.col("gene_rank").cast(pl.Utf8)).alias("rank_col")
    )
    rank_order = ["#1", "#2", "#3", "#4"]
    deg_pd = deg_plot_df.to_pandas()

    rect_a = (
        alt.Chart(deg_pd)
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.0)
        .encode(
            x=alt.X(
                "rank_col:N",
                sort=rank_order,
                title="Top DEG Marker Rank",
                axis=alt.Axis(labelAngle=0, labelFontSize=FONT_PX, titleFontSize=FONT_PX, titleFontWeight="bold"),
            ),
            y=alt.Y(
                "state_label:N",
                sort=ordered_state_labels,
                title="Selected Consensus-Meta42 State",
                axis=alt.Axis(
                    labelFontSize=FONT_PX,
                    titleFontSize=FONT_PX,
                    titleFontWeight="bold",
                    titlePadding=5,
                    labelLimit=215,
                ),
            ),
            color=alt.Color(
                "specificity_log2fc:Q",
                title="Marker Specificity (log2FC)",
                scale=alt.Scale(domain=[0.7, 3.5], scheme="purples"),
                legend=alt.Legend(
                    orient="top",
                    direction="horizontal",
                    gradientLength=115,
                    gradientThickness=8,
                    labelFontSize=FONT_PX,
                    titleFontSize=FONT_PX,
                ),
            ),
        )
    )
    text_a = (
        alt.Chart(deg_pd)
        .mark_text(fontSize=10.0, fontWeight="bold")
        .encode(
            x=alt.X("rank_col:N", sort=rank_order),
            y=alt.Y("state_label:N", sort=ordered_state_labels),
            text="gene:N",
            color=alt.condition(
                alt.datum.specificity_log2fc > 2.05,
                alt.value("#FFFFFF"),
                alt.value("#0F172A"),
            ),
        )
    )
    left_mat = (rect_a + text_a + panel_letter_layer(letter, x_px=-155.0, y_px=-58.0)).properties(
        width=mat_w,
        height=height,
    )

    state_gene_summary = (
        deg_short_df.group_by(["cell_state", "state_label", "enrichment_group", "state_bulk_beta"])
        .agg(
            pl.col("bulk_responder_tstat").mean().alias("mean_tstat"),
            pl.col("specificity_log2fc").mean().alias("mean_spec_lfc"),
        )
        .with_columns(
            pl.col("enrichment_group")
            .str.replace("-Enriched", "-|Enriched")
            .alias("enrichment_group_2l")
        )
        .sort("state_bulk_beta", descending=True)
    )
    sg_pd = state_gene_summary.to_pandas()
    grp_order = ["Responder-|Enriched (β > 0)", "Non-Responder-|Enriched (β < 0)"]
    grp_colors = [OKABE_VERMILION, OKABE_BLUE]

    zero_t = (
        alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
        .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
        .encode(x="x:Q")
    )
    bars_b = (
        alt.Chart(sg_pd)
        .mark_bar(stroke="#0F172A", strokeWidth=0.5)
        .encode(
            y=alt.Y("state_label:N", sort=ordered_state_labels, title=None, axis=alt.Axis(labels=False, ticks=False)),
            x=alt.X(
                "mean_tstat:Q",
                title="Bulk R vs. NR t-stat",
                axis=alt.Axis(
                    labelFontSize=FONT_PX,
                    titleFontSize=FONT_PX,
                    titleFontWeight="bold",
                    grid=True,
                    gridColor="#E2E8F0",
                ),
            ),
            color=alt.Color(
                "enrichment_group_2l:N",
                scale=alt.Scale(domain=grp_order, range=grp_colors),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="vertical",
                    offset=6,
                    labelExpr="split(datum.label, '|')",
                    rowPadding=14,
                    labelFontSize=FONT_PX - 1.0,
                    labelLimit=165,
                    symbolSize=65,
                ),
            ),
        )
        .properties(width=bar_w, height=height)
    )
    return independent(alt.hconcat(left_mat, bars_b + zero_t, spacing=8))


def build_observed_recist_radar_panel(
    f4m: ModuleType,
    proj_fractions_df: pl.DataFrame,
    states_df: pl.DataFrame,
    letter: str,
    side_px: int,
) -> Result[Any, str]:
    """(b) Heading-free Observed Clinical RECIST R vs. NR 12-spoke radar immunogram with 1:1 square aspect ratio."""
    compact_radar_labels: dict[str, tuple[str, str]] = {
        "M13": ("CXCL13+ Tex", "CD8+ T"),
        "M15": ("TCF7+ Stem CD8", "CD8+ T"),
        "M34": ("GZMK+ Tem CD8", "CD8+ T"),
        "M14": ("CXCL9+ M1-TAM", "Myeloid"),
        "M37": ("MS4A1+ TLS B", "B Cell"),
        "M29": ("JCHAIN+ Plasma", "Plasma"),
        "M23": ("FAP+ Mat. CAF", "Stromal"),
        "M17": ("ISLR+ myCAF", "Stromal"),
        "M19": ("S100A8+ MDSC", "Myeloid"),
        "M21": ("FOXP3+ Treg", "CD4+ Treg"),
        "M26": ("VIM+ZEB1+ EMT", "Tumor"),
        "M12": ("CA9+ Hypox.", "Tumor"),
    }
    f4m.META42_COMPACT_LABELS.update(compact_radar_labels)
    cfg = f4m.RadarPlotConfig()
    spokes_res: Result[Any, str] = f4m.select_radar_marker_spokes(states_df, cfg.preset)
    comp_res: Result[pl.DataFrame, str] = spokes_res.bind(
        lambda sp: f4m.compute_compositional_radar_features(proj_fractions_df, sp)
    )
    class_res: Result[pl.DataFrame, str] = comp_res.bind(
        lambda cdf: f4m.compute_oof_meta42_classifications(cdf, mode=cfg.classifier_mode, seed=cfg.seed)
    )
    match (spokes_res, class_res):
        case (Success(spokes), Success(classified_df)):
            obs_r_lbl = next(
                g for g in classified_df["obs_group"].unique().to_list() if str(g).startswith("Observed R ")
            )
            obs_nr_lbl = next(
                g for g in classified_df["obs_group"].unique().to_list() if str(g).startswith("Observed NR ")
            )
            obs_bounds = (-0.30, 0.30)
            obs_rings = (-0.20, -0.10, 0.00, 0.10, 0.20)
            poly_res: Result[tuple[pl.DataFrame, pl.DataFrame], str] = f4m.summarize_radar_polygons(
                compositional_df=classified_df,
                group_col="obs_group",
                group_order=(obs_r_lbl, obs_nr_lbl),
                spokes=spokes,
                radial_mode=cfg.radial_mode,
                val_min=obs_bounds[0],
                val_max=obs_bounds[1],
            )
            return poly_res.map(
                lambda pair: _flatten_radar_layers(
                    f4m.build_radar_panel_chart(
                        closed_df=pair[0],
                        spoke_stats_df=pair[1],
                        spokes=spokes,
                        group_order=(obs_r_lbl, obs_nr_lbl),
                        group_colors=(OKABE_VERMILION, OKABE_BLUE),
                        group_shapes=("circle", "square"),
                        ring_values=obs_rings,
                        val_min=obs_bounds[0],
                        val_max=obs_bounds[1],
                        unit_suffix="CLR z",
                        panel_title=letter,
                        panel_subtitle="",
                        show_sig_stars=True,
                        width_px=side_px,
                        height_px=side_px,
                    ),
                    letter=letter,
                    side_px=side_px,
                )
            )
        case _:
            return Failure("Failed to compute observed RECIST radar features.")


def _flatten_radar_layers(chart: Any, letter: str, side_px: int) -> Any:
    """Flatten nested LayerChart and enforce strict 1:1 square dimensions so the radar is an exact circle."""
    first = chart.layer[0]
    sub_layers = list(first.layer) if hasattr(first, "layer") else [first]
    return alt.LayerChart(
        layer=[*sub_layers, *chart.layer[1:], panel_letter_layer(letter, x_px=-85.0, y_px=-58.0)],
        width=side_px,
        height=side_px,
        view=chart.view,
    )


# ---------------------------------------------------------------------------
# Row 2: 90-Degree Rotated Ward-Biclustered Proximity + Vertical Meta42 Bulk β (95% CI)
#        + Column & Row Ward Dendrograms for (c), (d), and (e) in a Single Row
# ---------------------------------------------------------------------------

def _render_rotated_bicluster_panel(
    sub_clean: pl.DataFrame,
    meta_side_df: pl.DataFrame,
    y_domain: list[str],
    x_domain: list[str],
    row_brackets: Sequence[DendrogramBracket],
    col_brackets: Sequence[DendrogramBracket],
    y_col: str,
    y_axis_title: str,
    x_axis_title: str,
    letter: str,
    heat_w: int,
    heat_h: int,
    forest_h: int,
    y_domain_forest: list[float],
    y_axis_extent: int = 96,
) -> Any:
    """Render 90°-rotated panel:
    top:    Ward column dendrogram over Consensus-Meta42 states (with panel letter)
    middle: Vertical Meta42 Bulk β (95% CI) forest plot over Consensus-Meta42 states
    bottom: Biclustered proximity heatmap (Meta42 on x, External atlas on y, no badges) + right Ward row dendrogram
    """
    n_rows = len(y_domain)
    n_cols = len(x_domain)

    col_dendro_df = _brackets_to_polyline_df(col_brackets)
    col_dendro_line = (
        alt.Chart(col_dendro_df, view=alt.ViewConfig(strokeOpacity=0))
        .mark_line(color=COLOR_TEXT_PRIMARY, strokeWidth=1.05)
        .encode(
            x=alt.X(
                "u:Q",
                scale=alt.Scale(domain=[0.0, float(n_cols)], nice=False, zero=False),
                axis=None,
            ),
            y=alt.Y(
                "h:Q",
                scale=alt.Scale(domain=[0.0, 1.08], nice=False, zero=False),
                axis=None,
            ),
            detail="bracket_id:N",
            order="pt_ord:Q",
        )
    )
    col_dendro = (col_dendro_line + panel_letter_layer(letter, x_px=-34.0, y_px=-4.0)).properties(
        width=heat_w,
        height=34,
        view=alt.ViewConfig(strokeOpacity=0),
    )

    side_df = meta_side_df.with_columns(
        (pl.col("meta42_bulk_beta") - 1.96 * pl.col("meta42_bulk_se")).alias("ci_low"),
        (pl.col("meta42_bulk_beta") + 1.96 * pl.col("meta42_bulk_se")).alias("ci_high"),
        pl.when(pl.col("meta42_bulk_beta") >= 0.0)
        .then(pl.lit("Responder (β > 0)"))
        .otherwise(pl.lit("Non-Responder (β < 0)"))
        .alias("meta_dir"),
    )
    zero_h = (
        alt.Chart(pl.DataFrame({"y0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color="#64748B", strokeWidth=0.75)
        .encode(y="y0:Q")
    )
    vert_ci = (
        alt.Chart(side_df)
        .mark_rule(strokeWidth=1.5)
        .encode(
            x=alt.X(
                "meta_state_label:N",
                sort=x_domain,
                scale=alt.Scale(paddingInner=0.2, paddingOuter=0.1),
                axis=None,
            ),
            y=alt.Y(
                "ci_low:Q",
                title=["Meta42 Bulk β", "(95% CI)"],
                scale=alt.Scale(domain=y_domain_forest),
                axis=alt.Axis(
                    values=[v for v in range(-4, 5, 2) if y_domain_forest[0] + 0.3 <= v <= y_domain_forest[1] - 0.3],
                    labelFontSize=FONT_PX,
                    titleFontSize=FONT_PX - 0.8,
                    titleFontWeight="bold",
                    titlePadding=3,
                    minExtent=22,
                    maxExtent=22,
                ),
            ),
            y2="ci_high:Q",
            color=alt.Color(
                "meta_dir:N",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_VERMILION, OKABE_BLUE],
                ),
                legend=None,
            ),
        )
    )
    vert_pts = (
        alt.Chart(side_df)
        .mark_circle(size=50, opacity=1.0)
        .encode(
            x=alt.X("meta_state_label:N", sort=x_domain, axis=None),
            y="meta42_bulk_beta:Q",
            color=alt.Color(
                "meta_dir:N",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_VERMILION, OKABE_BLUE],
                ),
                legend=None,
            ),
        )
    )
    vert_stars_pos = (
        alt.Chart(side_df.filter((pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") >= 0.0)))
        .mark_text(align="center", baseline="bottom", dy=-3, fontSize=FONT_PX - 1.0, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(x=alt.X("meta_state_label:N", sort=x_domain, axis=None), y=alt.Y("ci_high:Q"), text="meta_star:N")
    )
    vert_stars_neg = (
        alt.Chart(side_df.filter((pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") < 0.0)))
        .mark_text(align="center", baseline="top", dy=4, fontSize=FONT_PX - 1.0, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(x=alt.X("meta_state_label:N", sort=x_domain, axis=None), y=alt.Y("ci_low:Q"), text="meta_star:N")
    )
    top_forest = (zero_h + vert_ci + vert_pts + vert_stars_pos + vert_stars_neg).properties(
        width=heat_w,
        height=forest_h,
    )

    heat_rects = (
        alt.Chart(sub_clean)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.5)
        .encode(
            x=alt.X(
                "meta_state_label:N",
                sort=x_domain,
                scale=alt.Scale(paddingInner=0.0, paddingOuter=0.0),
                title=None,
                axis=alt.Axis(
                    labelAngle=-40,
                    labelFontSize=FONT_PX - 0.5,
                    labelLimit=185,
                    minExtent=95,
                    maxExtent=95,
                ),
            ),
            y=alt.Y(
                f"{y_col}:N",
                sort=y_domain,
                scale=alt.Scale(paddingInner=0.0, paddingOuter=0.0),
                title=y_axis_title,
                axis=alt.Axis(
                    labelFontSize=FONT_PX - 1.0,
                    titleFontSize=FONT_PX - 0.6,
                    titleFontWeight="bold",
                    titlePadding=6,
                    minExtent=y_axis_extent,
                    maxExtent=y_axis_extent,
                    labelLimit=165,
                ),
            ),
            color=alt.Color(
                "composite_substate_score:Q",
                title="Substate Proximity",
                scale=alt.Scale(
                    domain=[-0.55, 0.0, 0.55],
                    range=["#7E22CE", "#FFFFFF", OKABE_BLUISH_GREEN],
                    clamp=True,
                ),
                legend=None,
            ),
        )
        .properties(width=heat_w, height=heat_h)
    )

    row_dendro_df = _brackets_to_polyline_df(row_brackets)
    row_dendro = (
        alt.Chart(row_dendro_df, view=alt.ViewConfig(strokeOpacity=0))
        .mark_line(color=COLOR_TEXT_PRIMARY, strokeWidth=1.05)
        .encode(
            x=alt.X(
                "h:Q",
                scale=alt.Scale(domain=[0.0, 1.08], nice=False, zero=False),
                axis=None,
            ),
            y=alt.Y(
                "u:Q",
                scale=alt.Scale(domain=[float(n_rows), 0.0], nice=False, zero=False),
                axis=None,
            ),
            detail="bracket_id:N",
            order="pt_ord:Q",
        )
        .properties(width=26, height=heat_h)
    )

    bottom_heat_row = independent(alt.hconcat(heat_rects, row_dendro, spacing=1))
    x_title_strip = (
        alt.Chart(
            pl.DataFrame({"x": [0.5 * float(heat_w)], "y": [0.0], "t": [x_axis_title]}),
            view=alt.ViewConfig(strokeOpacity=0),
        )
        .mark_text(
            align="center",
            baseline="top",
            dy=-4,
            fontSize=FONT_PX,
            fontWeight="bold",
            color=COLOR_TEXT_PRIMARY,
            clip=False,
        )
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=[0.0, float(heat_w)], nice=False, zero=False), axis=None),
            y=alt.Y("y:Q", scale=alt.Scale(domain=[0.0, 12.0], nice=False, zero=False), axis=None),
            text="t:N",
        )
        .properties(width=heat_w, height=12)
    )
    return independent(alt.vconcat(col_dendro, top_forest, bottom_heat_row, x_title_strip, spacing=8))


def build_centered_bottom_legends(total_w: int, center_x: float) -> Any:
    """Render a centered, vertically aligned horizontal legend strip for Substate Proximity and Meta42 Outcome."""
    cx = center_x
    bar_x0 = cx - 245.0
    bar_x1 = cx - 85.0
    n_steps = 48
    step_w = (bar_x1 - bar_x0) / float(n_steps)
    grad_rows = [
        {
            "x": bar_x0 + i * step_w,
            "x2": bar_x0 + (i + 1) * step_w + 0.35,
            "y": 15.0,
            "y2": 24.0,
            "val": -0.55 + 1.10 * ((i + 0.5) / float(n_steps)),
        }
        for i in range(n_steps)
    ]
    x_scale = alt.Scale(domain=[0.0, float(total_w)], nice=False, zero=False)
    y_scale = alt.Scale(domain=[36.0, 0.0], nice=False, zero=False)

    grad_rects = (
        alt.Chart(pl.DataFrame(grad_rows), view=alt.ViewConfig(strokeOpacity=0))
        .mark_rect()
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            x2="x2:Q",
            y=alt.Y("y:Q", scale=y_scale, axis=None),
            y2="y2:Q",
            color=alt.Color(
                "val:Q",
                scale=alt.Scale(
                    domain=[-0.55, 0.0, 0.55],
                    range=["#7E22CE", "#FFFFFF", OKABE_BLUISH_GREEN],
                    clamp=True,
                ),
                legend=None,
            ),
        )
    )
    grad_border = (
        alt.Chart(pl.DataFrame([{"x": bar_x0, "x2": bar_x1, "y": 15.0, "y2": 24.0}]))
        .mark_rect(fillOpacity=0.0, stroke=COLOR_HAIRLINE, strokeWidth=0.6)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            x2="x2:Q",
            y=alt.Y("y:Q", scale=y_scale, axis=None),
            y2="y2:Q",
        )
    )
    titles_layer = (
        alt.Chart(
            pl.DataFrame(
                [
                    {"x": bar_x0, "y": 1.0, "t": "Substate Proximity"},
                    {"x": cx - 20.0, "y": 1.0, "t": "Meta42 Outcome"},
                ]
            )
        )
        .mark_text(align="left", baseline="top", fontSize=FONT_PX, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.Y("y:Q", scale=y_scale, axis=None),
            text="t:N",
        )
    )
    ticks_left = (
        alt.Chart(pl.DataFrame([{"x": bar_x0, "y": 26.5, "t": "−0.5"}]))
        .mark_text(align="left", baseline="top", fontSize=FONT_PX - 0.5, color=COLOR_TEXT_SECONDARY)
        .encode(x=alt.X("x:Q", scale=x_scale, axis=None), y=alt.Y("y:Q", scale=y_scale, axis=None), text="t:N")
    )
    ticks_mid = (
        alt.Chart(pl.DataFrame([{"x": 0.5 * (bar_x0 + bar_x1), "y": 26.5, "t": "0.0"}]))
        .mark_text(align="center", baseline="top", fontSize=FONT_PX - 0.5, color=COLOR_TEXT_SECONDARY)
        .encode(x=alt.X("x:Q", scale=x_scale, axis=None), y=alt.Y("y:Q", scale=y_scale, axis=None), text="t:N")
    )
    ticks_right = (
        alt.Chart(pl.DataFrame([{"x": bar_x1, "y": 26.5, "t": "0.5"}]))
        .mark_text(align="right", baseline="top", fontSize=FONT_PX - 0.5, color=COLOR_TEXT_SECONDARY)
        .encode(x=alt.X("x:Q", scale=x_scale, axis=None), y=alt.Y("y:Q", scale=y_scale, axis=None), text="t:N")
    )
    outcome_pts = (
        alt.Chart(
            pl.DataFrame(
                [
                    {"x": cx - 14.0, "y": 20.0, "cat": "Responder (β > 0)"},
                    {"x": cx + 122.0, "y": 20.0, "cat": "Non-Responder (β < 0)"},
                ]
            )
        )
        .mark_circle(size=70, opacity=1.0)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.Y("y:Q", scale=y_scale, axis=None),
            color=alt.Color(
                "cat:N",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_VERMILION, OKABE_BLUE],
                ),
                legend=None,
            ),
        )
    )
    outcome_lbls = (
        alt.Chart(
            pl.DataFrame(
                [
                    {"x": cx - 4.0, "y": 20.5, "t": "Responder (β > 0)"},
                    {"x": cx + 132.0, "y": 20.5, "t": "Non-Responder (β < 0)"},
                ]
            )
        )
        .mark_text(align="left", baseline="middle", fontSize=FONT_PX, color=COLOR_TEXT_SECONDARY)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.Y("y:Q", scale=y_scale, axis=None),
            text="t:N",
        )
    )
    return independent(
        (grad_rects + grad_border + titles_layer + ticks_left + ticks_mid + ticks_right + outcome_pts + outcome_lbls)
        .properties(width=total_w, height=36)
    )


def build_coulton_myeloid_rotated_panel(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    letter: str,
    heat_w: int,
    heat_h: int,
    forest_h: int,
    y_axis_extent: int = 64,
) -> Result[Any, str]:
    """(c) 90°-rotated 4x4 Ward-biclustered Myeloid <-> Coulton TAM proximity + vertical Meta42 bulk beta."""
    match _extract_simplified_bicluster_tables(full_sim_df, nearest_df, coulton_ici_df):
        case Failure(err):
            return Failure(err)
        case Success(extracted):
            sub_clean_raw, meta_side_raw, _top_df, _y_dom, _x_dom, _rb, _cb = extracted
            sub_clean = sub_clean_raw.with_columns(
                pl.col("meta_state_label").map_elements(_shorten_state_label, return_dtype=pl.Utf8).alias("meta_state_label")
            )
            meta_side_filtered = meta_side_raw.with_columns(
                pl.col("meta_state_label").map_elements(_shorten_state_label, return_dtype=pl.Utf8).alias("meta_state_label")
            )
        case _:
            return Failure("Unexpected Result state in _extract_simplified_bicluster_tables")

    init_ext_rows = sorted(sub_clean["coulton_clean_label"].unique().to_list())
    init_meta_cols = sorted(sub_clean["meta_state_label"].unique().to_list())
    row_to_idx = {r: idx for idx, r in enumerate(init_ext_rows)}
    col_to_idx = {c: idx for idx, c in enumerate(init_meta_cols)}
    score_mat = np.zeros((len(init_ext_rows), len(init_meta_cols)), dtype=np.float64)
    for r_lbl, c_lbl, val in zip(
        sub_clean["coulton_clean_label"].to_list(),
        sub_clean["meta_state_label"].to_list(),
        sub_clean["composite_substate_score"].to_list(),
        strict=True,
    ):
        score_mat[row_to_idx[str(r_lbl)], col_to_idx[str(c_lbl)]] = float(val)

    match compute_ward_biclustering(score_mat, init_ext_rows, init_meta_cols):
        case Failure(b_err):
            return Failure(b_err)
        case Success(bic):
            y_domain = list(bic.ordered_rows)
            x_domain = list(bic.ordered_cols)
            row_brackets = bic.row_brackets
            col_brackets = bic.col_brackets
        case _:
            return Failure("Unexpected Result in compute_ward_biclustering for Coulton")

    return Success(
        _render_rotated_bicluster_panel(
            sub_clean=sub_clean,
            meta_side_df=meta_side_filtered,
            y_domain=y_domain,
            x_domain=x_domain,
            row_brackets=cast(Sequence[DendrogramBracket], row_brackets),
            col_brackets=cast(Sequence[DendrogramBracket], col_brackets),
            y_col="coulton_clean_label",
            y_axis_title="Coulton (2024) TAM Substates",
            x_axis_title="Significant Myeloid States (p < 0.01)",
            letter=letter,
            heat_w=heat_w,
            heat_h=heat_h,
            forest_h=forest_h,
            y_domain_forest=[-2.80, 3.30],
            y_axis_extent=y_axis_extent,
        )
    )


def build_simplified_external_rotated_panel(
    ext_sim_df: pl.DataFrame,
    ext_near_df: pl.DataFrame,
    states_df: pl.DataFrame,
    atlas_id: str,
    prefix_to_strip: str,
    letter: str,
    y_axis_title: str,
    x_axis_title: str,
    pval_cutoff: float,
    include_extra_states: tuple[str, ...],
    heat_w: int,
    heat_h: int,
    forest_h: int,
    y_domain_forest: list[float],
    y_axis_extent: int = 92,
) -> Result[Any, str]:
    """Build 90°-rotated Ward-biclustered proximity + vertical Meta42 Bulk β (95% CI) panel for Chu or Shi."""
    near_sub = ext_near_df.filter(
        (pl.col("atlas_id") == atlas_id)
        & ((pl.col("meta42_bulk_pval") < pval_cutoff) | pl.col("meta_state").is_in(list(include_extra_states)))
        & ((pl.col("rank_within_atlas") == 1) | pl.col("is_one_to_one_match"))
    )
    if near_sub.is_empty():
        return Failure(f"Zero retained pairs for simplified external panel {atlas_id}")

    retained_meta = near_sub["meta_state"].unique(maintain_order=True).to_list()
    retained_ext = near_sub["external_display"].unique(maintain_order=True).to_list()

    sim_sub = ext_sim_df.filter(
        (pl.col("atlas_id") == atlas_id)
        & pl.col("meta_state").is_in(retained_meta)
        & pl.col("external_display").is_in(retained_ext)
    )
    if sim_sub.is_empty():
        return Failure(f"Empty similarity submatrix for {atlas_id}")

    import re as _re

    clean_cols = [
        _re.sub(
            r"^S\d+_",
            "",
            str(c)[len(prefix_to_strip) :] if str(c).startswith(prefix_to_strip) else str(c),
        )
        for c in sim_sub["external_display"].to_list()
    ]
    sub_clean = sim_sub.with_columns(
        pl.Series("ext_clean_label", clean_cols),
        pl.col("meta_state_label").map_elements(_shorten_state_label, return_dtype=pl.Utf8).alias("meta_state_label"),
    )

    init_ext_rows = sorted(sub_clean["ext_clean_label"].unique().to_list())
    init_meta_cols = sorted(sub_clean["meta_state_label"].unique().to_list())
    row_to_idx = {r: idx for idx, r in enumerate(init_ext_rows)}
    col_to_idx = {c: idx for idx, c in enumerate(init_meta_cols)}
    score_mat = np.zeros((len(init_ext_rows), len(init_meta_cols)), dtype=np.float64)
    for r_lbl, c_lbl, val in zip(
        sub_clean["ext_clean_label"].to_list(),
        sub_clean["meta_state_label"].to_list(),
        sub_clean["composite_substate_score"].to_list(),
        strict=True,
    ):
        score_mat[row_to_idx[str(r_lbl)], col_to_idx[str(c_lbl)]] = float(val)

    match compute_ward_biclustering(score_mat, init_ext_rows, init_meta_cols):
        case Failure(b_err):
            return Failure(b_err)
        case Success(bic):
            y_domain = list(bic.ordered_rows)
            x_domain = list(bic.ordered_cols)
            row_brackets = bic.row_brackets
            col_brackets = bic.col_brackets
        case _:
            return Failure(f"Unexpected Result in compute_ward_biclustering for {atlas_id}")

    meta_eff = (
        states_df.filter(
            (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        )
        .select(
            pl.col("cell_state").alias("meta_state"),
            pl.col("bulk_deconv_beta").alias("meta42_bulk_beta"),
            pl.col("bulk_deconv_se").alias("meta42_bulk_se"),
            pl.col("bulk_deconv_pval").alias("meta42_bulk_pval"),
        )
    )
    meta_labels = sub_clean.select(["meta_state", "meta_state_label"]).unique()
    meta_side_base = meta_labels.join(meta_eff, on="meta_state", how="inner")
    meta_stars = [format_significance_stars(float(p)) for p in meta_side_base["meta42_bulk_pval"].to_list()]
    meta_side_df = meta_side_base.with_columns(pl.Series("meta_star", meta_stars))

    return Success(
        _render_rotated_bicluster_panel(
            sub_clean=sub_clean,
            meta_side_df=meta_side_df,
            y_domain=y_domain,
            x_domain=x_domain,
            row_brackets=row_brackets,
            col_brackets=col_brackets,
            y_col="ext_clean_label",
            y_axis_title=y_axis_title,
            x_axis_title=x_axis_title,
            letter=letter,
            heat_w=heat_w,
            heat_h=heat_h,
            forest_h=forest_h,
            y_domain_forest=y_domain_forest,
            y_axis_extent=y_axis_extent,
        )
    )


def build_combined_markers_radar_triad_typst_caption() -> str:
    """Generate a standalone importable Typst caption module for Figure 2 (panels a-e)."""
    return (
        "// Auto-generated by scripts/manuscript_pipeline/4j_plot_combined_markers_radar_coulton.py\n"
        "// Import in Typst via:\n"
        "// #import \"../captions/combined_markers_radar_triad_main_caption.typ\": caption_combined_markers_radar_triad_main\n\n"
        "#let caption_combined_markers_radar_triad_main = [\n"
        "  *Focused 14-state `Consensus-Meta42` top-4 protein-coding HGNC marker specificity and bulk responder association, observed clinical `RECIST` responder versus non-responder radar immunogram, and Ward-biclustered cross-atlas alignment across the `TME Triad` (`Coulton` `TAM`, `Chu` T-cell, and `Shi` `CAF` atlases).*\n"
        "  *(a)* Ranked protein-coding HGNC marker matrix (`Rank #1`–`Rank #4`) across the 7 top responder and 7 top non-responder `Consensus-Meta42` poles (`56` state-gene pairs; tile shading encodes one-vs-rest reference specificity $log_2 \"FC\"_(\"spec\")$ in $bold(Phi)_\"meta\"$) paired with pan-cancer bulk RNA-seq 4-gene signature responder vs. non-responder Welch $t$-statistics ($N = 1,015$ iAtlas patients).\n"
        "  *(b)* Empirical clinical `RECIST v1.1` Responder (`Observed R`, `CR/PR`, $n = 319$) versus Non-Responder (`Observed NR`, `SD/PD`, $n = 696$) 12-spoke radar immunogram across the 9 CRI-iAtlas trials for the 12 compartment-balanced `Consensus-Meta42` marker states; radial coordinates are within-trial standardized Aitchison centered log-ratios (`CLR z-score`, mean $plus.minus 95%$ CI; spoke asterisks denote Benjamini–Hochberg-adjusted two-sided Welch $t$-test false discovery rates: `*` $q < 0.05$, `**` $q < 0.01$, `***` $q < 0.001$).\n"
        "  *(c)* Simultaneously Ward-biclustered $4 times 4$ within-compartment Marker-Anchored Hybrid proximity matrix ($0.60$ whole-transcriptome `Phi_meta` specificity plus $0.40$ canonical top-4 HGNC marker enrichment $tanh(z_(\"DEG\") / 2)$) between the $4$ Coulton et al. @coultonUsingPancancerAtlas2024 true `TAM` substates achieving $q < 0.01$ in `CPI1000+` ($N = 1,446$) responder enrichment (`ECMMac`, `IFNGMac`, `IFNMac3`, `ICIMac2`; rows, with right-hand Ward dendrogram) and the $4$ `Consensus-Meta42` myeloid states achieving $p < 0.01$ in bulk deconvolution (`M19`, `M14`, `M41`, `M09`; columns), aligned with top Ward column dendrogram and `Pan-iAtlas-Meta` bulk deconvolution log-odds $hat(beta)_(\"meta\")$ ($plus.minus 95%$ CI, $N = 1,015$ patients; `*` $p < 0.05$, `**` $p < 0.01$, `***` $p < 10^(-3)$, `****` $p < 10^(-4)$).\n"
        "  *(d)* Simultaneously Ward-biclustered $8 times 8$ Marker-Anchored Hybrid proximity matrix between the $8$ matched Chu et al. pan-cancer T-cell substates @chuPancancerCellAtlas2023 (`CD4_Tfh`, `Tprolif_Tex`, `Tprolif_Treg`, `CD4_Tn`, `CD8_Tstr`, `CD8_t-Teff`, `NKT_FCGR3A`, `CD8_Teff`; rows) and the $8$ `Consensus-Meta42` `T/NK` states achieving $p < 0.01$ in bulk deconvolution (`M24`, `M13`, `M21`, `M15`, `M34`, `M05`, `M36`, `M06`; columns), paired with top Ward column dendrogram, `Pan-iAtlas-Meta` bulk deconvolution log-odds $hat(beta)_(\"meta\")$ ($plus.minus 95%$ CI), and right-hand Ward row dendrogram (full $16 times 39$ alignment in @fig-supp-atlas-chu-tcell).\n"
        "  *(e)* Simultaneously Ward-biclustered $5 times 4$ Marker-Anchored Hybrid proximity matrix between the $5$ matched Shi et al. pan-cancer `CAF` substates @shiCrosstissueMulticellularCoordination2025 (`S23_myCAF_ACTA2`, `S04_Fb_CD9`, `S07_Fb_TCF21`, `S03_Fb_LAMC1`, `S12_Fb_MMP11`; rows) and the $4$ `Consensus-Meta42` stromal/`CAF` states achieving $p < 0.01$ in bulk deconvolution (`M17`, `M23`, `M26`, `M16`; columns), paired with top Ward column dendrogram, `Pan-iAtlas-Meta` bulk deconvolution log-odds $hat(beta)_(\"meta\")$ ($plus.minus 95%$ CI), and right-hand Ward row dendrogram (full $6 times 19$ alignment in @fig-supp-atlas-shi-caf).\n"
        "]\n"
    )


def flatten_svg_nesting(svg_path: Path) -> None:
    """Unwrap passthrough <g> wrappers (no transform/clip-path/style) so Typst usvg never hits max nesting depth."""
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    root = ET.fromstring(svg_path.read_text(encoding="utf-8"))
    keep_attrs = {"transform", "clip-path", "style", "opacity", "fill", "stroke", "filter"}

    def simplify(elem: ET.Element) -> None:
        for child in list(elem):
            simplify(child)
        new_children: list[ET.Element] = []
        for child in list(elem):
            if child.tag.endswith("g") and not (set(child.attrib.keys()) & keep_attrs):
                new_children.extend(list(child))
            else:
                new_children.append(child)
        elem[:] = new_children

    simplify(root)
    svg_path.write_text(ET.tostring(root, encoding="unicode"), encoding="utf-8")


def configure_figure(chart: Any) -> Any:
    return (
        independent(chart)
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75, fill="#FFFFFF")
        .configure_axis(
            labelFont=FONT_STACK,
            titleFont=FONT_STACK,
            labelFontSize=FONT_PX,
            titleFontSize=FONT_PX,
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
            gridColor=COLOR_LIGHT_GREY,
            domainColor=COLOR_HAIRLINE,
            tickColor=COLOR_HAIRLINE,
        )
        .configure_legend(
            labelFont=FONT_STACK,
            titleFont=FONT_STACK,
            labelFontSize=FONT_PX,
            titleFontSize=FONT_PX,
            labelColor=COLOR_TEXT_SECONDARY,
            titleColor=COLOR_TEXT_PRIMARY,
        )
        .configure(font=FONT_STACK, background="#FFFFFF", padding={"left": 6, "right": 8, "top": 8, "bottom": 6})
    )


def run(config: CombinedMarkersConfig) -> Result[Path, str]:
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.caption_path.parent.mkdir(parents=True, exist_ok=True)
    f4 = load_sibling(config.scripts_dir / "manuscript_pipeline/4f_plot_meta42_responder_radar.py", "f4")
    match f4:
        case Success(f4m):
            setattr(f4m, "MIN_FONT_PX", FONT_PX)
            setattr(f4m, "TITLE_FONT_PX", LETTER_PX)
        case Failure(err):
            return Failure(err)

    paths = (
        config.benchmark_dir / "consensus_metastate_top_degs.parquet",
        config.benchmark_dir / "fractions/fractions_Consensus-Meta42_on_PanCancer-Combined.parquet",
        config.benchmark_dir / "milopy_iatlas_state_effects.parquet",
        config.coulton_dir / "meta42_vs_coulton_full_similarity.parquet",
        config.coulton_dir / "meta42_macrophage_nearest_coulton_substates.parquet",
        config.coulton_dir / "coulton_substate_ici_outcomes.parquet",
        config.pan_atlas_dir / "meta42_vs_external_atlases_similarity.parquet",
        config.pan_atlas_dir / "meta42_tcell_nearest_substates.parquet",
        config.pan_atlas_dir / "meta42_caf_nearest_substates.parquet",
    )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        return Failure("Missing inputs: " + ", ".join(missing))

    (
        deg_df,
        proj_df,
        states_df,
        full_sim_df,
        nearest_df,
        coulton_ici_df,
        ext_sim_df,
        tcell_near_df,
        caf_near_df,
    ) = (pl.read_parquet(p) for p in paths)

    panel_a = build_combined_deg_and_tstat_panel(deg_df, "a", 224, 102, config.row1_h)
    radar_res = build_observed_recist_radar_panel(f4m, proj_df, states_df, "b", config.row1_h)
    coulton_res = build_coulton_myeloid_rotated_panel(
        full_sim_df=full_sim_df,
        nearest_df=nearest_df,
        coulton_ici_df=coulton_ici_df,
        letter="c",
        heat_w=138,
        heat_h=168,
        forest_h=96,
        y_axis_extent=54,
    )
    chu_res = build_simplified_external_rotated_panel(
        ext_sim_df=ext_sim_df,
        ext_near_df=tcell_near_df,
        states_df=states_df,
        atlas_id="Chu_TCell_Atlas",
        prefix_to_strip="Chu ",
        letter="d",
        y_axis_title="Chu (2023) T-Cell Substates",
        x_axis_title="Significant T/NK States (p < 0.01)",
        pval_cutoff=0.01,
        include_extra_states=(),
        heat_w=218,
        heat_h=168,
        forest_h=96,
        y_domain_forest=[-3.20, 4.30],
        y_axis_extent=78,
    )
    shi_res = build_simplified_external_rotated_panel(
        ext_sim_df=ext_sim_df,
        ext_near_df=caf_near_df,
        states_df=states_df,
        atlas_id="Shi_PanCancer_Atlas_CAF",
        prefix_to_strip="Shi-CAF ",
        letter="e",
        y_axis_title="Shi (2025) CAF Substates",
        x_axis_title="Significant Stromal/CAF States (p < 0.01)",
        pval_cutoff=0.01,
        include_extra_states=(),
        heat_w=138,
        heat_h=168,
        forest_h=96,
        y_domain_forest=[-4.60, 1.50],
        y_axis_extent=112,
    )
    match (radar_res, coulton_res, chu_res, shi_res):
        case (Success(panel_b), Success(panel_c), Success(panel_d), Success(panel_e)):
            anchor = (
                alt.Chart(pl.DataFrame({"x": [0.0], "y": [0.0]}), view=alt.ViewConfig(strokeOpacity=0))
                .mark_point(opacity=0.0, size=0.0)
                .encode(x=alt.X("x:Q", axis=None), y=alt.Y("y:Q", axis=None))
                .properties(width=1, height=1)
            )
            row1_inner = independent(alt.hconcat(panel_a, panel_b, spacing=20))
            row2_inner = independent(alt.hconcat(panel_c, panel_d, panel_e, spacing=12))
            row1 = independent(alt.hconcat(anchor, row1_inner, spacing=0))
            row2 = independent(alt.hconcat(anchor, row2_inner, spacing=0))
            bottom_legends = build_centered_bottom_legends(total_w=940, center_x=470.0)
            figure = configure_figure(alt.vconcat(row1, row2, bottom_legends, spacing=12))
            svg_path, _ = export_altair_figure(
                figure, config.out_dir / "combined_markers_radar_coulton", max_vector_points=4000
            )
            flatten_svg_nesting(svg_path)
            config.caption_path.write_text(build_combined_markers_radar_triad_typst_caption(), encoding="utf-8")
            return Success(svg_path)
        case _:
            return Failure(
                "; ".join(
                    str(r.failure())
                    for r in (radar_res, coulton_res, chu_res, shi_res)
                    if isinstance(r, Failure)
                )
            )


def main() -> None:
    match run(CombinedMarkersConfig()):
        case Success(path):
            print(f"wrote {path}")
        case Failure(err):
            print(f"4j failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
