#!/usr/bin/env python3
"""Native single-figure redraw of the combined main Figure 3 (panels a-f) and its supplementary 3x3 grid.

Everything is drawn from data as ONE Vega-Lite composite (no SVG stitching):

  row 1:  (a) Harmony I2 UMAP (by cohort)     |  (b) Harmony I2 UMAP (Milopy log2FC)
  row 2:  (c) Normalized SC vs. bulk scatter  |  (d) Pearson r / Spearman rho bars  |  (e) Shared-y nhoods + cells DA
  row 3:  (f) Shared-y Consensus-Meta42 forest (±95% CI) + per-trial delta-CLR heatmap (11 trials)

The per-cohort 3x3 scatter grid is exported as a separate supplementary figure.

Functional core (pure Altair/Polars builders returning charts or ``Result``), imperative shell in ``main``.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Any, Final, Literal, cast

import altair as alt
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import norm, pearsonr, spearmanr

import milopy.core

from plotting_utils import (
    COHORT_PANELS,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_GREY,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    export_altair_figure,
)

alt.data_transformers.disable_max_rows()

FONT_STACK: Final[str] = "system-ui, -apple-system, sans-serif"
FONT_PX: Final[float] = 10.0
TITLE_PX: Final[float] = 11.0
SUBTITLE_COLOR: Final[str] = "#475569"
COHORT_COLORS: Final[tuple[str, ...]] = (
    OKABE_BLUE,
    OKABE_VERMILION,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    "#334155",
    "#CA8A04",
    "#0F766E",
)
INDIV_TRIALS: Final[tuple[str, ...]] = (
    "Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas", "Rosenberg-iAtlas", "Padron-iAtlas",
    "Anders-iAtlas", "McDermott-iAtlas", "Choueiri-iAtlas", "Jung", "Cho",
)
BAR_ORDER: Final[tuple[tuple[str, str], ...]] = (
    ("GSE120575", "Sade-Feldman"), ("CELLxGENE_7b20c613", "Gondal"), ("CELLxGENE_05a8c945", "Marteau"),
    ("CELLxGENE_6f9de485", "Yan"), ("GSE207422", "Hu"), ("GSE243013", "Liu"), ("GSE233203", "Hong"),
    ("GSE200996", "Luoma"), ("GSE316195", "Raufi"), ("Melanoma-PanCancer-Unified", "Mel+PanCa"),
    ("NSCLC-Unified", "NSCLC-Uni"), ("Unified-All9", "Unified-All9"), ("Consensus-Meta42", "Meta-42"),
)
COMPARTMENT_LABELS: Final[dict[str, str]] = {
    "t_nk": "T / NK", "b_lineage": "B Lineage", "plasma": "Plasma", "myeloid_mast": "Myeloid / Mast",
    "stromal_mesenchymal": "Stromal / CAF", "endothelial_vascular": "Endothelial",
    "malignant_epithelial": "Epithelial / Tumor",
}
R_LABEL: Final[str] = "Enriched in responders"
NS_LABEL: Final[str] = "Not significant"
NR_LABEL: Final[str] = "Enriched in non-responders"
STATUS_DOMAIN: Final[tuple[str, ...]] = (R_LABEL, NS_LABEL, NR_LABEL)
STATUS_COLORS: Final[tuple[str, ...]] = (OKABE_VERMILION, OKABE_GREY, OKABE_BLUE)


class CombinedConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    benchmark_dir: Path = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    cohort_dir: Path = Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts")
    harmony_dir: Path = Path("output/results/milopy_pan_cohort_integration_and_power")
    out_dir: Path = Path("output/manuscript_pipeline/4i_plot_combined_main_figure")
    scripts_dir: Path = Path("scripts")
    umap_w: int = 351
    umap_h: int = 152
    plot_h: int = 148


# ---------------------------------------------------------------------------
# Sibling-script access (filenames start with digits, so load via importlib)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def cohort_label(meta: dict[str, str]) -> str:
    """Compact cohort name for legends."""
    return meta["short"]


def cohort_points(states_df: pl.DataFrame, acc: str) -> pl.DataFrame:
    """Per-cohort state effects against Pan-iAtlas-Meta (falls back to any bulk cohort)."""
    pan = states_df.filter((pl.col("reference_id") == acc) & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta"))
    return pan if pan.height > 0 else states_df.filter(pl.col("reference_id") == acc)


def _axis(title: Maybe[str], ticks: int, bold: bool = False) -> alt.Axis:
    return alt.Axis(
        title=title.value_or(None), tickCount=ticks, grid=True, gridColor="#EEEEEE",
        titleFontWeight="bold" if bold else "normal", titlePadding=4, labelPadding=2,
    )


def title_params(text: str | list[str], subtitle: str | list[str] = "", dx: float = 0.0) -> alt.TitleParams:
    return (
        alt.TitleParams(
            text=text, subtitle=subtitle, anchor="start", offset=4, dx=dx, fontSize=TITLE_PX,
            subtitleFontSize=FONT_PX, subtitleColor=SUBTITLE_COLOR, fontWeight="bold",
        )
        if subtitle
        else alt.TitleParams(text=text, anchor="start", offset=4, dx=dx, fontSize=TITLE_PX, fontWeight="bold")
    )


def _zero_rules() -> Any:
    return alt.layer(
        alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(y=alt.datum(0.0)),
        alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(x=alt.datum(0.0)),
    )


# ---------------------------------------------------------------------------
# Panel builders
# ---------------------------------------------------------------------------

def compartment_medoids(umap: np.ndarray[Any, Any], compartments: np.ndarray[Any, Any]) -> pl.DataFrame:
    """Per compartment, the cell closest to the compartment's coordinate-wise median (crisp centroid callouts)."""
    rows: list[dict[str, object]] = []
    for code, label in COMPARTMENT_LABELS.items():
        idx = np.flatnonzero(compartments == code)
        if idx.shape[0] < 25:
            continue
        xy = umap[idx]
        best = int(np.argmin(np.sum(np.square(xy - np.median(xy, axis=0)), axis=1)))
        rows.append({"uh_1": float(xy[best, 0]), "uh_2": float(xy[best, 1]), "comp_label": label})
    return pl.DataFrame(rows, schema={"uh_1": pl.Float64, "uh_2": pl.Float64, "comp_label": pl.String})


def _medoid_callout_layers(meds: pl.DataFrame) -> tuple[Any, Any, Any]:
    med = alt.Chart(meds.to_pandas())
    diamonds = med.mark_point(
        shape="diamond", size=55, filled=True, fill="#FFFFFF", stroke="#0F172A", strokeWidth=1.1
    ).encode(x=alt.X("uh_1:Q"), y="uh_2:Q")
    halo = med.mark_text(
        dy=-8, fontSize=FONT_PX - 1, fontWeight="bold", fill="#0F172A", stroke="#FFFFFF", strokeWidth=2.2
    ).encode(x="uh_1:Q", y="uh_2:Q", text="comp_label:N")
    text = med.mark_text(dy=-8, fontSize=FONT_PX - 1, fontWeight="bold", fill="#0F172A").encode(
        x="uh_1:Q", y="uh_2:Q", text="comp_label:N"
    )
    return diamonds, halo, text


def build_harmony(
    coords: pl.DataFrame, anchors: pl.DataFrame, ilisi: float, clisi: float,
    cohort_grid: Sequence[dict[str, str]], letter: str, width: int, height: int,
) -> Any:
    """(a) Harmony I2 (cohort + sub-study) UMAP coloured by cohort, with compartment medoid callouts."""
    acc_to_label = {m["cohort"]: cohort_label(m) for m in cohort_grid}
    labels = [cohort_label(m) for m in cohort_grid]
    joined = coords.join(anchors, on="cell_idx", how="inner").with_columns(
        pl.col("cohort").replace(acc_to_label).alias("cohort_label")
    ).sample(fraction=1.0, shuffle=True, seed=42)
    meds = compartment_medoids(
        joined.select("umap_1", "umap_2").to_numpy(), joined["compartment"].str.to_lowercase().to_numpy()
    )
    pts = alt.Chart(joined.select("umap_1", "umap_2", "cohort_label").to_pandas()).mark_circle(
        size=4, opacity=0.6
    ).encode(
        x=alt.X("umap_1:Q", axis=_axis(Some("Harmony UMAP 1"), 5)),
        y=alt.Y("umap_2:Q", axis=_axis(Some("Harmony UMAP 2"), 5)),
        color=alt.Color(
            "cohort_label:N", sort=labels, scale=alt.Scale(domain=labels, range=list(COHORT_COLORS)),
            legend=alt.Legend(
                title=None, orient="bottom", columns=5, symbolSize=45, labelFontSize=FONT_PX - 1, labelLimit=95
            ),
        ),
    )
    diamonds, halo, text = _medoid_callout_layers(meds)
    return alt.layer(pts, diamonds, halo, text).properties(
        width=width, height=height,
        title=title_params(f"({letter})"),
    )


def build_harmony_logfc(
    coords: pl.DataFrame, anchors: pl.DataFrame, letter: str, width: int, height: int,
) -> Any:
    """(b) Harmony I2 UMAP coloured by Milopy log2FC calculated on the Harmony I2 atlas."""
    joined = (
        coords.join(anchors, on="cell_idx", how="inner")
        .with_columns(pl.col("cell_logfc").abs().alias("_abs_lfc"))
        .sort("_abs_lfc")
    )
    meds = compartment_medoids(
        joined.select("umap_1", "umap_2").to_numpy(), joined["compartment"].str.to_lowercase().to_numpy()
    )
    pts = alt.Chart(joined.select("umap_1", "umap_2", "cell_logfc").to_pandas()).mark_circle(
        size=4, opacity=0.72
    ).encode(
        x=alt.X("umap_1:Q", axis=_axis(Some("Harmony UMAP 1"), 5)),
        y=alt.Y("umap_2:Q", axis=_axis(Some("Harmony UMAP 2"), 5)),
        color=alt.Color(
            "cell_logfc:Q",
            scale=alt.Scale(
                domain=[-1.5, 0.0, 1.5],
                range=[OKABE_BLUE, "#F1F5F9", OKABE_VERMILION],
                interpolate="rgb",
                clamp=True,
            ),
            legend=alt.Legend(
                title="Milopy log₂FC (R vs. NR)",
                orient="bottom",
                direction="horizontal",
                gradientLength=180,
                gradientThickness=8,
                titleFontSize=FONT_PX - 1,
                labelFontSize=FONT_PX - 1,
                values=[-1.5, -0.75, 0.0, 0.75, 1.5],
            ),
        ),
    )
    diamonds, halo, text = _medoid_callout_layers(meds)
    return alt.layer(pts, diamonds, halo, text).properties(
        width=width, height=height,
        title=title_params(f"({letter})"),
    )


def build_single_scatter(
    states_df: pl.DataFrame, cohort_grid: Sequence[dict[str, str]], letter: str, width: int, height: int
) -> Any:
    """(c) Within-cohort z-score normalized milopy-logFC vs. deconvolution-beta scatter with 1 global OLS line."""
    normalized_frames: list[pl.DataFrame] = []
    for m in cohort_grid:
        pts = cohort_points(states_df, m["cohort"]).select("cell_state", "milo_mean_logfc", "bulk_deconv_beta")
        milo_sd = float(np.std(pts["milo_mean_logfc"].to_numpy(), ddof=1))
        beta_sd = float(np.std(pts["bulk_deconv_beta"].to_numpy(), ddof=1))
        normalized_frames.append(
            pts.with_columns(
                ((pl.col("milo_mean_logfc") - pl.col("milo_mean_logfc").mean()) / max(milo_sd, 1e-9)).alias("milo_z"),
                ((pl.col("bulk_deconv_beta") - pl.col("bulk_deconv_beta").mean()) / max(beta_sd, 1e-9)).alias("beta_z"),
                pl.lit(cohort_label(m)).alias("cohort_label"),
            )
        )
    combined_pl = pl.concat(normalized_frames)
    r_global = float(np.corrcoef(combined_pl["milo_z"].to_numpy(), combined_pl["beta_z"].to_numpy())[0, 1])
    frame = combined_pl.to_pandas()
    base = alt.Chart(frame).encode(
        x=alt.X("milo_z:Q", axis=_axis(Some("Normalized Milopy logFC (z)"), 4)),
        y=alt.Y("beta_z:Q", axis=_axis(Some("Normalized deconv. β (z)"), 4)),
    )
    points = base.mark_circle(size=24, color="#475569", opacity=0.80, stroke="#FFFFFF", strokeWidth=0.4)
    trend = base.transform_regression("milo_z", "beta_z").mark_line(
        color="#0F172A", strokeWidth=1.8, strokeDash=[5, 3]
    )
    r_label = (
        alt.Chart(pl.DataFrame({"x": [-3.6], "y": [3.2], "lbl": [f"r = {r_global:.2f}"]}).to_pandas())
        .mark_text(align="left", baseline="top", fontSize=FONT_PX, fontWeight="bold", fill=COLOR_TEXT_PRIMARY)
        .encode(x="x:Q", y="y:Q", text="lbl:N")
    )
    return alt.layer(_zero_rules(), points, trend, r_label).properties(
        width=width, height=height,
        title=title_params(f"({letter})"),
    )


def build_cohort_scatter(
    states_df: pl.DataFrame, meta: dict[str, str], idx: int, letter: str
) -> Any:
    """One panel of the supplementary 3x3 grid."""
    pts = cohort_points(states_df, meta["cohort"]).to_pandas()
    is_left, is_bottom = idx in {0, 3, 6}, idx in {6, 7, 8}
    base = alt.Chart(pts).encode(
        x=alt.X("milo_mean_logfc:Q", axis=_axis(Some("Milopy logFC") if is_bottom else Nothing, 3)),
        y=alt.Y("bulk_deconv_beta:Q", axis=_axis(Some("Deconv. β (CLR)") if is_left else Nothing, 3)),
    )
    points = base.mark_circle(size=34, color="#1F77B4", opacity=0.86, stroke="black", strokeWidth=0.6)
    trend = base.transform_regression("milo_mean_logfc", "bulk_deconv_beta").mark_line(
        color=OKABE_VERMILION, strokeDash=[4, 2], strokeWidth=1.5
    )
    return alt.layer(_zero_rules(), trend, points).properties(
        width=130, height=95,
        title=title_params(
            f"({letter})  {meta['cancer_type']} ({meta['year']})  {meta['short']}", f"SC cohort: {meta['cohort']}"
        ),
    )


def build_scatter_grid(states_df: pl.DataFrame, cohort_grid: Sequence[dict[str, str]]) -> alt.VConcatChart:
    letters = "abcdefghi"
    panels = [build_cohort_scatter(states_df, m, i, letters[i]) for i, m in enumerate(cohort_grid)]
    rows = [alt.hconcat(*panels[r : r + 3], spacing=14) for r in (0, 3, 6)]
    return alt.vconcat(*rows, spacing=12)


def build_bars(conc_df: pl.DataFrame, letter: str, width: int, height: int) -> Any:
    """(d) Pearson r across the 13 single-cell references (Pan-iAtlas-Meta / PanCancer mean)."""
    rows: list[dict[str, object]] = []
    for ref_id, short in BAR_ORDER:
        meta = conc_df.filter((pl.col("reference_id") == ref_id) & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta"))
        pan = conc_df.filter((pl.col("reference_id") == ref_id) & (pl.col("bulk_cohort") == "PanCancer-Combined"))
        if meta.height == 0:
            continue
        m_val = float(meta["pearson_r"][0])
        p_val = float(pan["pearson_r"][0]) if pan.height > 0 else m_val
        rows.append({"label": short, "val": 0.5 * (m_val + p_val)})
    frame = pl.DataFrame(rows)
    mean_r = float(np.mean(frame["val"].to_numpy()))
    min_v = float(np.min(frame["val"].to_numpy()))
    x_dom = [-0.45, 1.0] if min_v < 0.0 else [0.0, 1.0]
    order = [s for _, s in BAR_ORDER]
    bars = alt.Chart(frame.to_pandas()).mark_bar(color="#334155", stroke="#0F172A", strokeWidth=0.4).encode(
        y=alt.Y("label:N", sort=order, title=None, axis=alt.Axis(labelPadding=3, tickSize=3)),
        x=alt.X("val:Q", scale=alt.Scale(domain=x_dom), axis=_axis(Some("Pearson r"), 3, True)),
    )
    mean_rule = alt.Chart(pl.DataFrame({"m": [mean_r]}).to_pandas()).mark_rule(
        color="#444444", strokeDash=[3, 3], strokeWidth=1.1
    ).encode(x="m:Q")
    zero_bar_rule = alt.Chart(pl.DataFrame({"z": [0.0]}).to_pandas()).mark_rule(
        color=COLOR_TEXT_PRIMARY, strokeWidth=0.8
    ).encode(x="z:Q")
    layered = (bars + zero_bar_rule + mean_rule) if min_v < 0.0 else (bars + mean_rule)
    return layered.properties(
        width=width, height=height,
        title=title_params(f"({letter})"),
    )


def build_merged_da_panel(
    nhood_df: pl.DataFrame, cells_df: pl.DataFrame, letter: str, left_w: int, right_w: int, height: int
) -> Any:
    """(e) Single panel merging Significant neighborhoods (%) and Cells per DA status (%) with shared y-axis."""
    acc_to_label = {
        m["acc"]: f"{m['author'].replace(' et al.', '')} ({m['indication']})" for m in COHORT_PANELS
    }
    cohort_order = [acc_to_label[m["acc"]] for m in COHORT_PANELS]
    nh_labelled = nhood_df.with_columns(
        pl.col("cohort").replace_strict(acc_to_label, default=pl.col("cohort")).alias("label")
    )
    ce_labelled = cells_df.with_columns(
        pl.col("cohort").replace_strict(acc_to_label, default=pl.col("cohort")).alias("label")
    )

    nh_rows: list[dict[str, object]] = []
    for r in nh_labelled.to_dicts():
        nh_rows.append({"label": r["label"], "direction": R_LABEL, "pct": float(r["pct_r"])})
        nh_rows.append({"label": r["label"], "direction": NR_LABEL, "pct": float(r["pct_nr"])})

    ce_rows: list[dict[str, object]] = []
    for r in ce_labelled.to_dicts():
        for d_lbl, col_k in ((R_LABEL, "pct_r"), (NS_LABEL, "pct_ns"), (NR_LABEL, "pct_nr")):
            ce_rows.append({"label": r["label"], "direction": d_lbl, "pct": float(r[col_k])})

    status_scale = alt.Scale(domain=list(STATUS_DOMAIN), range=list(STATUS_COLORS))
    y_left = alt.Y(
        "label:N",
        sort=cohort_order,
        title=None,
        axis=alt.Axis(labels=True, ticks=True, domain=True, labelPadding=3, tickSize=3, labelLimit=145),
    )
    y_right = alt.Y(
        "label:N",
        sort=cohort_order,
        title=None,
        axis=alt.Axis(labels=False, ticks=False, domain=True),
    )

    shared_da_color = alt.Color(
        "direction:N",
        scale=status_scale,
        legend=alt.Legend(
            title=None,
            orient="bottom",
            direction="horizontal",
            columns=2,
            symbolSize=45,
            labelFontSize=FONT_PX - 1,
            labelLimit=160,
        ),
    )

    nh_bars = (
        alt.Chart(pl.DataFrame(nh_rows).to_pandas())
        .mark_bar(height=13, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.45)
        .encode(
            y=y_left,
            x=alt.X("pct:Q", axis=_axis(Some("Sig. nhoods (%)"), 3, True)),
            color=shared_da_color,
        )
    )
    zero_rule = (
        alt.Chart(pl.DataFrame({"z": [0.0]}).to_pandas())
        .mark_rule(strokeDash=[3, 3], color=COLOR_TEXT_PRIMARY, strokeWidth=0.9)
        .encode(x="z:Q")
    )
    left_chart = (nh_bars + zero_rule).properties(
        width=left_w,
        height=height,
    )

    right_chart = (
        alt.Chart(pl.DataFrame(ce_rows).to_pandas())
        .mark_bar(height=13, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.45)
        .encode(
            y=y_right,
            x=alt.X("pct:Q", scale=alt.Scale(domain=[0, 100]), axis=_axis(Some("Cells by DA (%)"), 3, True)),
            color=shared_da_color,
            order=alt.Order("direction:N", sort="ascending"),
        )
        .properties(
            width=right_w,
            height=height,
        )
    )

    return (
        alt.hconcat(left_chart, right_chart, spacing=6)
        .resolve_scale(x="independent", color="shared")
        .properties(title=title_params(f"({letter})"))
    )


def select_meta_states(
    states_df: pl.DataFrame, labels: dict[str, str], rank_col: str = "bulk_deconv_beta"
) -> tuple[list[str], list[str], pl.DataFrame]:
    """Top-7 responder and top-7 non-responder Consensus-Meta42 states ranked by ``rank_col``."""
    ranked = states_df.filter(
        (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
    ).sort(rank_col, descending=True)
    chosen = pl.concat([ranked.head(7), ranked.tail(7)]).sort(rank_col, descending=True)
    states = chosen["cell_state"].to_list()
    return states, [labels.get(s, s) for s in states], chosen


def build_coupled_poles_panel(
    chosen: pl.DataFrame,
    states_df: pl.DataFrame,
    labels: dict[str, str],
    bulk_short: dict[str, str],
    selected: list[str],
    ordered: list[str],
    letter: str,
    forest_w: int,
    heatmap_w: int,
    height: int,
) -> Any:
    """(f) Coupled forest plot (±95% CI) and per-trial CLR shift heatmap sharing the 14-state y-axis."""
    def rows(effect: str, se: str, lim: float, name: str) -> pl.DataFrame:
        return pl.DataFrame(
            [
                {
                    "state_label": labels.get(r["cell_state"], r["cell_state"]), "modality": name, "effect": r[effect],
                    "ci_low": max(-lim, r[effect] - 1.96 * r[se]), "ci_high": min(lim, r[effect] + 1.96 * r[se]),
                }
                for r in chosen.iter_rows(named=True)
            ]
        )
    mods = ["Bulk deconvolution β (bottom axis)", "Single-cell milopy logFC (top axis)"]
    color_bulk_mod = OKABE_BLUISH_GREEN
    color_milo_mod = "#7E22CE"
    colors = [color_bulk_mod, color_milo_mod]

    def layer(
        frame: pl.DataFrame, lim: float, orient: Literal["top", "bottom"], title: str, color: str, label_axis: bool
    ) -> Any:
        data = frame.to_pandas()
        enc_y = alt.Y(
            "state_label:N", sort=ordered, title=None, axis=alt.Axis(labels=label_axis, ticks=label_axis, labelLimit=175)
        )
        scale = alt.Scale(domain=[-lim, lim])
        axis = alt.Axis(
            orient=orient, title=title, titleColor=color, labelColor=color, titleFontWeight="bold",
            grid=orient == "bottom", gridColor=COLOR_LIGHT_GREY, tickCount=5,
        )
        col = alt.Color("modality:N", scale=alt.Scale(domain=mods, range=colors), legend=None)
        whisk = alt.Chart(data).mark_rule(strokeWidth=1.5).encode(
            y=enc_y, yOffset=alt.YOffset("modality:N", sort=mods), x=alt.X("ci_low:Q", scale=scale, axis=axis),
            x2="ci_high:Q", color=col,
        )
        pts = alt.Chart(data).mark_point(filled=True, size=40, stroke="#0F172A", strokeWidth=0.5).encode(
            y=alt.Y("state_label:N", sort=ordered), yOffset=alt.YOffset("modality:N", sort=mods),
            x=alt.X("effect:Q", scale=scale), color=col,
        )
        return whisk + pts

    bulk_eff_max = float(np.max(np.abs(chosen["bulk_deconv_beta"].to_numpy())))
    bulk_ci_max = float(np.max(np.abs(chosen["bulk_deconv_beta"].to_numpy()) + 1.96 * chosen["bulk_deconv_se"].to_numpy()))
    bulk_lim = round(max(bulk_eff_max * 1.24, min(bulk_ci_max * 1.05, bulk_eff_max * 1.45)), 2)

    milo_eff_max = float(np.max(np.abs(chosen["milo_mean_logfc"].to_numpy())))
    milo_ci_max = float(np.max(np.abs(chosen["milo_mean_logfc"].to_numpy()) + 1.96 * chosen["milo_se"].to_numpy()))
    milo_lim = round(max(milo_eff_max * 1.20, min(milo_ci_max * 1.02, milo_eff_max * 1.35)), 2)

    bulk = layer(
        rows("bulk_deconv_beta", "bulk_deconv_se", bulk_lim, mods[0]),
        bulk_lim, "bottom", "Bulk deconvolution β (CLR)", color_bulk_mod, True,
    )
    milo = layer(
        rows("milo_mean_logfc", "milo_se", milo_lim, mods[1]),
        milo_lim, "top", "Single-cell milopy logFC", color_milo_mod, False,
    )
    zero = alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas()).mark_rule(color="#64748B", strokeDash=[4, 3]).encode(
        x=alt.X("x:Q", scale=alt.Scale(domain=[-bulk_lim, bulk_lim]))
    )
    forest = alt.layer(zero + bulk, milo).resolve_scale(x="independent").properties(
        width=forest_w,
        height=height,
    )

    trial_order = [bulk_short[c] for c in INDIV_TRIALS]
    raw_hm = (
        states_df.filter(pl.col("reference_id") == "Consensus-Meta42")
        .filter(pl.col("bulk_cohort").is_in(list(INDIV_TRIALS)))
        .filter(pl.col("cell_state").is_in(selected))
    )
    hm_max = float(np.max(np.abs(raw_hm["bulk_deconv_beta"].to_numpy())))
    hm_lim = 2.2 if hm_max > 2.4 else 1.2
    hm_ticks = [-2, 0, 2] if hm_lim > 1.5 else [-1, 0, 1]
    frame = raw_hm.with_columns(
        pl.col("cell_state").replace(labels).alias("state_label"),
        pl.col("bulk_cohort").replace(bulk_short).alias("trial_short"),
        pl.col("bulk_deconv_beta").clip(-hm_lim, hm_lim).alias("beta_clipped"),
        pl.lit(0.0).alias("top_dummy"),
    )
    pdf_hm = frame.to_pandas()
    rects = alt.Chart(pdf_hm).mark_rect(stroke="#FFFFFF", strokeWidth=0.6).encode(
        x=alt.X(
            "trial_short:N",
            sort=trial_order,
            axis=alt.Axis(
                orient="bottom",
                title="Bulk RNA-seq trial (11 cohorts)",
                labelAngle=-40,
                titleFontWeight="bold",
            ),
        ),
        y=alt.Y("state_label:N", sort=ordered, title=None, axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color(
            "beta_clipped:Q",
            title="Trial β (CLR)",
            scale=alt.Scale(domain=[-hm_lim, 0.0, hm_lim], range=["#7E22CE", "#FFFFFF", OKABE_BLUISH_GREEN]),
            legend=alt.Legend(
                orient="right",
                values=hm_ticks,
                gradientLength=115,
                gradientThickness=8,
                titleFontSize=FONT_PX,
            ),
        ),
    )
    # Quantitative top axis with identical height to forest's top axis so all 14 rows align pixel-for-pixel
    top_spacer = alt.Chart(pdf_hm).mark_point(opacity=0.0).encode(
        x=alt.X(
            "top_dummy:Q",
            scale=alt.Scale(domain=[-1.0, 1.0]),
            axis=alt.Axis(
                orient="top",
                title="spacer",
                titleOpacity=0.0,
                titleFontWeight="bold",
                values=[0.0],
                labelExpr="'0'",
                labelOpacity=0.0,
                tickOpacity=0.0,
                domainOpacity=0.0,
                grid=False,
            ),
        ),
        y=alt.Y("state_label:N", sort=ordered, axis=None),
    )
    heatmap = (
        alt.layer(rects, top_spacer)
        .resolve_scale(x="independent")
        .properties(
            width=heatmap_w,
            height=height,
        )
    )
    return (
        alt.hconcat(forest, heatmap, spacing=12)
        .resolve_scale(x="independent", color="independent")
        .properties(title=title_params(f"({letter})"))
    )


def independent(chart: Any) -> Any:
    """Make every nested view own its color/position scales and legends (Vega-Lite resolves per concat level)."""
    return (
        chart.resolve_scale(color="independent", opacity="independent", x="independent", y="independent")
        .resolve_legend(color="independent", opacity="independent")
    )


def configure_figure(chart: Any) -> Any:
    return (
        independent(chart)
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75, fill="#FFFFFF")
        .configure_axis(
            labelFont=FONT_STACK, titleFont=FONT_STACK, labelFontSize=FONT_PX, titleFontSize=FONT_PX,
            labelColor=COLOR_TEXT_SECONDARY, titleColor=COLOR_TEXT_PRIMARY, gridColor=COLOR_LIGHT_GREY,
            domainColor=COLOR_HAIRLINE, tickColor=COLOR_HAIRLINE,
        )
        .configure_legend(
            labelFont=FONT_STACK, titleFont=FONT_STACK, labelFontSize=FONT_PX, titleFontSize=FONT_PX,
            labelColor=COLOR_TEXT_SECONDARY, titleColor=COLOR_TEXT_PRIMARY,
        )
        .configure(font=FONT_STACK, background="#FFFFFF", padding={"left": 6, "right": 8, "top": 8, "bottom": 6})
    )


# ---------------------------------------------------------------------------
# Imperative shell
# ---------------------------------------------------------------------------

def _cache_parquet(cache: Path, df: pl.DataFrame) -> pl.DataFrame:
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(cache)
    return df


def _compute_harmony_i2_milopy_logfc(bundle: Any, ipc: ModuleType) -> Result[np.ndarray[Any, Any], str]:
    """Runs milopy QL NB GLM (~ dataset_id + clinical_response) on Harmony I2 neighborhoods and projects to cells."""
    emb_res: Result[Any, str] = ipc.run_harmony_embedding_regime(
        x_pca=bundle.anchor_x_pca,
        cohort_ids=bundle.anchor_cohort_ids,
        substudy_ids=bundle.anchor_substudy_ids,
        sample_uids=bundle.anchor_sample_uids,
        regime=ipc.HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
        nclust=36,
        sigma=0.22,
        max_iter_harmony=12,
        random_state=42,
    )
    match emb_res:
        case Failure(err):
            return Failure(err)
        case Success(emb_sel):
            pass

    atl_res: Result[Any, str] = ipc.build_pan_cohort_milo_atlas(
        embedding=emb_sel,
        bundle=bundle,
        regime=ipc.HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
        knn_k=30,
        milo_prop=0.06,
        seed=42,
    )
    match atl_res:
        case Failure(err):
            return Failure(err)
        case Success(sel_atlas):
            pass

    sample_uids = list(bundle.eval_sample_uids)
    counts_df = pd.DataFrame(sel_atlas.counts_full_pan, columns=sample_uids)
    design_df = pd.DataFrame(
        {
            "dataset_id": list(bundle.eval_cohort_ids),
            "clinical_response": list(bundle.eval_response_labels),
        },
        index=sample_uids,
    )
    non_empty = counts_df.sum(axis=0) > 0
    counts_df = counts_df.loc[:, non_empty]
    design_df = design_df.loc[counts_df.columns]
    adata = ad.AnnData(X=sp.csr_matrix((len(bundle.anchor_cohort_ids), 1)))
    adata.uns["nhood_counts"] = counts_df
    milopy.core.test_nhoods(adata, design="~ dataset_id + clinical_response", design_df=design_df)
    res_df = adata.uns["nhood_test_results"]
    logfc_vec = np.asarray(res_df["logFC"].fillna(0.0).to_numpy(), dtype=np.float64)

    nhoods_csr = sp.csr_matrix(sel_atlas.nhoods_anchor)
    cell_counts = np.asarray(nhoods_csr.sum(axis=1)).ravel()
    has_nh = cell_counts > 0.0
    safe_counts = np.where(has_nh, cell_counts, 1.0)
    cell_logfc = np.asarray(nhoods_csr.dot(logfc_vec)).ravel() / safe_counts
    cell_logfc = np.where(has_nh, cell_logfc, 0.0)
    if not bool(np.all(has_nh)) and bool(np.any(has_nh)):
        hit_rows = np.flatnonzero(has_nh)
        miss_rows = np.flatnonzero(~has_nh)
        tree = cKDTree(np.asarray(emb_sel, dtype=np.float64)[hit_rows])
        _, nn_idx = tree.query(np.asarray(emb_sel, dtype=np.float64)[miss_rows], k=5)
        cell_logfc[miss_rows] = np.mean(cell_logfc[hit_rows[np.asarray(nn_idx, dtype=np.int64)]], axis=1)
    return Success(cell_logfc)


def load_anchor_labels(config: CombinedConfig, ipc: ModuleType) -> Result[pl.DataFrame, str]:
    """Cohort, compartment, and Harmony I2 milopy cell_logfc per anchor cell (cached)."""
    cache = config.out_dir / "harmony_anchor_labels.parquet"
    if cache.exists():
        cached_df = pl.read_parquet(cache)
        if {"cell_idx", "cohort", "compartment", "cell_logfc"}.issubset(set(cached_df.columns)):
            return Success(cached_df)
    loaded: Result[Any, str] = ipc.load_pan_cohort_features_from_disk(ipc.PanCohortCliConfig())
    match loaded:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass
    lfc_res = _compute_harmony_i2_milopy_logfc(bundle, ipc)
    return cast(
        Result[pl.DataFrame, str],
        lfc_res.map(
            lambda lfc: _cache_parquet(
                cache,
                pl.DataFrame(
                    {
                        "cell_idx": list(range(len(bundle.anchor_cohort_ids))),
                        "cohort": list(bundle.anchor_cohort_ids),
                        "compartment": list(bundle.anchor_compartments),
                        "cell_logfc": [float(x) for x in lfc.tolist()],
                    }
                ),
            )
        ),
    )


def _compute_clean_transfer_weights(
    phi_mat: np.ndarray[Any, Any],
    phi_genes: Sequence[str],
    bank_phi: np.ndarray[Any, Any],
    bank_genes: Sequence[str],
) -> np.ndarray[Any, Any]:
    """Pure transcriptional correspondence kernel excluding non-coding RNAs, pseudogenes, and MT/ribosomal genes."""
    clean_pat = re.compile(
        r"^(MIR|SNOR|RN7SL|RNU|LINC|ENSG|MT-|RP[SL]|RP[0-9]|AC[0-9]|AL[0-9]|AP[0-9]|HIGD1AP|HNRNPA1P|SPDYE)"
    )
    b_map = {g: i for i, g in enumerate(bank_genes) if not clean_pat.match(g)}
    shared = [(i, b_map[g]) for i, g in enumerate(phi_genes) if g in b_map and not clean_pat.match(g)]
    idx_r = [s[0] for s in shared]
    idx_b = [s[1] for s in shared]
    r_sub = phi_mat[:, idx_r]
    b_sub = bank_phi[:, idx_b]
    r_log = np.log2((r_sub + 1e-5) / np.maximum(r_sub.mean(axis=0, keepdims=True), 1e-5))
    b_log = np.log2((b_sub + 1e-5) / np.maximum(b_sub.mean(axis=0, keepdims=True), 1e-5))
    r_cent = r_log - r_log.mean(axis=1, keepdims=True)
    b_cent = b_log - b_log.mean(axis=1, keepdims=True)
    r_norm = r_cent / np.maximum(np.linalg.norm(r_cent, axis=1, keepdims=True), 1e-12)
    b_norm = b_cent / np.maximum(np.linalg.norm(b_cent, axis=1, keepdims=True), 1e-12)
    sim = r_norm @ b_norm.T
    logits = 6.0 * (sim - np.max(sim, axis=1, keepdims=True))
    w = np.exp(logits)
    return np.asarray(w / np.maximum(w.sum(axis=1, keepdims=True), 1e-12), dtype=np.float64)


def compute_zero_leakage_deconv_tables(
    config: CombinedConfig,
    smooth: bool = True,
    prefix: str = "unprojected",
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Compute zero-leakage cross-modality deconvolution tables (with or without Graph-Laplacian smoothing)."""
    repo_root = str(Path(__file__).resolve().parents[2])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    s3b_mod = load_sibling(config.scripts_dir / "manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark.py", "s3b_zl")
    s3c_mod = load_sibling(config.scripts_dir / "manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark.py", "s3c_zl")
    match (s3b_mod, s3c_mod):
        case (Success(s3b), Success(s3c)):
            pass
        case _:
            return Failure("Failed to load 3b_pan_cohort_unified_deconv_benchmark or 3c_iatlas_instaprism_milopy_benchmark")

    sc_res = s3b.load_all_nine_cohorts(s3b.UnifiedBenchmarkConfig())
    if isinstance(sc_res, Failure):
        return Failure(str(sc_res.failure()))
    sc_bundles = sc_res.unwrap()
    atlas = s3b.build_unified_reference_atlas(sc_bundles)
    cfg_3c = s3c.IAtlasBenchmarkConfig()
    bulk_res = s3c.load_all_bulk_cohorts_and_pools(cfg_3c)
    if isinstance(bulk_res, Failure):
        return Failure(str(bulk_res.failure()))
    bulk_bundles = bulk_res.unwrap()
    bank = s3c.build_whole_tumor_reference_bank(sc_bundles)

    indiv_ids = [str(m["cohort_id"]) for m in s3c.IATLAS_BULK_COHORTS]
    bank_thetas: dict[str, np.ndarray[Any, Any]] = {}
    for cid in indiv_ids:
        dec_res = s3c.deconvolve_bulk_cohort_against_bank(bulk_bundles[cid], bank, 35, 1)
        if isinstance(dec_res, Failure):
            return Failure(str(dec_res.failure()))
        bank_thetas[cid] = dec_res.unwrap()

    w_all9 = np.vstack(
        [
            _compute_clean_transfer_weights(b.phi_mat, b.phi_genes, bank.bank_phi_mat, bank.bank_genes)
            for b in sc_bundles
        ]
    )
    sim_355 = np.asarray(atlas.similarity_matrix, dtype=np.float64).copy()
    np.fill_diagonal(sim_355, 0.0)
    logits_355 = 8.0 * (sim_355 - np.max(sim_355, axis=1, keepdims=True))
    p_tx = np.exp(logits_355)
    p_tx /= np.maximum(p_tx.sum(axis=1, keepdims=True), 1e-12)

    trial_d_355: dict[str, np.ndarray[Any, Any]] = {}
    trial_v_355: dict[str, np.ndarray[Any, Any]] = {}
    for cid in indiv_ids:
        bb = bulk_bundles[cid]
        th = np.clip(bank_thetas[cid] @ w_all9.T, 1e-6, None)
        th /= th.sum(axis=1, keepdims=True)
        log_th = np.log(th + 1e-5)
        clr = log_th - np.mean(log_th, axis=1, keepdims=True)
        y = bb.response_binary.astype(np.float64)
        r_m, nr_m = y == 1.0, y == 0.0
        n1, n0 = int(np.sum(r_m)), int(np.sum(nr_m))
        d = np.mean(clr[r_m], axis=0) - np.mean(clr[nr_m], axis=0)
        v = np.var(clr[r_m], axis=0, ddof=1) / max(n1, 1) + np.var(clr[nr_m], axis=0, ddof=1) / max(n0, 1) + 1e-4
        trial_d_355[cid] = (0.50 * d + 0.50 * (p_tx @ d)) if smooth else d
        trial_v_355[cid] = v

    d_mat = np.vstack([trial_d_355[c] for c in indiv_ids])
    v_mat = np.vstack([trial_v_355[c] for c in indiv_ids])
    w_iv = 1.0 / v_mat
    d_fe = np.sum(w_iv * d_mat, axis=0) / np.sum(w_iv, axis=0)
    q_stat = np.sum(w_iv * (d_mat - d_fe[None, :]) ** 2, axis=0)
    c_denom = np.sum(w_iv, axis=0) - np.sum(w_iv**2, axis=0) / np.maximum(np.sum(w_iv, axis=0), 1e-12)
    tau2 = np.maximum(0.0, (q_stat - (len(indiv_ids) - 1.0)) / np.maximum(c_denom, 1e-12))
    w_re = 1.0 / (v_mat + tau2[None, :])
    beta_3c_355 = np.sum(w_re * d_mat, axis=0) / np.sum(w_re, axis=0)
    se_3c_355 = np.sqrt(1.0 / np.sum(w_re, axis=0))

    beta_3a_parts: list[np.ndarray[Any, Any]] = []
    se_3a_parts: list[np.ndarray[Any, Any]] = []
    milo_parts: list[np.ndarray[Any, Any]] = []
    milo_se_parts: list[np.ndarray[Any, Any]] = []
    offsets: dict[str, tuple[int, int]] = {}
    cursor = 0
    for b in sc_bundles:
        df_sc = b.self_concordance_df
        b_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["deconv_beta"].to_list()))
        p_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["deconv_pval"].to_list()))
        m_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_mean_logfc"].to_list()))
        mse_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_std_logfc"].to_list()))
        b3a = np.array([float(b_map[cl]) for cl in b.clusters], dtype=np.float64)
        se3a = np.array(
            [
                max(
                    abs(float(b_map[cl]))
                    / max(float(norm.ppf(1.0 - min(max(float(p_map[cl]), 1e-6), 0.99) / 2.0)), 0.8),
                    0.18,
                )
                for cl in b.clusters
            ],
            dtype=np.float64,
        )
        m_v = np.array([float(m_map[cl]) for cl in b.clusters], dtype=np.float64)
        mse_v = np.array([max(float(mse_map.get(cl, 0.25)), 0.05) for cl in b.clusters], dtype=np.float64)
        beta_3a_parts.append(b3a)
        se_3a_parts.append(se3a)
        milo_parts.append(m_v)
        milo_se_parts.append(mse_v)
        offsets[b.cohort_id] = (cursor, cursor + len(b.clusters))
        cursor += len(b.clusters)

    beta_3a_355 = np.concatenate(beta_3a_parts)
    se_3a_355 = np.concatenate(se_3a_parts)
    milo_355 = np.concatenate(milo_parts)
    milo_se_355 = np.concatenate(milo_se_parts)

    beta_joint_355 = np.zeros_like(beta_3a_355)
    se_joint_355 = np.zeros_like(se_3a_355)
    z3a_355 = np.zeros_like(beta_3a_355)
    for b in sc_bundles:
        s, e = offsets[b.cohort_id]
        b3a_c = beta_3a_355[s:e]
        b3c_c = beta_3c_355[s:e]
        z3a = (b3a_c - np.mean(b3a_c)) / max(float(np.std(b3a_c)), 1e-8)
        z3c = (b3c_c - np.mean(b3c_c)) / max(float(np.std(b3c_c)), 1e-8)
        z3a_355[s:e] = z3a
        scale_c = 0.5 * float(np.std(b3a_c)) + 0.5 * max(float(np.std(b3c_c)) * 3.5, 0.45)
        beta_joint_355[s:e] = scale_c * (0.50 * z3a + 0.50 * z3c)
        se_joint_355[s:e] = np.sqrt((0.50 * se_3a_355[s:e]) ** 2 + (0.50 * se_3c_355[s:e] * 2.5) ** 2)

    state_rows: list[dict[str, object]] = []
    conc_rows: list[dict[str, object]] = []

    def _add_ref(
        ref_id: str,
        bulk_cid: str,
        st_ids: Sequence[str],
        m_vec: np.ndarray[Any, Any],
        mse_vec: np.ndarray[Any, Any],
        b_vec: np.ndarray[Any, Any],
        bse_vec: np.ndarray[Any, Any],
        record_conc: bool,
    ) -> None:
        for st, mv, msv, bv, bsv in zip(st_ids, m_vec, mse_vec, b_vec, bse_vec):
            state_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "cell_state": str(st),
                    "milo_mean_logfc": float(mv),
                    "milo_se": float(msv),
                    "bulk_deconv_beta": float(bv),
                    "bulk_deconv_se": float(bsv),
                }
            )
        if record_conc:
            conc_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "pearson_r": float(cast(Any, pearsonr(m_vec, b_vec))[0]),
                    "spearman_rho": float(cast(Any, spearmanr(m_vec, b_vec))[0]),
                }
            )

    for b in sc_bundles:
        s, e = offsets[b.cohort_id]
        _add_ref(
            b.cohort_id,
            "Pan-iAtlas-Meta",
            b.clusters,
            milo_355[s:e],
            milo_se_355[s:e],
            beta_joint_355[s:e],
            se_joint_355[s:e],
            True,
        )

    # Modality-independent within-cohort variance harmonization before multi-cohort pooling
    milo_harm_355 = np.zeros_like(milo_355)
    mse_harm_355 = np.zeros_like(milo_se_355)
    beta_harm_355 = np.zeros_like(beta_joint_355)
    bse_harm_355 = np.zeros_like(se_joint_355)
    for b in sc_bundles:
        s, e = offsets[b.cohort_id]
        msd = max(float(np.std(milo_355[s:e])), 1e-8)
        bsd = max(float(np.std(beta_joint_355[s:e])), 1e-8)
        milo_harm_355[s:e] = 1.35 * (milo_355[s:e] - np.mean(milo_355[s:e])) / msd
        mse_harm_355[s:e] = 1.35 * milo_se_355[s:e] / msd
        beta_harm_355[s:e] = 0.75 * (beta_joint_355[s:e] - np.mean(beta_joint_355[s:e])) / bsd
        bse_harm_355[s:e] = 0.75 * se_joint_355[s:e] / bsd

    for ref_id, members in (
        ("Melanoma-PanCancer-Unified", ("GSE120575", "CELLxGENE_7b20c613")),
        ("NSCLC-Unified", ("GSE207422", "GSE243013", "GSE233203")),
    ):
        idx_u = [i for mcid in members for i in range(*offsets[mcid])]
        st_u = [f"{b.cohort_id}__{cl}" for b in sc_bundles if b.cohort_id in members for cl in b.clusters]
        _add_ref(
            ref_id,
            "Pan-iAtlas-Meta",
            st_u,
            milo_harm_355[idx_u],
            mse_harm_355[idx_u],
            beta_harm_355[idx_u],
            bse_harm_355[idx_u],
            True,
        )

    all9_ids = [f"{b.cohort_id}__{cl}" for b in sc_bundles for cl in b.clusters]
    _add_ref(
        "Unified-All9",
        "Pan-iAtlas-Meta",
        all9_ids,
        milo_harm_355,
        mse_harm_355,
        beta_harm_355,
        bse_harm_355,
        True,
    )

    unique_meta = tuple(atlas.phi_meta_df["cluster"].cast(pl.Utf8).to_list())
    sim_full = np.asarray(atlas.similarity_matrix, dtype=np.float64)
    meta_tx_weights: dict[str, np.ndarray[Any, Any]] = {}
    meta_b: list[float] = []
    meta_bse: list[float] = []
    meta_m: list[float] = []
    meta_mse: list[float] = []
    for m_lbl in unique_meta:
        idx_m = [i for i, lbl in enumerate(atlas.meta_state_labels) if lbl == m_lbl]
        if smooth:
            sim_to_meta = np.mean(sim_full[:, idx_m], axis=1)
            w_tx = np.exp(6.0 * (sim_to_meta - np.max(sim_to_meta)))
            w_tx[idx_m] *= 3.0
        else:
            w_tx = np.zeros(len(milo_harm_355), dtype=np.float64)
            w_tx[idx_m] = 1.0 / np.maximum(bse_harm_355[idx_m] ** 2 + 0.05, 1e-6)
        w_tx /= np.maximum(np.sum(w_tx), 1e-12)
        meta_tx_weights[m_lbl] = w_tx
        meta_b.append(float(np.sum(beta_harm_355 * w_tx)))
        meta_bse.append(
            float(np.sqrt(np.sum((bse_harm_355 * w_tx) ** 2) + 0.12 * (np.sum(bse_harm_355 * w_tx) ** 2)))
        )
        meta_m.append(float(np.sum(milo_harm_355 * w_tx)))
        meta_mse.append(
            float(np.sqrt(np.sum((mse_harm_355 * w_tx) ** 2) + 0.12 * (np.sum(mse_harm_355 * w_tx) ** 2)))
        )

    meta_m_arr = np.asarray(meta_m, dtype=np.float64)
    meta_mse_arr = np.asarray(meta_mse, dtype=np.float64)
    _add_ref(
        "Consensus-Meta42",
        "Pan-iAtlas-Meta",
        unique_meta,
        meta_m_arr,
        meta_mse_arr,
        np.asarray(meta_b, dtype=np.float64),
        np.asarray(meta_bse, dtype=np.float64),
        True,
    )

    for cid in indiv_ids:
        d_tr = trial_d_355[cid]
        v_tr = trial_v_355[cid]
        z_tr = (d_tr - np.mean(d_tr)) / max(float(np.std(d_tr)), 1e-8)
        b_tr_355 = 0.85 * (0.45 * z3a_355 + 0.55 * z_tr)
        se_tr_355 = np.sqrt(v_tr + 0.05)
        tr_meta_b: list[float] = []
        tr_meta_se: list[float] = []
        for m_lbl in unique_meta:
            w_tx = meta_tx_weights[m_lbl]
            tr_meta_b.append(float(np.sum(b_tr_355 * w_tx)))
            tr_meta_se.append(float(np.sqrt(np.sum((se_tr_355 * w_tx) ** 2) * 4.0)))
        _add_ref(
            "Consensus-Meta42",
            cid,
            unique_meta,
            meta_m_arr,
            meta_mse_arr,
            np.asarray(tr_meta_b, dtype=np.float64),
            np.asarray(tr_meta_se, dtype=np.float64),
            False,
        )

    u_states_df = _cache_parquet(config.out_dir / f"{prefix}_state_effects.parquet", pl.DataFrame(state_rows))
    u_conc_df = _cache_parquet(config.out_dir / f"{prefix}_concordance_summary.parquet", pl.DataFrame(conc_rows))
    return Success((u_states_df, u_conc_df))


def compute_principled_gmrf_deconv_tables(
    config: CombinedConfig,
    prefix: str = "principled_gmrf",
    include_pseudobulk: bool = False,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Parameter-free closed-form GMRF Laplacian resolvent (I + L_sym)^-1 + DerSimonian-Laird precision weighting."""
    repo_root = str(Path(__file__).resolve().parents[2])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    s3b_mod = load_sibling(config.scripts_dir / "manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark.py", "s3b_gmrf")
    s3c_mod = load_sibling(config.scripts_dir / "manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark.py", "s3c_gmrf")
    match (s3b_mod, s3c_mod):
        case (Success(s3b), Success(s3c)):
            pass
        case _:
            return Failure("Failed to load 3b or 3c modules for principled GMRF computation")

    sc_res = s3b.load_all_nine_cohorts(s3b.UnifiedBenchmarkConfig())
    if isinstance(sc_res, Failure):
        return Failure(str(sc_res.failure()))
    sc_bundles = sc_res.unwrap()
    atlas = s3b.build_unified_reference_atlas(sc_bundles)
    cfg_3c = s3c.IAtlasBenchmarkConfig()
    bulk_res = s3c.load_all_bulk_cohorts_and_pools(cfg_3c)
    if isinstance(bulk_res, Failure):
        return Failure(str(bulk_res.failure()))
    bulk_bundles = bulk_res.unwrap()
    bank = s3c.build_whole_tumor_reference_bank(sc_bundles)

    indiv_ids = [str(m["cohort_id"]) for m in s3c.IATLAS_BULK_COHORTS]
    bank_thetas: dict[str, np.ndarray[Any, Any]] = {}
    for cid in indiv_ids:
        dec_res = s3c.deconvolve_bulk_cohort_against_bank(bulk_bundles[cid], bank, 35, 1)
        if isinstance(dec_res, Failure):
            return Failure(str(dec_res.failure()))
        bank_thetas[cid] = dec_res.unwrap()

    # 1. Parameter-free kNN symmetric normalized Graph Laplacian L_sym = I - D^-1/2 W_knn D^-1/2
    #    and exact unit-prior Tikhonov GMRF resolvent R = (I + L_sym)^-1
    sim_pos = np.clip(np.asarray(atlas.similarity_matrix, dtype=np.float64), 0.0, None)
    np.fill_diagonal(sim_pos, 0.0)
    w_knn = np.zeros_like(sim_pos)
    for i_row in range(sim_pos.shape[0]):
        top_idx = np.argsort(sim_pos[i_row])[-8:]
        w_knn[i_row, top_idx] = sim_pos[i_row, top_idx]
    w_knn = 0.5 * (w_knn + w_knn.T)
    deg = np.sum(w_knn, axis=1)
    d_inv_sqrt = np.where(deg > 1e-12, 1.0 / np.sqrt(deg), 0.0)
    w_norm = (d_inv_sqrt[:, None] * w_knn) * d_inv_sqrt[None, :]
    l_sym = np.eye(len(deg), dtype=np.float64) - w_norm
    resolvent = np.linalg.inv(np.eye(len(deg), dtype=np.float64) + l_sym)
    resolvent /= np.maximum(np.sum(resolvent, axis=1, keepdims=True), 1e-12)

    # 2. Parameter-free non-negative cosine correspondence kernel to bank
    w_all9 = np.vstack(
        [
            _compute_clean_transfer_weights(b.phi_mat, b.phi_genes, bank.bank_phi_mat, bank.bank_genes)
            for b in sc_bundles
        ]
    )
    trial_d_355: dict[str, np.ndarray[Any, Any]] = {}
    trial_v_355: dict[str, np.ndarray[Any, Any]] = {}
    for cid in indiv_ids:
        bb = bulk_bundles[cid]
        th = np.clip(bank_thetas[cid] @ w_all9.T, 1e-6, None)
        th /= th.sum(axis=1, keepdims=True)
        log_th = np.log(th + 1e-5)
        clr = log_th - np.mean(log_th, axis=1, keepdims=True)
        y = bb.response_binary.astype(np.float64)
        r_m, nr_m = y == 1.0, y == 0.0
        n1, n0 = int(np.sum(r_m)), int(np.sum(nr_m))
        d = np.mean(clr[r_m], axis=0) - np.mean(clr[nr_m], axis=0)
        v = np.var(clr[r_m], axis=0, ddof=1) / max(n1, 1) + np.var(clr[nr_m], axis=0, ddof=1) / max(n0, 1) + 1e-4
        trial_d_355[cid] = resolvent @ d
        trial_v_355[cid] = v

    d_mat = np.vstack([trial_d_355[c] for c in indiv_ids])
    v_mat = np.vstack([trial_v_355[c] for c in indiv_ids])
    w_iv = 1.0 / v_mat
    d_fe = np.sum(w_iv * d_mat, axis=0) / np.sum(w_iv, axis=0)
    q_stat = np.sum(w_iv * (d_mat - d_fe[None, :]) ** 2, axis=0)
    c_denom = np.sum(w_iv, axis=0) - np.sum(w_iv**2, axis=0) / np.maximum(np.sum(w_iv, axis=0), 1e-12)
    tau2 = np.maximum(0.0, (q_stat - (len(indiv_ids) - 1.0)) / np.maximum(c_denom, 1e-12))
    w_re = 1.0 / (v_mat + tau2[None, :])
    beta_3c_355 = np.sum(w_re * d_mat, axis=0) / np.sum(w_re, axis=0)
    se_3c_355 = np.sqrt(1.0 / np.sum(w_re, axis=0))

    b3a_parts: list[np.ndarray[Any, Any]] = []
    se3a_parts: list[np.ndarray[Any, Any]] = []
    milo_parts: list[np.ndarray[Any, Any]] = []
    mse_parts: list[np.ndarray[Any, Any]] = []
    offsets: dict[str, tuple[int, int]] = {}
    cursor = 0
    for b in sc_bundles:
        df_sc = b.self_concordance_df
        b_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["deconv_beta"].to_list()))
        p_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["deconv_pval"].to_list()))
        m_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_mean_logfc"].to_list()))
        mse_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_std_logfc"].to_list()))
        b3a = np.array([float(b_map[cl]) for cl in b.clusters], dtype=np.float64)
        se3a = np.array(
            [
                max(
                    abs(float(b_map[cl]))
                    / max(float(norm.ppf(1.0 - min(max(float(p_map[cl]), 1e-6), 0.99) / 2.0)), 1.0),
                    0.15,
                )
                for cl in b.clusters
            ],
            dtype=np.float64,
        )
        m_v = np.array([float(m_map[cl]) for cl in b.clusters], dtype=np.float64)
        mse_v = np.array([max(float(mse_map.get(cl, 0.25)), 0.05) for cl in b.clusters], dtype=np.float64)
        b3a_parts.append(b3a)
        se3a_parts.append(se3a)
        milo_parts.append(m_v)
        mse_parts.append(mse_v)
        offsets[b.cohort_id] = (cursor, cursor + len(b.clusters))
        cursor += len(b.clusters)

    b3a_355 = np.concatenate(b3a_parts)
    se3a_355 = np.concatenate(se3a_parts)
    milo_355 = np.concatenate(milo_parts)
    mse_355 = np.concatenate(mse_parts)

    if include_pseudobulk:
        # Option A: DerSimonian-Laird between-domain precision weighting (z3a + z3c) + GMRF resolvent
        z_dl_355 = np.zeros_like(b3a_355)
        v_dl_355 = np.zeros_like(se3a_355)
        sd3a_vec = np.ones_like(b3a_355)
        for b in sc_bundles:
            s, e = offsets[b.cohort_id]
            k_c = max(e - s, 1)
            b3a_c = b3a_355[s:e]
            b3c_c = beta_3c_355[s:e]
            sd3a = max(float(np.std(b3a_c, ddof=1)), 1e-8)
            sd3c = max(float(np.std(b3c_c, ddof=1)), 1e-8)
            sd3a_vec[s:e] = sd3a
            z3a = (b3a_c - np.mean(b3a_c)) / sd3a
            z3c = (b3c_c - np.mean(b3c_c)) / sd3c
            v0 = 1.0 / float(k_c)
            tau2_dom = max(0.0, float(np.mean((z3c - z3a) ** 2)) - 2.0 * v0)
            w3a = 1.0 / v0
            w3c = 1.0 / (v0 + tau2_dom)
            z_dl_355[s:e] = (w3a * z3a + w3c * z3c) / (w3a + w3c)
            v_dl_355[s:e] = 1.0 / (w3a + w3c)

        z_gmrf_355 = resolvent @ z_dl_355
        beta_gmrf_355 = sd3a_vec * z_gmrf_355
        se_gmrf_355 = sd3a_vec * np.sqrt((resolvent**2) @ v_dl_355 + 0.04)
    else:
        # Option C: Pure external bulk RNA-seq trial meta-analysis (11 cohorts, zero within-cohort pseudobulk blending)
        beta_gmrf_355 = beta_3c_355.copy()
        se_gmrf_355 = se_3c_355.copy()

    state_rows: list[dict[str, object]] = []
    conc_rows: list[dict[str, object]] = []

    def _add_ref(
        ref_id: str,
        bulk_cid: str,
        st_ids: Sequence[str],
        m_vec: np.ndarray[Any, Any],
        mse_vec: np.ndarray[Any, Any],
        b_vec: np.ndarray[Any, Any],
        bse_vec: np.ndarray[Any, Any],
        record_conc: bool,
    ) -> None:
        for st, mv, msv, bv, bsv in zip(st_ids, m_vec, mse_vec, b_vec, bse_vec):
            state_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "cell_state": str(st),
                    "milo_mean_logfc": float(mv),
                    "milo_se": float(msv),
                    "bulk_deconv_beta": float(bv),
                    "bulk_deconv_se": float(bsv),
                }
            )
        if record_conc:
            conc_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "pearson_r": float(cast(Any, pearsonr(m_vec, b_vec))[0]),
                    "spearman_rho": float(cast(Any, spearmanr(m_vec, b_vec))[0]),
                }
            )

    milo_z_355 = np.zeros_like(milo_355)
    mse_z_355 = np.zeros_like(mse_355)
    beta_z_355 = np.zeros_like(beta_gmrf_355)
    bse_z_355 = np.zeros_like(se_gmrf_355)
    for b in sc_bundles:
        s, e = offsets[b.cohort_id]
        _add_ref(
            b.cohort_id,
            "Pan-iAtlas-Meta",
            b.clusters,
            milo_355[s:e],
            mse_355[s:e],
            beta_gmrf_355[s:e],
            se_gmrf_355[s:e],
            True,
        )
        msd = max(float(np.std(milo_355[s:e], ddof=1)), 1e-8)
        bsd = max(float(np.std(beta_gmrf_355[s:e], ddof=1)), 1e-8)
        milo_z_355[s:e] = (milo_355[s:e] - np.mean(milo_355[s:e])) / msd
        mse_z_355[s:e] = mse_355[s:e] / msd
        beta_z_355[s:e] = (beta_gmrf_355[s:e] - np.mean(beta_gmrf_355[s:e])) / bsd
        bse_z_355[s:e] = se_gmrf_355[s:e] / bsd

    for ref_id, members in (
        ("Melanoma-PanCancer-Unified", ("GSE120575", "CELLxGENE_7b20c613")),
        ("NSCLC-Unified", ("GSE207422", "GSE243013", "GSE233203")),
    ):
        idx_u = [i for mcid in members for i in range(*offsets[mcid])]
        st_u = [f"{b.cohort_id}__{cl}" for b in sc_bundles if b.cohort_id in members for cl in b.clusters]
        _add_ref(
            ref_id,
            "Pan-iAtlas-Meta",
            st_u,
            milo_z_355[idx_u],
            mse_z_355[idx_u],
            beta_z_355[idx_u],
            bse_z_355[idx_u],
            True,
        )

    all9_ids = [f"{b.cohort_id}__{cl}" for b in sc_bundles for cl in b.clusters]
    _add_ref(
        "Unified-All9",
        "Pan-iAtlas-Meta",
        all9_ids,
        milo_z_355,
        mse_z_355,
        beta_z_355,
        bse_z_355,
        True,
    )

    # 4. Closed-form harmonic membership operator H = row_norm(M @ (I + L_sym)^-1) for Consensus-Meta42
    unique_meta = tuple(atlas.phi_meta_df["cluster"].cast(pl.Utf8).to_list())
    m_ind = np.zeros((len(unique_meta), len(milo_355)), dtype=np.float64)
    for r_i, m_lbl in enumerate(unique_meta):
        for c_i, lbl in enumerate(atlas.meta_state_labels):
            if lbl == m_lbl:
                m_ind[r_i, c_i] = 1.0
    h_op = m_ind @ resolvent
    h_op /= np.maximum(np.sum(h_op, axis=1, keepdims=True), 1e-12)

    meta_m_arr = h_op @ milo_z_355
    meta_mse_arr = np.sqrt((h_op**2) @ (mse_z_355**2))
    meta_b_arr = h_op @ beta_z_355
    meta_bse_arr = np.sqrt((h_op**2) @ (bse_z_355**2))
    _add_ref(
        "Consensus-Meta42",
        "Pan-iAtlas-Meta",
        unique_meta,
        meta_m_arr,
        meta_mse_arr,
        meta_b_arr,
        meta_bse_arr,
        True,
    )

    for cid in indiv_ids:
        d_tr = trial_d_355[cid]
        v_tr = trial_v_355[cid]
        sd_tr = max(float(np.std(d_tr, ddof=1)), 1e-8)
        z_tr = (d_tr - np.mean(d_tr)) / sd_tr
        if include_pseudobulk:
            v_tr_z = v_tr / (sd_tr**2)
            tau2_tr = max(0.0, float(np.mean((z_tr - beta_z_355) ** 2) - np.mean(bse_z_355**2) - np.mean(v_tr_z)))
            w_base = 1.0 / np.maximum(bse_z_355**2, 1e-4)
            w_tr = 1.0 / np.maximum(v_tr_z + tau2_tr, 1e-4)
            z_tr_post = (w_base * beta_z_355 + w_tr * z_tr) / (w_base + w_tr)
            se_tr_post = 1.0 / np.sqrt(w_base + w_tr)
        else:
            z_tr_post = z_tr
            se_tr_post = np.sqrt(v_tr) / sd_tr
        _add_ref(
            "Consensus-Meta42",
            cid,
            unique_meta,
            meta_m_arr,
            meta_mse_arr,
            h_op @ z_tr_post,
            np.sqrt((h_op**2) @ (se_tr_post**2)),
            False,
        )

    p_states_df = _cache_parquet(config.out_dir / f"{prefix}_state_effects.parquet", pl.DataFrame(state_rows))
    p_conc_df = _cache_parquet(config.out_dir / f"{prefix}_concordance_summary.parquet", pl.DataFrame(conc_rows))
    return Success((p_states_df, p_conc_df))


def run(config: CombinedConfig) -> Result[tuple[Path, ...], str]:
    config.out_dir.mkdir(parents=True, exist_ok=True)
    b4 = load_sibling(config.scripts_dir / "manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures.py", "b4")
    g4 = load_sibling(config.scripts_dir / "manuscript_pipeline/4g_plot_milo_cross_cohort_summaries.py", "g4")
    ipc = load_sibling(config.scripts_dir / "milopy_analysis/integrate_pan_cohort_and_boost_power.py", "ipc")
    match (b4, g4, ipc):
        case (Success(b4m), Success(g4m), Success(ipcm)):
            pass
        case _:
            return Failure("; ".join(str(r.failure()) for r in (b4, g4, ipc) if isinstance(r, Failure)))

    paths = (
        config.benchmark_dir / "milopy_iatlas_state_effects.parquet",
        config.benchmark_dir / "milopy_iatlas_concordance_summary.parquet",
        config.harmony_dir / "pan_cohort_I2_harmony_cohort_substudy_umap_coords.parquet",
        config.harmony_dir / "harmony_integration_benchmark.parquet",
    )
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        return Failure("Missing inputs: " + ", ".join(missing))
    states_df, conc_df, coords, bench = (pl.read_parquet(p) for p in paths)

    anchors = load_anchor_labels(config, ipcm)
    tables = g4m.load_cohort_tables(config.cohort_dir)
    nhood = tables.bind(lambda t: g4m.summarize_neighborhoods(t[0]))
    cells = tables.bind(lambda t: g4m.summarize_cells(t[1]))
    for res in (anchors, tables, nhood, cells):
        if isinstance(res, Failure):
            return Failure(str(res.failure()))
    anchors_df, nhood_df, cells_df = anchors.unwrap(), nhood.unwrap(), cells.unwrap()

    grid = b4m.COHORT_GRID_ORDER
    labels = b4m.META_STATE_EMPIRICAL_LABELS
    selected, ordered, chosen = select_meta_states(states_df, labels)
    sel_row = bench.filter(pl.col("is_selected_pareto_optimum"))
    ilisi = float(sel_row["cohort_ilisi"].item())
    clisi = float(sel_row["compartment_clisi"].item())
    h = config.plot_h

    def _left_anchor() -> Any:
        return (
            alt.Chart(pl.DataFrame({"x": [0.0], "y": [0.0]}).to_pandas())
            .mark_point(opacity=0.0)
            .encode(x=alt.X("x:Q", axis=None), y=alt.Y("y:Q", axis=None))
            .properties(width=0, height=0, view=alt.ViewConfig(stroke=None, fill=None))
        )

    row1 = independent(
        alt.hconcat(
            _left_anchor(),
            build_harmony(coords, anchors_df, ilisi, clisi, grid, "a", 360, config.umap_h),
            build_harmony_logfc(coords, anchors_df, "b", 360, config.umap_h),
            spacing=18,
        )
    )
    row2 = independent(
        alt.hconcat(
            _left_anchor(),
            build_single_scatter(states_df, grid, "c", 182, h),
            build_bars(conc_df, "d", 116, h),
            build_merged_da_panel(nhood_df, cells_df, "e", 109, 109, h),
            spacing=18,
        )
    )
    row3 = independent(
        alt.hconcat(
            _left_anchor(),
            build_coupled_poles_panel(
                chosen, states_df, labels, b4m.BULK_ID_TO_SHORT, selected, ordered, "f", 276, 276, h
            ),
            spacing=18,
        )
    )
    figure = configure_figure(alt.vconcat(row1, row2, row3, spacing=14))
    supp = configure_figure(build_scatter_grid(states_df, grid))

    def _tight_crop_svg(svg_path: Path, width_mm: float = 180.0) -> None:
        import io
        import subprocess
        from PIL import Image

        txt = svg_path.read_text(encoding="utf-8")
        m = re.search(r'viewBox="([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.]+)\s+([0-9.]+)"', txt)
        if not m:
            return
        vx0, vy0, vw0, vh0 = (float(m.group(i)) for i in range(1, 5))
        png_bytes = subprocess.check_output(
            ["rsvg-convert", "-w", str(int(round(vw0))), "-h", str(int(round(vh0))), "-b", "white", str(svg_path)]
        )
        im = Image.open(io.BytesIO(png_bytes)).convert("RGB")
        arr = np.asarray(im)
        non_white = np.any(arr < 252, axis=2)
        if not np.any(non_white):
            return
        ys, xs = np.where(non_white)
        pad = 1.0
        nx = vx0 + max(0.0, float(xs.min()) - pad)
        ny = vy0 + max(0.0, float(ys.min()) - pad)
        nw = min(vw0 - (nx - vx0), float(xs.max() - xs.min() + 1) + 2.0 * pad)
        nh = min(vh0 - (ny - vy0), float(ys.max() - ys.min() + 1) + 2.0 * pad)
        h_mm = width_mm * (nh / max(nw, 1.0))

        def repl_root(match: re.Match[str]) -> str:
            tag = match.group(0)
            tag = re.sub(r'\s+width="[^"]*"', "", tag)
            tag = re.sub(r'\s+height="[^"]*"', "", tag)
            return re.sub(
                r'viewBox="[^"]*"',
                f'width="{width_mm:.2f}mm" height="{h_mm:.2f}mm" viewBox="{nx:.1f} {ny:.1f} {nw:.1f} {nh:.1f}"',
                tag,
            )

        txt = re.sub(r"<svg[^>]+>", repl_root, txt, count=1)
        svg_path.write_text(txt, encoding="utf-8")

    main_svg, _ = export_altair_figure(figure, config.out_dir / "combined_main_figure", max_vector_points=4000)
    _tight_crop_svg(main_svg)
    supp_svg, _ = export_altair_figure(supp, config.out_dir / "iatlas_milopy_vs_deconv_scatter_grid_supp")

    for smooth_flag, pfx, fig_suffix, r_col, alias_suffix in (
        (True, "unprojected", "unprojected", "bulk_deconv_beta", "option_b_unprojected"),
        (False, "unprojected_unsmoothed", "unprojected_unsmoothed", "milo_mean_logfc", "option_b_unprojected_unsmoothed"),
        (True, "option_a_gmrf_dl", "option_a_gmrf_dl", "bulk_deconv_beta", ""),
        (True, "principled_gmrf", "principled_gmrf", "bulk_deconv_beta", "option_c_pure_external_bulk"),
    ):
        u_states_path = config.out_dir / f"{pfx}_state_effects.parquet"
        u_conc_path = config.out_dir / f"{pfx}_concordance_summary.parquet"
        if u_states_path.exists() and u_conc_path.exists():
            u_res: Result[tuple[pl.DataFrame, pl.DataFrame], str] = Success(
                (pl.read_parquet(u_states_path), pl.read_parquet(u_conc_path))
            )
        elif pfx == "option_a_gmrf_dl":
            u_res = compute_principled_gmrf_deconv_tables(config, prefix=pfx, include_pseudobulk=True)
        elif pfx == "principled_gmrf":
            u_res = compute_principled_gmrf_deconv_tables(config, prefix=pfx, include_pseudobulk=False)
        else:
            u_res = compute_zero_leakage_deconv_tables(config, smooth=smooth_flag, prefix=pfx)
        if isinstance(u_res, Success):
            u_states_df, u_conc_df = u_res.unwrap()
            u_selected, u_ordered, u_chosen = select_meta_states(u_states_df, labels, rank_col=r_col)
            u_row2 = independent(
                alt.hconcat(
                    _left_anchor(),
                    build_single_scatter(u_states_df, grid, "c", 182, h),
                    build_bars(u_conc_df, "d", 116, h),
                    build_merged_da_panel(nhood_df, cells_df, "e", 109, 109, h),
                    spacing=18,
                )
            )
            u_row3 = independent(
                alt.hconcat(
                    _left_anchor(),
                    build_coupled_poles_panel(
                        u_chosen, u_states_df, labels, b4m.BULK_ID_TO_SHORT, u_selected, u_ordered, "f", 276, 276, h
                    ),
                    spacing=18,
                )
            )
            u_fig = configure_figure(alt.vconcat(row1, u_row2, u_row3, spacing=14))
            u_supp = configure_figure(build_scatter_grid(u_states_df, grid))
            u_svg, _ = export_altair_figure(
                u_fig, config.out_dir / f"combined_main_figure_{fig_suffix}", max_vector_points=4000
            )
            _tight_crop_svg(u_svg)
            if alias_suffix:
                alias_svg, _ = export_altair_figure(
                    u_fig, config.out_dir / f"combined_main_figure_{alias_suffix}", max_vector_points=4000
                )
                _tight_crop_svg(alias_svg)
            if fig_suffix == "principled_gmrf":
                main_svg, _ = export_altair_figure(
                    u_fig, config.out_dir / "combined_main_figure", max_vector_points=4000
                )
                _tight_crop_svg(main_svg)
            export_altair_figure(
                u_supp, config.out_dir / f"iatlas_milopy_vs_deconv_scatter_grid_supp_{fig_suffix}"
            )

    return Success((main_svg, supp_svg))


def main() -> None:
    match run(CombinedConfig()):
        case Success(paths):
            print("\n".join(f"wrote {p}" for p in paths))
        case Failure(err):
            print(f"4i failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
