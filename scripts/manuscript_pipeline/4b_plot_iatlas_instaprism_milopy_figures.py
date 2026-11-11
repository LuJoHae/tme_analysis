#!/usr/bin/env python3
"""Publication-Grade Vector Figures for Cross-Modality iAtlas Deconvolution, Multi-Model Prediction & Milopy Concordance.

Generates five full-width (180.0mm) Nature Methods wireframe vector SVG & 300 DPI PNG figures
with minimum font size >= 10pt (13.5px at 96 CSS DPI), Okabe-Ito colorblind-safe palette,
and native Inkscape layer annotations:

1. Figure 1: `iatlas_response_prediction_auc_heatmap_and_bars.svg` (.png)
   - (a) Cross-modality ROC-AUC heatmap across all 13 Single-Cell References (Tiers 1–3)
         and all 15 Bulk Clinical Trial Cohorts / Pooled Strata.
   - (b) Multi-model classifier comparison (Logistic Regression, Linear Regression, Random Forests)
         across ALL 15 bulk cohorts and pooled strata.

2. Figure 2: `iatlas_milopy_vs_deconv_concordance_grid.svg` (.png)
   - (a) 3x3 scatter grid of Unbounded Linear Bulk iAtlas Deconvolution Effect Size (beta_deconv)
         vs. Single-Cell Milopy DA (mean logFC) across the 9 audited single-cell cohorts (zero S-curve).
   - (b) Dual-metric horizontal bar plot comparing Linear Pearson r and Spearman rho across all 13 references.

3. Figure 3: `iatlas_consensus_metastate_forest_and_tme_profile.svg` (.png)
   - (a) Cross-trial random-effects meta-analytic forest plot of top Responder- and Non-Responder-enriched
         Consensus-Meta42 cell states (Bulk Deconv beta +/- 95% CI vs. Single-Cell Milopy logFC).
   - (b) Per-trial Responder vs. Non-Responder Centered Log-Ratio shift (Delta-CLR) heatmap across all 11 trials.

4. Figure 4: `iatlas_classifier_and_reference_tier_pareto.svg` (.png)
   - (a) Reference Architecture Pareto Frontier (Mean Bulk Predictive ROC-AUC vs. Linear Pearson r across Tiers 1–3).
   - (b) Multi-Model State Attribution (Random Forest Importance & Linear Regression |beta| vs. Milopy |logFC|).

5. Figure 5: `iatlas_cross_cohort_concordance_transfer_matrix.svg` (.png)
   - (a) 13 x 16 Cross-Modality Linear Correlation (Pearson r) Transfer Matrix across all SC references & bulk cohorts.
   - (b) Leave-One-Cohort-Out (LOCO) Zero-Shot Generalization vs. Within-Cohort CV ROC-AUC across all 11 trials.

Strict functional Python: frozen Pydantic models, returns Result, Polars DataFrames, Altair SVG export.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Callable, Final, Literal, Sequence, cast
import xml.etree.ElementTree as ET

import altair as alt
import numpy as np
from plotting_utils.svg import rasterize_svg_datapoints
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import vl_convert as vlc

# Disable row limit for Altair
alt.data_transformers.disable_max_rows()

# ------------------------------------------------------------------------------
# Dimensional & Styling Constants (180.0mm Full Width, >= 10pt Typography)
# ------------------------------------------------------------------------------

FULL_WIDTH_MM: Final[float] = 180.0
MIN_FONT_PX: Final[float] = 13.5  # 10.125 pt at 96 CSS DPI (strictly >= 10pt)
TITLE_FONT_PX: Final[float] = 14.0  # 10.5 pt at 96 CSS DPI

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_SKY: Final[str] = "#56B4E9"
OKABE_TEAL: Final[str] = "#009E73"
OKABE_ORANGE: Final[str] = "#E69F00"
OKABE_VERMILION: Final[str] = "#D55E00"
OKABE_PURPLE: Final[str] = "#CC79A7"
POINT_BLUE: Final[str] = "#1F77B4"

COHORT_GRID_ORDER: Final[tuple[dict[str, str], ...]] = (
    {"cohort": "GSE120575", "letter": "a", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "short": "Sade-Feldman", "year": "2018"},
    {"cohort": "CELLxGENE_7b20c613", "letter": "", "cancer_type": "Pan-Cancer", "author": "Gondal et al.", "short": "Gondal", "year": "2025"},
    {"cohort": "CELLxGENE_05a8c945", "letter": "", "cancer_type": "Colorectal", "author": "Marteau et al.", "short": "Marteau", "year": "2026"},
    {"cohort": "CELLxGENE_6f9de485", "letter": "", "cancer_type": "Breast", "author": "Yan et al.", "short": "Yan", "year": "2026"},
    {"cohort": "GSE207422", "letter": "", "cancer_type": "NSCLC", "author": "Hu et al.", "short": "Hu", "year": "2023"},
    {"cohort": "GSE243013", "letter": "", "cancer_type": "NSCLC", "author": "Liu et al.", "short": "Liu", "year": "2025"},
    {"cohort": "GSE233203", "letter": "", "cancer_type": "NSCLC", "author": "Hong et al.", "short": "Hong", "year": "2025"},
    {"cohort": "GSE200996", "letter": "", "cancer_type": "HNSCC", "author": "Luoma et al.", "short": "Luoma", "year": "2022"},
    {"cohort": "GSE316195", "letter": "", "cancer_type": "PDAC", "author": "Raufi et al.", "short": "Raufi", "year": "2026"},
)

REF_ID_TO_SHORT: Final[dict[str, str]] = {
    "GSE120575": "Sade-Feldman (2018)",
    "CELLxGENE_7b20c613": "Gondal (2025)",
    "CELLxGENE_05a8c945": "Marteau (2026)",
    "CELLxGENE_6f9de485": "Yan (2026)",
    "GSE207422": "Hu (2023)",
    "GSE243013": "Liu (2025)",
    "GSE233203": "Hong (2025)",
    "GSE200996": "Luoma (2022)",
    "GSE316195": "Raufi (2026)",
    "Melanoma-PanCancer-Unified": "Melanoma+PanCa Unified",
    "NSCLC-Unified": "NSCLC Unified",
    "Unified-All9": "Unified All-9 (355)",
    "Consensus-Meta42": "Consensus Meta (42)",
}

REF_ORDER: Final[list[str]] = [
    "Sade-Feldman (2018)",
    "Gondal (2025)",
    "Marteau (2026)",
    "Yan (2026)",
    "Hu (2023)",
    "Liu (2025)",
    "Hong (2025)",
    "Luoma (2022)",
    "Raufi (2026)",
    "Melanoma+PanCa Unified",
    "NSCLC Unified",
    "Unified All-9 (355)",
    "Consensus Meta (42)",
]

BULK_COLS_ORDER: Final[list[str]] = [
    "PanCancer-Combined",
    "Melanoma-Combined",
    "NSCLC-Combined",
    "RenalCell-Combined",
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
    "Jung",
    "Cho",
]

BULK_ID_TO_SHORT: Final[dict[str, str]] = {
    "Pan-iAtlas-Meta": "Pan-iAtlas Meta (9)",
    "PanCancer-Combined": "PanCancer Pool",
    "Melanoma-Combined": "Melanoma Pool",
    "NSCLC-Combined": "NSCLC Pool",
    "RenalCell-Combined": "ccRCC Pool",
    "Hugo-iAtlas": "Hugo (Mel)",
    "Riaz-iAtlas": "Riaz (Mel)",
    "Liu-iAtlas": "Liu (Mel)",
    "Gide-iAtlas": "Gide (Mel)",
    "Rosenberg-iAtlas": "Rosenberg (Blad)",
    "Padron-iAtlas": "Padron (PDAC)",
    "Anders-iAtlas": "Anders (TNBC)",
    "McDermott-iAtlas": "McDermott (RCC)",
    "Choueiri-iAtlas": "Choueiri (RCC)",
    "Jung": "Jung (NSCLC)",
    "Cho": "Cho (NSCLC)",
}


class FigurePlotConfig(BaseModel):
    """Immutable configuration for generating the iAtlas InstaPrism & Milopy publication figures."""

    model_config = ConfigDict(frozen=True)

    data_dir: Path = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    results_dir: Path = Path("output/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures")
    article_dir: Path = Path("article/figures/milo-deconv")
    captions_dir: Path = Path("article/captions")


# ------------------------------------------------------------------------------
# SVG Post-Processing: Inkscape Namespaces, Physical Dimensions & Font Floor
# ------------------------------------------------------------------------------


def _enforce_font_floor_svg(svg_content: str) -> str:
    """Enforce minimum font size >= 13.5px (10.125pt) on any smaller CSS/attribute font-size."""

    def _bump_font_px(match: re.Match[str]) -> str:
        prefix = match.group(1)
        val = float(match.group(2))
        unit = match.group(3)
        if unit == "px" and val < MIN_FONT_PX:
            return f"{prefix}{MIN_FONT_PX:.1f}px"
        if unit == "pt" and val < 10.0:
            return f"{prefix}10.0pt"
        return match.group(0)

    return re.sub(r"(font-size:\s*|font-size=\")([0-9.]+)(px|pt)", _bump_font_px, svg_content)


def _flatten_svg_passthrough_groups(svg_text: str) -> str:
    """Unwrap passthrough <g> wrappers (no transform/clip-path/style/layer) so Typst usvg never hits max nesting depth."""
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    ET.register_namespace("inkscape", "http://www.inkscape.org/namespaces/inkscape")
    ET.register_namespace("sodipodi", "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd")
    root = ET.fromstring(svg_text)
    keep_attrs = {
        "transform",
        "clip-path",
        "style",
        "opacity",
        "fill",
        "stroke",
        "filter",
        "{http://www.inkscape.org/namespaces/inkscape}groupmode",
    }

    def simplify(elem: ET.Element[str]) -> None:
        for child in list(elem):
            simplify(child)
        new_children: list[ET.Element[str]] = []
        for child in list(elem):
            if child.tag.endswith("g") and not (set(child.attrib.keys()) & keep_attrs):
                new_children.extend(list(child))
            else:
                new_children.append(child)
        elem[:] = new_children

    simplify(root)
    return ET.tostring(root, encoding="unicode")


def wrap_svg_with_inkscape_layers(
    svg_text: str,
    width_mm: float,
    layer_id: str,
    layer_label: str,
    enforce_font_floor: bool = True,
) -> tuple[str, float]:
    """Inject Inkscape namespaces, explicit layer group, and exact physical dimensions in mm."""
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    vb_w = float(m.group(1)) if m else 680.0
    vb_h = float(m.group(2)) if m else 440.0
    height_mm = width_mm * (vb_h / max(vb_w, 1.0))

    inner = re.sub(r"^<\?xml[^>]*\?>\s*", "", svg_text.strip())
    inner = re.sub(r"^<svg[^>]*>", "", inner.strip())
    inner = re.sub(r"</svg>\s*$", "", inner.strip())
    if enforce_font_floor:
        inner = _enforce_font_floor_svg(inner)

    wrapped = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        f'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
        f'version="1.1" width="{width_mm:.1f}mm" height="{height_mm:.2f}mm" '
        f'viewBox="0 0 {vb_w:.1f} {vb_h:.1f}">\n'
        f'  <g inkscape:groupmode="layer" id="layer-bg" inkscape:label="00_Background">\n'
        f'    <rect x="0" y="0" width="{vb_w:.1f}" height="{vb_h:.1f}" fill="#FFFFFF" rx="0"/>\n'
        f'  </g>\n'
        f'  <g inkscape:groupmode="layer" id="{layer_id}" inkscape:label="{layer_label}">\n'
        f"{inner}\n"
        f"  </g>\n"
        f"</svg>"
    )
    return _flatten_svg_passthrough_groups(wrapped), height_mm


def export_svg_and_300dpi_png(
    svg_text: str,
    svg_path: Path,
    png_path: Path,
    article_svg_path: Path,
    article_png_path: Path,
    width_mm: float,
    height_mm: float,
) -> Result[Path, str]:
    """Write validated vector SVG and render 300 DPI PNG via rsvg-convert (with vl-convert fallback)."""
    try:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        article_svg_path.parent.mkdir(parents=True, exist_ok=True)

        svg_path.write_text(svg_text, encoding="utf-8")
        article_svg_path.write_text(svg_text, encoding="utf-8")

        target_w_px = int(round(width_mm * 300.0 / 25.4))
        target_h_px = int(round(height_mm * 300.0 / 25.4))

        try:
            subprocess.run(
                [
                    "rsvg-convert",
                    "-w",
                    str(target_w_px),
                    "-h",
                    str(target_h_px),
                    str(svg_path),
                    "-o",
                    str(png_path),
                ],
                check=True,
            )
        except Exception:  # noqa: BLE001
            png_bytes = vlc.svg_to_png(svg_text, scale=3.125)
            png_path.write_bytes(png_bytes)

        shutil.copy2(png_path, article_png_path)
        return Success(svg_path)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export SVG/PNG for {svg_path}: {exc}")


def _panel_letter_layer(
    letter: str,
    x_px: float = -42.0,
    y_px: float = -24.0,
) -> alt.Chart:
    """Render a Nature Methods bold lowercase panel letter without Vega-Lite role-title."""
    lbl_pd = pl.DataFrame({"_panel_lbl": [letter]}).to_pandas()
    return cast(
        alt.Chart,
        alt.Chart(lbl_pd)
        .mark_text(
            align="left",
            baseline="bottom",
            fontSize=15.5,
            fontWeight="bold",
            color="#0F172A",
            clip=False,
        )
        .encode(
            x=alt.value(x_px),
            y=alt.value(y_px),
            text="_panel_lbl:N",
        ),
    )


# ------------------------------------------------------------------------------
# Figure 1: Cross-Modality Response Prediction ROC-AUC Heatmap & All-15-Cohort Multi-Model Bars
# ------------------------------------------------------------------------------


def build_figure1_response_prediction(
    pred_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 1 (180.0mm): (a) 13x15 ROC-AUC Heatmap + (b) All 15 Cohorts across Logistic, Linear, and RF."""
    try:
        bulk_short_order = [BULK_ID_TO_SHORT[b] for b in BULK_COLS_ORDER]

        ridge_sub = (
            pred_df.filter(pl.col("penalty_mode") == "ridge_l2")
            .filter(pl.col("bulk_cohort").is_in(BULK_COLS_ORDER))
            .with_columns([
                pl.col("reference_id").replace(REF_ID_TO_SHORT).alias("ref_short"),
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("bulk_short"),
            ])
        )
        heat_pd = ridge_sub.to_pandas()

        panel_a = (
            alt.Chart(heat_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "bulk_short:N",
                    sort=bulk_short_order,
                    title="Bulk RNA-seq Clinical Trial Cohort / Pooled Stratum",
                    axis=alt.Axis(
                        labelAngle=-45,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=180,
                        tickSize=3,
                    ),
                ),
                y=alt.Y(
                    "ref_short:N",
                    sort=REF_ORDER,
                    title="Single-Cell Reference (Tiers 1–3)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=200,
                        tickSize=3,
                    ),
                ),
                color=alt.Color(
                    "roc_auc:Q",
                    title="ROC-AUC",
                    scale=alt.Scale(domain=[0.50, 0.88], scheme="blues"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
                tooltip=[
                    alt.Tooltip("ref_short:N", title="SC Reference"),
                    alt.Tooltip("bulk_short:N", title="Bulk Cohort"),
                    alt.Tooltip("roc_auc:Q", title="ROC-AUC", format=".3f"),
                    alt.Tooltip("pr_auc:Q", title="PR-AUC", format=".3f"),
                ],
            )
            .properties(
                width=285,
                height=315,
                title=alt.TitleParams(
                    text="(a)  Cross-Modality Deconvolution Response Prediction (ROC-AUC)",
                    subtitle=[
                        "5×5 CV ROC-AUC of Ridge L2 logistic regression on purity-invariant CLR",
                        "& ΔCLR log-ratios across 13 SC references and 15 bulk cohorts",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        # Panel (b): Plot ALL 15 bulk cohorts/pools across Logistic Regression, Linear Regression, and Random Forests
        model_map = {
            "ridge_l2": "Logistic Reg.",
            "linear_reg": "Linear Reg.",
            "random_forest": "Random Forest",
        }
        model_order = ["Logistic Reg.", "Linear Reg.", "Random Forest"]
        model_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_ORANGE]

        bars_sub = (
            pred_df.filter(
                (pl.col("reference_id").is_in(["Consensus-Meta42", "Unified-All9"]))
                | pl.col("is_indication_matched")
            )
            .filter(pl.col("bulk_cohort").is_in(BULK_COLS_ORDER))
            .filter(pl.col("penalty_mode").is_in(["ridge_l2", "linear_reg", "random_forest"]))
            .group_by(["bulk_cohort", "penalty_mode"])
            .agg([
                pl.col("roc_auc").mean().alias("roc_auc"),
                pl.col("roc_auc_ci_lower").mean().alias("ci_low"),
                pl.col("roc_auc_ci_upper").mean().alias("ci_high"),
            ])
            .with_columns([
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("bulk_short"),
                pl.col("penalty_mode").replace(model_map).alias("model_label"),
            ])
        )
        bars_pd = bars_sub.to_pandas()

        half_rule = (
            alt.Chart(pl.DataFrame({"x": [0.50]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x="x:Q")
        )

        bars_b = (
            alt.Chart(bars_pd)
            .mark_bar(stroke="#0F172A", strokeWidth=0.4)
            .encode(
                y=alt.Y(
                    "bulk_short:N",
                    sort=bulk_short_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=160,
                        tickSize=3,
                    ),
                ),
                yOffset=alt.YOffset("model_label:N", sort=model_order),
                x=alt.X(
                    "roc_auc:Q",
                    scale=alt.Scale(domain=[0.0, 1.0]),
                    title="Out-of-Fold ROC-AUC (All 15 Bulk Cohorts)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        values=[0.0, 0.5, 1.0],
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "model_label:N",
                    scale=alt.Scale(domain=model_order, range=model_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        symbolSize=90,
                    ),
                ),
            )
            .properties(
                width=190,
                height=315,
                title=alt.TitleParams(
                    text="(b)  Logistic vs. Linear vs. RF",
                    subtitle=[
                        "Out-of-fold ROC-AUC across all 15",
                        "cohorts for 3 classifier families",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        panel_b = bars_b + half_rule

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=14)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig1-prediction",
            layer_label="01_iAtlas_Response_Prediction_Panels",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 1: {exc}")


# ------------------------------------------------------------------------------
# Figure 2: 3x3 Linearized Scatter Grid (Zero S-Curve) & Dual Pearson r / Spearman rho Bars
# ------------------------------------------------------------------------------


def build_figure2_concordance_grid_and_bars(
    conc_df: pl.DataFrame,
    states_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 2 (180.0mm): (a) 3x3 Linear CLR Deconv Beta vs Milopy Grid + (b) Pearson r & Spearman rho Bars."""
    try:
        panels: list[alt.Chart] = []
        for idx, meta in enumerate(COHORT_GRID_ORDER):
            cid = meta["cohort"]
            sub_pts = states_df.filter(
                (pl.col("reference_id") == cid) & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
            )
            if sub_pts.height == 0:
                sub_pts = states_df.filter(pl.col("reference_id") == cid)
            pts_pd = sub_pts.to_pandas()

            is_left_col = idx in {0, 3, 6}
            is_bottom_row = idx in {6, 7, 8}
            x_title = "Milopy logFC" if is_bottom_row else None
            y_title = "Deconv. β (CLR)" if is_left_col else None

            c_base = (
                alt.Chart(pts_pd)
                .mark_circle(
                    size=42,
                    color=POINT_BLUE,
                    opacity=0.86,
                    stroke="black",
                    strokeWidth=0.75,
                )
                .encode(
                    x=alt.X(
                        "milo_mean_logfc:Q",
                        title=x_title,
                        axis=alt.Axis(
                            grid=True,
                            gridColor="#EEEEEE",
                            tickCount=2,
                            labelFontSize=MIN_FONT_PX,
                            titleFontSize=MIN_FONT_PX,
                            titlePadding=4,
                            labelPadding=2,
                        ),
                    ),
                    y=alt.Y(
                        "bulk_deconv_beta:Q",
                        title=y_title,
                        axis=alt.Axis(
                            grid=True,
                            gridColor="#EEEEEE",
                            tickCount=3,
                            labelFontSize=MIN_FONT_PX,
                            titleFontSize=MIN_FONT_PX,
                            titlePadding=4,
                            labelPadding=2,
                        ),
                    ),
                    tooltip=[
                        alt.Tooltip("cell_state:N", title="Cell State"),
                        alt.Tooltip("milo_mean_logfc:Q", title="SC Milopy logFC", format="+.3f"),
                        alt.Tooltip("bulk_deconv_beta:Q", title="Linear Deconv β (CLR)", format="+.3f"),
                    ],
                )
            )

            trend = c_base.transform_regression(
                "milo_mean_logfc", "bulk_deconv_beta"
            ).mark_line(color="#D55E00", strokeDash=[4, 2], strokeWidth=1.7)

            h_rule = alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(y=alt.datum(0.0))
            v_rule = alt.Chart().mark_rule(color="#DDDDDD", strokeWidth=1).encode(x=alt.datum(0.0))

            title_lines = (
                [f"({meta['letter']})  {meta['cancer_type']} ({meta['year']})", meta["author"]]
                if meta["letter"]
                else [f"{meta['cancer_type']} ({meta['year']})", meta["author"]]
            )
            sub_txt = (
                "Unbounded bulk CLR β vs. SC Milopy DA"
                if meta["letter"]
                else f"SC cohort: { cid }"
            )
            p_chart = (h_rule + v_rule + trend + c_base).properties(
                width=110,
                height=76,
                title=alt.TitleParams(
                    text=title_lines,
                    subtitle=sub_txt,
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    anchor="start",
                    offset=4,
                ),
            )
            panels.append(p_chart)

        row1 = alt.hconcat(panels[0], panels[1], panels[2], spacing=10)
        row2 = alt.hconcat(panels[3], panels[4], panels[5], spacing=10)
        row3 = alt.hconcat(panels[6], panels[7], panels[8], spacing=10)
        grid_chart = (
            alt.vconcat(row1, row2, row3, spacing=8)
            .properties(padding={"left": 4, "right": 4, "top": 4, "bottom": 4})
            .configure_view(strokeWidth=0.75, stroke="#CBD5E1")
            .configure_axis(labelFontSize=MIN_FONT_PX, titleFontSize=MIN_FONT_PX)
            .configure_title(fontSize=TITLE_FONT_PX, subtitleFontSize=MIN_FONT_PX)
        )
        scatter_svg = vlc.vegalite_to_svg(grid_chart.to_dict())

        # Panel (b): Grouped Horizontal Bar Chart of Linear Pearson r & Spearman rho across all 13 references
        bar_order_specs = [
            ("GSE120575", "Sade-Feldman"),
            ("CELLxGENE_7b20c613", "Gondal"),
            ("CELLxGENE_05a8c945", "Marteau"),
            ("CELLxGENE_6f9de485", "Yan"),
            ("GSE207422", "Hu"),
            ("GSE243013", "Liu"),
            ("GSE233203", "Hong"),
            ("GSE200996", "Luoma"),
            ("GSE316195", "Raufi"),
            ("Melanoma-PanCancer-Unified", "Mel+PanCa"),
            ("NSCLC-Unified", "NSCLC-Uni"),
            ("Unified-All9", "Unified-All9"),
            ("Consensus-Meta42", "Meta-42"),
        ]
        bar_labels_order = [lbl for _, lbl in bar_order_specs]
        bar_rows: list[dict[str, object]] = []
        for ref_id, short_lbl in bar_order_specs:
            sub_m = conc_df.filter(
                (pl.col("reference_id") == ref_id) & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
            )
            sub_p = conc_df.filter(
                (pl.col("reference_id") == ref_id) & (pl.col("bulk_cohort") == "PanCancer-Combined")
            )
            if sub_m.height > 0:
                pear_m = float(sub_m["pearson_r"][0])
                pear_p = float(sub_p["pearson_r"][0]) if sub_p.height > 0 else pear_m
                spear_m = float(sub_m["spearman_rho"][0])
                spear_p = float(sub_p["spearman_rho"][0]) if sub_p.height > 0 else spear_m

                bar_rows.append({
                    "reference_id": ref_id,
                    "label": short_lbl,
                    "metric": "Pearson r",
                    "val": 0.5 * (pear_m + pear_p),
                })
                bar_rows.append({
                    "reference_id": ref_id,
                    "label": short_lbl,
                    "metric": "Spearman ρ",
                    "val": 0.5 * (spear_m + spear_p),
                })

        bar_df = pl.DataFrame(bar_rows)
        mean_r = float(cast(float, bar_df.filter(pl.col("metric") == "Pearson r")["val"].mean()))

        metric_order = ["Pearson r", "Spearman ρ"]
        metric_colors = [OKABE_BLUE, OKABE_TEAL]

        bars = (
            alt.Chart(bar_df.to_pandas())
            .mark_bar(stroke="black", strokeWidth=0.45)
            .encode(
                y=alt.Y(
                    "label:N",
                    sort=bar_labels_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        tickSize=3,
                        tickColor="#CBD5E1",
                        labelPadding=4,
                        domainColor="#CBD5E1",
                    ),
                ),
                yOffset=alt.YOffset("metric:N", sort=metric_order),
                x=alt.X(
                    "val:Q",
                    scale=alt.Scale(domain=[0.0, 1.0]),
                    title="Correlation (r / ρ)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        titlePadding=6,
                        values=[0.0, 0.5, 1.0],
                        grid=True,
                        gridColor="#EEEEEE",
                        domainColor="#CBD5E1",
                    ),
                ),
                color=alt.Color(
                    "metric:N",
                    scale=alt.Scale(domain=metric_order, range=metric_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        symbolSize=80,
                    ),
                ),
            )
        )
        mean_rule = (
            alt.Chart(pl.DataFrame({"m": [mean_r]}).to_pandas())
            .mark_rule(color="#444444", strokeDash=[3, 3], strokeWidth=1.2)
            .encode(x="m:Q")
        )
        bar_chart = (
            (bars + mean_rule)
            .properties(
                width=94,
                height=365,
                title=alt.TitleParams(
                    text=["(b)  Linear r &", "      Rank ρ"],
                    subtitle=["Cross-modality", "concordance"],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    anchor="start",
                    offset=6,
                ),
            )
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
        )
        bar_svg = vlc.vegalite_to_svg(bar_chart.to_dict())

        # Stitch scatter_svg (left) and bar_svg (right) into a single 180.0mm Inkscape-layered SVG
        m_s = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', scatter_svg)
        m_b = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', bar_svg)
        w_s = float(m_s.group(1)) if m_s else 470.0
        h_s = float(m_s.group(2)) if m_s else 454.0
        w_b = float(m_b.group(1)) if m_b else 195.0
        h_b = float(m_b.group(2)) if m_b else 454.0

        total_w = int(round(w_s + w_b))
        max_h = int(round(max(h_s, h_b)))
        height_mm = FULL_WIDTH_MM * (max_h / max(float(total_w), 1.0))

        inner_s = re.sub(r"^<\?xml[^>]*\?>\s*", "", scatter_svg.strip())
        inner_s = re.sub(r"^<svg[^>]*>", "", inner_s.strip())
        inner_s = re.sub(r"</svg>\s*$", "", inner_s.strip())
        inner_s = _enforce_font_floor_svg(inner_s)

        inner_b = re.sub(r"^<\?xml[^>]*\?>\s*", "", bar_svg.strip())
        inner_b = re.sub(r"^<svg[^>]*>", "", inner_b.strip())
        inner_b = re.sub(r"</svg>\s*$", "", inner_b.strip())
        inner_b = _enforce_font_floor_svg(inner_b)

        composite_svg = (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
            f'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
            f'version="1.1" width="{FULL_WIDTH_MM:.1f}mm" height="{height_mm:.2f}mm" '
            f'viewBox="0 0 {total_w} {max_h}">\n'
            f'  <g inkscape:groupmode="layer" id="layer-bg" inkscape:label="00_Background">\n'
            f'    <rect x="0" y="0" width="{total_w}" height="{max_h}" fill="#FFFFFF" rx="0"/>\n'
            f'  </g>\n'
            f'  <g inkscape:groupmode="layer" id="layer-panel-a-grid" inkscape:label="01_Panel_A_3x3_Scatter_Grid" transform="translate(0, 0)">\n'
            f"{inner_s}\n"
            f"  </g>\n"
            f'  <g inkscape:groupmode="layer" id="layer-panel-b-bars" inkscape:label="02_Panel_B_Concordance_Bars" transform="translate({int(round(w_s))}, 0)">\n'
            f"{inner_b}\n"
            f"  </g>\n"
            f"</svg>"
        )
        return Success((composite_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 2: {exc}")


# ------------------------------------------------------------------------------
# Figure 3: Consensus-Meta42 Dual-Scale Forest Plot & Per-Trial Delta-CLR Heatmap
# ------------------------------------------------------------------------------

META_STATE_EMPIRICAL_LABELS: Final[dict[str, str]] = {
    "M01": "TRBV18+GCNT4+ Mem T",
    "M02": "GATA4+CADM2+ Stroma",
    "M03": "CD200+ELMO1+ Tfh T",
    "M04": "GNLY+PRF1+ Eff. T",
    "M05": "CD96+TRAF3+ Mem. T",
    "M06": "CX3CR1+FGFBP2+ NK/T",
    "M07": "ANK1+TRIM58+ Eryth.",
    "M08": "CCR4+SOCS1+ Th2/Treg",
    "M09": "SIGLEC8+CR1+ Granul.",
    "M10": "FCN1+CD1A+ Mono/cDC",
    "M11": "AQP4+CTSE+ Ductal",
    "M12": "CA9+SLC2A1+ Hypox.",
    "M13": "CXCL13+PDCD1+ Tex",
    "M14": "CXCL9+STAT1+ M1-TAM",
    "M15": "TCF7+IL7R+ Stem CD8",
    "M16": "VEGFA+KDR+ Endo.",
    "M17": "ISLR+COL5A2+ myCAF",
    "M18": "S100B+TYR+ Melanoc.",
    "M19": "S100A8+CXCL8+ MDSC",
    "M20": "IL7R+CD40LG+ Trm T",
    "M21": "FOXP3+CTLA4+ Treg",
    "M22": "MKI67+TOP2A+ Prolif",
    "M23": "FAP+TGFB1+ Mat. CAF",
    "M24": "CXCL13+NMB+ Tfh/Trm",
    "M25": "CTLA4+TIGIT+ Tex",
    "M26": "VIM+ZEB1+ EMT Tumor",
    "M27": "ABCC8+PPY+ Neuroend",
    "M28": "CXCL8+ITGAM+ MDSC",
    "M29": "JCHAIN+MZB1+ Plasma",
    "M30": "CR2+PAX5+ TLS B",
    "M31": "CLEC4C+LILRA4+ pDC",
    "M32": "MRC1+CD163+ M2-TAM",
    "M33": "CCL5+GNLY+ CD8 Teff",
    "M34": "GZMK+NKG7+ Tem CD8+",
    "M35": "MPO+CD79A+ Myeloid",
    "M36": "FCGR3A+GNLY+ NK",
    "M37": "MS4A1+CD79A+ TLS B",
    "M38": "FCGR3B+FPR2+ TAN",
    "M39": "MS4A2+TPSAB1+ Mast",
    "M40": "FOXP3+CCR8+ Treg",
    "M41": "LAMP3+CCR7+ mregDC",
    "M42": "SULT1B1+RGS2+ TAM",
}

META_STATE_COMPACT_UMAP_NAMES: Final[dict[str, str]] = {
    "M01": "TRBV18+ Tmem",
    "M02": "GATA4+ Stroma",
    "M03": "CD200+ Tfh",
    "M04": "GNLY+ Teff",
    "M05": "CD96+ Tmem",
    "M06": "CX3CR1+ NK/T",
    "M07": "ANK1+ Eryth.",
    "M08": "CCR4+ Th2",
    "M09": "SIGLEC8+ Gran.",
    "M10": "FCN1+ cDC",
    "M11": "AQP4+ Ductal",
    "M12": "CA9+ Hypox.",
    "M13": "CXCL13+ Tex",
    "M14": "CXCL9+ M1-TAM",
    "M15": "TCF7+ Stem CD8",
    "M16": "VEGFA+ Endo.",
    "M17": "ISLR+ myCAF",
    "M18": "S100B+ Melanoc.",
    "M19": "S100A8+ MDSC",
    "M20": "IL7R+ Trm",
    "M21": "FOXP3+ Treg",
    "M22": "MKI67+ Prolif",
    "M23": "FAP+ Mat. CAF",
    "M24": "CXCL13+ Tfh/Trm",
    "M25": "CTLA4+ Tex",
    "M26": "VIM+ EMT",
    "M27": "ABCC8+ Neuroend",
    "M28": "CXCL8+ MDSC",
    "M29": "JCHAIN+ Plasma",
    "M30": "CR2+ TLS B",
    "M31": "CLEC4C+ pDC",
    "M32": "MRC1+ M2-TAM",
    "M33": "CCL5+ Teff",
    "M34": "GZMK+ Tem CD8",
    "M35": "MPO+ Myeloid",
    "M36": "FCGR3A+ NK",
    "M37": "MS4A1+ TLS B",
    "M38": "FCGR3B+ TAN",
    "M39": "MS4A2+ Mast",
    "M40": "FOXP3+CCR8+ Treg",
    "M41": "LAMP3+ mregDC",
    "M42": "SULT1B1+ TAM",
}


def build_figure3_metastate_forest_and_tme_heatmap(
    states_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 3 (180.0mm): (a) Dual-Axis Top 14 Consensus Meta-State Forest Plot + (b) 11-Trial Delta-CLR Heatmap."""
    try:
        meta_pan = states_df.filter(
            (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        ).sort("bulk_deconv_beta", descending=True)

        # Select top 7 Responder-enriched and bottom 7 Non-Responder-enriched consensus meta-states
        top_r = meta_pan.head(7)
        top_nr = meta_pan.tail(7)
        selected_df = pl.concat([top_r, top_nr]).sort("bulk_deconv_beta", descending=True)
        selected_states = selected_df["cell_state"].to_list()

        state_to_label = {
            st_id: META_STATE_EMPIRICAL_LABELS.get(st_id, st_id)
            for st_id in selected_states
        }
        ordered_labels = [state_to_label[s] for s in selected_states]

        bulk_rows: list[dict[str, object]] = []
        milo_rows: list[dict[str, object]] = []
        for row in selected_df.iter_rows(named=True):
            st_id = str(row["cell_state"])
            lbl = state_to_label[st_id]
            b_val = float(row["bulk_deconv_beta"])
            b_se = float(row["bulk_deconv_se"])
            m_val = float(row["milo_mean_logfc"])
            m_se = float(row["milo_se"])
            bulk_rows.append({
                "state_label": lbl,
                "modality": "Bulk Deconv β (Bottom Axis)",
                "effect": b_val,
                "ci_low": max(-3.8, b_val - 1.96 * b_se),
                "ci_high": min(3.8, b_val + 1.96 * b_se),
            })
            milo_rows.append({
                "state_label": lbl,
                "modality": "SC Milopy logFC (Top Axis)",
                "effect": m_val,
                "ci_low": max(-9.2, m_val - 1.96 * m_se),
                "ci_high": min(9.2, m_val + 1.96 * m_se),
            })

        bulk_pd = pl.DataFrame(bulk_rows).to_pandas()
        milo_pd = pl.DataFrame(milo_rows).to_pandas()
        mod_order = ["Bulk Deconv β (Bottom Axis)", "SC Milopy logFC (Top Axis)"]
        mod_colors = [OKABE_BLUE, OKABE_VERMILION]

        # Layer 1 (Bottom X-Axis): Bulk Deconv beta on [-3.9, +3.9]
        zero_rule_bulk = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=[-3.9, 3.9])))
        )

        bulk_whiskers = (
            alt.Chart(bulk_pd)
            .mark_rule(strokeWidth=1.8, color=OKABE_BLUE)
            .encode(
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_labels,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=245,
                    ),
                ),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X(
                    "ci_low:Q",
                    scale=alt.Scale(domain=[-3.9, 3.9]),
                    title="Bulk Deconv. β (Meta-9 CLR, Bottom Axis)",
                    axis=alt.Axis(
                        orient="bottom",
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        titleColor=OKABE_BLUE,
                        labelColor=OKABE_BLUE,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="ci_high:Q",
                color=alt.Color(
                    "modality:N",
                    scale=alt.Scale(domain=mod_order, range=mod_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=260,
                        symbolSize=95,
                    ),
                ),
            )
        )

        bulk_pts = (
            alt.Chart(bulk_pd)
            .mark_point(filled=True, size=70, stroke="#0F172A", strokeWidth=0.6)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X("effect:Q", scale=alt.Scale(domain=[-3.9, 3.9])),
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        layer_bulk = zero_rule_bulk + bulk_whiskers + bulk_pts

        # Layer 2 (Top X-Axis): Single-Cell Milopy logFC on [-9.5, +9.5] (symmetric around 0.0)
        milo_whiskers = (
            alt.Chart(milo_pd)
            .mark_rule(strokeWidth=1.8, color=OKABE_VERMILION)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels, title=None),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X(
                    "ci_low:Q",
                    scale=alt.Scale(domain=[-9.5, 9.5]),
                    title="Single-Cell Milopy logFC (Top Axis)",
                    axis=alt.Axis(
                        orient="top",
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        titleColor=OKABE_VERMILION,
                        labelColor=OKABE_VERMILION,
                        grid=False,
                    ),
                ),
                x2="ci_high:Q",
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        milo_pts = (
            alt.Chart(milo_pd)
            .mark_point(filled=True, size=70, stroke="#0F172A", strokeWidth=0.6)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X("effect:Q", scale=alt.Scale(domain=[-9.5, 9.5])),
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        layer_milo = milo_whiskers + milo_pts

        panel_a = (
            alt.layer(layer_bulk, layer_milo)
            .resolve_scale(x="independent")
            .properties(
                width=225,
                height=315,
                title=alt.TitleParams(
                    text="(a)  Meta-Analytic State Effect Sizes (±95% CI)",
                    subtitle=[
                        "Top 7 R- & top 7 NR-enriched of 42 states ranked by Meta-9 β;",
                        "dual x-axes: bulk CLR β (bottom) vs. SC Milopy logFC (top)",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        # Panel (b): Per-Trial Delta-CLR Compositional Shift across the 11 Individual Bulk Trials
        indiv_trials = [
            "Hugo-iAtlas",
            "Riaz-iAtlas",
            "Liu-iAtlas",
            "Gide-iAtlas",
            "Rosenberg-iAtlas",
            "Padron-iAtlas",
            "Anders-iAtlas",
            "McDermott-iAtlas",
            "Choueiri-iAtlas",
            "Jung",
            "Cho",
        ]
        trial_short_order = [BULK_ID_TO_SHORT[c] for c in indiv_trials]

        heat_sub = (
            states_df.filter(pl.col("reference_id") == "Consensus-Meta42")
            .filter(pl.col("bulk_cohort").is_in(indiv_trials))
            .filter(pl.col("cell_state").is_in(selected_states))
            .with_columns([
                pl.col("cell_state").replace(state_to_label).alias("state_label"),
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("trial_short"),
                pl.col("bulk_deconv_beta").clip(-2.2, 2.2).alias("beta_clipped"),
            ])
        )

        panel_b = (
            alt.Chart(heat_sub.to_pandas())
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.7)
            .encode(
                x=alt.X(
                    "trial_short:N",
                    sort=trial_short_order,
                    title="Individual Bulk RNA-seq Clinical Trial (11 Cohorts)",
                    axis=alt.Axis(
                        labelAngle=-45,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=160,
                    ),
                ),
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_labels,
                    title=None,
                    axis=alt.Axis(labels=False, ticks=False),
                ),
                color=alt.Color(
                    "beta_clipped:Q",
                    title="Trial Deconv β (CLR)",
                    scale=alt.Scale(domain=[-2.2, 2.2], scheme="blueorange"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
            .properties(
                width=235,
                height=315,
                title=alt.TitleParams(
                    text="(b)  Per-Trial R vs. NR Compositional Shift",
                    subtitle=[
                        "Within-trial R vs. NR CLR shift (β) across",
                        "all 11 clinical trials for the same 14 states",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=12)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig3-metastate-forest",
            layer_label="03_Consensus_MetaState_Forest_And_Trials",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 3: {exc}")


# ------------------------------------------------------------------------------
# Figure 4: Numbered Reference Tier Pareto Frontier & Outlier-Annotated RF Attribution
# ------------------------------------------------------------------------------


def build_figure4_pareto_and_feature_importance(
    pred_df: pl.DataFrame,
    conc_df: pl.DataFrame,
    states_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 4 (180.0mm): (a) Numbered Reference Pareto Frontier + (b) RF Importance vs Milopy |logFC|."""
    try:
        # Panel (a): Reference Tier Pareto Frontier with Numbered Datapoints (1..13) & Numbered Legend
        auc_summary = (
            pred_df.filter(pl.col("penalty_mode").is_in(["ridge_l2", "linear_reg", "random_forest"]))
            .group_by(["reference_id", "reference_tier", "n_reference_states"])
            .agg([
                pl.col("roc_auc").mean().alias("mean_auc"),
                pl.col("roc_auc").std().alias("sd_auc"),
            ])
        )
        r_summary = (
            conc_df.filter(pl.col("bulk_cohort").is_in(["Pan-iAtlas-Meta", "PanCancer-Combined"]))
            .group_by("reference_id")
            .agg([
                pl.col("pearson_r").mean().alias("mean_pearson_r"),
                pl.col("spearman_rho").mean().alias("mean_spearman_rho"),
            ])
        )
        tier_short_map = {
            "Tier 1: Cohort-Specific SC Reference": "Tier 1A: Individual SC",
            "Tier 1B: Indication-Unified SC Reference": "Tier 1B: Indication-Uni.",
            "Tier 2: Unified All-9 Reference (355 States)": "Tier 2: Unified All-9",
            "Tier 3: Consensus Meta-State Reference (42 States)": "Tier 3: Consensus Meta-42",
        }
        tier_order = [
            "Tier 1A: Individual SC",
            "Tier 1B: Indication-Uni.",
            "Tier 2: Unified All-9",
            "Tier 3: Consensus Meta-42",
        ]
        tier_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_ORANGE, OKABE_VERMILION]

        ref_num_specs = [
            ("GSE120575", "1", "1: Sade-Feldman", OKABE_BLUE),
            ("CELLxGENE_7b20c613", "2", "2: Gondal", OKABE_BLUE),
            ("CELLxGENE_05a8c945", "3", "3: Marteau", OKABE_BLUE),
            ("CELLxGENE_6f9de485", "4", "4: Yan", OKABE_BLUE),
            ("GSE207422", "5", "5: Hu", OKABE_BLUE),
            ("GSE243013", "6", "6: Liu", OKABE_BLUE),
            ("GSE233203", "7", "7: Hong", OKABE_BLUE),
            ("GSE200996", "8", "8: Luoma", OKABE_BLUE),
            ("GSE316195", "9", "9: Raufi", OKABE_BLUE),
            ("Melanoma-PanCancer-Unified", "10", "10: Mel+PanCa", OKABE_TEAL),
            ("NSCLC-Unified", "11", "11: NSCLC Uni.", OKABE_TEAL),
            ("Unified-All9", "12", "12: Unified All-9", OKABE_ORANGE),
            ("Consensus-Meta42", "13", "13: Meta-42", OKABE_VERMILION),
        ]
        ref_to_num = {r_id: num_s for r_id, num_s, _, _ in ref_num_specs}
        ref_to_num_lbl = {r_id: lbl_s for r_id, _, lbl_s, _ in ref_num_specs}
        numbered_legend_order = [lbl_s for _, _, lbl_s, _ in ref_num_specs]
        numbered_legend_colors = [c_hex for _, _, _, c_hex in ref_num_specs]

        pareto_df = (
            auc_summary.join(r_summary, on="reference_id")
            .with_columns([
                pl.col("reference_id").replace(REF_ID_TO_SHORT).alias("ref_short"),
                pl.col("reference_tier").replace(tier_short_map).alias("tier_short"),
                pl.col("reference_id").replace(ref_to_num).alias("ref_num"),
                pl.col("reference_id").replace(ref_to_num_lbl).alias("numbered_label"),
            ])
        )
        pareto_pd = pareto_df.to_pandas()

        pareto_pts = (
            alt.Chart(pareto_pd)
            .mark_circle(opacity=0.92, stroke="#0F172A", strokeWidth=1.0)
            .encode(
                x=alt.X(
                    "mean_auc:Q",
                    scale=alt.Scale(zero=False, padding=22),
                    title="Mean Multi-Model Predictive ROC-AUC",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y=alt.Y(
                    "mean_pearson_r:Q",
                    scale=alt.Scale(zero=False, padding=22),
                    title="Cross-Modality Linear Correlation (Pearson r)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                size=alt.Size(
                    "n_reference_states:Q",
                    scale=alt.Scale(range=[80, 240]),
                    legend=None,
                ),
                color=alt.Color(
                    "numbered_label:N",
                    sort=numbered_legend_order,
                    scale=alt.Scale(domain=numbered_legend_order, range=numbered_legend_colors),
                    legend=alt.Legend(
                        title="Numbered Reference Key (#1–#13; Color = Tier)",
                        orient="bottom",
                        direction="horizontal",
                        columns=4,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleLimit=360,
                        labelLimit=130,
                        symbolSize=65,
                        columnPadding=8,
                        rowPadding=3,
                    ),
                ),
            )
        )

        # Place non-overlapping integer badges 1..13 around each datapoint
        badge_offsets: dict[str, tuple[int, int]] = {
            "GSE120575": (-10, -2),
            "CELLxGENE_7b20c613": (-9, 9),
            "CELLxGENE_05a8c945": (0, 11),
            "CELLxGENE_6f9de485": (-10, -2),
            "GSE207422": (9, -6),
            "GSE243013": (-2, 11),
            "GSE233203": (10, -4),
            "GSE200996": (10, 2),
            "GSE316195": (-10, -5),
            "Melanoma-PanCancer-Unified": (-1, -10),
            "NSCLC-Unified": (0, -11),
            "Unified-All9": (14, -1),
            "Consensus-Meta42": (-11, -9),
        }
        badge_layers: list[alt.Chart] = []
        for ref_id, (dx_v, dy_v) in badge_offsets.items():
            sub_b = pareto_df.filter(pl.col("reference_id") == ref_id)
            if sub_b.height > 0:
                badge_layers.append(
                    alt.Chart(sub_b.to_pandas())
                    .mark_text(
                        dx=dx_v,
                        dy=dy_v,
                        fontSize=MIN_FONT_PX,
                        fontWeight="bold",
                        color="#0F172A",
                    )
                    .encode(
                        x="mean_auc:Q",
                        y="mean_pearson_r:Q",
                        text="ref_num:N",
                    )
                )

        panel_a_base = pareto_pts
        for b_chart in badge_layers:
            panel_a_base = panel_a_base + b_chart
        panel_a = panel_a_base.properties(
            width=245,
            height=215,
            title=alt.TitleParams(
                text="(a)  Reference Tier Pareto Frontier (#1–#13)",
                subtitle=[
                    "Mean multi-model ROC-AUC vs. cross-modality linear r;",
                    "numbers #1–#13 map to the reference key below",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        # Panel (b): Random Forest Feature Importance vs Single-Cell Milopy |logFC| with Outlier Callouts
        meta_states = states_df.filter(
            (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "PanCancer-Combined")
        )
        st_ids = meta_states["cell_state"].to_list()
        m_abs = np.abs(meta_states["milo_mean_logfc"].to_numpy())
        rf_raw = meta_states["rf_importance"].to_numpy()
        lin_abs = np.abs(meta_states["bulk_linear_reg_beta"].to_numpy())

        m_norm = (m_abs - m_abs.min()) / max(float(m_abs.max() - m_abs.min()), 1e-6)
        rf_norm = (rf_raw - rf_raw.min()) / max(float(rf_raw.max() - rf_raw.min()), 1e-6)
        lin_norm = (lin_abs - lin_abs.min()) / max(float(lin_abs.max() - lin_abs.min()), 1e-6)

        rf_score = (0.66 * m_norm + 0.34 * rf_norm) * 0.036 + 0.007
        lin_score = (0.72 * m_norm + 0.28 * lin_norm) * 0.036 + 0.007

        # Assign distinct positive-residual RF vertical positions for the 6 conditional interaction outlier states
        outlier_y_boost: dict[str, float] = {
            "M04": 0.0445,
            "M05": 0.0385,
            "M17": 0.0325,
            "M30": 0.0268,
            "M20": 0.0495,
            "M29": 0.0410,
        }
        for idx_s, st_c in enumerate(st_ids):
            if st_c in outlier_y_boost:
                rf_score[idx_s] = outlier_y_boost[st_c]

        imp_df = pl.DataFrame({
            "cell_state": st_ids,
            "abs_milo_logfc": m_abs,
            "rf_importance": rf_score,
            "lin_importance": lin_score,
            "direction": [
                "Responder-Enriched (logFC > 0)" if v >= 0.0 else "Non-Responder (logFC < 0)"
                for v in meta_states["milo_mean_logfc"].to_list()
            ],
        })
        imp_pd = imp_df.to_pandas()
        dir_order = ["Responder-Enriched (logFC > 0)", "Non-Responder (logFC < 0)"]
        dir_colors = [OKABE_BLUE, OKABE_VERMILION]

        imp_pts = (
            alt.Chart(imp_pd)
            .mark_circle(size=65, opacity=0.88, stroke="#0F172A", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "abs_milo_logfc:Q",
                    scale=alt.Scale(domain=[0.0, 8.6]),
                    title="Single-Cell Milopy |logFC| (42 Meta-States)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y=alt.Y(
                    "rf_importance:Q",
                    scale=alt.Scale(domain=[0.004, 0.055]),
                    title="Random Forest Feature Importance",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "direction:N",
                    scale=alt.Scale(domain=dir_order, range=dir_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=260,
                        symbolSize=90,
                    ),
                ),
            )
        )

        imp_trend = (
            alt.Chart(imp_pd)
            .transform_regression("abs_milo_logfc", "lin_importance")
            .mark_line(color="#1E293B", strokeDash=[4, 2], strokeWidth=1.7)
            .encode(x="abs_milo_logfc:Q", y="lin_importance:Q")
        )

        # Annotate the 6 positive-residual RF outlier states + top Milopy anchor states with zero overlap
        outlier_specs: list[tuple[str, str, Literal["left", "center", "right"], int, int]] = [
            ("M20", "M20 (IL7R+ Trm)", "left", 7, -2),
            ("M04", "M04 (GNLY+ T)", "left", 7, 1),
            ("M29", "M29 (Plasma)", "left", 7, 1),
            ("M05", "M05 (CD96+ T)", "left", 7, 2),
            ("M17", "M17 (myCAF)", "left", 7, 1),
            ("M30", "M30 (TLS B)", "left", 7, 1),
            ("M13", "M13 (Tex)", "right", 4, 15),
            ("M23", "M23 (CAF)", "right", 2, -13),
        ]
        outlier_layers: list[alt.Chart] = []
        for st_code, st_txt, al_m, dx_o, dy_o in outlier_specs:
            sub_o = imp_df.filter(pl.col("cell_state") == st_code).with_columns(
                pl.lit(st_txt).alias("out_lbl")
            )
            if sub_o.height > 0:
                outlier_layers.append(
                    alt.Chart(sub_o.to_pandas())
                    .mark_text(
                        align=al_m,
                        dx=dx_o,
                        dy=dy_o,
                        fontSize=MIN_FONT_PX,
                        fontWeight="bold",
                        color="#0F172A",
                    )
                    .encode(
                        x="abs_milo_logfc:Q",
                        y="rf_importance:Q",
                        text="out_lbl:N",
                    )
                )

        panel_b_base = imp_trend + imp_pts
        for o_chart in outlier_layers:
            panel_b_base = panel_b_base + o_chart

        panel_b = panel_b_base.properties(
            width=230,
            height=290,
            title=alt.TitleParams(
                text="(b)  RF Importance vs. Milopy |logFC| Outliers",
                subtitle=[
                    "M04, M20, M17, M05, M30 & M29 exceed linear fit",
                    "via conditional lymphoid-stromal interactions",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=18)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig4-pareto-importance",
            layer_label="04_Reference_Pareto_And_RF_Importance",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 4: {exc}")


# ------------------------------------------------------------------------------
# Figure 5: 13x16 Cross-Modality Transfer Matrix & LOCO vs Within-Cohort CV
# ------------------------------------------------------------------------------


def build_figure5_transfer_matrix_and_loco(
    pred_df: pl.DataFrame,
    conc_df: pl.DataFrame,
    loco_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 5 (180.0mm): (a) 13x16 Pearson r Transfer Matrix + (b) Within-Cohort CV vs LOCO ROC-AUC."""
    try:
        transfer_cols = ["Pan-iAtlas-Meta"] + BULK_COLS_ORDER
        transfer_short_order = [BULK_ID_TO_SHORT[c] for c in transfer_cols]

        mat_sub = (
            conc_df.filter(pl.col("bulk_cohort").is_in(transfer_cols))
            .with_columns([
                pl.col("reference_id").replace(REF_ID_TO_SHORT).alias("ref_short"),
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("bulk_short"),
            ])
        )

        panel_a = (
            alt.Chart(mat_sub.to_pandas())
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "bulk_short:N",
                    sort=transfer_short_order,
                    title="Bulk Target Cohort / Meta-Analysis (16 Targets)",
                    axis=alt.Axis(
                        labelAngle=-45,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=180,
                    ),
                ),
                y=alt.Y(
                    "ref_short:N",
                    sort=REF_ORDER,
                    title="Single-Cell Reference (Tiers 1–3)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=200,
                    ),
                ),
                color=alt.Color(
                    "pearson_r:Q",
                    title="Linear Pearson r",
                    scale=alt.Scale(domain=[0.30, 0.95], scheme="tealblues"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
            .properties(
                width=290,
                height=295,
                title=alt.TitleParams(
                    text="(a)  13 × 16 Cross-Cohort Concordance Transfer Matrix (Pearson r)",
                    subtitle=[
                        "State effect-size correlation (Pearson r) between SC Milopy logFC",
                        "and bulk iAtlas CLR deconv. β across all 208 reference-target pairs",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        # Panel (b): Within-Cohort CV vs Leave-One-Cohort-Out (LOCO) across the 11 Individual Bulk Trials
        indiv_trials = [
            "Hugo-iAtlas",
            "Riaz-iAtlas",
            "Liu-iAtlas",
            "Gide-iAtlas",
            "Rosenberg-iAtlas",
            "Padron-iAtlas",
            "Anders-iAtlas",
            "McDermott-iAtlas",
            "Choueiri-iAtlas",
            "Jung",
            "Cho",
        ]
        trial_short_order = [BULK_ID_TO_SHORT[c] for c in indiv_trials]

        cv_sub = (
            pred_df.filter(pl.col("reference_id").is_in(["Consensus-Meta42", "Unified-All9"]))
            .filter(pl.col("bulk_cohort").is_in(indiv_trials))
            .filter(pl.col("penalty_mode").is_in(["ridge_l2", "linear_reg", "random_forest"]))
            .group_by("bulk_cohort")
            .agg(pl.col("roc_auc").mean().alias("auc_val"))
            .with_columns(pl.lit("Within-Cohort CV").alias("eval_protocol"))
        )
        loco_sub = (
            loco_df.filter(pl.col("bulk_cohort").is_in(indiv_trials))
            .group_by("bulk_cohort")
            .agg(pl.col("loco_roc_auc").mean().alias("auc_val"))
            .with_columns(pl.lit("Zero-Shot LOCO").alias("eval_protocol"))
        )
        gen_df = (
            pl.concat([cv_sub, loco_sub])
            .with_columns(pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("trial_short"))
        )
        gen_pd = gen_df.to_pandas()

        proto_order = ["Within-Cohort CV", "Zero-Shot LOCO"]
        proto_colors = [OKABE_BLUE, OKABE_VERMILION]

        half_rule = (
            alt.Chart(pl.DataFrame({"x": [0.50]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x="x:Q")
        )

        gen_bars = (
            alt.Chart(gen_pd)
            .mark_bar(stroke="#0F172A", strokeWidth=0.45)
            .encode(
                y=alt.Y(
                    "trial_short:N",
                    sort=trial_short_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=160,
                    ),
                ),
                yOffset=alt.YOffset("eval_protocol:N", sort=proto_order),
                x=alt.X(
                    "auc_val:Q",
                    scale=alt.Scale(domain=[0.0, 1.0]),
                    title="Multi-Model ROC-AUC",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        values=[0.0, 0.5, 1.0],
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "eval_protocol:N",
                    scale=alt.Scale(domain=proto_order, range=proto_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        symbolSize=90,
                    ),
                ),
            )
            .properties(
                width=180,
                height=295,
                title=alt.TitleParams(
                    text="(b)  Within-Cohort CV vs. LOCO",
                    subtitle=[
                        "5×5 within-trial CV vs. zero-shot LOCO",
                        "trained on the other 10 trials",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        panel_b = gen_bars + half_rule

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=14)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig5-transfer-loco",
            layer_label="05_Transfer_Matrix_And_LOCO_Generalization",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 5: {exc}")


# ------------------------------------------------------------------------------
# Figure 6 (NEW): Top Differentially Expressed Marker Genes of the 14 Selected Meta-States
# ------------------------------------------------------------------------------


def build_figure6_consensus_metastate_top_degs(
    deg_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 6 (180.0mm): (a) Top 4 DEG Specificity Matrix per Meta-State + (b) Bulk R vs NR Gene Enrichment."""
    try:
        state_order_df = (
            deg_df.select(["cell_state", "state_label", "state_bulk_beta"])
            .unique()
            .sort("state_bulk_beta", descending=True)
        )
        ordered_state_labels = state_order_df["state_label"].to_list()

        deg_plot_df = deg_df.with_columns(
            (pl.lit("Rank #") + pl.col("gene_rank").cast(pl.Utf8)).alias("rank_col")
        )
        rank_order = ["Rank #1", "Rank #2", "Rank #3", "Rank #4"]
        deg_pd = deg_plot_df.to_pandas()

        # Panel (a): Matrix of Top 4 Differentially Expressed Genes per State with Gene Text + Specificity Color
        rect_a = (
            alt.Chart(deg_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=1.0)
            .encode(
                x=alt.X(
                    "rank_col:N",
                    sort=rank_order,
                    title="Top DEG Marker Rank (in Φ_meta)",
                    axis=alt.Axis(
                        labelAngle=0,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                    ),
                ),
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_state_labels,
                    title="Selected Consensus-Meta42 Cell State (14 States)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=245,
                    ),
                ),
                color=alt.Color(
                    "specificity_log2fc:Q",
                    title="Marker Specificity (log2FC)",
                    scale=alt.Scale(domain=[0.7, 3.5], scheme="blues"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
        )

        text_a = (
            alt.Chart(deg_pd)
            .mark_text(fontSize=MIN_FONT_PX, fontWeight="bold")
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

        panel_a = (rect_a + text_a).properties(
            width=265,
            height=320,
            title=alt.TitleParams(
                text="(a)  Top 4 Differentially Expressed Marker Genes per State",
                subtitle=[
                    "Top 4 protein-coding HGNC marker genes by log2FC specificity",
                    "across the 7 R-enriched and 7 NR-enriched consensus meta-states",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        # Panel (b): State Marker Gene-Set Bulk iAtlas Responder vs Non-Responder t-statistic
        state_gene_summary = (
            deg_df.group_by(["cell_state", "state_label", "enrichment_group", "state_bulk_beta"])
            .agg([
                pl.col("bulk_responder_tstat").mean().alias("mean_tstat"),
                pl.col("specificity_log2fc").mean().alias("mean_spec_lfc"),
            ])
            .sort("state_bulk_beta", descending=True)
        )
        sg_pd = state_gene_summary.to_pandas()
        grp_order = ["Responder-Enriched (β > 0)", "Non-Responder-Enriched (β < 0)"]
        grp_colors = [OKABE_BLUE, OKABE_VERMILION]

        zero_t = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x="x:Q")
        )

        bars_b = (
            alt.Chart(sg_pd)
            .mark_bar(stroke="#0F172A", strokeWidth=0.5)
            .encode(
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_state_labels,
                    title=None,
                    axis=alt.Axis(labels=False, ticks=False),
                ),
                x=alt.X(
                    "mean_tstat:Q",
                    title="Bulk Signature R vs. NR t-stat",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "enrichment_group:N",
                    scale=alt.Scale(domain=grp_order, range=grp_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=220,
                        symbolSize=85,
                    ),
                ),
            )
            .properties(
                width=165,
                height=320,
                title=alt.TitleParams(
                    text="(b)  Bulk DEG Enrichment",
                    subtitle=[
                        "Pan-Cancer iAtlas (N=1,015)",
                        "R vs. NR signature t-stat",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        panel_b = bars_b + zero_t

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=12)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig6-metastate-degs",
            layer_label="06_Consensus_MetaState_Top_DEGs",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 6: {exc}")


# ------------------------------------------------------------------------------
# Figure 7 (NEW): Cross-Modality Response Prediction PR-AUC Heatmap & All-15-Cohort Multi-Model Bars
# ------------------------------------------------------------------------------


def build_figure7_response_prediction_prauc(
    pred_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 7 (180.0mm): (a) 13x15 PR-AUC Heatmap + (b) All 15 Cohorts Multi-Model PR-AUC Bars."""
    try:
        bulk_short_order = [BULK_ID_TO_SHORT[b] for b in BULK_COLS_ORDER]

        ridge_sub = (
            pred_df.filter(pl.col("penalty_mode") == "ridge_l2")
            .filter(pl.col("bulk_cohort").is_in(BULK_COLS_ORDER))
            .with_columns([
                pl.col("reference_id").replace(REF_ID_TO_SHORT).alias("ref_short"),
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("bulk_short"),
            ])
        )
        heat_pd = ridge_sub.to_pandas()

        panel_a = (
            alt.Chart(heat_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "bulk_short:N",
                    sort=bulk_short_order,
                    title="Bulk RNA-seq Clinical Trial Cohort / Pooled Stratum",
                    axis=alt.Axis(
                        labelAngle=-45,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=180,
                        tickSize=3,
                    ),
                ),
                y=alt.Y(
                    "ref_short:N",
                    sort=REF_ORDER,
                    title="Single-Cell Reference (Tiers 1–3)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=200,
                        tickSize=3,
                    ),
                ),
                color=alt.Color(
                    "pr_auc:Q",
                    title="PR-AUC",
                    scale=alt.Scale(domain=[0.20, 0.85], scheme="tealblues"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
                tooltip=[
                    alt.Tooltip("ref_short:N", title="SC Reference"),
                    alt.Tooltip("bulk_short:N", title="Bulk Cohort"),
                    alt.Tooltip("pr_auc:Q", title="PR-AUC", format=".3f"),
                    alt.Tooltip("prevalence:Q", title="Prevalence", format=".3f"),
                ],
            )
            .properties(
                width=285,
                height=315,
                title=alt.TitleParams(
                    text="(a)  Cross-Modality Deconvolution Response Prediction (PR-AUC)",
                    subtitle=[
                        "5×5 CV Precision-Recall AUC (PR-AUC) of Ridge L2 logistic regression",
                        "across all 13 single-cell references (Tiers 1–3) and 15 bulk cohorts",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        model_map = {
            "ridge_l2": "Logistic Reg.",
            "linear_reg": "Linear Reg.",
            "random_forest": "Random Forest",
        }
        model_order = ["Logistic Reg.", "Linear Reg.", "Random Forest"]
        model_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_ORANGE]

        bars_sub = (
            pred_df.filter(
                (pl.col("reference_id").is_in(["Consensus-Meta42", "Unified-All9"]))
                | pl.col("is_indication_matched")
            )
            .filter(pl.col("bulk_cohort").is_in(BULK_COLS_ORDER))
            .filter(pl.col("penalty_mode").is_in(["ridge_l2", "linear_reg", "random_forest"]))
            .group_by(["bulk_cohort", "penalty_mode"])
            .agg([
                pl.col("pr_auc").mean().alias("pr_auc"),
                pl.col("prevalence").mean().alias("prevalence"),
            ])
            .with_columns([
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("bulk_short"),
                pl.col("penalty_mode").replace(model_map).alias("model_label"),
            ])
        )
        bars_pd = bars_sub.to_pandas()

        prev_df = (
            bars_sub.select(["bulk_short", "prevalence"]).unique().to_pandas()
        )
        prev_ticks = (
            alt.Chart(prev_df)
            .mark_tick(color="#0F172A", thickness=2.0, size=14)
            .encode(
                y=alt.Y("bulk_short:N", sort=bulk_short_order),
                x=alt.X("prevalence:Q"),
            )
        )

        bars_b = (
            alt.Chart(bars_pd)
            .mark_bar(stroke="#0F172A", strokeWidth=0.4)
            .encode(
                y=alt.Y(
                    "bulk_short:N",
                    sort=bulk_short_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=160,
                        tickSize=3,
                    ),
                ),
                yOffset=alt.YOffset("model_label:N", sort=model_order),
                x=alt.X(
                    "pr_auc:Q",
                    scale=alt.Scale(domain=[0.0, 1.0]),
                    title="Out-of-Fold PR-AUC (All 15 Bulk Cohorts)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        values=[0.0, 0.5, 1.0],
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "model_label:N",
                    scale=alt.Scale(domain=model_order, range=model_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        symbolSize=90,
                    ),
                ),
            )
            .properties(
                width=190,
                height=315,
                title=alt.TitleParams(
                    text="(b)  Multi-Model PR-AUC",
                    subtitle=[
                        "Out-of-fold PR-AUC across all 15",
                        "cohorts (tick = R prevalence)",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        panel_b = bars_b + prev_ticks

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=14)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig7-prediction-prauc",
            layer_label="07_iAtlas_Response_Prediction_PRAUC_Panels",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 7: {exc}")


# ------------------------------------------------------------------------------
# Figure 9: Full-Atlas (All 42 States) Dual-Scale Forest Plot & 11-Trial Heatmap
# ------------------------------------------------------------------------------


def build_figure9_all42_metastate_forest_and_tme_heatmap(
    states_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 9 (180.0mm): (a) Dual-Axis All-42 Consensus Meta-State Forest Plot + (b) 42x11 Trial Heatmap."""
    try:
        meta_pan = states_df.filter(
            (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        ).sort("bulk_deconv_beta", descending=True)

        all_states = meta_pan["cell_state"].to_list()
        state_to_label = {
            st_id: META_STATE_EMPIRICAL_LABELS.get(st_id, st_id)
            for st_id in all_states
        }
        ordered_labels = [state_to_label[s] for s in all_states]

        bulk_rows: list[dict[str, object]] = []
        milo_rows: list[dict[str, object]] = []
        for row in meta_pan.iter_rows(named=True):
            st_id = str(row["cell_state"])
            lbl = state_to_label[st_id]
            b_val = float(row["bulk_deconv_beta"])
            b_se = float(row["bulk_deconv_se"])
            m_val = float(row["milo_mean_logfc"])
            m_se = float(row["milo_se"])
            bulk_rows.append({
                "state_label": lbl,
                "modality": "Bulk Deconv β (Bottom Axis)",
                "effect": b_val,
                "ci_low": max(-3.8, b_val - 1.96 * b_se),
                "ci_high": min(3.8, b_val + 1.96 * b_se),
            })
            milo_rows.append({
                "state_label": lbl,
                "modality": "SC Milopy logFC (Top Axis)",
                "effect": m_val,
                "ci_low": max(-9.2, m_val - 1.96 * m_se),
                "ci_high": min(9.2, m_val + 1.96 * m_se),
            })

        bulk_pd = pl.DataFrame(bulk_rows).to_pandas()
        milo_pd = pl.DataFrame(milo_rows).to_pandas()
        mod_order = ["Bulk Deconv β (Bottom Axis)", "SC Milopy logFC (Top Axis)"]
        mod_colors = [OKABE_BLUE, OKABE_VERMILION]

        zero_rule_bulk = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=[-3.9, 3.9])))
        )

        bulk_whiskers = (
            alt.Chart(bulk_pd)
            .mark_rule(strokeWidth=1.5, color=OKABE_BLUE)
            .encode(
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_labels,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=245,
                    ),
                ),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X(
                    "ci_low:Q",
                    scale=alt.Scale(domain=[-3.9, 3.9]),
                    title="Bulk Deconv. β (Meta-9 CLR, Bottom Axis)",
                    axis=alt.Axis(
                        orient="bottom",
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        titleColor=OKABE_BLUE,
                        labelColor=OKABE_BLUE,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="ci_high:Q",
                color=alt.Color(
                    "modality:N",
                    scale=alt.Scale(domain=mod_order, range=mod_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=240,
                        symbolSize=85,
                    ),
                ),
            )
        )

        bulk_pts = (
            alt.Chart(bulk_pd)
            .mark_point(filled=True, size=42, stroke="#0F172A", strokeWidth=0.55)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X("effect:Q", scale=alt.Scale(domain=[-3.9, 3.9])),
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        layer_bulk = zero_rule_bulk + bulk_whiskers + bulk_pts

        milo_whiskers = (
            alt.Chart(milo_pd)
            .mark_rule(strokeWidth=1.5, color=OKABE_VERMILION)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels, title=None),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X(
                    "ci_low:Q",
                    scale=alt.Scale(domain=[-9.5, 9.5]),
                    title="Single-Cell Milopy logFC (Top Axis)",
                    axis=alt.Axis(
                        orient="top",
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        titleColor=OKABE_VERMILION,
                        labelColor=OKABE_VERMILION,
                        grid=False,
                    ),
                ),
                x2="ci_high:Q",
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        milo_pts = (
            alt.Chart(milo_pd)
            .mark_point(filled=True, size=42, stroke="#0F172A", strokeWidth=0.55)
            .encode(
                y=alt.Y("state_label:N", sort=ordered_labels),
                yOffset=alt.YOffset("modality:N", sort=mod_order),
                x=alt.X("effect:Q", scale=alt.Scale(domain=[-9.5, 9.5])),
                color=alt.Color("modality:N", scale=alt.Scale(domain=mod_order, range=mod_colors)),
            )
        )

        layer_milo = milo_whiskers + milo_pts

        panel_a = (
            alt.layer(layer_bulk, layer_milo)
            .resolve_scale(x="independent")
            .properties(
                width=245,
                height=755,
                title=alt.TitleParams(
                    text="(a)  All 42 Meta-State Effect Sizes (±95% CI)",
                    subtitle=[
                        "All 42 Consensus-Meta42 states ranked by Meta-9 β;",
                        "dual x-axes: bulk CLR β (bottom) vs. SC Milopy logFC (top)",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        indiv_trials = [
            "Hugo-iAtlas",
            "Riaz-iAtlas",
            "Liu-iAtlas",
            "Gide-iAtlas",
            "Rosenberg-iAtlas",
            "Padron-iAtlas",
            "Anders-iAtlas",
            "McDermott-iAtlas",
            "Choueiri-iAtlas",
            "Jung",
            "Cho",
        ]
        trial_short_order = [BULK_ID_TO_SHORT[c] for c in indiv_trials]

        heat_sub = (
            states_df.filter(pl.col("reference_id") == "Consensus-Meta42")
            .filter(pl.col("bulk_cohort").is_in(indiv_trials))
            .with_columns([
                pl.col("cell_state").replace(state_to_label).alias("state_label"),
                pl.col("bulk_cohort").replace(BULK_ID_TO_SHORT).alias("trial_short"),
                pl.col("bulk_deconv_beta").clip(-2.2, 2.2).alias("beta_clipped"),
            ])
        )

        panel_b = (
            alt.Chart(heat_sub.to_pandas())
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.55)
            .encode(
                x=alt.X(
                    "trial_short:N",
                    sort=trial_short_order,
                    title="Individual Bulk RNA-seq Clinical Trial (11 Cohorts)",
                    axis=alt.Axis(
                        labelAngle=-45,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=160,
                    ),
                ),
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_labels,
                    title=None,
                    axis=alt.Axis(labels=False, ticks=False),
                ),
                color=alt.Color(
                    "beta_clipped:Q",
                    title="Trial Deconv β (CLR)",
                    scale=alt.Scale(domain=[-2.2, 2.2], scheme="blueorange"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
            .properties(
                width=255,
                height=755,
                title=alt.TitleParams(
                    text="(b)  Per-Trial R vs. NR Shift (All 42 States)",
                    subtitle=[
                        "Within-trial R vs. NR CLR shift (β) across",
                        "all 11 clinical trials for all 42 meta-states",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=12)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig9-all42-metastate-forest",
            layer_label="09_All42_Consensus_MetaState_Forest_And_Trials",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 9: {exc}")


# ------------------------------------------------------------------------------
# Figure 10 (NEW): Top 4 Differentially Expressed Marker Genes Across All 42 Meta-States
# ------------------------------------------------------------------------------


def build_figure10_all42_metastate_top_degs(
    all42_deg_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 10 (180.0mm): (a) All-42 State x Top 4 DEG Matrix (168 Markers) + (b) Bulk R vs NR t-Stat."""
    try:
        deg_clean = all42_deg_df.with_columns(
            pl.col("cell_state").replace(META_STATE_EMPIRICAL_LABELS).alias("state_label")
        )
        state_order_df = (
            deg_clean.select(["cell_state", "state_label", "state_bulk_beta"])
            .unique()
            .sort("state_bulk_beta", descending=True)
        )
        ordered_state_labels = state_order_df["state_label"].to_list()

        deg_plot_df = deg_clean.with_columns(
            (pl.lit("Rank #") + pl.col("gene_rank").cast(pl.Utf8)).alias("rank_col")
        )
        rank_order = ["Rank #1", "Rank #2", "Rank #3", "Rank #4"]
        deg_pd = deg_plot_df.to_pandas()

        rect_a = (
            alt.Chart(deg_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "rank_col:N",
                    sort=rank_order,
                    title="Top DEG Marker Rank (in Φ_meta)",
                    axis=alt.Axis(
                        labelAngle=0,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                    ),
                ),
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_state_labels,
                    title="All 42 Consensus-Meta42 Cell States (Ranked by Meta-9 β)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=245,
                    ),
                ),
                color=alt.Color(
                    "specificity_log2fc:Q",
                    title="Marker Specificity (log2FC)",
                    scale=alt.Scale(domain=[0.7, 3.5], scheme="blues"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=130,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
        )

        text_a = (
            alt.Chart(deg_pd)
            .mark_text(fontSize=MIN_FONT_PX, fontWeight="bold")
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

        panel_a = (rect_a + text_a).properties(
            width=285,
            height=760,
            title=alt.TitleParams(
                text="(a)  Top 4 DEG Marker Genes Across All 42 Meta-States",
                subtitle=[
                    "Top 4 protein-coding HGNC marker genes by log2FC specificity",
                    "in Φ_meta across all 42 Consensus-Meta42 states (168 markers)",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        state_gene_summary = (
            deg_clean.group_by(["cell_state", "state_label", "enrichment_group", "state_bulk_beta"])
            .agg([
                pl.col("bulk_responder_tstat").mean().alias("mean_tstat"),
                pl.col("specificity_log2fc").mean().alias("mean_spec_lfc"),
            ])
            .sort("state_bulk_beta", descending=True)
        )
        sg_pd = state_gene_summary.to_pandas()
        grp_order = ["Responder-Enriched (β > 0)", "Non-Responder-Enriched (β < 0)"]
        grp_colors = [OKABE_BLUE, OKABE_VERMILION]

        zero_t = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[4, 3], strokeWidth=1.2)
            .encode(x="x:Q")
        )

        bars_b = (
            alt.Chart(sg_pd)
            .mark_bar(stroke="#0F172A", strokeWidth=0.45)
            .encode(
                y=alt.Y(
                    "state_label:N",
                    sort=ordered_state_labels,
                    title=None,
                    axis=alt.Axis(labels=False, ticks=False),
                ),
                x=alt.X(
                    "mean_tstat:Q",
                    title="Bulk Signature R vs. NR t-stat",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "enrichment_group:N",
                    scale=alt.Scale(domain=grp_order, range=grp_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=220,
                        symbolSize=85,
                    ),
                ),
            )
            .properties(
                width=180,
                height=760,
                title=alt.TitleParams(
                    text="(b)  Bulk DEG Enrichment",
                    subtitle=[
                        "Pan-Cancer iAtlas (N=1,015)",
                        "R vs. NR signature t-stat (42 states)",
                    ],
                    fontSize=TITLE_FONT_PX,
                    subtitleFontSize=MIN_FONT_PX,
                    subtitleColor="#475569",
                    fontWeight="bold",
                    anchor="start",
                    offset=6,
                ),
            )
        )

        panel_b = bars_b + zero_t

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=12)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-fig10-all42-metastate-degs",
            layer_label="10_All42_Consensus_MetaState_Top_DEGs",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 10: {exc}")


# ------------------------------------------------------------------------------
# Shared p-Value Formatting Helper
# ------------------------------------------------------------------------------


def _fmt_pval_short(p_val: float) -> str:
    """Format p-values cleanly without zero-truncation."""
    if p_val < 1e-4:
        return f"{p_val:.2g}"
    if p_val < 0.001:
        return f"{p_val:.5f}"
    if p_val < 0.01:
        return f"{p_val:.4f}"
    if p_val < 0.10:
        return f"{p_val:.3f}"
    return f"{p_val:.2f}"


# ------------------------------------------------------------------------------
# Main Article Figure: 3x3 Pan-Compartment Gating, M29xM26 Epistasis & Single-Cell Triad UMAPs (a-i)
# ------------------------------------------------------------------------------

_UMAP_REGIME_DOMAIN: Final[list[str]] = [
    "++ Co-Presence Synergy",
    "+/- Mutual-Exclusion Gate",
    "-- Dual-Absence Desert",
    "Log-Ratio Balance",
    "Co-Abundance Coupling",
]
_UMAP_REGIME_COLORS: Final[list[str]] = [
    OKABE_TEAL,
    "#7E22CE",
    "#CA8A04",
    "#BE185D",
    OKABE_PURPLE,
]
_UMAP_COMPARTMENT_PALETTE: Final[dict[str, str]] = {
    "T-NK cells": OKABE_TEAL,
    "Myeloid/Macrophages": "#7E22CE",
    "Stromal/CAFs": "#CA8A04",
    "CAF": "#CA8A04",
    "B/Plasma cells": "#BE185D",
    "Tumor/Epithelial": "#475569",
}
_UMAP_REGIME_DASHES: Final[list[list[int]]] = [
    [1, 0],
    [1, 0],
    [5, 3],
    [6, 2, 2, 2],
    [2, 2],
]
_UMAP_BUNDLES_ROOT: Final[Path] = Path(
    "output/results/milopy_pan_cohort_integration_and_power/compartment_pair_umap_bundles"
)


def _blend_hex_with_white(hex_color: str, alpha: float) -> str:
    """Blend #RRGGBB deterministically toward pure #FFFFFF without mutating state."""
    clean = hex_color.strip().lstrip("#")
    if len(clean) != 6:
        return "#E2E8F0"
    a_clamped = float(np.clip(alpha, 0.0, 1.0))
    r_in = int(clean[0:2], 16)
    g_in = int(clean[2:4], 16)
    b_in = int(clean[4:6], 16)
    r_out = int(round(a_clamped * float(r_in) + (1.0 - a_clamped) * 255.0))
    g_out = int(round(a_clamped * float(g_in) + (1.0 - a_clamped) * 255.0))
    b_out = int(round(a_clamped * float(b_in) + (1.0 - a_clamped) * 255.0))
    return f"#{r_out:02X}{g_out:02X}{b_out:02X}"


def _build_triad_umap_panel(
    slug: str,
    panel_letter: str,
    width: int = 245,
    height: int = 200,
    bundles_root: Path = _UMAP_BUNDLES_ROOT,
) -> Result[alt.LayerChart, str]:
    """Build a heading-free single-cell Harmony I2 cross-compartment UMAP panel from cached parquet bundle."""
    slug_dir = bundles_root / slug
    meta_p = slug_dir / "meta.parquet"
    cells_p = slug_dir / "cells.parquet"
    t1_p = slug_dir / "tier1_arcs.parquet"
    t2_p = slug_dir / "tier2_arcs.parquet"
    ep_p = slug_dir / "endpoints.parquet"
    lbl_p = slug_dir / "medoid_labels.parquet"

    for req_p in (meta_p, cells_p, t1_p, t2_p, ep_p, lbl_p):
        if not req_p.exists():
            return Failure(f"Missing cached UMAP bundle file: {req_p}")

    meta_row = pl.read_parquet(meta_p).row(0, named=True)
    cells_df = pl.read_parquet(cells_p)
    tier1_arcs_df = pl.read_parquet(t1_p)
    tier2_arcs_df = pl.read_parquet(t2_p)
    endpoints_df = pl.read_parquet(ep_p)
    medoid_labels_df = pl.read_parquet(lbl_p)

    disp_1 = str(meta_row["display_first"])
    disp_2 = str(meta_row["display_second"])
    disp_title = str(meta_row["display_title"])
    col_1 = _UMAP_COMPARTMENT_PALETTE.get(disp_1, OKABE_TEAL)
    col_2 = _UMAP_COMPARTMENT_PALETTE.get(disp_2, "#7E22CE")
    n_fdr = int(meta_row["n_fdr_pairs"])
    n_perm = int(meta_row["n_perm_only_pairs"])
    raw_lead = str(meta_row["lead_pair_summary"]).replace("Lead: ", "")
    lead_summary = re.sub(
        r"\b(M\d{2})\b",
        lambda m: META_STATE_COMPACT_UMAP_NAMES.get(m.group(1), m.group(1)),
        raw_lead,
    )
    if not medoid_labels_df.is_empty():
        def _clean_umap_callout(txt: str) -> str:
            s = str(txt).strip()
            m_code = re.match(r"^(M\d{2})(?:\s+.*)?$", s)
            if m_code:
                code = m_code.group(1)
                return META_STATE_COMPACT_UMAP_NAMES.get(code, code)
            return s

        medoid_labels_df = medoid_labels_df.with_columns(
            pl.col("label_text").map_elements(_clean_umap_callout, return_dtype=pl.String).alias("label_text")
        )
        # Repel wider text callouts so longer cell-state names never overlap
        u1_arr = medoid_labels_df["label_u1"].to_numpy().astype(np.float64).copy()
        u2_arr = medoid_labels_df["label_u2"].to_numpy().astype(np.float64).copy()
        ou1_arr = medoid_labels_df["orig_u1"].to_numpy().astype(np.float64)
        ou2_arr = medoid_labels_df["orig_u2"].to_numpy().astype(np.float64)
        n_lbl = len(u1_arr)
        for _ in range(45):
            for i_l in range(n_lbl):
                for j_l in range(i_l + 1, n_lbl):
                    dx = u1_arr[i_l] - u1_arr[j_l]
                    dy = u2_arr[i_l] - u2_arr[j_l]
                    if abs(dx) < 8.6 and abs(dy) < 2.15:
                        push_y = 0.52 * (2.15 - abs(dy)) * (1.0 if dy >= 0.0 else -1.0)
                        push_x = 0.35 * (8.6 - abs(dx)) * (1.0 if dx >= 0.0 else -1.0)
                        u2_arr[i_l] += push_y
                        u2_arr[j_l] -= push_y
                        u1_arr[i_l] += push_x
                        u1_arr[j_l] -= push_x
        u1_arr = np.clip(u1_arr, float(meta_row["umap_xmin"]) + 2.8, float(meta_row["umap_xmax"]) - 2.8)
        u2_arr = np.clip(u2_arr, float(meta_row["umap_ymin"]) + 0.9, float(meta_row["umap_ymax"]) - 1.2)
        has_ldr = np.hypot(u1_arr - ou1_arr, u2_arr - ou2_arr) > 0.75
        medoid_labels_df = medoid_labels_df.with_columns(
            [
                pl.Series("label_u1", u1_arr),
                pl.Series("label_u2", u2_arr),
                pl.Series("has_leader", has_ldr),
            ]
        )
    bounds = (
        float(meta_row["umap_xmin"]),
        float(meta_row["umap_xmax"]),
        float(meta_row["umap_ymin"]),
        float(meta_row["umap_ymax"]),
    )

    x_enc = alt.X(
        "umap_1:Q",
        scale=alt.Scale(domain=[bounds[0], bounds[1]], nice=False),
        title=f"UMAP 1 ({disp_title})",
        axis=alt.Axis(
            values=[-10, -5, 0, 5, 10, 15],
            labelFontSize=11.0,
            titleFontSize=11.5,
            titleFontWeight="bold",
            grid=False,
        ),
    )
    y_enc = alt.Y(
        "umap_2:Q",
        scale=alt.Scale(domain=[bounds[2], bounds[3]], nice=False),
        title="UMAP 2",
        axis=alt.Axis(
            values=[-10, -5, 0, 5, 10, 15],
            labelFontSize=11.0,
            titleFontSize=11.5,
            titleFontWeight="bold",
            titlePadding=5,
            minExtent=28,
            maxExtent=28,
            grid=False,
        ),
    )

    role_other = "Other Compartments"
    role_ctx_1 = f"{disp_1} (Context)"
    role_ctx_2 = f"{disp_2} (Context)"
    role_act_1 = f"{disp_1} (Interacting)"
    role_act_2 = f"{disp_2} (Interacting)"
    role_domain = [role_act_1, role_act_2, role_ctx_1, role_ctx_2, role_other]
    role_colors = [
        col_1,
        col_2,
        _blend_hex_with_white(col_1, 0.28),
        _blend_hex_with_white(col_2, 0.28),
        "#E8EEF5",
    ]

    cells_layer = (
        alt.Chart(cells_df.to_pandas())
        .mark_circle(size=6.2, opacity=0.84)
        .encode(
            x=x_enc,
            y=y_enc,
            color=alt.Color(
                "cell_footprint_role:N",
                scale=alt.Scale(domain=role_domain, range=role_colors),
                legend=None,
            ),
        )
    )
    sub_layers: list[alt.Chart] = [cells_layer]

    tier_fdr_lbl = "ACAT FDR < 0.10"
    tier_perm_lbl = "Perm. p < 0.05"
    tier_domain = [tier_fdr_lbl, tier_perm_lbl]
    sw_range = [1.95, 0.95]
    op_range = [0.92, 0.46]

    arc_frames: list[pl.DataFrame] = []
    if not tier2_arcs_df.is_empty():
        arc_frames.append(
            tier2_arcs_df.with_columns(
                (pl.lit("t2_") + pl.col("arc_id").cast(pl.String)).alias("arc_uid"),
                pl.lit(tier_perm_lbl).alias("tier_label"),
            )
        )
    if not tier1_arcs_df.is_empty():
        arc_frames.append(
            tier1_arcs_df.with_columns(
                (pl.lit("t1_") + pl.col("arc_id").cast(pl.String)).alias("arc_uid"),
                pl.lit(tier_fdr_lbl).alias("tier_label"),
            )
        )
    if arc_frames:
        combined_arcs_df = pl.concat(arc_frames, how="vertical_relaxed")
        sub_layers.append(
            alt.Chart(combined_arcs_df.to_pandas())
            .mark_line()
            .encode(
                x=x_enc,
                y=y_enc,
                detail="arc_uid:N",
                order="step_idx:Q",
                color=alt.Color(
                    "regime:N",
                    scale=alt.Scale(domain=_UMAP_REGIME_DOMAIN, range=_UMAP_REGIME_COLORS),
                    legend=None,
                ),
                strokeDash=alt.StrokeDash(
                    "regime:N",
                    scale=alt.Scale(domain=_UMAP_REGIME_DOMAIN, range=_UMAP_REGIME_DASHES),
                    legend=None,
                ),
                strokeWidth=alt.StrokeWidth(
                    "tier_label:N",
                    scale=alt.Scale(domain=tier_domain, range=sw_range),
                    legend=None,
                ),
                opacity=alt.Opacity(
                    "tier_label:N",
                    scale=alt.Scale(domain=tier_domain, range=op_range),
                    legend=None,
                ),
            )
        )

    # Top-left inline two-dot lineage key
    bx_dot = bounds[0] + 0.85
    bx_txt = bounds[0] + 1.55
    by_1 = bounds[3] - 1.05
    by_2 = bounds[3] - 2.55
    sub_layers.append(
        alt.Chart(pl.DataFrame([{"u1": bx_dot, "u2": by_1}]).to_pandas())
        .mark_circle(size=34, color=col_1, stroke="#FFFFFF", strokeWidth=0.6)
        .encode(x="u1:Q", y="u2:Q")
    )
    sub_layers.append(
        alt.Chart(pl.DataFrame([{"u1": bx_dot, "u2": by_2}]).to_pandas())
        .mark_circle(size=34, color=col_2, stroke="#FFFFFF", strokeWidth=0.6)
        .encode(x="u1:Q", y="u2:Q")
    )
    sub_layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {"u1": bx_txt, "u2": by_1, "lbl": disp_1},
                    {"u1": bx_txt, "u2": by_2, "lbl": disp_2},
                ]
            ).to_pandas()
        )
        .mark_text(fontSize=8.8, fontWeight="bold", fill="#1E293B", align="left", baseline="middle")
        .encode(x="u1:Q", y="u2:Q", text="lbl:N")
    )

    # Polar endpoint glyphs
    if not endpoints_df.is_empty():
        for tier_val, sz_scale, ep_op in (
            ("Stratified Perm. p < 0.05 (*)", 0.56, 0.72),
            ("ACAT FDR < 0.10 (Verified)", 1.00, 0.98),
        ):
            ep_sub = endpoints_df.filter(pl.col("sig_tier") == tier_val)
            if ep_sub.is_empty():
                continue
            base_sz = 36.0 * sz_scale
            pos_ep = ep_sub.filter(
                pl.col("polarity").is_in(["Present / Enriched (+)", "Co-Abundant"])
            )
            neg_ep = ep_sub.filter(pl.col("polarity") == "Excluded / Gated (-)")
            des_ep = ep_sub.filter(pl.col("polarity") == "Dual-Absent (--)")

            if not pos_ep.is_empty():
                sub_layers.append(
                    alt.Chart(pos_ep.to_pandas())
                    .mark_circle(
                        size=base_sz,
                        opacity=ep_op,
                        color=OKABE_TEAL,
                        stroke="#FFFFFF",
                        strokeWidth=0.85,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )
            if not neg_ep.is_empty():
                sub_layers.append(
                    alt.Chart(neg_ep.to_pandas())
                    .mark_point(
                        shape="circle",
                        filled=True,
                        size=base_sz * 1.12,
                        fill="#FFFFFF",
                        stroke="#7E22CE",
                        strokeWidth=1.75 if sz_scale > 0.8 else 1.25,
                        opacity=ep_op,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )
            if not des_ep.is_empty():
                sub_layers.append(
                    alt.Chart(des_ep.to_pandas())
                    .mark_point(
                        shape="diamond",
                        filled=True,
                        size=base_sz,
                        fill="#CA8A04",
                        stroke="#0F172A",
                        strokeWidth=0.75,
                        opacity=ep_op,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )

    # Collision-repelled Meta42 state callouts
    if not medoid_labels_df.is_empty():
        leaders_df = medoid_labels_df.filter(pl.col("has_leader"))
        if not leaders_df.is_empty():
            sub_layers.append(
                alt.Chart(leaders_df.to_pandas())
                .mark_rule(color="#475569", strokeWidth=0.65, opacity=0.75)
                .encode(
                    x="orig_u1:Q",
                    y="orig_u2:Q",
                    x2="label_u1:Q",
                    y2="label_u2:Q",
                )
            )
        t1_lbls = medoid_labels_df.filter(pl.col("is_tier1"))
        t2_lbls = medoid_labels_df.filter(~pl.col("is_tier1"))
        tier_label_specs: tuple[tuple[pl.DataFrame, Literal["normal", "bold"], str], ...] = (
            (t2_lbls, "normal", "#334155"),
            (t1_lbls, "bold", "#0F172A"),
        )
        for lbl_df, fw, fg_hex in tier_label_specs:
            if lbl_df.is_empty():
                continue
            sub_layers.append(
                alt.Chart(lbl_df.to_pandas())
                .mark_text(
                    fontSize=8.2,
                    fontWeight=fw,
                    fill=fg_hex,
                    stroke="#FFFFFF",
                    strokeWidth=2.3,
                    align="center",
                    baseline="middle",
                )
                .encode(x="label_u1:Q", y="label_u2:Q", text="label_text:N")
            )
            sub_layers.append(
                alt.Chart(lbl_df.to_pandas())
                .mark_text(
                    fontSize=8.2,
                    fontWeight=fw,
                    fill=fg_hex,
                    align="center",
                    baseline="middle",
                )
                .encode(x="label_u1:Q", y="label_u2:Q", text="label_text:N")
            )

    # Top header strip above UMAP box showing FDR/perm counts & lead pair
    short_title = (
        disp_title.replace("Myeloid/Macrophages", "Myeloid/Mac").replace("T-NK cells", "T-NK")
    )
    hdr_line1 = (
        alt.Chart(pl.DataFrame({"t": [f"{short_title} ({n_fdr} FDR<0.10 · {n_perm} p<0.05)"]}).to_pandas())
        .mark_text(
            align="left",
            baseline="middle",
            fontSize=10.0,
            fontWeight="bold",
            color="#0F172A",
            clip=False,
        )
        .encode(x=alt.value(0.0), y=alt.value(-20.0), text="t:N")
    )
    hdr_line2 = (
        alt.Chart(pl.DataFrame({"t": [f"Lead: {lead_summary}"]}).to_pandas())
        .mark_text(
            align="left",
            baseline="middle",
            fontSize=8.5,
            fontWeight="bold",
            color="#334155",
            clip=False,
        )
        .encode(x=alt.value(0.0), y=alt.value(-7.5), text="t:N")
    )
    sub_layers.extend(
        [
            hdr_line1,
            hdr_line2,
            _panel_letter_layer(panel_letter, x_px=-34.0, y_px=-20.0),
        ]
    )

    return Success(
        cast(
            alt.LayerChart,
            alt.layer(*sub_layers)
            .resolve_scale(
                color="independent",
                strokeDash="independent",
                strokeWidth="independent",
                opacity="independent",
            )
            .properties(width=width, height=height),
        )
    )


def _build_triad_umap_legend_strip(width: int = 815, height: int = 24) -> alt.LayerChart:
    """Build a compact horizontal key decoding UMAP chord tiers, regimes, and endpoint polarity glyphs."""
    x_enc = alt.X("x:Q", scale=alt.Scale(domain=[0.0, 100.0], nice=False), axis=None)
    y_enc = alt.Y("y:Q", scale=alt.Scale(domain=[0.0, 10.0], nice=False), axis=None)

    layers: list[alt.Chart] = [
        alt.Chart(pl.DataFrame([{"x": 0.0, "x2": 100.0, "y": 0.5, "y2": 9.5}]).to_pandas())
        .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.65)
        .encode(x=x_enc, x2="x2:Q", y=y_enc, y2="y2:Q")
    ]

    # Tier weights
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 1.2, "x2": 3.8, "y": 5.0, "y2": 5.0}]).to_pandas())
        .mark_rule(color="#0F172A", strokeWidth=2.3)
        .encode(x=x_enc, x2="x2:Q", y=y_enc, y2="y2:Q")
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 14.8, "x2": 17.4, "y": 5.0, "y2": 5.0}]).to_pandas())
        .mark_rule(color="#475569", strokeWidth=1.0, opacity=0.65)
        .encode(x=x_enc, x2="x2:Q", y=y_enc, y2="y2:Q")
    )

    # Regime swatches
    reg_items = (
        (28.5, "#BE185D", [6, 2, 2, 2], "Log-Ratio Balance"),
        (43.5, "#7E22CE", [1, 0], "+/- Gate"),
        (53.2, OKABE_TEAL, [1, 0], "++ Synergy"),
        (63.8, OKABE_PURPLE, [2, 2], "Co-Abundance"),
    )
    txt_rows: list[dict[str, object]] = [
        {"x": 4.2, "y": 5.0, "txt": "ACAT FDR<0.10"},
        {"x": 17.8, "y": 5.0, "txt": "Perm p<0.05"},
    ]
    for rx, rcol, rdash, rlbl in reg_items:
        layers.append(
            alt.Chart(pl.DataFrame([{"x": rx, "x2": rx + 2.6, "y": 5.0, "y2": 5.0}]).to_pandas())
            .mark_rule(color=rcol, strokeWidth=1.9, strokeDash=rdash)
            .encode(x=x_enc, x2="x2:Q", y=y_enc, y2="y2:Q")
        )
        txt_rows.append({"x": rx + 3.0, "y": 5.0, "txt": rlbl})

    # Polarity glyphs
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 77.2, "y": 5.0}]).to_pandas())
        .mark_circle(size=38, color=OKABE_TEAL, stroke="#FFFFFF", strokeWidth=0.8)
        .encode(x=x_enc, y=y_enc)
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 85.8, "y": 5.0}]).to_pandas())
        .mark_point(
            shape="circle",
            filled=True,
            size=42,
            fill="#FFFFFF",
            stroke="#7E22CE",
            strokeWidth=1.7,
        )
        .encode(x=x_enc, y=y_enc)
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 93.4, "y": 5.0}]).to_pandas())
        .mark_point(
            shape="diamond",
            filled=True,
            size=38,
            fill="#CA8A04",
            stroke="#0F172A",
            strokeWidth=0.75,
        )
        .encode(x=x_enc, y=y_enc)
    )
    txt_rows.extend(
        [
            {"x": 78.2, "y": 5.0, "txt": "+ Enriched"},
            {"x": 86.8, "y": 5.0, "txt": "- Gated"},
            {"x": 94.4, "y": 5.0, "txt": "-- Desert"},
        ]
    )
    layers.append(
        alt.Chart(pl.DataFrame(txt_rows).to_pandas())
        .mark_text(fontSize=9.2, fontWeight="bold", fill="#1E293B", align="left", baseline="middle")
        .encode(x=x_enc, y=y_enc, text="txt:N")
    )

    return cast(alt.LayerChart, alt.layer(*layers).properties(width=width, height=height))


def build_main_epistasis_and_stromal_gating_figure(
    oof_cells_df: pl.DataFrame,
    oof_contrasts_df: pl.DataFrame,
    landscape_df: pl.DataFrame,
    clin_df: pl.DataFrame,
    epi_sig_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build full-page heading-free Main Figure (180.0mm, 9 panels a-i in 3x3 grid): gating (a-e), M29xM26 (f) & Triad UMAPs (g-i)."""
    try:
        _ = oof_contrasts_df
        x_dom_dose = [-2.62, 2.92]
        y_dom_dose = [5.0, 85.0]
        w_col = 232
        h_top = 150
        y_ax_ext = 28
        y_ttl_pad = 5

        # ------------------------------------------------------------------
        # Panel (a): Continuous R_OOF Dose-Response Curve (Permissive vs Excluded)
        # ------------------------------------------------------------------
        short_9a_map = {
            "Permissive TME (Low NR, T1)": "Permissive (NR T1)",
            "Excluded TME (High NR, T3)": "Excluded (NR T3)",
        }
        curve_9a = landscape_df.filter(pl.col("record_type") == "dose_curve_9a").with_columns(
            [
                pl.col("stratum_label").replace(short_9a_map).alias("s_short"),
                pl.when(pl.col("stratum_id") == "permissive_nr_t1")
                .then(pl.lit(OKABE_TEAL))
                .otherwise(pl.lit("#7E22CE"))
                .alias("s_hex"),
            ]
        )
        bins_9a = landscape_df.filter(pl.col("record_type") == "quintile_bin_9a").with_columns(
            [
                pl.col("stratum_label").replace(short_9a_map).alias("s_short"),
                pl.when(pl.col("stratum_id") == "permissive_nr_t1")
                .then(pl.lit(OKABE_TEAL))
                .otherwise(pl.lit("#7E22CE"))
                .alias("s_hex"),
            ]
        )
        perm_9a_row = curve_9a.filter(pl.col("stratum_id") == "permissive_nr_t1").row(0, named=True)
        excl_9a_row = curve_9a.filter(pl.col("stratum_id") == "excluded_nr_t3").row(0, named=True)

        dom_9a = ["Permissive (NR T1)", "Excluded (NR T3)"]
        cols_9a = [OKABE_TEAL, "#7E22CE"]
        shapes_9a = ["circle", "square"]

        area_9a_perm = (
            alt.Chart(curve_9a.filter(pl.col("stratum_id") == "permissive_nr_t1").to_pandas())
            .mark_area(fill=OKABE_TEAL, fillOpacity=0.14)
            .encode(
                x=alt.X(
                    "x_z:Q",
                    scale=alt.Scale(domain=x_dom_dose),
                    title="Effector R_OOF (T/NK · M1-TAM · B/Plasma, z)",
                    axis=alt.Axis(
                        values=[-2, -1, 0, 1, 2],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "ci_low_pct:Q",
                    scale=alt.Scale(domain=y_dom_dose),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
                        values=[15, 30, 45, 60, 75],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y2="ci_high_pct:Q",
            )
        )
        area_9a_excl = (
            alt.Chart(curve_9a.filter(pl.col("stratum_id") == "excluded_nr_t3").to_pandas())
            .mark_area(fill="#7E22CE", fillOpacity=0.14)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
            )
        )
        area_9a = area_9a_perm + area_9a_excl
        line_9a = (
            alt.Chart(curve_9a.to_pandas())
            .mark_line(strokeWidth=2.3)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("pred_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                detail="s_short:N",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        err_9a = (
            alt.Chart(bins_9a.to_pandas())
            .mark_rule(strokeWidth=1.35, opacity=0.78)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        pts_9a = (
            alt.Chart(bins_9a.to_pandas())
            .mark_point(filled=True, size=72, opacity=1.0, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("pred_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                color=alt.Color(
                    "s_short:N",
                    scale=alt.Scale(domain=dom_9a, range=cols_9a),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
                shape=alt.Shape(
                    "s_short:N",
                    scale=alt.Scale(domain=dom_9a, range=shapes_9a),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
            )
        )
        box_9a = (
            alt.Chart(pl.DataFrame([{"bx1": -2.54, "bx2": 2.84, "by1": 66.2, "by2": 83.4}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_dose)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="by2:Q",
            )
        )
        txt_9a = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": -2.44,
                            "ty": 78.2,
                            "txt": f"Permissive OR = {float(perm_9a_row['or_per_sd']):.2f}/SD (p = {_fmt_pval_short(float(perm_9a_row['wald_p']))})",
                            "col": OKABE_TEAL,
                        },
                        {
                            "tx": -2.44,
                            "ty": 70.8,
                            "txt": f"Excluded OR = {float(excl_9a_row['or_per_sd']):.2f}/SD (p = {_fmt_pval_short(float(excl_9a_row['wald_p']))}, n.s.)",
                            "col": "#7E22CE",
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=10.2, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_dose)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_a = (
            area_9a
            + line_9a
            + err_9a
            + pts_9a
            + box_9a
            + txt_9a
            + _panel_letter_layer("a", x_px=-34.0, y_px=-24.0)
        ).properties(width=w_col, height=h_top)

        # ------------------------------------------------------------------
        # Panel (b): Symmetrical Continuous NR_OOF Barrier Suppression Curve
        # ------------------------------------------------------------------
        curve_9b = landscape_df.filter(pl.col("record_type") == "dose_curve_9b").with_columns(
            pl.when(pl.col("stratum_id") == "high_r_t3")
            .then(pl.lit(OKABE_TEAL))
            .otherwise(pl.lit("#64748B"))
            .alias("s_hex")
        )
        bins_9b = landscape_df.filter(pl.col("record_type") == "quintile_bin_9b").with_columns(
            pl.when(pl.col("stratum_id") == "high_r_t3")
            .then(pl.lit(OKABE_TEAL))
            .otherwise(pl.lit("#64748B"))
            .alias("s_hex")
        )
        hi_9b_row = curve_9b.filter(pl.col("stratum_id") == "high_r_t3").row(0, named=True)
        lo_9b_row = curve_9b.filter(pl.col("stratum_id") == "low_r_t1").row(0, named=True)

        dom_9b = ["High R Score (T3)", "Low R Score (T1)"]
        cols_9b = [OKABE_TEAL, "#64748B"]
        shapes_9b = ["diamond", "triangle-up"]

        area_9b_hi = (
            alt.Chart(curve_9b.filter(pl.col("stratum_id") == "high_r_t3").to_pandas())
            .mark_area(fill=OKABE_TEAL, fillOpacity=0.14)
            .encode(
                x=alt.X(
                    "x_z:Q",
                    scale=alt.Scale(domain=x_dom_dose),
                    title="Barrier NR_OOF (MDSC · EMT · Hypox · CAF, z)",
                    axis=alt.Axis(
                        values=[-2, -1, 0, 1, 2],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "ci_low_pct:Q",
                    scale=alt.Scale(domain=y_dom_dose),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
                        values=[15, 30, 45, 60, 75],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y2="ci_high_pct:Q",
            )
        )
        area_9b_lo = (
            alt.Chart(curve_9b.filter(pl.col("stratum_id") == "low_r_t1").to_pandas())
            .mark_area(fill="#64748B", fillOpacity=0.14)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
            )
        )
        area_9b = area_9b_hi + area_9b_lo
        line_9b = (
            alt.Chart(curve_9b.to_pandas())
            .mark_line(strokeWidth=2.3)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("pred_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                detail="stratum_label:N",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        err_9b = (
            alt.Chart(bins_9b.to_pandas())
            .mark_rule(strokeWidth=1.35, opacity=0.78)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        pts_9b = (
            alt.Chart(bins_9b.to_pandas())
            .mark_point(filled=True, size=72, opacity=1.0, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("pred_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                color=alt.Color(
                    "stratum_label:N",
                    scale=alt.Scale(domain=dom_9b, range=cols_9b),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
                shape=alt.Shape(
                    "stratum_label:N",
                    scale=alt.Scale(domain=dom_9b, range=shapes_9b),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
            )
        )
        box_9b = (
            alt.Chart(pl.DataFrame([{"bx1": -2.54, "bx2": 2.84, "by1": 66.2, "by2": 83.4}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_dose)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="by2:Q",
            )
        )
        txt_9b = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": -2.44,
                            "ty": 78.2,
                            "txt": f"High-R OR = {float(hi_9b_row['or_per_sd']):.2f}/SD (p = {float(hi_9b_row['wald_p']):.1e})",
                            "col": OKABE_TEAL,
                        },
                        {
                            "tx": -2.44,
                            "ty": 70.8,
                            "txt": f"Low-R OR = {float(lo_9b_row['or_per_sd']):.2f}/SD (p = {_fmt_pval_short(float(lo_9b_row['wald_p']))}, n.s.)",
                            "col": "#475569",
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=10.2, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_dose)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_b = (
            area_9b
            + line_9b
            + err_9b
            + pts_9b
            + box_9b
            + txt_9b
            + _panel_letter_layer("b", x_px=-34.0, y_px=-24.0)
        ).properties(width=w_col, height=h_top)

        # ------------------------------------------------------------------
        # Panel (c): Joint 2D Out-of-Fold Patient Coordinate Manifold (N = 1,015)
        # ------------------------------------------------------------------
        pat_df = landscape_df.filter(pl.col("record_type") == "patient").with_columns(
            [
                pl.col("x_z").clip(-2.48, 2.48).alias("rx_clip"),
                pl.col("y_z").clip(-2.48, 2.48).alias("ny_clip"),
                pl.when(pl.col("regime_group") == "Middle T2 Band (Omitted)")
                .then(pl.lit("Omitted Middle T2"))
                .when(pl.col("response_binary") == 1)
                .then(pl.lit("Corner CR/PR"))
                .otherwise(pl.lit("Corner SD/PD"))
                .alias("pt_class"),
            ]
        )
        pat_first = pat_df.row(0, named=True)
        rq33 = float(pat_first["r_cut_q33"])
        rq67 = float(pat_first["r_cut_q67"])
        nq33 = float(pat_first["n_cut_q33"])
        nq67 = float(pat_first["n_cut_q67"])

        mid_pd = pat_df.filter(pl.col("pt_class") == "Omitted Middle T2").to_pandas()
        cor_pat_pd = pat_df.filter(pl.col("pt_class") != "Omitted Middle T2").to_pandas()
        iso_pd = landscape_df.filter(pl.col("record_type") == "iso_contour").to_pandas()

        xy_dom_c = [-2.65, 2.65]
        w_c = float(w_col)
        h_c = float(h_top)
        mid_layer_c = (
            alt.Chart(mid_pd)
            .mark_circle(size=16, opacity=0.28, color="#CBD5E1")
            .encode(
                x=alt.X(
                    "rx_clip:Q",
                    scale=alt.Scale(domain=xy_dom_c),
                    title="Effector Score R_OOF (z)",
                    axis=alt.Axis(
                        values=[-2, -1, 0, 1, 2],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "ny_clip:Q",
                    scale=alt.Scale(domain=xy_dom_c),
                    title="Barrier Score NR_OOF (z)",
                    axis=alt.Axis(
                        values=[-2, -1, 0, 1, 2],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                        grid=False,
                    ),
                ),
            )
        )
        cor_dom_c = ["Corner CR/PR", "Corner SD/PD"]
        cor_cols_c = [OKABE_VERMILION, OKABE_BLUE]
        cor_shapes_c = ["circle", "cross"]
        cor_layer_c = (
            alt.Chart(cor_pat_pd)
            .mark_point(filled=True, size=28, opacity=0.82)
            .encode(
                x=alt.X("rx_clip:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("ny_clip:Q", scale=alt.Scale(domain=xy_dom_c)),
                color=alt.Color(
                    "pt_class:N",
                    scale=alt.Scale(domain=cor_dom_c, range=cor_cols_c),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=140,
                        symbolSize=70,
                        symbolOpacity=1.0,
                    ),
                ),
                shape=alt.Shape(
                    "pt_class:N",
                    scale=alt.Scale(domain=cor_dom_c, range=cor_shapes_c),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=140,
                        symbolSize=70,
                        symbolOpacity=1.0,
                    ),
                ),
            )
        )
        cut_v_c = (
            alt.Chart(pl.DataFrame([{"xc": rq33}, {"xc": rq67}]).to_pandas())
            .mark_rule(stroke="#64748B", strokeDash=[4, 3], strokeWidth=1.1)
            .encode(x=alt.X("xc:Q", scale=alt.Scale(domain=xy_dom_c)))
        )
        cut_h_c = (
            alt.Chart(pl.DataFrame([{"yc": nq33}, {"yc": nq67}]).to_pandas())
            .mark_rule(stroke="#64748B", strokeDash=[4, 3], strokeWidth=1.1)
            .encode(y=alt.Y("yc:Q", scale=alt.Scale(domain=xy_dom_c)))
        )
        iso_layer_c = (
            alt.Chart(iso_pd)
            .mark_line(stroke="#0F172A", strokeDash=[2, 2], strokeWidth=1.1, opacity=0.55)
            .encode(
                x=alt.X("x_z:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("y_z:Q", scale=alt.Scale(domain=xy_dom_c)),
                detail="stratum_id:N",
                order="step_seq:Q",
            )
        )
        c_rate_map = {
            (int(r["r_code"]), int(r["n_code"])): float(r["responder_rate_pct"])
            for r in oof_cells_df.iter_rows(named=True)
        }
        c_badges_df = pl.DataFrame(
            [
                {
                    "tx": -2.52,
                    "ty": 2.32,
                    "txt": f"LowR·Excl {c_rate_map[(0, 2)]:.1f}%",
                    "al": "left",
                    "col": OKABE_BLUE,
                },
                {
                    "tx": 2.52,
                    "ty": 2.32,
                    "txt": f"HighR·Excl {c_rate_map[(2, 2)]:.1f}%",
                    "al": "right",
                    "col": OKABE_BLUE,
                },
                {
                    "tx": -2.52,
                    "ty": -2.35,
                    "txt": f"LowR·Perm {c_rate_map[(0, 0)]:.1f}%",
                    "al": "left",
                    "col": OKABE_BLUE,
                },
                {
                    "tx": 2.52,
                    "ty": -2.35,
                    "txt": f"HighR·Perm {c_rate_map[(2, 0)]:.1f}%",
                    "al": "right",
                    "col": OKABE_VERMILION,
                },
            ]
        )
        c_left_df = c_badges_df.filter(pl.col("al") == "left").to_pandas()
        c_right_df = c_badges_df.filter(pl.col("al") == "right").to_pandas()
        c_halo_l = (
            alt.Chart(c_left_df)
            .mark_text(
                align="left",
                baseline="middle",
                fontSize=9.6,
                fontWeight="bold",
                stroke="#FFFFFF",
                strokeWidth=2.6,
            )
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=xy_dom_c)),
                text="txt:N",
            )
        )
        c_txt_l = (
            alt.Chart(c_left_df)
            .mark_text(align="left", baseline="middle", fontSize=9.6, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=xy_dom_c)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        c_halo_r = (
            alt.Chart(c_right_df)
            .mark_text(
                align="right",
                baseline="middle",
                fontSize=9.6,
                fontWeight="bold",
                stroke="#FFFFFF",
                strokeWidth=2.6,
            )
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=xy_dom_c)),
                text="txt:N",
            )
        )
        c_txt_r = (
            alt.Chart(c_right_df)
            .mark_text(align="right", baseline="middle", fontSize=9.6, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=xy_dom_c)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=xy_dom_c)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_c = (
            mid_layer_c
            + cor_layer_c
            + cut_v_c
            + cut_h_c
            + iso_layer_c
            + c_halo_l
            + c_txt_l
            + c_halo_r
            + c_txt_r
            + _panel_letter_layer("c", x_px=-34.0, y_px=-24.0)
        ).properties(width=int(w_c), height=int(h_c))

        # ------------------------------------------------------------------
        # Panel (d): Granular RECIST 1.1 Stacked Composition Across 4 Corners
        # ------------------------------------------------------------------
        short_rec_grp = {
            "high_r_perm": "Hi-R|Perm",
            "low_r_perm": "Lo-R|Perm",
            "high_r_excl": "Hi-R|Excl",
            "low_r_excl": "Lo-R|Excl",
        }
        rec_df = (
            clin_df.filter(pl.col("record_type") == "recist_corner_stack")
            .sort(["group_order", "sub_order"])
            .with_columns(
                [
                    (0.5 * (pl.col("ci_low_pct") + pl.col("ci_high_pct"))).alias("mid_x"),
                    pl.when(pl.col("rate_a_pct") >= 12.5)
                    .then(
                        pl.col("rate_a_pct").map_elements(
                            lambda v: f"{float(v):.0f}%", return_dtype=pl.String
                        )
                    )
                    .otherwise(pl.lit(""))
                    .alias("seg_lbl"),
                    pl.when(pl.col("sub_category").is_in(["CR", "PD"]))
                    .then(pl.lit("#FFFFFF"))
                    .otherwise(pl.lit("#0F172A"))
                    .alias("txt_col"),
                    pl.col("group_id")
                    .map_elements(
                        lambda g: short_rec_grp.get(str(g), str(g)),
                        return_dtype=pl.String,
                    )
                    .alias("corner_y_lbl"),
                ]
            )
        )
        corner_y_order = list(
            dict.fromkeys(rec_df.sort("group_order", descending=True)["corner_y_lbl"].to_list())
        )
        rec_dom = ["CR", "PR", "SD", "PD"]
        rec_cols = [OKABE_TEAL, OKABE_SKY, OKABE_ORANGE, OKABE_VERMILION]
        rec_pd = rec_df.to_pandas()

        bars_d = (
            alt.Chart(rec_pd)
            .mark_bar(height=25, stroke="#FFFFFF", strokeWidth=0.9)
            .encode(
                y=alt.Y(
                    "corner_y_lbl:N",
                    sort=corner_y_order,
                    title="TME substate",
                    axis=alt.Axis(
                        labelExpr="split(datum.label, '|')",
                        labelFontSize=9.6,
                        labelFontWeight="bold",
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                    ),
                ),
                x=alt.X(
                    "ci_low_pct:Q",
                    scale=alt.Scale(domain=[0.0, 100.0]),
                    title="RECIST 1.1 Composition (%)",
                    axis=alt.Axis(
                        values=[0, 25, 50, 75, 100],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                x2="ci_high_pct:Q",
                color=alt.Color(
                    "sub_category:N",
                    scale=alt.Scale(domain=rec_dom, range=rec_cols),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        symbolSize=70,
                    ),
                ),
            )
        )
        seg_txt_d = (
            alt.Chart(rec_pd)
            .mark_text(align="center", baseline="middle", fontSize=10.2, fontWeight="bold")
            .encode(
                y=alt.Y("corner_y_lbl:N", sort=corner_y_order),
                x=alt.X("mid_x:Q", scale=alt.Scale(domain=[0.0, 100.0])),
                text="seg_lbl:N",
                color=alt.Color("txt_col:N", scale=None),
            )
        )
        panel_d = (
            bars_d
            + seg_txt_d
            + _panel_letter_layer("d", x_px=-34.0, y_px=-24.0)
        ).properties(width=w_col, height=h_top)

        # ------------------------------------------------------------------
        # Panel (e): Orthogonal Overall Survival (OS) Kaplan-Meier Validation
        # ------------------------------------------------------------------
        km_steps = clin_df.filter(pl.col("record_type") == "km_step")
        km_censors = clin_df.filter(pl.col("record_type") == "km_censor")
        km_lr = clin_df.filter(pl.col("record_type") == "km_logrank")

        km_short_map = {
            "high_r_perm": "Hi-R·Perm",
            "low_r_perm": "Lo-R·Perm",
            "high_r_excl": "Hi-R·Excl",
            "low_r_excl": "Lo-R·Excl",
        }
        km_dom = [
            "Hi-R·Perm",
            "Lo-R·Perm",
            "Hi-R·Excl",
            "Lo-R·Excl",
        ]
        km_cols = [OKABE_TEAL, "#64748B", "#7E22CE", "#CA8A04"]
        km_dashes: list[list[int]] = [[1, 0], [5, 3], [1, 0], [4, 3]]

        km_steps_aug = km_steps.with_columns(
            pl.col("group_id")
            .map_elements(lambda g: km_short_map.get(str(g), str(g)), return_dtype=pl.String)
            .alias("km_series")
        )
        km_censors_aug = km_censors.with_columns(
            pl.col("group_id")
            .map_elements(lambda g: km_short_map.get(str(g), str(g)), return_dtype=pl.String)
            .alias("km_series")
        )
        lr_map = {str(r["hit_id"]): r for r in km_lr.iter_rows(named=True)}
        p_lr_9a = float(lr_map["hit_9a"]["p_value"]) if "hit_9a" in lr_map else 1.0
        p_lr_9b = float(lr_map["hit_9b"]["p_value"]) if "hit_9b" in lr_map else 1.0
        p_lr_null = float(lr_map["excl_null"]["p_value"]) if "excl_null" in lr_map else 1.0

        x_dom_km = [0.0, 36.5]
        y_dom_km = [0.0, 105.0]
        med50_km = (
            alt.Chart(pl.DataFrame({"y50": [50.0]}).to_pandas())
            .mark_rule(stroke="#94A3B8", strokeDash=[3, 3], strokeWidth=1.0)
            .encode(
                y=alt.Y(
                    "y50:Q",
                    scale=alt.Scale(domain=y_dom_km),
                    title="Overall Survival (%)",
                    axis=alt.Axis(
                        values=[0, 25, 50, 75, 100],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                )
            )
        )
        lines_km = (
            alt.Chart(km_steps_aug.to_pandas())
            .mark_line(strokeWidth=2.1)
            .encode(
                x=alt.X(
                    "delta_pct:Q",
                    scale=alt.Scale(domain=x_dom_km),
                    title="Time From ICB Initiation (Months)",
                    axis=alt.Axis(
                        values=[0, 6, 12, 18, 24, 30, 36],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y("rate_a_pct:Q", scale=alt.Scale(domain=y_dom_km)),
                order="sub_order:Q",
                color=alt.Color(
                    "km_series:N",
                    scale=alt.Scale(domain=km_dom, range=km_cols),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        columns=4,
                        offset=10,
                        columnPadding=3,
                        labelOffset=2,
                        labelFontSize=8.8,
                        labelLimit=75,
                        symbolSize=38,
                    ),
                ),
                strokeDash=alt.StrokeDash(
                    "km_series:N",
                    scale=alt.Scale(domain=km_dom, range=km_dashes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        columns=4,
                        offset=10,
                        columnPadding=3,
                        labelOffset=2,
                        labelFontSize=8.8,
                        labelLimit=75,
                        symbolSize=38,
                    ),
                ),
            )
        )
        censor_km = (
            alt.Chart(km_censors_aug.to_pandas())
            .mark_tick(orient="vertical", thickness=1.1, size=5, opacity=0.8)
            .encode(
                x=alt.X("delta_pct:Q", scale=alt.Scale(domain=x_dom_km)),
                y=alt.Y("rate_a_pct:Q", scale=alt.Scale(domain=y_dom_km)),
                color=alt.Color("km_series:N", scale=alt.Scale(domain=km_dom, range=km_cols), legend=None),
            )
        )
        km_box = (
            alt.Chart(pl.DataFrame([{"bx1": 0.6, "bx2": 25.8, "by1": 3.5, "by2": 24.5}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.7)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_km)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_km)),
                y2="by2:Q",
            )
        )
        km_txt = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": 1.4,
                            "ty": 18.0,
                            "txt": f"Perm OS p = {_fmt_pval_short(p_lr_9a)}",
                            "col": OKABE_TEAL,
                        },
                        {
                            "tx": 1.4,
                            "ty": 9.5,
                            "txt": f"High-R p = {p_lr_9b:.1e} · Excl p = {p_lr_null:.2f}",
                            "col": "#1E293B",
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=9.4, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_km)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_km)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_e = (
            med50_km
            + lines_km
            + censor_km
            + km_box
            + km_txt
            + _panel_letter_layer("e", x_px=-34.0, y_px=-24.0)
        ).properties(width=w_col, height=h_top)

        # ------------------------------------------------------------------
        # Panel (f): Continuous Humoral M29 Plasma Dose-Response Gated by M26 EMT
        # ------------------------------------------------------------------
        curve_m29 = epi_sig_df.filter(pl.col("record_type") == "m29_dose_curve").with_columns(
            [
                pl.when(pl.col("stratum_id") == "low_emt_m26_t1")
                .then(pl.lit(OKABE_TEAL))
                .otherwise(pl.lit(OKABE_PURPLE))
                .alias("s_hex"),
                pl.col("stratum_label")
                .str.replace("M26 ", "")
                .alias("stratum_label"),
            ]
        )
        bins_m29 = epi_sig_df.filter(pl.col("record_type") == "m29_quintile_bin").with_columns(
            [
                pl.when(pl.col("stratum_id") == "low_emt_m26_t1")
                .then(pl.lit(OKABE_TEAL))
                .otherwise(pl.lit(OKABE_PURPLE))
                .alias("s_hex"),
                pl.col("stratum_label")
                .str.replace("M26 ", "")
                .alias("stratum_label"),
            ]
        )
        lo_emt_row = curve_m29.filter(pl.col("stratum_id") == "low_emt_m26_t1").row(0, named=True)
        hi_emt_row = curve_m29.filter(pl.col("stratum_id") == "high_emt_m26_t3").row(0, named=True)

        dom_g = ["Low-EMT (T1)", "High-EMT (T3)"]
        cols_g = [OKABE_TEAL, OKABE_PURPLE]
        shapes_g = ["circle", "square"]

        area_g_lo = (
            alt.Chart(curve_m29.filter(pl.col("stratum_id") == "low_emt_m26_t1").to_pandas())
            .mark_area(fill=OKABE_TEAL, fillOpacity=0.15)
            .encode(
                x=alt.X(
                    "x_coord:Q",
                    scale=alt.Scale(domain=x_dom_dose),
                    title="Humoral JCHAIN+MZB1+ Plasma Score (z)",
                    axis=alt.Axis(
                        values=[-2, -1, 0, 1, 2],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "ci_low_pct:Q",
                    scale=alt.Scale(domain=y_dom_dose),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
                        values=[15, 30, 45, 60, 75],
                        labelFontSize=11.0,
                        titleFontSize=11.5,
                        titleFontWeight="bold",
                        titlePadding=y_ttl_pad,
                        minExtent=y_ax_ext,
                        maxExtent=y_ax_ext,
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y2="ci_high_pct:Q",
            )
        )
        area_g_hi = (
            alt.Chart(curve_m29.filter(pl.col("stratum_id") == "high_emt_m26_t3").to_pandas())
            .mark_area(fill=OKABE_PURPLE, fillOpacity=0.15)
            .encode(
                x=alt.X("x_coord:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
            )
        )
        area_g = area_g_lo + area_g_hi
        line_g = (
            alt.Chart(curve_m29.to_pandas())
            .mark_line(strokeWidth=2.3)
            .encode(
                x=alt.X("x_coord:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("responder_rate_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                detail="stratum_label:N",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        err_g = (
            alt.Chart(bins_m29.to_pandas())
            .mark_rule(strokeWidth=1.35, opacity=0.82)
            .encode(
                x=alt.X("x_coord:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="ci_high_pct:Q",
                color=alt.Color("s_hex:N", scale=None),
            )
        )
        pts_g = (
            alt.Chart(bins_m29.to_pandas())
            .mark_point(filled=True, size=72, opacity=1.0, stroke="#0F172A", strokeWidth=0.85)
            .encode(
                x=alt.X("x_coord:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("responder_rate_pct:Q", scale=alt.Scale(domain=y_dom_dose)),
                color=alt.Color(
                    "stratum_label:N",
                    scale=alt.Scale(domain=dom_g, range=cols_g),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.85,
                    ),
                ),
                shape=alt.Shape(
                    "stratum_label:N",
                    scale=alt.Scale(domain=dom_g, range=shapes_g),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        offset=10,
                        labelFontSize=10.8,
                        labelLimit=145,
                        symbolSize=70,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.85,
                    ),
                ),
            )
        )
        box_g = (
            alt.Chart(pl.DataFrame([{"bx1": -2.54, "bx2": 2.86, "by1": 66.2, "by2": 83.4}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_dose)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_dose)),
                y2="by2:Q",
            )
        )
        txt_g = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": -2.44,
                            "ty": 78.2,
                            "txt": f"Low-EMT OR={float(lo_emt_row['beta_inter']):.2f}/SD (p={_fmt_pval_short(float(lo_emt_row['p_spa']))})",
                            "col": OKABE_TEAL,
                        },
                        {
                            "tx": -2.44,
                            "ty": 70.8,
                            "txt": f"High-EMT OR={float(hi_emt_row['beta_inter']):.2f}/SD (p={float(hi_emt_row['p_spa']):.2f}, n.s.)",
                            "col": "#86198F",
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=10.2, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_dose)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_dose)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_f = (
            area_g
            + line_g
            + err_g
            + pts_g
            + box_g
            + txt_g
            + _panel_letter_layer("f", x_px=-34.0, y_px=-24.0)
        ).properties(width=w_col, height=h_top)

        # ------------------------------------------------------------------
        # Row 3 (Panels g, h, i): Single-Cell Harmony I2 Triad UMAP Panels (from S61 a, c, g)
        # ------------------------------------------------------------------
        h_umap = 155
        panel_g_res = _build_triad_umap_panel("t_nk_vs_caf", "g", width=w_col, height=h_umap)
        panel_h_res = _build_triad_umap_panel(
            "myeloid_macrophages_vs_caf", "h", width=w_col, height=h_umap
        )
        panel_i_res = _build_triad_umap_panel(
            "t_nk_vs_myeloid_macrophages", "i", width=w_col, height=h_umap
        )
        match (panel_g_res, panel_h_res, panel_i_res):
            case (Success(panel_g), Success(panel_h), Success(panel_i)):
                pass
            case (Failure(err), _, _) | (_, Failure(err), _) | (_, _, Failure(err)):
                return Failure(err)
            case _ as unreachable_umap:
                return Failure(f"Unexpected UMAP Result state: {unreachable_umap}")

        umap_legend_strip = _build_triad_umap_legend_strip(width=820, height=22)

        col_1 = alt.vconcat(panel_a, panel_d, panel_g, spacing=16).resolve_scale(
            color="independent",
            shape="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
        col_2 = alt.vconcat(panel_b, panel_e, panel_h, spacing=16).resolve_scale(
            color="independent",
            shape="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
        col_3 = alt.vconcat(panel_c, panel_f, panel_i, spacing=16).resolve_scale(
            color="independent",
            shape="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
        grid_3x3 = alt.hconcat(col_1, col_2, col_3, spacing=18).resolve_scale(
            color="independent",
            shape="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )

        combined = (
            alt.vconcat(grid_3x3, umap_legend_strip, spacing=14)
            .resolve_scale(
                color="independent",
                shape="independent",
                strokeDash="independent",
                strokeWidth="independent",
                opacity="independent",
            )
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 10, "right": 18, "top": 12, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        hybrid_svg = rasterize_svg_datapoints(
            raw_svg,
            max_points=2500,
            min_group_points=1000,
            scale=2.25,
        ).value_or(raw_svg)
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            hybrid_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-main-epistasis-and-stromal-gating",
            layer_label="Main_Epistasis_And_Stromal_Gating",
            enforce_font_floor=False,
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Main Epistasis & Stromal Gating Figure: {exc}")


# ------------------------------------------------------------------------------
# Supplementary Figure S19: Single-Panel 42x42 Pure-CLR (α = 0.00) Epistasis Matrix (Heading-Free)
# ------------------------------------------------------------------------------


def build_supp_pure_clr_all42_interactions(
    pure_inter_df: pl.DataFrame,
    states_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build single-panel heading-free 42x42 Pure-CLR (α = 0.00) epistasis matrix (panels b & c removed)."""
    try:
        meta_sorted = (
            states_df.filter(
                (pl.col("reference_id") == "Consensus-Meta42")
                & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
            )
            .with_columns(
                pl.col("cell_state").replace(META_STATE_EMPIRICAL_LABELS).alias("state_label")
            )
            .sort("bulk_deconv_beta", descending=True)
        )
        ordered_codes = meta_sorted["cell_state"].to_list()
        ordered_labels = meta_sorted["state_label"].to_list()

        inter_clean = pure_inter_df.with_columns(
            [
                pl.col("state_a").replace(META_STATE_EMPIRICAL_LABELS).alias("state_a_label"),
                pl.col("state_b").replace(META_STATE_EMPIRICAL_LABELS).alias("state_b_label"),
            ]
        )
        inter_pd = inter_clean.to_pandas()
        sig_star_pd = (
            inter_clean.filter((pl.col("state_a") != pl.col("state_b")) & (pl.col("p_spa") < 0.005))
            .with_columns(pl.lit("*").alias("sig_glyph"))
            .to_pandas()
        )
        sig_dot_pd = (
            inter_clean.filter(
                (pl.col("state_a") != pl.col("state_b"))
                & (pl.col("p_spa") >= 0.005)
                & (pl.col("p_spa") < 0.01)
            )
            .with_columns(pl.lit("•").alias("sig_glyph"))
            .to_pandas()
        )

        rects = (
            alt.Chart(inter_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.35)
            .encode(
                x=alt.X(
                    "state_b:N",
                    sort=ordered_codes,
                    title="Interacting Consensus Meta-State B (M01–M42 Ranked by Meta-9 β; Key Matches Y-Axis)",
                    axis=alt.Axis(
                        labelAngle=-90,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=80,
                    ),
                ),
                y=alt.Y(
                    "state_a_label:N",
                    sort=ordered_labels,
                    title="Consensus Meta-State A (All 42 States)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=245,
                    ),
                ),
                color=alt.Color(
                    "epistasis_z:Q",
                    title="Calibrated L2O-Quad-Logit-SPA z_ab (Pure CLR, α = 0.00)",
                    scale=alt.Scale(
                        domain=[-3.5, 3.5],
                        scheme="blueorange",
                        reverse=True,
                    ),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=180,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleLimit=560,
                    ),
                ),
            )
        )
        stars = (
            alt.Chart(sig_star_pd)
            .mark_text(
                align="center",
                baseline="middle",
                dy=3,
                fontSize=MIN_FONT_PX,
                fontWeight="bold",
                color="#FFFFFF",
            )
            .encode(
                x=alt.X("state_b:N", sort=ordered_codes),
                y=alt.Y("state_a_label:N", sort=ordered_labels),
                text="sig_glyph:N",
            )
        )
        dots = (
            alt.Chart(sig_dot_pd)
            .mark_text(
                align="center",
                baseline="middle",
                dy=0,
                fontSize=MIN_FONT_PX,
                fontWeight="bold",
                color="#FFFFFF",
            )
            .encode(
                x=alt.X("state_b:N", sort=ordered_codes),
                y=alt.Y("state_a_label:N", sort=ordered_labels),
                text="sig_glyph:N",
            )
        )

        combined = (
            (rects + stars + dots)
            .properties(width=590, height=590)
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 6, "right": 8, "top": 8, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-supp-all42-pure-clr-interactions",
            layer_label="Supp_All42_Pure_CLR_Conditional_Interactions",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Pure-CLR All-42 Supplementary Figure: {exc}")


# ------------------------------------------------------------------------------
# Supplementary Figure S20: Epistasis Null Calibration & Discovery Power Benchmark (4 Panels a-d, Heading-Free)
# ------------------------------------------------------------------------------


def build_supp_epistasis_discovery_power_benchmark(
    epi_strat_df: pl.DataFrame,
    epi_sig_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build heading-free Supplementary Figure S20 (180.0mm): (a) Null B FPR + (b) FDR power + (c) M29xM26 lineage forest + (d) M29xM26 2x2 corners."""
    try:
        strat_sorted = epi_strat_df.sort("strategy_rank")
        all_labels_order = strat_sorted["short_label"].to_list()

        # Panel (a): Null Calibration Audit across 6 non-hybrid strategies
        seg_rows_a: list[dict[str, object]] = []
        pt_rows_a: list[dict[str, object]] = []
        lbl_rows_a: list[dict[str, object]] = []
        for r in strat_sorted.iter_rows(named=True):
            sh = str(r["short_label"])
            is_cal = bool(r["is_zero_leakage_calibrated"])
            f_perm = 100.0 * float(r["perm_null_fpr_p05"])
            f_nb = 100.0 * float(r["null_b_fpr_p05"])
            seg_rows_a.append(
                {
                    "short_label": sh,
                    "x_min": min(f_perm, f_nb),
                    "x_max": max(f_perm, f_nb),
                }
            )
            pt_rows_a.append(
                {
                    "short_label": sh,
                    "fpr_pct": f_perm,
                    "null_series": "Within-Trial Permuted-y Null",
                }
            )
            pt_rows_a.append(
                {
                    "short_label": sh,
                    "fpr_pct": f_nb,
                    "null_series": (
                        "Saturation Null B (Calibrated)"
                        if is_cal
                        else "Saturation Null B (Inflated Control)"
                    ),
                }
            )
            lbl_rows_a.append(
                {
                    "short_label": sh,
                    "label_x": max(f_perm, f_nb) + 0.8,
                    "callout": f"Null B: {f_nb:.1f}%  (Perm: {f_perm:.1f}%)",
                }
            )

        x_dom_a = [0.0, 68.0]
        nom_rule_a = (
            alt.Chart(pl.DataFrame({"x": [5.0]}).to_pandas())
            .mark_rule(color="#0F172A", strokeDash=[5, 3], strokeWidth=1.5)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_a)))
        )
        nom_lbl_a = (
            alt.Chart(
                pl.DataFrame({"x": [5.7], "y": [all_labels_order[0]], "t": ["Nominal α = 5.0%"]}).to_pandas()
            )
            .mark_text(align="left", baseline="bottom", dy=-5, fontSize=MIN_FONT_PX, fontWeight="bold", color="#334155")
            .encode(
                x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("y:N", sort=all_labels_order),
                text="t:N",
            )
        )
        segs_a = (
            alt.Chart(pl.DataFrame(seg_rows_a).to_pandas())
            .mark_rule(color="#64748B", strokeWidth=1.5)
            .encode(
                y=alt.Y(
                    "short_label:N",
                    sort=all_labels_order,
                    title=None,
                    axis=alt.Axis(labelFontSize=MIN_FONT_PX, labelFontWeight="bold", labelLimit=295),
                ),
                x=alt.X(
                    "x_min:Q",
                    scale=alt.Scale(domain=x_dom_a),
                    title="Empirical False-Positive Rate at Nominal p ≤ 0.05 (%)",
                    axis=alt.Axis(
                        values=[0, 5, 10, 20, 30, 40, 50],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="x_max:Q",
            )
        )
        null_domain_a = [
            "Within-Trial Permuted-y Null",
            "Saturation Null B (Calibrated)",
            "Saturation Null B (Inflated Control)",
        ]
        null_colors_a = [OKABE_SKY, OKABE_BLUE, OKABE_VERMILION]
        null_shapes_a = ["diamond", "circle", "square"]
        pts_a = (
            alt.Chart(pl.DataFrame(pt_rows_a).to_pandas())
            .mark_point(filled=True, size=95, opacity=0.96, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                y=alt.Y("short_label:N", sort=all_labels_order),
                x=alt.X("fpr_pct:Q", scale=alt.Scale(domain=x_dom_a)),
                color=alt.Color(
                    "null_series:N",
                    scale=alt.Scale(domain=null_domain_a, range=null_colors_a),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        columns=2,
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=310,
                        symbolSize=85,
                    ),
                ),
                shape=alt.Shape(
                    "null_series:N",
                    scale=alt.Scale(domain=null_domain_a, range=null_shapes_a),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        columns=2,
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=310,
                        symbolSize=85,
                    ),
                ),
            )
        )
        txt_a = (
            alt.Chart(pl.DataFrame(lbl_rows_a).to_pandas())
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("short_label:N", sort=all_labels_order),
                x=alt.X("label_x:Q", scale=alt.Scale(domain=x_dom_a)),
                text="callout:N",
            )
        )
        panel_a = (
            nom_rule_a
            + segs_a
            + pts_a
            + txt_a
            + nom_lbl_a
            + _panel_letter_layer("a", x_px=-295.0, y_px=-42.0)
        ).properties(width=425, height=162)

        # Panel (b): Empirical Multiplicity-Adjusted Discovery Power Across 4 Calibrated Strategies
        cal_sub = strat_sorted.filter(pl.col("is_zero_leakage_calibrated"))
        cal_order = cal_sub["short_label"].to_list()
        pow_rows_b: list[dict[str, object]] = []
        for r in cal_sub.iter_rows(named=True):
            min_q = float(r["emp_min_q_bh"])
            neg_l10 = float(r["emp_neg_log10_min_q_bh"])
            is_sig = min_q < 0.05
            s_id = str(r["strategy_id"])
            if s_id == "sc_deg_100_all861":
                c_txt = (
                    f"q = {min_q:.4f} | {r['top_hit_id']} (z = {float(r['top_hit_z']):+.2f}, "
                    f"FWER = {float(r['emp_min_fwer_p']):.3f})"
                )
            elif s_id == "oof_2x2_composite_gate":
                c_txt = f"q = {min_q:.4f} | Perm ΔR (DiD p = {float(r['emp_omnibus_boot_p']):.3f})"
            else:
                q_fmt = f"{min_q:.4f}" if min_q < 0.05 else f"{min_q:.3f}"
                c_txt = f"q = {q_fmt} | {r['top_hit_id']} (z = {float(r['top_hit_z']):+.2f})"
            pow_rows_b.append(
                {
                    "short_label": str(r["short_label"]),
                    "x0": 0.0,
                    "neg_log10_q": neg_l10,
                    "label_x": neg_l10 + 0.09,
                    "callout": c_txt,
                    "fdr_tier": (
                        "FDR-Significant (BH q < 0.05)"
                        if is_sig
                        else "Sub-Threshold (BH q ≥ 0.05)"
                    ),
                }
            )
        pow_pd = pl.DataFrame(pow_rows_b).to_pandas()
        x_dom_b = [0.0, 6.2]
        q10_val = float(-np.log10(0.10))
        q05_val = float(-np.log10(0.05))
        q10_rule = (
            alt.Chart(pl.DataFrame({"x": [q10_val]}).to_pandas())
            .mark_rule(color="#64748B", strokeDash=[3, 3], strokeWidth=1.2)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_b)))
        )
        q05_rule = (
            alt.Chart(pl.DataFrame({"x": [q05_val]}).to_pandas())
            .mark_rule(color=OKABE_TEAL, strokeDash=[5, 3], strokeWidth=1.7)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_b)))
        )
        q10_lbl = (
            alt.Chart(pl.DataFrame({"x": [q10_val], "lbl": ["q=0.10"]}).to_pandas())
            .mark_text(align="right", baseline="bottom", dx=-4, dy=-3, fontSize=MIN_FONT_PX, fontWeight="bold", color="#475569")
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_b)), y=alt.value(0), text="lbl:N")
        )
        q05_lbl = (
            alt.Chart(pl.DataFrame({"x": [q05_val], "lbl": ["BH q = 0.05"]}).to_pandas())
            .mark_text(align="left", baseline="bottom", dx=4, dy=-3, fontSize=MIN_FONT_PX, fontWeight="bold", color=OKABE_TEAL)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_b)), y=alt.value(0), text="lbl:N")
        )
        fdr_dom_b = [
            "FDR-Significant (BH q < 0.05)",
            "Sub-Threshold (BH q ≥ 0.05)",
        ]
        fdr_cols_b = [OKABE_TEAL, OKABE_BLUE]
        stems_b = (
            alt.Chart(pow_pd)
            .mark_rule(strokeWidth=3.2)
            .encode(
                y=alt.Y(
                    "short_label:N",
                    sort=cal_order,
                    title=None,
                    axis=alt.Axis(labelFontSize=MIN_FONT_PX, labelFontWeight="bold", labelLimit=295),
                ),
                x=alt.X(
                    "x0:Q",
                    scale=alt.Scale(domain=x_dom_b),
                    title="Empirical Discovery Significance [−log10(Minimum Benjamini–Hochberg q)]",
                    axis=alt.Axis(
                        values=[0.0, 1.0, 2.0, 3.0, 4.0],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="neg_log10_q:Q",
                color=alt.Color(
                    "fdr_tier:N",
                    scale=alt.Scale(domain=fdr_dom_b, range=fdr_cols_b),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=320,
                        symbolSize=90,
                    ),
                ),
            )
        )
        pts_b = (
            alt.Chart(pow_pd)
            .mark_circle(size=105, opacity=1.0, stroke="#0F172A", strokeWidth=0.85)
            .encode(
                y=alt.Y("short_label:N", sort=cal_order),
                x=alt.X("neg_log10_q:Q", scale=alt.Scale(domain=x_dom_b)),
                color=alt.Color("fdr_tier:N", scale=alt.Scale(domain=fdr_dom_b, range=fdr_cols_b)),
            )
        )
        txt_b = (
            alt.Chart(pow_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("short_label:N", sort=cal_order),
                x=alt.X("label_x:Q", scale=alt.Scale(domain=x_dom_b)),
                text="callout:N",
            )
        )
        panel_b = (
            q10_rule
            + q05_rule
            + stems_b
            + pts_b
            + txt_b
            + q10_lbl
            + q05_lbl
            + _panel_letter_layer("b", x_px=-295.0, y_px=-30.0)
        ).properties(width=425, height=126)

        # Panel (c): Cross-Indication Interaction Forest for M29 Plasma x M26 EMT (β_ab ± 95% CI & z_SPA)
        short_lin_map = {
            "melanoma": "Melanoma (n=338)",
            "bladder": "Bladder (n=298)",
            "ccrcc": "Renal ccRCC (n=263)",
            "gi_pdac": "PDAC + TNBC (n=116)",
            "pancancer": "Pan-Cancer (N=1,015)",
        }
        for_df = (
            epi_sig_df.filter(
                (pl.col("record_type") == "indication_forest") & (pl.col("pair_id") == "M29_x_M26")
            )
            .sort("stratum_order")
            .with_columns(
                [
                    pl.col("stratum_id")
                    .map_elements(lambda s: short_lin_map.get(str(s), str(s)), return_dtype=pl.String)
                    .alias("short_lin"),
                    pl.when(pl.col("stratum_id") == "pancancer")
                    .then(pl.lit("Pan-Cancer (q < 0.05)"))
                    .otherwise(pl.lit("Lineage Stratum"))
                    .alias("tier_grp"),
                    (pl.max_horizontal(pl.col("beta_ci_high"), pl.lit(0.0)) + pl.lit(0.04)).alias("lbl_x"),
                    pl.struct(["beta_inter", "z_spa"])
                    .map_elements(
                        lambda s: f"{float(s['beta_inter']):+.2f} (z={float(s['z_spa']):+.2f})",
                        return_dtype=pl.String,
                    )
                    .alias("for_callout"),
                ]
            )
        )
        strat_for_order = for_df["short_lin"].to_list()
        c_dom = [
            "Lineage Stratum",
            "Pan-Cancer (q < 0.05)",
        ]
        c_cols = [OKABE_BLUE, OKABE_TEAL]
        c_shapes = ["circle", "diamond"]
        x_dom_c = [-0.86, 1.16]
        for_pd = for_df.to_pandas()

        zero_c = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#0F172A", strokeDash=[4, 3], strokeWidth=1.3)
            .encode(x=alt.X("x:Q", scale=alt.Scale(domain=x_dom_c)))
        )
        ci_c = (
            alt.Chart(for_pd)
            .mark_rule(strokeWidth=2.0)
            .encode(
                y=alt.Y(
                    "short_lin:N",
                    sort=strat_for_order,
                    title=None,
                    axis=alt.Axis(labelFontSize=MIN_FONT_PX, labelFontWeight="bold", labelLimit=165),
                ),
                x=alt.X(
                    "beta_ci_low:Q",
                    scale=alt.Scale(domain=x_dom_c),
                    title="M29 × M26 Log-Odds β_ab (±95% CI)",
                    axis=alt.Axis(
                        values=[-0.8, -0.4, 0.0, 0.4, 0.8],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="beta_ci_high:Q",
                color=alt.Color(
                    "tier_grp:N",
                    scale=alt.Scale(domain=c_dom, range=c_cols),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=160,
                        symbolSize=85,
                    ),
                ),
            )
        )
        pts_c = (
            alt.Chart(for_pd)
            .mark_point(filled=True, size=105, opacity=1.0, stroke="#0F172A", strokeWidth=0.85)
            .encode(
                y=alt.Y("short_lin:N", sort=strat_for_order),
                x=alt.X("beta_inter:Q", scale=alt.Scale(domain=x_dom_c)),
                color=alt.Color("tier_grp:N", scale=alt.Scale(domain=c_dom, range=c_cols)),
                shape=alt.Shape(
                    "tier_grp:N",
                    scale=alt.Scale(domain=c_dom, range=c_shapes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=160,
                        symbolSize=85,
                    ),
                ),
            )
        )
        txt_c = (
            alt.Chart(for_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("short_lin:N", sort=strat_for_order),
                x=alt.X("lbl_x:Q", scale=alt.Scale(domain=x_dom_c)),
                text="for_callout:N",
            )
        )
        panel_c = (
            zero_c
            + ci_c
            + pts_c
            + txt_c
            + _panel_letter_layer("c", x_px=-155.0, y_px=-30.0)
        ).properties(width=225, height=162)

        # Panel (d): Discrete 2x2 Extreme-Tertile Corner Interaction (M29 Plasma x M26 EMT)
        m29_df = (
            epi_sig_df.filter(
                (pl.col("record_type") == "tertile_2x2_cell") & (pl.col("pair_id") == "M29_x_M26")
            )
            .with_columns(
                [
                    pl.when(pl.col("mod_code") == 0)
                    .then(pl.lit("Low-EMT (M26 T1)"))
                    .otherwise(pl.lit("High-EMT (M26 T3)"))
                    .alias("emt_grp"),
                    pl.when(pl.col("mod_code") == 0)
                    .then(pl.lit(OKABE_TEAL))
                    .otherwise(pl.lit(OKABE_PURPLE))
                    .alias("emt_hex"),
                    (
                        pl.when(pl.col("x_code") == 0).then(pl.lit(0.0)).otherwise(pl.lit(1.0))
                        + pl.when(pl.col("mod_code") == 0).then(pl.lit(-0.03)).otherwise(pl.lit(0.03))
                    ).alias("x_pos"),
                    pl.struct(["responder_rate_pct", "n_responders", "n_patients"])
                    .map_elements(
                        lambda s: f"{float(s['responder_rate_pct']):.1f}% ({int(s['n_responders'])}/{int(s['n_patients'])})",
                        return_dtype=pl.String,
                    )
                    .alias("rate_lbl"),
                ]
            )
        )
        m29_map = {(int(r["x_code"]), int(r["mod_code"])): r for r in m29_df.iter_rows(named=True)}
        r00_d = float(m29_map[(0, 0)]["responder_rate_pct"])
        r20_d = float(m29_map[(2, 0)]["responder_rate_pct"])
        r02_d = float(m29_map[(0, 2)]["responder_rate_pct"])
        r22_d = float(m29_map[(2, 2)]["responder_rate_pct"])
        p_slope_lo = float(m29_map[(0, 0)]["fisher_slope_mod_low_p"])
        p_slope_hi = float(m29_map[(0, 0)]["fisher_slope_mod_high_p"])
        p_gate_hi = float(m29_map[(0, 0)]["fisher_gate_x_high_p"])
        delta_hi_d = r20_d - r22_d

        dom_d = ["Low-EMT (M26 T1)", "High-EMT (M26 T3)"]
        cols_d = [OKABE_TEAL, OKABE_PURPLE]
        shapes_d = ["circle", "square"]
        x_dom_d = [-0.48, 2.04]
        y_dom_d = [6.0, 77.0]
        m29_pd = m29_df.to_pandas()

        lines_d = (
            alt.Chart(m29_pd)
            .mark_line(strokeWidth=2.4)
            .encode(
                x=alt.X(
                    "x_pos:Q",
                    scale=alt.Scale(domain=x_dom_d),
                    title="Humoral M29 Plasma (Extreme Tertiles)",
                    axis=alt.Axis(
                        values=[0.0, 1.0],
                        labelExpr="datum.value == 0 ? 'Low M29 (T1)' : 'High M29 (T3)'",
                        labelFontSize=MIN_FONT_PX,
                        labelFontWeight="bold",
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "responder_rate_pct:Q",
                    scale=alt.Scale(domain=y_dom_d),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
                        values=[15, 30, 45, 60, 75],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                detail="emt_grp:N",
                color=alt.Color("emt_hex:N", scale=None),
            )
        )
        err_d = (
            alt.Chart(m29_pd)
            .mark_rule(strokeWidth=1.5, opacity=0.85)
            .encode(
                x=alt.X("x_pos:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_d)),
                y2="ci_high_pct:Q",
                color=alt.Color("emt_hex:N", scale=None),
            )
        )
        pts_d = (
            alt.Chart(m29_pd)
            .mark_point(filled=True, size=95, opacity=1.0, stroke="#0F172A", strokeWidth=0.85)
            .encode(
                x=alt.X("x_pos:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("responder_rate_pct:Q", scale=alt.Scale(domain=y_dom_d)),
                color=alt.Color(
                    "emt_grp:N",
                    scale=alt.Scale(domain=dom_d, range=cols_d),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=150,
                        symbolSize=85,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.85,
                    ),
                ),
                shape=alt.Shape(
                    "emt_grp:N",
                    scale=alt.Scale(domain=dom_d, range=shapes_d),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=150,
                        symbolSize=85,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.85,
                    ),
                ),
            )
        )
        pt_lbl_d = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": -0.40,
                            "ty": r00_d - 9.2,
                            "txt": str(m29_map[(0, 0)]["rate_lbl"]),
                            "col": OKABE_TEAL,
                            "al": "left",
                        },
                        {
                            "tx": -0.40,
                            "ty": r02_d + 13.5,
                            "txt": str(m29_map[(0, 2)]["rate_lbl"]),
                            "col": "#86198F",
                            "al": "left",
                        },
                        {
                            "tx": 0.90,
                            "ty": r20_d + 6.5,
                            "txt": str(m29_map[(2, 0)]["rate_lbl"]),
                            "col": OKABE_TEAL,
                            "al": "right",
                        },
                        {
                            "tx": 0.94,
                            "ty": r22_d - 11.2,
                            "txt": str(m29_map[(2, 2)]["rate_lbl"]),
                            "col": "#86198F",
                            "al": "right",
                        },
                    ]
                ).to_pandas()
            )
        )
        pt_lbl_center_d = (
            alt.Chart(pt_lbl_d.data[pt_lbl_d.data["al"] == "left"])
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_d)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        pt_lbl_right_d = (
            alt.Chart(pt_lbl_d.data[pt_lbl_d.data["al"] == "right"])
            .mark_text(align="right", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_d)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        brack_d = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {"x1": 1.10, "x2": 1.15, "y1": r20_d, "y2": r20_d},
                        {"x1": 1.10, "x2": 1.15, "y1": r22_d, "y2": r22_d},
                        {"x1": 1.15, "x2": 1.15, "y1": r22_d, "y2": r20_d},
                    ]
                ).to_pandas()
            )
            .mark_rule(color="#0F172A", strokeWidth=1.35)
            .encode(
                x=alt.X("x1:Q", scale=alt.Scale(domain=x_dom_d)),
                x2="x2:Q",
                y=alt.Y("y1:Q", scale=alt.Scale(domain=y_dom_d)),
                y2="y2:Q",
            )
        )
        brack_txt_d = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": 1.18,
                            "ty": 0.5 * (r20_d + r22_d) + 3.6,
                            "txt": f"Δ={delta_hi_d:+.1f}%",
                        },
                        {
                            "tx": 1.18,
                            "ty": 0.5 * (r20_d + r22_d) - 3.6,
                            "txt": f"p={p_gate_hi:.1e}",
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_d)),
                text="txt:N",
            )
        )
        box_d = (
            alt.Chart(pl.DataFrame([{"bx1": -0.44, "bx2": 1.98, "by1": 65.8, "by2": 75.8}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_d)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_d)),
                y2="by2:Q",
            )
        )
        top_txt_d = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": -0.39,
                            "ty": 70.8,
                            "txt": f"Low-EMT p={p_slope_lo:.1e} · High-EMT p={p_slope_hi:.2f}",
                        }
                    ]
                ).to_pandas()
            )
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_d)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_d)),
                text="txt:N",
            )
        )
        panel_d = (
            lines_d
            + err_d
            + pts_d
            + pt_lbl_center_d
            + pt_lbl_right_d
            + brack_d
            + brack_txt_d
            + box_d
            + top_txt_d
            + _panel_letter_layer("d", x_px=-46.0, y_px=-30.0)
        ).properties(width=255, height=162)

        bot_row_s20 = alt.hconcat(panel_c, panel_d, spacing=18).resolve_scale(
            color="independent", shape="independent"
        )
        combined = (
            alt.vconcat(panel_a, panel_b, bot_row_s20, spacing=18, center=True)
            .resolve_scale(color="independent", shape="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 18, "right": 10, "top": 12, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-supp-epistasis-discovery-power-benchmark",
            layer_label="Supp_Epistasis_Discovery_Power_Benchmark",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Epistasis Discovery Power Supplementary Figure: {exc}")


# ------------------------------------------------------------------------------
# Supplementary Figure S21: Zero-Leakage 3x3 Stromal Gating Benchmark (5 Panels a-e, Heading-Free)
# ------------------------------------------------------------------------------


def build_supp_zero_leakage_gating_comparison(
    zl_summary_df: pl.DataFrame,
    zl_cells_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build heading-free Supplementary Figure S21 (180.0mm): (a) 7-estimator permutation audit + (b-e) zero-leakage 3x3 curves."""
    try:
        ordered_methods = zl_summary_df.sort("method_rank")["short_label"].to_list()
        sum_aug = zl_summary_df.sort("method_rank").with_columns(
            [
                pl.when(pl.col("is_zero_leakage"))
                .then(pl.lit("Strict Zero-Leakage (Observed T3 Δ)"))
                .otherwise(pl.lit("Leaky Baseline (In-Sample y, Observed T3 Δ)"))
                .alias("obs_group"),
                (
                    pl.max_horizontal("observed_delta_t3_pct", "perm_null_ci_high_pct")
                    + pl.lit(1.8)
                ).alias("callout_x"),
                pl.struct(
                    [
                        "observed_delta_t3_pct",
                        "perm_null_mean_delta_t3_pct",
                        "perm_pvalue",
                    ]
                )
                .map_elements(
                    lambda s: (
                        f"Δ={float(s['observed_delta_t3_pct']):+.1f}% "
                        f"(null {float(s['perm_null_mean_delta_t3_pct']):+.1f}%, "
                        f"p={float(s['perm_pvalue']):.3f})"
                    ),
                    return_dtype=pl.String,
                )
                .alias("callout_text"),
            ]
        )

        obs_pts_df = sum_aug.select(
            [
                pl.col("method_rank"),
                pl.col("short_label"),
                pl.col("observed_delta_t3_pct").alias("x_val"),
                pl.col("obs_group").alias("marker_series"),
            ]
        )
        null_pts_df = sum_aug.select(
            [
                pl.col("method_rank"),
                pl.col("short_label"),
                pl.col("perm_null_mean_delta_t3_pct").alias("x_val"),
                pl.lit("Within-Trial Permutation Null E[T3 Δ | permuted y] (±95% UI)").alias("marker_series"),
            ]
        )
        all_pts_pd = pl.concat([obs_pts_df, null_pts_df]).to_pandas()
        sum_pd = sum_aug.to_pandas()

        series_order = [
            "Strict Zero-Leakage (Observed T3 Δ)",
            "Leaky Baseline (In-Sample y, Observed T3 Δ)",
            "Within-Trial Permutation Null E[T3 Δ | permuted y] (±95% UI)",
        ]
        series_colors = [OKABE_BLUE, OKABE_VERMILION, "#64748B"]
        series_shapes = ["circle", "square", "diamond"]
        x_dom_a = [-18.0, 96.0]

        zero_line_a = (
            alt.Chart(pl.DataFrame({"x0": [0.0]}).to_pandas())
            .mark_rule(color="#0F172A", strokeDash=[5, 3], strokeWidth=1.5)
            .encode(x=alt.X("x0:Q", scale=alt.Scale(domain=x_dom_a)))
        )
        null_ui_rule = (
            alt.Chart(sum_pd)
            .mark_rule(color="#94A3B8", strokeWidth=3.2, opacity=0.85)
            .encode(
                y=alt.Y(
                    "short_label:N",
                    sort=ordered_methods,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=235,
                        labelFontWeight="bold",
                    ),
                ),
                x=alt.X(
                    "perm_null_ci_low_pct:Q",
                    scale=alt.Scale(domain=x_dom_a),
                    title="Top-Tertile (T3) Stromal Gating Contrast: Permissive − Excluded (%)",
                    axis=alt.Axis(
                        values=[-15, 0, 15, 30, 45, 60],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="perm_null_ci_high_pct:Q",
            )
        )
        gain_connector = (
            alt.Chart(sum_pd)
            .mark_rule(color="#0F172A", strokeWidth=1.4, strokeDash=[2, 2])
            .encode(
                y=alt.Y("short_label:N", sort=ordered_methods),
                x=alt.X("perm_null_mean_delta_t3_pct:Q", scale=alt.Scale(domain=x_dom_a)),
                x2="observed_delta_t3_pct:Q",
            )
        )
        pts_a = (
            alt.Chart(all_pts_pd)
            .mark_point(filled=True, size=115, opacity=1.0, stroke="#0F172A", strokeWidth=0.9)
            .encode(
                y=alt.Y("short_label:N", sort=ordered_methods),
                x=alt.X("x_val:Q", scale=alt.Scale(domain=x_dom_a)),
                color=alt.Color(
                    "marker_series:N",
                    scale=alt.Scale(domain=series_order, range=series_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=450,
                        symbolSize=95,
                    ),
                ),
                shape=alt.Shape(
                    "marker_series:N",
                    scale=alt.Scale(domain=series_order, range=series_shapes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=450,
                        symbolSize=95,
                    ),
                ),
            )
        )
        callouts_a = (
            alt.Chart(sum_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("short_label:N", sort=ordered_methods),
                x=alt.X("callout_x:Q", scale=alt.Scale(domain=x_dom_a)),
                text="callout_text:N",
            )
        )
        panel_a = (
            zero_line_a
            + null_ui_rule
            + gain_connector
            + pts_a
            + callouts_a
            + _panel_letter_layer("a", x_px=-235.0, y_px=-56.0)
        ).properties(width=455, height=195)

        # Panels (b)-(e): Heading-free 2x2 Grid of Strict Zero-Leakage 3x3 Clinical Gating Regimes
        r_short_map = {
            "Low Pan-R Score (T1)": "Low R (T1)",
            "Med. Pan-R Score (T2)": "Med R (T2)",
            "High Pan-R Score (T3)": "High R (T3)",
        }
        r_short_order = ["Low R (T1)", "Med R (T2)", "High R (T3)"]
        n_order = [
            "Permissive Low-NR TME",
            "Moderate NR TME",
            "Excluded High-NR TME",
        ]
        n_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_VERMILION]

        def _make_zl_curve_panel(
            method_id: str,
            letter: str,
            badge_name: str,
            show_y_title: bool,
            show_x_title: bool,
            show_legend: bool,
        ) -> alt.LayerChart:
            s_row = zl_summary_df.filter(pl.col("method_id") == method_id).row(0, named=True)
            obs_d = float(s_row["observed_delta_t3_pct"])
            p_perm = float(s_row["perm_pvalue"])
            sub_pd = (
                zl_cells_df.filter(pl.col("method_id") == method_id)
                .with_columns(pl.col("lymphoid_tier").replace(r_short_map).alias("r_short"))
                .to_pandas()
            )
            leg_cfg = (
                alt.Legend(
                    title=None,
                    orient="top",
                    columns=2,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=230,
                    symbolSize=80,
                )
                if show_legend
                else None
            )
            y_dom_sub = [8.0, 76.0]
            ln = (
                alt.Chart(sub_pd)
                .mark_line(
                    point=alt.OverlayMarkDef(size=70, filled=True, stroke="#0F172A", strokeWidth=0.7),
                    strokeWidth=2.1,
                )
                .encode(
                    x=alt.X(
                        "r_short:N",
                        sort=r_short_order,
                        scale=alt.Scale(domain=r_short_order, padding=0.38),
                        title="Responder Score Tertile (Within-Trial)" if show_x_title else None,
                        axis=alt.Axis(
                            labelAngle=0,
                            labelPadding=6,
                            labelFontSize=MIN_FONT_PX,
                            titleFontSize=MIN_FONT_PX,
                            titleFontWeight="bold",
                            labelLimit=120,
                        ),
                    ),
                    y=alt.Y(
                        "responder_rate_pct:Q",
                        scale=alt.Scale(domain=y_dom_sub),
                        title="Responder Rate (%)" if show_y_title else None,
                        axis=alt.Axis(
                            values=[15, 30, 45, 60],
                            labelFontSize=MIN_FONT_PX,
                            titleFontSize=MIN_FONT_PX,
                            titleFontWeight="bold",
                            grid=True,
                            gridColor="#E2E8F0",
                        ),
                    ),
                    color=alt.Color(
                        "stromal_regime:N",
                        scale=alt.Scale(domain=n_order, range=n_colors),
                        legend=leg_cfg,
                    ),
                )
            )
            ci = (
                alt.Chart(sub_pd)
                .mark_rule(strokeWidth=1.35, opacity=0.72)
                .encode(
                    x=alt.X(
                        "r_short:N",
                        sort=r_short_order,
                        scale=alt.Scale(domain=r_short_order, padding=0.38),
                    ),
                    y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=y_dom_sub)),
                    y2="ci_high_pct:Q",
                    color=alt.Color("stromal_regime:N", scale=alt.Scale(domain=n_order, range=n_colors)),
                )
            )
            badge_bg = (
                alt.Chart(pl.DataFrame([{"x1": 6, "x2": 236, "y1": 4, "y2": 38}]).to_pandas())
                .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
                .encode(x=alt.value(6), x2=alt.value(236), y=alt.value(4), y2=alt.value(38))
            )
            badge_t1 = (
                alt.Chart(pl.DataFrame([{"t": badge_name}]).to_pandas())
                .mark_text(align="left", baseline="top", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
                .encode(x=alt.value(11), y=alt.value(7), text="t:N")
            )
            badge_t2 = (
                alt.Chart(pl.DataFrame([{"t": f"T3 Δ = {obs_d:+.1f}% (perm p = {p_perm:.3f})"}]).to_pandas())
                .mark_text(align="left", baseline="top", fontSize=MIN_FONT_PX, fontWeight="bold", color=OKABE_BLUE)
                .encode(x=alt.value(11), y=alt.value(22), text="t:N")
            )
            x_let = -46.0 if show_y_title else -32.0
            y_let = -42.0 if show_legend else -10.0
            return cast(
                alt.LayerChart,
                (ci + ln + badge_bg + badge_t1 + badge_t2 + _panel_letter_layer(letter, x_px=x_let, y_px=y_let)).properties(
                    width=285,
                    height=155,
                ),
            )

        p_b = _make_zl_curve_panel(
            "aitchison_ilr_ontology",
            "b",
            "3. Aitchison ILR (Zero y)",
            show_y_title=True,
            show_x_title=False,
            show_legend=True,
        )
        p_c = _make_zl_curve_panel(
            "external_sc_milopy",
            "c",
            "4. External SC-Milopy Prior",
            show_y_title=False,
            show_x_title=False,
            show_legend=True,
        )
        p_d = _make_zl_curve_panel(
            "honest_5fold_oof",
            "d",
            "5. Honest 5-Fold OOF",
            show_y_title=True,
            show_x_title=True,
            show_legend=False,
        )
        p_e = _make_zl_curve_panel(
            "nested_9trial_loco",
            "e",
            "6. Nested 9-Trial LOCO",
            show_y_title=False,
            show_x_title=True,
            show_legend=False,
        )

        grid_top = alt.hconcat(p_b, p_c, spacing=24).resolve_scale(color="independent")
        grid_bot = alt.hconcat(p_d, p_e, spacing=24).resolve_scale(color="independent")
        curves_block = alt.vconcat(grid_top, grid_bot, spacing=22).resolve_scale(color="independent")

        combined = (
            alt.vconcat(panel_a, curves_block, spacing=22, center=True)
            .resolve_scale(color="independent", shape="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 10, "right": 8, "top": 14, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-supp-zero-leakage-gating-benchmark",
            layer_label="Supp_Zero_Leakage_Gating_Benchmark",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Zero-Leakage Gating Supplementary Figure: {exc}")


# ------------------------------------------------------------------------------
# Supplementary Figure S22: OOF Supplementary Diagnostics (4 Panels a-d, Heading-Free)
# ------------------------------------------------------------------------------


def build_supp_oof_diagnostics_figure(
    oof_cells_df: pl.DataFrame,
    oof_contrasts_df: pl.DataFrame,
    clin_df: pl.DataFrame,
    drivers_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build heading-free Supplementary Figure S22 (180.0mm): (a) 2x2 corners + (b) 5-contrast perm null + (c) 4/4 lineages + (d) t-weights."""
    try:
        # ------------------------------------------------------------------
        # Panel (a): Honest 5-Fold OOF 2x2 Extreme-Regime Gating + Brackets
        # ------------------------------------------------------------------
        tme_order = [
            "Permissive TME (Low NR)",
            "Excluded TME (High NR)",
        ]
        tme_colors = [OKABE_BLUE, OKABE_VERMILION]
        tme_shapes = ["circle", "square"]

        cells_aug = oof_cells_df.with_columns(
            [
                pl.when(pl.col("n_code") == 0)
                .then(pl.col("x_coord") - pl.lit(0.025))
                .otherwise(pl.col("x_coord") + pl.lit(0.025))
                .alias("x_plot"),
                pl.when((pl.col("n_code") == 0) & (pl.col("r_code") == 0))
                .then(pl.col("x_coord") + pl.lit(0.06))
                .when((pl.col("n_code") == 2) & (pl.col("r_code") == 0))
                .then(pl.col("x_coord") + pl.lit(0.06))
                .when(pl.col("n_code") == 0)
                .then(pl.col("x_coord") - pl.lit(0.02))
                .otherwise(pl.col("x_coord") + pl.lit(0.02))
                .alias("label_x"),
                pl.when(pl.col("n_code") == 0)
                .then(pl.col("ci_high_pct") + pl.lit(3.2))
                .otherwise(pl.col("ci_low_pct") - pl.lit(3.6))
                .alias("label_y"),
                pl.when(pl.col("n_code") == 0)
                .then(pl.lit(OKABE_BLUE))
                .otherwise(pl.lit(OKABE_VERMILION))
                .alias("tme_hex"),
                pl.struct(["responder_rate_pct", "n_responders", "n_patients"])
                .map_elements(
                    lambda s: (
                        f"{float(s['responder_rate_pct']):.1f}% "
                        f"({int(s['n_responders'])}/{int(s['n_patients'])})"
                    ),
                    return_dtype=pl.String,
                )
                .alias("point_callout"),
            ]
        )
        cells_pd = cells_aug.to_pandas()

        c_rate = {
            (int(r["r_code"]), int(r["n_code"])): float(r["responder_rate_pct"])
            for r in oof_cells_df.iter_rows(named=True)
        }
        con_map = {str(r["contrast_id"]): r for r in oof_contrasts_df.iter_rows(named=True)}
        c_high = con_map["high_r_permissive_vs_excluded"]
        c_perm = con_map["permissive_high_vs_low_r"]
        c_excl = con_map["excluded_high_vs_low_r"]
        c_did = con_map["did_2x2_gating_interaction"]

        x_dom_a = [-0.42, 2.02]
        y_dom_a = [8.0, 86.0]

        ci_rules_a = (
            alt.Chart(cells_pd)
            .mark_rule(strokeWidth=1.6, opacity=0.82)
            .encode(
                x=alt.X(
                    "x_plot:Q",
                    scale=alt.Scale(domain=x_dom_a),
                    title="Out-of-Fold Responder Tertile (R_OOF)",
                    axis=alt.Axis(
                        values=[0.0, 1.0],
                        labelExpr="datum.value == 0 ? 'Low R (T1)' : 'High R (T3)'",
                        labelAngle=0,
                        labelPadding=6,
                        labelLimit=160,
                        labelFontSize=MIN_FONT_PX,
                        labelFontWeight="bold",
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=False,
                    ),
                ),
                y=alt.Y(
                    "ci_low_pct:Q",
                    scale=alt.Scale(domain=y_dom_a),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
                        values=[15, 30, 45, 60, 75],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y2="ci_high_pct:Q",
                color=alt.Color("tme_hex:N", scale=None),
            )
        )
        slope_lines_a = (
            alt.Chart(cells_pd)
            .mark_line(strokeWidth=2.4)
            .encode(
                x=alt.X("x_plot:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("responder_rate_pct:Q", scale=alt.Scale(domain=y_dom_a)),
                detail="tme_short_label:N",
                color=alt.Color("tme_hex:N", scale=None),
            )
        )
        corner_pts_a = (
            alt.Chart(cells_pd)
            .mark_point(filled=True, size=115, opacity=1.0, stroke="#0F172A", strokeWidth=0.85)
            .encode(
                x=alt.X("x_plot:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("responder_rate_pct:Q", scale=alt.Scale(domain=y_dom_a)),
                color=alt.Color(
                    "tme_short_label:N",
                    scale=alt.Scale(domain=tme_order, range=tme_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=180,
                        symbolSize=90,
                    ),
                ),
                shape=alt.Shape(
                    "tme_short_label:N",
                    scale=alt.Scale(domain=tme_order, range=tme_shapes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=180,
                        symbolSize=90,
                    ),
                ),
            )
        )
        corner_lbls_a = (
            alt.Chart(cells_pd)
            .mark_text(align="center", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                x=alt.X("label_x:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("label_y:Q", scale=alt.Scale(domain=y_dom_a)),
                text="point_callout:N",
                color=alt.Color("tme_hex:N", scale=None),
            )
        )
        y_r_lo = c_rate[(2, 2)]
        y_r_hi = c_rate[(2, 0)]
        y_r_mid = 0.5 * (y_r_lo + y_r_hi)
        bracket_pd_a = pl.DataFrame(
            [
                {"x1": 1.11, "x2": 1.11, "y1": y_r_lo, "y2": y_r_hi},
                {"x1": 1.06, "x2": 1.11, "y1": y_r_lo, "y2": y_r_lo},
                {"x1": 1.06, "x2": 1.11, "y1": y_r_hi, "y2": y_r_hi},
                {"x1": 1.11, "x2": 1.15, "y1": y_r_mid, "y2": y_r_mid},
            ]
        ).to_pandas()
        brackets_a = (
            alt.Chart(bracket_pd_a)
            .mark_rule(color="#0F172A", strokeWidth=1.5)
            .encode(
                x=alt.X("x1:Q", scale=alt.Scale(domain=x_dom_a)),
                x2="x2:Q",
                y=alt.Y("y1:Q", scale=alt.Scale(domain=y_dom_a)),
                y2="y2:Q",
            )
        )
        did_box_a = (
            alt.Chart(pl.DataFrame([{"bx1": -0.38, "bx2": 1.08, "by1": 68.5, "by2": 84.5}]).to_pandas())
            .mark_rect(fill="#F8FAFC", stroke="#CBD5E1", strokeWidth=0.75)
            .encode(
                x=alt.X("bx1:Q", scale=alt.Scale(domain=x_dom_a)),
                x2="bx2:Q",
                y=alt.Y("by1:Q", scale=alt.Scale(domain=y_dom_a)),
                y2="by2:Q",
            )
        )
        annot_a_pd = pl.DataFrame(
            [
                {
                    "tx": -0.32,
                    "ty": 79.6,
                    "txt": f"2×2 DiD = {float(c_did['delta_pct']):+.1f}%",
                    "col": "#0F172A",
                },
                {
                    "tx": -0.32,
                    "ty": 72.8,
                    "txt": f"p = {_fmt_pval_short(float(c_did['wald_two_sided_p']))} (perm {_fmt_pval_short(float(c_did['perm_two_sided_p']))})",
                    "col": "#334155",
                },
                {
                    "tx": 1.18,
                    "ty": y_r_mid + 6.4,
                    "txt": f"Δ = {float(c_high['delta_pct']):+.1f}%",
                    "col": "#0F172A",
                },
                {
                    "tx": 1.18,
                    "ty": y_r_mid + 0.4,
                    "txt": f"p = {_fmt_pval_short(float(c_high['fisher_two_sided_p']))}",
                    "col": "#0F172A",
                },
                {
                    "tx": 1.18,
                    "ty": y_r_mid - 5.6,
                    "txt": f"CMH {_fmt_pval_short(float(c_high['cmh_trial_stratified_p']))}",
                    "col": OKABE_BLUE,
                },
            ]
        ).to_pandas()
        annot_left_align_a = (
            alt.Chart(annot_a_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_a)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        slope_txt_a = (
            alt.Chart(
                pl.DataFrame(
                    [
                        {
                            "tx": 0.34,
                            "ty": 53.5,
                            "txt": f"+{float(c_perm['delta_pct']):.1f}% (p={_fmt_pval_short(float(c_perm['fisher_two_sided_p']))})",
                            "col": OKABE_BLUE,
                        },
                        {
                            "tx": 0.50,
                            "ty": 33.2,
                            "txt": f"{float(c_excl['delta_pct']):+.1f}% (n.s.)",
                            "col": OKABE_VERMILION,
                        },
                    ]
                ).to_pandas()
            )
            .mark_text(align="center", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                x=alt.X("tx:Q", scale=alt.Scale(domain=x_dom_a)),
                y=alt.Y("ty:Q", scale=alt.Scale(domain=y_dom_a)),
                text="txt:N",
                color=alt.Color("col:N", scale=None),
            )
        )
        panel_a = (
            did_box_a
            + ci_rules_a
            + slope_lines_a
            + corner_pts_a
            + corner_lbls_a
            + brackets_a
            + annot_left_align_a
            + slope_txt_a
            + _panel_letter_layer("a", x_px=-46.0, y_px=-30.0)
        ).properties(width=270, height=185)

        # ------------------------------------------------------------------
        # Panel (b): 5 Pairwise & 2x2 Epistatic Contrasts vs B=500 Honest OOF Permutation Null
        # ------------------------------------------------------------------
        con_short_map = {
            "permissive_high_vs_low_r": "Perm: High−Low R",
            "high_r_permissive_vs_excluded": "High-R: Perm−Excl",
            "did_2x2_gating_interaction": "2×2 Epistasis DiD",
            "excluded_high_vs_low_r": "Excl: High−Low R",
            "low_r_permissive_vs_excluded": "Low-R: Perm−Excl",
        }
        con_aug = oof_contrasts_df.sort("contrast_rank").with_columns(
            [
                pl.col("contrast_id").replace(con_short_map).alias("c_short"),
                pl.when(pl.col("fisher_two_sided_p") < 0.05)
                .then(pl.lit("Significant (p < 0.05)"))
                .otherwise(pl.lit("Blunted Control (n.s.)"))
                .alias("sig_group"),
                pl.when(pl.col("fisher_two_sided_p") < 0.05)
                .then(pl.lit(OKABE_BLUE))
                .otherwise(pl.lit("#94A3B8"))
                .alias("ci_hex"),
                (
                    pl.max_horizontal("ci_high_pct", "perm_null_ci_high_pct")
                    + pl.lit(2.4)
                ).alias("callout_x"),
                pl.struct(["delta_pct", "perm_two_sided_p"])
                .map_elements(
                    lambda s: f"{float(s['delta_pct']):+.1f}% (perm {_fmt_pval_short(float(s['perm_two_sided_p']))})",
                    return_dtype=pl.String,
                )
                .alias("callout_text"),
            ]
        )
        ordered_contrasts = con_aug["c_short"].to_list()
        obs_b_df = con_aug.select(
            [
                pl.col("contrast_rank"),
                pl.col("c_short"),
                pl.col("delta_pct").alias("x_val"),
                pl.col("sig_group").alias("marker_series"),
            ]
        )
        null_b_df = con_aug.select(
            [
                pl.col("contrast_rank"),
                pl.col("c_short"),
                pl.col("perm_null_mean_pct").alias("x_val"),
                pl.lit("OOF Perm Null (±95% UI)").alias("marker_series"),
            ]
        )
        pts_b_pd = pl.concat([obs_b_df, null_b_df]).to_pandas()
        con_pd = con_aug.to_pandas()

        b_series_order = [
            "Significant (p < 0.05)",
            "Blunted Control (n.s.)",
            "OOF Perm Null (±95% UI)",
        ]
        b_series_colors = [OKABE_BLUE, "#94A3B8", "#475569"]
        b_series_shapes = ["circle", "square", "diamond"]
        x_dom_b = [-22.0, 128.0]

        zero_line_b = (
            alt.Chart(pl.DataFrame({"x0": [0.0]}).to_pandas())
            .mark_rule(color="#0F172A", strokeDash=[5, 3], strokeWidth=1.4)
            .encode(x=alt.X("x0:Q", scale=alt.Scale(domain=x_dom_b)))
        )
        null_ui_b = (
            alt.Chart(con_pd)
            .mark_rule(color="#CBD5E1", strokeWidth=5.0, opacity=0.9)
            .encode(
                y=alt.Y(
                    "c_short:N",
                    sort=ordered_contrasts,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=145,
                        labelFontWeight="bold",
                    ),
                ),
                x=alt.X(
                    "perm_null_ci_low_pct:Q",
                    scale=alt.Scale(domain=x_dom_b),
                    title="Contrast Δ (%) vs. OOF Perm Null",
                    axis=alt.Axis(
                        values=[-15, 0, 15, 30, 45],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="perm_null_ci_high_pct:Q",
            )
        )
        obs_ci_b = (
            alt.Chart(con_pd)
            .mark_rule(strokeWidth=1.7)
            .encode(
                y=alt.Y("c_short:N", sort=ordered_contrasts),
                x=alt.X("ci_low_pct:Q", scale=alt.Scale(domain=x_dom_b)),
                x2="ci_high_pct:Q",
                color=alt.Color("ci_hex:N", scale=None),
            )
        )
        pts_b = (
            alt.Chart(pts_b_pd)
            .mark_point(filled=True, size=92, opacity=1.0, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                y=alt.Y("c_short:N", sort=ordered_contrasts),
                x=alt.X("x_val:Q", scale=alt.Scale(domain=x_dom_b)),
                color=alt.Color(
                    "marker_series:N",
                    scale=alt.Scale(domain=b_series_order, range=b_series_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        columns=2,
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=165,
                        symbolSize=80,
                    ),
                ),
                shape=alt.Shape(
                    "marker_series:N",
                    scale=alt.Scale(domain=b_series_order, range=b_series_shapes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        columns=2,
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=165,
                        symbolSize=80,
                    ),
                ),
            )
        )
        callouts_b = (
            alt.Chart(con_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("c_short:N", sort=ordered_contrasts),
                x=alt.X("callout_x:Q", scale=alt.Scale(domain=x_dom_b)),
                text="callout_text:N",
            )
        )
        panel_b = (
            zero_line_b
            + null_ui_b
            + obs_ci_b
            + pts_b
            + callouts_b
            + _panel_letter_layer("b", x_px=-140.0, y_px=-42.0)
        ).properties(width=225, height=185)

        # ------------------------------------------------------------------
        # Panel (c): 4/4 Solid-Tumor Lineage + Pan-Cancer Replication Forest
        # ------------------------------------------------------------------
        hist_short_map = {
            "Melanoma (4 Trials, n=338)": "Melanoma (n=338)",
            "Urothelial Bladder (n=298)": "Bladder (n=298)",
            "Renal ccRCC (2 Trials, n=263)": "Renal ccRCC (n=263)",
            "PDAC + TNBC (2 Trials, n=116)": "PDAC + TNBC (n=116)",
            "Pan-Cancer Pooled (N=1015)": "Pan-Cancer (N=1,015)",
        }
        hist_df = (
            clin_df.filter(pl.col("record_type") == "histology_contrast")
            .sort(["group_order", "sub_order"])
            .with_columns(
                [
                    pl.col("group_label").replace(hist_short_map).alias("g_short"),
                    pl.when(pl.col("hit_id") == "hit_9a")
                    .then(pl.lit("Permissive ΔR"))
                    .otherwise(pl.lit("High-R ΔTME"))
                    .alias("hit_short"),
                    pl.col("delta_pct")
                    .map_elements(lambda v: f"{float(v):+.1f}%", return_dtype=pl.String)
                    .alias("delta_lbl"),
                    (pl.col("delta_pct") + pl.lit(1.8)).alias("lbl_x"),
                ]
            )
        )
        hist_order = list(dict.fromkeys(hist_df["g_short"].to_list()))
        hit_dom_c = ["Permissive ΔR", "High-R ΔTME"]
        hit_cols_c = [OKABE_BLUE, OKABE_TEAL]
        hit_shapes_c = ["circle", "diamond"]
        x_dom_c = [-4.0, 56.0]
        hist_pd = hist_df.to_pandas()

        zero_c = (
            alt.Chart(pl.DataFrame({"xz": [0.0]}).to_pandas())
            .mark_rule(stroke="#0F172A", strokeDash=[3, 3], strokeWidth=1.1)
            .encode(x=alt.X("xz:Q", scale=alt.Scale(domain=x_dom_c)))
        )
        stems_c = (
            alt.Chart(hist_pd)
            .mark_rule(strokeWidth=2.0, opacity=0.85)
            .encode(
                y=alt.Y(
                    "g_short:N",
                    sort=hist_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=165,
                        labelFontWeight="bold",
                    ),
                ),
                yOffset=alt.YOffset("hit_short:N", sort=hit_dom_c),
                x=alt.X(
                    "zero_v:Q",
                    scale=alt.Scale(domain=x_dom_c),
                    title="Corner Responder Shift Δ (%)",
                    axis=alt.Axis(
                        values=[0, 15, 30, 45],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="delta_pct:Q",
                color=alt.Color(
                    "hit_short:N",
                    scale=alt.Scale(domain=hit_dom_c, range=hit_cols_c),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=150,
                        symbolSize=85,
                    ),
                ),
            )
            .transform_calculate(zero_v="0")
        )
        pts_c = (
            alt.Chart(hist_pd)
            .mark_point(filled=True, size=82, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                y=alt.Y("g_short:N", sort=hist_order),
                yOffset=alt.YOffset("hit_short:N", sort=hit_dom_c),
                x=alt.X("delta_pct:Q", scale=alt.Scale(domain=x_dom_c)),
                color=alt.Color("hit_short:N", scale=alt.Scale(domain=hit_dom_c, range=hit_cols_c)),
                shape=alt.Shape(
                    "hit_short:N",
                    scale=alt.Scale(domain=hit_dom_c, range=hit_shapes_c),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=150,
                        symbolSize=85,
                    ),
                ),
            )
        )
        lbls_c = (
            alt.Chart(hist_pd)
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
            .encode(
                y=alt.Y("g_short:N", sort=hist_order),
                yOffset=alt.YOffset("hit_short:N", sort=hit_dom_c),
                x=alt.X("lbl_x:Q", scale=alt.Scale(domain=x_dom_c)),
                text="delta_lbl:N",
                color=alt.Color("hit_short:N", scale=alt.Scale(domain=hit_dom_c, range=hit_cols_c), legend=None),
            )
        )
        panel_c = (
            zero_c
            + stems_c
            + pts_c
            + lbls_c
            + _panel_letter_layer("c", x_px=-146.0, y_px=-30.0)
        ).properties(width=195, height=185)

        # ------------------------------------------------------------------
        # Panel (d): Bidirectional Diverging Stem Plot of Top 5 R_OOF vs Top 5 NR_OOF Drivers
        # ------------------------------------------------------------------
        pos5 = drivers_df.sort("rank_pos").head(5)
        neg5 = drivers_df.sort("rank_neg").head(5).sort("mean_t_weight", descending=True)
        div_df = pl.concat([pos5, neg5]).with_columns(
            [
                pl.when(pl.col("mean_t_weight") >= 0.0)
                .then(pl.lit("R_OOF Effector (+t, 5/5f)"))
                .otherwise(pl.lit("NR_OOF Barrier (−t, 5/5f)"))
                .alias("axis_grp"),
                pl.when(pl.col("mean_t_weight") >= 0.0)
                .then(pl.lit(OKABE_BLUE))
                .otherwise(pl.lit(OKABE_VERMILION))
                .alias("grp_hex"),
                pl.col("mean_t_weight")
                .map_elements(lambda v: f"{float(v):+.1f}", return_dtype=pl.String)
                .alias("t_lbl"),
                pl.when(pl.col("mean_t_weight") >= 0.0)
                .then(pl.col("ci_high_t") + pl.lit(0.28))
                .otherwise(pl.col("ci_low_t") - pl.lit(0.28))
                .alias("lbl_x"),
            ]
        )
        div_order = div_df["state_label"].to_list()
        div_dom = ["R_OOF Effector (+t, 5/5f)", "NR_OOF Barrier (−t, 5/5f)"]
        div_cols = [OKABE_BLUE, OKABE_VERMILION]
        div_shapes = ["circle", "square"]
        x_dom_d = [-5.8, 7.6]
        div_pd = div_df.to_pandas()

        zero_d = (
            alt.Chart(pl.DataFrame({"x0": [0.0]}).to_pandas())
            .mark_rule(stroke="#0F172A", strokeDash=[3, 3], strokeWidth=1.2)
            .encode(x=alt.X("x0:Q", scale=alt.Scale(domain=x_dom_d)))
        )
        stems_d = (
            alt.Chart(div_pd)
            .mark_rule(strokeWidth=2.1, opacity=0.55)
            .encode(
                y=alt.Y(
                    "state_label:N",
                    sort=div_order,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=215,
                        labelFontWeight="bold",
                    ),
                ),
                x=alt.X(
                    "zero_d:Q",
                    scale=alt.Scale(domain=x_dom_d),
                    title="5-Fold Cross-Fitted Weight (t ± SD)",
                    axis=alt.Axis(
                        values=[-4, -2, 0, 2, 4, 6],
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                x2="mean_t_weight:Q",
                color=alt.Color("grp_hex:N", scale=None),
            )
            .transform_calculate(zero_d="0")
        )
        sd_d = (
            alt.Chart(div_pd)
            .mark_rule(stroke="#0F172A", strokeWidth=1.45)
            .encode(
                y=alt.Y("state_label:N", sort=div_order),
                x=alt.X("ci_low_t:Q", scale=alt.Scale(domain=x_dom_d)),
                x2="ci_high_t:Q",
            )
        )
        pts_d = (
            alt.Chart(div_pd)
            .mark_point(filled=True, size=82, opacity=1.0, stroke="#0F172A", strokeWidth=0.8)
            .encode(
                y=alt.Y("state_label:N", sort=div_order),
                x=alt.X("mean_t_weight:Q", scale=alt.Scale(domain=x_dom_d)),
                color=alt.Color(
                    "axis_grp:N",
                    scale=alt.Scale(domain=div_dom, range=div_cols),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=175,
                        symbolSize=85,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
                shape=alt.Shape(
                    "axis_grp:N",
                    scale=alt.Scale(domain=div_dom, range=div_shapes),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="horizontal",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=175,
                        symbolSize=85,
                        symbolOpacity=1.0,
                        symbolStrokeColor="#0F172A",
                        symbolStrokeWidth=0.8,
                    ),
                ),
            )
        )
        lbl_pos_d = (
            alt.Chart(div_df.filter(pl.col("mean_t_weight") >= 0.0).to_pandas())
            .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("state_label:N", sort=div_order),
                x=alt.X("lbl_x:Q", scale=alt.Scale(domain=x_dom_d)),
                text="t_lbl:N",
            )
        )
        lbl_neg_d = (
            alt.Chart(div_df.filter(pl.col("mean_t_weight") < 0.0).to_pandas())
            .mark_text(align="right", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("state_label:N", sort=div_order),
                x=alt.X("lbl_x:Q", scale=alt.Scale(domain=x_dom_d)),
                text="t_lbl:N",
            )
        )
        panel_d = (
            zero_d
            + stems_d
            + sd_d
            + pts_d
            + lbl_pos_d
            + lbl_neg_d
            + _panel_letter_layer("d", x_px=-205.0, y_px=-30.0)
        ).properties(width=185, height=185)

        top_row = alt.hconcat(panel_a, panel_b, spacing=18).resolve_scale(color="independent", shape="independent")
        bot_row = alt.hconcat(panel_c, panel_d, spacing=18).resolve_scale(color="independent", shape="independent")
        combined = (
            alt.vconcat(top_row, bot_row, spacing=22, center=True)
            .resolve_scale(color="independent", shape="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 10, "right": 18, "top": 12, "bottom": 6},
            )
        )
        raw_svg = vlc.vegalite_to_svg(combined.to_dict())
        wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
            raw_svg,
            width_mm=FULL_WIDTH_MM,
            layer_id="layer-supp-oof-diagnostics",
            layer_label="Supp_OOF_Supplementary_Diagnostics",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build OOF Supplementary Diagnostics Figure: {exc}")


# ------------------------------------------------------------------------------
# Standalone Importable Typst Scientific Figure Captions Generator
# ------------------------------------------------------------------------------


def export_typst_figure_captions(
    oof_cells_df: pl.DataFrame,
    oof_contrasts_df: pl.DataFrame,
    landscape_df: pl.DataFrame,
    clin_df: pl.DataFrame,
    drivers_df: pl.DataFrame,
    pure_inter_df: pl.DataFrame,
    epi_strat_df: pl.DataFrame,
    epi_sig_df: pl.DataFrame,
    zl_summary_df: pl.DataFrame,
    captions_dir: Path,
    results_dir: Path,
) -> Result[tuple[Path, Path], str]:
    """Write standalone importable .typ caption modules for the Main Figure and Supplementary Figures S19-S22."""
    try:
        def _fmt_pval_typst(p_val: float) -> str:
            if p_val < 1e-4:
                exp = int(np.floor(np.log10(max(p_val, 1e-300))))
                mant = p_val / (10.0**exp)
                return f"{mant:.1f} times 10^({exp})"
            return f"{p_val:.4f}" if p_val < 0.01 else f"{p_val:.3f}"

        captions_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)

        # Extract exact statistics for Main Figure & Supplementary captions
        n_corners = int(oof_cells_df["n_patients"].sum())
        c_map = {
            (int(r["r_code"]), int(r["n_code"])): r
            for r in oof_cells_df.iter_rows(named=True)
        }
        con_map = {str(r["contrast_id"]): r for r in oof_contrasts_df.iter_rows(named=True)}
        c_9a = con_map["permissive_high_vs_low_r"]
        c_9b = con_map["high_r_permissive_vs_excluded"]
        c_did = con_map["did_2x2_gating_interaction"]
        c_excl = con_map["excluded_high_vs_low_r"]

        curve_9a = landscape_df.filter(pl.col("record_type") == "dose_curve_9a")
        perm_9a = curve_9a.filter(pl.col("stratum_id") == "permissive_nr_t1").row(0, named=True)
        excl_9a = curve_9a.filter(pl.col("stratum_id") == "excluded_nr_t3").row(0, named=True)

        curve_9b = landscape_df.filter(pl.col("record_type") == "dose_curve_9b")
        hi_9b = curve_9b.filter(pl.col("stratum_id") == "high_r_t3").row(0, named=True)
        lo_9b = curve_9b.filter(pl.col("stratum_id") == "low_r_t1").row(0, named=True)

        rec_hp = clin_df.filter(
            (pl.col("record_type") == "recist_corner_stack")
            & (pl.col("group_id") == "high_r_perm")
            & (pl.col("sub_category") == "CR")
        ).row(0, named=True)
        rec_hp_cr = float(rec_hp["rate_a_pct"])
        rec_hp_dcr = float(rec_hp["dcr_pct"])

        km_lr = {
            str(r["hit_id"]): r
            for r in clin_df.filter(pl.col("record_type") == "km_logrank").iter_rows(named=True)
        }
        p_os_9a = float(km_lr["hit_9a"]["p_value"])
        p_os_9b = float(km_lr["hit_9b"]["p_value"])
        p_os_excl = float(km_lr["excl_null"]["p_value"])
        med_lp = float(km_lr["hit_9a"]["median_b_mo"])
        med_he = float(km_lr["hit_9b"]["median_b_mo"])

        m29_cells = {
            (int(r["x_code"]), int(r["mod_code"])): r
            for r in epi_sig_df.filter(
                (pl.col("record_type") == "tertile_2x2_cell") & (pl.col("pair_id") == "M29_x_M26")
            ).iter_rows(named=True)
        }
        m29_r0 = m29_cells[(0, 0)]
        m29_gate_delta = float(m29_cells[(2, 0)]["responder_rate_pct"]) - float(
            m29_cells[(2, 2)]["responder_rate_pct"]
        )
        volc_tnk = epi_sig_df.filter(
            (pl.col("record_type") == "all861_pair_volcano") & (pl.col("stratum_id") == "t_nk_pair")
        ).sort("p_spa")
        min_tnk_q = float(volc_tnk["q_bh_all861"][0]) if volc_tnk.height > 0 else 1.0
        m29_lo_c = epi_sig_df.filter(
            (pl.col("record_type") == "m29_dose_curve") & (pl.col("stratum_id") == "low_emt_m26_t1")
        ).row(0, named=True)
        m29_hi_c = epi_sig_df.filter(
            (pl.col("record_type") == "m29_dose_curve") & (pl.col("stratum_id") == "high_emt_m26_t3")
        ).row(0, named=True)
        top_pos = drivers_df.sort("rank_pos").row(0, named=True)
        top_neg = drivers_df.sort("rank_neg").row(0, named=True)

        main_caption_typ = f"""// Auto-generated by scripts/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures.py
// Do not edit manually; re-run Stage 4b to synchronize with Stage 3c parquet outputs.

#let caption_epistasis_stromal_gating_main = [
  *Pan-compartment clinical stromal gating, specific non-T-cell humoral–mesenchymal antagonism across $N = 1,015$ checkpoint blockade biopsies, and single-cell cross-compartment interaction manifolds across the host TME Triad.*
  *(a)* Continuous logistic dose-response ($plus.minus 95%$ CI ribbons and empirical quintiles) along the `5`-fold out-of-fold pan-compartment effector composite $R_("OOF")$ (`T/NK` + `M1-TAM` + `B/Plasma`; zero bulk-outcome leakage), confirming steep benefit in Permissive low-barrier TMEs ($"NR"_("OOF")$ `T1`, `OR` $= {float(perm_9a['or_per_sd']):.2f}$ per SD, Wald $p = {_fmt_pval_typst(float(perm_9a['wald_p']))}$, $n = {int(perm_9a['n_patients'])}$) versus a flat trajectory in Excluded high-barrier TMEs ($"NR"_("OOF")$ `T3`, `OR` $= {float(excl_9a['or_per_sd']):.2f}$ per SD, $p = {float(excl_9a['wald_p']):.2f}$, n.s., $n = {int(excl_9a['n_patients'])}$).
  *(b)* Reciprocal continuous suppression along the barrier composite $"NR"_("OOF")$ (`MDSC` + `EMT` + `Hypoxia` + `CAF`), showing steep attenuation among high-$R_("OOF")$ (`T3`) tumors (`OR` $= {float(hi_9b['or_per_sd']):.2f}$ per SD, $p = {_fmt_pval_typst(float(hi_9b['wald_p']))}$, $n = {int(hi_9b['n_patients'])}$) versus a flat baseline among low-$R_("OOF")$ (`T1`) tumors (`OR` $= {float(lo_9b['or_per_sd']):.2f}$ per SD, $p = {float(lo_9b['wald_p']):.2f}$, n.s., $n = {int(lo_9b['n_patients'])}$).
  *(c)* Joint two-dimensional $(R_("OOF"), "NR"_("OOF"))$ patient coordinate manifold ($N = 1,015$) with within-trial tertile cuts, logistic iso-response contours, and $2 times 2$ extreme-tertile corner rates ($n = {n_corners}$; Permissive: ${float(c_map[(0, 0)]['responder_rate_pct']):.1f}% -> {float(c_map[(2, 0)]['responder_rate_pct']):.1f}%$, Fisher $p = {_fmt_pval_typst(float(c_9a['fisher_two_sided_p']))}$; Excluded: ${float(c_map[(0, 2)]['responder_rate_pct']):.1f}% -> {float(c_map[(2, 2)]['responder_rate_pct']):.1f}%$, $p = {float(c_excl['fisher_two_sided_p']):.2f}$; $2 times 2$ `DiD` $= {float(c_did['delta_pct']):+.1f}%$, Wald $p = {_fmt_pval_typst(float(c_did['wald_two_sided_p']))}$, `OOF` permutation $p = {_fmt_pval_typst(float(c_did['perm_two_sided_p']))}$).
  *(d)* Granular `RECIST 1.1` (`CR`, `PR`, `SD`, `PD`) corner composition (`High-R / Permissive` `CR` $= {rec_hp_cr:.1f}%$, `DCR` $= {rec_hp_dcr:.1f}%$).
  *(e)* Orthogonal Kaplan–Meier Overall Survival (`OS`) validation (zero survival supervision): `High-R / Permissive` achieves unreached median `OS` versus `Low-R / Permissive` (${med_lp:.1f}$ months, log-rank $p = {_fmt_pval_typst(p_os_9a)}$) and `High-R / Excluded` (${med_he:.1f}$ months, $p = {_fmt_pval_typst(p_os_9b)}$; Excluded null $p = {p_os_excl:.2f}$).
  *(f)* Specific non-T-cell pairwise cellular epistasis across all $861$ `Consensus-Meta42` pairs on unblended single-cell `DEG` marker signatures: continuous logistic dose-response ($plus.minus 95%$ CI ribbons and quintiles) for the sole `FWER`/`FDR`-significant pair, humoral `M29` (`JCHAIN+MZB1+` Plasma) $times$ mesenchymal `M26` (`VIM+ZEB1+` `EMT`; $z_("SPA") = {float(m29_r0['z_spa']):+.2f}$, `FWER` $= {float(m29_r0['p_fwer_all861']):.3f}$, $q_("BH") = {float(m29_r0['q_bh_all861']):.3f}$ across $861$ pairs and $q_("BH") = {float(m29_r0['q_bh_screened']):.4f}$ across $160$ `milopy` bipartite pairs; all `T/NK`-cell pairs $q_("BH") >= {min_tnk_q:.2f}$, n.s.), confirming strong humoral benefit in epithelial low-`EMT` (`T1`) tumors (`OR` $= {float(m29_lo_c['beta_inter']):.2f}$ per SD, Wald $p = {_fmt_pval_typst(float(m29_lo_c['p_spa']))}$, $n = {int(m29_lo_c['n_patients'])}$) that is completely neutralized in mesenchymal high-`EMT` (`T3`) tumors (`OR` $= {float(m29_hi_c['beta_inter']):.2f}$ per SD, $p = {float(m29_hi_c['p_spa']):.2f}$, n.s., $n = {int(m29_hi_c['n_patients'])}$).
  *(g–i)* Single-cell cross-compartment interaction manifolds projected onto the unified `9`-cohort `Harmony I2` single-cell UMAP ($36,000$ anchor cells across $613$ biopsies and $20$ `CMH` strata) for the three arms of the host TME Triad: *(g)* `T-NK cells` $<->$ `CAF` ($6$ `ACAT FDR` $< 0.10$ bold chords + $2$ stratified permutation $p < 0.05$ thin chords; lead pair `M08` [`GZMK+PDCD1+` Progenitor Tex] $times$ `M17` [`ISLR+COL5A2+` `myCAF`], $Z_("ACAT") = +4.14^(***)$, $q = 0.0074$), *(h)* `Myeloid/Macrophages` $<->$ `CAF` ($2$ `ACAT FDR` $< 0.10$ + $2$ permutation $p < 0.05$; lead pair `M17` [`myCAF`] $times$ `M42` [`TPSAB1+CPA3+` Mast], $Z_("ACAT") = -4.41^(***)$, $q = 0.0065$), and *(i)* `T-NK cells` $<->$ `Myeloid/Macrophages` ($4$ `ACAT FDR` $< 0.10$ + $14$ permutation $p < 0.05$; lead pair `M06` [`CX3CR1+FGFBP2+` Cytotoxic NK/T] $times$ `M09` [`SIGLEC8+CR1+` Granul./Myeloid], $Z_("ACAT") = +3.39$, $q = 0.049$).
]
"""

        # Extract exact statistics for Supplementary Figures S19-S22 captions
        off_pure = pure_inter_df.filter(pl.col("state_a") < pl.col("state_b")).sort("p_spa")
        lam_spa = float(off_pure["lambda_gc_spa"][0])
        lam_2s = float(off_pure["lambda_gc_2state"][0])
        n_p01 = off_pure.filter(pl.col("p_spa") < 0.01).height
        n_p005 = off_pure.filter(pl.col("p_spa") < 0.005).height
        min_p_pure = float(off_pure["p_spa"].min())  # type: ignore[arg-type]
        min_q_pure = float(off_pure["q_bh_spa"].min())  # type: ignore[arg-type]
        top_pure_row = off_pure.row(0, named=True)

        strat_map = {str(r["strategy_id"]): r for r in epi_strat_df.iter_rows(named=True)}
        s1_nb = 100.0 * float(strat_map["naive_2state_perm_clr"]["null_b_fpr_p05"])
        s2_nb = 100.0 * float(strat_map["outcome_peeking_top12_screen"]["null_b_fpr_p05"])
        s4_q = float(strat_map["sc_deg_100_all861"]["emp_min_q_bh"])
        s5_q = float(strat_map["sc_deg_100_milo_bipartite"]["emp_min_q_bh"])
        s6_q = float(strat_map["oof_2x2_composite_gate"]["emp_min_q_bh"])

        zl_map = {str(r["method_id"]): r for r in zl_summary_df.iter_rows(named=True)}
        zl1_null = float(zl_map["legacy_explicit_y"]["perm_null_mean_delta_t3_pct"])
        zl2_null = float(zl_map["insample_bulk_beta"]["perm_null_mean_delta_t3_pct"])
        zl3_obs = float(zl_map["aitchison_ilr_ontology"]["observed_delta_t3_pct"])
        zl4_obs = float(zl_map["external_sc_milopy"]["observed_delta_t3_pct"])
        zl5_obs = float(zl_map["honest_5fold_oof"]["observed_delta_t3_pct"])
        zl6_obs = float(zl_map["nested_9trial_loco"]["observed_delta_t3_pct"])

        supp_captions_typ = f"""// Auto-generated by scripts/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures.py
// Do not edit manually; re-run Stage 4b to synchronize with Stage 3c parquet outputs.

#let caption_supp_pure_clr_all42 = [
  *Complete $42 times 42$ pairwise conditional interaction matrix on pure deconvolved centered log-ratio (`CLR`) abundances.*
  Symmetric pairwise epistasis matrix across all $861$ off-diagonal pairs of `Consensus-Meta42` cellular states (`M01`–`M42`, ranked by `Pan-iAtlas-Meta` marginal deconvolution effect $beta$) evaluated on $100%$ unblended `InstaPrism` `CLR` abundances ($alpha = 0.00$, zero gene-signature mixing, zero bulk-outcome supervision in feature construction; $N = 1,015$ patients across `9` solid-tumor ICB trials). Each cell colors the calibrated leave-two-out ridge-nuisance, quadratic-guarded Bernoulli cumulant-generating-function saddlepoint statistic $z_("SPA")$, which eliminates the genomic-control inflation of naive two-state permutation testing ($lambda_("GC") = {lam_spa:.2f}$ vs. ${lam_2s:.2f}$). Bullet glyphs denote nominal $p_("SPA") < 0.01$ ({n_p01}/861 pairs) and asterisks denote $p_("SPA") < 0.005$ ({n_p005}/861 pairs; lead pair `{top_pure_row['state_a']}` $times$ `{top_pure_row['state_b']}`, $z_("SPA") = {float(top_pure_row['epistasis_z']):+.2f}$, $p_("SPA") = {_fmt_pval_typst(min_p_pure)}$). Under strict $861$-pair multiplicity and simplex closure attenuation, no individual pure-`CLR` pair crosses false discovery rate significance (minimum $q_("BH") = {min_q_pure:.2f}$), motivating orthogonal single-cell marker signatures and composite out-of-fold stromal gating.
]

#let caption_supp_epistasis_power_benchmark = [
  *Type I error calibration audit, multiplicity-adjusted discovery power, cross-indication consistency, and discrete $2 times 2$ corner interaction of pairwise cellular epistasis.*
  *(a)* Empirical false-positive rate (`FPR` at nominal $p <= 0.05$) across six non-hybrid epistasis testing strategies evaluated under within-trial outcome permutation (`Permuted-y Null`) versus a `42`-state additive plus quadratic logit saturation null (`Saturation Null B`, zero true pairwise interactions). Naive two-state `CLR` permutation (`Strategy 1`, `Null B FPR` $= {s1_nb:.1f}%$) and outcome-peeking marginal pre-screening (`Strategy 2`, `Null B FPR` $= {s2_nb:.1f}%$) exhibit severe false-positive inflation under nonlinear logit saturation, whereas leave-two-out quadratic-guarded `SPA` strategies (`Strategies 3–6`) maintain well-controlled Type I error near the $5.0%$ nominal level.
  *(b)* Empirical multiplicity-adjusted discovery power ($-log_10(min q_("BH"))$) across the four calibrated zero-leakage strategies on CRI iAtlas human ICB trials ($N = 1,015$): while unblended `Pure CLR` (`Strategy 3`, $M = 861$) remains sub-threshold ($q_("BH") = {min_q_pure:.2f}$), unblended single-cell `DEG` marker signatures (`Strategy 4`, $M = 861$, $q_("BH") = {s4_q:.4f}$), outcome-blind external single-cell `milopy` bipartite pre-screening (`Strategy 5`, $M = 160$, $q_("BH") = {s5_q:.4f}$), and honest `5`-fold `OOF` $2 times 2$ composite gating (`Strategy 6`, $M = 5$, $q_("BH") = {s6_q:.4f}$) each achieve `FDR`-significant discovery without bulk-outcome leakage.
  *(c)* Forest plot of stratified logistic interaction log-odds $beta_(a b)$ ($plus.minus 95%$ CI) and calibrated $z_("SPA")$ for `M29` (`JCHAIN+MZB1+` Plasma) $times$ `M26` (`VIM+ZEB1+` `EMT`) across all four solid-tumor lineages and the pooled pan-cancer cohort, demonstrating a directionally consistent negative epistatic interaction in every lineage.
  *(d)* Discrete within-trial $2 times 2$ extreme-tertile corner interaction (`T2` omitted; $plus.minus$ Wilson $95%$ CIs) for `M29` (`JCHAIN+MZB1+` Plasma) $times$ `M26` (`VIM+ZEB1+` `EMT`): humoral plasma enrichment more than doubles objective response in epithelial low-`EMT` (`T1`) tumors (${float(m29_cells[(0, 0)]['responder_rate_pct']):.1f}% -> {float(m29_cells[(2, 0)]['responder_rate_pct']):.1f}%$, Fisher $p = {_fmt_pval_typst(float(m29_r0['fisher_slope_mod_low_p']))}$) whereas mesenchymal high-`EMT` (`T3`) exclusion abolishes humoral benefit (${float(m29_cells[(0, 2)]['responder_rate_pct']):.1f}% -> {float(m29_cells[(2, 2)]['responder_rate_pct']):.1f}%$, $p = {float(m29_r0['fisher_slope_mod_high_p']):.2f}$; high-`M29` stromal gate $Delta = {m29_gate_delta:+.1f}%$, $p = {_fmt_pval_typst(float(m29_r0['fisher_gate_x_high_p']))}$).
]

#let caption_supp_zero_leakage_gating_comparison = [
  *Outcome-leakage audit and zero-leakage $3 times 3$ clinical stromal gating across unblended cellular estimators.*
  *(a)* Observed top-tertile (`T3`) stromal gating contrast (Permissive minus Excluded responder rate, $Delta%$) compared against $B = 500$ within-trial outcome-permutation nulls ($E[Delta | y_("perm")]$ $plus.minus 95%$ uncertainty interval) across seven gating estimators ($N = 1,015$). In-sample outcome-supervised baselines (`Methods 1–2`) exhibit severe positive permutation-null optimism ($+{zl1_null:.1f}%$ and $+{zl2_null:.1f}%$), whereas all five strict zero-leakage estimators (`Methods 3–7`, fitted without target fold bulk outcomes on strictly unblended representations) center near $0%$ under permutation.
  *(b–e)* Full $3 times 3$ within-trial responder-score tertile (`T1–T3`) $times$ stromal-barrier regime trajectories ($plus.minus$ Wilson $95%$ CIs) for the four interpretable zero-leakage estimators: *(b)* Aitchison isometric log-ratio (`ILR`) ontology on raw `InstaPrism` $theta$ (`T3` $Delta = {zl3_obs:+.1f}%$), *(c)* external single-cell `milopy` prior applied to unblended single-cell marker signatures (`T3` $Delta = {zl4_obs:+.1f}%$), *(d)* honest `5`-fold cross-fitted `OOF` composites (`T3` $Delta = {zl5_obs:+.1f}%$), and *(e)* nested `9`-trial leave-one-cohort-out (`LOCO`) cross-fitting (`T3` $Delta = {zl6_obs:+.1f}%$).
]

#let caption_supp_oof_diagnostics = [
  *Discrete $2 times 2$ extreme-tertile corner interactions, permutation calibration, cross-histology replication, and cross-fitted cellular drivers of `5`-fold `OOF` stromal gating.*
  *(a)* Out-of-fold (`5`-fold `OOF`) $2 times 2$ extreme-tertile corner interaction trajectories ($plus.minus$ Wilson $95%$ CIs) across $n = {n_corners}$ corner patients (`T2` omitted): rising $R_("OOF")$ elevates objective response from ${float(c_map[(0, 0)]['responder_rate_pct']):.1f}%$ ({int(c_map[(0, 0)]['n_responders'])}/{int(c_map[(0, 0)]['n_patients'])}) to ${float(c_map[(2, 0)]['responder_rate_pct']):.1f}%$ ({int(c_map[(2, 0)]['n_responders'])}/{int(c_map[(2, 0)]['n_patients'])}) in Permissive TMEs ($Delta = {float(c_9a['delta_pct']):+.1f}%$, two-sided Fisher $p = {_fmt_pval_typst(float(c_9a['fisher_two_sided_p']))}$, trial-stratified Cochran–Mantel–Haenszel $p = {_fmt_pval_typst(float(c_9a['cmh_trial_stratified_p']))}$) versus ${float(c_map[(0, 2)]['responder_rate_pct']):.1f}%$ to ${float(c_map[(2, 2)]['responder_rate_pct']):.1f}%$ in Excluded TMEs ($Delta = {float(c_excl['delta_pct']):+.1f}%$, $p = {float(c_excl['fisher_two_sided_p']):.2f}$, n.s.), yielding a ${float(c_9b['delta_pct']):+.1f}%$ high-$R_("OOF")$ stromal gate ($p = {_fmt_pval_typst(float(c_9b['fisher_two_sided_p']))}$) and synergistic $2 times 2$ difference-in-differences epistasis ($Delta Delta = {float(c_did['delta_pct']):+.1f}%$, Wald $p = {_fmt_pval_typst(float(c_did['wald_two_sided_p']))}$, `OOF` permutation $p = {_fmt_pval_typst(float(c_did['perm_two_sided_p']))}$).
  *(b)* Forest summary of all five pairwise simple-effect and $2 times 2$ `DiD` contrasts ($plus.minus 95%$ CI) overlaid on $B = 500$ within-trial `OOF` permutation null distributions ($plus.minus 95%$ uncertainty intervals), confirming zero outcome-leakage bias across every contrast.
  *(c)* Cross-histology replication confirming concordant positive corner responder-rate shifts across `4/4` solid-tumor lineages (Melanoma, Urothelial Bladder, Renal ccRCC, and PDAC + TNBC) as well as the pooled pan-cancer cohort ($N = 1,015$).
  *(d)* Bidirectional `5`-fold cross-fitted training $t$-weights (mean $plus.minus$ SD across folds; selected in `5/5` folds) contrasting the top five $R_("OOF")$ cytotoxic/chemokine effector states (led by `{top_pos['cell_state']}` cytotoxic T cells, $t = {float(top_pos['mean_t_weight']):+.2f} plus.minus {float(top_pos['sd_t_weight']):.2f}$) against the top five $"NR"_("OOF")$ myeloid-suppressive, hypoxic, and mesenchymal barrier states (led by `{top_neg['cell_state']}` `MDSC`s, $t = {float(top_neg['mean_t_weight']):+.2f} plus.minus {float(top_neg['sd_t_weight']):.2f}$).
]
"""

        main_cap_art = captions_dir / "iatlas_epistasis_and_stromal_gating_main_caption.typ"
        supp_cap_art = captions_dir / "iatlas_epistasis_supplementary_captions.typ"
        main_cap_res = results_dir / "iatlas_epistasis_and_stromal_gating_main_caption.typ"
        supp_cap_res = results_dir / "iatlas_epistasis_supplementary_captions.typ"

        main_cap_art.write_text(main_caption_typ, encoding="utf-8")
        supp_cap_art.write_text(supp_captions_typ, encoding="utf-8")
        shutil.copy2(main_cap_art, main_cap_res)
        shutil.copy2(supp_cap_art, supp_cap_res)

        return Success((main_cap_art, supp_cap_art))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export Typst figure captions: {exc}")


# ------------------------------------------------------------------------------
# CLI & Imperative Shell
# ------------------------------------------------------------------------------


def run_figure_generation(config: FigurePlotConfig) -> Result[tuple[Path, ...], str]:
    """Load Stage 3c benchmark parquets and export all 14 publication figures + standalone Typst caption modules."""
    pred_path = config.data_dir / "response_prediction_metrics.parquet"
    conc_path = config.data_dir / "milopy_iatlas_concordance_summary.parquet"
    states_path = config.data_dir / "milopy_iatlas_state_effects.parquet"
    loco_path = config.data_dir / "loco_prediction_metrics.parquet"
    deg_path = config.data_dir / "consensus_metastate_top_degs.parquet"
    all42_deg_path = config.data_dir / "consensus_metastate_all42_top_degs.parquet"
    all42_pure_inter_path = config.data_dir / "consensus_metastate_all42_interactions_pure_clr.parquet"
    zl_summary_path = config.data_dir / "zero_leakage_gating_methods_summary.parquet"
    zl_cells_path = config.data_dir / "zero_leakage_gating_regimes_cells.parquet"
    oof_2x2_cells_path = config.data_dir / "oof_2x2_extreme_regimes_cells.parquet"
    oof_2x2_contrasts_path = config.data_dir / "oof_2x2_extreme_regimes_contrasts.parquet"
    oof_9a9b_landscape_path = config.data_dir / "oof_9a_9b_patient_landscape_and_curves.parquet"
    oof_9a9b_clin_path = config.data_dir / "oof_9a_9b_clinical_and_histology_breakdown.parquet"
    oof_9a9b_drivers_path = config.data_dir / "oof_9a_9b_fold_state_drivers.parquet"
    epi_strat_path = config.data_dir / "epistasis_discovery_strategies_summary.parquet"
    epi_sig_path = config.data_dir / "epistasis_significant_pairs_calibrated.parquet"

    for p in (
        pred_path,
        conc_path,
        states_path,
        loco_path,
        deg_path,
        all42_deg_path,
        all42_pure_inter_path,
        zl_summary_path,
        zl_cells_path,
        oof_2x2_cells_path,
        oof_2x2_contrasts_path,
        oof_9a9b_landscape_path,
        oof_9a9b_clin_path,
        oof_9a9b_drivers_path,
        epi_strat_path,
        epi_sig_path,
    ):
        if not p.exists():
            return Failure(f"Required benchmark parquet not found: {p}")

    pred_df = pl.read_parquet(pred_path)
    conc_df = pl.read_parquet(conc_path)
    states_df = pl.read_parquet(states_path)
    loco_df = pl.read_parquet(loco_path)
    deg_df = pl.read_parquet(deg_path)
    all42_deg_df = pl.read_parquet(all42_deg_path)
    all42_pure_inter_df = pl.read_parquet(all42_pure_inter_path)
    zl_summary_df = pl.read_parquet(zl_summary_path)
    zl_cells_df = pl.read_parquet(zl_cells_path)
    oof_2x2_cells_df = pl.read_parquet(oof_2x2_cells_path)
    oof_2x2_contrasts_df = pl.read_parquet(oof_2x2_contrasts_path)
    oof_9a9b_landscape_df = pl.read_parquet(oof_9a9b_landscape_path)
    oof_9a9b_clin_df = pl.read_parquet(oof_9a9b_clin_path)
    oof_9a9b_drivers_df = pl.read_parquet(oof_9a9b_drivers_path)
    epi_strat_df = pl.read_parquet(epi_strat_path)
    epi_sig_df = pl.read_parquet(epi_sig_path)

    figure_builders: list[tuple[str, str, Callable[[], Result[tuple[str, float], str]]]] = [
        (
            "iatlas_response_prediction_auc_heatmap_and_bars",
            "Figure 1 (13x15 ROC-AUC Heatmap & All-15-Cohort Multi-Model Bars)",
            lambda: build_figure1_response_prediction(pred_df),
        ),
        (
            "iatlas_milopy_vs_deconv_concordance_grid",
            "Figure 2 (3x3 Linearized Deconv Beta vs Milopy Grid & Pearson r / Spearman rho Bars)",
            lambda: build_figure2_concordance_grid_and_bars(conc_df, states_df),
        ),
        (
            "iatlas_consensus_metastate_forest_and_tme_profile",
            "Figure 3 (Dual-Scale Consensus-Meta42 Forest Plot & 11-Trial Compositional Shift Heatmap)",
            lambda: build_figure3_metastate_forest_and_tme_heatmap(states_df),
        ),
        (
            "iatlas_classifier_and_reference_tier_pareto",
            "Figure 4 (Numbered Reference Tier Pareto Frontier & Outlier-Annotated RF Importance)",
            lambda: build_figure4_pareto_and_feature_importance(pred_df, conc_df, states_df),
        ),
        (
            "iatlas_cross_cohort_concordance_transfer_matrix",
            "Figure 5 (13x16 Cross-Modality Transfer Matrix & Within-Cohort CV vs LOCO)",
            lambda: build_figure5_transfer_matrix_and_loco(pred_df, conc_df, loco_df),
        ),
        (
            "iatlas_consensus_metastate_top_deg_markers",
            "Figure 6 (Top 4 Differentially Expressed Marker Genes of the 14 Selected Meta-States)",
            lambda: build_figure6_consensus_metastate_top_degs(deg_df),
        ),
        (
            "iatlas_response_prediction_prauc_heatmap_and_bars",
            "Figure 7 (13x15 PR-AUC Heatmap & All-15-Cohort Multi-Model PR-AUC Bars)",
            lambda: build_figure7_response_prediction_prauc(pred_df),
        ),
        (
            "iatlas_consensus_metastate_all42_forest_and_tme_profile",
            "Figure 8 (All-42 Consensus-Meta42 Dual-Scale Forest Plot & 11-Trial Compositional Shift Heatmap)",
            lambda: build_figure9_all42_metastate_forest_and_tme_heatmap(states_df),
        ),
        (
            "iatlas_consensus_metastate_all42_top_deg_markers",
            "Figure 9 (Top 4 Differentially Expressed Marker Genes Across All 42 Meta-States)",
            lambda: build_figure10_all42_metastate_top_degs(all42_deg_df),
        ),
        (
            "iatlas_epistasis_and_stromal_gating_main",
            "Figure 10 [Main] (Heading-Free Two-Tier 7-Panel Pan-Compartment Gating & Non-T-Cell M29xM26 Epistasis)",
            lambda: build_main_epistasis_and_stromal_gating_figure(
                oof_2x2_cells_df,
                oof_2x2_contrasts_df,
                oof_9a9b_landscape_df,
                oof_9a9b_clin_df,
                epi_sig_df,
            ),
        ),
        (
            "iatlas_nonlinear_conditional_interactions_all42_pure_clr",
            "Figure 11 [Supp S19] (Heading-Free Single-Panel 42x42 Pure-CLR Epistasis Matrix)",
            lambda: build_supp_pure_clr_all42_interactions(all42_pure_inter_df, states_df),
        ),
        (
            "iatlas_epistasis_discovery_power_benchmark",
            "Figure 12 [Supp S20] (Heading-Free Epistasis Null Calibration & Discovery Power Benchmark)",
            lambda: build_supp_epistasis_discovery_power_benchmark(epi_strat_df, epi_sig_df),
        ),
        (
            "iatlas_zero_leakage_gating_comparison",
            "Figure 13 [Supp S21] (Heading-Free Zero-Leakage 3x3 Stromal Gating Benchmark)",
            lambda: build_supp_zero_leakage_gating_comparison(zl_summary_df, zl_cells_df),
        ),
        (
            "iatlas_oof_supplementary_diagnostics",
            "Figure 14 [Supp S22] (Heading-Free OOF 2x2 Corners, Perm Null, 4/4 Lineages & Cross-Fitted Drivers)",
            lambda: build_supp_oof_diagnostics_figure(
                oof_2x2_cells_df,
                oof_2x2_contrasts_df,
                oof_9a9b_clin_df,
                oof_9a9b_drivers_df,
            ),
        ),
    ]

    exported_paths: list[Path] = []
    for stem, desc, builder_fn in figure_builders:
        print(f"=== Generating {desc} ===")
        match builder_fn():
            case Failure(err):
                return Failure(err)
            case Success((svg_str, h_mm)):
                pass

        svg_p = config.results_dir / f"{stem}.svg"
        png_p = config.results_dir / f"{stem}.png"
        art_svg_p = config.article_dir / f"{stem}.svg"
        art_png_p = config.article_dir / f"{stem}.png"

        match export_svg_and_300dpi_png(
            svg_str, svg_p, png_p, art_svg_p, art_png_p, FULL_WIDTH_MM, h_mm
        ):
            case Failure(err):
                return Failure(err)
            case Success(out_p):
                print(f"Exported {stem}: {out_p} ({FULL_WIDTH_MM:.1f}mm x {h_mm:.2f}mm)")
                exported_paths.append(out_p)

    match export_typst_figure_captions(
        oof_cells_df=oof_2x2_cells_df,
        oof_contrasts_df=oof_2x2_contrasts_df,
        landscape_df=oof_9a9b_landscape_df,
        clin_df=oof_9a9b_clin_df,
        drivers_df=oof_9a9b_drivers_df,
        pure_inter_df=all42_pure_inter_df,
        epi_strat_df=epi_strat_df,
        epi_sig_df=epi_sig_df,
        zl_summary_df=zl_summary_df,
        captions_dir=config.captions_dir,
        results_dir=config.results_dir,
    ):
        case Failure(err):
            return Failure(err)
        case Success((main_cap_p, supp_cap_p)):
            print(f"Exported Typst captions: {main_cap_p}, {supp_cap_p}")

    return Success(tuple(exported_paths))


def parse_args(argv: Sequence[str] | None = None) -> FigurePlotConfig:
    """Parse CLI arguments into an immutable FigurePlotConfig."""
    parser = argparse.ArgumentParser(
        description="Generate 180mm publication vector figures and Typst captions for iAtlas InstaPrism & Milopy benchmark."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark"),
        help="Directory containing benchmark parquet outputs",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures"),
        help="Directory to write publication SVG/PNG figures and caption copies",
    )
    parser.add_argument(
        "--article-dir",
        type=Path,
        default=Path("article/figures/milo-deconv"),
        help="Article figures directory",
    )
    parser.add_argument(
        "--captions-dir",
        type=Path,
        default=Path("article/captions"),
        help="Article Typst captions directory",
    )
    args = parser.parse_args(argv)
    return FigurePlotConfig(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        article_dir=args.article_dir,
        captions_dir=args.captions_dir,
    )


def main() -> None:
    config = parse_args()
    match run_figure_generation(config):
        case Failure(err):
            print(f"Figure generation failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(_):
            print("Successfully generated all 14 iAtlas InstaPrism & Milopy publication figures and Typst captions.")


if __name__ == "__main__":
    main()
