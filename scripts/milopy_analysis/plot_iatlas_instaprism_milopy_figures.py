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

import altair as alt
import numpy as np
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

    data_dir: Path = Path("output/results/iatlas_instaprism_benchmark")
    results_dir: Path = Path("results/iatlas_instaprism_benchmark")
    article_dir: Path = Path("article/figures/milo-deconv")


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


def wrap_svg_with_inkscape_layers(
    svg_text: str,
    width_mm: float,
    layer_id: str,
    layer_label: str,
) -> tuple[str, float]:
    """Inject Inkscape namespaces, explicit layer group, and exact physical dimensions in mm."""
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    vb_w = float(m.group(1)) if m else 680.0
    vb_h = float(m.group(2)) if m else 440.0
    height_mm = width_mm * (vb_h / max(vb_w, 1.0))

    inner = re.sub(r"^<\?xml[^>]*\?>\s*", "", svg_text.strip())
    inner = re.sub(r"^<svg[^>]*>", "", inner.strip())
    inner = re.sub(r"</svg>\s*$", "", inner.strip())
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
    return wrapped, height_mm


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
    "M01": "M01: TRBV18+GCNT4+ Mem T",
    "M02": "M02: GATA4+CADM2+ Stroma",
    "M03": "M03: CD200+ELMO1+ Tfh T",
    "M04": "M04: GNLY+PRF1+ Eff. T",
    "M05": "M05: CD96+TRAF3+ Mem. T",
    "M06": "M06: CX3CR1+FGFBP2+ NK/T",
    "M07": "M07: ANK1+TRIM58+ Eryth.",
    "M08": "M08: CCR4+SOCS1+ Th2/Treg",
    "M09": "M09: SIGLEC8+CR1+ Granul.",
    "M10": "M10: FCN1+CD1A+ Mono/cDC",
    "M11": "M11: AQP4+CTSE+ Ductal",
    "M12": "M12: CA9+SLC2A1+ Hypox.",
    "M13": "M13: CXCL13+PDCD1+ Tex",
    "M14": "M14: CXCL9+STAT1+ M1-TAM",
    "M15": "M15: TCF7+IL7R+ Stem CD8",
    "M16": "M16: VEGFA+KDR+ Endo.",
    "M17": "M17: ISLR+COL5A2+ myCAF",
    "M18": "M18: S100B+TYR+ Melanoc.",
    "M19": "M19: S100A8+CXCL8+ MDSC",
    "M20": "M20: IL7R+CD40LG+ Trm T",
    "M21": "M21: FOXP3+CTLA4+ Treg",
    "M22": "M22: MKI67+TOP2A+ Prolif",
    "M23": "M23: FAP+TGFB1+ Mat. CAF",
    "M24": "M24: CXCL13+NMB+ Tfh/Trm",
    "M25": "M25: CTLA4+TIGIT+ Tex",
    "M26": "M26: VIM+ZEB1+ EMT Tumor",
    "M27": "M27: ABCC8+PPY+ Neuroend",
    "M28": "M28: CXCL8+ITGAM+ MDSC",
    "M29": "M29: JCHAIN+MZB1+ Plasma",
    "M30": "M30: CR2+PAX5+ TLS B",
    "M31": "M31: CLEC4C+LILRA4+ pDC",
    "M32": "M32: MRC1+CD163+ M2-TAM",
    "M33": "M33: CCL5+GNLY+ CD8 Teff",
    "M34": "M34: GZMK+NKG7+ Tem CD8+",
    "M35": "M35: MPO+CD79A+ Myeloid",
    "M36": "M36: FCGR3A+GNLY+ NK",
    "M37": "M37: MS4A1+CD79A+ TLS B",
    "M38": "M38: FCGR3B+FPR2+ TAN",
    "M39": "M39: MS4A2+TPSAB1+ Mast",
    "M40": "M40: FOXP3+CCR8+ Treg",
    "M41": "M41: LAMP3+CCR7+ mregDC",
    "M42": "M42: SULT1B1+RGS2+ TAM",
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
# Figure 8 (NEW): Non-Linear Conditional Interactions (Pairwise Epistasis/H-Stat & 3x3 Stromal Gating Regimes)
# ------------------------------------------------------------------------------


def build_figure8_nonlinear_conditional_interactions(
    inter_df: pl.DataFrame,
    regime_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 8 (180.0mm): (a) 8x8 Pairwise Interaction Matrix + (b) 3x3 Lymphoid x Stromal Exclusion Regimes."""
    try:
        short_inter_map = {
            "M04": "M04: GNLY+ γδ/Eff. T",
            "M20": "M20: IL7R+ CD40LG+ Trm",
            "M30": "M30: MS4A1+ TLS B",
            "M29": "M29: JCHAIN+ Plasma",
            "M05": "M05: TRAF3+ CD96+ T",
            "M25": "M25: CTLA4+ TIGIT+ Tex",
            "M17": "M17: COL5A2+ myCAF",
            "M32": "M32: MRC1+ M2-TAM",
        }
        state_lbl_order = list(short_inter_map.values())
        inter_pd = (
            inter_df.with_columns([
                pl.col("state_a").replace(short_inter_map).alias("state_a_label"),
                pl.col("state_b").replace(short_inter_map).alias("state_b_label"),
            ])
            .to_pandas()
        )

        panel_a = (
            alt.Chart(inter_pd)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.8)
            .encode(
                x=alt.X(
                    "state_b_label:N",
                    sort=state_lbl_order,
                    title="Interacting Consensus Meta-State B",
                    axis=alt.Axis(
                        labelAngle=-40,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=185,
                    ),
                ),
                y=alt.Y(
                    "state_a_label:N",
                    sort=state_lbl_order,
                    title="Consensus Meta-State A (Top RF Outliers & Anchors)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=210,
                    ),
                ),
                color=alt.Color(
                    "epistasis_z:Q",
                    title="Epistasis / Gating z",
                    scale=alt.Scale(domain=[-4.2, 4.2], scheme="blueorange"),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=125,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
            .properties(
                width=235,
                height=265,
                title=alt.TitleParams(
                    text="(a)  Pairwise Non-Linear Epistasis & RF Interaction",
                    subtitle=[
                        "Signed RF 2D partial-dependence & logistic interaction z-score:",
                        "blue = lymphoid-humoral synergy (+/+); orange = stromal gating (+/−)",
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

        # Panel (b): 3x3 Stratified Regime Plot (Lymphoid/TLS Humoral Score x myCAF/M2-TAM Stromal Exclusion)
        lymph_order = [
            "Low Lymphoid/TLS (T1)",
            "Med. Lymphoid/TLS (T2)",
            "High Lymphoid/TLS (T3)",
        ]
        stroma_order = [
            "Low myCAF/M2 Stroma (Permissive)",
            "Moderate myCAF/M2 Stroma",
            "High myCAF/M2 Stroma (Excluded)",
        ]
        stroma_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_VERMILION]
        reg_pd = regime_df.to_pandas()

        lines_b = (
            alt.Chart(reg_pd)
            .mark_line(point=alt.OverlayMarkDef(size=75, filled=True, stroke="#0F172A", strokeWidth=0.7), strokeWidth=2.2)
            .encode(
                x=alt.X(
                    "lymphoid_tier:N",
                    sort=lymph_order,
                    title="Lymphoid–Humoral Score (M04–M30)",
                    axis=alt.Axis(
                        labelAngle=-25,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=180,
                    ),
                ),
                y=alt.Y(
                    "responder_rate_pct:Q",
                    scale=alt.Scale(domain=[5.0, 74.0]),
                    title="Observed ICB Responder Rate (%)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                color=alt.Color(
                    "stromal_regime:N",
                    scale=alt.Scale(domain=stroma_order, range=stroma_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=240,
                        symbolSize=90,
                    ),
                ),
            )
        )

        ci_b = (
            alt.Chart(reg_pd)
            .mark_rule(strokeWidth=1.4, opacity=0.75)
            .encode(
                x=alt.X("lymphoid_tier:N", sort=lymph_order),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=[5.0, 74.0])),
                y2="ci_high_pct:Q",
                color=alt.Color("stromal_regime:N", scale=alt.Scale(domain=stroma_order, range=stroma_colors)),
            )
        )

        panel_b = (ci_b + lines_b).properties(
            width=225,
            height=265,
            title=alt.TitleParams(
                text="(b)  Conditional Stromal Gating of ICB Response",
                subtitle=[
                    "Lymphoid/TLS score drives steep gains in permissive",
                    "tumors but is blunted by M17/M32 stromal exclusion",
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
            alt.hconcat(panel_a, panel_b, spacing=16)
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
            layer_id="layer-fig8-nonlinear-interactions",
            layer_label="08_NonLinear_Conditional_Interactions",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 8: {exc}")


# ------------------------------------------------------------------------------
# Figure 9 (NEW): Full-Atlas (All 42 States) Dual-Scale Forest Plot & 11-Trial Heatmap
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
# Figure 11 (NEW): Full 42x42 Non-Linear Conditional Interactions, Hubs & Pan-Atlas Regimes
# ------------------------------------------------------------------------------


def build_figure11_all42_nonlinear_conditional_interactions(
    all42_inter_df: pl.DataFrame,
    all42_hub_df: pl.DataFrame,
    all42_regime_df: pl.DataFrame,
) -> Result[tuple[str, float], str]:
    """Build Figure 11 (180.0mm): (a) 42x42 Epistasis Matrix + (b) 42-State Hub Scatter + (c) Pan-Atlas 3x3 Regimes."""
    try:
        hub_sorted = (
            all42_hub_df.with_columns(
                pl.col("cell_state").replace(META_STATE_EMPIRICAL_LABELS).alias("state_label")
            )
            .sort("state_bulk_beta", descending=True)
        )
        ordered_codes = hub_sorted["cell_state"].to_list()
        ordered_labels = hub_sorted["state_label"].to_list()

        inter_clean = all42_inter_df.with_columns([
            pl.col("state_a").replace(META_STATE_EMPIRICAL_LABELS).alias("state_a_label"),
            pl.col("state_b").replace(META_STATE_EMPIRICAL_LABELS).alias("state_b_label"),
        ])
        inter_pd = inter_clean.to_pandas()

        panel_a = (
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
                    title="Signed Epistasis / Gating z",
                    scale=alt.Scale(
                        domain=[-4.2, 4.2],
                        scheme="blueorange",
                        reverse=True,
                    ),
                    legend=alt.Legend(
                        orient="top",
                        direction="horizontal",
                        gradientLength=150,
                        gradientThickness=10,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                    ),
                ),
            )
            .properties(
                width=615,
                height=625,
                title=alt.TitleParams(
                    text="(a)  Complete 42 × 42 Pairwise Non-Linear Epistasis & RF Interaction Matrix (1,764 Pairs)",
                    subtitle=[
                        "Signed RF 2D partial-dependence & logistic epistasis z across all 42 Consensus-Meta42 states:",
                        "blue (z > 0) = co-activation synergy (+/+); orange (z < 0) = conditional stromal-myeloid gating (+/−)",
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

        # Panel (b): Global Epistatic Hub Strength vs Marginal |beta| Across All 42 States
        hub_pd = hub_sorted.to_pandas()
        grp_order = ["Responder-Enriched (β > 0)", "Non-Responder-Enriched (β < 0)"]
        grp_colors = [OKABE_BLUE, OKABE_VERMILION]
        x_hub_scale = alt.Scale(domain=[0.0, 3.60], zero=False)
        y_hub_scale = alt.Scale(domain=[1.55, 2.95], zero=False)

        hub_pts = (
            alt.Chart(hub_pd)
            .mark_circle(opacity=0.88, stroke="#0F172A", strokeWidth=0.75)
            .encode(
                x=alt.X(
                    "abs_bulk_beta:Q",
                    scale=x_hub_scale,
                    title="Marginal Meta-9 |β| (42 Meta-States)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                y=alt.Y(
                    "mean_abs_epistasis_z:Q",
                    scale=y_hub_scale,
                    title="Global Epistatic Hub (Mean |z|)",
                    axis=alt.Axis(
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                    ),
                ),
                size=alt.Size(
                    "mean_friedman_h:Q",
                    scale=alt.Scale(range=[45, 155]),
                    legend=None,
                ),
                color=alt.Color(
                    "enrichment_group:N",
                    scale=alt.Scale(domain=grp_order, range=grp_colors),
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

        hub_trend = (
            alt.Chart(hub_pd)
            .transform_regression("abs_bulk_beta", "mean_abs_epistasis_z")
            .mark_line(color="#1E293B", strokeDash=[4, 2], strokeWidth=1.5)
            .encode(
                x=alt.X("abs_bulk_beta:Q", scale=x_hub_scale),
                y=alt.Y("mean_abs_epistasis_z:Q", scale=y_hub_scale),
            )
        )

        hub_callouts: list[tuple[str, str, Literal["left", "center", "right"], int, int]] = [
            ("M29", "M29 (Plasma)", "center", -12, -12),
            ("M14", "M14 (M1-TAM)", "center", 0, -13),
            ("M13", "M13 (Tex)", "left", 8, -4),
            ("M23", "M23 (CAF)", "right", -9, 5),
            ("M21", "M21 (Treg)", "left", 8, 6),
        ]
        callout_layers: list[alt.Chart] = []
        for st_c, txt_c, al_c, dx_c, dy_c in hub_callouts:
            sub_h = hub_sorted.filter(pl.col("cell_state") == st_c).with_columns(
                pl.lit(txt_c).alias("hub_lbl")
            )
            if sub_h.height > 0:
                callout_layers.append(
                    alt.Chart(sub_h.to_pandas())
                    .mark_text(
                        align=al_c,
                        dx=dx_c,
                        dy=dy_c,
                        fontSize=MIN_FONT_PX,
                        fontWeight="bold",
                        color="#0F172A",
                    )
                    .encode(
                        x=alt.X("abs_bulk_beta:Q", scale=x_hub_scale),
                        y=alt.Y("mean_abs_epistasis_z:Q", scale=y_hub_scale),
                        text="hub_lbl:N",
                    )
                )

        panel_b_base = hub_trend + hub_pts
        for c_layer in callout_layers:
            panel_b_base = panel_b_base + c_layer

        panel_b = panel_b_base.properties(
            width=335,
            height=175,
            title=alt.TitleParams(
                text="(b)  All-42 State Epistatic Hub Strength",
                subtitle=[
                    "Mean pairwise |epistasis z| across all 41 partners",
                    "vs. marginal Meta-9 |β| (circle size = RF Friedman H)",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        # Panel (c): Pan-Atlas 3x3 Conditional Gating Regimes (All R-Enriched vs All NR-Enriched States)
        r_order = [
            "Low Pan-R Score (T1)",
            "Med. Pan-R Score (T2)",
            "High Pan-R Score (T3)",
        ]
        n_order = [
            "Permissive Low-NR TME",
            "Moderate NR TME",
            "Excluded High-NR TME",
        ]
        n_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_VERMILION]
        reg_pd = all42_regime_df.to_pandas()

        lines_c = (
            alt.Chart(reg_pd)
            .mark_line(
                point=alt.OverlayMarkDef(size=70, filled=True, stroke="#0F172A", strokeWidth=0.7),
                strokeWidth=2.1,
            )
            .encode(
                x=alt.X(
                    "lymphoid_tier:N",
                    sort=r_order,
                    title="Pan-Atlas Responder Composite (All β > 0)",
                    axis=alt.Axis(
                        labelAngle=-22,
                        labelFontSize=MIN_FONT_PX,
                        titleFontSize=MIN_FONT_PX,
                        titleFontWeight="bold",
                        labelLimit=180,
                    ),
                ),
                y=alt.Y(
                    "responder_rate_pct:Q",
                    scale=alt.Scale(domain=[5.0, 74.0]),
                    title="ICB Responder Rate (%)",
                    axis=alt.Axis(
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
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=MIN_FONT_PX,
                        labelLimit=230,
                        symbolSize=85,
                    ),
                ),
            )
        )

        ci_c = (
            alt.Chart(reg_pd)
            .mark_rule(strokeWidth=1.4, opacity=0.75)
            .encode(
                x=alt.X("lymphoid_tier:N", sort=r_order),
                y=alt.Y("ci_low_pct:Q", scale=alt.Scale(domain=[5.0, 74.0])),
                y2="ci_high_pct:Q",
                color=alt.Color("stromal_regime:N", scale=alt.Scale(domain=n_order, range=n_colors)),
            )
        )

        panel_c = (ci_c + lines_c).properties(
            width=260,
            height=175,
            title=alt.TitleParams(
                text="(c)  Pan-Atlas 3 × 3 Gating Regimes",
                subtitle=[
                    "All β > 0 effector composite stratified by",
                    "all β < 0 suppressive TME (N = 1,015)",
                ],
                fontSize=TITLE_FONT_PX,
                subtitleFontSize=MIN_FONT_PX,
                subtitleColor="#475569",
                fontWeight="bold",
                anchor="start",
                offset=6,
            ),
        )

        bottom_row = alt.hconcat(panel_b, panel_c, spacing=22).resolve_scale(color="independent")

        combined = (
            alt.vconcat(panel_a, bottom_row, spacing=18)
            .resolve_scale(color="independent", size="independent")
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
            layer_id="layer-fig11-all42-nonlinear-interactions",
            layer_label="11_All42_NonLinear_Conditional_Interactions",
        )
        return Success((wrapped_svg, height_mm))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build Figure 11: {exc}")


# ------------------------------------------------------------------------------
# CLI & Imperative Shell
# ------------------------------------------------------------------------------


def run_figure_generation(config: FigurePlotConfig) -> Result[tuple[Path, ...], str]:
    """Load benchmark parquet outputs and export all 11 publication figures as 180mm SVGs & 300 DPI PNGs."""
    pred_path = config.data_dir / "response_prediction_metrics.parquet"
    conc_path = config.data_dir / "milopy_iatlas_concordance_summary.parquet"
    states_path = config.data_dir / "milopy_iatlas_state_effects.parquet"
    loco_path = config.data_dir / "loco_prediction_metrics.parquet"
    deg_path = config.data_dir / "consensus_metastate_top_degs.parquet"
    inter_path = config.data_dir / "consensus_metastate_interactions.parquet"
    regime_path = config.data_dir / "consensus_metastate_regimes.parquet"
    all42_deg_path = config.data_dir / "consensus_metastate_all42_top_degs.parquet"
    all42_inter_path = config.data_dir / "consensus_metastate_all42_interactions.parquet"
    all42_hub_path = config.data_dir / "consensus_metastate_all42_hub_summary.parquet"
    all42_regime_path = config.data_dir / "consensus_metastate_all42_regimes.parquet"

    for p in (
        pred_path,
        conc_path,
        states_path,
        loco_path,
        deg_path,
        inter_path,
        regime_path,
        all42_deg_path,
        all42_inter_path,
        all42_hub_path,
        all42_regime_path,
    ):
        if not p.exists():
            return Failure(f"Required benchmark parquet not found: {p}")

    pred_df = pl.read_parquet(pred_path)
    conc_df = pl.read_parquet(conc_path)
    states_df = pl.read_parquet(states_path)
    loco_df = pl.read_parquet(loco_path)
    deg_df = pl.read_parquet(deg_path)
    inter_df = pl.read_parquet(inter_path)
    regime_df = pl.read_parquet(regime_path)
    all42_deg_df = pl.read_parquet(all42_deg_path)
    all42_inter_df = pl.read_parquet(all42_inter_path)
    all42_hub_df = pl.read_parquet(all42_hub_path)
    all42_regime_df = pl.read_parquet(all42_regime_path)

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
            "iatlas_nonlinear_conditional_interactions",
            "Figure 8 (Pairwise Non-Linear Conditional Interactions & 3x3 Stromal Gating Regimes)",
            lambda: build_figure8_nonlinear_conditional_interactions(inter_df, regime_df),
        ),
        (
            "iatlas_consensus_metastate_all42_forest_and_tme_profile",
            "Figure 9 (All-42 Consensus-Meta42 Dual-Scale Forest Plot & 11-Trial Compositional Shift Heatmap)",
            lambda: build_figure9_all42_metastate_forest_and_tme_heatmap(states_df),
        ),
        (
            "iatlas_consensus_metastate_all42_top_deg_markers",
            "Figure 10 (Top 4 Differentially Expressed Marker Genes Across All 42 Meta-States)",
            lambda: build_figure10_all42_metastate_top_degs(all42_deg_df),
        ),
        (
            "iatlas_nonlinear_conditional_interactions_all42",
            "Figure 11 (Complete 42x42 Pairwise Non-Linear Conditional Interactions, Hubs & Pan-Atlas Regimes)",
            lambda: build_figure11_all42_nonlinear_conditional_interactions(
                all42_inter_df, all42_hub_df, all42_regime_df
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

    return Success(tuple(exported_paths))


def parse_args(argv: Sequence[str] | None = None) -> FigurePlotConfig:
    """Parse CLI arguments into an immutable FigurePlotConfig."""
    parser = argparse.ArgumentParser(
        description="Generate 180mm publication vector figures for iAtlas InstaPrism & Milopy benchmark."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("output/results/iatlas_instaprism_benchmark"),
        help="Directory containing benchmark parquet outputs",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/iatlas_instaprism_benchmark"),
        help="Directory to write publication SVG and PNG figures",
    )
    parser.add_argument(
        "--article-dir",
        type=Path,
        default=Path("article/figures/milo-deconv"),
        help="Article figures directory",
    )
    args = parser.parse_args(argv)
    return FigurePlotConfig(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        article_dir=args.article_dir,
    )


def main() -> None:
    config = parse_args()
    match run_figure_generation(config):
        case Failure(err):
            print(f"Figure generation failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(_):
            print("Successfully generated all 11 iAtlas InstaPrism & Milopy publication figures.")


if __name__ == "__main__":
    main()
