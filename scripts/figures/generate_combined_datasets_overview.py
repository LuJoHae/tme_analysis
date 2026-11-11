#!/usr/bin/env python3
"""
Generate a 3-panel publication-grade vector SVG compendium (180mm x 230mm)
with clean Nature Methods typography (6.0–6.4 pt body/tables, 7.1 pt titles, 8.25 pt badges):
  - Panel a: 9 Bulk RNA-seq ICB Clinical Trial Cohorts (n = 1,097 patients)
  - Panel b: 9 Single-Cell ICB Response Cohorts used in Milo (milopy) Analysis (n = 613 evaluations, 3,597,092 cells)
  - Panel c: 19 Single-Cell Reference Datasets used in the Multi-Resolution Deconvolution & Integration Benchmark (1,001,078 cells)

Coordinate system:
  - Physical dimensions: 180mm x 230mm
  - Native 96 DPI user units: viewBox="0 0 680 870" (1 pt = 1.3333 px, 6.38pt = 8.5 px)
"""

import html
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Final

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# ==============================================================================
# CONSTANTS & NATURE METHODS PALETTE (OKABE-ITO COLORBLIND-SAFE)
# ==============================================================================
SVG_WIDTH_MM: Final[str] = "180mm"
SVG_HEIGHT_MM: Final[str] = "230mm"
VIEWBOX_W: Final[int] = 680
VIEWBOX_H: Final[int] = 870

COLOR_CANVAS_BG: Final[str] = "#FFFFFF"
COLOR_PANEL_BG: Final[str] = "#FFFFFF"
COLOR_SUBTLE_FILL: Final[str] = "#F8FAFC"
COLOR_BORDER_HAIRLINE: Final[str] = "#CBD5E1"
COLOR_DIVIDER_RULE: Final[str] = "#E2E8F0"

COLOR_TEXT_PRIMARY: Final[str] = "#0F172A"
COLOR_TEXT_SECONDARY: Final[str] = "#334155"
COLOR_TEXT_MUTED: Final[str] = "#475569"

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_SKY_BLUE: Final[str] = "#56B4E9"
OKABE_BLUISH_GREEN: Final[str] = "#009E73"
OKABE_ORANGE: Final[str] = "#E69F00"
OKABE_VERMILION: Final[str] = "#D55E00"
OKABE_REDDISH_PURPLE: Final[str] = "#CC79A7"

FONT_SANS: Final[str] = "system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"


# ==============================================================================
# IMMUTABLE DATA MODELS (FROZEN PYDANTIC, NO NONE)
# ==============================================================================
class BulkTrialCohort(BaseModel):
    """Immutable model for a bulk RNA-seq ICB clinical trial cohort (Panel a)."""

    model_config = ConfigDict(frozen=True)

    study: str
    cancer_type: str
    cancer_color: str
    therapy: str
    n_total: int
    n_resp: int
    n_non_resp: int
    orr_pct: float


class SingleCellMiloCohort(BaseModel):
    """Immutable model for a single-cell ICB cohort analyzed with Milo (Panel b)."""

    model_config = ConfigDict(frozen=True)

    publication: str
    cancer_type: str
    cancer_color: str
    therapy: str
    cells_raw: int
    n_resp: int
    n_non_resp: int
    orr_pct: float


class DeconvReferenceDataset(BaseModel):
    """Immutable model for a single-cell reference dataset in the benchmark plot (Panel c)."""

    model_config = ConfigDict(frozen=True)

    index_str: str
    study: str
    cancer_type: str
    cancer_color: str
    cells_full: int


# ==============================================================================
# DATASET COMPENDIUM DATA DEFINITIONS
# ==============================================================================
BULK_TRIAL_COHORTS: Final[tuple[BulkTrialCohort, ...]] = (
    BulkTrialCohort(
        study="Hugo (2016)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        therapy="Pembrolizumab (PD-1)",
        n_total=27,
        n_resp=14,
        n_non_resp=13,
        orr_pct=51.9,
    ),
    BulkTrialCohort(
        study="Riaz (2017)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        therapy="Nivolumab (PD-1)",
        n_total=107,
        n_resp=34,
        n_non_resp=73,
        orr_pct=31.8,
    ),
    BulkTrialCohort(
        study="Liu (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        therapy="Nivolumab / Pembro (PD-1)",
        n_total=122,
        n_resp=48,
        n_non_resp=74,
        orr_pct=39.3,
    ),
    BulkTrialCohort(
        study="Gide (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        therapy="Anti-PD-1 ± Ipi (Combo)",
        n_total=91,
        n_resp=63,
        n_non_resp=28,
        orr_pct=69.2,
    ),
    BulkTrialCohort(
        study="Rosenberg (2016)",
        cancer_type="Bladder",
        cancer_color=OKABE_SKY_BLUE,
        therapy="Atezolizumab (PD-L1)",
        n_total=347,
        n_resp=68,
        n_non_resp=279,
        orr_pct=19.6,
    ),
    BulkTrialCohort(
        study="McDermott (2018)",
        cancer_type="Renal (ccRCC)",
        cancer_color=OKABE_BLUISH_GREEN,
        therapy="Atezo ± Bevacizumab",
        n_total=263,
        n_resp=97,
        n_non_resp=166,
        orr_pct=36.9,
    ),
    BulkTrialCohort(
        study="Choueiri (2016)",
        cancer_type="Renal (ccRCC)",
        cancer_color=OKABE_BLUISH_GREEN,
        therapy="Nivolumab (PD-1)",
        n_total=16,
        n_resp=6,
        n_non_resp=10,
        orr_pct=37.5,
    ),
    BulkTrialCohort(
        study="Padron (2022)",
        cancer_type="Pancreatic",
        cancer_color=OKABE_BLUE,
        therapy="Nivo + Sotigalimab",
        n_total=93,
        n_resp=16,
        n_non_resp=77,
        orr_pct=17.2,
    ),
    BulkTrialCohort(
        study="Anders (2021)",
        cancer_type="Breast (TNBC)",
        cancer_color=OKABE_REDDISH_PURPLE,
        therapy="Atezo + Nab-paclitaxel",
        n_total=31,
        n_resp=11,
        n_non_resp=20,
        orr_pct=35.5,
    ),
)

SC_MILO_COHORTS: Final[tuple[SingleCellMiloCohort, ...]] = (
    SingleCellMiloCohort(
        publication="Sade-Feldman (2018)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        therapy="Anti-PD-1 ± CTLA-4",
        cells_raw=16291,
        n_resp=14,
        n_non_resp=25,
        orr_pct=35.9,
    ),
    SingleCellMiloCohort(
        publication="Gondal (2025)",
        cancer_type="Pan-Cancer",
        cancer_color=OKABE_ORANGE,
        therapy="Nivolumab / Pembro",
        cells_raw=332972,
        n_resp=41,
        n_non_resp=82,
        orr_pct=33.3,
    ),
    SingleCellMiloCohort(
        publication="Marteau (2026)",
        cancer_type="Colorectal",
        cancer_color=OKABE_BLUISH_GREEN,
        therapy="Pembro / Nivolumab",
        cells_raw=1107923,
        n_resp=48,
        n_non_resp=29,
        orr_pct=62.3,
    ),
    SingleCellMiloCohort(
        publication="Yan (2026)",
        cancer_type="Breast",
        cancer_color=OKABE_REDDISH_PURPLE,
        therapy="Neoadjuvant AC Chemo",
        cells_raw=424114,
        n_resp=45,
        n_non_resp=37,
        orr_pct=54.9,
    ),
    SingleCellMiloCohort(
        publication="Hu (2023)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        therapy="Anti-PD-1 + Chemo",
        cells_raw=91997,
        n_resp=4,
        n_non_resp=10,
        orr_pct=28.6,
    ),
    SingleCellMiloCohort(
        publication="Liu (2025)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        therapy="PD-(L)1 + Chemo",
        cells_raw=1254749,
        n_resp=130,
        n_non_resp=112,
        orr_pct=53.7,
    ),
    SingleCellMiloCohort(
        publication="Hong (2025)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        therapy="Atezo + Bev + Chemo",
        cells_raw=40142,
        n_resp=3,
        n_non_resp=4,
        orr_pct=42.9,
    ),
    SingleCellMiloCohort(
        publication="Luoma (2022)",
        cancer_type="HNSCC",
        cancer_color=OKABE_VERMILION,
        therapy="Nivolumab ± Ipi",
        cells_raw=245253,
        n_resp=4,
        n_non_resp=8,
        orr_pct=33.3,
    ),
    SingleCellMiloCohort(
        publication="Raufi (2026)",
        cancer_type="Pancreatic",
        cancer_color=OKABE_BLUE,
        therapy="Cemiplimab + CXCR4i",
        cells_raw=83651,
        n_resp=15,
        n_non_resp=2,
        orr_pct=88.2,
    ),
)

DECONV_REF_DATASETS: Final[tuple[DeconvReferenceDataset, ...]] = (
    DeconvReferenceDataset(
        index_str="01",
        study="Sade-Feldman (2018)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        cells_full=16288,
    ),
    DeconvReferenceDataset(
        index_str="02",
        study="Jerby-Arnon (2018)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        cells_full=7186,
    ),
    DeconvReferenceDataset(
        index_str="03",
        study="Ma (2019)",
        cancer_type="Liver (HCC)",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=5115,
    ),
    DeconvReferenceDataset(
        index_str="04",
        study="Yost (2019)",
        cancer_type="Skin (BCC)",
        cancer_color=OKABE_ORANGE,
        cells_full=15000,
    ),
    DeconvReferenceDataset(
        index_str="05",
        study="Maynard (2020)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        cells_full=27489,
    ),
    DeconvReferenceDataset(
        index_str="06",
        study="Pelka (2021)",
        cancer_type="Colorectal",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=65000,
    ),
    DeconvReferenceDataset(
        index_str="07",
        study="Becker (2022)",
        cancer_type="Colorectal",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=30000,
    ),
    DeconvReferenceDataset(
        index_str="08",
        study="Khaliq (2022)",
        cancer_type="Colorectal",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=25000,
    ),
    DeconvReferenceDataset(
        index_str="09",
        study="Azizi (2018)",
        cancer_type="Breast",
        cancer_color=OKABE_REDDISH_PURPLE,
        cells_full=45000,
    ),
    DeconvReferenceDataset(
        index_str="10",
        study="Leader (2021)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        cells_full=35000,
    ),
    DeconvReferenceDataset(
        index_str="11",
        study="Kim (2020)",
        cancer_type="LUAD",
        cancer_color=OKABE_SKY_BLUE,
        cells_full=40000,
    ),
    DeconvReferenceDataset(
        index_str="12",
        study="Borcherding (2021)",
        cancer_type="ccRCC",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=25000,
    ),
    DeconvReferenceDataset(
        index_str="13",
        study="Pu (2021)",
        cancer_type="Thyroid",
        cancer_color=OKABE_VERMILION,
        cells_full=20000,
    ),
    DeconvReferenceDataset(
        index_str="14",
        study="Durante (2020)",
        cancer_type="Uveal Mel.",
        cancer_color=OKABE_ORANGE,
        cells_full=10000,
    ),
    DeconvReferenceDataset(
        index_str="15",
        study="Biermann (2022)",
        cancer_type="Mel. Brain",
        cancer_color=OKABE_ORANGE,
        cells_full=12000,
    ),
    DeconvReferenceDataset(
        index_str="16",
        study="Cheng (2021)",
        cancer_type="Pan-Cancer",
        cancer_color=OKABE_BLUE,
        cells_full=390000,
    ),
    DeconvReferenceDataset(
        index_str="17",
        study="Sharma (2020)",
        cancer_type="Liver (HCC)",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=15000,
    ),
    DeconvReferenceDataset(
        index_str="18",
        study="Lu (2022)",
        cancer_type="Liver (HCC)",
        cancer_color=OKABE_BLUISH_GREEN,
        cells_full=18000,
    ),
    DeconvReferenceDataset(
        index_str="19",
        study="Qian (2020)",
        cancer_type="Pan-Cancer",
        cancer_color=OKABE_BLUE,
        cells_full=200000,
    ),
)


# ==============================================================================
# PURE FUNCTIONAL SVG RENDERERS (FONTS: 8.0px–8.5px / 6.0–6.4pt TABLES)
# ==============================================================================
def render_bulk_row(idx: int, item: BulkTrialCohort, y_offset: int) -> str:
    """Render a single bulk RNA-seq trial row in Panel a (height = 22px, font-size = 8.5px = 6.38pt)."""
    row_h = 22
    bar_w = 114
    raw_w_r = int(round((item.n_resp / item.n_total) * bar_w))
    w_r = max(36, min(bar_w - 44, raw_w_r))
    w_nr = bar_w - w_r
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL
    r_center = w_r // 2
    nr_center = w_r + (w_nr // 2)

    return f"""
      <!-- Bulk Row {idx + 1}: {html.escape(item.study)} -->
      <g id="bulk-row-{idx + 1:02d}" transform="translate(0, {y_offset})">
        <rect x="0" y="0" width="668" height="{row_h}" fill="{row_bg}" />
        <line x1="0" y1="{row_h}" x2="668" y2="{row_h}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        <text x="5" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.study)}</text>
        <circle cx="150" cy="11.0" r="2.8" fill="{item.cancer_color}" />
        <text x="158" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.cancer_type)}</text>
        <text x="256" y="14.0" font-size="8.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">{html.escape(item.therapy)}</text>
        <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{item.n_total}</text>
        <g transform="translate(495, 3)">
          <rect x="0" y="0" width="{w_r}" height="16" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{w_r}" y="0" width="{w_nr}" height="16" fill="{OKABE_VERMILION}" />
          <text x="{r_center}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{item.n_resp}R</text>
          <text x="{nr_center}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{item.n_non_resp}NR</text>
        </g>
        <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{item.orr_pct:.1f}%</text>
      </g>"""


def render_sc_milo_row(idx: int, item: SingleCellMiloCohort, y_offset: int) -> str:
    """Render a single-cell ICB cohort row in Panel b (height = 22px, font-size = 8.5px = 6.38pt)."""
    row_h = 22
    bar_w = 114
    n_eval = item.n_resp + item.n_non_resp
    raw_w_r = int(round((item.n_resp / n_eval) * bar_w))
    w_r = max(36, min(bar_w - 42, raw_w_r))
    w_nr = bar_w - w_r
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL
    r_center = w_r // 2
    nr_center = w_r + (w_nr // 2)

    return f"""
      <!-- SC Milo Row {idx + 1}: {html.escape(item.publication)} -->
      <g id="sc-milo-row-{idx + 1:02d}" transform="translate(0, {y_offset})">
        <rect x="0" y="0" width="668" height="{row_h}" fill="{row_bg}" />
        <line x1="0" y1="{row_h}" x2="668" y2="{row_h}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
        <text x="5" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.publication)}</text>
        <circle cx="150" cy="11.0" r="2.8" fill="{item.cancer_color}" />
        <text x="158" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.cancer_type)}</text>
        <text x="256" y="14.0" font-size="8.5" font-weight="400" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">{html.escape(item.therapy)}</text>
        <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">{item.cells_raw:,}</text>
        <g transform="translate(495, 3)">
          <rect x="0" y="0" width="{w_r}" height="16" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{w_r}" y="0" width="{w_nr}" height="16" fill="{OKABE_VERMILION}" />
          <text x="{r_center}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{item.n_resp}R</text>
          <text x="{nr_center}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{item.n_non_resp}NR</text>
        </g>
        <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{item.orr_pct:.1f}%</text>
      </g>"""


def render_deconv_ref_row(idx: int, item: DeconvReferenceDataset, y_offset: int) -> str:
    """Render a single reference dataset row in Panel c sub-table (width = 331px, height = 22px, font-size = 8.5px = 6.38pt)."""
    row_h = 22
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL

    return f"""
        <!-- Deconv Ref Row {item.index_str}: {html.escape(item.study)} -->
        <g id="deconv-ref-row-{item.index_str}" transform="translate(0, {y_offset})">
          <rect x="0" y="0" width="331" height="{row_h}" fill="{row_bg}" />
          <line x1="0" y1="{row_h}" x2="331" y2="{row_h}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.5" />
          <text x="4" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_MUTED}" letter-spacing="-0.2">{html.escape(item.index_str)}</text>
          <text x="22" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.study)}</text>
          <circle cx="165" cy="11.0" r="2.8" fill="{item.cancer_color}" />
          <text x="172" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">{html.escape(item.cancer_type)}</text>
          <text x="326" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">{item.cells_full:,}</text>
        </g>"""


def render_svg() -> str:
    """Declaratively render the complete 180mm x 230mm 3-panel dataset compendium SVG with >=10pt fonts."""
    # Aggregates for Panel a
    bulk_total_n = sum(c.n_total for c in BULK_TRIAL_COHORTS)
    bulk_total_r = sum(c.n_resp for c in BULK_TRIAL_COHORTS)
    bulk_total_nr = sum(c.n_non_resp for c in BULK_TRIAL_COHORTS)
    bulk_orr = (bulk_total_r / bulk_total_n) * 100.0
    bulk_bar_w = 114
    bulk_w_r = max(44, int(round((bulk_total_r / bulk_total_n) * bulk_bar_w)))
    bulk_w_nr = bulk_bar_w - bulk_w_r

    # Aggregates for Panel b
    sc_total_cells = sum(c.cells_raw for c in SC_MILO_COHORTS)
    sc_total_r = sum(c.n_resp for c in SC_MILO_COHORTS)
    sc_total_nr = sum(c.n_non_resp for c in SC_MILO_COHORTS)
    sc_total_eval = sc_total_r + sc_total_nr
    sc_orr = (sc_total_r / sc_total_eval) * 100.0
    sc_bar_w = 114
    sc_w_r = int(round((sc_total_r / sc_total_eval) * sc_bar_w))
    sc_w_nr = sc_bar_w - sc_w_r

    # Aggregates for Panel c
    ref_total_cells = sum(c.cells_full for c in DECONV_REF_DATASETS)

    bulk_rows_svg = "".join(
        render_bulk_row(idx, item, 22 + idx * 22)
        for idx, item in enumerate(BULK_TRIAL_COHORTS)
    )

    sc_milo_rows_svg = "".join(
        render_sc_milo_row(idx, item, 22 + idx * 22)
        for idx, item in enumerate(SC_MILO_COHORTS)
    )

    left_ref_rows_svg = "".join(
        render_deconv_ref_row(idx, item, 22 + idx * 22)
        for idx, item in enumerate(DECONV_REF_DATASETS[:10])
    )

    right_ref_rows_svg = "".join(
        render_deconv_ref_row(idx, item, 22 + idx * 22)
        for idx, item in enumerate(DECONV_REF_DATASETS[10:])
    )

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{SVG_WIDTH_MM}" height="{SVG_HEIGHT_MM}" viewBox="0 0 {VIEWBOX_W} {VIEWBOX_H}"
     version="1.1" id="figure-combined-datasets-compendium"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <!-- ====================================================================== -->
  <!-- LAYER 00: CANVAS BACKGROUND & OUTER PANEL CONTAINERS                   -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{VIEWBOX_W}" height="{VIEWBOX_H}" fill="{COLOR_CANVAS_BG}" />

    <!-- Panel A Container Card (Bulk ICB Trials) -->
    <rect id="card-panel-a" x="3" y="4" width="674" height="301"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel B Container Card (Single-Cell Milo Response Cohorts) -->
    <rect id="card-panel-b" x="3" y="309" width="674" height="274"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />

    <!-- Panel C Container Card (Deconvolution Reference & Integration Datasets) -->
    <rect id="card-panel-c" x="3" y="587" width="674" height="274"
          fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 01: PANEL BADGES                                                 -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Panel_Headers">
    <!-- Panel Badge [a] -->
    <g id="badge-panel-a" transform="translate(8, 8)">
      <text x="0" y="13" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(a)</text>
      <text x="20" y="13" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.15">Bulk RNA-seq Clinical Trial Cohorts: Regimens &amp; RECIST Response Rates (n = 1,097)</text>
    </g>

    <!-- Panel Badge [b] -->
    <g id="badge-panel-b" transform="translate(8, 313)">
      <text x="0" y="13" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(b)</text>
      <text x="20" y="13" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.15">Single-Cell ICB Response Cohorts: Milo Compendium (9 Cohorts, n = 613, {sc_total_cells:,} Cells)</text>
    </g>

    <!-- Panel Badge [c] -->
    <g id="badge-panel-c" transform="translate(8, 591)">
      <text x="0" y="13" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(c)</text>
      <text x="20" y="13" font-size="9.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.15">Single-Cell Reference Datasets in Deconvolution &amp; Integration Benchmark (19 Studies)</text>
    </g>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 02: PANEL A — BULK RNA-SEQ CLINICAL TRIAL COHORTS                -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_Bulk_ICB_Trials">
    <g id="grp-bulk-table" transform="translate(6, 29)">
      <!-- Outer Table Border -->
      <rect x="0" y="0" width="668" height="242" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />

      <!-- Header Bar -->
      <rect x="0" y="0" width="668" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
      <text x="5" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Trial Cohort</text>
      <text x="158" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Cancer</text>
      <text x="256" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Regimen &amp; Target</text>
      <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">N</text>
      <text x="495" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">RECIST (R / NR)</text>
      <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">ORR</text>
      {bulk_rows_svg}

      <!-- Overall Summary Row (y = 220, height = 22) -->
      <g id="bulk-row-summary" transform="translate(0, 220)">
        <rect x="0" y="0" width="668" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <text x="5" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Overall (9 Cohorts)</text>
        <text x="158" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">5 Solid Types</text>
        <text x="256" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">Anti-PD-(L)1 ± Combo</text>
        <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{bulk_total_n:,}</text>
        <g transform="translate(495, 3)">
          <rect x="0" y="0" width="{bulk_w_r}" height="16" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{bulk_w_r}" y="0" width="{bulk_w_nr}" height="16" fill="{OKABE_VERMILION}" />
          <text x="{bulk_w_r // 2}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{bulk_total_r}R</text>
          <text x="{bulk_w_r + (bulk_w_nr // 2)}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{bulk_total_nr}NR</text>
        </g>
        <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{bulk_orr:.1f}%</text>
      </g>

      <!-- RECIST Response Legend Row (y = 247, height = 23) -->
      <g id="grp-bulk-strips" transform="translate(0, 247)">
        <rect x="0" y="0" width="668" height="23" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.5" />
        <rect x="10" y="7.0" width="9" height="9" fill="{OKABE_BLUISH_GREEN}" />
        <text x="24" y="14.5" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">R: Responder (CR / PR, n = 357, 32.5%)</text>
        <rect x="344" y="7.0" width="9" height="9" fill="{OKABE_VERMILION}" />
        <text x="358" y="14.5" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">NR: Non-Responder (SD / PD, n = 740, 67.5%)</text>
      </g>
    </g>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 03: PANEL B — SINGLE-CELL MILO RESPONSE COHORTS                  -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Single_Cell_Milo_Cohorts">
    <g id="grp-sc-milo-table" transform="translate(6, 335)">
      <!-- Outer Table Border -->
      <rect x="0" y="0" width="668" height="242" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />

      <!-- Header Bar -->
      <rect x="0" y="0" width="668" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
      <text x="5" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Primary Study</text>
      <text x="158" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Cancer</text>
      <text x="256" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Regimen</text>
      <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">Cell Number</text>
      <text x="495" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Response (R / NR)</text>
      <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">ORR</text>
      {sc_milo_rows_svg}

      <!-- Single-Cell Milo Summary Row (y = 220, height = 22) -->
      <g id="sc-milo-row-summary" transform="translate(0, 220)">
        <rect x="0" y="0" width="668" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <text x="5" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Milo Total (9 Cohorts)</text>
        <text x="158" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">6 Solid Types</text>
        <text x="256" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">PD-(L)1 ± Combo</text>
        <text x="486" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">{sc_total_cells:,}</text>
        <g transform="translate(495, 3)">
          <rect x="0" y="0" width="{sc_w_r}" height="16" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{sc_w_r}" y="0" width="{sc_w_nr}" height="16" fill="{OKABE_VERMILION}" />
          <text x="{sc_w_r // 2}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{sc_total_r}R</text>
          <text x="{sc_w_r + (sc_w_nr // 2)}" y="10.8" font-size="8.0" font-weight="700" fill="#FFFFFF" text-anchor="middle" letter-spacing="-0.2">{sc_total_nr}NR</text>
        </g>
        <text x="663" y="14.0" font-size="8.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end" letter-spacing="-0.2">{sc_orr:.1f}%</text>
      </g>
    </g>
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 04: PANEL C — DECONVOLUTION REFERENCE & INTEGRATION DATASETS     -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Deconv_Reference_Datasets">
    <g id="grp-deconv-ref-content" transform="translate(6, 613)">

      <!-- Left Sub-Table: Datasets 01 to 10 (x = 0, width = 331, height = 242) -->
      <g id="grp-ref-table-left" transform="translate(0, 0)">
        <rect x="0" y="0" width="331" height="242" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <rect x="0" y="0" width="331" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <text x="4" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">#</text>
        <text x="22" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Study</text>
        <text x="168" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Cancer</text>
        <text x="326" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">Cell Number</text>
        {left_ref_rows_svg}
      </g>

      <!-- Right Sub-Table: Datasets 11 to 19 + Total Row (x = 337, width = 331, height = 242) -->
      <g id="grp-ref-table-right" transform="translate(337, 0)">
        <rect x="0" y="0" width="331" height="242" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <rect x="0" y="0" width="331" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
        <text x="4" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">#</text>
        <text x="22" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Study</text>
        <text x="168" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">Cancer</text>
        <text x="326" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">Cell Number</text>
        {right_ref_rows_svg}

        <!-- Row 10 of Right Table: Aggregate Summary of 19 Datasets (y = 220, height = 22) -->
        <g id="deconv-ref-row-total" transform="translate(0, 220)">
          <rect x="0" y="0" width="331" height="22" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.6" />
          <text x="4" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Σ</text>
          <text x="22" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" letter-spacing="-0.2">19 Reference Studies</text>
          <text x="168" y="14.0" font-size="8.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}" letter-spacing="-0.2">12 Types</text>
          <text x="326" y="14.0" font-size="8.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end" letter-spacing="-0.2">{ref_total_cells:,}</text>
        </g>
      </g>
    </g>
  </g>
</svg>
"""


def validate_and_export(svg_content: str, out_paths: tuple[Path, ...]) -> Result[tuple[Path, ...], str]:
    """Validate XML structure and export SVG and 300 DPI PNG files."""
    try:
        root = ET.fromstring(svg_content)
    except ET.ParseError as exc:
        return Failure(f"Invalid SVG XML generated: {exc}")

    ns_inkscape = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
    layers = tuple(elem for elem in root.iter() if elem.attrib.get(ns_inkscape) == "layer")
    if len(layers) < 5:
        return Failure(f"Expected 5 Inkscape layers, found {len(layers)}")

    # Verify every <text> element has physical font-size >= 7.5px (>= 5.6 pt)
    for elem in root.iter("{http://www.w3.org/2000/svg}text"):
        fs_str = elem.attrib.get("font-size", "0")
        fs_val = float(fs_str.rstrip("ptx"))
        if fs_val < 7.5:
            return Failure(f"Found text element with font-size {fs_str} < 7.5")

    written: list[Path] = []
    for svg_path in out_paths:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        svg_path.write_text(svg_content, encoding="utf-8")
        written.append(svg_path)

        png_path = svg_path.with_suffix(".png")
        from plotting_utils import rasterize_svg_to_png

        match rasterize_svg_to_png(svg_path, png_path, width_px=2126, background="white"):
            case Success(_):
                written.append(png_path)
            case Failure(err):
                return Failure(err)

    return Success(tuple(written))


def main() -> int:
    """Imperative shell entry point."""
    repo_root = Path(__file__).resolve().parents[2]
    targets = (
        repo_root / "article" / "figures" / "dataset_overview" / "figure_combined_datasets_overview.svg",
        repo_root / "output" / "reports" / "figure_combined_datasets_overview.svg",
    )

    res = validate_and_export(render_svg(), targets)
    if isinstance(res, Failure):
        sys.stderr.write(f"[ERROR] {res.failure()}\n")
        return 1
    paths = res.unwrap()
    for p in paths:
        size_bytes = p.stat().st_size
        sys.stdout.write(f"[OK] Generated {p.suffix.upper().lstrip('.')}: {p.relative_to(repo_root)} ({size_bytes:,} bytes)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
