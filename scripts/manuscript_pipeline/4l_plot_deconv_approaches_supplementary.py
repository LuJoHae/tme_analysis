#!/usr/bin/env python3
"""Generate publication-grade supplementary comparison figures across all Figure 4 deconvolution approaches.

Synthesizes results from:
  - Pseudobulk-calibrated cross-modality frameworks:
      * Option A (GMRF + DerSimonian-Laird precision weighting, r = 0.83)
      * Option B Unsmoothed (Unprojected 55% external bulk + 45% SC pseudobulk, r = 0.68)
      * Option B GMRF (Unprojected + GMRF resolvent smoothing, r = 0.67)
  - 100% Pure external bulk RNA-seq frameworks (11 clinical trials, N = 965 patients, 0% pseudobulk):
      * Baseline Option C (D0 GondalBank indirect transfer + GMRF, r = 0.20)
      * Native Unsupervised HVG / Markers (G0_Filtered_HVG, G1_Unsup_Markers)
      * Option C Alt 1 / Alt 1b (G2b Strictly Held-Out Cross-Atlas LOCO Milopy + NuSVR / Ridge)
      * Option C Alt 2 (G2a State-Blind Cohort-Shared Milopy + DWLS)
      * Option C Alt 3 / Alt 4 (G2d State-Adaptive Milopy +/- GMRF + DWLS)

Renders 4 full-page Nature Methods wireframe vector figures (.svg + 300-DPI .png):
  1. supp_fig_deconv_factorial_solver_gene_matrix.{svg,png}
  2. supp_fig_deconv_four_axis_validation_suite.{svg,png}
  3. supp_fig_deconv_approaches_head_to_head_comparison.{svg,png}
  4. supp_fig_deconv_gene_selection_mechanisms.{svg,png}
"""

from __future__ import annotations

import io
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final, Literal, assert_never, cast

import altair as alt
import numpy as np
import polars as pl
from PIL import Image
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.stats import pearsonr, spearmanr

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
for _p_str in (str(REPO_ROOT), str(SCRIPT_DIR)):
    if _p_str not in sys.path:
        sys.path.insert(0, _p_str)

from plotting_utils import (  # noqa: E402
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLACK,
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

# ------------------------------------------------------------------------------
# Canonical Labels (Strictly Stripped of M01-M42 Numeric State Prefixes)
# ------------------------------------------------------------------------------

META_STATE_CLEAN_LABELS: Final[dict[str, str]] = {
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

CANONICAL_POLE_STATES: Final[tuple[str, ...]] = (
    "M30",  # CR2+PAX5+ TLS B
    "M37",  # MS4A1+CD79A+ TLS B
    "M29",  # JCHAIN+MZB1+ Plasma
    "M24",  # CXCL13+NMB+ Tfh/Trm
    "M13",  # CXCL13+PDCD1+ Tex
    "M34",  # GZMK+NKG7+ Tem CD8+
    "M14",  # CXCL9+STAT1+ M1-TAM
    "M32",  # MRC1+CD163+ M2-TAM
    "M16",  # VEGFA+KDR+ Endo.
    "M19",  # S100A8+CXCL8+ MDSC
    "M12",  # CA9+SLC2A1+ Hypox.
    "M26",  # VIM+ZEB1+ EMT Tumor
    "M17",  # ISLR+COL5A2+ myCAF
    "M23",  # FAP+TGFB1+ Mat. CAF
)

SC_COHORT_ORDER: Final[tuple[str, ...]] = (
    "GSE120575",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE207422",
    "GSE243013",
    "GSE233203",
    "GSE200996",
    "GSE316195",
)

SC_COHORT_SHORT: Final[dict[str, str]] = {
    "GSE120575": "Sade-Feldman",
    "CELLxGENE_7b20c613": "Gondal",
    "CELLxGENE_05a8c945": "Marteau",
    "CELLxGENE_6f9de485": "Yan",
    "GSE207422": "Hu",
    "GSE243013": "Liu",
    "GSE233203": "Hong",
    "GSE200996": "Luoma",
    "GSE316195": "Raufi",
}

SC_COHORT_COLORS: Final[list[str]] = [
    OKABE_BLUE,
    OKABE_VERMILION,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_SKY_BLUE,
    OKABE_REDDISH_PURPLE,
    "#332288",
    "#882255",
    OKABE_GREY,
]

REF_13_ORDER: Final[list[str]] = [
    "Sade-Feldman",
    "Gondal",
    "Marteau",
    "Yan",
    "Hu",
    "Liu",
    "Hong",
    "Luoma",
    "Raufi",
    "Mel+PanCa",
    "NSCLC-Uni",
    "Unified-All9",
    "Meta-42",
]

REF_ID_TO_COMPACT: Final[dict[str, str]] = {
    "GSE120575": "Sade-Feldman",
    "CELLxGENE_7b20c613": "Gondal",
    "CELLxGENE_05a8c945": "Marteau",
    "CELLxGENE_6f9de485": "Yan",
    "GSE207422": "Hu",
    "GSE243013": "Liu",
    "GSE233203": "Hong",
    "GSE200996": "Luoma",
    "GSE316195": "Raufi",
    "Melanoma-PanCancer-Unified": "Mel+PanCa",
    "NSCLC-Unified": "NSCLC-Uni",
    "Unified-All9": "Unified-All9",
    "Consensus-Meta42": "Meta-42",
}

SOLVER_ID_TO_LABEL: Final[dict[str, str]] = {
    "D0_Baseline_GondalBank": "D0: GondalBank",
    "D1_InstaPrism_EM": "InstaPrism EM",
    "D2_DWLS": "DWLS",
    "D3_CIBERSORT_NuSVR": "CIBERSORT NuSVR",
    "D4_Ridge_LogContrast": "Ridge Log-Contrast",
    "D5_Ensemble": "4-Solver Ensemble",
}

SOLVER_5_ORDER: Final[list[str]] = [
    "InstaPrism EM",
    "DWLS",
    "CIBERSORT NuSVR",
    "Ridge Log-Contrast",
    "4-Solver Ensemble",
]

STRATEGY_ID_TO_LABEL: Final[dict[str, str]] = {
    "Baseline_D0_GondalBank": "I. D0: GondalBank + GMRF",
    "G0_Filtered_HVG": "I. Filtered HVG (G0)",
    "G1_Unsup_Markers": "I. Unsup. Specificity (G1)",
    "G2b_LOCO_Milo": "II. LOCO Milo (G2b)",
    "G2b_LOCO_Milo_GMRF": "II. LOCO Milo + GMRF (Alt 1)",
    "G2a_Cohort_Milo": "III. Cohort Milo (G2a)",
    "G2a_Cohort_Milo_GMRF": "III. Cohort Milo + GMRF (Alt 2)",
    "G2d_State_Milo_GMRF": "IV. State Milo + GMRF (Alt 3)",
    "G2d_State_Milo": "IV. State Milo Direct (Alt 4)",
}

STRATEGY_8_ORDER: Final[list[str]] = [
    "I. Filtered HVG (G0)",
    "I. Unsup. Specificity (G1)",
    "II. LOCO Milo (G2b)",
    "II. LOCO Milo + GMRF (Alt 1)",
    "III. Cohort Milo (G2a)",
    "III. Cohort Milo + GMRF (Alt 2)",
    "IV. State Milo + GMRF (Alt 3)",
    "IV. State Milo Direct (Alt 4)",
]

HygieneTier = Literal["Unsupervised", "LOCO Cross-Atlas", "Cohort-Shared", "State-Adaptive", "Calibrated"]


# ------------------------------------------------------------------------------
# Immutable Configuration & Loaded Bundle Models
# ------------------------------------------------------------------------------


class SuppComparisonConfig(BaseModel):
    """Immutable configuration for generating the 4 deconvolution supplementary comparison figures."""

    model_config = ConfigDict(frozen=True)

    explore_4k_dir: Path = Path("output/manuscript_pipeline/4k_explore_option_c_deconv_alternatives")
    fig4_cache_dir: Path = Path("output/manuscript_pipeline/4i_plot_combined_main_figure")
    stage3b_dir: Path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark")
    out_dir: Path = Path("output/manuscript_pipeline/4l_plot_deconv_approaches_supplementary")
    article_fig_dir: Path = Path("article/figures/milo-deconv")
    width_mm: float = 180.0


class ParadigmSpec(BaseModel):
    """Immutable descriptor for one flagship cross-modality deconvolution paradigm."""

    model_config = ConfigDict(frozen=True)

    key: str
    cache_prefix: str
    row_label: str
    compact_label: str
    scatter_title: str
    regime_group: Literal["Pseudobulk-Calibrated", "Pure External Bulk"]
    hygiene_tier: HygieneTier
    include_in_scatter8: bool


FLAGSHIP_PARADIGMS: Final[tuple[ParadigmSpec, ...]] = (
    ParadigmSpec(
        key="opt_a",
        cache_prefix="option_a_gmrf_dl",
        row_label="[Calibrated] Option A: GMRF + DL",
        compact_label="Option A (GMRF+DL)",
        scatter_title="(a1) Option A: GMRF + DL",
        regime_group="Pseudobulk-Calibrated",
        hygiene_tier="Calibrated",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="opt_b_raw",
        cache_prefix="unprojected_unsmoothed",
        row_label="[Calibrated] Option B: Unsmoothed",
        compact_label="Option B (Unsmoothed)",
        scatter_title="(a2) Option B: Unsmoothed",
        regime_group="Pseudobulk-Calibrated",
        hygiene_tier="Calibrated",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="opt_b_gmrf",
        cache_prefix="unprojected",
        row_label="[Calibrated] Option B: GMRF Smooth",
        compact_label="Option B (GMRF)",
        scatter_title="Option B: Unprojected GMRF",
        regime_group="Pseudobulk-Calibrated",
        hygiene_tier="Calibrated",
        include_in_scatter8=False,
    ),
    ParadigmSpec(
        key="d0_base",
        cache_prefix="principled_gmrf",
        row_label="[Pure Bulk] Baseline D0: GondalBank",
        compact_label="Option C: D0 Baseline",
        scatter_title="(a3) Option C: D0 Baseline",
        regime_group="Pure External Bulk",
        hygiene_tier="Unsupervised",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="g1_unsup",
        cache_prefix="option_c_unsup_g1_instaprism",
        row_label="[Pure Bulk] Unsup. G1 + InstaPrism",
        compact_label="Unsup. G1 + InstaPrism",
        scatter_title="(a4) Unsup. G1 + InstaPrism",
        regime_group="Pure External Bulk",
        hygiene_tier="Unsupervised",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="alt1_nusvr",
        cache_prefix="option_c_alt1_loco_milo_nusvr",
        row_label="[Pure Bulk] Alt 1: LOCO + NuSVR",
        compact_label="Alt 1: LOCO + NuSVR",
        scatter_title="(a5) Alt 1: LOCO + NuSVR",
        regime_group="Pure External Bulk",
        hygiene_tier="LOCO Cross-Atlas",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="alt1b_ridge",
        cache_prefix="option_c_alt1b_loco_milo_ridge",
        row_label="[Pure Bulk] Alt 1b: LOCO + Ridge",
        compact_label="Alt 1b: LOCO + Ridge",
        scatter_title="Option C Alt 1b: LOCO + Ridge",
        regime_group="Pure External Bulk",
        hygiene_tier="LOCO Cross-Atlas",
        include_in_scatter8=False,
    ),
    ParadigmSpec(
        key="alt2_cohort",
        cache_prefix="option_c_alt2_cohort_milo_dwls",
        row_label="[Pure Bulk] Alt 2: Cohort Milo + DWLS",
        compact_label="Alt 2: Cohort + DWLS",
        scatter_title="(a6) Alt 2: Cohort + DWLS",
        regime_group="Pure External Bulk",
        hygiene_tier="Cohort-Shared",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="alt3_state_gmrf",
        cache_prefix="option_c_alt3_state_milo_gmrf_dwls",
        row_label="[Pure Bulk] Alt 3: State Milo + GMRF",
        compact_label="Alt 3: State + GMRF",
        scatter_title="(a7) Alt 3: State + GMRF",
        regime_group="Pure External Bulk",
        hygiene_tier="State-Adaptive",
        include_in_scatter8=True,
    ),
    ParadigmSpec(
        key="alt4_state_raw",
        cache_prefix="option_c_alt4_state_milo_dwls",
        row_label="[Pure Bulk] Alt 4: State Milo + DWLS",
        compact_label="Alt 4: State + DWLS",
        scatter_title="(a8) Alt 4: State + DWLS",
        regime_group="Pure External Bulk",
        hygiene_tier="State-Adaptive",
        include_in_scatter8=True,
    ),
)


class LoadedComparisonBundle(BaseModel):
    """Immutable bundle of all validated Polars DataFrames required for the 4 supplementary figures."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    grid_bench_df: pl.DataFrame
    loto_folds_df: pl.DataFrame
    perm_null_df: pl.DataFrame
    gene_space_df: pl.DataFrame
    cluster_meta_df: pl.DataFrame
    paradigm_states: dict[str, pl.DataFrame]
    paradigm_conc: dict[str, pl.DataFrame]


# ------------------------------------------------------------------------------
# Pure Data Loading & Validation Helpers
# ------------------------------------------------------------------------------


def classify_strategy_tier(strategy_id: str) -> HygieneTier:
    """Exhaustively classify a strategy identifier into its information-hygiene tier."""
    match strategy_id:
        case "Baseline_D0_GondalBank" | "G0_Filtered_HVG" | "G1_Unsup_Markers":
            return "Unsupervised"
        case "G2b_LOCO_Milo" | "G2b_LOCO_Milo_GMRF":
            return "LOCO Cross-Atlas"
        case "G2a_Cohort_Milo" | "G2a_Cohort_Milo_GMRF":
            return "Cohort-Shared"
        case "G2d_State_Milo_GMRF" | "G2d_State_Milo":
            return "State-Adaptive"
        case _:
            return "Unsupervised"


def tier_color(tier: HygieneTier) -> str:
    """Exhaustively map an information-hygiene tier to its Okabe-Ito colorblind-safe color."""
    match tier:
        case "Unsupervised":
            return OKABE_GREY
        case "LOCO Cross-Atlas":
            return OKABE_BLUISH_GREEN
        case "Cohort-Shared":
            return OKABE_BLUE
        case "State-Adaptive":
            return OKABE_ORANGE
        case "Calibrated":
            return OKABE_REDDISH_PURPLE
        case _ as unreachable:
            assert_never(unreachable)


def strip_numeric_state_prefix(raw_label: str) -> str:
    """Remove any M01-M42 or L0-L99 numeric state prefix while preserving biological markers."""
    if raw_label in META_STATE_CLEAN_LABELS:
        return META_STATE_CLEAN_LABELS[raw_label]
    cleaned = re.sub(r"^M\d{2}\s*[:·\-]\s*", "", raw_label.strip())
    return re.sub(r"^M\d{2}\s+", "", cleaned)


def load_all_comparison_data(config: SuppComparisonConfig) -> Result[LoadedComparisonBundle, str]:
    """Load and validate all granular benchmark and paradigm Parquet tables."""
    required_4k = (
        config.explore_4k_dir / "option_c_alternatives_grid_benchmark.parquet",
        config.explore_4k_dir / "option_c_loto_trial_folds.parquet",
        config.explore_4k_dir / "option_c_permutation_null_samples.parquet",
        config.explore_4k_dir / "option_c_gene_space_concordance.parquet",
        config.stage3b_dir / "unified_all_states_detailed_concordance.parquet",
    )
    for req_p in required_4k:
        if not req_p.exists():
            return Failure(f"Missing required benchmark table: {req_p}")

    paradigm_states: dict[str, pl.DataFrame] = {}
    paradigm_conc: dict[str, pl.DataFrame] = {}
    for spec in FLAGSHIP_PARADIGMS:
        st_p = config.fig4_cache_dir / f"{spec.cache_prefix}_state_effects.parquet"
        co_p = config.fig4_cache_dir / f"{spec.cache_prefix}_concordance_summary.parquet"
        if not st_p.exists() or not co_p.exists():
            return Failure(f"Missing cached Figure 4 paradigm tables for {spec.cache_prefix}: {st_p}")
        paradigm_states[spec.key] = pl.read_parquet(st_p)
        paradigm_conc[spec.key] = pl.read_parquet(co_p)

    return Success(
        LoadedComparisonBundle(
            grid_bench_df=pl.read_parquet(required_4k[0]),
            loto_folds_df=pl.read_parquet(required_4k[1]),
            perm_null_df=pl.read_parquet(required_4k[2]),
            gene_space_df=pl.read_parquet(required_4k[3]),
            cluster_meta_df=pl.read_parquet(required_4k[4]),
            paradigm_states=paradigm_states,
            paradigm_conc=paradigm_conc,
        )
    )


# ------------------------------------------------------------------------------
# Nature Methods Styling & Inkscape Vector Export Helpers
# ------------------------------------------------------------------------------


def style_publication_chart(chart: Any) -> Any:
    """Apply Nature Methods minimalist wireframe styling with outward ticks and zero drop shadows."""
    return (
        chart.configure(
            background=COLOR_BG_WHITE,
            font="system-ui, -apple-system, sans-serif",
            padding={"left": 28, "right": 28, "top": 14, "bottom": 14},
        )
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75, fill=COLOR_BG_WHITE)
        .configure_axis(
            labelColor=COLOR_TEXT_PRIMARY,
            titleColor=COLOR_TEXT_PRIMARY,
            labelFontSize=9.5,
            titleFontSize=10.5,
            titleFontWeight="bold",
            tickColor="#1E293B",
            domainColor="#1E293B",
            domainWidth=0.75,
            tickWidth=0.75,
            tickSize=4,
            gridColor="#F1F5F9",
            labelLimit=420,
            titlePadding=10,
        )
        .configure_legend(
            labelColor=COLOR_TEXT_PRIMARY,
            titleColor=COLOR_TEXT_PRIMARY,
            labelFontSize=9.2,
            titleFontSize=9.8,
            titleFontWeight="bold",
            strokeColor=COLOR_HAIRLINE,
            fillColor=COLOR_BG_WHITE,
            padding=6,
            columnPadding=16,
            labelLimit=360,
        )
        .configure_title(
            color=COLOR_TEXT_PRIMARY,
            fontSize=11.0,
            fontWeight="bold",
            anchor="start",
        )
    )


def _alt_data(*, values: Sequence[dict[str, Any]]) -> Any:
    """Typed wrapper around Altair's inline Data constructor for strict mypy compliance."""
    return getattr(alt, "Data")(values=list(values))


def finalize_and_verify_svg(
    chart: Any,
    out_dir: Path,
    article_dir: Path,
    stem: str,
    width_mm: float = 180.0,
    max_vector_points: int = 5000,
) -> Result[tuple[Path, Path], str]:
    """Export Altair chart to SVG, inject Inkscape layers, tight-crop, validate XML, and rasterize 300-DPI PNG."""
    svg_path, _ = export_altair_figure(chart, out_dir / stem, scale=3.0, max_vector_points=max_vector_points)
    raw_txt = svg_path.read_text(encoding="utf-8")

    m = re.search(r'viewBox="([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.]+)\s+([0-9.]+)"', raw_txt)
    if m is not None:
        vx0, vy0, vw0, vh0 = (float(m.group(i)) for i in range(1, 5))
        png_probe = subprocess.check_output(
            ["rsvg-convert", "-w", str(int(round(vw0))), "-h", str(int(round(vh0))), "-b", "white", str(svg_path)]
        )
        arr = np.asarray(Image.open(io.BytesIO(png_probe)).convert("RGB"))
        non_white = np.any(arr < 252, axis=2)
        if np.any(non_white):
            ys, xs = np.where(non_white)
            pad = 8.0
            nx = vx0 + max(0.0, float(xs.min()) - pad)
            ny = vy0 + max(0.0, float(ys.min()) - pad)
            nw = min(vw0 - (nx - vx0), float(xs.max() - xs.min() + 1) + 2.0 * pad)
            nh = min(vh0 - (ny - vy0), float(ys.max() - ys.min() + 1) + 2.0 * pad)
            h_mm = width_mm * (nh / max(nw, 1.0))

            def _inject_root_attrs(match: re.Match[str]) -> str:
                tag = match.group(0)
                tag = re.sub(r'\s+width="[^"]*"', "", tag)
                tag = re.sub(r'\s+height="[^"]*"', "", tag)
                if "xmlns:inkscape=" not in tag:
                    tag = tag.replace(
                        "<svg",
                        '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
                        'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"',
                        1,
                    )
                return re.sub(
                    r'viewBox="[^"]*"',
                    f'width="{width_mm:.2f}mm" height="{h_mm:.2f}mm" viewBox="{nx:.1f} {ny:.1f} {nw:.1f} {nh:.1f}"',
                    tag,
                )

            raw_txt = re.sub(r"<svg[^>]+>", _inject_root_attrs, raw_txt, count=1)
            if 'inkscape:groupmode="layer"' not in raw_txt:
                raw_txt = raw_txt.replace(
                    '<g transform="translate(',
                    '<g inkscape:groupmode="layer" id="layer-01-panels" inkscape:label="01_Publication_Panels" transform="translate(',
                    1,
                )
            svg_path.write_text(raw_txt, encoding="utf-8")

    try:
        ET.parse(svg_path)
    except ET.ParseError as exc:
        return Failure(f"Malformed XML in {svg_path}: {exc}")

    if "feDropShadow" in raw_txt:
        return Failure(f"Prohibited feDropShadow detected in {svg_path}")

    png_path = svg_path.with_suffix(".png")
    subprocess.run(
        ["rsvg-convert", "-d", "300", "-p", "300", "-b", "white", str(svg_path), "-o", str(png_path)],
        check=True,
    )

    article_dir.mkdir(parents=True, exist_ok=True)
    dst_svg = article_dir / svg_path.name
    dst_png = article_dir / png_path.name
    shutil.copy2(svg_path, dst_svg)
    shutil.copy2(png_path, dst_png)
    return Success((dst_svg, dst_png))


# ------------------------------------------------------------------------------
# Supplementary Figure 1: Factorial Solver x Gene Matrix & 10x13 Cohort Heatmap
# ------------------------------------------------------------------------------


def build_supp_fig1_factorial_matrix(bundle: LoadedComparisonBundle) -> Any:
    """Build Supplementary Figure 1: (a) 8x5 factorial matrix, (b) 355 vs Meta-42 slope, (c) 10x13 cohort heatmap."""
    d0_row = bundle.grid_bench_df.filter(pl.col("strategy_id") == "Baseline_D0_GondalBank")
    d0_r = float(d0_row["global_pearson_r"][0])
    d0_rho = float(d0_row["global_spearman_rho"][0])

    # Panel (a): 8 strategies x 5 direct deconvolution solvers
    mat_df = (
        bundle.grid_bench_df.lazy()
        .filter(pl.col("strategy_id") != "Baseline_D0_GondalBank")
        .with_columns(
            pl.col("strategy_id").replace(STRATEGY_ID_TO_LABEL).alias("strategy_label"),
            pl.col("deconv_tool").replace(SOLVER_ID_TO_LABEL).alias("solver_label"),
            (
                pl.col("global_pearson_r").map_elements(lambda v: f"{v:.2f}", return_dtype=pl.String)
                + " ("
                + pl.col("positive_cohorts").cast(pl.String)
                + "/9)"
            ).alias("cell_text"),
            pl.when(pl.col("global_pearson_r") >= 0.42)
            .then(pl.lit("#FFFFFF"))
            .otherwise(pl.lit("#0F172A"))
            .alias("text_color"),
        )
        .collect()
    )
    mat_data = _alt_data(values=mat_df.to_dicts())

    pa_rect = (
        alt.Chart(mat_data)
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.2)
        .encode(
            x=alt.X(
                "solver_label:N",
                sort=SOLVER_5_ORDER,
                title="Direct Bulk Deconvolution Solver",
                axis=alt.Axis(labelAngle=-22),
            ),
            y=alt.Y("strategy_label:N", sort=STRATEGY_8_ORDER, title="Single-Cell Gene-Selection & GMRF Tier"),
            color=alt.Color(
                "global_pearson_r:Q",
                title="Global r",
                scale=alt.Scale(
                    domain=[0.08, 0.25, 0.42, 0.65],
                    range=["#F8FAFC", "#BAE6FD", OKABE_SKY_BLUE, OKABE_BLUE],
                ),
                legend=alt.Legend(orient="right", gradientLength=120),
            ),
        )
    )
    pa_txt_dark = (
        alt.Chart(_alt_data(values=mat_df.filter(pl.col("global_pearson_r") < 0.42).to_dicts()))
        .mark_text(fontSize=8.5, fontWeight="bold", color="#0F172A")
        .encode(
            x=alt.X("solver_label:N", sort=SOLVER_5_ORDER),
            y=alt.Y("strategy_label:N", sort=STRATEGY_8_ORDER),
            text=alt.Text("cell_text:N"),
        )
    )
    pa_txt_light = (
        alt.Chart(_alt_data(values=mat_df.filter(pl.col("global_pearson_r") >= 0.42).to_dicts()))
        .mark_text(fontSize=8.5, fontWeight="bold", color="#FFFFFF")
        .encode(
            x=alt.X("solver_label:N", sort=SOLVER_5_ORDER),
            y=alt.Y("strategy_label:N", sort=STRATEGY_8_ORDER),
            text=alt.Text("cell_text:N"),
        )
    )
    panel_a = alt.layer(pa_rect, pa_txt_dark, pa_txt_light).properties(
        width=310,
        height=225,
        title=alt.TitleParams(
            text="(a) Pure External Bulk Factorial Grid (8 Gene Tiers × 5 Solvers)",
            subtitle=f"Cell text: Global Pearson r (Pos. SC cohorts/9) | Indirect D0 Baseline: r = {d0_r:+.2f}, ρ = {d0_rho:+.2f}",
            subtitleFontSize=9.0,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (b): Granularity scaling Unified-All9 (K=355) vs Consensus-Meta42 (K=42)
    gran_order = ["I. D0: GondalBank + GMRF", *STRATEGY_8_ORDER]
    dwls_or_d0 = bundle.grid_bench_df.filter(
        (pl.col("strategy_id") == "Baseline_D0_GondalBank") | (pl.col("deconv_tool") == "D2_DWLS")
    ).with_columns(
        pl.col("strategy_id").replace(STRATEGY_ID_TO_LABEL).alias("strategy_label"),
        (pl.col("meta42_r") - pl.col("global_pearson_r")).alias("delta_r"),
    )
    gran_pts_rows: list[dict[str, object]] = []
    for row in dwls_or_d0.iter_rows(named=True):
        gran_pts_rows.append(
            {
                "strategy_label": row["strategy_label"],
                "granularity": "Fine-Grained All-9 (355 states)",
                "r_val": float(row["global_pearson_r"]),
            }
        )
        gran_pts_rows.append(
            {
                "strategy_label": row["strategy_label"],
                "granularity": "Consensus Meta-42 (42 states)",
                "r_val": float(row["meta42_r"]),
            }
        )
    pb_rules = (
        alt.Chart(_alt_data(values=dwls_or_d0.to_dicts()))
        .mark_rule(stroke="#64748B", strokeWidth=1.5)
        .encode(
            y=alt.Y("strategy_label:N", sort=gran_order, title=""),
            x=alt.X("global_pearson_r:Q", title="Cross-Modality Pearson r", scale=alt.Scale(domain=[0.0, 0.96])),
            x2=alt.X2("meta42_r:Q"),
        )
    )
    pb_pts = (
        alt.Chart(_alt_data(values=gran_pts_rows))
        .mark_point(filled=True, size=68, stroke="#FFFFFF", strokeWidth=0.8)
        .encode(
            y=alt.Y("strategy_label:N", sort=gran_order, axis=alt.Axis(labels=True)),
            x=alt.X("r_val:Q"),
            color=alt.Color(
                "granularity:N",
                title="Reference Granularity",
                scale=alt.Scale(
                    domain=["Fine-Grained All-9 (355 states)", "Consensus Meta-42 (42 states)"],
                    range=[OKABE_SKY_BLUE, OKABE_BLUE],
                ),
                legend=alt.Legend(orient="bottom", columns=1),
            ),
            shape=alt.Shape(
                "granularity:N",
                scale=alt.Scale(
                    domain=["Fine-Grained All-9 (355 states)", "Consensus Meta-42 (42 states)"],
                    range=["circle", "diamond"],
                ),
                legend=alt.Legend(orient="bottom", columns=1),
            ),
        )
    )
    pb_delta = (
        alt.Chart(
            _alt_data(
                values=dwls_or_d0.with_columns(
                    (pl.lit("Δ") + pl.col("delta_r").map_elements(lambda v: f"{v:+.2f}", return_dtype=pl.String)).alias(
                        "delta_txt"
                    )
                ).to_dicts()
            )
        )
        .mark_text(align="left", dx=6, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_SECONDARY)
        .encode(
            y=alt.Y("strategy_label:N", sort=gran_order),
            x=alt.X("meta42_r:Q"),
            text=alt.Text("delta_txt:N"),
        )
    )
    panel_b = alt.layer(pb_rules, pb_pts, pb_delta).properties(
        width=225,
        height=225,
        title=alt.TitleParams(
            text="(b) Consensus Aggregation Lift (355 → 42 States)",
            subtitle="Held-constant DWLS solver (+ D0 baseline) across tiers",
            subtitleFontSize=9.0,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (c): 10 Flagship Paradigms x 13 Single-Cell References Concordance Matrix
    conc_rows: list[dict[str, object]] = []
    paradigm_row_order = [s.row_label for s in FLAGSHIP_PARADIGMS]
    for spec in FLAGSHIP_PARADIGMS:
        c_df = bundle.paradigm_conc[spec.key]
        for r_row in c_df.iter_rows(named=True):
            ref_raw = str(r_row["reference_id"])
            if ref_raw not in REF_ID_TO_COMPACT:
                continue
            r_val = float(r_row["pearson_r"])
            conc_rows.append(
                {
                    "paradigm_label": spec.row_label,
                    "regime_group": spec.regime_group,
                    "ref_short": REF_ID_TO_COMPACT[ref_raw],
                    "pearson_r": r_val,
                    "r_txt": f"{r_val:+.2f}",
                }
            )
    conc_all_df = pl.DataFrame(conc_rows)
    pc_rect = (
        alt.Chart(_alt_data(values=conc_all_df.to_dicts()))
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.0)
        .encode(
            x=alt.X(
                "ref_short:N",
                sort=REF_13_ORDER,
                title="Single-Cell Cohort & Hierarchical Consensus Reference (13 Tiers)",
                axis=alt.Axis(labelAngle=-28),
            ),
            y=alt.Y("paradigm_label:N", sort=paradigm_row_order, title=""),
            color=alt.Color(
                "pearson_r:Q",
                title="Pearson r",
                scale=alt.Scale(
                    domain=[-0.35, 0.0, 0.45, 0.90],
                    range=[OKABE_VERMILION, "#F8FAFC", OKABE_SKY_BLUE, OKABE_BLUE],
                ),
                legend=alt.Legend(orient="right", gradientLength=145),
            ),
        )
    )
    pc_txt_dark = (
        alt.Chart(_alt_data(values=conc_all_df.filter(pl.col("pearson_r").abs() < 0.52).to_dicts()))
        .mark_text(fontSize=8.3, fontWeight="bold", color="#0F172A")
        .encode(
            x=alt.X("ref_short:N", sort=REF_13_ORDER),
            y=alt.Y("paradigm_label:N", sort=paradigm_row_order),
            text=alt.Text("r_txt:N"),
        )
    )
    pc_txt_light = (
        alt.Chart(_alt_data(values=conc_all_df.filter(pl.col("pearson_r").abs() >= 0.52).to_dicts()))
        .mark_text(fontSize=8.3, fontWeight="bold", color="#FFFFFF")
        .encode(
            x=alt.X("ref_short:N", sort=REF_13_ORDER),
            y=alt.Y("paradigm_label:N", sort=paradigm_row_order),
            text=alt.Text("r_txt:N"),
        )
    )
    panel_c = alt.layer(pc_rect, pc_txt_dark, pc_txt_light).properties(
        width=635,
        height=245,
        title=alt.TitleParams(
            text="(c) Cross-Paradigm × 13-Reference Concordance Matrix: Pseudobulk-Calibrated (Options A/B) vs. Pure External Bulk (Option C & Alts 1–4)",
            subtitle="Resolves negative correlations in FACS-sorted cohorts (Sade-Feldman CD45+, Luoma CD3+/CD45+) and PDAC (Raufi) without bulk label leakage",
            subtitleFontSize=9.0,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    top_row = alt.hconcat(panel_a, panel_b, spacing=24).resolve_scale(color="independent", shape="independent")
    return style_publication_chart(alt.vconcat(top_row, panel_c, spacing=24).resolve_scale(color="independent"))


# ------------------------------------------------------------------------------
# Supplementary Figure 2: Four-Axis Methodological Validation Suite (V1-V4)
# ------------------------------------------------------------------------------


def build_supp_fig2_validation_suite(bundle: LoadedComparisonBundle) -> Any:
    """Build Supplementary Figure 2: V1 dual permutation nulls, V3 LOTO jackknife, V2 PB CLR recovery, V4 LOTO AUC."""
    flagship_order = [
        "Option C Baseline (D0 GondalBank)",
        "Filtered HVG (G0 + DWLS)",
        "Unsupervised Markers (G1 + InstaPrism)",
        "Alt 1: LOCO Milo + NuSVR",
        "Alt 1b: LOCO Milo + Ridge",
        "Alt 2: Cohort Milo + DWLS",
        "Alt 3: State Milo + GMRF + DWLS",
        "Alt 4: State Milo + DWLS",
    ]

    # Panel (a): Validation V1 — Bulk patient-label permutation null (B=1,000) + Active state-label null
    v1_summary = (
        bundle.perm_null_df.lazy()
        .filter(pl.col("null_type") == "bulk_patient_perm")
        .group_by("paradigm_label")
        .agg(
            pl.col("null_r").quantile(0.025).alias("q025"),
            pl.col("null_r").quantile(0.975).alias("q975"),
            pl.col("null_r").mean().alias("null_mean"),
            pl.col("observed_r").first().alias("obs_r"),
            pl.col("perm_p_value").first().alias("perm_p"),
            pl.col("state_perm_mean_r").first().alias("state_null_mean"),
        )
        .with_columns(
            pl.col("perm_p")
            .map_elements(lambda p: f"p={p:.3f}", return_dtype=pl.String)
            .alias("p_label")
        )
        .collect()
    )
    v1_pts_rows: list[dict[str, object]] = []
    for row in v1_summary.iter_rows(named=True):
        v1_pts_rows.append(
            {
                "paradigm_label": row["paradigm_label"],
                "metric_type": "Bulk Patient-Label Null Mean (V1a)",
                "r_val": float(row["null_mean"]),
            }
        )
        v1_pts_rows.append(
            {
                "paradigm_label": row["paradigm_label"],
                "metric_type": "Active Target State-Shuffle Null (V1b)",
                "r_val": float(row["state_null_mean"]),
            }
        )
        v1_pts_rows.append(
            {
                "paradigm_label": row["paradigm_label"],
                "metric_type": "Observed External Bulk Concordance",
                "r_val": float(row["obs_r"]),
            }
        )

    v1_zero = (
        alt.Chart(_alt_data(values=[{"x0": 0.0}]))
        .mark_rule(stroke="#94A3B8", strokeDash=[3, 3], strokeWidth=1.0)
        .encode(x="x0:Q")
    )
    v1_ci = (
        alt.Chart(_alt_data(values=v1_summary.to_dicts()))
        .mark_rule(strokeWidth=6.5, color="#CBD5E1", opacity=0.9)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order, title=""),
            x=alt.X(
                "q025:Q",
                title="Global Cross-Modality Pearson r",
                scale=alt.Scale(domain=[-0.66, 0.88]),
                axis=alt.Axis(values=[-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]),
            ),
            x2=alt.X2("q975:Q"),
        )
    )
    v1_marks = (
        alt.Chart(_alt_data(values=v1_pts_rows))
        .mark_point(filled=True, size=68, stroke="#FFFFFF", strokeWidth=0.7)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order),
            x=alt.X("r_val:Q"),
            color=alt.Color(
                "metric_type:N",
                title="Validation V1 Distribution",
                scale=alt.Scale(
                    domain=[
                        "Bulk Patient-Label Null Mean (V1a)",
                        "Active Target State-Shuffle Null (V1b)",
                        "Observed External Bulk Concordance",
                    ],
                    range=[OKABE_GREY, OKABE_ORANGE, OKABE_BLUE],
                ),
                legend=alt.Legend(orient="top", columns=1),
            ),
            shape=alt.Shape(
                "metric_type:N",
                scale=alt.Scale(
                    domain=[
                        "Bulk Patient-Label Null Mean (V1a)",
                        "Active Target State-Shuffle Null (V1b)",
                        "Observed External Bulk Concordance",
                    ],
                    range=["square", "diamond", "circle"],
                ),
                legend=alt.Legend(orient="top", columns=1),
            ),
        )
    )
    v1_plbl = (
        alt.Chart(_alt_data(values=v1_summary.to_dicts()))
        .mark_text(align="left", dx=7, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order),
            x=alt.X("obs_r:Q"),
            text=alt.Text("p_label:N"),
        )
    )
    panel_a = alt.layer(v1_zero, v1_ci, v1_marks, v1_plbl).properties(
        width=285,
        height=210,
        title=alt.TitleParams(
            text="(a) Validation V1: Bulk & State Permutation Nulls",
            subtitle="Grey bar: 95% bulk label null (B = 1,000) | Diamond: active state-shuffle null",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Map flagship pairs for Panels (b) and (d)
    key_to_lbl: dict[tuple[str, str], str] = {
        ("Baseline_D0_GondalBank", "D0_Baseline_GondalBank"): "Option C Baseline (D0 GondalBank)",
        ("G0_Filtered_HVG", "D2_DWLS"): "Filtered HVG (G0 + DWLS)",
        ("G1_Unsup_Markers", "D1_InstaPrism_EM"): "Unsupervised Markers (G1 + InstaPrism)",
        ("G2b_LOCO_Milo_GMRF", "D3_CIBERSORT_NuSVR"): "Alt 1: LOCO Milo + NuSVR",
        ("G2b_LOCO_Milo_GMRF", "D4_Ridge_LogContrast"): "Alt 1b: LOCO Milo + Ridge",
        ("G2a_Cohort_Milo_GMRF", "D2_DWLS"): "Alt 2: Cohort Milo + DWLS",
        ("G2d_State_Milo_GMRF", "D2_DWLS"): "Alt 3: State Milo + GMRF + DWLS",
        ("G2d_State_Milo", "D2_DWLS"): "Alt 4: State Milo + DWLS",
    }
    folds_sub_rows: list[dict[str, object]] = []
    for row in bundle.loto_folds_df.iter_rows(named=True):
        pair = (str(row["strategy_id"]), str(row["deconv_tool"]))
        if pair in key_to_lbl:
            folds_sub_rows.append(
                {
                    **row,
                    "paradigm_label": key_to_lbl[pair],
                    "hygiene_tier": classify_strategy_tier(str(row["strategy_id"])),
                }
            )
    folds_sub_df = pl.DataFrame(folds_sub_rows)
    folds_range_df = (
        folds_sub_df.group_by("paradigm_label", "hygiene_tier")
        .agg(
            pl.col("loto_holdout_global_r").min().alias("min_r"),
            pl.col("loto_holdout_global_r").max().alias("max_r"),
            pl.col("full_global_r").first().alias("full_r"),
            pl.col("weighted_mean_auc").first().alias("mean_auc"),
        )
        .with_columns(
            pl.col("mean_auc").map_elements(lambda v: f"{v:.3f}", return_dtype=pl.String).alias("auc_txt")
        )
    )

    # Panel (b): Validation V3 — Leave-One-Trial-Out (LOTO) Bulk Jackknife Across 11 Trials
    v3_rule = (
        alt.Chart(_alt_data(values=folds_range_df.to_dicts()))
        .mark_rule(stroke="#334155", strokeWidth=2.0)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order, title=""),
            x=alt.X(
                "min_r:Q",
                title="Leave-One-Trial-Out (LOTO) Global Pearson r (11 Folds)",
                scale=alt.Scale(domain=[0.05, 0.72]),
            ),
            x2=alt.X2("max_r:Q"),
        )
    )
    v3_dots = (
        alt.Chart(_alt_data(values=folds_sub_df.to_dicts()))
        .mark_circle(size=42, opacity=0.78, stroke="#FFFFFF", strokeWidth=0.5)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order, axis=alt.Axis(labels=False)),
            x=alt.X("loto_holdout_global_r:Q"),
            color=alt.Color(
                "hygiene_tier:N",
                title="Information-Hygiene Tier",
                scale=alt.Scale(
                    domain=["Unsupervised", "LOCO Cross-Atlas", "Cohort-Shared", "State-Adaptive"],
                    range=[OKABE_GREY, OKABE_BLUISH_GREEN, OKABE_BLUE, OKABE_ORANGE],
                ),
                legend=alt.Legend(orient="top", columns=2),
            ),
        )
    )
    v3_full = (
        alt.Chart(_alt_data(values=folds_range_df.to_dicts()))
        .mark_point(shape="diamond", filled=True, size=65, color="#0F172A")
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order),
            x=alt.X("full_r:Q"),
        )
    )
    panel_b = alt.layer(v3_rule, v3_dots, v3_full).properties(
        width=285,
        height=210,
        title=alt.TitleParams(
            text="(b) Validation V3: 11-Trial LOTO Bulk Jackknife Stability",
            subtitle="Points: each dropped clinical trial (Hugo–Cho) | Diamond: all-11-trial r",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (c): Validation V2 — 613-Sample SC Pseudobulk CLR Recovery vs. External Bulk Concordance
    v2_df = bundle.grid_bench_df.with_columns(
        pl.col("deconv_tool").replace(SOLVER_ID_TO_LABEL).alias("solver_label"),
        pl.col("strategy_id")
        .map_elements(classify_strategy_tier, return_dtype=pl.String)
        .alias("hygiene_tier"),
    )
    callouts_map: dict[tuple[str, str], tuple[str, str, int, int]] = {
        ("Baseline_D0_GondalBank", "D0_Baseline_GondalBank"): ("D0 GondalBank", "right", -8, 1),
        ("G0_Filtered_HVG", "D2_DWLS"): ("G0 Filtered HVG", "left", 7, -7),
        ("G1_Unsup_Markers", "D3_CIBERSORT_NuSVR"): ("G1 Unsup. Markers", "left", 14, 11),
        ("G2b_LOCO_Milo_GMRF", "D3_CIBERSORT_NuSVR"): ("Alt 1 (LOCO+NuSVR)", "right", -8, -8),
        ("G2a_Cohort_Milo_GMRF", "D2_DWLS"): ("Alt 2 (Cohort+DWLS)", "right", -8, -7),
        ("G2d_State_Milo_GMRF", "D2_DWLS"): ("Alt 3 (State+GMRF)", "left", 8, 7),
        ("G2d_State_Milo", "D2_DWLS"): ("Alt 4 (State+DWLS)", "left", 8, -8),
    }
    callout_rows: list[dict[str, object]] = []
    for row in v2_df.iter_rows(named=True):
        pair = (str(row["strategy_id"]), str(row["deconv_tool"]))
        if pair in callouts_map:
            lbl_c, al_c, dx_c, dy_c = callouts_map[pair]
            callout_rows.append(
                {
                    "v2_pseudobulk_clr_r": float(row["v2_pseudobulk_clr_r"]),
                    "global_pearson_r": float(row["global_pearson_r"]),
                    "callout": lbl_c,
                    "align": al_c,
                    "dx": dx_c,
                    "dy": dy_c,
                }
            )

    solver_dom = ["D0: GondalBank", *SOLVER_5_ORDER]
    solver_cols = [OKABE_GREY, OKABE_SKY_BLUE, OKABE_BLUE, OKABE_BLUISH_GREEN, OKABE_ORANGE, OKABE_VERMILION]
    v2_pts = (
        alt.Chart(_alt_data(values=v2_df.to_dicts()))
        .mark_point(filled=True, size=62, opacity=0.9, stroke="#FFFFFF", strokeWidth=0.6)
        .encode(
            x=alt.X(
                "v2_pseudobulk_clr_r:Q",
                title="SC Pseudobulk Ground-Truth CLR Recovery (r, N = 613 Samples)",
                scale=alt.Scale(domain=[0.20, 0.60]),
            ),
            y=alt.Y(
                "global_pearson_r:Q",
                title="External Bulk Cross-Modality r (N = 355)",
                scale=alt.Scale(domain=[0.08, 0.70]),
            ),
            color=alt.Color(
                "solver_label:N",
                title="Solver",
                scale=alt.Scale(domain=solver_dom, range=solver_cols),
                legend=alt.Legend(orient="top", columns=3),
            ),
            shape=alt.Shape(
                "hygiene_tier:N",
                title="Tier",
                scale=alt.Scale(
                    domain=["Unsupervised", "LOCO Cross-Atlas", "Cohort-Shared", "State-Adaptive"],
                    range=["circle", "triangle-up", "square", "diamond"],
                ),
                legend=alt.Legend(orient="right"),
            ),
        )
    )
    # Render callouts per offset group so dx/dy stay crisp
    callout_layers: list[Any] = [v2_pts]
    for c_item in callout_rows:
        callout_layers.append(
            alt.Chart(_alt_data(values=[c_item]))
            .mark_text(
                align=cast(Any, str(c_item["align"])),
                dx=int(cast(int, c_item["dx"])),
                dy=int(cast(int, c_item["dy"])),
                fontSize=8.3,
                fontWeight="bold",
                color="#0F172A",
            )
            .encode(
                x=alt.X("v2_pseudobulk_clr_r:Q"),
                y=alt.Y("global_pearson_r:Q"),
                text=alt.Text("callout:N"),
            )
        )
    panel_c = alt.layer(*callout_layers).properties(
        width=285,
        height=215,
        title=alt.TitleParams(
            text="(c) Validation V2: Pseudobulk CLR Recovery vs. Bulk Concordance",
            subtitle="Native cohort operators lift pseudobulk CLR recovery over D0 (up to 0.50 vs. 0.35)",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (d): Validation V4 — Out-of-Fold (LOTO) External Bulk Clinical Response ROC-AUC (N = 965)
    folds_v4_df = folds_range_df.with_columns(
        pl.lit(0.42).alias("base_x"),
        pl.lit(0.97).alias("txt_x"),
        pl.col("mean_auc").map_elements(lambda v: f"Mean: {v:.3f}", return_dtype=pl.String).alias("auc_summary"),
    )
    v4_rand = (
        alt.Chart(_alt_data(values=[{"x0": 0.50}]))
        .mark_rule(stroke="#94A3B8", strokeDash=[3, 3], strokeWidth=1.1)
        .encode(x="x0:Q")
    )
    v4_bars = (
        alt.Chart(_alt_data(values=folds_v4_df.to_dicts()))
        .mark_bar(height=13, opacity=0.88)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order, title=""),
            x=alt.X(
                "base_x:Q",
                title="Out-of-Trial LOTO Clinical Response ROC-AUC",
                scale=alt.Scale(domain=[0.42, 1.16]),
                axis=alt.Axis(values=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
            ),
            x2=alt.X2("mean_auc:Q"),
            color=alt.Color(
                "hygiene_tier:N",
                scale=alt.Scale(
                    domain=["Unsupervised", "LOCO Cross-Atlas", "Cohort-Shared", "State-Adaptive"],
                    range=[OKABE_GREY, OKABE_BLUISH_GREEN, OKABE_BLUE, OKABE_ORANGE],
                ),
                legend=None,
            ),
        )
    )
    v4_trial_dots = (
        alt.Chart(_alt_data(values=folds_sub_df.to_dicts()))
        .mark_circle(color="#0F172A", opacity=0.45)
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order, axis=alt.Axis(labels=False)),
            x=alt.X("loto_test_roc_auc:Q", scale=alt.Scale(domain=[0.42, 1.16])),
            size=alt.Size("n_patients:Q", title="Trial N", scale=alt.Scale(range=[18, 84]), legend=None),
        )
    )
    v4_txt = (
        alt.Chart(_alt_data(values=folds_v4_df.to_dicts()))
        .mark_text(align="left", fontSize=8.3, fontWeight="bold", color="#0F172A")
        .encode(
            y=alt.Y("paradigm_label:N", sort=flagship_order),
            x=alt.X("txt_x:Q"),
            text=alt.Text("auc_summary:N"),
        )
    )
    panel_d = alt.layer(v4_rand, v4_bars, v4_trial_dots, v4_txt).properties(
        width=295,
        height=215,
        title=alt.TitleParams(
            text="(d) Validation V4: Out-of-Trial Response ROC-AUC",
            subtitle="Bars: mean LOTO AUC | Dots: 11 held-out trials (by N)",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    row_top = alt.hconcat(panel_a, panel_b, spacing=24).resolve_scale(color="independent", shape="independent")
    row_bot = alt.hconcat(panel_c, panel_d, spacing=24).resolve_scale(color="independent", shape="independent")
    return style_publication_chart(
        alt.vconcat(row_top, row_bot, spacing=24).resolve_scale(color="independent", shape="independent")
    )


# ------------------------------------------------------------------------------
# Supplementary Figure 3: Head-to-Head 8-Paradigm Scatters, Agreement & Poles
# ------------------------------------------------------------------------------


def _build_normalized_355_df(states_df: pl.DataFrame) -> tuple[pl.DataFrame, float, float, int]:
    """Compute within-cohort z-normalized 355-state scatter DataFrame and exact summary correlations."""
    parts: list[pl.DataFrame] = []
    pos_cnt = 0
    for cid in SC_COHORT_ORDER:
        sub = states_df.filter((pl.col("reference_id") == cid) & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta"))
        m_arr = sub["milo_mean_logfc"].to_numpy()
        b_arr = sub["bulk_deconv_beta"].to_numpy()
        m_sd = max(float(np.std(m_arr, ddof=1)), 1e-8)
        b_sd = max(float(np.std(b_arr, ddof=1)), 1e-8)
        mz = (m_arr - float(np.mean(m_arr))) / m_sd
        bz = (b_arr - float(np.mean(b_arr))) / b_sd
        r_c = float(cast(Any, pearsonr(mz, bz))[0])
        if r_c > 0.0:
            pos_cnt += 1
        parts.append(
            pl.DataFrame(
                {
                    "cohort_short": [SC_COHORT_SHORT[cid]] * len(mz),
                    "milo_z": mz.tolist(),
                    "bulk_z": bz.tolist(),
                }
            )
        )
    full_df = pl.concat(parts)
    mz_all = full_df["milo_z"].to_numpy()
    bz_all = full_df["bulk_z"].to_numpy()
    r_glob = float(cast(Any, pearsonr(mz_all, bz_all))[0])
    rho_glob = float(cast(Any, spearmanr(mz_all, bz_all))[0])
    return full_df, r_glob, rho_glob, pos_cnt


def _make_mini_scatter(
    spec: ParadigmSpec,
    states_df: pl.DataFrame,
    show_x_title: bool,
    show_y_title: bool,
    show_legend: bool,
) -> Any:
    """Create one locked-scale [-3.8, +3.8] cross-modality scatter panel for the 2x4 comparison grid."""
    norm_df, r_glob, rho_glob, pos_cnt = _build_normalized_355_df(states_df)
    mz = norm_df["milo_z"].to_numpy()
    bz = norm_df["bulk_z"].to_numpy()
    slope, intercept = np.polyfit(mz, bz, 1)
    xs = np.linspace(-3.4, 3.4, 25)
    ys = slope * xs + intercept

    cohort_shorts = [SC_COHORT_SHORT[c] for c in SC_COHORT_ORDER]
    z_rules_x = alt.Chart(_alt_data(values=[{"v": 0.0}])).mark_rule(stroke="#CBD5E1", strokeDash=[3, 3]).encode(x="v:Q")
    z_rules_y = alt.Chart(_alt_data(values=[{"v": 0.0}])).mark_rule(stroke="#CBD5E1", strokeDash=[3, 3]).encode(y="v:Q")

    pts = (
        alt.Chart(_alt_data(values=norm_df.to_dicts()))
        .mark_circle(size=22, opacity=0.78)
        .encode(
            x=alt.X(
                "milo_z:Q",
                title="Single-cell milopy DA (z)" if show_x_title else "",
                scale=alt.Scale(domain=[-3.8, 3.8]),
            ),
            y=alt.Y(
                "bulk_z:Q",
                title="Bulk deconv β (z)" if show_y_title else "",
                scale=alt.Scale(domain=[-3.8, 3.8]),
            ),
            color=alt.Color(
                "cohort_short:N",
                title="Single-Cell Cohort",
                scale=alt.Scale(domain=cohort_shorts, range=SC_COHORT_COLORS),
                legend=alt.Legend(orient="top", columns=9) if show_legend else None,
            ),
        )
    )
    ols_line = (
        alt.Chart(_alt_data(values=[{"x": float(xv), "y": float(yv)} for xv, yv in zip(xs, ys)]))
        .mark_line(color="#0F172A", strokeWidth=1.6)
        .encode(x="x:Q", y="y:Q")
    )
    stat_badge = (
        alt.Chart(
            _alt_data(
                values=[
                    {
                        "bx": -3.55,
                        "by": 3.35,
                        "btxt": f"r = {r_glob:+.2f} | ρ = {rho_glob:+.2f} ({pos_cnt}/9+)",
                    }
                ]
            )
        )
        .mark_text(align="left", baseline="top", fontSize=8.5, fontWeight="bold", color="#0F172A")
        .encode(x="bx:Q", y="by:Q", text="btxt:N")
    )
    return alt.layer(z_rules_x, z_rules_y, pts, ols_line, stat_badge).properties(
        width=136,
        height=124,
        title=alt.TitleParams(text=spec.scatter_title, fontSize=9.0, fontWeight="bold"),
    )


def build_supp_fig3_head_to_head(bundle: LoadedComparisonBundle) -> Any:
    """Build Supplementary Figure 3: (a) 2x4 scatter grid, (b) 9x9 Spearman agreement, (c) Meta-42 pole stability."""
    scatter_specs = [s for s in FLAGSHIP_PARADIGMS if s.include_in_scatter8]
    row1_charts = [
        _make_mini_scatter(
            scatter_specs[i],
            bundle.paradigm_states[scatter_specs[i].key],
            show_x_title=False,
            show_y_title=(i == 0),
            show_legend=(i == 0),
        )
        for i in range(4)
    ]
    row2_charts = [
        _make_mini_scatter(
            scatter_specs[i],
            bundle.paradigm_states[scatter_specs[i].key],
            show_x_title=True,
            show_y_title=(i == 4),
            show_legend=False,
        )
        for i in range(4, 8)
    ]
    panel_a = alt.vconcat(
        alt.hconcat(*row1_charts, spacing=18),
        alt.hconcat(*row2_charts, spacing=18),
        spacing=12,
    )

    # Extract Consensus-Meta42 vectors across Single-Cell Milopy + the 8 scatter paradigms
    meta_ids = [f"M{i:02d}" for i in range(1, 43)]
    ref_st_0 = bundle.paradigm_states[scatter_specs[0].key].filter(
        (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
    )
    sc_map = dict(zip(ref_st_0["cell_state"].to_list(), ref_st_0["milo_mean_logfc"].to_list()))
    vec_dict: dict[str, np.ndarray[Any, Any]] = {
        "SC Milopy": np.array([float(sc_map[m]) for m in meta_ids], dtype=np.float64)
    }
    method_order = ["SC Milopy", *[s.compact_label for s in scatter_specs]]
    for spec in scatter_specs:
        m_df = bundle.paradigm_states[spec.key].filter(
            (pl.col("reference_id") == "Consensus-Meta42") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        )
        b_map = dict(zip(m_df["cell_state"].to_list(), m_df["bulk_deconv_beta"].to_list()))
        vec_dict[spec.compact_label] = np.array([float(b_map[m]) for m in meta_ids], dtype=np.float64)

    # Panel (b): 9x9 Pairwise Spearman Rank-Agreement Matrix across 42 Consensus-Meta42 States
    agree_rows: list[dict[str, object]] = []
    for m_a in method_order:
        for m_b in method_order:
            rho_ab = float(cast(Any, spearmanr(vec_dict[m_a], vec_dict[m_b]))[0])
            agree_rows.append(
                {
                    "method_row": m_a,
                    "method_col": m_b,
                    "rho": rho_ab,
                    "rho_txt": f"{rho_ab:.2f}",
                }
            )
    agree_df = pl.DataFrame(agree_rows)
    pb_rect = (
        alt.Chart(_alt_data(values=agree_df.to_dicts()))
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.0)
        .encode(
            x=alt.X("method_col:N", sort=method_order, title="", axis=alt.Axis(labelAngle=-35)),
            y=alt.Y("method_row:N", sort=method_order, title=""),
            color=alt.Color(
                "rho:Q",
                title="Spearman ρ",
                scale=alt.Scale(domain=[0.15, 0.55, 1.0], range=["#F8FAFC", OKABE_SKY_BLUE, OKABE_BLUE]),
                legend=alt.Legend(orient="right", gradientLength=115),
            ),
        )
    )
    pb_dark = (
        alt.Chart(_alt_data(values=agree_df.filter(pl.col("rho") < 0.62).to_dicts()))
        .mark_text(fontSize=7.8, fontWeight="bold", color="#0F172A")
        .encode(x=alt.X("method_col:N", sort=method_order), y=alt.Y("method_row:N", sort=method_order), text="rho_txt:N")
    )
    pb_light = (
        alt.Chart(_alt_data(values=agree_df.filter(pl.col("rho") >= 0.62).to_dicts()))
        .mark_text(fontSize=7.8, fontWeight="bold", color="#FFFFFF")
        .encode(x=alt.X("method_col:N", sort=method_order), y=alt.Y("method_row:N", sort=method_order), text="rho_txt:N")
    )
    panel_b = alt.layer(pb_rect, pb_dark, pb_light).properties(
        width=235,
        height=235,
        title=alt.TitleParams(
            text="(b) 9×9 Inter-Paradigm Rank Consensus",
            subtitle="Pairwise Spearman ρ on Consensus-Meta42 (K = 42)",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (c): Standardized Effect-Size & Polarity Stability across 14 Canonical Responder/Non-Responder Poles
    pole_labels_order = [strip_numeric_state_prefix(m) for m in CANONICAL_POLE_STATES]
    pole_rows: list[dict[str, object]] = []
    for m_name in method_order:
        v_raw = vec_dict[m_name]
        v_z = (v_raw - float(np.mean(v_raw))) / max(float(np.std(v_raw, ddof=1)), 1e-8)
        z_by_id = dict(zip(meta_ids, v_z.tolist()))
        for st_code in CANONICAL_POLE_STATES:
            zv = float(np.clip(z_by_id[st_code], -2.6, 2.6))
            pole_rows.append(
                {
                    "state_clean": strip_numeric_state_prefix(st_code),
                    "method_col": m_name,
                    "effect_z": zv,
                    "z_txt": f"{zv:+.1f}",
                }
            )
    pole_df = pl.DataFrame(pole_rows)
    pc_rect = (
        alt.Chart(_alt_data(values=pole_df.to_dicts()))
        .mark_rect(stroke="#FFFFFF", strokeWidth=0.9)
        .encode(
            x=alt.X("method_col:N", sort=method_order, title="", axis=alt.Axis(labelAngle=-35)),
            y=alt.Y("state_clean:N", sort=pole_labels_order, title="Canonical Consensus-Meta42 Poles"),
            color=alt.Color(
                "effect_z:Q",
                title="Effect z",
                scale=alt.Scale(
                    domain=[-2.4, 0.0, 2.4],
                    range=[OKABE_VERMILION, "#F8FAFC", OKABE_BLUE],
                ),
                legend=alt.Legend(orient="right", gradientLength=115),
            ),
        )
    )
    pc_dark = (
        alt.Chart(_alt_data(values=pole_df.filter(pl.col("effect_z").abs() < 1.35).to_dicts()))
        .mark_text(fontSize=7.5, fontWeight="bold", color="#0F172A")
        .encode(x=alt.X("method_col:N", sort=method_order), y=alt.Y("state_clean:N", sort=pole_labels_order), text="z_txt:N")
    )
    pc_light = (
        alt.Chart(_alt_data(values=pole_df.filter(pl.col("effect_z").abs() >= 1.35).to_dicts()))
        .mark_text(fontSize=7.5, fontWeight="bold", color="#FFFFFF")
        .encode(x=alt.X("method_col:N", sort=method_order), y=alt.Y("state_clean:N", sort=pole_labels_order), text="z_txt:N")
    )
    panel_c = alt.layer(pc_rect, pc_dark, pc_light).properties(
        width=310,
        height=235,
        title=alt.TitleParams(
            text="(c) Canonical Biomarker Polarity Stability",
            subtitle="TLS B, Plasma & CXCL13+ T vs. myCAF, Hypoxia & MDSC",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    bot_row = alt.hconcat(panel_b, panel_c, spacing=24).resolve_scale(color="independent")
    return style_publication_chart(alt.vconcat(panel_a, bot_row, spacing=22).resolve_scale(color="independent"))


# ------------------------------------------------------------------------------
# Supplementary Figure 4: Gene-Space Mechanisms, FACS Rescue & Ablation Waterfall
# ------------------------------------------------------------------------------


def build_supp_fig4_mechanisms(bundle: LoadedComparisonBundle) -> Any:
    """Build Supplementary Figure 4: (a) gene-space concordance, (b) FACS melanoma rescue, (c) ablation waterfall."""
    # Panel (a): Gene-Space Concordance (U = 4,529 union genes)
    g_df = (
        bundle.gene_space_df.lazy()
        .filter(pl.col("is_measured_sc") & (pl.col("n_bulk_trials_measured") >= 3))
        .with_columns(
            pl.when(pl.col("is_cell_cycle_or_dissoc") | (~pl.col("is_clean_gene")))
            .then(pl.lit("Filtered Confounder (Cycle / Ribo)"))
            .when((pl.col("mean_sc_milo_r") >= 0.12) & (pl.col("bulk_meta_z") >= 1.0))
            .then(pl.lit("Conserved Responder (Immune / TLS)"))
            .when((pl.col("mean_sc_milo_r") <= -0.10) & (pl.col("bulk_meta_z") <= -0.8))
            .then(pl.lit("Conserved Non-Responder (CAF)"))
            .otherwise(pl.lit("Background / Context HVG"))
            .alias("gene_class")
        )
        .collect()
    )

    highlight_genes: dict[str, tuple[str, int, int]] = {
        "CXCL13": ("left", 6, -6),
        "CXCL9": ("right", -6, -6),
        "MS4A1": ("left", 6, 3),
        "JCHAIN": ("right", -6, 6),
        "CR2": ("left", 6, 0),
        "FAP": ("left", 7, 8),
        "COL5A2": ("right", -8, 9),
        "MKI67": ("right", -6, -9),
        "TOP2A": ("left", 6, 9),
        "UBE2C": ("right", -8, -8),
        "FOS": ("left", 7, -5),
    }
    g_Sub_bg = g_df.filter(pl.col("gene_class") == "Background / Context HVG")
    g_Sub_fg = g_df.filter(pl.col("gene_class") != "Background / Context HVG")

    gene_dom = [
        "Conserved Responder (Immune / TLS)",
        "Conserved Non-Responder (CAF)",
        "Filtered Confounder (Cycle / Ribo)",
        "Background / Context HVG",
    ]
    gene_cols = [OKABE_BLUE, OKABE_BLUISH_GREEN, OKABE_VERMILION, "#CBD5E1"]

    pa_x0 = alt.Chart(_alt_data(values=[{"x": 0.0}])).mark_rule(stroke="#94A3B8", strokeDash=[3, 3]).encode(x="x:Q")
    pa_y0 = alt.Chart(_alt_data(values=[{"y": 0.0}])).mark_rule(stroke="#94A3B8", strokeDash=[3, 3]).encode(y="y:Q")
    pa_bg = (
        alt.Chart(_alt_data(values=g_Sub_bg.to_dicts()))
        .mark_circle(size=14, opacity=0.32, color="#CBD5E1")
        .encode(
            x=alt.X(
                "mean_sc_milo_r:Q",
                title="Pan-Single-Cell Milopy Specificity Alignment (Mean r_SC)",
                scale=alt.Scale(domain=[-0.45, 0.48]),
                axis=alt.Axis(values=[-0.4, -0.2, 0.0, 0.2, 0.4]),
            ),
            y=alt.Y(
                "bulk_meta_z:Q",
                title="11-Trial External Bulk Responder Meta-Z (z_bulk)",
                scale=alt.Scale(domain=[-4.4, 5.4]),
            ),
        )
    )
    pa_fg = (
        alt.Chart(_alt_data(values=g_Sub_fg.to_dicts()))
        .mark_circle(size=32, opacity=0.82)
        .encode(
            x=alt.X("mean_sc_milo_r:Q"),
            y=alt.Y("bulk_meta_z:Q"),
            color=alt.Color(
                "gene_class:N",
                title="Gene-Space Concordance Class",
                scale=alt.Scale(domain=gene_dom, range=gene_cols),
                legend=alt.Legend(orient="top", columns=2),
            ),
        )
    )
    gene_layers: list[Any] = [pa_x0, pa_y0, pa_bg, pa_fg]
    for row in g_df.filter(pl.col("gene").is_in(list(highlight_genes.keys()))).iter_rows(named=True):
        gsym = str(row["gene"])
        al_g, dx_g, dy_g = highlight_genes[gsym]
        gene_layers.append(
            alt.Chart(_alt_data(values=[row]))
            .mark_text(align=cast(Any, al_g), dx=dx_g, dy=dy_g, fontSize=8.2, fontWeight="bold", color="#0F172A")
            .encode(x="mean_sc_milo_r:Q", y="bulk_meta_z:Q", text="gene:N")
        )
    panel_a = alt.layer(*gene_layers).properties(
        width=295,
        height=225,
        title=alt.TitleParams(
            text="(a) Cross-Platform Gene-Space Concordance (U = 4,529 Genes)",
            subtitle="Separates conserved TLS/T/CAF drivers from discordant proliferation confounders",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (b): Mechanistic Rescue in CD45+ FACS-Sorted Melanoma (Sade-Feldman GSE120575, K = 38 states)
    sf_meta_map: dict[str, str] = {}
    for row in bundle.cluster_meta_df.filter(pl.col("cohort") == "GSE120575").iter_rows(named=True):
        m_code = str(row["meta_state"])
        sf_meta_map[str(row["cell_state"])] = strip_numeric_state_prefix(m_code)

    sf_specs: tuple[tuple[str, str, str], ...] = (
        ("d0_base", "D0 GondalBank", OKABE_VERMILION),
        ("alt1_nusvr", "Alt 1: LOCO + NuSVR", OKABE_BLUISH_GREEN),
        ("alt2_cohort", "Alt 2: Cohort Milo + DWLS", OKABE_BLUE),
    )
    sf_pts_rows: list[dict[str, object]] = []
    sf_line_rows: list[dict[str, object]] = []
    sf_z_by_key: dict[str, dict[str, tuple[float, float]]] = {}
    sf_dom: list[str] = []
    sf_cols: list[str] = []
    sf_rs: dict[str, float] = {}

    for p_key, p_base_lbl, p_col in sf_specs:
        sub = bundle.paradigm_states[p_key].filter(
            (pl.col("reference_id") == "GSE120575") & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        )
        sts = sub["cell_state"].cast(pl.String).to_list()
        mv = sub["milo_mean_logfc"].to_numpy()
        bv = sub["bulk_deconv_beta"].to_numpy()
        mz = (mv - float(np.mean(mv))) / max(float(np.std(mv, ddof=1)), 1e-8)
        bz = (bv - float(np.mean(bv))) / max(float(np.std(bv, ddof=1)), 1e-8)
        r_sf = float(cast(Any, pearsonr(mz, bz))[0])
        sf_rs[p_key] = r_sf
        p_lbl = f"{p_base_lbl} (r = {r_sf:+.2f})"
        sf_dom.append(p_lbl)
        sf_cols.append(p_col)
        sf_z_by_key[p_key] = {st: (float(m_i), float(b_i)) for st, m_i, b_i in zip(sts, mz, bz)}
        for st, m_i, b_i in zip(sts, mz, bz):
            sf_pts_rows.append(
                {
                    "method_lbl": p_lbl,
                    "state_bio": sf_meta_map.get(st, "Immune Subset"),
                    "milo_z": float(m_i),
                    "bulk_z": float(b_i),
                }
            )
        slp, icp = np.polyfit(mz, bz, 1)
        for xv in np.linspace(-2.6, 3.2, 20):
            sf_line_rows.append({"method_lbl": p_lbl, "x": float(xv), "y": float(slp * xv + icp)})

    # Top rescued responder states (deduplicated by biological phenotype label to guarantee zero overlap)
    raw_candidates: list[dict[str, object]] = []
    for st, (mz0, bz0) in sf_z_by_key["d0_base"].items():
        _, bz2 = sf_z_by_key["alt2_cohort"][st]
        if mz0 > 0.65 and (bz2 - bz0) > 0.75:
            raw_candidates.append(
                {
                    "milo_z": mz0,
                    "bz_d0": bz0,
                    "bz_alt2": bz2,
                    "bio_lbl": sf_meta_map.get(st, "TLS / Stem Lymphoid"),
                }
            )
    raw_candidates = sorted(
        raw_candidates,
        key=lambda d: float(cast(float, d["bz_alt2"])) - float(cast(float, d["bz_d0"])),
        reverse=True,
    )
    rescue_arrows: list[dict[str, object]] = []
    seen_bio: set[str] = set()
    for rc in raw_candidates:
        b_name = str(rc["bio_lbl"])
        if b_name not in seen_bio:
            seen_bio.add(b_name)
            rescue_arrows.append(rc)
        if len(rescue_arrows) >= 2:
            break

    pb_zero_y = alt.Chart(_alt_data(values=[{"y": 0.0}])).mark_rule(stroke="#94A3B8", strokeDash=[3, 3]).encode(y="y:Q")
    pb_arr = (
        alt.Chart(_alt_data(values=rescue_arrows))
        .mark_rule(stroke="#475569", strokeDash=[2, 2], strokeWidth=1.2)
        .encode(x="milo_z:Q", y="bz_d0:Q", y2="bz_alt2:Q")
    )
    pb_pts = (
        alt.Chart(_alt_data(values=sf_pts_rows))
        .mark_point(filled=True, size=44, opacity=0.82)
        .encode(
            x=alt.X(
                "milo_z:Q",
                title="Sade-Feldman Milopy DA (z, K = 38 States)",
                scale=alt.Scale(domain=[-2.9, 3.8]),
            ),
            y=alt.Y(
                "bulk_z:Q",
                title="External Bulk Deconvolution β (z)",
                scale=alt.Scale(domain=[-2.8, 3.1]),
            ),
            color=alt.Color(
                "method_lbl:N",
                title="Melanoma CD45+ Deconvolution",
                scale=alt.Scale(domain=sf_dom, range=sf_cols),
                legend=alt.Legend(orient="top", columns=1),
            ),
        )
    )
    pb_lines = (
        alt.Chart(_alt_data(values=sf_line_rows))
        .mark_line(strokeWidth=1.8)
        .encode(
            x="x:Q",
            y="y:Q",
            color=alt.Color("method_lbl:N", scale=alt.Scale(domain=sf_dom, range=sf_cols), legend=None),
        )
    )
    pb_lbls = (
        alt.Chart(_alt_data(values=rescue_arrows))
        .mark_text(align="right", dx=-6, dy=-6, fontSize=8.0, fontWeight="bold", color="#0F172A")
        .encode(x="milo_z:Q", y="bz_alt2:Q", text="bio_lbl:N")
    )
    panel_b = alt.layer(pb_zero_y, pb_arr, pb_pts, pb_lines, pb_lbls).properties(
        width=305,
        height=225,
        title=alt.TitleParams(
            text="(b) CD45+ FACS Cohort Rescue (Sade-Feldman Melanoma)",
            subtitle=f"Cohort-native operators reverse T/B penalty ({sf_rs['d0_base']:+.2f} → {sf_rs['alt2_cohort']:+.2f})",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    # Panel (c): Cumulative Architectural Progression & Hygiene Spectrum across All 9 Milestones
    waterfall_specs: tuple[tuple[str, str, HygieneTier], ...] = (
        ("g1_unsup", "1. Pure Bulk: Unsup. Native Markers (G1 + InstaPrism)", "Unsupervised"),
        ("d0_base", "2. Pure Bulk: Baseline Option C (D0 GondalBank + GMRF)", "Unsupervised"),
        ("alt1_nusvr", "3. Pure Bulk Alt 1: LOCO Milo + NuSVR + GMRF", "LOCO Cross-Atlas"),
        ("alt1b_ridge", "4. Pure Bulk Alt 1b: LOCO Milo + Ridge + GMRF", "LOCO Cross-Atlas"),
        ("alt2_cohort", "5. Pure Bulk Alt 2: Cohort Milo + DWLS + GMRF", "Cohort-Shared"),
        ("alt3_state_gmrf", "6. Pure Bulk Alt 3: State Milo + GMRF + DWLS", "State-Adaptive"),
        ("alt4_state_raw", "7. Pure Bulk Alt 4: State Milo + Direct DWLS", "State-Adaptive"),
        ("opt_b_raw", "8. [Calibrated] Option B: Unprojected (55% Bulk + 45% PB)", "Calibrated"),
        ("opt_a", "9. [Calibrated] Option A: GMRF + DL Precision", "Calibrated"),
    )
    wf_rows: list[dict[str, object]] = []
    wf_order = [w[1] for w in waterfall_specs]
    for p_key, step_lbl, h_tier in waterfall_specs:
        _, r_g, _, pos_n = _build_normalized_355_df(bundle.paradigm_states[p_key])
        c_df = bundle.paradigm_conc[p_key].filter(pl.col("reference_id") == "Consensus-Meta42")
        m42_r = float(c_df["pearson_r"][0])
        wf_rows.append(
            {
                "step_label": step_lbl,
                "hygiene_tier": h_tier,
                "global_r": r_g,
                "meta42_r": m42_r,
                "label_x": max(r_g, m42_r),
                "summary_txt": f"r = {r_g:+.2f} | Meta-42 = {m42_r:+.2f} ({pos_n}/9+)",
            }
        )
    wf_df = pl.DataFrame(wf_rows)
    tier_dom = ["Unsupervised", "LOCO Cross-Atlas", "Cohort-Shared", "State-Adaptive", "Calibrated"]
    tier_cols = [tier_color(cast(HygieneTier, t)) for t in tier_dom]

    pc_bars = (
        alt.Chart(_alt_data(values=wf_df.to_dicts()))
        .mark_bar(height=14, opacity=0.9)
        .encode(
            y=alt.Y("step_label:N", sort=wf_order, title=""),
            x=alt.X(
                "global_r:Q",
                title="Cross-Modality Pearson Concordance (Bars: Unified-All9 N = 355 | Diamonds: Consensus-Meta42 K = 42)",
                scale=alt.Scale(domain=[0.0, 1.25]),
            ),
            color=alt.Color(
                "hygiene_tier:N",
                title="Methodological Regime",
                scale=alt.Scale(domain=tier_dom, range=tier_cols),
                legend=alt.Legend(orient="top", columns=5),
            ),
        )
    )
    pc_m42_rule = (
        alt.Chart(_alt_data(values=wf_df.to_dicts()))
        .mark_rule(stroke="#334155", strokeDash=[2, 2], strokeWidth=1.2)
        .encode(y=alt.Y("step_label:N", sort=wf_order), x="global_r:Q", x2="meta42_r:Q")
    )
    pc_m42_pt = (
        alt.Chart(_alt_data(values=wf_df.to_dicts()))
        .mark_point(shape="diamond", filled=True, size=65, color="#0F172A")
        .encode(y=alt.Y("step_label:N", sort=wf_order), x="meta42_r:Q")
    )
    pc_lbl = (
        alt.Chart(_alt_data(values=wf_df.to_dicts()))
        .mark_text(align="left", dx=8, fontSize=8.4, fontWeight="bold", color="#0F172A")
        .encode(
            y=alt.Y("step_label:N", sort=wf_order),
            x=alt.X("label_x:Q"),
            text=alt.Text("summary_txt:N"),
        )
    )
    panel_c = alt.layer(pc_bars, pc_m42_rule, pc_m42_pt, pc_lbl).properties(
        width=605,
        height=205,
        title=alt.TitleParams(
            text="(c) Stepwise Architectural Ablation & Methodological Hygiene Spectrum Across All Explored Approaches",
            subtitle="Isolates marginal concordance gains across pure external bulk tiers (Steps 1–7, N = 965 patients) vs. pseudobulk-calibrated anchors (Steps 8–9)",
            subtitleFontSize=8.8,
            subtitleColor=COLOR_TEXT_SECONDARY,
        ),
    )

    top_row = alt.hconcat(panel_a, panel_b, spacing=24).resolve_scale(color="independent")
    return style_publication_chart(alt.vconcat(top_row, panel_c, spacing=24).resolve_scale(color="independent"))


# ------------------------------------------------------------------------------
# Imperative Shell Entrypoint
# ------------------------------------------------------------------------------


def run(config: SuppComparisonConfig) -> Result[tuple[Path, ...], str]:
    """Load validated benchmark bundles and export all 4 supplementary comparison figures."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.article_fig_dir.mkdir(parents=True, exist_ok=True)

    bundle_res = load_all_comparison_data(config)
    if isinstance(bundle_res, Failure):
        return Failure(str(bundle_res.failure()))
    bundle = bundle_res.unwrap()

    specs_to_render: tuple[tuple[str, Any, int], ...] = (
        ("supp_fig_deconv_factorial_solver_gene_matrix", build_supp_fig1_factorial_matrix(bundle), 5000),
        ("supp_fig_deconv_four_axis_validation_suite", build_supp_fig2_validation_suite(bundle), 5000),
        ("supp_fig_deconv_approaches_head_to_head_comparison", build_supp_fig3_head_to_head(bundle), 6000),
        ("supp_fig_deconv_gene_selection_mechanisms", build_supp_fig4_mechanisms(bundle), 2500),
    )

    written: list[Path] = []
    for stem, chart_obj, max_pts in specs_to_render:
        exp_res = finalize_and_verify_svg(
            chart_obj,
            out_dir=config.out_dir,
            article_dir=config.article_fig_dir,
            stem=stem,
            width_mm=config.width_mm,
            max_vector_points=max_pts,
        )
        if isinstance(exp_res, Failure):
            return Failure(str(exp_res.failure()))
        svg_p, png_p = exp_res.unwrap()
        written.extend((svg_p, png_p))

    return Success(tuple(written))


def main() -> None:
    """CLI entrypoint for generating the 4 deconvolution approaches supplementary figures."""
    match run(SuppComparisonConfig()):
        case Success(paths):
            print("\n".join(f"wrote {p}" for p in paths))
        case Failure(err):
            print(f"4l failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
