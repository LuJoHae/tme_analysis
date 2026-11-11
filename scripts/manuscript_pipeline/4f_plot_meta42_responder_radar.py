#!/usr/bin/env python3
"""Publication-Grade Radar Immunogram for Consensus-Meta42 Responder vs. Non-Responder Fractions.

Constructs a 4-panel (180.0 mm full-width) Nature Methods wireframe vector SVG and 300 DPI PNG
visualizing deconvolved Consensus-Meta42 cell-state fractions across Responder and Non-Responder
groups in the Pan-iAtlas clinical trial compendium (N = 1,015 across 9 solid-tumor trials):

1. Panel (a): Consensus-Meta42 Classifier-Predicted Responder (y_hat = 1) vs. Non-Responder (y_hat = 0)
   12-spoke radar immunogram (5x5 repeated stratified CV L2-penalized logistic regression, within-trial
   standardized compositional CLR z-score +/- 95% CI around the z = 0 neutral reference ring).
2. Panel (b): Observed Clinical RECIST 1.1 Responder (CR/PR, n = 319) vs. Non-Responder (SD/PD, n = 696)
   12-spoke radar immunogram with Benjamini-Hochberg FDR significance annotations per spoke.
3. Panel (c): Zero-Leakage 2x2 Out-of-Fold Stromal-Gating Regime Radar on raw unprojected InstaPrism
   fractions (High-R / Permissive vs. High-R / Excluded vs. Low-R / Excluded), illustrating how
   myofibroblastic CAF (M17) co-expansion gates lymphoid and humoral effector infiltration.
4. Panel (d): Direct 12-Marker Sub-Simplex Fraction Share (theta_tilde_12 %, summing to 100% across
   the 12 marker states, parity = 8.33%) paired dumbbell & 95% CI readout with Meta-9 beta callouts.

Strict functional Python: frozen Pydantic models, returns Result/Maybe (zero None/Optional),
Polars-only DataFrames passed natively to Altair, Inkscape layer wrapping, and 300 DPI rsvg-convert verification.
"""

from __future__ import annotations

import argparse
from enum import StrEnum
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Final, Sequence, assert_never

import altair as alt
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import stats
from sklearn.ensemble import RandomForestClassifier  # type: ignore[import-untyped]
from sklearn.linear_model import LogisticRegression, Ridge  # type: ignore[import-untyped]
from sklearn.model_selection import RepeatedStratifiedKFold  # type: ignore[import-untyped]
from sklearn.preprocessing import StandardScaler  # type: ignore[import-untyped]
import vl_convert as vlc

# Disable row limit for Altair vector serialization
alt.data_transformers.disable_max_rows()

# ------------------------------------------------------------------------------
# Nature Methods Dimensional & Okabe-Ito Colorblind-Safe Styling Constants
# ------------------------------------------------------------------------------

FULL_WIDTH_MM: Final[float] = 180.0
MIN_FONT_PX: Final[float] = 13.5  # 10.125 pt at 96 CSS DPI (strictly >= 10 pt)
TITLE_FONT_PX: Final[float] = 14.0  # 10.5 pt at 96 CSS DPI

OKABE_BLUE: Final[str] = "#0072B2"       # Responder / Effector pole
OKABE_VERMILION: Final[str] = "#D55E00"  # Non-Responder / Barrier pole
OKABE_TEAL: Final[str] = "#009E73"       # High-R / Permissive regime
OKABE_PURPLE: Final[str] = "#CC79A7"     # High-R / Excluded gated regime
OKABE_ORANGE: Final[str] = "#E69F00"     # Secondary highlight
OKABE_SKY: Final[str] = "#56B4E9"        # Auxiliary reference
COLOR_SLATE_DARK: Final[str] = "#0F172A"
COLOR_SLATE_MID: Final[str] = "#334155"
COLOR_SLATE_MUTED: Final[str] = "#64748B"
COLOR_HAIRLINE: Final[str] = "#CBD5E1"
COLOR_GRID_LIGHT: Final[str] = "#E2E8F0"
COLOR_BG_WHITE: Final[str] = "#FFFFFF"

ALL_META42_COLS: Final[tuple[str, ...]] = tuple(f"M{i:02d}" for i in range(1, 43))

# Canonical compact spoke labels & compartment tags for all 42 Consensus-Meta42 states
META42_COMPACT_LABELS: Final[dict[str, tuple[str, str]]] = {
    "M01": ("TRBV18+ Mem T", "T/NK"),
    "M02": ("GATA4+ Stroma", "Stromal"),
    "M03": ("CD200+ Tfh T", "T/NK"),
    "M04": ("GNLY+PRF1+ Teff", "T/NK"),
    "M05": ("CD96+ Mem T", "T/NK"),
    "M06": ("CX3CR1+ NK/T", "T/NK"),
    "M07": ("ANK1+ Eryth.", "Other"),
    "M08": ("CCR4+ Th2/Treg", "T/NK"),
    "M09": ("SIGLEC8+ Gran.", "Myeloid"),
    "M10": ("FCN1+CD1A+ cDC", "Myeloid"),
    "M11": ("AQP4+ Ductal", "Tumor"),
    "M12": ("CA9+SLC2A1+ Hypox.", "Tumor"),
    "M13": ("CXCL13+PDCD1+ Tex", "CD8+ T"),
    "M14": ("CXCL9+STAT1+ M1-TAM", "Myeloid"),
    "M15": ("TCF7+IL7R+ Stem CD8", "CD8+ T"),
    "M16": ("VEGFA+KDR+ Endo.", "Endo."),
    "M17": ("ISLR+COL5A2+ myCAF", "Stromal"),
    "M18": ("S100B+TYR+ Melan.", "Tumor"),
    "M19": ("S100A8+CXCL8+ MDSC", "Myeloid"),
    "M20": ("IL7R+CD40LG+ Trm", "T/NK"),
    "M21": ("FOXP3+CTLA4+ Treg", "CD4+ Treg"),
    "M22": ("MKI67+TOP2A+ Prolif", "Tumor"),
    "M23": ("FAP+TGFB1+ Mat. CAF", "Stromal"),
    "M24": ("CXCL13+NMB+ Tfh", "T/NK"),
    "M25": ("CTLA4+TIGIT+ Tex", "CD8+ T"),
    "M26": ("VIM+ZEB1+ EMT", "Tumor"),
    "M27": ("ABCC8+ Neuroend.", "Tumor"),
    "M28": ("CXCL8+ITGAM+ MDSC", "Myeloid"),
    "M29": ("JCHAIN+MZB1+ Plasma", "Plasma"),
    "M30": ("CR2+PAX5+ TLS B", "B Cell"),
    "M31": ("CLEC4C+LILRA4+ pDC", "Myeloid"),
    "M32": ("MRC1+CD163+ M2-TAM", "Myeloid"),
    "M33": ("CCL5+GNLY+ CD8 Teff", "CD8+ T"),
    "M34": ("GZMK+NKG7+ Tem CD8", "CD8+ T"),
    "M35": ("MPO+CD79A+ Myeloid", "Myeloid"),
    "M36": ("FCGR3A+GNLY+ NK", "NK"),
    "M37": ("MS4A1+CD79A+ TLS B", "B Cell"),
    "M38": ("FCGR3B+FPR2+ TAN", "Myeloid"),
    "M39": ("MS4A2+TPSAB1+ Mast", "Myeloid"),
    "M40": ("FOXP3+CCR8+ Treg", "CD4+ Treg"),
    "M41": ("LAMP3+CCR7+ mregDC", "Myeloid"),
    "M42": ("SULT1B1+RGS2+ TAM", "Myeloid"),
}

# Curated 12-state pan-cancer immunogram ordered clockwise from 12:00:
#   0..5  (Right Hemisphere, 0°..150°):  6 Responder Effector-Humoral-APC Markers
#   6..11 (Left Hemisphere,  180°..330°): 6 Non-Responder Stromal-Suppressive-Mesenchymal Barriers
CURATED_12_MARKER_ORDER: Final[tuple[str, ...]] = (
    "M13",  # 0°   (12:00) CXCL13+PDCD1+ Tex       [CD8+ T Effector]
    "M15",  # 30°  (1:00)  TCF7+IL7R+ Stem CD8     [CD8+ T Renewal]
    "M34",  # 60°  (2:00)  GZMK+NKG7+ Tem CD8      [CD8+ T Cytotoxic]
    "M14",  # 90°  (3:00)  CXCL9+STAT1+ M1-TAM     [IFN-gamma M1 Macrophage]
    "M37",  # 120° (4:00)  MS4A1+CD79A+ TLS B      [Follicular TLS B]
    "M29",  # 150° (5:00)  JCHAIN+MZB1+ Plasma     [Humoral Plasma / Epistasis Hub]
    "M23",  # 180° (6:00)  FAP+TGFB1+ Mat. CAF     [Matrix-Remodeling CAF]
    "M17",  # 210° (7:00)  ISLR+COL5A2+ myCAF      [Myofibroblastic CAF / #1 RF Gini]
    "M19",  # 240° (8:00)  S100A8+CXCL8+ MDSC      [PMN-MDSC / ECM-Myeloid]
    "M21",  # 270° (9:00)  FOXP3+CTLA4+ Treg       [Effector Regulatory T]
    "M26",  # 300° (10:00) VIM+ZEB1+ EMT Tumor     [Mesenchymal Dedifferentiation]
    "M12",  # 330° (11:00) CA9+SLC2A1+ Hypox.      [Hypoxic Metabolic Barrier]
)


# ------------------------------------------------------------------------------
# 1. Immutable Domain Models & Enums (Zero None / Optional)
# ------------------------------------------------------------------------------


class MarkerSelectionPreset(StrEnum):
    """Preset strategies for selecting Consensus-Meta42 radar marker spokes."""

    CURATED_12_IMMUNOGRAM = "curated_12"
    TOP_14_META_POLES = "top_14"


class RadialScaleMode(StrEnum):
    """Radial coordinate scaling mode for compositional cell-state fractions."""

    WITHIN_TRIAL_CLR_Z = "within_trial_clr_z"
    REL_FOLD_PCT = "rel_fold_pct"


class ClassifierMode(StrEnum):
    """Supervised classifier family used to classify samples on Consensus-Meta42."""

    RIDGE_L2 = "ridge_l2"
    LINEAR_REG = "linear_reg"
    RANDOM_FOREST = "random_forest"


class SpokePole(StrEnum):
    """Biological hemisphere pole of a radar marker state."""

    RESPONDER_NICHE = "Responder Effector-Humoral Niche"
    NON_RESPONDER_BARRIER = "Non-Responder Stromal-Suppressive Barrier"


class TextAnchorGroup(StrEnum):
    """Horizontal text-anchor bucket for collision-free perimeter spoke labels."""

    TOP_CENTER = "top_center"
    BOTTOM_CENTER = "bottom_center"
    RIGHT_SIDE = "right_side"
    LEFT_SIDE = "left_side"


class MarkerSpokeSpec(BaseModel):
    """Immutable specification of a single radar spoke around the clock face."""

    model_config = ConfigDict(frozen=True)

    spoke_index: int
    state_id: str
    short_label: str
    compartment: str
    pole: SpokePole
    angle_rad: float
    angle_deg: float
    sin_a: float
    cos_a: float
    anchor_group: TextAnchorGroup
    bulk_deconv_beta: float
    milo_mean_logfc: float


class ExemplarPatientSpec(BaseModel):
    """Immutable specification of a single exemplar patient for individual radar profiling."""

    model_config = ConfigDict(frozen=True)

    panel_tag: str
    sample_id: str
    short_patient_label: str
    pole_category: str
    stroke_hex: str
    marker_shape: str


EXEMPLAR_PATIENTS: Final[tuple[ExemplarPatientSpec, ...]] = (
    ExemplarPatientSpec(
        panel_tag="(a)",
        sample_id="Liu-iAtlas__Liu_Sample137",
        short_patient_label="R1: Liu_Sample137 (Melanoma CR)",
        pole_category="Responder",
        stroke_hex=OKABE_BLUE,
        marker_shape="circle",
    ),
    ExemplarPatientSpec(
        panel_tag="(b)",
        sample_id="Rosenberg-iAtlas__10239",
        short_patient_label="R2: IMvigor210 #10239 (Bladder PR)",
        pole_category="Responder",
        stroke_hex=OKABE_TEAL,
        marker_shape="diamond",
    ),
    ExemplarPatientSpec(
        panel_tag="(c)",
        sample_id="Hugo-iAtlas__Pt10",
        short_patient_label="NR1: Hugo_Pt10 (Melanoma PD)",
        pole_category="Non-Responder",
        stroke_hex=OKABE_VERMILION,
        marker_shape="square",
    ),
    ExemplarPatientSpec(
        panel_tag="(d)",
        sample_id="McDermott-iAtlas__p01492",
        short_patient_label="NR2: IMmotion150 #p01492 (ccRCC PD)",
        pole_category="Non-Responder",
        stroke_hex=OKABE_ORANGE,
        marker_shape="triangle-down",
    ),
)


class RadarPlotConfig(BaseModel):
    """Immutable configuration for Stage 4f Consensus-Meta42 radar figure generation."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    data_dir: Path = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    results_dir: Path = Path("output/manuscript_pipeline/4f_plot_meta42_responder_radar")
    article_dir: Maybe[Path] = Some(Path("article/figures/milo-deconv"))
    preset: MarkerSelectionPreset = MarkerSelectionPreset.CURATED_12_IMMUNOGRAM
    radial_mode: RadialScaleMode = RadialScaleMode.WITHIN_TRIAL_CLR_Z
    classifier_mode: ClassifierMode = ClassifierMode.RIDGE_L2
    seed: int = 42


# ------------------------------------------------------------------------------
# 2. Pure Functional Core: Marker Selection, CLR-Z Transform & OOF Classification
# ------------------------------------------------------------------------------


def _classify_text_anchor(sin_a: float, cos_a: float) -> TextAnchorGroup:
    """Purely assign a spoke angle to one of four non-overlapping text alignment buckets."""
    if abs(sin_a) <= 0.18:
        return TextAnchorGroup.TOP_CENTER if cos_a >= 0.0 else TextAnchorGroup.BOTTOM_CENTER
    return TextAnchorGroup.RIGHT_SIDE if sin_a > 0.18 else TextAnchorGroup.LEFT_SIDE


def select_radar_marker_spokes(
    states_df: pl.DataFrame,
    preset: MarkerSelectionPreset,
) -> Result[tuple[MarkerSpokeSpec, ...], str]:
    """Select and order K Consensus-Meta42 marker states clockwise around the radar perimeter."""
    meta_pan = (
        states_df.lazy()
        .filter(
            (pl.col("reference_id") == "Consensus-Meta42")
            & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
        )
        .sort("bulk_deconv_beta", descending=True)
        .collect()
    )
    if meta_pan.is_empty():
        return Failure("No Consensus-Meta42 Pan-iAtlas-Meta rows found in states_df.")

    beta_map: dict[str, float] = {
        str(r["cell_state"]): float(r["bulk_deconv_beta"])
        for r in meta_pan.iter_rows(named=True)
    }
    milo_map: dict[str, float] = {
        str(r["cell_state"]): float(r["milo_mean_logfc"])
        for r in meta_pan.iter_rows(named=True)
    }

    match preset:
        case MarkerSelectionPreset.CURATED_12_IMMUNOGRAM:
            ordered_ids = list(CURATED_12_MARKER_ORDER)
        case MarkerSelectionPreset.TOP_14_META_POLES:
            top_r = meta_pan.head(7)["cell_state"].cast(pl.Utf8).to_list()
            # Order NR states from most negative (6:00) clockwise back toward 11:00
            top_nr = (
                meta_pan.tail(7)
                .sort("bulk_deconv_beta", descending=False)["cell_state"]
                .cast(pl.Utf8)
                .to_list()
            )
            ordered_ids = top_r + top_nr
        case _ as unreachable:
            assert_never(unreachable)

    missing = [s for s in ordered_ids if s not in beta_map]
    if missing:
        return Failure(f"Selected marker states missing from Pan-iAtlas-Meta table: {missing}")

    n_spokes = len(ordered_ids)
    half_k = n_spokes // 2
    spokes: list[MarkerSpokeSpec] = []
    for idx, st_id in enumerate(ordered_ids):
        angle_rad = 2.0 * math.pi * float(idx) / float(n_spokes)
        angle_deg = 360.0 * float(idx) / float(n_spokes)
        sin_a = math.sin(angle_rad)
        cos_a = math.cos(angle_rad)
        short_lbl, comp = META42_COMPACT_LABELS.get(st_id, (st_id, "TME"))
        pole = (
            SpokePole.RESPONDER_NICHE
            if idx < half_k
            else SpokePole.NON_RESPONDER_BARRIER
        )
        spokes.append(
            MarkerSpokeSpec(
                spoke_index=idx,
                state_id=st_id,
                short_label=short_lbl,
                compartment=comp,
                pole=pole,
                angle_rad=angle_rad,
                angle_deg=angle_deg,
                sin_a=sin_a,
                cos_a=cos_a,
                anchor_group=_classify_text_anchor(sin_a, cos_a),
                bulk_deconv_beta=beta_map[st_id],
                milo_mean_logfc=milo_map[st_id],
            )
        )
    return Success(tuple(spokes))


def compute_compositional_radar_features(
    fractions_df: pl.DataFrame,
    spokes: Sequence[MarkerSpokeSpec],
    eps: float = 1e-4,
) -> Result[pl.DataFrame, str]:
    """Compute sample-wise Aitchison CLR, within-trial standardized CLR-Z, relative fold %, and sub-simplex %."""
    required_cols = {"sample_id", "cohort", "response_binary", *ALL_META42_COLS}
    missing = required_cols - set(fractions_df.columns)
    if missing:
        return Failure(f"Fractions DataFrame missing required columns: {sorted(missing)}")
    if fractions_df.is_empty():
        return Failure("Fractions DataFrame is empty.")

    spoke_ids = [s.state_id for s in spokes]
    # 1. Sample-wise log2(theta + eps) across all 42 states and row-mean centering (Aitchison CLR)
    log_exprs = [
        (pl.col(c).clip(lower_bound=0.0) + eps).log(base=2.0).alias(f"_log_{c}")
        for c in ALL_META42_COLS
    ]
    with_logs = fractions_df.lazy().with_columns(log_exprs)
    mean_log_expr = pl.mean_horizontal([pl.col(f"_log_{c}") for c in ALL_META42_COLS]).alias("_row_mean_log")
    with_clr = with_logs.with_columns(mean_log_expr).with_columns([
        (pl.col(f"_log_{c}") - pl.col("_row_mean_log")).alias(f"clr_{c}")
        for c in ALL_META42_COLS
    ])

    # 2. Unsupervised within-cohort standardization of CLR (zero access to response_binary)
    #    plus relative fold-change % vs cohort geometric mean: 100 * (2^(clr - mean_clr_cohort) - 1)
    z_and_fold_exprs: list[pl.Expr] = []
    for c in ALL_META42_COLS:
        c_col = pl.col(f"clr_{c}")
        c_mean = c_col.mean().over("cohort")
        c_std = pl.when(c_col.std(ddof=0).over("cohort") > 1e-6).then(
            c_col.std(ddof=0).over("cohort")
        ).otherwise(1.0)
        z_and_fold_exprs.append(((c_col - c_mean) / c_std).alias(f"z_{c}"))
        z_and_fold_exprs.append(
            (100.0 * ((2.0 ** (c_col - c_mean)) - 1.0)).alias(f"fold_pct_{c}")
        )

    # 3. K-marker sub-simplex percentage share (sums to exactly 100.0% across the selected spokes per sample)
    sub_sum_expr = pl.sum_horizontal([pl.col(c).clip(lower_bound=1e-12) for c in spoke_ids]).alias("_sub_sum")
    sub_share_exprs = [
        (100.0 * pl.col(c).clip(lower_bound=1e-12) / pl.col("_sub_sum")).alias(f"sub_pct_{c}")
        for c in spoke_ids
    ]

    keep_cols = (
        ["sample_id", "cohort", "response_binary"]
        + list(ALL_META42_COLS)
        + [f"clr_{c}" for c in ALL_META42_COLS]
        + [f"z_{c}" for c in ALL_META42_COLS]
        + [f"fold_pct_{c}" for c in ALL_META42_COLS]
        + [f"sub_pct_{c}" for c in spoke_ids]
    )
    out_df = (
        with_clr.with_columns(z_and_fold_exprs)
        .with_columns(sub_sum_expr)
        .with_columns(sub_share_exprs)
        .select(keep_cols)
        .collect()
    )
    return Success(out_df)


def compute_oof_meta42_classifications(
    compositional_df: pl.DataFrame,
    mode: ClassifierMode = ClassifierMode.RIDGE_L2,
    n_splits: int = 5,
    n_repeats: int = 5,
    seed: int = 42,
) -> Result[pl.DataFrame, str]:
    """Compute honest 5x5 repeated stratified CV out-of-fold Responder classifications on Consensus-Meta42.

    Strictly fits StandardScaler and supervised classifier weights inside each training fold
    (zero test-fold outcome leakage), matching Stage 3c Task 1 (`evaluate_regularized_logistic_cv`).
    """
    if compositional_df.height < 10:
        return Failure("Too few samples to run cross-validated classification.")

    y_arr = compositional_df["response_binary"].to_numpy().astype(np.int32)
    n_pos = int(np.sum(y_arr == 1))
    n_neg = int(np.sum(y_arr == 0))
    if n_pos < 2 or n_neg < 2:
        return Failure("Both Responder and Non-Responder classes require >= 2 samples.")

    raw_mat = compositional_df.select(list(ALL_META42_COLS)).to_numpy().astype(np.float64)
    vst_mat = np.arcsin(np.sqrt(np.clip(raw_mat, 0.0, 1.0)))
    clr_z_mat = (
        compositional_df.select([f"z_{c}" for c in ALL_META42_COLS])
        .to_numpy()
        .astype(np.float64)
    )
    cohorts_arr = compositional_df["cohort"].to_numpy()
    vst_std = vst_mat.copy()
    for u_cid in np.unique(cohorts_arr):
        mask = cohorts_arr == u_cid
        if int(np.sum(mask)) >= 4:
            sub_v = vst_mat[mask, :]
            vst_std[mask, :] = (sub_v - sub_v.mean(axis=0, keepdims=True)) / np.maximum(
                sub_v.std(axis=0, keepdims=True), 1e-6
            )

    x_features = np.column_stack([vst_std, clr_z_mat])
    actual_splits = min(n_splits, n_pos, n_neg)
    rskf = RepeatedStratifiedKFold(
        n_splits=actual_splits,
        n_repeats=n_repeats,
        random_state=seed,
    )

    oof_prob_accum = np.zeros(len(y_arr), dtype=np.float64)
    oof_counts = np.zeros(len(y_arr), dtype=np.float64)

    for fold_i, (train_idx, test_idx) in enumerate(rskf.split(x_features, y_arr)):
        scaler = StandardScaler()
        x_tr = scaler.fit_transform(x_features[train_idx])
        x_te = scaler.transform(x_features[test_idx])
        y_tr = y_arr[train_idx]

        match mode:
            case ClassifierMode.RIDGE_L2:
                clf_log = LogisticRegression(
                    C=0.22,
                    l1_ratio=0.0,
                    solver="lbfgs",
                    class_weight="balanced",
                    max_iter=500,
                    random_state=seed,
                )
                clf_log.fit(x_tr, y_tr)
                fold_preds = np.asarray(clf_log.predict_proba(x_te)[:, 1], dtype=np.float64)
            case ClassifierMode.LINEAR_REG:
                reg_lin = Ridge(alpha=4.2, random_state=seed)
                w_pos = float(len(y_tr)) / max(2.0 * float(np.sum(y_tr == 1)), 1.0)
                w_neg = float(len(y_tr)) / max(2.0 * float(np.sum(y_tr == 0)), 1.0)
                sw = np.where(y_tr == 1, w_pos, w_neg)
                reg_lin.fit(x_tr, y_tr.astype(np.float64), sample_weight=sw)
                fold_preds = np.clip(
                    np.asarray(reg_lin.predict(x_te), dtype=np.float64),
                    0.01,
                    0.99,
                )
            case ClassifierMode.RANDOM_FOREST:
                rf_clf = RandomForestClassifier(
                    n_estimators=65,
                    max_depth=4,
                    min_samples_leaf=max(2, int(len(train_idx) // 18)),
                    max_features="sqrt",
                    class_weight="balanced_subsample",
                    random_state=seed + fold_i,
                    n_jobs=1,
                )
                rf_clf.fit(x_tr, y_tr)
                fold_preds = np.asarray(rf_clf.predict_proba(x_te)[:, 1], dtype=np.float64)
            case _ as unreachable:
                assert_never(unreachable)

        oof_prob_accum[test_idx] += fold_preds
        oof_counts[test_idx] += 1.0

    oof_probs = oof_prob_accum / np.maximum(oof_counts, 1.0)
    pred_bin = (oof_probs >= 0.50).astype(np.int32)
    n_pred_r = int(np.sum(pred_bin == 1))
    n_pred_nr = int(np.sum(pred_bin == 0))
    n_obs_r = int(np.sum(y_arr == 1))
    n_obs_nr = int(np.sum(y_arr == 0))

    pred_labels = [
        f"Predicted R (n={n_pred_r})" if b == 1 else f"Predicted NR (n={n_pred_nr})"
        for b in pred_bin.tolist()
    ]
    obs_labels = [
        f"Observed R (CR/PR, n={n_obs_r})" if b == 1 else f"Observed NR (SD/PD, n={n_obs_nr})"
        for b in y_arr.tolist()
    ]

    enriched_df = compositional_df.with_columns([
        pl.Series("oof_prob", oof_probs, dtype=pl.Float64),
        pl.Series("pred_binary", pred_bin, dtype=pl.Int32),
        pl.Series("pred_group", pred_labels, dtype=pl.Utf8),
        pl.Series("obs_group", obs_labels, dtype=pl.Utf8),
    ])
    return Success(enriched_df)


def _benjamini_hochberg_qvalues(p_values: Sequence[float]) -> tuple[float, ...]:
    """Purely compute Benjamini-Hochberg FDR q-values from a sequence of two-sided p-values."""
    m = len(p_values)
    if m == 0:
        return ()
    arr = np.asarray(p_values, dtype=np.float64)
    order = np.argsort(arr)
    ranked_p = arr[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    raw_q = ranked_p * float(m) / ranks
    monotone_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    clipped_q = np.clip(monotone_q, 0.0, 1.0)
    out = np.empty(m, dtype=np.float64)
    out[order] = clipped_q
    return tuple(float(v) for v in out.tolist())


def _map_value_to_radius(
    val: float,
    val_min: float,
    val_max: float,
    r_inner: float = 0.18,
    r_outer: float = 0.88,
) -> float:
    """Affinely map a scalar metric in [val_min, val_max] onto polar radius [r_inner, r_outer]."""
    span = max(val_max - val_min, 1e-9)
    frac = min(1.0, max(0.0, (val - val_min) / span))
    return r_inner + (r_outer - r_inner) * frac


def summarize_radar_polygons(
    compositional_df: pl.DataFrame,
    group_col: str,
    group_order: Sequence[str],
    spokes: Sequence[MarkerSpokeSpec],
    radial_mode: RadialScaleMode,
    val_min: float,
    val_max: float,
    r_inner: float = 0.18,
    r_outer: float = 0.88,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Summarize group means, 95% CIs, BH-FDR q-values, and closed polar-to-Cartesian polygon coordinates.

    Returns:
        (closed_polygon_df, spoke_stats_df)
        where `closed_polygon_df` has K + 1 rows per group (vertex_order 0..K, with K repeating 0)
        and `spoke_stats_df` has K rows per group for vertices, CI whiskers, and FDR callouts.
    """
    if group_col not in compositional_df.columns:
        return Failure(f"Grouping column '{group_col}' not found in DataFrame.")

    col_prefix = (
        "z_" if radial_mode == RadialScaleMode.WITHIN_TRIAL_CLR_Z else "fold_pct_"
    )

    # Compute two-sided Welch t-test p-values between the first two groups in group_order across K spokes
    p_vals: list[float] = []
    if len(group_order) >= 2:
        g0_df = compositional_df.filter(pl.col(group_col) == group_order[0])
        g1_df = compositional_df.filter(pl.col(group_col) == group_order[1])
        for sp in spokes:
            v0 = g0_df[f"{col_prefix}{sp.state_id}"].to_numpy()
            v1 = g1_df[f"{col_prefix}{sp.state_id}"].to_numpy()
            if len(v0) >= 2 and len(v1) >= 2:
                t_res = stats.ttest_ind(v0, v1, equal_var=False)
                p_raw = float(np.asarray(t_res.pvalue).item())
                p_vals.append(p_raw if math.isfinite(p_raw) else 1.0)
            else:
                p_vals.append(1.0)
    else:
        p_vals = [1.0 for _ in spokes]

    q_vals = _benjamini_hochberg_qvalues(p_vals)

    spoke_rows: list[dict[str, object]] = []
    closed_rows: list[dict[str, object]] = []

    for g_idx, g_name in enumerate(group_order):
        sub_df = compositional_df.filter(pl.col(group_col) == g_name)
        if sub_df.is_empty():
            return Failure(f"Group '{g_name}' has zero rows in '{group_col}'.")

        n_g = sub_df.height
        group_vertices: list[dict[str, object]] = []
        # Small angular jitter (+/- 0.018 rad) for CI whiskers when >2 groups overlap
        whisker_angle_offset = (float(g_idx) - 0.5 * float(len(group_order) - 1)) * 0.018

        for k_i, sp in enumerate(spokes):
            metric_vec = sub_df[f"{col_prefix}{sp.state_id}"].to_numpy().astype(np.float64)
            sub_pct_vec = sub_df[f"sub_pct_{sp.state_id}"].to_numpy().astype(np.float64)
            raw_pct_vec = 100.0 * sub_df[sp.state_id].to_numpy().astype(np.float64)

            mean_val = float(np.mean(metric_vec))
            sem_val = float(np.std(metric_vec, ddof=1) / math.sqrt(max(n_g, 1))) if n_g > 1 else 0.0
            ci_low = mean_val - 1.96 * sem_val
            ci_high = mean_val + 1.96 * sem_val

            sub_mean = float(np.mean(sub_pct_vec))
            sub_sem = float(np.std(sub_pct_vec, ddof=1) / math.sqrt(max(n_g, 1))) if n_g > 1 else 0.0

            r_mean = _map_value_to_radius(mean_val, val_min, val_max, r_inner, r_outer)
            r_low = _map_value_to_radius(ci_low, val_min, val_max, r_inner, r_outer)
            r_high = _map_value_to_radius(ci_high, val_min, val_max, r_inner, r_outer)

            ang_w = sp.angle_rad + whisker_angle_offset
            sin_w = math.sin(ang_w)
            cos_w = math.cos(ang_w)

            q_v = q_vals[k_i]
            sig_stars = (
                "***"
                if q_v < 0.001
                else ("**" if q_v < 0.01 else ("*" if q_v < 0.05 else ""))
            )

            row_dict: dict[str, object] = {
                "group": g_name,
                "group_order_idx": g_idx,
                "n_samples": n_g,
                "spoke_index": sp.spoke_index,
                "vertex_order": sp.spoke_index,
                "state_id": sp.state_id,
                "short_label": sp.short_label,
                "compartment": sp.compartment,
                "pole": sp.pole.value,
                "angle_deg": sp.angle_deg,
                "mean_val": mean_val,
                "sem_val": sem_val,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "sub_simplex_pct": sub_mean,
                "sub_simplex_ci_low": sub_mean - 1.96 * sub_sem,
                "sub_simplex_ci_high": sub_mean + 1.96 * sub_sem,
                "raw_fraction_pct": float(np.mean(raw_pct_vec)),
                "bulk_deconv_beta": sp.bulk_deconv_beta,
                "milo_mean_logfc": sp.milo_mean_logfc,
                "p_value": p_vals[k_i],
                "q_bh": q_v,
                "sig_stars": sig_stars,
                "radius": r_mean,
                "r_low": r_low,
                "r_high": r_high,
                "x": r_mean * sp.sin_a,
                "y": r_mean * sp.cos_a,
                "x_low": r_low * sin_w,
                "y_low": r_low * cos_w,
                "x_high": r_high * sin_w,
                "y_high": r_high * cos_w,
            }
            spoke_rows.append(row_dict)
            group_vertices.append(row_dict)
            closed_rows.append(row_dict)

        # Append closure vertex (vertex_order = K) identical to spoke_index = 0
        first_vertex = dict(group_vertices[0])
        first_vertex["vertex_order"] = len(spokes)
        closed_rows.append(first_vertex)

    return Success((pl.DataFrame(closed_rows), pl.DataFrame(spoke_rows)))


# ------------------------------------------------------------------------------
# 3. Declarative Altair Polar-in-Cartesian Radar Chart Builders
# ------------------------------------------------------------------------------


def _build_radar_grid_and_spokes_layers(
    spokes: Sequence[MarkerSpokeSpec],
    ring_values: Sequence[float],
    val_min: float,
    val_max: float,
    unit_suffix: str,
    show_sig_stars_map: Maybe[dict[str, str]],
    x_dom: list[float],
    y_dom: list[float],
    r_inner: float = 0.18,
    r_outer: float = 0.88,
) -> alt.LayerChart:
    """Construct concentric grid rings, highlighted z=0 reference ring, radial spokes, and perimeter labels."""
    no_axis = alt.Axis(labels=False, ticks=False, domain=False, grid=False, title="")
    n_circle_pts = 72
    angles = np.linspace(0.0, 2.0 * math.pi, n_circle_pts + 1)

    ring_rows: list[dict[str, object]] = []
    zero_ring_rows: list[dict[str, object]] = []
    tick_rows: list[dict[str, object]] = []

    for idx, r_val in enumerate(ring_values):
        rad = _map_value_to_radius(r_val, val_min, val_max, r_inner, r_outer)
        is_zero = abs(r_val) < 1e-6
        target_list = zero_ring_rows if is_zero else ring_rows
        for p_i, ang in enumerate(angles):
            target_list.append({
                "ring_val": f"{r_val:+.2f}",
                "pt_order": p_i,
                "x": float(rad * math.sin(float(ang))),
                "y": float(rad * math.cos(float(ang))),
            })
        # Show 3 crisp radial ticks (inner minimum, dashed z=0 baseline, outer maximum) along the 15° gap
        if is_zero or idx == 0 or idx == len(ring_values) - 1:
            lbl_str = "0" if is_zero else f"{r_val:+.2g}"
            tick_ang = math.radians(15.0)
            tick_rows.append({
                "x": float(rad * math.sin(tick_ang)),
                "y": float(rad * math.cos(tick_ang)),
                "tick_label": lbl_str,
                "is_zero": "zero" if is_zero else "grid",
            })

    # Hemisphere outer arc indicators (Right = Responder Effector Niche, Left = Non-Responder Barrier)
    arc_r = r_outer + 0.045
    right_arc_angles = np.linspace(math.radians(-8.0), math.radians(158.0), 40)
    left_arc_angles = np.linspace(math.radians(172.0), math.radians(338.0), 40)
    arc_rows: list[dict[str, object]] = [
        {
            "hemisphere": "Responder Pole",
            "pt_order": i,
            "x": float(arc_r * math.sin(float(a))),
            "y": float(arc_r * math.cos(float(a))),
            "hex_col": OKABE_VERMILION,
        }
        for i, a in enumerate(right_arc_angles)
    ] + [
        {
            "hemisphere": "Non-Responder Pole",
            "pt_order": i,
            "x": float(arc_r * math.sin(float(a))),
            "y": float(arc_r * math.cos(float(a))),
            "hex_col": OKABE_BLUE,
        }
        for i, a in enumerate(left_arc_angles)
    ]

    # Radial spokes from inner ring to outer ring
    spoke_line_rows: list[dict[str, object]] = []
    label_rows: list[dict[str, object]] = []
    stars_lookup = show_sig_stars_map.value_or({})

    for sp in spokes:
        spoke_line_rows.append({
            "state_id": sp.state_id,
            "x0": float(r_inner * 0.55 * sp.sin_a),
            "y0": float(r_inner * 0.55 * sp.cos_a),
            "x1": float(r_outer * sp.sin_a),
            "y1": float(r_outer * sp.cos_a),
        })
        r_lbl = r_outer + 0.115 if sp.anchor_group in (TextAnchorGroup.TOP_CENTER, TextAnchorGroup.BOTTOM_CENTER) else r_outer + 0.095
        star_s = stars_lookup.get(sp.state_id, "")
        display_lbl = f"{sp.short_label} {star_s}".strip()
        pole_col = OKABE_VERMILION if sp.pole == SpokePole.RESPONDER_NICHE else OKABE_BLUE
        label_rows.append({
            "state_id": sp.state_id,
            "display_label": display_lbl,
            "anchor_group": sp.anchor_group.value,
            "pole_hex": pole_col,
            "x": float(r_lbl * sp.sin_a),
            "y": float(r_lbl * sp.cos_a),
        })

    rings_chart = (
        alt.Chart(pl.DataFrame(ring_rows))
        .mark_line(stroke=COLOR_GRID_LIGHT, strokeWidth=0.85, opacity=0.95)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            detail="ring_val:N",
            order="pt_order:Q",
        )
    )

    zero_ring_chart = (
        alt.Chart(pl.DataFrame(zero_ring_rows))
        .mark_line(stroke=COLOR_SLATE_MID, strokeWidth=1.35, strokeDash=[4, 3], opacity=0.95)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            detail="ring_val:N",
            order="pt_order:Q",
        )
    )

    arcs_chart = (
        alt.Chart(pl.DataFrame(arc_rows))
        .mark_line(strokeWidth=2.6, opacity=0.42)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            detail="hemisphere:N",
            color=alt.Color("hex_col:N", scale=None),
            order="pt_order:Q",
        )
    )

    spokes_chart = (
        alt.Chart(pl.DataFrame(spoke_line_rows))
        .mark_rule(stroke=COLOR_HAIRLINE, strokeWidth=0.85)
        .encode(
            x=alt.X("x0:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y0:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            x2="x1:Q",
            y2="y1:Q",
        )
    )

    ticks_df = pl.DataFrame(tick_rows)
    ticks_halo = (
        alt.Chart(ticks_df)
        .mark_text(
            align="center",
            baseline="middle",
            fontSize=MIN_FONT_PX,
            fontWeight="bold",
            stroke=COLOR_BG_WHITE,
            strokeWidth=3.2,
            opacity=0.92,
        )
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="tick_label:N",
        )
    )
    ticks_chart = (
        alt.Chart(ticks_df)
        .mark_text(
            align="center",
            baseline="middle",
            fontSize=MIN_FONT_PX,
            fontWeight="bold",
            color=COLOR_SLATE_DARK,
        )
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="tick_label:N",
        )
    )

    lbl_df = pl.DataFrame(label_rows)
    lbl_top = (
        alt.Chart(lbl_df.filter(pl.col("anchor_group") == TextAnchorGroup.TOP_CENTER.value))
        .mark_text(align="center", baseline="bottom", fontSize=MIN_FONT_PX, fontWeight="bold")
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="display_label:N",
            color=alt.Color("pole_hex:N", scale=None),
        )
    )
    lbl_bot = (
        alt.Chart(lbl_df.filter(pl.col("anchor_group") == TextAnchorGroup.BOTTOM_CENTER.value))
        .mark_text(align="center", baseline="top", fontSize=MIN_FONT_PX, fontWeight="bold")
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="display_label:N",
            color=alt.Color("pole_hex:N", scale=None),
        )
    )
    lbl_right = (
        alt.Chart(lbl_df.filter(pl.col("anchor_group") == TextAnchorGroup.RIGHT_SIDE.value))
        .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="display_label:N",
            color=alt.Color("pole_hex:N", scale=None),
        )
    )
    lbl_left = (
        alt.Chart(lbl_df.filter(pl.col("anchor_group") == TextAnchorGroup.LEFT_SIDE.value))
        .mark_text(align="right", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            text="display_label:N",
            color=alt.Color("pole_hex:N", scale=None),
        )
    )

    return alt.LayerChart(
        layer=[
            rings_chart,
            spokes_chart,
            zero_ring_chart,
            arcs_chart,
            ticks_halo,
            ticks_chart,
            lbl_top,
            lbl_bot,
            lbl_right,
            lbl_left,
        ]
    )


def build_radar_panel_chart(
    closed_df: pl.DataFrame,
    spoke_stats_df: pl.DataFrame,
    spokes: Sequence[MarkerSpokeSpec],
    group_order: Sequence[str],
    group_colors: Sequence[str],
    group_shapes: Sequence[str],
    ring_values: Sequence[float],
    val_min: float,
    val_max: float,
    unit_suffix: str,
    panel_title: str,
    panel_subtitle: Sequence[str],
    show_sig_stars: bool = False,
    width_px: int = 325,
    height_px: int = 325,
) -> alt.LayerChart:
    """Build a single publication-grade radar chart panel with closed polygons and radial 95% CIs."""
    no_axis = alt.Axis(labels=False, ticks=False, domain=False, grid=False, title="")
    x_dom = [-1.42, 1.42]
    y_dom = [-1.42, 1.42]

    stars_map: Maybe[dict[str, str]] = (
        Some({
            str(r["state_id"]): str(r["sig_stars"])
            for r in spoke_stats_df.filter(pl.col("group") == group_order[0]).iter_rows(named=True)
        })
        if show_sig_stars
        else Nothing
    )

    base_grid = _build_radar_grid_and_spokes_layers(
        spokes=spokes,
        ring_values=ring_values,
        val_min=val_min,
        val_max=val_max,
        unit_suffix=unit_suffix,
        show_sig_stars_map=stars_map,
        x_dom=x_dom,
        y_dom=y_dom,
    )

    color_map = dict(zip(group_order, group_colors))
    closed_with_hex = closed_df.with_columns(
        pl.col("group").replace(color_map).alias("grp_hex")
    )
    stats_with_hex = spoke_stats_df.with_columns(
        pl.col("group").replace(color_map).alias("grp_hex")
    )

    polygons = (
        alt.Chart(closed_with_hex)
        .mark_line(fillOpacity=0.14, strokeWidth=2.3)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            detail="group:N",
            color=alt.Color("grp_hex:N", scale=None),
            fill=alt.Color("grp_hex:N", scale=None),
            order="vertex_order:Q",
        )
    )

    ci_whiskers = (
        alt.Chart(stats_with_hex)
        .mark_rule(strokeWidth=1.55, opacity=0.88)
        .encode(
            x=alt.X("x_low:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y_low:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            x2="x_high:Q",
            y2="y_high:Q",
            color=alt.Color("grp_hex:N", scale=None),
        )
    )

    vertices = (
        alt.Chart(stats_with_hex)
        .mark_point(filled=True, size=72, opacity=1.0, stroke=COLOR_SLATE_DARK, strokeWidth=0.8)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            color=alt.Color(
                "group:N",
                scale=alt.Scale(domain=list(group_order), range=list(group_colors)),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2 if len(group_order) <= 2 else 1,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=260,
                    symbolSize=85,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.8,
                ),
            ),
            shape=alt.Shape(
                "group:N",
                scale=alt.Scale(domain=list(group_order), range=list(group_shapes)),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2 if len(group_order) <= 2 else 1,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=260,
                    symbolSize=85,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.8,
                ),
            ),
        )
    )

    return alt.LayerChart(
        layer=[base_grid, polygons, ci_whiskers, vertices],
        width=width_px,
        height=height_px,
        view=alt.ViewConfig(stroke=None, fill=COLOR_BG_WHITE),
        title=alt.TitleParams(
            text=panel_title,
            subtitle=list(panel_subtitle),
            fontSize=TITLE_FONT_PX,
            subtitleFontSize=MIN_FONT_PX,
            subtitleColor=COLOR_SLATE_MID,
            fontWeight="bold",
            anchor="start",
            offset=6,
        ),
    )


def build_subsimplex_dumbbell_panel(
    pred_stats_df: pl.DataFrame,
    obs_stats_df: pl.DataFrame,
    spokes: Sequence[MarkerSpokeSpec],
    pred_group_order: Sequence[str],
    obs_group_order: Sequence[str],
    width_px: int = 325,
    height_px: int = 325,
) -> alt.LayerChart:
    """Build Panel (d): 12-Marker Sub-Simplex Share (%) paired dumbbell & 95% CI readout."""
    ordered_labels = [sp.short_label for sp in spokes]
    parity_pct = 100.0 / float(len(spokes))

    # Extract Predicted R/NR and Observed R/NR sub-simplex percentage shares
    p_r = pred_stats_df.filter(pl.col("group") == pred_group_order[0])
    p_nr = pred_stats_df.filter(pl.col("group") == pred_group_order[1])
    o_r = obs_stats_df.filter(pl.col("group") == obs_group_order[0])
    o_nr = obs_stats_df.filter(pl.col("group") == obs_group_order[1])

    p_r_map = {str(r["state_id"]): r for r in p_r.iter_rows(named=True)}
    p_nr_map = {str(r["state_id"]): r for r in p_nr.iter_rows(named=True)}
    o_r_map = {str(r["state_id"]): r for r in o_r.iter_rows(named=True)}
    o_nr_map = {str(r["state_id"]): r for r in o_nr.iter_rows(named=True)}

    pts_rows: list[dict[str, object]] = []
    link_rows: list[dict[str, object]] = []
    annot_rows: list[dict[str, object]] = []

    for sp in spokes:
        pr_row = p_r_map[sp.state_id]
        pnr_row = p_nr_map[sp.state_id]
        or_row = o_r_map[sp.state_id]
        onr_row = o_nr_map[sp.state_id]

        pr_v = float(pr_row["sub_simplex_pct"])
        pnr_v = float(pnr_row["sub_simplex_pct"])
        or_v = float(or_row["sub_simplex_pct"])
        onr_v = float(onr_row["sub_simplex_pct"])

        link_rows.append({
            "short_label": sp.short_label,
            "x_r": pr_v,
            "x_nr": pnr_v,
            "link_hex": OKABE_BLUE if pr_v >= pnr_v else OKABE_VERMILION,
        })

        pts_rows.extend([
            {
                "short_label": sp.short_label,
                "stratum": "Pred. R (Meta42 CV)",
                "sub_pct": pr_v,
                "ci_low": float(pr_row["sub_simplex_ci_low"]),
                "ci_high": float(pr_row["sub_simplex_ci_high"]),
                "pt_hex": OKABE_BLUE,
            },
            {
                "short_label": sp.short_label,
                "stratum": "Pred. NR (Meta42 CV)",
                "sub_pct": pnr_v,
                "ci_low": float(pnr_row["sub_simplex_ci_low"]),
                "ci_high": float(pnr_row["sub_simplex_ci_high"]),
                "pt_hex": OKABE_VERMILION,
            },
            {
                "short_label": sp.short_label,
                "stratum": "Obs. R (CR/PR)",
                "sub_pct": or_v,
                "ci_low": float(or_row["sub_simplex_ci_low"]),
                "ci_high": float(or_row["sub_simplex_ci_high"]),
                "pt_hex": OKABE_TEAL,
            },
            {
                "short_label": sp.short_label,
                "stratum": "Obs. NR (SD/PD)",
                "sub_pct": onr_v,
                "ci_low": float(onr_row["sub_simplex_ci_low"]),
                "ci_high": float(onr_row["sub_simplex_ci_high"]),
                "pt_hex": OKABE_ORANGE,
            },
        ])

        annot_rows.append({
            "short_label": sp.short_label,
            "x_pos": 10.22,
            "beta_txt": f"β={sp.bulk_deconv_beta:+.2f}",
            "pole_hex": OKABE_BLUE if sp.bulk_deconv_beta >= 0.0 else OKABE_VERMILION,
        })

    x_dom = [7.35, 10.85]
    stratum_order = [
        "Pred. R (Meta42 CV)",
        "Pred. NR (Meta42 CV)",
        "Obs. R (CR/PR)",
        "Obs. NR (SD/PD)",
    ]
    stratum_cols = [OKABE_BLUE, OKABE_VERMILION, OKABE_TEAL, OKABE_ORANGE]
    stratum_shapes = ["circle", "square", "triangle-up", "diamond"]

    parity_rule = (
        alt.Chart(pl.DataFrame({"x0": [parity_pct]}))
        .mark_rule(stroke=COLOR_SLATE_MID, strokeDash=[4, 3], strokeWidth=1.25)
        .encode(x=alt.X("x0:Q", scale=alt.Scale(domain=x_dom)))
    )

    links = (
        alt.Chart(pl.DataFrame(link_rows))
        .mark_rule(strokeWidth=2.0, opacity=0.55)
        .encode(
            y=alt.Y(
                "short_label:N",
                sort=ordered_labels,
                title=None,
                axis=alt.Axis(
                    labelFontSize=MIN_FONT_PX,
                    labelFontWeight="bold",
                    labelLimit=225,
                ),
            ),
            x=alt.X(
                "x_nr:Q",
                scale=alt.Scale(domain=x_dom),
                title=f"12-Marker Sub-Simplex Share (%, Parity = {parity_pct:.2f}%)",
                axis=alt.Axis(
                    values=[7.5, 8.0, round(parity_pct, 2), 9.0, 9.5, 10.0],
                    labelFontSize=MIN_FONT_PX,
                    titleFontSize=MIN_FONT_PX,
                    titleFontWeight="bold",
                    grid=True,
                    gridColor=COLOR_GRID_LIGHT,
                ),
            ),
            x2="x_r:Q",
            color=alt.Color("link_hex:N", scale=None),
        )
    )

    pts_df = pl.DataFrame(pts_rows)
    ci_bars = (
        alt.Chart(pts_df.filter(pl.col("stratum").is_in(["Pred. R (Meta42 CV)", "Pred. NR (Meta42 CV)"])))
        .mark_rule(strokeWidth=1.35, opacity=0.85)
        .encode(
            y=alt.Y("short_label:N", sort=ordered_labels),
            x=alt.X("ci_low:Q", scale=alt.Scale(domain=x_dom)),
            x2="ci_high:Q",
            color=alt.Color("pt_hex:N", scale=None),
        )
    )

    dots = (
        alt.Chart(pts_df)
        .mark_point(filled=True, size=68, opacity=0.95, stroke=COLOR_SLATE_DARK, strokeWidth=0.75)
        .encode(
            y=alt.Y("short_label:N", sort=ordered_labels),
            x=alt.X("sub_pct:Q", scale=alt.Scale(domain=x_dom)),
            color=alt.Color(
                "stratum:N",
                scale=alt.Scale(domain=stratum_order, range=stratum_cols),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=165,
                    symbolSize=75,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.75,
                ),
            ),
            shape=alt.Shape(
                "stratum:N",
                scale=alt.Scale(domain=stratum_order, range=stratum_shapes),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=165,
                    symbolSize=75,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.75,
                ),
            ),
        )
    )

    beta_labels = (
        alt.Chart(pl.DataFrame(annot_rows))
        .mark_text(align="left", baseline="middle", fontSize=MIN_FONT_PX, fontWeight="bold")
        .encode(
            y=alt.Y("short_label:N", sort=ordered_labels),
            x=alt.X("x_pos:Q", scale=alt.Scale(domain=x_dom)),
            text="beta_txt:N",
            color=alt.Color("pole_hex:N", scale=None),
        )
    )

    return alt.LayerChart(
        layer=[parity_rule, links, ci_bars, dots, beta_labels],
        width=width_px,
        height=height_px,
        title=alt.TitleParams(
            text="(d)  12-Marker Sub-Simplex Share (%) & Meta-9 β",
            subtitle=[
                "Re-closed 12-marker simplex share (sums to 100%, parity = 8.33%)",
                "across Pred. vs. Obs. R/NR groups (+/- 95% CI) & Meta-9 CLR β",
            ],
            fontSize=TITLE_FONT_PX,
            subtitleFontSize=MIN_FONT_PX,
            subtitleColor=COLOR_SLATE_MID,
            fontWeight="bold",
            anchor="start",
            offset=6,
        ),
    )


# ------------------------------------------------------------------------------
# 4. SVG Inkscape Layer Wrapping & 300 DPI Verification
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
    """Inject Inkscape namespaces, explicit background & figure layer groups, and exact dimensions in mm."""
    m = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg_text)
    vb_w = float(m.group(1)) if m else 760.0
    vb_h = float(m.group(2)) if m else 680.0
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
        f'  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">\n'
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
    article_dir: Maybe[Path],
    width_mm: float,
    height_mm: float,
) -> Result[tuple[Path, Path], str]:
    """Write validated Inkscape SVG and render 300 DPI PNG via rsvg-convert (with vl_convert fallback)."""
    try:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        svg_path.write_text(svg_text, encoding="utf-8")

        target_w_px = int(round(width_mm * 300.0 / 25.4))
        target_h_px = int(round(height_mm * 300.0 / 25.4))

        rsvg_bin = shutil.which("rsvg-convert") or "/opt/homebrew/bin/rsvg-convert"
        if Path(rsvg_bin).exists():
            subprocess.run(
                [
                    rsvg_bin,
                    "-w",
                    str(target_w_px),
                    "-h",
                    str(target_h_px),
                    "-b",
                    "white",
                    str(svg_path),
                    "-o",
                    str(png_path),
                ],
                check=True,
                capture_output=True,
            )
        else:
            png_bytes = vlc.svg_to_png(svg_text, scale=3.125)
            png_path.write_bytes(png_bytes)

        match article_dir:
            case Some(art_dir):
                art_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(svg_path, art_dir / svg_path.name)
                shutil.copy2(png_path, art_dir / png_path.name)
            case Nothing:
                pass

        return Success((svg_path, png_path))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export SVG/PNG for {svg_path}: {exc}")


# ------------------------------------------------------------------------------
# 5. End-to-End 4-Panel Radar Figure Assembly & Pipeline Execution
# ------------------------------------------------------------------------------


def build_meta42_responder_radar_figure(
    proj_fractions_df: pl.DataFrame,
    raw_fractions_df: pl.DataFrame,
    states_df: pl.DataFrame,
    oof_landscape_df: pl.DataFrame,
    config: RadarPlotConfig,
) -> Result[tuple[str, float, pl.DataFrame], str]:
    """Assemble the complete 4-panel Consensus-Meta42 Responder vs. Non-Responder radar figure."""
    match select_radar_marker_spokes(states_df, config.preset):
        case Failure(err):
            return Failure(err)
        case Success(spokes):
            pass

    match compute_compositional_radar_features(proj_fractions_df, spokes):
        case Failure(err):
            return Failure(err)
        case Success(proj_comp_df):
            pass

    match compute_oof_meta42_classifications(
        proj_comp_df,
        mode=config.classifier_mode,
        seed=config.seed,
    ):
        case Failure(err):
            return Failure(err)
        case Success(classified_df):
            pass

    match compute_compositional_radar_features(raw_fractions_df, spokes):
        case Failure(err):
            return Failure(err)
        case Success(raw_comp_df):
            pass

    # Determine group labels and radial bounds based on RadialScaleMode
    pred_r_lbl = next(
        g for g in classified_df["pred_group"].unique().to_list() if g.startswith("Predicted R ")
    )
    pred_nr_lbl = next(
        g for g in classified_df["pred_group"].unique().to_list() if g.startswith("Predicted NR ")
    )
    obs_r_lbl = next(
        g for g in classified_df["obs_group"].unique().to_list() if g.startswith("Observed R ")
    )
    obs_nr_lbl = next(
        g for g in classified_df["obs_group"].unique().to_list() if g.startswith("Observed NR ")
    )

    match config.radial_mode:
        case RadialScaleMode.WITHIN_TRIAL_CLR_Z:
            pred_bounds = (-0.65, 0.65)
            pred_rings = (-0.50, -0.25, 0.00, 0.25, 0.50)
            obs_bounds = (-0.30, 0.30)
            obs_rings = (-0.20, -0.10, 0.00, 0.10, 0.20)
            regime_bounds = (-1.25, 1.25)
            regime_rings = (-1.00, -0.50, 0.00, 0.50, 1.00)
            unit_lbl = "CLR z"
        case RadialScaleMode.REL_FOLD_PCT:
            pred_bounds = (-12.0, 12.0)
            pred_rings = (-10.0, -5.0, 0.0, 5.0, 10.0)
            obs_bounds = (-6.0, 6.0)
            obs_rings = (-4.0, -2.0, 0.0, 2.0, 4.0)
            regime_bounds = (-20.0, 20.0)
            regime_rings = (-15.0, -7.5, 0.0, 7.5, 15.0)
            unit_lbl = "Fold %"
        case _ as unreachable:
            assert_never(unreachable)

    # Panel (a): Classifier-Predicted R vs. NR on Consensus-Meta42
    match summarize_radar_polygons(
        compositional_df=classified_df,
        group_col="pred_group",
        group_order=(pred_r_lbl, pred_nr_lbl),
        spokes=spokes,
        radial_mode=config.radial_mode,
        val_min=pred_bounds[0],
        val_max=pred_bounds[1],
    ):
        case Failure(err):
            return Failure(err)
        case Success((pred_closed_df, pred_stats_df)):
            pass

    # Panel (b): Observed Clinical RECIST R (CR/PR) vs. NR (SD/PD)
    match summarize_radar_polygons(
        compositional_df=classified_df,
        group_col="obs_group",
        group_order=(obs_r_lbl, obs_nr_lbl),
        spokes=spokes,
        radial_mode=config.radial_mode,
        val_min=obs_bounds[0],
        val_max=obs_bounds[1],
    ):
        case Failure(err):
            return Failure(err)
        case Success((obs_closed_df, obs_stats_df)):
            pass

    # Panel (c): Zero-Leakage 2x2 Out-of-Fold Stromal-Gating Regimes on Raw InstaPrism Fractions
    patient_oof = (
        oof_landscape_df.lazy()
        .filter(pl.col("record_type") == "patient")
        .select(["sample_id", "regime_group"])
        .collect()
    )
    raw_joined = raw_comp_df.join(patient_oof, on="sample_id", how="inner")
    regime_map = {
        "High-R / Permissive (T3, T1)": "High-R / Permissive (ORR 51.6%, n=128)",
        "High-R / Excluded (T3, T3)": "High-R / Excluded (ORR 27.6%, n=98)",
        "Low-R / Excluded (T1, T3)": "Low-R / Excluded (ORR 25.0%, n=92)",
    }
    raw_regimes_df = (
        raw_joined.filter(pl.col("regime_group").is_in(list(regime_map.keys())))
        .with_columns(pl.col("regime_group").replace(regime_map).alias("regime_display"))
    )
    regime_order = (
        "High-R / Permissive (ORR 51.6%, n=128)",
        "High-R / Excluded (ORR 27.6%, n=98)",
        "Low-R / Excluded (ORR 25.0%, n=92)",
    )
    match summarize_radar_polygons(
        compositional_df=raw_regimes_df,
        group_col="regime_display",
        group_order=regime_order,
        spokes=spokes,
        radial_mode=config.radial_mode,
        val_min=regime_bounds[0],
        val_max=regime_bounds[1],
    ):
        case Failure(err):
            return Failure(err)
        case Success((reg_closed_df, reg_stats_df)):
            pass

    panel_a = build_radar_panel_chart(
        closed_df=pred_closed_df,
        spoke_stats_df=pred_stats_df,
        spokes=spokes,
        group_order=(pred_r_lbl, pred_nr_lbl),
        group_colors=(OKABE_BLUE, OKABE_VERMILION),
        group_shapes=("circle", "square"),
        ring_values=pred_rings,
        val_min=pred_bounds[0],
        val_max=pred_bounds[1],
        unit_suffix=unit_lbl,
        panel_title="(a)  Meta42 Classifier-Predicted R vs. NR Radar",
        panel_subtitle=(
            "5x5 OOF CV L2-logistic classification on Consensus-Meta42 (N = 1,015);",
            "within-trial CLR z-score (+/- 95% CI) around dashed z = 0 baseline ring",
        ),
        show_sig_stars=False,
    )

    panel_b = build_radar_panel_chart(
        closed_df=obs_closed_df,
        spoke_stats_df=obs_stats_df,
        spokes=spokes,
        group_order=(obs_r_lbl, obs_nr_lbl),
        group_colors=(OKABE_BLUE, OKABE_VERMILION),
        group_shapes=("circle", "square"),
        ring_values=obs_rings,
        val_min=obs_bounds[0],
        val_max=obs_bounds[1],
        unit_suffix=unit_lbl,
        panel_title="(b)  Observed Clinical RECIST R vs. NR Radar",
        panel_subtitle=(
            "Empirical CR/PR (n = 319) vs. SD/PD (n = 696) across 9 iAtlas trials;",
            "spoke asterisks denote BH-adjusted Welch t-test FDR (* <0.05, ** <0.01, *** <0.001)",
        ),
        show_sig_stars=True,
    )

    panel_c = build_radar_panel_chart(
        closed_df=reg_closed_df,
        spoke_stats_df=reg_stats_df,
        spokes=spokes,
        group_order=regime_order,
        group_colors=(OKABE_TEAL, OKABE_PURPLE, OKABE_VERMILION),
        group_shapes=("circle", "triangle-up", "square"),
        ring_values=regime_rings,
        val_min=regime_bounds[0],
        val_max=regime_bounds[1],
        unit_suffix=unit_lbl,
        panel_title="(c)  2x2 OOF Stromal-Gating Regimes (Raw InstaPrism)",
        panel_subtitle=(
            "Unblended raw InstaPrism Meta42 fractions across 5-fold OOF corners:",
            "High-R / Excluded co-expands M17 myCAF (+0.25 SD) alongside TLS B / Tex",
        ),
        show_sig_stars=False,
    )

    panel_d = build_subsimplex_dumbbell_panel(
        pred_stats_df=pred_stats_df,
        obs_stats_df=obs_stats_df,
        spokes=spokes,
        pred_group_order=(pred_r_lbl, pred_nr_lbl),
        obs_group_order=(obs_r_lbl, obs_nr_lbl),
    )

    top_row = alt.hconcat(panel_a, panel_b, spacing=135).resolve_scale(
        color="independent",
        shape="independent",
    )
    bot_row = alt.hconcat(panel_c, panel_d, spacing=135).resolve_scale(
        color="independent",
        shape="independent",
    )
    combined = (
        alt.vconcat(top_row, bot_row, spacing=22, center=False)
        .resolve_scale(color="independent", shape="independent")
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75, fill=COLOR_BG_WHITE)
        .configure(
            font="system-ui, -apple-system, sans-serif",
            background=COLOR_BG_WHITE,
            padding={"left": 115, "right": 105, "top": 12, "bottom": 12},
        )
    )

    try:
        raw_svg = str(vlc.vegalite_to_svg(combined.to_dict()))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to serialize Altair radar figure to SVG: {exc}")

    wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
        raw_svg,
        width_mm=FULL_WIDTH_MM,
        layer_id="layer-01-meta42-responder-radar",
        layer_label="01_Meta42_Responder_Radar_Panels",
    )

    summary_df = pl.concat([
        pred_stats_df.with_columns(pl.lit("panel_a_predicted_r_vs_nr").alias("panel_id")),
        obs_stats_df.with_columns(pl.lit("panel_b_observed_recist_r_vs_nr").alias("panel_id")),
        reg_stats_df.with_columns(pl.lit("panel_c_oof_2x2_stromal_regimes").alias("panel_id")),
    ])
    return Success((wrapped_svg, height_mm, summary_df))


def summarize_exemplar_patient_polygons(
    classified_df: pl.DataFrame,
    oof_landscape_df: pl.DataFrame,
    spokes: Sequence[MarkerSpokeSpec],
    exemplars: Sequence[ExemplarPatientSpec] = EXEMPLAR_PATIENTS,
    val_min: float = -2.50,
    val_max: float = 2.50,
    r_inner: float = 0.18,
    r_outer: float = 0.88,
) -> Result[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame], str]:
    """Summarize individual exemplar patient polygons and cohort R/NR reference centroid rings.

    Returns:
        (patient_closed_df, patient_spoke_df, centroid_closed_df)
    """
    patient_meta = (
        oof_landscape_df.lazy()
        .filter(pl.col("record_type") == "patient")
        .select([
            "sample_id",
            "histology",
            "response_recist",
            "regime_group",
            "os_months",
            "os_event",
        ])
        .collect()
    )
    joined_df = classified_df.join(patient_meta, on="sample_id", how="left")

    # Compute Observed CR/PR and SD/PD centroid polygons on [-2.50, +2.50] for background context
    centroid_closed_rows: list[dict[str, object]] = []
    for obs_val, ref_label in ((1, "Observed R Mean (CR/PR)"), (0, "Observed NR Mean (SD/PD)")):
        sub_obs = joined_df.filter(pl.col("response_binary") == obs_val)
        grp_pts: list[dict[str, object]] = []
        for sp in spokes:
            m_val = float(sub_obs[f"z_{sp.state_id}"].mean())  # type: ignore[arg-type]
            rad = _map_value_to_radius(m_val, val_min, val_max, r_inner, r_outer)
            pt_dict: dict[str, object] = {
                "ref_group": ref_label,
                "vertex_order": sp.spoke_index,
                "state_id": sp.state_id,
                "mean_z": m_val,
                "radius": rad,
                "x": float(rad * sp.sin_a),
                "y": float(rad * sp.cos_a),
            }
            grp_pts.append(pt_dict)
            centroid_closed_rows.append(pt_dict)
        first_pt = dict(grp_pts[0])
        first_pt["vertex_order"] = len(spokes)
        centroid_closed_rows.append(first_pt)

    patient_spoke_rows: list[dict[str, object]] = []
    patient_closed_rows: list[dict[str, object]] = []

    for ex_idx, ex in enumerate(exemplars):
        pt_df = joined_df.filter(pl.col("sample_id") == ex.sample_id)
        if pt_df.is_empty():
            return Failure(f"Exemplar patient '{ex.sample_id}' not found in classified DataFrame.")
        pt_row = pt_df.row(0, named=True)
        oof_p = float(pt_row["oof_prob"])
        recist_str = str(pt_row.get("response_recist") or ("CR/PR" if int(pt_row["response_binary"]) == 1 else "SD/PD"))
        hist_str = str(pt_row.get("histology") or str(pt_row["cohort"]))
        regime_str = str(pt_row.get("regime_group") or "")
        os_m = float(pt_row.get("os_months") or -1.0)
        os_ev = int(pt_row.get("os_event") or -1)

        r_sub_sum = sum(float(pt_row[f"sub_pct_{sp.state_id}"]) for sp in spokes[:6])
        nr_sub_sum = sum(float(pt_row[f"sub_pct_{sp.state_id}"]) for sp in spokes[6:])

        ex_vertices: list[dict[str, object]] = []
        for sp in spokes:
            z_val = float(pt_row[f"z_{sp.state_id}"])
            sub_pct = float(pt_row[f"sub_pct_{sp.state_id}"])
            raw_pct = 100.0 * float(pt_row[sp.state_id])
            rad = _map_value_to_radius(z_val, val_min, val_max, r_inner, r_outer)
            row_d: dict[str, object] = {
                "exemplar_idx": ex_idx,
                "panel_tag": ex.panel_tag,
                "sample_id": ex.sample_id,
                "short_patient_label": ex.short_patient_label,
                "pole_category": ex.pole_category,
                "stroke_hex": ex.stroke_hex,
                "marker_shape": ex.marker_shape,
                "cohort": str(pt_row["cohort"]),
                "histology": hist_str,
                "recist": recist_str,
                "oof_prob": oof_p,
                "regime_group": regime_str,
                "os_months": os_m,
                "os_event": os_ev,
                "r_pole_sub_pct": r_sub_sum,
                "nr_pole_sub_pct": nr_sub_sum,
                "spoke_index": sp.spoke_index,
                "vertex_order": sp.spoke_index,
                "state_id": sp.state_id,
                "short_label": sp.short_label,
                "compartment": sp.compartment,
                "pole": sp.pole.value,
                "angle_deg": sp.angle_deg,
                "z_val": z_val,
                "sub_simplex_pct": sub_pct,
                "raw_fraction_pct": raw_pct,
                "radius": rad,
                "x": float(rad * sp.sin_a),
                "y": float(rad * sp.cos_a),
            }
            patient_spoke_rows.append(row_d)
            ex_vertices.append(row_d)
            patient_closed_rows.append(row_d)

        close_d = dict(ex_vertices[0])
        close_d["vertex_order"] = len(spokes)
        patient_closed_rows.append(close_d)

    return Success((
        pl.DataFrame(patient_closed_rows),
        pl.DataFrame(patient_spoke_rows),
        pl.DataFrame(centroid_closed_rows),
    ))


def build_single_patient_radar_panel(
    pt_closed_df: pl.DataFrame,
    pt_spoke_df: pl.DataFrame,
    centroid_closed_df: pl.DataFrame,
    spokes: Sequence[MarkerSpokeSpec],
    exemplar: ExemplarPatientSpec,
    val_min: float = -2.50,
    val_max: float = 2.50,
    width_px: int = 325,
    height_px: int = 325,
) -> alt.LayerChart:
    """Build an individual patient radar panel overlaid on cohort R/NR reference centroid rings."""
    no_axis = alt.Axis(labels=False, ticks=False, domain=False, grid=False, title="")
    x_dom = [-1.42, 1.42]
    y_dom = [-1.42, 1.42]
    ring_values = (-2.0, -1.0, 0.0, 1.0, 2.0)

    base_grid = _build_radar_grid_and_spokes_layers(
        spokes=spokes,
        ring_values=ring_values,
        val_min=val_min,
        val_max=val_max,
        unit_suffix="CLR z",
        show_sig_stars_map=Nothing,
        x_dom=x_dom,
        y_dom=y_dom,
    )

    first_row = pt_spoke_df.row(0, named=True)
    oof_p = float(first_row["oof_prob"])
    r_share = float(first_row["r_pole_sub_pct"])
    nr_share = float(first_row["nr_pole_sub_pct"])
    os_m = float(first_row["os_months"])
    os_ev = int(first_row["os_event"])
    os_txt = (
        f"OS = {os_m:.1f} mo ({'event' if os_ev == 1 else 'censored'})"
        if os_m > 0.0
        else "OS = N/A (IMmotion150)"
    )

    is_resp = exemplar.pole_category == "Responder"
    target_centroid = "Observed R Mean (CR/PR)" if is_resp else "Observed NR Mean (SD/PD)"
    centroid_sub = centroid_closed_df.filter(pl.col("ref_group") == target_centroid)

    pt_legend_label = (
        f"Patient Profile (P_OOF(R) = {oof_p:.2f})"
    )
    combined_poly_df = pl.concat([
        centroid_sub.select([
            pl.lit(target_centroid).alias("series_label"),
            pl.col("vertex_order"),
            pl.col("x"),
            pl.col("y"),
        ]),
        pt_closed_df.select([
            pl.lit(pt_legend_label).alias("series_label"),
            pl.col("vertex_order"),
            pl.col("x"),
            pl.col("y"),
        ]),
    ])

    centroid_layer = (
        alt.Chart(centroid_sub)
        .mark_line(
            stroke=COLOR_SLATE_MID,
            strokeWidth=1.55,
            strokeDash=[5, 3],
            fill=COLOR_SLATE_MID,
            fillOpacity=0.06,
        )
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            order="vertex_order:Q",
        )
    )

    patient_poly = (
        alt.Chart(pt_closed_df)
        .mark_line(
            stroke=exemplar.stroke_hex,
            fill=exemplar.stroke_hex,
            fillOpacity=0.20,
            strokeWidth=2.5,
        )
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            order="vertex_order:Q",
        )
    )

    # Create legend-anchoring point layer with patient vertices and a reference centroid key
    centroid_vertices = centroid_sub.filter(pl.col("vertex_order") < len(spokes)).select([
        pl.lit(target_centroid).alias("series_label"),
        pl.col("x"),
        pl.col("y"),
    ])
    patient_vertices = pt_spoke_df.select([
        pl.lit(pt_legend_label).alias("series_label"),
        pl.col("x"),
        pl.col("y"),
    ])
    legend_pts_df = pl.concat([patient_vertices, centroid_vertices])
    series_order = [pt_legend_label, target_centroid]
    series_colors = [exemplar.stroke_hex, COLOR_SLATE_MID]
    series_shapes = [exemplar.marker_shape, "circle"]

    vertices_layer = (
        alt.Chart(legend_pts_df)
        .mark_point(filled=True, size=68, opacity=0.95, stroke=COLOR_SLATE_DARK, strokeWidth=0.75)
        .encode(
            x=alt.X("x:Q", scale=alt.Scale(domain=x_dom), axis=no_axis),
            y=alt.Y("y:Q", scale=alt.Scale(domain=y_dom), axis=no_axis),
            color=alt.Color(
                "series_label:N",
                scale=alt.Scale(domain=series_order, range=series_colors),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=240,
                    symbolSize=80,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.75,
                ),
            ),
            shape=alt.Shape(
                "series_label:N",
                scale=alt.Scale(domain=series_order, range=series_shapes),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="horizontal",
                    columns=2,
                    labelFontSize=MIN_FONT_PX,
                    labelLimit=240,
                    symbolSize=80,
                    symbolOpacity=1.0,
                    symbolStrokeColor=COLOR_SLATE_DARK,
                    symbolStrokeWidth=0.75,
                ),
            ),
        )
    )

    _ = combined_poly_df
    return alt.LayerChart(
        layer=[base_grid, centroid_layer, patient_poly, vertices_layer],
        width=width_px,
        height=height_px,
        view=alt.ViewConfig(stroke=None, fill=COLOR_BG_WHITE),
        title=alt.TitleParams(
            text=f"{exemplar.panel_tag}  {exemplar.short_patient_label}",
            subtitle=[
                f"OOF P(R) = {oof_p:.2f} | {os_txt}",
                f"12-Marker Share: Responder Pole = {r_share:.1f}% vs. Barrier Pole = {nr_share:.1f}%",
            ],
            fontSize=TITLE_FONT_PX,
            subtitleFontSize=MIN_FONT_PX,
            subtitleColor=COLOR_SLATE_MID,
            fontWeight="bold",
            anchor="start",
            offset=6,
        ),
    )


def build_exemplar_patients_radar_figure(
    proj_fractions_df: pl.DataFrame,
    states_df: pl.DataFrame,
    oof_landscape_df: pl.DataFrame,
    config: RadarPlotConfig,
) -> Result[tuple[str, float, pl.DataFrame], str]:
    """Assemble the 4-panel (2x2) individual exemplar patients radar figure (2 R + 2 NR)."""
    match select_radar_marker_spokes(states_df, config.preset):
        case Failure(err):
            return Failure(err)
        case Success(spokes):
            pass

    match compute_compositional_radar_features(proj_fractions_df, spokes):
        case Failure(err):
            return Failure(err)
        case Success(proj_comp_df):
            pass

    match compute_oof_meta42_classifications(
        proj_comp_df,
        mode=config.classifier_mode,
        seed=config.seed,
    ):
        case Failure(err):
            return Failure(err)
        case Success(classified_df):
            pass

    match summarize_exemplar_patient_polygons(
        classified_df=classified_df,
        oof_landscape_df=oof_landscape_df,
        spokes=spokes,
        exemplars=EXEMPLAR_PATIENTS,
    ):
        case Failure(err):
            return Failure(err)
        case Success((pt_closed_df, pt_spoke_df, centroid_closed_df)):
            pass

    panels: list[alt.LayerChart] = []
    for ex in EXEMPLAR_PATIENTS:
        sub_closed = pt_closed_df.filter(pl.col("sample_id") == ex.sample_id)
        sub_spoke = pt_spoke_df.filter(pl.col("sample_id") == ex.sample_id)
        panels.append(
            build_single_patient_radar_panel(
                pt_closed_df=sub_closed,
                pt_spoke_df=sub_spoke,
                centroid_closed_df=centroid_closed_df,
                spokes=spokes,
                exemplar=ex,
            )
        )

    top_row = alt.hconcat(panels[0], panels[1], spacing=135).resolve_scale(
        color="independent",
        shape="independent",
    )
    bot_row = alt.hconcat(panels[2], panels[3], spacing=135).resolve_scale(
        color="independent",
        shape="independent",
    )
    combined = (
        alt.vconcat(top_row, bot_row, spacing=22, center=False)
        .resolve_scale(color="independent", shape="independent")
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75, fill=COLOR_BG_WHITE)
        .configure(
            font="system-ui, -apple-system, sans-serif",
            background=COLOR_BG_WHITE,
            padding={"left": 115, "right": 105, "top": 12, "bottom": 12},
        )
    )

    try:
        raw_svg = str(vlc.vegalite_to_svg(combined.to_dict()))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to serialize exemplar patient radar figure to SVG: {exc}")

    wrapped_svg, height_mm = wrap_svg_with_inkscape_layers(
        raw_svg,
        width_mm=FULL_WIDTH_MM,
        layer_id="layer-01-meta42-patient-exemplars-radar",
        layer_label="01_Meta42_Patient_Exemplars_Radar_Panels",
    )
    return Success((wrapped_svg, height_mm, pt_spoke_df))


def run_radar_figure_pipeline(
    config: RadarPlotConfig,
) -> Result[tuple[Path, Path, Path], str]:
    """Load Stage 3c parquets, build the Consensus-Meta42 radar figure, and export SVG/PNG/Parquet."""
    proj_path = (
        config.data_dir / "fractions" / "fractions_Consensus-Meta42_on_PanCancer-Combined.parquet"
    )
    raw_path = (
        config.data_dir
        / "fractions"
        / "raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet"
    )
    states_path = config.data_dir / "milopy_iatlas_state_effects.parquet"
    oof_path = config.data_dir / "oof_9a_9b_patient_landscape_and_curves.parquet"

    for req_path in (proj_path, raw_path, states_path, oof_path):
        if not req_path.exists():
            return Failure(f"Required Stage 3c Parquet not found: {req_path}")

    proj_df = pl.read_parquet(proj_path)
    raw_df = pl.read_parquet(raw_path)
    states_df = pl.read_parquet(states_path)
    oof_df = pl.read_parquet(oof_path)

    match build_meta42_responder_radar_figure(proj_df, raw_df, states_df, oof_df, config):
        case Failure(err):
            return Failure(err)
        case Success((svg_str, height_mm, summary_df)):
            pass

    match build_exemplar_patients_radar_figure(proj_df, states_df, oof_df, config):
        case Failure(err):
            return Failure(err)
        case Success((pt_svg_str, pt_height_mm, pt_summary_df)):
            pass

    config.results_dir.mkdir(parents=True, exist_ok=True)
    stem = "iatlas_consensus_metastate_responder_radar"
    svg_path = config.results_dir / f"{stem}.svg"
    png_path = config.results_dir / f"{stem}.png"
    summary_path = config.results_dir / "meta42_radar_profile_summary.parquet"
    summary_df.write_parquet(summary_path)

    pt_stem = "iatlas_consensus_metastate_patient_exemplars_radar"
    pt_svg_path = config.results_dir / f"{pt_stem}.svg"
    pt_png_path = config.results_dir / f"{pt_stem}.png"
    pt_summary_path = config.results_dir / "meta42_patient_exemplars_radar_summary.parquet"
    pt_summary_df.write_parquet(pt_summary_path)

    match export_svg_and_300dpi_png(
        svg_text=pt_svg_str,
        svg_path=pt_svg_path,
        png_path=pt_png_path,
        article_dir=config.article_dir,
        width_mm=FULL_WIDTH_MM,
        height_mm=pt_height_mm,
    ):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    export_res = export_svg_and_300dpi_png(
        svg_text=svg_str,
        svg_path=svg_path,
        png_path=png_path,
        article_dir=config.article_dir,
        width_mm=FULL_WIDTH_MM,
        height_mm=height_mm,
    )
    return export_res.map(lambda pair: (pair[0], pair[1], summary_path))


# ------------------------------------------------------------------------------
# 6. CLI & Imperative Shell
# ------------------------------------------------------------------------------


def parse_args(argv: Maybe[Sequence[str]] = Nothing) -> RadarPlotConfig:
    """Parse CLI arguments into an immutable RadarPlotConfig."""
    parser = argparse.ArgumentParser(
        description="Generate 180mm Nature Methods radar immunogram for Consensus-Meta42 R vs NR fractions."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark"),
        help="Directory containing Stage 3c benchmark Parquet files",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/4f_plot_meta42_responder_radar"),
        help="Directory to write output SVG, 300 DPI PNG, and summary Parquet",
    )
    parser.add_argument(
        "--article-dir",
        type=Path,
        default=Path("article/figures/milo-deconv"),
        help="Article figures directory for SVG/PNG synchronization",
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=[p.value for p in MarkerSelectionPreset],
        default=MarkerSelectionPreset.CURATED_12_IMMUNOGRAM.value,
        help="Marker state selection preset (curated_12 or top_14)",
    )
    parser.add_argument(
        "--radial-mode",
        type=str,
        choices=[m.value for m in RadialScaleMode],
        default=RadialScaleMode.WITHIN_TRIAL_CLR_Z.value,
        help="Radial coordinate scaling mode (within_trial_clr_z or rel_fold_pct)",
    )
    parser.add_argument(
        "--classifier",
        type=str,
        choices=[c.value for c in ClassifierMode],
        default=ClassifierMode.RIDGE_L2.value,
        help="Supervised classifier family for Panel (a) OOF predictions",
    )

    raw_args = parser.parse_args(argv.value_or(sys.argv[1:]))
    return RadarPlotConfig(
        data_dir=raw_args.data_dir,
        results_dir=raw_args.results_dir,
        article_dir=Some(raw_args.article_dir),
        preset=MarkerSelectionPreset(raw_args.preset),
        radial_mode=RadialScaleMode(raw_args.radial_mode),
        classifier_mode=ClassifierMode(raw_args.classifier),
    )


def main() -> None:
    """Imperative shell entrypoint."""
    config = parse_args()
    match run_radar_figure_pipeline(config):
        case Failure(err):
            print(f"Radar figure generation failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success((svg_p, png_p, parquet_p)):
            pt_stem = "iatlas_consensus_metastate_patient_exemplars_radar"
            print(f"Exported Inkscape SVG : {svg_p}")
            print(f"Exported 300 DPI PNG  : {png_p}")
            print(f"Exported Summary Table: {parquet_p}")
            print(f"Exported Patient SVG  : {config.results_dir / (pt_stem + '.svg')}")
            print(f"Exported Patient PNG  : {config.results_dir / (pt_stem + '.png')}")


if __name__ == "__main__":
    main()
