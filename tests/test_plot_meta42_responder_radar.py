#!/usr/bin/env python3
"""Property-based (Hypothesis) and unit tests for Stage 4f Consensus-Meta42 radar chart."""

from __future__ import annotations

import importlib.util
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pytest
from hypothesis import HealthCheck, given, settings, strategies as st
from returns.maybe import Nothing
from returns.pipeline import is_successful
from returns.result import Failure, Success

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "manuscript_pipeline"
    / "4f_plot_meta42_responder_radar.py"
)
SPEC = importlib.util.spec_from_file_location("plot_meta42_responder_radar", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
radar_mod: Any = importlib.util.module_from_spec(SPEC)
sys.modules["plot_meta42_responder_radar"] = radar_mod
SPEC.loader.exec_module(radar_mod)

ALL_META42_COLS: tuple[str, ...] = radar_mod.ALL_META42_COLS
CURATED_12_MARKER_ORDER: tuple[str, ...] = radar_mod.CURATED_12_MARKER_ORDER
MarkerSelectionPreset = radar_mod.MarkerSelectionPreset
RadialScaleMode = radar_mod.RadialScaleMode
ClassifierMode = radar_mod.ClassifierMode
SpokePole = radar_mod.SpokePole
RadarPlotConfig = radar_mod.RadarPlotConfig
select_radar_marker_spokes = radar_mod.select_radar_marker_spokes
compute_compositional_radar_features = radar_mod.compute_compositional_radar_features
compute_oof_meta42_classifications = radar_mod.compute_oof_meta42_classifications
summarize_radar_polygons = radar_mod.summarize_radar_polygons
wrap_svg_with_inkscape_layers = radar_mod.wrap_svg_with_inkscape_layers
run_radar_figure_pipeline = radar_mod.run_radar_figure_pipeline


def _make_synthetic_states_df() -> pl.DataFrame:
    """Create a deterministic 42-state Pan-iAtlas-Meta summary table."""
    rows: list[dict[str, object]] = []
    # Give positive betas to the 6 curated responder markers and negative betas to the 6 barrier markers
    curated_r = set(CURATED_12_MARKER_ORDER[:6])
    curated_nr = set(CURATED_12_MARKER_ORDER[6:])
    for idx, st_id in enumerate(ALL_META42_COLS):
        if st_id in curated_r:
            beta = 2.5 - 0.2 * float(CURATED_12_MARKER_ORDER.index(st_id))
            lfc = 0.5
        elif st_id in curated_nr:
            beta = -3.0 + 0.2 * float(CURATED_12_MARKER_ORDER.index(st_id) - 6)
            lfc = -0.5
        else:
            beta = 0.05 * math.sin(float(idx))
            lfc = 0.01 * math.cos(float(idx))
        rows.append({
            "reference_id": "Consensus-Meta42",
            "bulk_cohort": "Pan-iAtlas-Meta",
            "cell_state": st_id,
            "bulk_deconv_beta": beta,
            "milo_mean_logfc": lfc,
        })
    return pl.DataFrame(rows)


def _make_synthetic_fractions_df(
    raw_weights: np.ndarray,
    cohorts: list[str],
    responses: list[int],
) -> pl.DataFrame:
    """Construct a valid compositional 42-state fractions Polars DataFrame."""
    row_sums = np.maximum(raw_weights.sum(axis=1, keepdims=True), 1e-12)
    norm_weights = raw_weights / row_sums
    data: dict[str, Any] = {
        "sample_id": [f"S_{i:03d}" for i in range(len(cohorts))],
        "cohort": cohorts,
        "response_binary": responses,
    }
    for j, col_name in enumerate(ALL_META42_COLS):
        data[col_name] = norm_weights[:, j].tolist()
    return pl.DataFrame(data)


def test_select_radar_marker_spokes_curated_12_and_top_14() -> None:
    """Verify spoke selection geometry, pole balance, and angular spacing for both presets."""
    states_df = _make_synthetic_states_df()

    res_12 = select_radar_marker_spokes(states_df, MarkerSelectionPreset.CURATED_12_IMMUNOGRAM)
    assert is_successful(res_12)
    spokes_12 = res_12.unwrap()
    assert len(spokes_12) == 12
    assert tuple(s.state_id for s in spokes_12) == CURATED_12_MARKER_ORDER
    # First 6 are Responder Niche (0° to 150°), next 6 are Non-Responder Barrier (180° to 330°)
    for idx, sp in enumerate(spokes_12):
        assert math.isclose(sp.angle_deg, 30.0 * float(idx), abs_tol=1e-9)
        expected_pole = (
            SpokePole.RESPONDER_NICHE if idx < 6 else SpokePole.NON_RESPONDER_BARRIER
        )
        assert sp.pole == expected_pole

    res_14 = select_radar_marker_spokes(states_df, MarkerSelectionPreset.TOP_14_META_POLES)
    assert is_successful(res_14)
    spokes_14 = res_14.unwrap()
    assert len(spokes_14) == 14
    assert sum(1 for s in spokes_14 if s.pole == SpokePole.RESPONDER_NICHE) == 7
    assert sum(1 for s in spokes_14 if s.pole == SpokePole.NON_RESPONDER_BARRIER) == 7


@settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
@given(
    n_per_cohort=st.integers(min_value=6, max_value=14),
    seed=st.integers(min_value=1, max_value=10_000),
    scale_factor=st.floats(min_value=0.25, max_value=4.0, allow_nan=False, allow_infinity=False),
)
def test_property_compositional_clr_and_subsimplex_invariants(
    n_per_cohort: int,
    seed: int,
    scale_factor: float,
) -> None:
    """Property test:
    1. Row-wise CLR values sum to 0 across all 42 states.
    2. Within-trial CLR-Z values have mean 0 and unit SD within each cohort.
    3. Sub-simplex percentages sum to 100.0% across the K selected spokes for every patient.
    4. Positive multiplicative scaling of raw weights leaves sub-simplex % invariant.
    """
    rng = np.random.default_rng(seed)
    n_total = 2 * n_per_cohort
    raw_weights = rng.uniform(0.005, 0.10, size=(n_total, 42))
    cohorts = ["Trial_A"] * n_per_cohort + ["Trial_B"] * n_per_cohort
    responses = [i % 2 for i in range(n_total)]

    states_df = _make_synthetic_states_df()
    spokes = select_radar_marker_spokes(
        states_df, MarkerSelectionPreset.CURATED_12_IMMUNOGRAM
    ).unwrap()

    df_base = _make_synthetic_fractions_df(raw_weights, cohorts, responses)
    df_scaled = _make_synthetic_fractions_df(raw_weights * scale_factor, cohorts, responses)

    res_base = compute_compositional_radar_features(df_base, spokes)
    res_scaled = compute_compositional_radar_features(df_scaled, spokes)
    assert is_successful(res_base)
    assert is_successful(res_scaled)

    comp_base = res_base.unwrap()
    comp_scaled = res_scaled.unwrap()

    # 1. Row-wise CLR sum == 0
    clr_mat = comp_base.select([f"clr_{c}" for c in ALL_META42_COLS]).to_numpy()
    assert np.allclose(clr_mat.sum(axis=1), 0.0, atol=1e-6)

    # 2. Within-cohort standardized CLR-Z has zero mean and unit variance per cohort
    for cid in ("Trial_A", "Trial_B"):
        sub = comp_base.filter(pl.col("cohort") == cid)
        z_mat = sub.select([f"z_{c}" for c in ALL_META42_COLS]).to_numpy()
        assert np.allclose(z_mat.mean(axis=0), 0.0, atol=1e-6)
        assert np.allclose(z_mat.std(axis=0, ddof=0), 1.0, atol=1e-5)

    # 3. Sub-simplex percentages sum to 100% per sample
    spoke_sub_cols = [f"sub_pct_{s.state_id}" for s in spokes]
    sub_mat_base = comp_base.select(spoke_sub_cols).to_numpy()
    sub_mat_scaled = comp_scaled.select(spoke_sub_cols).to_numpy()
    assert np.allclose(sub_mat_base.sum(axis=1), 100.0, atol=1e-6)

    # 4. Scale invariance of normalized fractions and sub-simplex %
    assert np.allclose(sub_mat_base, sub_mat_scaled, atol=1e-6)


@settings(
    max_examples=20,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
@given(
    seed=st.integers(min_value=1, max_value=5_000),
    val_bound=st.floats(min_value=0.35, max_value=1.50, allow_nan=False, allow_infinity=False),
)
def test_property_radar_polygon_closure_and_radial_bounds(
    seed: int,
    val_bound: float,
) -> None:
    """Property test:
    1. Every group polygon has K + 1 rows where vertex_order=0 and vertex_order=K have identical (x, y).
    2. All mapped radii lie strictly within [r_inner, r_outer].
    3. BH-adjusted q-values satisfy 0 <= p_value <= q_bh <= 1.0.
    """
    rng = np.random.default_rng(seed)
    n_samples = 24
    raw_weights = rng.dirichlet(np.ones(42), size=n_samples)
    cohorts = ["Trial_A"] * 12 + ["Trial_B"] * 12
    responses = [1 if i < 6 or (12 <= i < 18) else 0 for i in range(n_samples)]

    states_df = _make_synthetic_states_df()
    spokes = select_radar_marker_spokes(
        states_df, MarkerSelectionPreset.CURATED_12_IMMUNOGRAM
    ).unwrap()
    frac_df = _make_synthetic_fractions_df(raw_weights, cohorts, responses)
    comp_df = compute_compositional_radar_features(frac_df, spokes).unwrap().with_columns(
        pl.when(pl.col("response_binary") == 1)
        .then(pl.lit("Responder"))
        .otherwise(pl.lit("Non-Responder"))
        .alias("grp")
    )

    r_inner = 0.18
    r_outer = 0.88
    poly_res = summarize_radar_polygons(
        compositional_df=comp_df,
        group_col="grp",
        group_order=("Responder", "Non-Responder"),
        spokes=spokes,
        radial_mode=RadialScaleMode.WITHIN_TRIAL_CLR_Z,
        val_min=-val_bound,
        val_max=val_bound,
        r_inner=r_inner,
        r_outer=r_outer,
    )
    assert is_successful(poly_res)
    closed_df, stats_df = poly_res.unwrap()

    k_spokes = len(spokes)
    assert closed_df.height == 2 * (k_spokes + 1)
    assert stats_df.height == 2 * k_spokes

    for g_name in ("Responder", "Non-Responder"):
        g_poly = closed_df.filter(pl.col("group") == g_name).sort("vertex_order")
        assert g_poly.height == k_spokes + 1
        first_row = g_poly.row(0, named=True)
        last_row = g_poly.row(k_spokes, named=True)
        assert first_row["vertex_order"] == 0
        assert last_row["vertex_order"] == k_spokes
        assert math.isclose(float(first_row["x"]), float(last_row["x"]), abs_tol=1e-12)
        assert math.isclose(float(first_row["y"]), float(last_row["y"]), abs_tol=1e-12)

    radii = stats_df["radius"].to_numpy()
    r_lows = stats_df["r_low"].to_numpy()
    r_highs = stats_df["r_high"].to_numpy()
    assert np.all((radii >= r_inner - 1e-9) & (radii <= r_outer + 1e-9))
    assert np.all((r_lows >= r_inner - 1e-9) & (r_lows <= r_outer + 1e-9))
    assert np.all((r_highs >= r_inner - 1e-9) & (r_highs <= r_outer + 1e-9))

    p_vals = stats_df["p_value"].to_numpy()
    q_vals = stats_df["q_bh"].to_numpy()
    assert np.all((p_vals >= 0.0) & (p_vals <= 1.0))
    assert np.all((q_vals >= p_vals - 1e-12) & (q_vals <= 1.0 + 1e-12))


def test_end_to_end_svg_xml_layers_and_parquet_integrity(tmp_path: Path) -> None:
    """Verify full Stage 4f pipeline execution on cached Stage 3c artifacts and validate SVG/XML."""
    data_dir = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    if not data_dir.exists():
        pytest.skip("Stage 3c benchmark directory not present.")

    config = RadarPlotConfig(
        data_dir=data_dir,
        results_dir=tmp_path,
        article_dir=Nothing,
        preset=MarkerSelectionPreset.CURATED_12_IMMUNOGRAM,
        radial_mode=RadialScaleMode.WITHIN_TRIAL_CLR_Z,
        classifier_mode=ClassifierMode.RIDGE_L2,
    )
    res = run_radar_figure_pipeline(config)
    assert is_successful(res)
    svg_path, png_path, parquet_path = res.unwrap()

    assert svg_path.exists() and svg_path.stat().st_size > 10_000
    assert png_path.exists() and png_path.stat().st_size > 50_000
    assert parquet_path.exists()

    svg_text = svg_path.read_text(encoding="utf-8")
    assert 'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"' in svg_text
    assert 'inkscape:groupmode="layer"' in svg_text
    assert "feDropShadow" not in svg_text

    # Parse XML to guarantee strict well-formedness
    root = ET.fromstring(svg_text)
    assert root.tag.endswith("svg")
    assert root.attrib.get("width") == "180.0mm"

    summary_df = pl.read_parquet(parquet_path)
    assert summary_df.height == 84  # (2 + 2 + 3 groups) * 12 spokes
    assert int(summary_df.null_count().sum_horizontal().item()) == 0

    pt_svg_path = tmp_path / "iatlas_consensus_metastate_patient_exemplars_radar.svg"
    pt_png_path = tmp_path / "iatlas_consensus_metastate_patient_exemplars_radar.png"
    pt_parquet_path = tmp_path / "meta42_patient_exemplars_radar_summary.parquet"

    assert pt_svg_path.exists() and pt_svg_path.stat().st_size > 10_000
    assert pt_png_path.exists() and pt_png_path.stat().st_size > 50_000
    assert pt_parquet_path.exists()

    pt_svg_text = pt_svg_path.read_text(encoding="utf-8")
    assert 'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"' in pt_svg_text
    assert 'inkscape:groupmode="layer"' in pt_svg_text
    assert "feDropShadow" not in pt_svg_text
    pt_root = ET.fromstring(pt_svg_text)
    assert pt_root.tag.endswith("svg")
    assert pt_root.attrib.get("width") == "180.0mm"

    pt_summary_df = pl.read_parquet(pt_parquet_path)
    assert pt_summary_df.height == 48  # 4 exemplar patients * 12 spokes
    assert int(pt_summary_df.null_count().sum_horizontal().item()) == 0
    for sample_id in pt_summary_df["sample_id"].unique().to_list():
        sub = pt_summary_df.filter(pl.col("sample_id") == sample_id)
        assert math.isclose(float(sub["sub_simplex_pct"].sum()), 100.0, abs_tol=1e-5)
        r_share = float(sub["r_pole_sub_pct"][0])
        nr_share = float(sub["nr_pole_sub_pct"][0])
        assert math.isclose(r_share + nr_share, 100.0, abs_tol=1e-5)
