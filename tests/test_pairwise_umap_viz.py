"""Unit, Hypothesis property-based, and hybrid Inkscape SVG tests for pairwise UMAP visualization.

Verifies:
1. Manifold-faithful geometric medoids (`compute_nhood_umap_medoids`): every medoid strictly
   coincides with an actual observed member cell `x_{i*} in C_u` (even across torn bimodal UMAP
   islands where centroids float in empty space) and transforms equivariantly under translation
   and positive isotropic scaling.
2. Outward-bowing quadratic Bezier chords (`discretize_outward_bezier_arc`): exact endpoint
   interpolation `B(0) = m_u`, `B(1) = m_v`, outward bow away from global UMAP median, and exact
   `(u, v) <-> (v, u)` swap invariance `B_{vu}(1 - t) == B_{uv}(t)` + polarity swap symmetry.
3. De novo pairwise hub isolation (`summarize_nhood_pairwise_hubness`) and signed 1-to-many
   seed-conditioned partner field (`compute_anchor_partner_field`).
4. End-to-end 4-panel composite & cross-cohort grid export (`run_pairwise_umap_pipeline`): valid
   Inkscape layer XML (`01_Background`, `02_Pairwise_UMAP_Panels`), 300 DPI rasterized dense cell
   clouds (`<image>`), preserved pure vector Bezier `<path>` arcs & `<text>` labels, and zero
   drop shadows.
"""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Result, Success
from scipy import sparse as sp

from milopy.pairwise import (
    ALL_QUADRANT_STATES,
    EnrichmentDirection,
    HubNoveltyStatus,
    NhoodCoordinateMethod,
    PairwiseAnalysisResult,
    PairwiseUmapBundle,
    QuadrantState,
    UmapPairVizConfig,
    compute_anchor_partner_field,
    compute_nhood_umap_medoids,
    discretize_outward_bezier_arc,
    summarize_nhood_pairwise_hubness,
    swap_quadrant,
)
from milopy.pairwise.umap_viz import classify_quadrant_visual_regime_and_polarity
from plotting_utils.svg import INKSCAPE_NS, SVG_NS  # type: ignore[import-untyped]
from scripts.milopy_analysis.visualize_pairwise_umap import (
    CohortRawVizArtifacts,
    PairRepresentationMode,
    PairwiseUmapCliConfig,
    run_pairwise_umap_pipeline,
)


@settings(max_examples=35)
@given(
    seed=st.integers(min_value=1, max_value=50000),
    tx=st.floats(min_value=-25.0, max_value=25.0, allow_nan=False, allow_infinity=False),
    ty=st.floats(min_value=-25.0, max_value=25.0, allow_nan=False, allow_infinity=False),
    scale=st.floats(min_value=0.25, max_value=6.0, allow_nan=False, allow_infinity=False),
)
def test_prop_geometric_medoid_manifold_membership_and_affine_equivariance(
    seed: int,
    tx: float,
    ty: float,
    scale: float,
) -> None:
    """Proves every medoid is an exact observed member cell on the UMAP manifold & affine-equivariant."""
    rng = np.random.default_rng(seed)
    # Two torn UMAP islands centered at (-10, 0) and (+10, 0) with empty void around (0, 0)
    island_left = rng.normal(loc=[-10.0, 0.0], scale=0.6, size=(30, 2))
    island_right = rng.normal(loc=[10.0, 0.0], scale=0.6, size=(30, 2))
    coords = np.vstack([island_left, island_right])

    # Neighborhood 0 spans left island; Neighborhood 1 spans right island;
    # Neighborhood 2 is split across the tear (17 left + 13 right) -> arithmetic mean lands in void!
    rows_0 = list(range(0, 15))
    rows_1 = list(range(30, 45))
    rows_2 = list(range(10, 27)) + list(range(45, 58))
    all_rows = rows_0 + rows_1 + rows_2
    all_cols = [0] * len(rows_0) + [1] * len(rows_1) + [2] * len(rows_2)
    nhoods_mat = sp.csc_matrix(
        (np.ones(len(all_rows), dtype=np.float64), (all_rows, all_cols)),
        shape=(60, 3),
    )

    med_orig = compute_nhood_umap_medoids(
        nhoods_mat=nhoods_mat,
        umap_coords=coords,
        method=NhoodCoordinateMethod.GEOMETRIC_MEDOID,
    ).unwrap()

    member_sets = (set(rows_0), set(rows_1), set(rows_2))
    for u_idx in range(3):
        row_u = med_orig.row(u_idx, named=True)
        chosen_cell = int(row_u["medoid_cell_idx"])
        assert chosen_cell in member_sets[u_idx]
        mx = float(row_u["umap_1"])
        my = float(row_u["umap_2"])
        assert np.isclose(mx, coords[chosen_cell, 0], atol=1e-12)
        assert np.isclose(my, coords[chosen_cell, 1], atol=1e-12)
        # Never lands in the empty tear [-6, +6]
        assert abs(mx) > 6.0

    transformed = coords * scale + np.array([[tx, ty]], dtype=np.float64)
    med_trans = compute_nhood_umap_medoids(
        nhoods_mat=nhoods_mat,
        umap_coords=transformed,
        method=NhoodCoordinateMethod.GEOMETRIC_MEDOID,
    ).unwrap()

    assert med_trans["medoid_cell_idx"].to_list() == med_orig["medoid_cell_idx"].to_list()
    np.testing.assert_allclose(
        med_trans["umap_1"].to_numpy(),
        med_orig["umap_1"].to_numpy() * scale + tx,
        rtol=1e-10,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        med_trans["umap_2"].to_numpy(),
        med_orig["umap_2"].to_numpy() * scale + ty,
        rtol=1e-10,
        atol=1e-10,
    )


@settings(max_examples=45)
@given(
    xu=st.floats(min_value=-15.0, max_value=15.0, allow_nan=False, allow_infinity=False),
    yu=st.floats(min_value=-15.0, max_value=15.0, allow_nan=False, allow_infinity=False),
    dx=st.floats(min_value=0.5, max_value=12.0, allow_nan=False, allow_infinity=False),
    dy=st.floats(min_value=-12.0, max_value=12.0, allow_nan=False, allow_infinity=False),
    gx=st.floats(min_value=-10.0, max_value=10.0, allow_nan=False, allow_infinity=False),
    gy=st.floats(min_value=-10.0, max_value=10.0, allow_nan=False, allow_infinity=False),
    curvature=st.floats(min_value=0.10, max_value=0.45, allow_nan=False, allow_infinity=False),
    q_idx=st.integers(min_value=0, max_value=3),
    is_resp=st.booleans(),
)
def test_prop_bezier_chord_endpoints_outward_bow_and_index_swap_invariance(
    xu: float,
    yu: float,
    dx: float,
    dy: float,
    gx: float,
    gy: float,
    curvature: float,
    q_idx: int,
    is_resp: bool,
) -> None:
    """Proves B(0)=m_u, B(1)=m_v, outward bow from global center, and exact (u,v)<->(v,u) invariance."""
    xv = xu + dx
    yv = yu + dy
    bx_uv, by_uv = discretize_outward_bezier_arc(
        p_u=(xu, yu),
        p_v=(xv, yv),
        global_center=(gx, gy),
        curvature=curvature,
        n_steps=21,
    )
    bx_vu, by_vu = discretize_outward_bezier_arc(
        p_u=(xv, yv),
        p_v=(xu, yu),
        global_center=(gx, gy),
        curvature=curvature,
        n_steps=21,
    )

    # 1. Exact endpoint interpolation
    assert np.isclose(bx_uv[0], xu, atol=1e-12) and np.isclose(by_uv[0], yu, atol=1e-12)
    assert np.isclose(bx_uv[-1], xv, atol=1e-12) and np.isclose(by_uv[-1], yv, atol=1e-12)

    # 2. Outward bow: Bezier midpoint t=0.5 sits strictly further from global_center than straight midpoint
    mid_straight_x = 0.5 * (xu + xv)
    mid_straight_y = 0.5 * (yu + yv)
    d_straight = float(np.hypot(mid_straight_x - gx, mid_straight_y - gy))
    d_bezier = float(np.hypot(bx_uv[10] - gx, by_uv[10] - gy))
    assert d_bezier + 1e-10 >= d_straight

    # 3. Exact geometric reversal invariance under (u, v) <-> (v, u)
    np.testing.assert_allclose(bx_uv, bx_vu[::-1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(by_uv, by_vu[::-1], rtol=1e-12, atol=1e-12)

    # 4. Exact regime & endpoint polarity swap invariance under swap_quadrant(q)
    q_state = ALL_QUADRANT_STATES[q_idx]
    direction = EnrichmentDirection.RESPONDER if is_resp else EnrichmentDirection.NON_RESPONDER
    reg_uv, pol_u, pol_v = classify_quadrant_visual_regime_and_polarity(q_state, direction)
    reg_vu, swp_v, swp_u = classify_quadrant_visual_regime_and_polarity(
        swap_quadrant(q_state), direction
    )
    assert reg_uv == reg_vu
    assert pol_u == swp_u
    assert pol_v == swp_v


def test_denovo_hub_and_anchor_partner_projection() -> None:
    """Verifies DE_NOVO_HUB isolation (univariate FDR >= 0.10) and signed anchor partner Z field."""
    n_cells = 60
    n_nhoods = 6
    rows: list[int] = []
    cols: list[int] = []
    coords = np.zeros((n_cells, 2), dtype=np.float64)
    island_centers = (
        (-8.0, 6.0),
        (8.0, 6.0),
        (-8.0, -6.0),
        (8.0, -6.0),
        (0.0, 9.0),
        (0.0, -9.0),
    )
    for u in range(n_nhoods):
        cx, cy = island_centers[u]
        for loc in range(10):
            c_idx = u * 10 + loc
            coords[c_idx, 0] = cx + 0.1 * float(loc)
            coords[c_idx, 1] = cy + 0.05 * float(loc % 3)
            rows.append(c_idx)
            cols.append(u)

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, n_nhoods),
    )
    state_labels = ("M29_Plasma", "M30_TLS_B", "M04_Eff_T", "M17_myCAF", "M05_Mem_T", "M32_TAM")
    # Univ FDR: 0, 1, 3, 5 are NS (>= 0.10); 2 and 4 are Univ Sig (< 0.10)
    univ_fdr = (0.45, 0.55, 0.01, 0.65, 0.02, 0.85)

    pairs_df = pl.DataFrame(
        {
            "nhood_u": [0, 2],
            "nhood_v": [1, 3],
            "quadrant": [QuadrantState.BOTH_PRESENT.value, QuadrantState.FIRST_ONLY.value],
            "direction": [EnrichmentDirection.RESPONDER.value, EnrichmentDirection.RESPONDER.value],
            "log_or_quadrant": [3.2, 2.9],
            "pvalue_iut": [0.0004, 0.0008],
            "fdr_iut": [0.008, 0.012],
            "is_significant": [True, True],
        }
    )
    cont_df = pl.DataFrame(
        schema={
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "axis": pl.String,
            "direction": pl.String,
            "effect_size": pl.Float64,
            "pvalue_iut": pl.Float64,
            "fdr_iut": pl.Float64,
            "is_significant": pl.Boolean,
        }
    )
    modules_df = pl.DataFrame(
        {
            "pair_module_id": [0, 1],
            "lead_nhood_u": [0, 2],
            "lead_nhood_v": [1, 3],
            "quadrant": [QuadrantState.BOTH_PRESENT.value, QuadrantState.FIRST_ONLY.value],
            "direction": [EnrichmentDirection.RESPONDER.value, EnrichmentDirection.RESPONDER.value],
            "best_pvalue_iut": [0.0004, 0.0008],
            "best_fdr_iut": [0.008, 0.012],
        }
    )
    pw_res = PairwiseAnalysisResult(
        pairs_df=pairs_df,
        continuous_pairs_df=cont_df,
        modules_df=modules_df,
        n_total_nhoods=n_nhoods,
        n_active_nhoods=n_nhoods,
        n_evaluated_pairs=15,
        n_significant_quadrant_pairs=2,
        n_significant_continuous_pairs=0,
        n_modules=2,
    )

    medoids_df = compute_nhood_umap_medoids(
        nhoods_mat=nhoods_mat,
        umap_coords=coords,
        state_labels=Some(state_labels),
    ).unwrap()
    nodes_df, _, bg_df = summarize_nhood_pairwise_hubness(
        nhoods_mat=nhoods_mat,
        umap_coords=coords,
        medoids_df=medoids_df,
        pairwise_res=pw_res,
        univariate_fdr=Some(univ_fdr),
    ).unwrap()

    statuses = nodes_df["hub_novelty_status"].to_list()
    assert statuses[0] == HubNoveltyStatus.DE_NOVO_HUB.value
    assert statuses[1] == HubNoveltyStatus.DE_NOVO_HUB.value
    assert statuses[2] == HubNoveltyStatus.SHARED_HUB.value
    assert statuses[3] == HubNoveltyStatus.DE_NOVO_HUB.value
    assert statuses[4] == HubNoveltyStatus.UNIVARIATE_ONLY.value
    assert statuses[5] == HubNoveltyStatus.INACTIVE.value

    # Anchor on nhood 0 (synergy with 1) vs nhood 2 (gating antagonism against 3)
    anch_lbl_syn, anch_bg_syn, _, anch_arcs_syn = compute_anchor_partner_field(
        nhoods_mat=nhoods_mat,
        nhood_nodes_df=nodes_df,
        bg_cells_df=bg_df,
        pairwise_res=pw_res,
        anchor_nhood=Some(0),
    ).unwrap()
    assert "M29_Plasma" in anch_lbl_syn
    assert anch_arcs_syn.height > 0
    cells_nh1 = anch_bg_syn.filter(pl.col("cell_idx").is_between(10, 19))["anchor_partner_z"].to_numpy()
    assert bool(np.all(cells_nh1 > 2.0))

    _, anch_bg_gate, _, _ = compute_anchor_partner_field(
        nhoods_mat=nhoods_mat,
        nhood_nodes_df=nodes_df,
        bg_cells_df=bg_df,
        pairwise_res=pw_res,
        anchor_nhood=Some(2),
    ).unwrap()
    cells_nh3 = anch_bg_gate.filter(pl.col("cell_idx").is_between(30, 39))["anchor_partner_z"].to_numpy()
    assert bool(np.all(cells_nh3 < -2.0))


def test_end_to_end_composite_svg_hybrid_raster_vector_integrity(tmp_path: Path) -> None:
    """Verifies 4-panel composite & grid export valid Inkscape XML with hybrid <image> + vector <path>."""
    rng = np.random.default_rng(2026)
    n_donors = 24
    cells_per_donor = 65  # 1,560 cells -> exceeds max_points=1000 & min_group_points=250!
    n_cells = n_donors * cells_per_donor
    n_nhoods = 6

    island_xy = np.array(
        [
            [-9.0, 6.0],
            [9.0, 6.0],
            [-9.0, -6.0],
            [9.0, -6.0],
            [0.0, 8.5],
            [0.0, -8.5],
        ],
        dtype=np.float64,
    )
    umap_1 = np.empty(n_cells, dtype=np.float64)
    umap_2 = np.empty(n_cells, dtype=np.float64)
    pt_ids: list[str] = []
    resps: list[str] = []
    cts: list[str] = []
    rows: list[int] = []
    cols: list[int] = []

    ct_names = ("Plasma", "TLS_B", "Eff_T", "myCAF", "Mem_T", "TAM")
    for d_i in range(n_donors):
        is_r = d_i < 12
        r_lbl = "responder" if is_r else "non-responder"
        for loc in range(cells_per_donor):
            c_idx = d_i * cells_per_donor + loc
            nh_id = loc % n_nhoods
            umap_1[c_idx] = float(island_xy[nh_id, 0] + rng.normal(0.0, 0.45))
            umap_2[c_idx] = float(island_xy[nh_id, 1] + rng.normal(0.0, 0.45))
            pt_ids.append(f"PT_{d_i:02d}")
            resps.append(r_lbl)
            cts.append(ct_names[nh_id])
            # Donor quadrant structure: R donors have nhood 0 & 1 present and nhood 3 absent
            keep_in_nh = True
            if nh_id in (0, 1) and not is_r and (d_i % 3 != 0):
                keep_in_nh = False
            if nh_id == 3 and is_r and (d_i % 4 != 0):
                keep_in_nh = False
            if keep_in_nh:
                rows.append(c_idx)
                cols.append(nh_id)

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, n_nhoods),
    )
    cells_df = pl.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_cells)],
            "patient_id": pt_ids,
            "clinical_response": resps,
            "cell_type": cts,
            "UMAP1": umap_1,
            "UMAP2": umap_2,
        }
    )
    da_df = pl.DataFrame(
        {
            "Nhood": list(range(n_nhoods)),
            "logFC": [0.9, 0.8, 1.2, -1.3, 0.1, -0.2],
            "FDR": [0.22, 0.30, 0.03, 0.18, 0.70, 0.80],
            "Nhood_CellType": list(ct_names),
        }
    )
    q_df = pl.DataFrame(
        {
            "nhood_u": [0, 2],
            "nhood_v": [1, 3],
            "quadrant": [QuadrantState.BOTH_PRESENT.value, QuadrantState.FIRST_ONLY.value],
            "direction": [EnrichmentDirection.RESPONDER.value, EnrichmentDirection.RESPONDER.value],
            "log_or_quadrant": [3.1, 2.8],
            "pvalue_iut": [0.0005, 0.0012],
            "fdr_iut": [0.009, 0.018],
            "is_significant": [True, True],
        }
    )
    c_df = pl.DataFrame(
        schema={
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "axis": pl.String,
            "direction": pl.String,
            "effect_size": pl.Float64,
            "pvalue_iut": pl.Float64,
            "fdr_iut": pl.Float64,
            "is_significant": pl.Boolean,
        }
    )
    m_df = pl.DataFrame(
        {
            "pair_module_id": [0, 1],
            "lead_nhood_u": [0, 2],
            "lead_nhood_v": [1, 3],
            "quadrant": [QuadrantState.BOTH_PRESENT.value, QuadrantState.FIRST_ONLY.value],
            "direction": [EnrichmentDirection.RESPONDER.value, EnrichmentDirection.RESPONDER.value],
            "best_pvalue_iut": [0.0005, 0.0012],
            "best_fdr_iut": [0.009, 0.018],
        }
    )
    meta_sub = pl.DataFrame(
        {
            "accession": ["SYNTH_UMAP_01"] * n_nhoods,
            "nhood_idx": list(range(n_nhoods)),
            "meta_state": ["M29", "M30", "M04", "M17", "M05", "M32"],
            "meta_state_label": [
                "M29: JCHAIN+MZB1+ Plasma",
                "M30: CR2+PAX5+ TLS B",
                "M04: GNLY+PRF1+ Eff. T",
                "M17: ISLR+COL5A2+ myCAF",
                "M05: CD96+TRAF3+ Mem. T",
                "M32: MRC1+CD163+ M2-TAM",
            ],
        }
    )

    def _stub_loader(
        accession: str,
        config: PairwiseUmapCliConfig,
        all_meta42_df: pl.DataFrame,
    ) -> Result[CohortRawVizArtifacts, str]:
        _ = (config, all_meta42_df)
        return Success(
            CohortRawVizArtifacts(
                accession=accession,
                nhoods_mat=nhoods_mat,
                cells_df=cells_df,
                da_results_df=da_df,
                t0_quadrant_df=q_df,
                t0_continuous_df=c_df,
                t0_modules_df=m_df,
                meta42_sub_df=meta_sub,
            )
        )

    out_dir = tmp_path / "umap_out"
    cli_cfg = PairwiseUmapCliConfig(
        output_dir=out_dir,
        pair_mode=PairRepresentationMode.RAW_MODULES,
        cohort_filter=Some(("SYNTH_UMAP_01",)),
        max_chords=16,
        max_bg_cells=1560,
    )
    summaries = run_pairwise_umap_pipeline(cli_cfg=cli_cfg, loader=_stub_loader).unwrap()
    assert len(summaries) == 1
    assert summaries[0].n_displayed_chords == 2

    comp_svg = out_dir / "SYNTH_UMAP_01" / "SYNTH_UMAP_01_raw_modules_umap_composite.svg"
    comp_png = out_dir / "SYNTH_UMAP_01" / "SYNTH_UMAP_01_raw_modules_umap_composite.png"
    grid_svg = out_dir / "figures" / "nine_cohorts_pairwise_umap_grid.svg"
    grid_png = out_dir / "figures" / "nine_cohorts_pairwise_umap_grid.png"

    for p in (comp_svg, comp_png, grid_svg, grid_png):
        assert p.exists() and p.stat().st_size > 1000

    svg_text = comp_svg.read_text(encoding="utf-8")
    assert "feDropShadow" not in svg_text
    assert "xmlns:inkscape=" in svg_text
    assert "xmlns:sodipodi=" in svg_text

    root = ET.fromstring(svg_text)
    ink_layers = [
        el.attrib.get(f"{{{INKSCAPE_NS}}}label", "")
        for el in root.iter()
        if el.attrib.get(f"{{{INKSCAPE_NS}}}groupmode") == "layer"
    ]
    assert ink_layers == ["01_Background", "02_Pairwise_UMAP_Panels"]

    # Dense single-cell UMAP clouds (>= 250 marks) rasterized into <image>; Bezier chords & text stay vector
    image_elems = [el for el in root.iter() if el.tag == f"{{{SVG_NS}}}image"]
    path_elems = [el for el in root.iter() if el.tag == f"{{{SVG_NS}}}path"]
    text_elems = [el for el in root.iter() if el.tag == f"{{{SVG_NS}}}text"]
    assert len(image_elems) >= 4
    assert len(path_elems) >= 10
    assert len(text_elems) >= 10


def test_super_meta42_default_cross_compartment_mixing_qc_and_cmh_strata(tmp_path: Path) -> None:
    """Verifies SUPER_META42 defaults enforce cross-compartment gating, donor mixing QC, CMH strata & 1/state-pair."""
    rng = np.random.default_rng(777)
    # 60 donors across 2 <PMID>ggg<Patient> sub-studies with imbalanced R rates (20% vs 80%):
    #   S1 (27124452): 6 R + 24 NR = 30 donors
    #   S2 (33958794): 24 R + 6 NR = 30 donors
    donor_ids: list[str] = []
    donor_resps: list[str] = []
    donor_studies: list[str] = []
    for i in range(30):
        donor_ids.append(f"27124452gggPtS1_{i:02d}")
        donor_resps.append("responder" if i < 6 else "non-responder")
        donor_studies.append("27124452")
    for i in range(30):
        donor_ids.append(f"33958794gggPtS2_{i:02d}")
        donor_resps.append("responder" if i < 24 else "non-responder")
        donor_studies.append("33958794")

    # 7 cell-disjoint neighborhoods on distinct UMAP islands:
    #   0: M04 (Eff_T, T_NK) — primary true R+-gated partner with 1 (M17 myCAF)
    #   1: M17 (myCAF, STROMAL_MESENCHYMAL) — excluded in R, enriched in NR
    #   2: M04 (Eff_T, T_NK) — parallel sub-patch of M04 with slightly weaker R+- gate vs 1 -> tests 1/state-pair!
    #   3: M03 (Tfh_T, T_NK) & 4: M05 (Mem_T, T_NK) — intra-T_NK AND purely S2-study-confounded!
    #   5: M30 (TLS_B, B_LINEAGE) — single-donor private clone (all cells in donor 30) -> tests MixingQCConfig!
    #   6: M10 (Mono/cDC, MYELOID_MAST) — uniform background across all 60 donors
    n_nhoods = 7
    island_xy = np.array(
        [
            [-6.5, 4.5],
            [6.5, -4.5],
            [-5.2, 5.4],
            [-6.8, 2.0],
            [-4.8, 1.8],
            [5.8, 4.8],
            [0.0, -5.5],
        ],
        dtype=np.float64,
    )
    ct_names = ("Eff_T", "myCAF", "Eff_T", "Tfh_T", "Mem_T", "TLS_B", "Mono")

    umap_1_list: list[float] = []
    umap_2_list: list[float] = []
    pt_col: list[str] = []
    resp_col: list[str] = []
    ct_col: list[str] = []
    mat_rows: list[int] = []
    mat_cols: list[int] = []

    def _emit_cells(d_idx: int, nh_idx: int, n_emit: int) -> None:
        for _ in range(n_emit):
            c_id = len(pt_col)
            umap_1_list.append(float(island_xy[nh_idx, 0] + rng.normal(0.0, 0.35)))
            umap_2_list.append(float(island_xy[nh_idx, 1] + rng.normal(0.0, 0.35)))
            pt_col.append(donor_ids[d_idx])
            resp_col.append(donor_resps[d_idx])
            ct_col.append(ct_names[nh_idx])
            mat_rows.append(c_id)
            mat_cols.append(nh_idx)

    for d_idx in range(len(donor_ids)):
        is_r = donor_resps[d_idx] == "responder"
        st_id = donor_studies[d_idx]
        # Uniform baseline neighborhood 6 so every donor has library_size >= 30
        _emit_cells(d_idx, 6, 32)

        # True within-stratum cross-compartment +- gate between (0/2: M04) and (1: M17):
        #   Responders: 28/30 are (u+, v-); Non-responders: split evenly across (u+, v+), (u-, v+), (u-, v-)
        if is_r:
            if d_idx % 15 != 0:
                _emit_cells(d_idx, 0, 8)
                if d_idx % 10 != 0:
                    _emit_cells(d_idx, 2, 7)
            else:
                _emit_cells(d_idx, 1, 8)
        else:
            bucket = d_idx % 3
            if bucket == 0:
                _emit_cells(d_idx, 0, 8)
                _emit_cells(d_idx, 2, 7)
                _emit_cells(d_idx, 1, 8)
            elif bucket == 1:
                _emit_cells(d_idx, 1, 8)

        # Confounded intra-T_NK pair (3: M03, 4: M05): present in 90% of S2 regardless of R/NR!
        if st_id == "33958794" and (d_idx % 10 != 0):
            _emit_cells(d_idx, 3, 7)
            _emit_cells(d_idx, 4, 7)

    # Single-donor clonal spike in neighborhood 5 (M30 TLS_B) inside 1 responder donor
    _emit_cells(30, 5, 55)

    n_cells = len(pt_col)
    nhoods_mat = sp.csc_matrix(
        (np.ones(len(mat_rows), dtype=np.float64), (mat_rows, mat_cols)),
        shape=(n_cells, n_nhoods),
    )
    cells_df = pl.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_cells)],
            "patient_id": pt_col,
            "clinical_response": resp_col,
            "cell_type": ct_col,
            "UMAP1": np.asarray(umap_1_list, dtype=np.float64),
            "UMAP2": np.asarray(umap_2_list, dtype=np.float64),
        }
    )
    da_df = pl.DataFrame(
        {
            "Nhood": list(range(n_nhoods)),
            "logFC": [0.8, -0.9, 0.7, 1.4, 1.3, 2.0, 0.0],
            "FDR": [0.18, 0.22, 0.25, 0.02, 0.03, 0.40, 0.95],
            "Nhood_CellType": list(ct_names),
        }
    )
    meta_sub = pl.DataFrame(
        {
            "accession": ["SYNTH_CMH_01"] * n_nhoods,
            "nhood_idx": list(range(n_nhoods)),
            "meta_state": ["M04", "M17", "M04", "M03", "M05", "M30", "M10"],
            "meta_state_label": [
                "M04: GNLY+PRF1+ Eff. T",
                "M17: ISLR+COL5A2+ myCAF",
                "M04: GNLY+PRF1+ Eff. T",
                "M03: CD200+ELMO1+ Tfh T",
                "M05: CD96+TRAF3+ Mem. T",
                "M30: CR2+PAX5+ TLS B",
                "M10: FCN1+VCAN+ Mono/cDC",
            ],
        }
    )

    def _stub_loader(
        accession: str,
        config: PairwiseUmapCliConfig,
        all_meta42_df: pl.DataFrame,
    ) -> Result[CohortRawVizArtifacts, str]:
        _ = (config, all_meta42_df)
        return Success(
            CohortRawVizArtifacts(
                accession=accession,
                nhoods_mat=nhoods_mat,
                cells_df=cells_df,
                da_results_df=da_df,
                t0_quadrant_df=pl.DataFrame(),
                t0_continuous_df=pl.DataFrame(),
                t0_modules_df=pl.DataFrame(),
                meta42_sub_df=meta_sub,
            )
        )

    out_dir = tmp_path / "umap_cmh_out"
    cli_cfg = PairwiseUmapCliConfig(
        output_dir=out_dir,
        pair_mode=PairRepresentationMode.SUPER_META42,
        cohort_filter=Some(("SYNTH_CMH_01",)),
        max_chords=16,
        max_bg_cells=1500,
    )
    summaries = run_pairwise_umap_pipeline(cli_cfg=cli_cfg, loader=_stub_loader).unwrap()
    assert len(summaries) == 1
    s0 = summaries[0]
    # Duplicate (0, 1) and (2, 1) M04<->M17 sub-patches consolidate to strictly 1 discrete quadrant chord,
    # while intra-T_NK M03<->M05 and single-donor clonal M30 are completely excluded!
    assert s0.n_sig_quadrant == 1
    assert s0.n_displayed_chords >= 1
    assert "M04" in s0.lead_pair_desc and "M17" in s0.lead_pair_desc

    arcs_df = pl.read_parquet(
        out_dir / "SYNTH_CMH_01" / "SYNTH_CMH_01_super_meta42_bundled_arcs.parquet"
    )
    nodes_df = pl.read_parquet(
        out_dir / "SYNTH_CMH_01" / "SYNTH_CMH_01_super_meta42_nhood_umap_nodes.parquet"
    )
    idx_to_lbl = dict(zip(nodes_df["nhood_idx"].to_list(), nodes_df["state_label"].to_list()))
    unique_arcs = arcs_df.select(["arc_id", "nhood_u", "nhood_v"]).unique()
    for r in unique_arcs.to_dicts():
        lu = str(idx_to_lbl[int(r["nhood_u"])])
        lv = str(idx_to_lbl[int(r["nhood_v"])])
        assert "M30" not in lu and "M30" not in lv
        assert not ("M03" in lu and "M05" in lv)
        assert not ("M05" in lu and "M03" in lv)

    comp_svg = out_dir / "SYNTH_CMH_01" / "SYNTH_CMH_01_super_meta42_umap_composite.svg"
    comp_png = out_dir / "SYNTH_CMH_01" / "SYNTH_CMH_01_super_meta42_umap_composite.png"
    assert comp_svg.exists() and comp_png.exists()
    svg_text = comp_svg.read_text(encoding="utf-8")
    assert "feDropShadow" not in svg_text
    assert "Sub-Study CMH" in svg_text
    assert isinstance(PairwiseUmapBundle, type)
