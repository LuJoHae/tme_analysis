"""Unit and Hypothesis property-based tests for scripts/milopy_analysis/prefilter_and_deconv_bridge.py.

Verifies:
1. Exact outcome-blindness of `map_cohort_nhoods_to_meta42_outcome_blind` under arbitrary clinical
   response, univariate DA logFC, and bulk beta permutations via Hypothesis.
2. Exact binomial bounds & complement symmetry of `wilson_95_ci` via Hypothesis.
3. Hub-preserving Mantel node-permutation calibration under structured row-hub nulls vs. planted
   synergistic/antagonistic interaction block matrices.
4. Pure end-to-end execution of `run_full_bridge_pipeline` via dependency-injected Protocol stubs,
   myCAF gating regime recovery, and XML well-formedness / Nature Methods wireframe compliance of
   all 4 vector SVG figures via `vl-convert`.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Result, Success
from scipy import sparse as sp
import vl_convert as vlc

from scripts.milopy_analysis.prefilter_and_deconv_bridge import (
    AltairChart,
    BridgeConfig,
    CohortInputsBundle,
    FOCUSED_8_META_STATES,
    GlobalDeconvInputs,
    H5adValidationMetrics,
    NhoodMeta42Mapping,
    build_deconv_vs_pairwise_milo_chart,
    build_meta42_hub_and_pseudobulk_chart,
    build_prefilter_power_ablation_chart,
    build_singlecell_mycaf_gating_chart,
    map_cohort_nhoods_to_meta42_outcome_blind,
    parse_cli_args,
    run_full_bridge_pipeline,
    run_mantel_node_permutation_test,
    wilson_95_ci,
)


def _build_synthetic_bridge_fixture(
    accession: str = "SYNTH_WT_01",
    seed: int = 123,
) -> tuple[CohortInputsBundle, GlobalDeconvInputs]:
    """Constructs a whole-tumor synthetic cohort harboring M04/M29/M30 effectors and M17 myCAF."""
    rng = np.random.default_rng(seed)
    n_donors = 48  # 12 per quadrant of (Lymphoid, M17_myCAF)
    cells_per_donor = 20
    n_cells = n_donors * cells_per_donor
    n_nhoods = 10

    # States: L0->M04 (Eff T), L1->M29 (Plasma), L2->M30 (TLS B), L3->M17 (myCAF), L4->M32 (M2-TAM)
    lk_to_meta = {
        "L0": "M04",
        "L1": "M29",
        "L2": "M30",
        "L3": "M17",
        "L4": "M32",
    }
    meta_labels = {
        "M04": "M04: GNLY+PRF1+ Eff. T",
        "M05": "M05: CD96+TRAF3+ Mem. T",
        "M17": "M17: ISLR+COL5A2+ myCAF",
        "M20": "M20: IL7R+CD40LG+ Trm T",
        "M25": "M25: CTLA4+TIGIT+ Tex",
        "M29": "M29: JCHAIN+MZB1+ Plasma",
        "M30": "M30: CR2+PAX5+ TLS B",
        "M32": "M32: MRC1+CD163+ M2-TAM",
    }

    # Donors 0..11: Lymphoid+ / myCAF- (10/12 R = 83%)
    # Donors 12..23: Lymphoid+ / myCAF+ (2/12 R = 17% -> M17 stromal exclusion gate!)
    # Donors 24..35: Lymphoid- / myCAF- (3/12 R = 25%)
    # Donors 36..47: Lymphoid- / myCAF+ (1/12 R = 8%)
    donor_resps = (
        ["responder"] * 10 + ["non-responder"] * 2
        + ["responder"] * 2 + ["non-responder"] * 10
        + ["responder"] * 3 + ["non-responder"] * 9
        + ["responder"] * 1 + ["non-responder"] * 11
    )

    pt_ids: list[str] = []
    cell_resps: list[str] = []
    cell_types: list[str] = []
    rows: list[int] = []
    cols: list[int] = []
    pb_rows: list[dict[str, object]] = []

    # Pair of nhoods per state (0,1 -> L0/M04; 2,3 -> L1/M29; 4,5 -> L2/M30; 6,7 -> L3/M17; 8,9 -> L4/M32)
    nhood_cts = [
        "T cell",
        "T cell",
        "plasma cell",
        "plasma cell",
        "B cell",
        "B cell",
        "fibroblast",
        "fibroblast",
        "macrophage",
        "macrophage",
    ]

    for d_idx in range(n_donors):
        sid = f"DONOR_{d_idx:02d}"
        resp = donor_resps[d_idx]
        lym_on = d_idx < 24
        caf_on = (12 <= d_idx < 24) or (d_idx >= 36)
        tam_on = bool(rng.random() < (0.65 if caf_on else 0.30))

        state_counts = {
            "L0": 5 if lym_on else 0,
            "L1": 4 if (lym_on and (resp == "responder" or d_idx % 3 == 0)) else 0,
            "L2": 4 if (lym_on and (resp == "responder" or d_idx % 4 == 0)) else 0,
            "L3": 5 if caf_on else 0,
            "L4": 4 if tam_on else 2,
        }
        # Fill remaining cells up to cells_per_donor in L4
        used = sum(state_counts.values())
        if used < cells_per_donor:
            state_counts["L4"] += cells_per_donor - used

        base = d_idx * cells_per_donor
        cursor = 0
        for st_idx, lk in enumerate(("L0", "L1", "L2", "L3", "L4")):
            cnt = min(state_counts[lk], cells_per_donor - cursor)
            ct_name = nhood_cts[2 * st_idx]
            for c_off in range(cnt):
                cid = base + cursor + c_off
                pt_ids.append(sid)
                cell_resps.append(resp)
                cell_types.append(ct_name)
                if c_off < max(cnt - 1, 1):
                    rows.append(cid)
                    cols.append(2 * st_idx)
                if c_off >= 1 or cnt == 1:
                    rows.append(cid)
                    cols.append(2 * st_idx + 1)
            cursor += cnt

            tf = float(cnt) / float(cells_per_donor)
            inf_f = float(np.clip(tf + rng.normal(0.0, 0.015), 0.0, 1.0))
            pb_rows.append(
                {
                    "sample_id": sid,
                    "cell_state": lk,
                    "inferred_fraction": inf_f,
                    "n_cells": cnt,
                    "n_total": cells_per_donor,
                    "response": resp,
                    "treatment_status": "Pre-treatment",
                    "true_fraction": tf,
                }
            )

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, n_nhoods),
    )
    cells_df = pl.DataFrame(
        {
            "cell_id": [f"c_{i}" for i in range(n_cells)],
            "patient_id": pt_ids,
            "clinical_response": cell_resps,
            "cell_type": cell_types,
        }
    )
    da_df = pl.DataFrame(
        {
            "Nhood": list(range(n_nhoods)),
            "logFC": [0.8, 0.7, 1.4, 1.3, 1.2, 1.1, -1.5, -1.4, -0.4, -0.3],
            "FDR": [0.20, 0.25, 0.04, 0.05, 0.06, 0.08, 0.03, 0.04, 0.50, 0.60],
            "Nhood_CellType": nhood_cts,
        }
    )
    pb_df = pl.DataFrame(pb_rows)
    phi_df = pl.DataFrame(
        {
            "cluster": ["L0", "L1", "L2", "L3", "L4"],
            "GNLY": [0.02, 0.0001, 0.0001, 0.0, 0.0],
            "JCHAIN": [0.0, 0.03, 0.001, 0.0, 0.0],
            "MS4A1": [0.0, 0.001, 0.025, 0.0, 0.0],
            "COL5A2": [0.0, 0.0, 0.0, 0.03, 0.0002],
            "MRC1": [0.0, 0.0, 0.0, 0.0002, 0.02],
        }
    )

    cohort_bundle = CohortInputsBundle(
        accession=accession,
        nhoods_mat=nhoods_mat,
        cells_df=cells_df,
        da_results_df=da_df,
        pseudobulk_fractions_df=pb_df,
        reference_phi_df=phi_df,
    )

    conc_3b_df = pl.DataFrame(
        {
            "cohort": [accession] * 5,
            "cell_state": ["L0", "L1", "L2", "L3", "L4"],
            "meta_state": ["M04", "M29", "M30", "M17", "M32"],
            "milo_mean_logfc": [0.75, 1.35, 1.15, -1.45, -0.35],
            "deconv_beta": [0.5, 1.1, 0.9, -1.2, -0.3],
        }
    )

    # Build 8x8 interactions table covering FOCUSED_8_META_STATES
    int_rows: list[dict[str, object]] = []
    for sa in FOCUSED_8_META_STATES:
        for sb in FOCUSED_8_META_STATES:
            if sa == sb:
                ez = 0.0
                itype = "Self (Main Effect)"
            elif {sa, sb} in ({"M29", "M30"}, {"M04", "M29"}, {"M04", "M30"}):
                ez = 3.4
                itype = "Co-Activation Synergy (+/+)"
            elif "M17" in {sa, sb} and ({sa, sb} & {"M04", "M29", "M30"}):
                ez = -3.3
                itype = "Antagonistic Gate (+/−)"
            else:
                ez = -0.5
                itype = "Weak Background"
            int_rows.append(
                {
                    "state_a": sa,
                    "state_b": sb,
                    "state_a_label": meta_labels[sa],
                    "state_b_label": meta_labels[sb],
                    "state_a_beta": 1.0 if sa in ("M04", "M05", "M29", "M30") else -1.0,
                    "state_b_beta": 1.0 if sb in ("M04", "M05", "M29", "M30") else -1.0,
                    "interaction_type": itype,
                    "friedman_h": abs(ez) * 0.02,
                    "epistasis_z": ez,
                    "signed_interaction": ez * 0.02,
                }
            )
    int_df = pl.DataFrame(int_rows)

    hub_df = pl.DataFrame(
        {
            "cell_state": list(FOCUSED_8_META_STATES),
            "state_label": [meta_labels[s] for s in FOCUSED_8_META_STATES],
            "enrichment_group": [
                "Responder-Enriched (beta > 0)"
                if s in ("M04", "M05", "M29", "M30")
                else "Non-Responder-Enriched (beta < 0)"
                for s in FOCUSED_8_META_STATES
            ],
            "state_bulk_beta": [1.2, 0.8, -1.8, -0.4, -0.5, 1.9, 1.5, -0.9],
            "abs_bulk_beta": [1.2, 0.8, 1.8, 0.4, 0.5, 1.9, 1.5, 0.9],
            "state_milo_logfc": [1.0, 0.6, -1.6, -0.3, -0.4, 1.7, 1.3, -0.7],
            "mean_abs_epistasis_z": [2.4, 1.5, 2.6, 1.1, 1.2, 2.7, 2.5, 1.8],
            "mean_friedman_h": [0.06] * 8,
            "max_synergy_z": [3.4, 2.0, 0.5, 0.4, 0.4, 3.6, 3.4, 0.8],
            "min_gating_z": [-3.3, -1.8, -3.5, -1.2, -1.3, -3.3, -3.2, -2.5],
        }
    )
    deg_df = pl.DataFrame(
        {
            "cell_state": list(FOCUSED_8_META_STATES),
            "gene": ["GNLY", "CD96", "COL5A2", "IL7R", "TIGIT", "JCHAIN", "CR2", "MRC1"],
        }
    )

    global_deconv = GlobalDeconvInputs(
        concordance_3b_df=conc_3b_df,
        iatlas_all42_interactions_df=int_df,
        iatlas_focused8_interactions_df=int_df,
        iatlas_hub_summary_df=hub_df,
        iatlas_top_degs_df=deg_df,
    )
    return cohort_bundle, global_deconv


@settings(max_examples=20, deadline=None)
@given(perm_seed=st.integers(min_value=1, max_value=10000))
def test_prop_nhood_meta42_mapper_strictly_outcome_blind(perm_seed: int) -> None:
    """Proves neighborhood->L_k->Meta42 mapping is 100% invariant to clinical outcomes & DA logFCs."""
    bundle_orig, global_orig = _build_synthetic_bridge_fixture(seed=77)
    map_orig = map_cohort_nhoods_to_meta42_outcome_blind(
        inputs=bundle_orig,
        concordance_3b_df=global_orig.concordance_3b_df,
        hub_summary_df=global_orig.iatlas_hub_summary_df,
    ).unwrap()

    rng = np.random.default_rng(perm_seed)
    flipped_cell_resps = [
        "responder" if bool(rng.random() < 0.5) else "non-responder"
        for _ in range(bundle_orig.cells_df.height)
    ]
    flipped_pb_resps = [
        "responder" if bool(rng.random() < 0.5) else "non-responder"
        for _ in range(bundle_orig.pseudobulk_fractions_df.height)
    ]
    scrambled_da = bundle_orig.da_results_df.with_columns(
        pl.Series("logFC", rng.normal(0.0, 5.0, size=bundle_orig.da_results_df.height)),
        pl.Series("FDR", rng.uniform(0.0, 1.0, size=bundle_orig.da_results_df.height)),
    )
    bundle_perm = CohortInputsBundle(
        accession=bundle_orig.accession,
        nhoods_mat=bundle_orig.nhoods_mat,
        cells_df=bundle_orig.cells_df.with_columns(
            pl.Series("clinical_response", flipped_cell_resps)
        ),
        da_results_df=scrambled_da,
        pseudobulk_fractions_df=bundle_orig.pseudobulk_fractions_df.with_columns(
            pl.Series("response", flipped_pb_resps)
        ),
        reference_phi_df=bundle_orig.reference_phi_df,
    )

    map_perm = map_cohort_nhoods_to_meta42_outcome_blind(
        inputs=bundle_perm,
        concordance_3b_df=global_orig.concordance_3b_df,
        hub_summary_df=global_orig.iatlas_hub_summary_df,
    ).unwrap()

    assert map_orig.nhood_meta_states == map_perm.nhood_meta_states
    assert (
        map_orig.annotations_df["mapped_cell_state"].to_list()
        == map_perm.annotations_df["mapped_cell_state"].to_list()
    )
    np.testing.assert_allclose(
        map_orig.annotations_df["composite_affinity"].to_numpy(),
        map_perm.annotations_df["composite_affinity"].to_numpy(),
        atol=1e-12,
    )
    # Also verify exact recovery of planted L0..L4 -> M04, M29, M30, M17, M32
    assert map_orig.nhood_meta_states == (
        "M04",
        "M04",
        "M29",
        "M29",
        "M30",
        "M30",
        "M17",
        "M17",
        "M32",
        "M32",
    )


@settings(max_examples=50, deadline=None)
@given(
    n_total=st.integers(min_value=1, max_value=500),
    frac=st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
)
def test_prop_wilson_95_ci_invariants(n_total: int, frac: float) -> None:
    """Proves Wilson 95% CI obeys [0, 1] bounds, brackets p_hat, and satisfies exact complement symmetry."""
    k = int(round(frac * n_total))
    p_hat = float(k) / float(n_total)
    lo, hi = wilson_95_ci(k, n_total)
    assert 0.0 <= lo <= p_hat + 1e-12
    assert p_hat - 1e-12 <= hi <= 1.0

    lo_comp, hi_comp = wilson_95_ci(n_total - k, n_total)
    np.testing.assert_allclose(lo, 1.0 - hi_comp, atol=1e-12)
    np.testing.assert_allclose(hi, 1.0 - lo_comp, atol=1e-12)


def test_mantel_node_permutation_null_and_planted_block_power() -> None:
    """Verifies Mantel node-permutation test detects planted block concordance and stays calm under H_0."""
    states = [f"M{i:02d}" for i in range(1, 13)]
    rng = np.random.default_rng(999)

    # Planted block matrix: M01..M04 synergize (+3), M05..M08 gate M01..M04 (-3)
    pair_x: dict[tuple[str, str], float] = {}
    pair_y_planted: dict[tuple[str, str], float] = {}
    pair_y_null: dict[tuple[str, str], float] = {}

    hub_effects = rng.normal(0.0, 1.5, size=len(states))
    for i, sa in enumerate(states):
        for j in range(i + 1, len(states)):
            sb = states[j]
            base = 3.0 if (i < 4 and j < 4) else (-3.0 if (i < 4 and 4 <= j < 8) else 0.0)
            pair_x[(sa, sb)] = base + float(rng.normal(0.0, 0.25))
            pair_y_planted[(sa, sb)] = 0.8 * base + float(rng.normal(0.0, 0.35))
            # Independent additive hub-structured null matrix
            pair_y_null[(sa, sb)] = float(hub_effects[i] + hub_effects[j] + rng.normal(0.0, 0.5))

    res_planted = run_mantel_node_permutation_test(
        states=states,
        pair_to_x=pair_x,
        pair_to_y=pair_y_planted,
        panel_name="Planted",
        n_perms=1000,
        seed=42,
    )
    assert res_planted.pearson_r > 0.85
    assert res_planted.pvalue_two_sided < 0.005
    assert res_planted.directional_concordance_strong > 0.95

    res_null = run_mantel_node_permutation_test(
        states=states,
        pair_to_x=pair_x,
        pair_to_y=pair_y_null,
        panel_name="Null",
        n_perms=1000,
        seed=42,
    )
    assert res_null.pvalue_two_sided > 0.05


def test_end_to_end_injected_bridge_pipeline_and_svg_xml_validity() -> None:
    """Verifies pure end-to-end execution via Protocol stubs and XML validity of all 4 SVG figures."""
    cohort_bundle, global_deconv = _build_synthetic_bridge_fixture()
    exported_svgs: list[tuple[str, str]] = []

    def stub_cohort_loader(
        accession: str,
        milo_dir: Path,
        deconv_3a_dir: Path,
    ) -> Result[CohortInputsBundle, str]:
        _ = (accession, milo_dir, deconv_3a_dir)
        return Success(cohort_bundle)

    def stub_global_loader(
        deconv_3b_path: Path,
        iatlas_dir: Path,
    ) -> Result[GlobalDeconvInputs, str]:
        _ = (deconv_3b_path, iatlas_dir)
        return Success(global_deconv)

    def stub_h5ad_validator(
        local_h5ad_dir: Path,
        cohort_inputs: CohortInputsBundle,
        mapping: NhoodMeta42Mapping,
        lk_to_meta: Mapping[str, str],
    ) -> Maybe[H5adValidationMetrics]:
        _ = (local_h5ad_dir, mapping, lk_to_meta)
        return Some(
            H5adValidationMetrics(
                accession=cohort_inputs.accession,
                n_nhoods=10,
                n_shared_genes=5,
                exact_lk_agreement=1.0,
                exact_meta42_agreement=1.0,
                exact_compartment_agreement=1.0,
            )
        )

    def stub_exporter(chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
        svg_str = vlc.vegalite_to_svg(chart.to_json())
        exported_svgs.append((base_path.name, svg_str))
        return Success((base_path.with_suffix(".svg"), base_path.with_suffix(".png")))

    cfg = parse_cli_args(
        ["--cohorts", "SYNTH_WT_01", "--n-mantel-perms", "200", "--fdr-threshold", "0.10"]
    )
    res = run_full_bridge_pipeline(
        config=cfg,
        cohort_loader=stub_cohort_loader,
        global_loader=stub_global_loader,
        h5ad_validator=stub_h5ad_validator,
        exporter=stub_exporter,
        persist_outputs=False,
    )
    assert isinstance(res, Success)
    out = res.unwrap()

    assert out.ablation_benchmark_df.height == 6
    assert out.nhood_annotations_df.height == 10
    assert isinstance(out.h5ad_validation, Some)

    # Verify M17 myCAF gating regime recovery: Lymphoid+ / myCAF- beats Lymphoid+ / myCAF+
    reg_map = {
        str(r["gating_regime"]): float(r["responder_pct"])
        for r in out.mycaf_gating_regimes_df.to_dicts()
    }
    assert reg_map["Lymphoid+ / myCAF-"] > reg_map["Lymphoid+ / myCAF+"] + 30.0

    # Compile and verify XML & Nature Methods wireframe standard on all 4 charts
    charts = (
        build_prefilter_power_ablation_chart(out.ablation_benchmark_df).unwrap(),
        build_deconv_vs_pairwise_milo_chart(
            out.meta42_concordance_df, out.mantel_whole_tumor, out.mantel_focused8
        ).unwrap(),
        build_singlecell_mycaf_gating_chart(
            out.mycaf_gating_regimes_df, out.mycaf_cohort_forest_df
        ).unwrap(),
        build_meta42_hub_and_pseudobulk_chart(
            out.hub_concordance_df, out.within_cohort_pseudobulk_df
        ).unwrap(),
    )
    for ch in charts:
        svg_xml = vlc.vegalite_to_svg(ch.to_json())
        root = ET.fromstring(svg_xml)
        assert root.tag.endswith("svg")
        assert "feDropShadow" not in svg_xml
        assert len(svg_xml) > 1000
