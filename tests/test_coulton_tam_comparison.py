"""Property-based (`hypothesis`), unit, SVG, and AST exclusivity tests for `coulton_tam_comparison`."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Final
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Some
from returns.pipeline import is_successful
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from scripts.coulton_tam_comparison.compare_meta42_coulton_tam import (
    execute_comparison_pipeline,
    linearize_and_standardize_bulk_matrix,
    parse_cli_config,
)
from scripts.coulton_tam_comparison.core import (
    accumulate_coulton_batch,
    benjamini_hochberg_fdr,
    compute_centered_log2_specificity,
    compute_mann_whitney_auroc,
    compute_stratum_similarity_df,
    compute_ward_biclustering,
    create_empty_coulton_accumulator,
    deduplicate_gene_columns,
    deserialize_coulton_centroids_from_polars,
    finalize_coulton_centroids,
    fit_penalized_logistic_newton,
    inverse_variance_meta_analysis,
    row_cosine_similarity,
    row_spearman_correlation,
    serialize_coulton_centroids_to_polars,
    strip_coulton_cluster_number,
)
from scripts.coulton_tam_comparison.models import (
    AltairChart,
    AssaySource,
    BulkIcbCohortBundle,
    COULTON_ORDERED_CLUSTERS,
    ComparisonConfig,
    CoultonCentroidBundle,
    META42_BIO_LABELS,
    Meta42ReferenceBundle,
    SignificanceTier,
    TissueStratum,
)
from scripts.coulton_tam_comparison.plotting import (
    build_simplified_bicluster_typst_caption,
    build_simplified_biclustered_myeloid_coulton_chart,
    classify_significance_tier,
    format_significance_stars,
)
from scripts.coulton_tam_comparison.render_figures import inject_inkscape_layers_pure
from tests.test_tme_datasets_api_exclusivity import audit_python_script

REPO_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
COMPARISON_PKG_DIR: Final[Path] = REPO_ROOT / "scripts" / "coulton_tam_comparison"


@given(
    n_cells=st.integers(min_value=8, max_value=45),
    n_genes=st.integers(min_value=4, max_value=14),
    batch_size=st.integers(min_value=2, max_value=13),
    seed=st.integers(min_value=0, max_value=10_000),
)
@settings(max_examples=30)
def test_sparse_csr_batch_accumulator_matches_dense_exact(
    n_cells: int,
    n_genes: int,
    batch_size: int,
    seed: int,
) -> None:
    """Property: multi-batch sparse CSR reduction equals single-shot exact dense `expm1` means."""
    rng = np.random.default_rng(seed)
    k_clusters = 4
    cluster_labels = tuple(f"C{i}" for i in range(k_clusters))
    gene_symbols = tuple(f"GENE{g}" for g in range(n_genes))

    raw_dense = np.where(
        rng.random((n_cells, n_genes)) < 0.45,
        rng.uniform(0.05, 4.0, size=(n_cells, n_genes)),
        0.0,
    ).astype(np.float64)
    cluster_codes = rng.integers(0, k_clusters, size=n_cells, dtype=np.int64)
    is_tumor = rng.random(n_cells) < 0.65

    acc = create_empty_coulton_accumulator(
        assay=AssaySource.RNA_LOG1P,
        cluster_labels=cluster_labels,
        gene_symbols=gene_symbols,
    )
    csr_full = sp.csr_matrix(raw_dense)
    for start in range(0, n_cells, batch_size):
        stop = min(start + batch_size, n_cells)
        step_res = accumulate_coulton_batch(
            acc=acc,
            batch_csr=csr_full[start:stop],
            cluster_codes=cluster_codes[start:stop],
            is_tumor_mask=is_tumor[start:stop],
        )
        assert is_successful(step_res)
        acc = step_res.unwrap()

    bundle = finalize_coulton_centroids(acc)
    lin_dense = np.expm1(raw_dense)
    det_dense = (raw_dense > 0.0).astype(np.float64)

    for k in range(k_clusters):
        all_mask = cluster_codes == k
        tum_mask = all_mask & is_tumor
        assert int(bundle.all_cell_counts[k]) == int(np.sum(all_mask))
        assert int(bundle.tumor_cell_counts[k]) == int(np.sum(tum_mask))

        if int(np.sum(all_mask)) > 0:
            np.testing.assert_allclose(
                bundle.all_mean_linear[k],
                lin_dense[all_mask].mean(axis=0),
                atol=1e-10,
            )
            np.testing.assert_allclose(
                bundle.all_mean_log[k],
                raw_dense[all_mask].mean(axis=0),
                atol=1e-10,
            )
            np.testing.assert_allclose(
                bundle.all_pct_expressing[k],
                det_dense[all_mask].mean(axis=0) * 100.0,
                atol=1e-10,
            )
        if int(np.sum(tum_mask)) > 0:
            np.testing.assert_allclose(
                bundle.tumor_mean_linear[k],
                lin_dense[tum_mask].mean(axis=0),
                atol=1e-10,
            )


@given(
    n_rows=st.integers(min_value=4, max_value=12),
    n_genes=st.integers(min_value=6, max_value=20),
    scale_factor=st.floats(min_value=0.1, max_value=25.0),
    seed=st.integers(min_value=0, max_value=10_000),
)
@settings(max_examples=30)
def test_centered_log2_specificity_invariants(
    n_rows: int,
    n_genes: int,
    scale_factor: float,
    seed: int,
) -> None:
    """Property: reference-centered specificity has zero column mean over ref rows & scale invariance."""
    rng = np.random.default_rng(seed)
    mat = rng.uniform(0.01, 10.0, size=(n_rows, n_genes)).astype(np.float64)
    ref_idx = list(range(max(2, n_rows // 2)))

    spec_a = compute_centered_log2_specificity(mat, ref_idx)
    spec_b = compute_centered_log2_specificity(mat * scale_factor, ref_idx)

    np.testing.assert_allclose(spec_a[ref_idx, :].mean(axis=0), 0.0, atol=1e-10)
    np.testing.assert_allclose(spec_a, spec_b, atol=1e-10)

    cos = row_cosine_similarity(spec_a, spec_b)
    assert np.all(cos >= -1.0) and np.all(cos <= 1.0)
    spear = row_spearman_correlation(spec_a, spec_b)
    assert np.all(spear >= -1.0) and np.all(spear <= 1.0)


def test_within_compartment_centering_recovers_substates_despite_shared_baseline() -> None:
    """Within-compartment specificity recovers 1-to-1 TAM substate identity under massive pan-mac baseline."""
    # 3 states sharing 90% pan-macrophage housekeeping density + small substate-specific markers
    base = np.array([50.0, 45.0, 40.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], dtype=np.float64)
    state_ifng = base + np.array([0.0, 0.0, 0.0, 4.0, 3.5, 0.0, 0.0, 0.0, 0.0])
    state_m2 = base + np.array([0.0, 0.0, 0.0, 0.0, 0.0, 4.2, 3.8, 0.0, 0.0])
    state_ecm = base + np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 4.0])

    meta_mat = np.vstack([state_ifng, state_m2, state_ecm])
    coulton_mat = np.vstack([state_ifng * 1.1, state_m2 * 0.9, state_ecm * 1.2])

    spec_m = compute_centered_log2_specificity(meta_mat, [0, 1, 2])
    spec_c = compute_centered_log2_specificity(coulton_mat, [0, 1, 2])
    cos = row_cosine_similarity(spec_m, spec_c)

    assert np.argmax(cos[0]) == 0
    assert np.argmax(cos[1]) == 1
    assert np.argmax(cos[2]) == 2
    assert float(np.min(np.diag(cos))) > 0.98


def test_marker_anchored_hybrid_similarity_recovers_canonical_m14_match() -> None:
    """Stage 3c canonical 4-DEG anchor recovers M14 -> 8_IFNGMac even under flat HVG baseline."""
    genes = tuple(
        ["CXCL9", "STAT1", "GBP1", "IDO1", "MRC1", "CD163", "FABP4", "C1QA"]
        + [f"G{i}" for i in range(1, 55)]
    )
    n_genes = len(genes)
    meta_ids = ("M14", "M32", "M42")
    flat_phi = np.full((3, n_genes), 1.0 / n_genes, dtype=np.float64)
    phi_df = pl.DataFrame(
        {"cluster": list(meta_ids), **{g: flat_phi[:, idx] for idx, g in enumerate(genes)}}
    )
    degs_df = pl.DataFrame(
        {
            "cell_state": ["M14", "M14", "M14", "M14", "M32", "M32", "M32", "M32"],
            "gene": ["CXCL9", "STAT1", "GBP1", "IDO1", "MRC1", "CD163", "FABP4", "C1QA"],
            "specificity_log2fc": [3.25, 2.52, 1.82, 1.22, 3.25, 2.52, 1.82, 1.22],
        }
    )

    coulton_lin = np.full((len(COULTON_ORDERED_CLUSTERS), n_genes), 0.10, dtype=np.float64)
    coulton_lin[COULTON_ORDERED_CLUSTERS.index("8_IFNGMac"), :4] = 6.0
    coulton_lin[COULTON_ORDERED_CLUSTERS.index("21_HemeMac"), 4:8] = 6.0
    counts = np.array(
        [2 if lbl == "23_NA" else 300 for lbl in COULTON_ORDERED_CLUSTERS],
        dtype=np.int64,
    )
    bundle = CoultonCentroidBundle(
        assay=AssaySource.RNA_LOG1P,
        cluster_labels=COULTON_ORDERED_CLUSTERS,
        gene_symbols=genes,
        all_cell_counts=counts,
        tumor_cell_counts=counts,
        all_mean_linear=coulton_lin,
        tumor_mean_linear=coulton_lin,
        all_mean_log=np.log1p(coulton_lin),
        tumor_mean_log=np.log1p(coulton_lin),
        all_pct_expressing=np.full_like(coulton_lin, 30.0),
        tumor_pct_expressing=np.full_like(coulton_lin, 30.0),
    )

    res = compute_stratum_similarity_df(
        phi_meta_df=phi_df,
        coulton_bundle=bundle,
        stratum=TissueStratum.TUMOR_ONLY,
        top_degs_df=Some(degs_df),
    )
    assert is_successful(res)
    sim_df = res.unwrap()
    m14_top1 = sim_df.filter(
        (pl.col("meta_state") == "M14") & (pl.col("rank_among_true_tams") == 1)
    ).row(0, named=True)
    m32_top1 = sim_df.filter(
        (pl.col("meta_state") == "M32") & (pl.col("rank_among_true_tams") == 1)
    ).row(0, named=True)

    assert m14_top1["coulton_cluster"] == "8_IFNGMac"
    assert float(m14_top1["canonical_deg_z"]) > 2.0
    assert "CXCL9" in str(m14_top1["shared_signature_genes"])
    assert m32_top1["coulton_cluster"] == "21_HemeMac"


@given(
    shift=st.floats(min_value=0.8, max_value=2.5),
    seed=st.integers(min_value=0, max_value=5_000),
)
@settings(max_examples=25)
def test_logistic_newton_and_meta_analysis_sign_symmetry(shift: float, seed: int) -> None:
    """Property: Newton-Raphson logistic slope & AUROC flip symmetrically under score negation."""
    rng = np.random.default_rng(seed)
    n_Half = 30
    y = np.concatenate([np.zeros(n_Half), np.ones(n_Half)]).astype(np.float64)
    noise = rng.normal(0.0, 0.6, size=2 * n_Half)
    scores_pos = noise + y * shift
    scores_neg = -scores_pos

    intercept = np.ones(2 * n_Half, dtype=np.float64)
    b_pos, se_pos = fit_penalized_logistic_newton(np.column_stack([intercept, scores_pos]), y)
    b_neg, se_neg = fit_penalized_logistic_newton(np.column_stack([intercept, scores_neg]), y)

    assert float(b_pos[1]) > 0.0
    assert abs(float(b_pos[1]) + float(b_neg[1])) < 1e-5
    np.testing.assert_allclose(se_pos[1], se_neg[1], atol=1e-5)

    auc_pos = compute_mann_whitney_auroc(scores_pos, y)
    auc_neg = compute_mann_whitney_auroc(scores_neg, y)
    assert auc_pos > 0.5
    assert abs((auc_pos + auc_neg) - 1.0) < 1e-10

    m_beta, m_se, m_p = inverse_variance_meta_analysis(
        np.array([b_pos[1], b_pos[1]], dtype=np.float64),
        np.array([se_pos[1], se_pos[1]], dtype=np.float64),
    )
    assert m_beta > 0.0
    assert m_se < float(se_pos[1])
    assert 0.0 <= m_p <= 1.0


@given(
    pvals=st.lists(
        st.floats(min_value=1e-12, max_value=1.0, allow_nan=False, allow_infinity=False),
        min_size=1,
        max_size=30,
    )
)
@settings(max_examples=30)
def test_benjamini_hochberg_monotonicity_and_bounds(pvals: list[float]) -> None:
    """Property: BH FDR values stay in [0, 1], dominate raw p-values, and preserve rank order."""
    p_arr = np.asarray(pvals, dtype=np.float64)
    q_arr = benjamini_hochberg_fdr(p_arr)
    assert np.all(q_arr >= 0.0) and np.all(q_arr <= 1.0)
    assert np.all(q_arr + 1e-12 >= p_arr)
    order = np.argsort(p_arr)
    assert np.all(np.diff(q_arr[order]) >= -1e-12)


def test_polars_centroid_roundtrip_and_gene_dedup() -> None:
    """Verify duplicate HGNC aggregation and lossless Polars DataFrame roundtrip."""
    acc = create_empty_coulton_accumulator(
        assay=AssaySource.RNA_LOG1P,
        cluster_labels=("0_AlvMac", "8_IFNGMac"),
        gene_symbols=("CXCL9", "CXCL9", "MRC1"),
    )
    csr = sp.csr_matrix(np.array([[1.0, 0.5, 0.0], [0.0, 2.0, 1.5]], dtype=np.float64))
    res = accumulate_coulton_batch(
        acc=acc,
        batch_csr=csr,
        cluster_codes=np.array([0, 1], dtype=np.int64),
        is_tumor_mask=np.array([True, True], dtype=np.bool_),
    )
    assert is_successful(res)
    dedup = deduplicate_gene_columns(
        raw_symbols=("CXCL9", "CXCL9", "MRC1"),
        bundle=finalize_coulton_centroids(res.unwrap()),
    )
    assert dedup.gene_symbols == ("CXCL9", "MRC1")
    df = serialize_coulton_centroids_to_polars((dedup,))
    reloaded = deserialize_coulton_centroids_from_polars(df)
    assert is_successful(reloaded)
    re_b = reloaded.unwrap()[0]
    assert re_b.gene_symbols == ("CXCL9", "MRC1")
    np.testing.assert_allclose(re_b.tumor_mean_linear, dedup.tumor_mean_linear, rtol=1e-5, atol=1e-6)


def test_inject_inkscape_layers_pure_wireframe_compliance() -> None:
    """Verify SVG post-processor injects Inkscape layers, #FFFFFF canvas, rx=0, and zero shadows."""
    minimal_svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="400" height="300">'
        '<defs><clipPath id="c1"/></defs>'
        '<g class="mark-rect"><rect x="10" y="10" width="50" height="20" fill="#0072B2"/></g>'
        "</svg>"
    )
    res = inject_inkscape_layers_pure(minimal_svg, figure_id="test_fig")
    assert is_successful(res)
    layered_svg = res.unwrap()
    assert "xmlns:inkscape=" in layered_svg
    assert "xmlns:sodipodi=" in layered_svg
    assert "feDropShadow" not in layered_svg

    root = ET.fromstring(layered_svg)
    ink_ns = "http://www.inkscape.org/namespaces/inkscape"
    svg_ns = "http://www.w3.org/2000/svg"
    layers = [
        el.attrib.get(f"{{{ink_ns}}}label", "")
        for el in root.findall(f"{{{svg_ns}}}g")
        if el.attrib.get(f"{{{ink_ns}}}groupmode") == "layer"
    ]
    assert layers == ["00_Canvas_Background", "01_Data_Panels_And_Axes"]


def test_end_to_end_pipeline_with_synthetic_injected_dependencies(tmp_path: Path) -> None:
    """Verify full `execute_comparison_pipeline` via pure Protocol dependency injection."""
    rng = np.random.default_rng(42)
    genes = tuple(
        [
            "CXCL9", "CXCL10", "STAT1", "GBP1", "MRC1", "CD163", "FOLR2", "SELENOP",
            "SULT1B1", "RGS2", "COL1A1", "COL1A2", "SPARC", "S100A8", "S100A9", "VCAN",
            "MKI67", "TOP2A", "CD68", "CD14", "CSF1R", "TYROBP", "C1QA", "C1QB", "C1QC",
        ]
        + [f"GENE{i}" for i in range(1, 55)]
    )
    n_genes = len(genes)
    meta_ids = tuple(sorted(META42_BIO_LABELS.keys()))

    phi_mat = rng.uniform(0.001, 0.02, size=(len(meta_ids), n_genes))
    # Give M14 strong CXCL9/CXCL10/STAT1, M32 strong MRC1/CD163/FOLR2
    phi_mat[meta_ids.index("M14"), :4] = 0.20
    phi_mat[meta_ids.index("M32"), 4:8] = 0.20
    phi_mat = phi_mat / phi_mat.sum(axis=1, keepdims=True)

    phi_df = pl.DataFrame(
        {"cluster": list(meta_ids), **{g: phi_mat[:, idx] for idx, g in enumerate(genes)}}
    )
    eff_df = pl.DataFrame(
        {
            "reference_id": ["Consensus-Meta42"] * len(meta_ids),
            "bulk_cohort": ["Pan-iAtlas-Meta"] * len(meta_ids),
            "bulk_short_name": ["Pan-iAtlas Meta (N=1,015)"] * len(meta_ids),
            "cell_state": list(meta_ids),
            "milo_mean_logfc": [0.8 if m == "M14" else (-0.4 if m == "M32" else 0.05) for m in meta_ids],
            "milo_se": [0.15] * len(meta_ids),
            "bulk_deconv_beta": [0.65 if m == "M14" else (-0.45 if m == "M32" else 0.02) for m in meta_ids],
            "bulk_deconv_se": [0.12] * len(meta_ids),
            "bulk_deconv_pval": [1e-5 if m in {"M14", "M32"} else 0.6 for m in meta_ids],
            "bulk_multivariate_beta": [0.3] * len(meta_ids),
            "delta_clr": [0.2] * len(meta_ids),
            "quadrant": ["Concordant Responder (+/+)" if m == "M14" else "Concordant Non-Responder (-/-)" for m in meta_ids],
        }
    )
    hub_df = pl.DataFrame(
        {
            "cell_state": list(meta_ids),
            "enrichment_group": ["Responder-Enriched" if m == "M14" else "Non-Responder-Enriched" for m in meta_ids],
            "mean_abs_epistasis_z": [1.5] * len(meta_ids),
            "max_synergy_z": [2.2] * len(meta_ids),
            "min_gating_z": [-2.1] * len(meta_ids),
        }
    )

    coulton_lin = rng.uniform(0.01, 0.2, size=(len(COULTON_ORDERED_CLUSTERS), n_genes))
    coulton_lin[COULTON_ORDERED_CLUSTERS.index("8_IFNGMac"), :4] = 5.0
    coulton_lin[COULTON_ORDERED_CLUSTERS.index("17_IFNMac3"), :3] = 2.5
    coulton_lin[COULTON_ORDERED_CLUSTERS.index("1_MetM2Mac"), 4:8] = 5.0
    coulton_pct = np.full_like(coulton_lin, 25.0)
    counts = np.array(
        [2 if lbl == "23_NA" else 500 for lbl in COULTON_ORDERED_CLUSTERS],
        dtype=np.int64,
    )

    def fake_streamer(assay: AssaySource, batch_size: int) -> Result[CoultonCentroidBundle, str]:
        assert batch_size >= 1000
        return Success(
            CoultonCentroidBundle(
                assay=assay,
                cluster_labels=COULTON_ORDERED_CLUSTERS,
                gene_symbols=genes,
                all_cell_counts=counts,
                tumor_cell_counts=counts,
                all_mean_linear=coulton_lin,
                tumor_mean_linear=coulton_lin,
                all_mean_log=np.log1p(coulton_lin),
                tumor_mean_log=np.log1p(coulton_lin),
                all_pct_expressing=coulton_pct,
                tumor_pct_expressing=coulton_pct,
            )
        )

    def fake_meta_loader(_cfg: ComparisonConfig) -> Result[Meta42ReferenceBundle, str]:
        return Success(
            Meta42ReferenceBundle(
                phi_meta_df=phi_df,
                concordance_3b_df=pl.DataFrame({"state": list(meta_ids)}),
                hub_summary_df=hub_df,
                state_effects_df=eff_df,
                top_degs_df=pl.DataFrame({"metastate_id": ["M14"], "gene_symbol": ["CXCL9"]}),
            )
        )

    def fake_bulk_loader() -> Result[tuple[BulkIcbCohortBundle, ...], str]:
        raw_bulk = rng.uniform(1.0, 12.0, size=(40, n_genes))
        z_mat = linearize_and_standardize_bulk_matrix(raw_bulk)
        y_vec = np.array([1.0 if i >= 20 else 0.0 for i in range(40)], dtype=np.float64)
        return Success(
            (
                BulkIcbCohortBundle(
                    cohort_id="Synthetic-ICB",
                    cancer_type="Melanoma",
                    short_name="Synthetic (Melanoma)",
                    n_samples=40,
                    n_responders=20,
                    response_binary=y_vec,
                    genes=genes,
                    gene_to_idx={g: idx for idx, g in enumerate(genes)},
                    log2_z_mat=z_mat,
                ),
            )
        )

    def fake_exporter(_chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
        svg_p = base_path.with_suffix(".svg")
        png_p = base_path.with_suffix(".png")
        svg_p.write_text("<svg/>", encoding="utf-8")
        png_p.write_bytes(b"PNG")
        return Success((svg_p, png_p))

    cfg = ComparisonConfig(
        output_dir=tmp_path / "out",
        include_sct_sensitivity=False,
        reuse_cached_centroids=False,
        render_plots=True,
    )
    res = execute_comparison_pipeline(
        config=cfg,
        streamer=fake_streamer,
        meta42_loader=fake_meta_loader,
        bulk_loader=fake_bulk_loader,
        exporter=fake_exporter,
    )
    match res:
        case Failure(err):
            raise AssertionError(f"Synthetic pipeline failed: {err}")
        case Success((outputs, written)):
            assert len(written) == 7 + 1 + 8 + 1
            cap_file = tmp_path / "out" / "coulton_tam_simplified_bicluster_main_caption.typ"
            assert cap_file.is_file()
            cap_text = cap_file.read_text(encoding="utf-8")
            assert "#let caption_coulton_simplified_bicluster_main = [" in cap_text
            assert "@coultonUsingPancancerAtlas2024" in cap_text
            assert "`****`" in cap_text

            chart_res = build_simplified_biclustered_myeloid_coulton_chart(
                full_sim_df=outputs.full_similarity_df,
                nearest_df=outputs.nearest_substates_df,
                coulton_ici_df=outputs.coulton_ici_outcomes_df,
            )
            assert is_successful(chart_res)
            import vl_convert as vlc

            svg_with_letters = vlc.vegalite_to_svg(chart_res.unwrap().to_json())
            assert "role-title" not in svg_with_letters
            assert ">a<" in svg_with_letters
            assert ">b<" in svg_with_letters
            assert ">c<" in svg_with_letters
            assert ">#1<" in svg_with_letters
            assert ">#2<" not in svg_with_letters
            assert ">#3<" not in svg_with_letters

            chart_no_letters = build_simplified_biclustered_myeloid_coulton_chart(
                full_sim_df=outputs.full_similarity_df,
                nearest_df=outputs.nearest_substates_df,
                coulton_ici_df=outputs.coulton_ici_outcomes_df,
                show_panel_letters=False,
            )
            assert is_successful(chart_no_letters)
            svg_no_letters = vlc.vegalite_to_svg(chart_no_letters.unwrap().to_json())
            assert "role-title" not in svg_no_letters
            assert ">a<" not in svg_no_letters

            cap_pure_res = build_simplified_bicluster_typst_caption(
                full_sim_df=outputs.full_similarity_df,
                nearest_df=outputs.nearest_substates_df,
                coulton_ici_df=outputs.coulton_ici_outcomes_df,
            )
            assert is_successful(cap_pure_res)

            m14_top1 = outputs.nearest_substates_df.filter(
                (pl.col("meta_state") == "M14") & (pl.col("rank_among_true_tams") == 1)
            )["coulton_cluster"][0]
            m32_top1 = outputs.nearest_substates_df.filter(
                (pl.col("meta_state") == "M32") & (pl.col("rank_among_true_tams") == 1)
            )["coulton_cluster"][0]
            assert m14_top1 == "8_IFNGMac"
            assert m32_top1 == "1_MetM2Mac"


def test_coulton_comparison_ast_purity_and_tme_datasets_exclusivity() -> None:
    """Enforce tme_datasets API exclusivity and zero `None`/`del`/`global`/`nonlocal` in package."""
    py_files = sorted(COMPARISON_PKG_DIR.glob("*.py"))
    assert len(py_files) == 6

    for script_path in py_files:
        audit = audit_python_script(script_path)
        assert len(audit.forbidden_calls) == 0, (
            f"{script_path.name} contains forbidden direct dataset calls: {audit.forbidden_calls}"
        )
        assert len(audit.h5ad_globs) == 0, (
            f"{script_path.name} globs .h5ad directly: {audit.h5ad_globs}"
        )

        tree = ast.parse(script_path.read_text(encoding="utf-8"), filename=str(script_path))
        for node in ast.walk(tree):
            assert not isinstance(node, (ast.Delete, ast.Global, ast.Nonlocal)), (
                f"Forbidden imperative statement {type(node).__name__} in {script_path.name}@L{getattr(node, 'lineno', 0)}"
            )
            if isinstance(node, ast.Constant) and node.value is None:
                raise AssertionError(
                    f"Literal `None` is banned by code-style-guide.md; found in {script_path.name}@L{node.lineno}"
                )


def test_cli_manuscript_out_dir_and_snakemake_makefile_wiring() -> None:
    """Verify CLI `--manuscript-out-dir` flag and Snakemake + Makefile `olm` rule registration."""
    cfg = parse_cli_config(
        [
            "--manuscript-out-dir",
            "/storage/halu/data-test/output/manuscript_pipeline",
            "--output-dir",
            "/storage/halu/data-test/results/coulton_tam_comparison",
            "--skip-plots",
        ]
    )
    assert cfg.phi_meta_path == Path(
        "/storage/halu/data-test/output/manuscript_pipeline/"
        "3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet"
    )
    assert cfg.state_effects_path == Path(
        "/storage/halu/data-test/output/manuscript_pipeline/"
        "3c_iatlas_instaprism_milopy_benchmark/milopy_iatlas_state_effects.parquet"
    )
    assert cfg.render_plots is False

    smk_path = REPO_ROOT / "workflow" / "rules" / "coulton_tam_comparison.smk"
    olm_prof_path = REPO_ROOT / "workflow" / "profiles" / "olm" / "config.v8+.yaml"
    snakefile_text = (REPO_ROOT / "Snakefile").read_text(encoding="utf-8")
    makefile_text = (REPO_ROOT / "Makefile").read_text(encoding="utf-8")

    assert smk_path.is_file()
    assert olm_prof_path.is_file()
    smk_text = smk_path.read_text(encoding="utf-8")
    assert "rule all_coulton_tam_comparison:" in smk_text
    assert "rule compare_meta42_coulton_tam_tables:" in smk_text
    assert "rule render_coulton_tam_comparison_figures:" in smk_text
    assert 'include: "workflow/rules/coulton_tam_comparison.smk"' in snakefile_text
    assert "coulton-tam-comparison-olm:" in makefile_text
    assert "SNAKEMAKE_OLM" in makefile_text


@given(
    p1=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    p2=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
)
@settings(max_examples=50)
def test_significance_stars_monotonicity_and_boundaries(p1: float, p2: float) -> None:
    """Property: Smaller p-values never produce fewer significance asterisks."""
    low_p, high_p = (p1, p2) if p1 <= p2 else (p2, p1)
    s_low = format_significance_stars(low_p)
    s_high = format_significance_stars(high_p)
    assert set(s_low) <= {"*"}
    assert set(s_high) <= {"*"}
    assert len(s_low) >= len(s_high)
    assert len(s_low) <= 4


def test_significance_tier_exact_cutoffs() -> None:
    """Verify exact 4-tier significance boundaries (* <0.05, ** <0.01, *** <1e-3, **** <1e-4) and NaN."""
    assert classify_significance_tier(float("nan")) == SignificanceTier.NS
    assert format_significance_stars(float("nan")) == ""
    assert format_significance_stars(0.05) == ""
    assert format_significance_stars(0.0499) == "*"
    assert format_significance_stars(0.01) == "*"
    assert format_significance_stars(0.0099) == "**"
    assert format_significance_stars(0.001) == "**"
    assert format_significance_stars(0.00099) == "***"
    assert format_significance_stars(0.0001) == "***"
    assert format_significance_stars(0.000099) == "****"
    assert format_significance_stars(1.37e-11) == "****"


@given(
    n_rows=st.integers(min_value=2, max_value=8),
    n_cols=st.integers(min_value=2, max_value=10),
    seed=st.integers(min_value=0, max_value=10_000),
)
@settings(max_examples=30)
def test_ward_biclustering_permutation_set_and_bounds(
    n_rows: int,
    n_cols: int,
    seed: int,
) -> None:
    """Property: Ward biclustering returns exact bijections and normalized leaf brackets in [0.5, K-0.5]."""
    rng = np.random.default_rng(seed)
    mat = rng.normal(0.0, 1.0, size=(n_rows, n_cols)).astype(np.float64)
    r_labels = [f"R{i}" for i in range(n_rows)]
    c_labels = [f"C{j}" for j in range(n_cols)]
    res = compute_ward_biclustering(mat, r_labels, c_labels)
    assert is_successful(res)
    bic = res.unwrap()
    assert set(bic.ordered_rows) == set(r_labels)
    assert len(bic.ordered_rows) == n_rows
    assert set(bic.ordered_cols) == set(c_labels)
    assert len(bic.ordered_cols) == n_cols
    assert len(bic.row_brackets) == n_rows - 1
    assert len(bic.col_brackets) == n_cols - 1
    for br in bic.row_brackets:
        assert 0.5 <= br.u_left <= float(n_rows) - 0.5
        assert 0.5 <= br.u_right <= float(n_rows) - 0.5
        assert 0.0 <= br.h_merge <= 1.0 + 1e-9
    for br in bic.col_brackets:
        assert 0.5 <= br.u_left <= float(n_cols) - 0.5
        assert 0.5 <= br.u_right <= float(n_cols) - 0.5
        assert 0.0 <= br.h_merge <= 1.0 + 1e-9


def test_strip_coulton_cluster_number_preserves_substate_identity() -> None:
    """Verify leading `<digits>_` cluster numbers are stripped while preserving unique substate suffixes."""
    assert strip_coulton_cluster_number("8_IFNGMac") == "IFNGMac"
    assert strip_coulton_cluster_number("18_ECMMac") == "ECMMac"
    assert strip_coulton_cluster_number("3_ICIMac1") == "ICIMac1"
    assert strip_coulton_cluster_number("4_ICIMac2") == "ICIMac2"
    assert strip_coulton_cluster_number("1_MetM2Mac") == "MetM2Mac"
    assert strip_coulton_cluster_number("IFNGMac") == "IFNGMac"
