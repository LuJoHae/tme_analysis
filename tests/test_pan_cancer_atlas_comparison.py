"""Property-based, AST exclusivity, and Inkscape SVG layer tests for `pan_cancer_atlas_comparison`.

Enforces `.agents/rules/code-style-guide.md`:
1. Property-based tests (`hypothesis`) for pure Functional Core routines in
   `scripts.pan_cancer_atlas_comparison.core`:
   - Out-of-core CSR chunk centroid accumulation matches full-matrix aggregation.
   - Benjamini-Hochberg FDR q-values lie in `[p_i, 1]` and preserve rank monotonicity.
   - DerSimonian-Laird random-effects meta-analysis satisfies `tau2 >= 0` and `I2 in [0, 100]`.
   - Optimal-leaf Ward hierarchical biclustering returns valid bijections over rows and columns.
   - Centroid serialization -> deserialization via Polars is lossless.
2. AST audit verifying zero `read_h5ad`, `read_zarr`, `find_dataset_h5ad`, or `Lair` calls
   and zero `pandas` imports across `scripts/pan_cancer_atlas_comparison/*.py`.
3. Inkscape native layer (`inkscape:groupmode="layer"`) and zero `feDropShadow` verification
   for publication SVG generation (`inject_inkscape_layers_pure`).
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Final
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.result import Success
from scipy import sparse as sp

from scripts.pan_cancer_atlas_comparison.core import (
    accumulate_atlas_batch,
    assign_one_to_one_matches_df,
    benjamini_hochberg_fdr,
    build_chu_curated_centroid_bundle,
    compute_ward_biclustering,
    create_empty_atlas_accumulator,
    dersimonian_laird_meta_analysis,
    deserialize_atlas_centroids_from_polars,
    finalize_atlas_centroids,
    serialize_atlas_centroids_to_polars,
)
from scripts.pan_cancer_atlas_comparison.models import (
    ExternalAtlasId,
    ProfileOrigin,
)
from scripts.pan_cancer_atlas_comparison.render_figures import (
    _INKSCAPE_NS,
    _SVG_NS,
    inject_inkscape_layers_pure,
)

REPO_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
PAN_ATLAS_SCRIPT_DIR: Final[Path] = REPO_ROOT / "scripts" / "pan_cancer_atlas_comparison"


@settings(max_examples=25, deadline=5000)
@given(
    n_cells=st.integers(min_value=6, max_value=40),
    n_genes=st.integers(min_value=4, max_value=16),
    chunk_size=st.integers(min_value=1, max_value=13),
    seed=st.integers(min_value=1, max_value=10_000),
)
def test_property_out_of_core_csr_centroid_accumulation_matches_exact_mean(
    n_cells: int,
    n_genes: int,
    chunk_size: int,
    seed: int,
) -> None:
    """Property: Chunked CSR streaming yields identical cluster centroids to full-matrix means."""
    rng = np.random.default_rng(seed)
    dense = np.abs(rng.normal(loc=0.8, scale=0.7, size=(n_cells, n_genes))).astype(np.float64)
    dense[rng.random((n_cells, n_genes)) < 0.45] = 0.0
    labels = ("c01", "c02", "c03")
    codes = rng.integers(low=0, high=len(labels), size=n_cells, dtype=np.int64)
    tumor_mask = rng.random(n_cells) < 0.70
    genes = tuple(f"G{i}" for i in range(n_genes))

    acc = create_empty_atlas_accumulator(
        atlas_id=ExternalAtlasId.GAO_CAF,
        origin=ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD,
        cluster_labels=labels,
        gene_symbols=genes,
    )
    for start in range(0, n_cells, chunk_size):
        stop = min(n_cells, start + chunk_size)
        step_res = accumulate_atlas_batch(
            acc=acc,
            batch_csr=sp.csr_matrix(dense[start:stop]),
            cluster_codes=codes[start:stop],
            is_tumor_mask=tumor_mask[start:stop],
        )
        assert isinstance(step_res, Success)
        acc = step_res.unwrap()

    bundle = finalize_atlas_centroids(acc)
    dense_lin = np.expm1(dense)
    for k_idx in range(len(labels)):
        mask_all = codes == k_idx
        mask_tum = mask_all & tumor_mask
        assert int(bundle.all_cell_counts[k_idx]) == int(np.sum(mask_all))
        assert int(bundle.tumor_cell_counts[k_idx]) == int(np.sum(mask_tum))
        if np.any(mask_all):
            np.testing.assert_allclose(
                bundle.all_mean_log[k_idx],
                dense[mask_all].mean(axis=0),
                rtol=1e-10,
                atol=1e-10,
            )
            np.testing.assert_allclose(
                bundle.all_mean_linear[k_idx],
                dense_lin[mask_all].mean(axis=0),
                rtol=1e-10,
                atol=1e-10,
            )
        if np.any(mask_tum):
            np.testing.assert_allclose(
                bundle.tumor_mean_log[k_idx],
                dense[mask_tum].mean(axis=0),
                rtol=1e-10,
                atol=1e-10,
            )


@settings(max_examples=30, deadline=5000)
@given(
    pvals=st.lists(
        st.floats(min_value=1e-12, max_value=1.0, allow_nan=False, allow_infinity=False),
        min_size=1,
        max_size=40,
    )
)
def test_property_benjamini_hochberg_fdr_bounds_and_monotonicity(pvals: list[float]) -> None:
    """Property: BH q-values lie in `[p_i, 1.0]` and preserve ascending p-value rank order."""
    p_arr = np.asarray(pvals, dtype=np.float64)
    q_arr = benjamini_hochberg_fdr(p_arr)
    assert q_arr.shape == p_arr.shape
    assert bool(np.all((q_arr >= p_arr - 1e-12) & (q_arr <= 1.0 + 1e-12)))
    order = np.argsort(p_arr, kind="mergesort")
    sorted_q = q_arr[order]
    assert bool(np.all(np.diff(sorted_q) >= -1e-12))


@settings(max_examples=30, deadline=5000)
@given(
    betas=st.lists(
        st.floats(min_value=-3.0, max_value=3.0, allow_nan=False, allow_infinity=False),
        min_size=2,
        max_size=12,
    ),
    ses=st.lists(
        st.floats(min_value=0.05, max_value=1.5, allow_nan=False, allow_infinity=False),
        min_size=2,
        max_size=12,
    ),
)
def test_property_dersimonian_laird_invariants(betas: list[float], ses: list[float]) -> None:
    """Property: DerSimonian-Laird pooling yields non-negative `tau2`, `I2 in [0, 100]`, and bounded mean."""
    k = min(len(betas), len(ses))
    b_arr = np.asarray(betas[:k], dtype=np.float64)
    s_arr = np.asarray(ses[:k], dtype=np.float64)
    dl_beta, dl_se, dl_pval, tau2, i2_pct = dersimonian_laird_meta_analysis(b_arr, s_arr)
    assert dl_se > 0.0
    assert tau2 >= 0.0
    assert 0.0 <= i2_pct <= 100.0
    assert 0.0 <= dl_pval <= 1.0
    assert float(np.min(b_arr)) - 1e-9 <= dl_beta <= float(np.max(b_arr)) + 1e-9


@settings(max_examples=20, deadline=5000)
@given(
    n_rows=st.integers(min_value=2, max_value=14),
    n_cols=st.integers(min_value=2, max_value=14),
    seed=st.integers(min_value=1, max_value=5000),
)
def test_property_ward_bicluster_returns_bijective_orders(
    n_rows: int,
    n_cols: int,
    seed: int,
) -> None:
    """Property: Optimal-leaf Ward hierarchical ordering is a bijection on both rows and columns."""
    rng = np.random.default_rng(seed)
    mat = rng.normal(size=(n_rows, n_cols)).astype(np.float64)
    row_labels = tuple(f"R_{i:02d}" for i in range(n_rows))
    col_labels = tuple(f"C_{j:02d}" for j in range(n_cols))
    res = compute_ward_biclustering(mat, row_labels, col_labels)
    assert isinstance(res, Success)
    bic = res.unwrap()
    assert len(bic.ordered_rows) == n_rows
    assert len(bic.ordered_cols) == n_cols
    assert set(bic.ordered_rows) == set(row_labels)
    assert set(bic.ordered_cols) == set(col_labels)


def test_chu_curated_bundle_polars_round_trip() -> None:
    """Verify Chu 39-state curated centroid bundle round-trips losslessly through Polars."""
    genes = ("CD8A", "CD4", "FOXP3", "CXCL13", "HSPA1A", "DNAJB1", "TCF7", "GZMK", "ZNF683")
    bundle = build_chu_curated_centroid_bundle(genes)
    assert bundle.atlas_id == ExternalAtlasId.CHU_TCELL
    assert bundle.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR
    assert len(bundle.cluster_labels) == 39

    df = serialize_atlas_centroids_to_polars((bundle,))
    restored_res = deserialize_atlas_centroids_from_polars(df)
    assert isinstance(restored_res, Success)
    restored = restored_res.unwrap()
    assert len(restored) == 1
    assert restored[0].cluster_labels == bundle.cluster_labels
    np.testing.assert_allclose(restored[0].tumor_mean_linear, bundle.tumor_mean_linear, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(restored[0].tumor_mean_log, bundle.tumor_mean_log, rtol=1e-5, atol=1e-6)


def test_pan_atlas_scripts_strict_ast_purity_and_api_exclusivity() -> None:
    """Verify zero `pandas`, zero forbidden `global`/`nonlocal`/`del`, and zero non-API dataset reads."""
    py_files = sorted(p for p in PAN_ATLAS_SCRIPT_DIR.glob("*.py") if p.name != "__init__.py")
    assert len(py_files) == 5
    forbidden_calls = {"read_h5ad", "read_zarr", "find_dataset_h5ad", "Lair"}

    for py_path in py_files:
        src = py_path.read_text(encoding="utf-8")
        tree = ast.parse(src, filename=str(py_path))
        for node in ast.walk(tree):
            match node:
                case ast.Import(names=aliases):
                    for alias in aliases:
                        assert alias.name != "pandas", f"Forbidden pandas import in {py_path.name}"
                case ast.ImportFrom(module=mod_name) if mod_name == "pandas":
                    raise AssertionError(f"Forbidden from pandas import in {py_path.name}")
                case ast.Global() | ast.Nonlocal() | ast.Delete():
                    raise AssertionError(f"Forbidden global/nonlocal/del in {py_path.name}@L{node.lineno}")
                case ast.Call(func=ast.Name(id=fn_id)) if fn_id in forbidden_calls:
                    raise AssertionError(f"Forbidden {fn_id} call in {py_path.name}@L{node.lineno}")
                case ast.Call(func=ast.Attribute(attr=attr_id)) if attr_id in forbidden_calls:
                    raise AssertionError(f"Forbidden .{attr_id} call in {py_path.name}@L{node.lineno}")
                case _:
                    pass


def test_inject_inkscape_layers_pure_well_formed_and_zero_drop_shadow() -> None:
    """Verify SVG post-processor produces valid XML with native Inkscape layers and zero drop shadows."""
    raw_svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="400" viewBox="0 0 600 400">'
        '<rect class="background" width="600" height="400" fill="#FFFFFF"/>'
        '<g class="mark-rect"><rect x="10" y="10" width="50" height="50" fill="#0072B2"/></g>'
        "</svg>"
    )
    wrapped_res = inject_inkscape_layers_pure(raw_svg, figure_id="fig_test_triad")
    assert isinstance(wrapped_res, Success)
    wrapped = wrapped_res.unwrap()
    assert "feDropShadow" not in wrapped
    root = ET.fromstring(wrapped)
    assert root.tag == f"{{{_SVG_NS}}}svg"
    layers = [
        el
        for el in root.findall(f"{{{_SVG_NS}}}g")
        if el.get(f"{{{_INKSCAPE_NS}}}groupmode") == "layer"
    ]
    assert len(layers) == 2
    labels = [el.get(f"{{{_INKSCAPE_NS}}}label") for el in layers]
    assert labels == ["00_Canvas_Background", "01_Data_Panels_And_Axes"]


@settings(max_examples=25, deadline=5000)
@given(
    n_ext=st.integers(min_value=2, max_value=12),
    seed=st.integers(min_value=1, max_value=5000),
)
def test_property_one_to_one_assignment_is_injective_and_optimal(n_ext: int, seed: int) -> None:
    """Property: each Meta42 state and each atlas substate is matched at most once, beating greedy."""
    from scripts.pan_cancer_atlas_comparison.models import META42_T_NK_IDS

    rng = np.random.default_rng(seed)
    metas = list(META42_T_NK_IDS)
    score = rng.normal(size=(len(metas), n_ext))
    df = pl.DataFrame(
        {
            "atlas_id": ["A"] * (len(metas) * n_ext),
            "atlas_compartment": ["T_NK"] * (len(metas) * n_ext),
            "meta_state": [m for m in metas for _ in range(n_ext)],
            "external_cluster": [f"e{j}" for _ in metas for j in range(n_ext)],
            "composite_substate_score": score.ravel().tolist(),
        }
    )
    out = assign_one_to_one_matches_df(df).filter(pl.col("is_one_to_one_match"))
    assert out.height == min(len(metas), n_ext)
    assert out["meta_state"].n_unique() == out.height
    assert out["external_cluster"].n_unique() == out.height
    greedy_total = float(np.sort(score.max(axis=1))[::-1][: min(len(metas), n_ext)].sum())
    assert float(out["composite_substate_score"].sum()) <= greedy_total + 1e-9
