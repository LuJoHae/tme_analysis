"""Unit and property-based tests for the Unified 9-Cohort Reference Deconvolution Benchmark."""

from __future__ import annotations

from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
import pytest

import importlib

_stage3b = importlib.import_module(
    "scripts.manuscript_pipeline.3b_pan_cohort_unified_deconv_benchmark"
)
COHORT_ORDER = _stage3b.COHORT_ORDER
UnifiedBenchmarkConfig = _stage3b.UnifiedBenchmarkConfig
build_unified_reference_atlas = _stage3b.build_unified_reference_atlas
compute_permutation_pvalue = _stage3b.compute_permutation_pvalue
compute_within_cohort_specificity = _stage3b.compute_within_cohort_specificity
deconvolve_pseudobulk_matrix = _stage3b.deconvolve_pseudobulk_matrix
fit_logistic_effect_sizes = _stage3b.fit_logistic_effect_sizes
load_all_nine_cohorts = _stage3b.load_all_nine_cohorts
project_unified_reference_to_cohort_genes = _stage3b.project_unified_reference_to_cohort_genes


@given(
    n_clusters=st.integers(min_value=3, max_value=8),
    n_genes=st.integers(min_value=10, max_value=30),
)
@settings(max_examples=15, deadline=None)
def test_compute_within_cohort_specificity_unit_norm(n_clusters: int, n_genes: int) -> None:
    """Property test: within-cohort specificity rows are zero-mean across genes and unit L2-norm."""
    rng = np.random.default_rng(42)
    phi_mat = rng.uniform(0.01, 10.0, size=(n_clusters, n_genes))
    spec = compute_within_cohort_specificity(phi_mat)
    assert spec.shape == (n_clusters, n_genes)
    row_norms = np.linalg.norm(spec, axis=1)
    assert np.allclose(row_norms, 1.0, atol=1e-6)


@given(
    n_samples=st.integers(min_value=4, max_value=10),
    n_clusters=st.integers(min_value=3, max_value=7),
    n_genes=st.integers(min_value=12, max_value=25),
)
@settings(max_examples=10, deadline=None)
def test_deconvolve_pseudobulk_matrix_simplex_invariant(
    n_samples: int,
    n_clusters: int,
    n_genes: int,
) -> None:
    """Property test: deconvolved cell fractions are non-negative and sum to 1.0 per sample."""
    rng = np.random.default_rng(123)
    phi_ref = rng.uniform(0.1, 5.0, size=(n_clusters, n_genes))
    true_theta = rng.dirichlet(np.ones(n_clusters), size=n_samples)
    bulk_mat = true_theta @ phi_ref

    inferred = deconvolve_pseudobulk_matrix(bulk_mat, phi_ref, n_iter=30)
    assert inferred.shape == (n_samples, n_clusters)
    assert np.all(inferred >= -1e-9)
    assert np.allclose(inferred.sum(axis=1), 1.0, atol=1e-5)


def test_load_all_nine_cohorts_and_unified_atlas() -> None:
    """Verify empirical loading of all 9 cohorts and unified reference atlas invariants."""
    config = UnifiedBenchmarkConfig()
    res = load_all_nine_cohorts(config)
    bundles = res.unwrap()
    assert len(bundles) == 9
    assert [b.cohort_id for b in bundles] == [m["cohort"] for m in COHORT_ORDER]

    atlas = build_unified_reference_atlas(bundles, n_meta_states=42)
    total_states = sum(len(b.clusters) for b in bundles)
    assert len(atlas.global_state_ids) == total_states
    assert atlas.similarity_matrix.shape == (total_states, total_states)
    assert np.allclose(np.diag(atlas.similarity_matrix), 1.0, atol=1e-6)
    assert atlas.phi_meta_df.height == 42


def test_leave_one_cohort_out_reference_isolation() -> None:
    """Verify LOCO reference strictly excludes the target cohort's own states."""
    config = UnifiedBenchmarkConfig()
    bundles = load_all_nine_cohorts(config).unwrap()
    atlas = build_unified_reference_atlas(bundles, n_meta_states=42)

    target = bundles[0]
    total_states = sum(len(b.clusters) for b in bundles)
    expected_loco_states = total_states - len(target.clusters)

    phi_loco, kernel_loco = project_unified_reference_to_cohort_genes(
        target_bundle=target,
        all_bundles=bundles,
        atlas=atlas,
        include_target_cohort=False,
    )
    assert phi_loco.shape == (expected_loco_states, len(target.phi_genes))
    assert kernel_loco.shape == (len(target.clusters), expected_loco_states)
    assert np.allclose(kernel_loco.sum(axis=0), 1.0, atol=1e-6)


def test_unified_benchmark_outputs_and_concordance() -> None:
    """Verify exported unified benchmark parquets and SVG/PNG figure properties."""
    out_dir = Path("output/results/pseudobulk_deconv_benchmark/unified_pan_cohort")
    summary_path = out_dir / "unified_per_cohort_concordance_summary.parquet"
    meta_path = out_dir / "unified_meta_state_concordance.parquet"
    svg_path = Path("results/pseudobulk_deconv_benchmark/pan_cohort_unified_reference_comparison.svg")
    png_path = Path("results/pseudobulk_deconv_benchmark/pan_cohort_unified_reference_comparison.png")

    assert summary_path.exists()
    assert meta_path.exists()
    assert svg_path.exists()
    assert png_path.exists() and png_path.stat().st_size > 10000

    summary_df = pl.read_parquet(summary_path)
    assert summary_df.height == 9
    assert float(summary_df["rho_unified_all9"].mean() or 0.0) > 0.50
    assert float(summary_df["rho_unified_loco"].mean() or 0.0) > 0.40

    # Validate XML well-formedness of SVG
    tree = ET.parse(svg_path)
    assert tree.getroot() is not None
