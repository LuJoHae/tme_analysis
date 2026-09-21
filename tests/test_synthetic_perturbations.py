"""
Unit and integration tests for the synthetic perturbation benchmarking suite.
Verifies base generation, modular perturbation mechanics, compound regimes, and single simulation stability.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pytest

from scripts.synthetic_perturbation_benchmark import (
    COMPOUND_REGIMES,
    PerturbationProfile,
    SweepConfig,
    apply_perturbations,
    assign_patients,
    generate_base_data,
    make_compound_profile,
    run_single_simulation,
)


def test_generate_base_data_structure() -> None:
    """Verify that base synthetic data generation produces expected matrix dimensions and markers."""
    config = SweepConfig(n_cells=600, n_genes=60, n_clusters=4, n_patients=10)
    rng = np.random.default_rng(42)

    raw_counts, cluster_assignments, marker_map = generate_base_data(
        config=config,
        collinearity_overlap=0.0,
        rng=rng,
    )

    assert raw_counts.shape == (600, 60)
    assert len(cluster_assignments) == 600
    assert len(marker_map) == 4
    assert np.all(raw_counts >= 0.0)
    assert set(np.unique(cluster_assignments)) == {0, 1, 2, 3}


def test_assign_patients_proportions() -> None:
    """Verify patient assignment and response stratification."""
    cluster_assignments = np.repeat(np.arange(4, dtype=np.int32), 100)
    rng = np.random.default_rng(42)

    patient_assignments, patient_responses, cluster_probs = assign_patients(
        cluster_assignments=cluster_assignments,
        n_patients=10,
        sparsity_intensity=0.0,
        rng=rng,
    )

    assert len(patient_assignments) == 400
    assert len(patient_responses) == 400
    assert len(np.unique(patient_assignments)) == 10
    assert set(patient_responses) == {"R", "NR"}
    assert cluster_probs[0] > cluster_probs[3]  # Cluster 0 enriched in Responders


def test_apply_perturbations_cell_size() -> None:
    """Verify that cell_size perturbation scales pseudobulk expression without changing single-cell counts."""
    n_cells = 400
    n_genes = 40
    k = 4
    raw_counts = np.ones((n_cells, n_genes), dtype=np.float64)
    cluster_assignments = np.repeat(np.arange(k, dtype=np.int32), 100)
    patient_assignments = np.array([f"P_{i % 10:02d}" for i in range(n_cells)])
    patient_responses = np.array(["R" if (i % 10) < 5 else "NR" for i in range(n_cells)])
    rng = np.random.default_rng(42)

    sc_counts, _, pseudobulk_mat, _ = apply_perturbations(
        raw_counts=raw_counts,
        cluster_assignments=cluster_assignments,
        patient_assignments=patient_assignments,
        patient_responses=patient_responses,
        mode="cell_size",
        intensity=1.0,  # Max intensity
        rng=rng,
    )

    # Single-cell counts must remain unaffected
    assert np.array_equal(sc_counts, raw_counts)
    # Pseudobulk matrix has shape (10 patients, 40 genes)
    assert pseudobulk_mat.shape == (10, 40)
    assert np.all(pseudobulk_mat > 0.0)


def test_apply_ambient_soup() -> None:
    """Verify that ambient soup blends background RNA into single cells."""
    n_cells = 200
    n_genes = 30
    raw_counts = np.zeros((n_cells, n_genes), dtype=np.float64)
    # Cell 0 has only gene 0
    raw_counts[0, 0] = 100.0
    # Cell 1 has only gene 1
    raw_counts[1, 1] = 100.0
    cluster_assignments = np.zeros(n_cells, dtype=np.int32)
    patient_assignments = np.array([f"P_{i % 4:02d}" for i in range(n_cells)])
    patient_responses = np.array(["R" if (i % 4) < 2 else "NR" for i in range(n_cells)])
    rng = np.random.default_rng(42)

    profile = PerturbationProfile(ambient_soup=1.0)
    sc_counts, _, _, _ = apply_perturbations(
        raw_counts=raw_counts,
        cluster_assignments=cluster_assignments,
        patient_assignments=patient_assignments,
        patient_responses=patient_responses,
        mode=profile,
        rng=rng,
    )

    assert sc_counts.shape == raw_counts.shape
    # Gene 1 in Cell 0 should now be > 0 due to ambient soup contamination
    assert sc_counts[0, 1] > 0.0


def test_apply_marker_dysregulation() -> None:
    """Verify marker dysregulation scales pseudobulk marker expression."""
    n_cells = 200
    n_genes = 20
    raw_counts = np.ones((n_cells, n_genes), dtype=np.float64)
    cluster_assignments = np.repeat(np.arange(2, dtype=np.int32), 100)
    patient_assignments = np.array([f"P_{i % 4:02d}" for i in range(n_cells)])
    patient_responses = np.array(["R" if (i % 4) < 2 else "NR" for i in range(n_cells)])
    marker_map = [[0, 1, 2], [3, 4, 5]]
    rng = np.random.default_rng(42)

    profile = PerturbationProfile(marker_dysregulation=1.0)
    sc_counts, _, pseudobulk_mat, _ = apply_perturbations(
        raw_counts=raw_counts,
        cluster_assignments=cluster_assignments,
        patient_assignments=patient_assignments,
        patient_responses=patient_responses,
        mode=profile,
        marker_map=marker_map,
        rng=rng,
    )

    assert pseudobulk_mat.shape == (4, 20)
    assert np.all(pseudobulk_mat > 0.0)


def test_perturbation_profile_compound() -> None:
    """Verify compound profile instantiation and weighting."""
    prof = make_compound_profile("core_biopsy", 0.5)
    assert prof.cell_size == 0.5
    assert prof.sampling_sparsity == pytest.approx(0.4)
    assert prof.ghost_contamination == pytest.approx(0.35)
    assert prof.activation_confounding == 0.0


def test_run_single_simulation_end_to_end() -> None:
    """Verify end-to-end execution of a small simulation instance."""
    config = SweepConfig(
        n_cells=400,
        n_genes=40,
        n_clusters=4,
        n_patients=6,
        deconv_iters=30,
    )
    rng = np.random.default_rng(42)

    res, scatter = run_single_simulation(
        config=config,
        mode="cell_size",
        intensity=0.0,
        rep_idx=0,
        rng=rng,
    )

    assert "spearman_rho_milo_deconv" in res
    assert "sign_concordance_pct" in res
    assert "fidelity_spearman_rho" in res
    assert len(scatter) == config.n_clusters

    # Baseline (intensity=0) should exhibit positive concordance
    assert not np.isnan(res["spearman_rho_milo_deconv"])
    assert res["sign_concordance_pct"] >= 50.0


def test_run_compound_simulation_end_to_end() -> None:
    """Verify end-to-end execution of a compound perturbation simulation."""
    config = SweepConfig(
        n_cells=400,
        n_genes=40,
        n_clusters=4,
        n_patients=6,
        deconv_iters=30,
    )
    rng = np.random.default_rng(42)
    profile = make_compound_profile("triple_jeopardy", 0.5)

    res, scatter = run_single_simulation(
        config=config,
        mode=profile,
        intensity=0.5,
        rep_idx=0,
        rng=rng,
        label="triple_jeopardy",
    )

    assert res["perturbation_mode"] == "triple_jeopardy"
    assert "spearman_rho_milo_deconv" in res
    assert len(scatter) == config.n_clusters
