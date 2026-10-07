"""Unit tests for multi-cohort single-cell sampling, gene alignment, PCA, and Leiden clustering."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from returns.maybe import Nothing, Some
from returns.result import Success

from tme_datasets import (
    ClusterAnalysisSpec,
    CohortSamplingMode,
    HarmonizeConfig,
    HarmonizeMode,
    SampledSingleCellResult,
    SingleCellSamplingSpec,
    compute_cohort_cell_allocations,
    find_dataset_h5ad,
    generate_random_cluster_spec,
    generate_random_sampling_spec,
    resolve_sampled_cohorts,
    run_pca_knn_leiden,
    sample_single_cell_cohorts,
)


def test_resolve_sampled_cohorts() -> None:
    """Verify deterministic cohort resolution and filtering."""
    # 1. Explicit cohorts
    spec_explicit = SingleCellSamplingSpec(
        cohort_ids=Some(("GSE120575", "Maynard_NSCLC")),
        require_cached_h5ad=False,
    )
    res_explicit = resolve_sampled_cohorts(spec_explicit)
    assert isinstance(res_explicit, Success)
    assert res_explicit.unwrap() == ("GSE120575", "Maynard_NSCLC")

    # 2. Cancer type filter (Melanoma)
    spec_melanoma = SingleCellSamplingSpec(
        cancer_types=Some(("Melanoma",)),
        require_cached_h5ad=False,
    )
    res_melanoma = resolve_sampled_cohorts(spec_melanoma)
    assert isinstance(res_melanoma, Success)
    melanoma_ids = res_melanoma.unwrap()
    assert "GSE120575" in melanoma_ids
    assert "Maynard_NSCLC" not in melanoma_ids

    # 3. Deterministic random cohort sampling with seed
    spec_random1 = SingleCellSamplingSpec(
        n_cohorts=Some(3),
        seed=Some(42),
        require_cached_h5ad=False,
    )
    spec_random2 = SingleCellSamplingSpec(
        n_cohorts=Some(3),
        seed=Some(42),
        require_cached_h5ad=False,
    )
    res1 = resolve_sampled_cohorts(spec_random1)
    res2 = resolve_sampled_cohorts(spec_random2)
    assert isinstance(res1, Success)
    assert isinstance(res2, Success)
    assert len(res1.unwrap()) == 3
    assert res1.unwrap() == res2.unwrap()  # Exactly reproducible with seed


def test_compute_cohort_cell_allocations() -> None:
    """Verify all 5 cell allocation modes across cohorts."""
    cohort_ids = ("cohort_A", "cohort_B", "cohort_C")
    cohort_sizes = {"cohort_A": 500, "cohort_B": 2000, "cohort_C": 5000}

    # 1. Fixed per cohort
    spec_fixed = SingleCellSamplingSpec(
        mode=CohortSamplingMode.FIXED_PER_COHORT,
        n_cells_per_cohort=1000,
    )
    alloc_fixed = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, spec_fixed)
    assert alloc_fixed["cohort_A"] == 500  # Clamped to total
    assert alloc_fixed["cohort_B"] == 1000
    assert alloc_fixed["cohort_C"] == 1000

    # 2. Fraction per cohort (10%)
    spec_frac = SingleCellSamplingSpec(
        mode=CohortSamplingMode.FRACTION_PER_COHORT,
        fraction_per_cohort=0.10,
    )
    alloc_frac = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, spec_frac)
    assert alloc_frac["cohort_A"] == 50
    assert alloc_frac["cohort_B"] == 200
    assert alloc_frac["cohort_C"] == 500

    # 3. Explicit counts
    spec_counts = SingleCellSamplingSpec(
        mode=CohortSamplingMode.EXPLICIT_COUNTS,
        explicit_cell_counts={"cohort_A": 150, "cohort_B": 800},
        n_cells_per_cohort=300,  # Fallback for cohort_C
    )
    alloc_counts = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, spec_counts)
    assert alloc_counts["cohort_A"] == 150
    assert alloc_counts["cohort_B"] == 800
    assert alloc_counts["cohort_C"] == 300

    # 4. Explicit fractions
    spec_fractions = SingleCellSamplingSpec(
        mode=CohortSamplingMode.EXPLICIT_FRACTIONS,
        explicit_cell_fractions={"cohort_A": 0.5, "cohort_B": 0.05},
        fraction_per_cohort=0.02,  # Fallback for cohort_C
    )
    alloc_fractions = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, spec_fractions)
    assert alloc_fractions["cohort_A"] == 250
    assert alloc_fractions["cohort_B"] == 100
    assert alloc_fractions["cohort_C"] == 100

    # 5. Global budget (Total 750 cells)
    spec_budget = SingleCellSamplingSpec(
        mode=CohortSamplingMode.GLOBAL_BUDGET,
        global_cell_budget=Some(750),
    )
    alloc_budget = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, spec_budget)
    # Total available: 7500. Proportions: A=500/7500=0.0667, B=2000/7500=0.2667, C=5000/7500=0.6667
    assert alloc_budget["cohort_A"] == 50
    assert alloc_budget["cohort_B"] == 200
    assert alloc_budget["cohort_C"] == 500
    assert sum(alloc_budget.values()) == 750


def test_run_pca_knn_leiden_synthetic() -> None:
    """Verify HVG, joint PCA, kNN graph, and Leiden clustering on synthetic multi-batch AnnData."""
    rng = np.random.default_rng(42)
    # 60 cells x 80 genes
    X_dense = rng.poisson(lam=3.0, size=(60, 80)).astype(np.float32)
    # Introduce cluster signal: cells 0-30 have high genes 0-10, cells 30-60 have high genes 11-20
    X_dense[:30, :10] += 15.0
    X_dense[30:, 10:20] += 15.0

    obs = pd.DataFrame({
        "dataset_id": ["cohort_1"] * 30 + ["cohort_2"] * 30,
    }, index=[f"cell_{i}" for i in range(60)])
    var = pd.DataFrame(index=[f"GENE_{j}" for j in range(80)])

    adata = ad.AnnData(X=sp.csr_matrix(X_dense), obs=obs, var=var)

    spec = ClusterAnalysisSpec(
        n_top_genes=40,
        n_pcs=10,
        n_neighbors=8,
        leiden_resolution=0.6,
        seed=Some(42),
    )

    res = run_pca_knn_leiden(adata, spec)
    assert isinstance(res, Success)
    clustered = res.unwrap()

    # PCA checks
    assert "X_pca" in clustered.obsm
    assert clustered.obsm["X_pca"].shape == (60, 10)
    assert "pca" in clustered.uns

    # kNN graph checks
    assert "neighbors" in clustered.uns
    assert "connectivities" in clustered.obsp
    assert "distances" in clustered.obsp

    # Leiden checks
    assert "leiden" in clustered.obs.columns
    unique_clusters = clustered.obs["leiden"].unique()
    assert len(unique_clusters) >= 2


@pytest.mark.skipif(
    not (isinstance(find_dataset_h5ad("GSE120575"), Some) and isinstance(find_dataset_h5ad("Maynard_NSCLC"), Some)),
    reason="Requires cached processed GSE120575 and Maynard_NSCLC datasets",
)
def test_sample_single_cell_cohorts_end_to_end() -> None:
    """Verify complete multi-cohort sampling, gene intersection, PCA, and Leiden on real cached datasets."""
    sampling_spec = SingleCellSamplingSpec(
        cohort_ids=Some(("GSE120575", "Maynard_NSCLC")),
        mode=CohortSamplingMode.FIXED_PER_COHORT,
        n_cells_per_cohort=100,
        seed=Some(42),
    )
    cluster_spec = ClusterAnalysisSpec(
        n_top_genes=500,
        n_pcs=15,
        n_neighbors=10,
        leiden_resolution=0.5,
        seed=Some(42),
    )
    harmonize_config = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        min_shared_genes=1000,
    )

    res = sample_single_cell_cohorts(
        sampling_spec=sampling_spec,
        cluster_spec=cluster_spec,
        harmonize_config=harmonize_config,
    )
    assert isinstance(res, Success)
    result: SampledSingleCellResult = res.unwrap()

    # Verify result model structure
    assert result.total_cells == 200
    assert result.cells_per_cohort["GSE120575"] == 100
    assert result.cells_per_cohort["Maynard_NSCLC"] == 100
    assert result.n_shared_genes > 1000
    assert result.n_clusters >= 2

    # Verify AnnData contents
    adata = result.adata
    assert adata.n_obs == 200
    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape[1] == 15
    assert "leiden" in adata.obs.columns
    assert "dataset_id" in adata.obs.columns
    assert "sampling_summary" in adata.uns

    # Ensure all processed single-cell features are canonical Ensembl IDs
    assert all(g.startswith("ENSG") for g in adata.var_names)


def test_generate_random_sampling_spec() -> None:
    """Verify randomized parameter generation across diverse sampling modes."""
    # 1. Randomized explicit counts mode
    spec_counts = generate_random_sampling_spec(
        available_cohorts=("GSE120575", "Maynard_NSCLC"),
        mode=CohortSamplingMode.EXPLICIT_COUNTS,
        min_cells_per_cohort=100,
        max_cells_per_cohort=500,
        seed=42,
    )
    assert spec_counts.mode == CohortSamplingMode.EXPLICIT_COUNTS
    assert len(spec_counts.explicit_cell_counts) >= 1
    for count in spec_counts.explicit_cell_counts.values():
        assert 100 <= count <= 500

    # 2. Randomized explicit fractions mode
    spec_fracs = generate_random_sampling_spec(
        available_cohorts=("GSE120575", "Maynard_NSCLC"),
        mode=CohortSamplingMode.EXPLICIT_FRACTIONS,
        min_fraction=0.05,
        max_fraction=0.20,
        seed=42,
    )
    assert spec_fracs.mode == CohortSamplingMode.EXPLICIT_FRACTIONS
    assert len(spec_fracs.explicit_cell_fractions) >= 1
    for frac in spec_fracs.explicit_cell_fractions.values():
        assert 0.05 <= frac <= 0.20

    # 3. Randomized without explicit mode (picks random mode)
    spec_auto = generate_random_sampling_spec(seed=123)
    assert isinstance(spec_auto.mode, CohortSamplingMode)
    assert isinstance(spec_auto.cohort_ids, Some)
    assert len(spec_auto.cohort_ids.unwrap()) >= 2

    # 4. Deterministic reproducibility with identical seed
    spec_a = generate_random_sampling_spec(seed=999)
    spec_b = generate_random_sampling_spec(seed=999)
    assert spec_a.mode == spec_b.mode
    assert spec_a.cohort_ids == spec_b.cohort_ids
    assert spec_a.n_cells_per_cohort == spec_b.n_cells_per_cohort
    assert spec_a.fraction_per_cohort == spec_b.fraction_per_cohort


def test_generate_random_cluster_spec() -> None:
    """Verify randomized cluster analysis parameter generation."""
    spec = generate_random_cluster_spec(
        min_pcs=15,
        max_pcs=30,
        min_neighbors=10,
        max_neighbors=25,
        seed=42,
    )
    assert 15 <= spec.n_pcs <= 30
    assert 10 <= spec.n_neighbors <= 25
    assert 0.3 <= spec.leiden_resolution <= 1.5
    assert spec.n_top_genes in [500, 1000, 1500, 2000, 3000]


def test_sampler_qc_filtering_cells_and_genes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify sampler strictly excludes cells that failed QC and only sees QC-passing genes."""
    from tme_datasets.sampling.cohort_sampler import sample_and_harmonize_cohorts

    # Create synthetic cohort: 10 cells (5 pass, 5 fail) x 6 genes (4 pass, 2 fail)
    # Cell 0-4: passing QC (is_retained_qc = True)
    # Cell 5-9: failing QC (is_retained_qc = False, low counts/droplets)
    # Gene 0-3: expressed in passing cells (pass QC)
    # Gene 4: 0 counts across all cells (fail QC)
    # Gene 5: explicitly flagged as is_retained_qc = False (fail QC)
    X = np.zeros((10, 6), dtype=np.float32)
    # Fill expression for passing cells on genes 0-3
    for c in range(5):
        X[c, :4] = [10.0, 15.0, 8.0, 12.0]
    # Failing cells have low or random counts
    X[5:, :2] = 1.0

    obs = pd.DataFrame(
        {
            "is_retained_qc": [True, True, True, True, True, False, False, False, False, False],
            "total_counts": [1000.0] * 5 + [20.0] * 5,
        },
        index=[f"cell_{i}" for i in range(10)],
    )
    var = pd.DataFrame(
        {
            "gene_name": [f"GENE_{i}" for i in range(6)],
            "is_retained_qc": [True, True, True, True, True, False],
        },
        index=[f"ENSG0000000000{i}" for i in range(6)],
    )
    mock_adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

    monkeypatch.setattr(
        "tme_datasets.sampling.cohort_sampler.load_dataset",
        lambda cid, **kwargs: Success(mock_adata.copy()),
    )

    sampling_spec = SingleCellSamplingSpec(
        cohort_ids=Some(("mock_cohort",)),
        mode=CohortSamplingMode.FIXED_PER_COHORT,
        n_cells_per_cohort=4,
        only_qc_passing_cells=True,
        only_qc_passing_genes=True,
        min_cells_per_gene=3,
        seed=Some(42),
    )
    harmonize_config = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        min_shared_genes=1,
    )

    res = sample_and_harmonize_cohorts(
        cohort_ids=("mock_cohort",),
        cell_allocations={"mock_cohort": 4},
        spec=sampling_spec,
        harmonize_config=harmonize_config,
    )
    assert isinstance(res, Success)
    sampled = res.unwrap()

    # 1. Assert only QC-passing cells were sampled
    assert sampled.n_obs == 4
    # Every sampled cell must have is_retained_qc == True
    assert all(sampled.obs["is_retained_qc"] == True)
    # The cell IDs must only be drawn from cell_0 through cell_4
    assert all(any(c in name for c in ["cell_0", "cell_1", "cell_2", "cell_3", "cell_4"]) for name in sampled.obs_names)
    assert not any(any(c in name for c in ["cell_5", "cell_6", "cell_7", "cell_8", "cell_9"]) for name in sampled.obs_names)

    # 2. Assert sampler only sees QC-passing genes
    # Gene 4 (0 counts) and Gene 5 (is_retained_qc=False) must be completely excluded
    assert sampled.n_vars == 4
    assert set(sampled.var_names) == {f"ENSG0000000000{i}" for i in range(4)}
    assert "ENSG00000000004" not in sampled.var_names
    assert "ENSG00000000005" not in sampled.var_names


def test_sampler_qc_filtering_adaptive_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify on-the-fly adaptive QC executes when raw AnnData lacks is_retained_qc annotations."""
    from tme_datasets.sampling.cohort_sampler import sample_and_harmonize_cohorts

    # AnnData with 50 cells x 20 genes (no is_retained_qc column)
    # 40 viable cells with high counts and genes, 10 dead cells / droplets with 5 counts
    rng = np.random.default_rng(42)
    X = np.zeros((50, 20), dtype=np.float32)
    X[:40, :] = rng.poisson(lam=10.0, size=(40, 20)).astype(np.float32)
    X[40:, :2] = 1.0  # 10 droplets with total counts = 2

    obs = pd.DataFrame(index=[f"c_{i}" for i in range(50)])
    var = pd.DataFrame(
        {"gene_name": [f"G_{i}" for i in range(20)]},
        index=[f"ENSG000000000{i:02d}" for i in range(20)],
    )
    raw_mock = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

    monkeypatch.setattr(
        "tme_datasets.sampling.cohort_sampler.load_dataset",
        lambda cid, **kwargs: Success(raw_mock.copy()),
    )

    from tme_datasets.models import QualityControlSpec

    spec = SingleCellSamplingSpec(
        cohort_ids=Some(("raw_cohort",)),
        mode=CohortSamplingMode.FIXED_PER_COHORT,
        n_cells_per_cohort=15,
        only_qc_passing_cells=True,
        only_qc_passing_genes=True,
        min_cells_per_gene=3,
        qc_spec=Some(QualityControlSpec(min_counts_per_cell=50, min_genes_per_cell=5, max_pct_mitochondrial=50.0)),
        seed=Some(42),
    )
    harmonize_config = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        min_shared_genes=1,
    )

    res = sample_and_harmonize_cohorts(
        cohort_ids=("raw_cohort",),
        cell_allocations={"raw_cohort": 15},
        spec=spec,
        harmonize_config=harmonize_config,
    )
    assert isinstance(res, Success)
    sampled = res.unwrap()
    assert sampled.n_obs == 15
    # The droplets (c_40 to c_49) must NOT be among the sampled cells
    for name in sampled.obs_names:
        assert not any(f"c_{d}" in name for d in range(40, 50))

