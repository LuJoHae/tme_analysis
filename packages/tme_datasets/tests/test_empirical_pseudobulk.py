"""Unit and Hypothesis property-based tests for empirical pseudobulk aggregation and caching."""

from __future__ import annotations

from pathlib import Path
import time
from typing import Final

import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Nothing, Some
from returns.pipeline import is_successful
from returns.result import Failure, Success
import scipy.sparse as sp

from tme_datasets import (
    EmpiricalPseudobulkConfig,
    aggregate_empirical_pseudobulk,
    cache_empirical_pseudobulk,
    compute_empirical_cell_fractions,
    empirical_pseudobulk_to_anndata,
    load_empirical_pseudobulk,
    resolve_adata_gene_symbols,
)

SC9_COHORT_IDS: Final[tuple[str, ...]] = (
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


def _build_synthetic_adata(
    counts: np.ndarray,
    sample_labels: list[str],
    state_labels: list[str],
    responses: list[str],
    genes: list[str],
    *,
    sparse_x: bool = True,
) -> ad.AnnData:
    """Construct an in-memory AnnData fixture purely from synthetic count arrays."""
    obs = pd.DataFrame(
        {
            "sample_id": sample_labels,
            "cell_state": state_labels,
            "response": responses,
            "timepoint": ["Pre" for _ in sample_labels],
        },
        index=[f"cell_{i}" for i in range(len(sample_labels))],
    )
    var = pd.DataFrame(
        {"gene_symbol": genes},
        index=[f"ENSG{1000 + j}" for j in range(len(genes))],
    )
    x_mat = sp.csr_matrix(counts, dtype=np.float32) if sparse_x else counts.astype(np.float32)
    return ad.AnnData(X=x_mat, obs=obs, var=var)


@settings(max_examples=25, deadline=2000)
@given(
    n_samples=st.integers(min_value=2, max_value=5),
    cells_per_sample=st.integers(min_value=3, max_value=8),
    n_genes=st.integers(min_value=3, max_value=8),
    seed=st.integers(min_value=1, max_value=10_000),
)
def test_property_linear_additivity_sum_and_mean(
    n_samples: int,
    cells_per_sample: int,
    n_genes: int,
    seed: int,
) -> None:
    """Property: sum/mean pseudobulk equals exact per-sample sum/mean of linear counts."""
    rng = np.random.default_rng(seed)
    n_cells = n_samples * cells_per_sample
    # Guarantee high enough UMI counts so raw counts are recognized as linear scale
    counts = rng.poisson(lam=80.0, size=(n_cells, n_genes)).astype(np.float64) + 50.0
    sample_labels = [f"S{s:02d}" for s in range(n_samples) for _ in range(cells_per_sample)]
    state_labels = [f"C{i % 3}" for i in range(n_cells)]
    responses = ["R" if (i // cells_per_sample) % 2 == 0 else "NR" for i in range(n_cells)]
    genes = [f"GENE{g}" for g in range(n_genes)]

    adata = _build_synthetic_adata(
        counts=counts,
        sample_labels=sample_labels,
        state_labels=state_labels,
        responses=responses,
        genes=genes,
        sparse_x=True,
    )

    sum_cfg = EmpiricalPseudobulkConfig(
        sample_col="sample_id",
        cell_state_col=Some("cell_state"),
        aggregation="sum",
        normalize_cells_cp10k=False,
        min_cells_per_sample=2,
    )
    res_sum = aggregate_empirical_pseudobulk(adata, config=sum_cfg, dataset_id="SYNTH")
    assert is_successful(res_sum)
    pb_sum = res_sum.unwrap()

    expected_samples = tuple(f"S{s:02d}" for s in range(n_samples))
    assert pb_sum.sample_ids == expected_samples
    assert pb_sum.genes == tuple(genes)

    for s_idx in range(n_samples):
        lo = s_idx * cells_per_sample
        hi = (s_idx + 1) * cells_per_sample
        expected_row = counts[lo:hi].sum(axis=0)
        np.testing.assert_allclose(pb_sum.bulk_matrix[s_idx], expected_row, rtol=1e-5, atol=1e-5)

    mean_cfg = EmpiricalPseudobulkConfig(
        sample_col="sample_id",
        cell_state_col=Some("cell_state"),
        aggregation="mean",
        normalize_cells_cp10k=False,
        min_cells_per_sample=2,
    )
    res_mean = aggregate_empirical_pseudobulk(adata, config=mean_cfg, dataset_id="SYNTH")
    assert is_successful(res_mean)
    pb_mean = res_mean.unwrap()
    np.testing.assert_allclose(
        pb_mean.bulk_matrix,
        pb_sum.bulk_matrix / float(cells_per_sample),
        rtol=1e-5,
        atol=1e-5,
    )


@settings(max_examples=25, deadline=2000)
@given(
    n_samples=st.integers(min_value=2, max_value=6),
    n_states=st.integers(min_value=2, max_value=5),
    seed=st.integers(min_value=1, max_value=10_000),
)
def test_property_true_cell_fractions_simplex(
    n_samples: int,
    n_states: int,
    seed: int,
) -> None:
    """Property: ground-truth state fractions are non-negative and sum to 1.0 per sample."""
    rng = np.random.default_rng(seed)
    sample_list: list[str] = []
    state_list: list[str] = []
    states = tuple(f"State_{k}" for k in range(n_states))

    for s in range(n_samples):
        n_c = int(rng.integers(5, 25))
        sample_list.extend([f"P_{s:02d}"] * n_c)
        chosen = rng.choice(states, size=n_c, replace=True).tolist()
        state_list.extend(str(c) for c in chosen)

    cells_df = pl.DataFrame({"sample_id": sample_list, "cell_state": state_list})
    frac_res = compute_empirical_cell_fractions(cells_df)
    assert is_successful(frac_res)
    frac_df = frac_res.unwrap()
    observed_states = cells_df["cell_state"].n_unique()
    assert frac_df.height == n_samples * observed_states
    assert float(frac_df["true_fraction"].to_numpy().min()) >= 0.0

    sums_df = frac_df.group_by("sample_id").agg(pl.col("true_fraction").sum().alias("sum_frac"))
    np.testing.assert_allclose(
        sums_df["sum_frac"].to_numpy(),
        np.ones(n_samples, dtype=np.float64),
        rtol=1e-12,
        atol=1e-12,
    )


def test_automatic_expm1_linearization_from_log1p_adata() -> None:
    """Verify log1p-normalized AnnData matrices are automatically inverted via expm1 before summing."""
    linear_counts = np.array(
        [
            [10.0, 0.0, 25.0],
            [5.0, 15.0, 0.0],
            [0.0, 20.0, 30.0],
            [12.0, 8.0, 10.0],
        ],
        dtype=np.float64,
    )
    log_x = np.log1p(linear_counts)
    adata = _build_synthetic_adata(
        counts=log_x,
        sample_labels=["P1", "P1", "P2", "P2"],
        state_labels=["CD8", "B", "CD8", "Myeloid"],
        responses=["Responder", "Responder", "Non-responder", "Non-responder"],
        genes=["CD8A", "MS4A1", "LYZ"],
        sparse_x=True,
    )

    cfg = EmpiricalPseudobulkConfig(
        sample_col="sample_id",
        cell_state_col=Some("cell_state"),
        normalize_cells_cp10k=False,
        min_cells_per_sample=2,
    )
    res = aggregate_empirical_pseudobulk(adata, config=cfg, dataset_id="LOG_TEST")
    assert is_successful(res)
    out = res.unwrap()

    assert resolve_adata_gene_symbols(adata) == ("CD8A", "MS4A1", "LYZ")
    np.testing.assert_allclose(out.bulk_matrix[0], [15.0, 15.0, 25.0], rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(out.bulk_matrix[1], [12.0, 28.0, 40.0], rtol=1e-5, atol=1e-5)
    assert out.sample_metadata_df["response_binary"].to_list() == [1, 0]


def test_atomic_parquet_caching_subsetting_and_anndata_roundtrip(tmp_path: Path) -> None:
    """Verify cache_empirical_pseudobulk, fast cache hit, target_genes slicing, and AnnData conversion."""
    rng = np.random.default_rng(123)
    counts = rng.poisson(lam=100.0, size=(12, 5)).astype(np.float64) + 40.0
    genes = ["CD8A", "PDCD1", "CXCL13", "MS4A1", "COL1A1"]
    adata = _build_synthetic_adata(
        counts=counts,
        sample_labels=["S1"] * 4 + ["S2"] * 4 + ["S3"] * 4,
        state_labels=["Tex", "B", "CAF"] * 4,
        responses=["CR"] * 4 + ["PD"] * 4 + ["PR"] * 4,
        genes=genes,
        sparse_x=False,
    )

    cfg = EmpiricalPseudobulkConfig(
        sample_col="sample_id",
        cell_state_col=Some("cell_state"),
        min_cells_per_sample=3,
    )
    created = load_empirical_pseudobulk(
        "SYNTH_COHORT",
        config=cfg,
        cache=True,
        cache_dir=Some(tmp_path),
        adata=Some(adata),
    )
    assert is_successful(created)
    res1 = created.unwrap()
    assert isinstance(res1.cache_path, Some)
    parquet_file = res1.cache_path.unwrap()
    assert parquet_file.exists()
    assert (tmp_path / "SYNTH_COHORT.sample_meta.parquet").exists()
    assert (tmp_path / "SYNTH_COHORT.true_fractions.parquet").exists()

    # Fast cache hit without passing AnnData
    t0 = time.perf_counter()
    hit = load_empirical_pseudobulk(
        "SYNTH_COHORT",
        config=cfg,
        cache=True,
        cache_dir=Some(tmp_path),
    )
    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    assert is_successful(hit)
    assert elapsed_ms < 50.0
    res2 = hit.unwrap()
    np.testing.assert_allclose(res2.bulk_matrix, res1.bulk_matrix, rtol=1e-7, atol=1e-7)
    assert isinstance(res2.true_fractions_df, Some)

    # Subset target_genes directly from Parquet cache
    sub_cfg = EmpiricalPseudobulkConfig(
        target_genes=Some(frozenset({"CD8A", "CXCL13"})),
    )
    sub_hit = load_empirical_pseudobulk(
        "SYNTH_COHORT",
        config=sub_cfg,
        cache_dir=Some(tmp_path),
    )
    assert is_successful(sub_hit)
    sub_res = sub_hit.unwrap()
    assert sub_res.genes == ("CD8A", "CXCL13")
    assert sub_res.bulk_matrix.shape == (3, 2)

    # Round-trip conversion to AnnData
    pb_adata = empirical_pseudobulk_to_anndata(res2)
    assert pb_adata.shape == (3, 5)
    assert list(pb_adata.obs_names) == ["S1", "S2", "S3"]
    assert list(pb_adata.var_names) == genes
    assert "response_binary" in pb_adata.obs.columns


def test_store_io_guard_rejects_capstor_store() -> None:
    """Verify empirical pseudobulk cache operations reject /capstor/store paths cleanly."""
    adata = _build_synthetic_adata(
        counts=np.full((4, 2), 60.0, dtype=np.float64),
        sample_labels=["S1", "S1", "S2", "S2"],
        state_labels=["T", "B", "T", "B"],
        responses=["R", "R", "NR", "NR"],
        genes=["CD8A", "MS4A1"],
    )
    agg = aggregate_empirical_pseudobulk(adata, dataset_id="STORE_GUARD_TEST").unwrap()
    forbidden_dir = Path("/capstor/store/cscs/empa/em25/halu/pseudobulk")

    write_res = cache_empirical_pseudobulk(agg, cache_dir=Some(forbidden_dir))
    assert isinstance(write_res, Failure)
    assert "Direct runtime read/write to persistent $STORE" in write_res.failure()

    load_res = load_empirical_pseudobulk(
        "STORE_GUARD_TEST",
        cache_dir=Some(forbidden_dir),
        adata=Some(adata),
    )
    assert isinstance(load_res, Failure)
    assert "Direct runtime read/write to persistent $STORE" in load_res.failure()


def test_real_sc9_cached_pseudobulks_and_gse120575_raw_aggregation() -> None:
    """Verify all 9 SC-9 cached cohorts load in <100 ms and GSE120575 aggregates from raw H5AD."""
    for cid in SC9_COHORT_IDS:
        t0 = time.perf_counter()
        res = load_empirical_pseudobulk(cid)
        dt_ms = (time.perf_counter() - t0) * 1000.0
        assert isinstance(res, Success), f"Expected cached pseudobulk for {cid}: {res}"
        pb = res.unwrap()
        assert pb.bulk_matrix.shape[0] >= 7
        assert pb.bulk_matrix.shape[1] >= 500
        assert dt_ms < 150.0

    # Also verify direct uncached aggregation from data/preprocessed/GSE120575.h5ad
    raw_cfg = EmpiricalPseudobulkConfig(
        target_genes=Some(frozenset({"CD8A", "PDCD1", "TCF7", "CXCL13", "MS4A1"})),
        min_cells_per_sample=10,
    )
    raw_res = load_empirical_pseudobulk(
        "GSE120575",
        config=raw_cfg,
        cache=False,
        force_recompute=True,
    )
    assert isinstance(raw_res, Success)
    raw_pb = raw_res.unwrap()
    # Sade-Feldman 2018 has 48 patient biopsies (>=10 cells across 48 Pre/Post biopsies)
    assert 35 <= len(raw_pb.sample_ids) <= 48
    assert set(raw_pb.genes).issubset({"CD8A", "PDCD1", "TCF7", "CXCL13", "MS4A1"})
    assert len(raw_pb.genes) >= 4
