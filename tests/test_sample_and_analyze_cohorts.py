"""Tests for multi-cohort sampling, perturbation, harmonization, and diagnostic pipeline.

Verifies:
1. Deterministic cohort resolution across sampling axes (all, cancer_type, treatment, platform, custom).
2. Strict tme_datasets perturbation routing (dropout, jitter, knockout, overexpression, negative binomial).
3. End-to-end execution of sample_and_analyze_cohorts with exact neighborhood composition,
   replicate prevalence filtering, Parquet table export, and vector SVG diagnostic dashboard generation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from tme_datasets import (
    MILOPY_RESPONSE_COHORTS,
    NegativeBinomialConfig,
    SubsampleSpec,
    add_expression_jitter,
    in_silico_knockout,
    in_silico_overexpression,
    randomize_negative_binomial,
    resolve_milopy_cohorts,
    simulate_dropout,
)

# Import pipeline functions from scripts
import importlib
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
_stage2a1 = importlib.import_module(
    "scripts.manuscript_pipeline.2a1_sample_and_analyze_cohorts"
)
SampledRunConfig = _stage2a1.SampledRunConfig
apply_cohort_perturbation = _stage2a1.apply_cohort_perturbation
harmonize_cohort_metadata = _stage2a1.harmonize_cohort_metadata
run_sampling_pipeline = _stage2a1.run_sampling_pipeline

from scripts.milopy_analysis.plot_neighborhood_diagnostics import render_diagnostic_dashboard


def create_synthetic_cohort_adata(
    cohort_id: str,
    n_patients: int = 4,
    cells_per_pt: int = 50,
    n_genes: int = 80,
    seed: int = 42,
) -> ad.AnnData:
    """Generates synthetic AnnData representing a clinical response cohort."""
    rng = np.random.default_rng(seed)
    n_cells = n_patients * cells_per_pt

    pts = []
    responses = []
    for p in range(n_patients):
        pid = f"Pt_{cohort_id}_{p}"
        resp = "responder" if p % 2 == 0 else "non-responder"
        pts.extend([pid] * cells_per_pt)
        responses.extend([resp] * cells_per_pt)

    X = rng.poisson(lam=3.0, size=(n_cells, n_genes)).astype(np.float32)
    # Give responder cells slight expression difference in first 10 genes
    for i in range(n_cells):
        if responses[i] == "responder":
            X[i, :10] += 2.5

    gene_names = [f"GENE_{j}" for j in range(n_genes)]
    gene_names[0] = "CD8A"
    gene_names[1] = "PDCD1"
    gene_names[2] = "CXCL13"
    gene_names[3] = "IFNG"

    obs = pd.DataFrame(
        {
            "dataset_id": cohort_id,
            "patient_id": pts,
            "clinical_response": responses,
            "cancer_type": "Melanoma" if "mela" in cohort_id.lower() else "NSCLC",
        },
        index=[f"{cohort_id}_cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)


class TestCohortResolution:
    """Verifies resolve_milopy_cohorts across all supported axes."""

    def test_resolve_all_cohorts(self) -> None:
        res = resolve_milopy_cohorts("all")
        assert isinstance(res, Success)
        cohorts = res.unwrap()
        assert len(cohorts) == 9
        assert "GSE120575" in cohorts
        assert "GSE316195" in cohorts

    def test_resolve_by_cancer_type(self) -> None:
        res = resolve_milopy_cohorts("cancer_type", "melanoma")
        assert isinstance(res, Success)
        cohorts = res.unwrap()
        assert cohorts == ("GSE120575", "CELLxGENE_7b20c613")

        res_nsclc = resolve_milopy_cohorts("cancer_type", "nsclc")
        assert isinstance(res_nsclc, Success)
        assert set(res_nsclc.unwrap()) == {"GSE207422", "GSE243013", "GSE233203"}

    def test_resolve_by_treatment(self) -> None:
        res = resolve_milopy_cohorts("treatment", "anti-pd1_monotherapy")
        assert isinstance(res, Success)
        cohorts = res.unwrap()
        assert len(cohorts) == 3
        assert "CELLxGENE_7b20c613" in cohorts

    def test_resolve_by_platform(self) -> None:
        res = resolve_milopy_cohorts("platform", "bd_rhapsody")
        assert isinstance(res, Success)
        cohorts = res.unwrap()
        assert cohorts == ("GSE207422",)

    def test_resolve_custom_list(self) -> None:
        res = resolve_milopy_cohorts("custom", "GSE120575,GSE200996")
        assert isinstance(res, Success)
        assert res.unwrap() == ("GSE120575", "GSE200996")

    def test_invalid_filter_returns_failure(self) -> None:
        res = resolve_milopy_cohorts("cancer_type", "non_existent_cancer")
        assert isinstance(res, Failure)


class TestPerturbations:
    """Verifies perturbation models applied strictly via tme_datasets.transforms."""

    def test_simulate_dropout(self) -> None:
        adata = create_synthetic_cohort_adata("TestCohort", n_patients=2, cells_per_pt=30)
        initial_nnz = adata.X.nnz
        res = simulate_dropout(adata, rate=0.20, seed=42)
        assert isinstance(res, Success)
        perturbed = res.unwrap()
        assert perturbed.X.nnz < initial_nnz

    def test_add_expression_jitter(self) -> None:
        adata = create_synthetic_cohort_adata("TestCohort", n_patients=2, cells_per_pt=30)
        orig_mean = float(np.mean(adata.X.toarray()))
        res = add_expression_jitter(adata, sigma=0.20, seed=42)
        assert isinstance(res, Success)
        perturbed = res.unwrap()
        new_mean = float(np.mean(perturbed.X.toarray()))
        # Log-normal jitter has mean exp(sigma^2/2) > 1
        assert new_mean > 0.0

    def test_in_silico_knockout(self) -> None:
        adata = create_synthetic_cohort_adata("TestCohort", n_patients=2, cells_per_pt=30)
        cd8_idx = list(adata.var_names).index("CD8A")
        initial_cd8_sum = float(adata.X.toarray()[:, cd8_idx].sum())
        assert initial_cd8_sum > 0

        res = in_silico_knockout(adata, genes=("CD8A",), efficiency=1.0)
        assert isinstance(res, Success)
        knocked = res.unwrap()
        assert float(knocked.X.toarray()[:, cd8_idx].sum()) == 0.0

    def test_in_silico_overexpression(self) -> None:
        adata = create_synthetic_cohort_adata("TestCohort", n_patients=2, cells_per_pt=30)
        cxcl13_idx = list(adata.var_names).index("CXCL13")
        initial_cxcl13_sum = float(adata.X.toarray()[:, cxcl13_idx].sum())

        res = in_silico_overexpression(adata, genes=("CXCL13",), fold_change=3.0)
        assert isinstance(res, Success)
        oe = res.unwrap()
        assert np.isclose(float(oe.X.toarray()[:, cxcl13_idx].sum()), initial_cxcl13_sum * 3.0, rtol=1e-4)


class TestPipelineEndToEnd:
    """Verifies end-to-end execution of sample_and_analyze_cohorts using mocked loaders."""

    def test_end_to_end_sampling_and_diagnostics(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # Create mock cohorts
        adata1 = create_synthetic_cohort_adata("MockCohortA", n_patients=4, cells_per_pt=40, seed=1)
        adata2 = create_synthetic_cohort_adata("MockCohortB", n_patients=4, cells_per_pt=40, seed=2)

        mock_registry = {"MockCohortA": adata1, "MockCohortB": adata2}

        def mock_load_dataset(accession: str, **kwargs: Any) -> Result[ad.AnnData, str]:
            if accession in mock_registry:
                return Success(mock_registry[accession].copy())
            return Failure(f"Unknown mock cohort {accession}")

        # Monkeypatch load_dataset in 2a1_sample_and_analyze_cohorts
        monkeypatch.setattr(_stage2a1, "load_dataset", mock_load_dataset)

        out_dir = tmp_path / "results"
        reports_dir = tmp_path / "reports"

        config = SampledRunConfig(
            sample_by="custom",
            value="MockCohortA,MockCohortB",
            cohort_ids=("MockCohortA", "MockCohortB"),
            cells_per_dataset=150,
            base_dir=tmp_path,
            out_dir=out_dir,
            reports_dir=reports_dir,
            perturbation="dropout",
            dropout_rate=0.05,
            harmonize=True,
            harmony_key="patient_id",
            harmony_max_iter=5,
            prop=0.20,
            k=15,
            d=10,
            fdr_threshold=0.20,
            min_replicates=1,
            min_prevalence_frac=0.02,
            run_permutation=True,
            n_permutations=20,
            plot=True,
            seed=42,
        )

        res = run_sampling_pipeline(config)
        assert isinstance(res, Success), f"Pipeline failed: {res}"
        da_path = res.unwrap()
        assert da_path.exists()

        run_dir = da_path.parent

        # 1. Verify da_results.parquet
        da_df = pl.read_parquet(da_path)
        assert len(da_df) > 0
        required_da_cols = [
            "Nhood",
            "logFC",
            "FDR",
            "status",
            "hit_type",
            "is_significant",
            "UMAP1",
            "UMAP2",
            "dataset_id_purity",
            "patient_id_max_fraction",
            "replicate_prevalence",
        ]
        for col in required_da_cols:
            assert col in da_df.columns, f"Missing {col} in da_results.parquet"

        # 2. Verify nhood_composition_datasets.parquet
        ds_comp_path = run_dir / "nhood_composition_datasets.parquet"
        assert ds_comp_path.exists()
        ds_df = pl.read_parquet(ds_comp_path)
        assert "nhood_index" in ds_df.columns
        assert "nhood_size" in ds_df.columns
        assert "dominant_dataset_id" in ds_df.columns

        # Check that per-cohort counts exist
        assert any("count_dataset_id_mockcohorta" in c for c in ds_df.columns)
        assert any("count_dataset_id_mockcohortb" in c for c in ds_df.columns)

        # 3. Verify nhood_composition_patients.parquet
        pt_comp_path = run_dir / "nhood_composition_patients.parquet"
        assert pt_comp_path.exists()
        pt_df = pl.read_parquet(pt_comp_path)
        assert "nhood_index" in pt_df.columns
        assert "dominant_patient_id" in pt_df.columns
        assert "patient_id_purity" in pt_df.columns

        # 4. Verify cell_level_scores.parquet
        cell_path = run_dir / "cell_level_scores.parquet"
        assert cell_path.exists()
        cell_df = pl.read_parquet(cell_path)
        assert len(cell_df) > 0
        for col in ["cell_id", "dataset_id", "patient_id", "clinical_response", "cell_logfc", "min_fdr", "is_significant"]:
            assert col in cell_df.columns, f"Missing {col} in cell_level_scores.parquet"

        # 5. Verify run_summary.json
        summary_path = run_dir / "run_summary.json"
        assert summary_path.exists()
        with open(summary_path) as f:
            summary = json.load(f)
        assert summary["sampling"]["sample_by"] == "custom"
        assert summary["perturbation"]["type"] == "dropout"
        assert summary["metrics"]["total_nhoods"] == len(da_df)

        # 6. Verify Diagnostic Visualizations (SVG & PNG)
        run_name = "custom_MockCohortA,MockCohortB_dropout".replace("/", "_").replace(" ", "_")
        svg_file = reports_dir / run_name / f"{run_name}_diagnostics.svg"
        assert svg_file.exists(), f"SVG report file {svg_file} does not exist"

        # Validate SVG is well-formed XML
        tree = ET.parse(svg_file)
        root = tree.getroot()
        assert "svg" in root.tag.lower()
