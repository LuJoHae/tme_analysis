"""Unit tests for single-cell deconvolution reference construction pipeline."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets.deconvolution import (
    DeconvolutionReferenceConfig,
    DeconvolutionReferenceResult,
    build_deconvolution_reference,
    calculate_cnv_proxy_scores,
    detect_malignant_cells,
    export_to_bayesprism,
    export_to_instaprism,
)


@pytest.fixture
def mock_sc_anndata() -> ad.AnnData:
    """Construct a synthetic multi-lineage single-cell dataset with known properties."""
    rng = np.random.default_rng(42)
    n_genes = 20
    gene_names = [f"GENE_{i}" for i in range(15)] + [
        "MT-CO1",    # Confounding: Mitochondrial
        "RPS18",     # Confounding: Ribosomal
        "TRAV1-1",   # Confounding: TCR
        "MLANA",     # Confounding: Melanoma marker
        "EPCAM",     # Epithelial marker
    ]

    # Create 5 distinct cell states across 3 lineages
    states_def = [
        ("CD8_Effector", "T_cell", False, 30, 2.0),
        ("CD8_Exhausted", "T_cell", False, 25, 2.1),   # Highly correlated with CD8_Effector
        ("Treg", "T_cell", False, 20, 1.8),
        ("Macrophage_M1", "Myeloid", False, 25, 5.0),  # Larger library size
        ("Macrophage_M2", "Myeloid", False, 25, 4.8),  # Correlated with M1
        ("Melanoma_Tumor", "Malignant", True, 40, 10.0), # High mRNA content
        ("Rare_Contaminant", "Unknown", False, 3, 1.0),  # < 15 cells, should be filtered
    ]

    all_counts = []
    obs_states = []
    obs_types = []
    obs_malignant = []

    for state, cell_type, is_mal, n_cells, lib_factor in states_def:
        # Base expression profile per state
        base_lambda = rng.uniform(0.5, 5.0, size=n_genes) * lib_factor
        if is_mal:
            # Overexpress tumor marker MLANA
            base_lambda[gene_names.index("MLANA")] = 50.0

        counts = rng.poisson(lam=base_lambda, size=(n_cells, n_genes)).astype(np.float32)
        all_counts.append(counts)
        obs_states.extend([state] * n_cells)
        obs_types.extend([cell_type] * n_cells)
        obs_malignant.extend([is_mal] * n_cells)

    X_mat = np.vstack(all_counts)
    total_cells = X_mat.shape[0]

    obs_df = pd.DataFrame(
        {
            "cell_state": obs_states,
            "cell_type": obs_types,
            "is_malignant": obs_malignant,
        },
        index=[f"cell_{i}" for i in range(total_cells)],
    )

    var_df = pd.DataFrame(
        {"gene_symbol": gene_names, "contig": ["1"] * 10 + ["2"] * 10},
        index=gene_names,
    )

    return ad.AnnData(X=X_mat, obs=obs_df, var=var_df)


def test_deconvolution_reference_config_defaults() -> None:
    """Verify frozen configuration defaults."""
    cfg = DeconvolutionReferenceConfig()
    assert cfg.cell_state_key == "cell_state"
    assert cfg.min_cells_per_state == 15
    assert cfg.collinearity_threshold == 0.85
    assert cfg.pseudo_min == 1e-8
    assert cfg.filter_confounding is True
    assert cfg.normalize_multinomial is True


def test_build_deconvolution_reference_pipeline(mock_sc_anndata: ad.AnnData) -> None:
    """Verify end-to-end deconvolution reference building from AnnData."""
    config = DeconvolutionReferenceConfig(
        cell_state_key="cell_state",
        cell_type_key=Some("cell_type"),
        malignant_key=Some("is_malignant"),
        min_cells_per_state=15,
        filter_confounding=True,
    )

    res = build_deconvolution_reference(mock_sc_anndata, config)
    assert isinstance(res, Success)
    ref = res.unwrap()
    assert isinstance(ref, DeconvolutionReferenceResult)

    # 1. Rare state (< 15 cells) should be filtered out
    assert "Rare_Contaminant" not in ref.cell_states
    assert len(ref.cell_states) == 6
    assert len(ref.cell_types) == 3  # T_cell, Myeloid, Malignant

    # 2. Confounding genes filtered (MT-CO1, RPS18, TRAV1-1, MLANA should be stripped)
    for banned in ["MT-CO1", "RPS18", "TRAV1-1", "MLANA"]:
        assert banned not in ref.gene_names
        assert banned not in ref.phi_state.columns

    # 3. Multinomial row normalization: each row must sum to 1.0
    mat_state, states, genes = ref.to_numpy(which="state")
    assert mat_state.shape == (6, len(ref.gene_names))
    row_sums = mat_state.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, rtol=1e-5)

    mat_type, types, _ = ref.to_numpy(which="type")
    assert mat_type.shape == (3, len(ref.gene_names))
    np.testing.assert_allclose(mat_type.sum(axis=1), 1.0, rtol=1e-5)

    # 4. Hierarchy table checks
    assert isinstance(ref.hierarchy_table, pl.DataFrame)
    assert "cell_state" in ref.hierarchy_table.columns
    assert "cell_type" in ref.hierarchy_table.columns
    assert "is_malignant" in ref.hierarchy_table.columns
    assert "n_cells" in ref.hierarchy_table.columns

    # Verify Malignant tagging
    mal_rows = ref.hierarchy_table.filter(pl.col("is_malignant"))
    assert mal_rows["cell_state"].to_list() == ["Melanoma_Tumor"]
    assert ref.malignant_states == ("Melanoma_Tumor",)

    # 5. Relative mRNA scaling factors (Tumor cells should have highest relative content)
    mrna_df = ref.mrna_scaling
    assert isinstance(mrna_df, pl.DataFrame)
    tumor_scaling = mrna_df.filter(pl.col("cell_state") == "Melanoma_Tumor")["relative_rna_content"][0]
    treg_scaling = mrna_df.filter(pl.col("cell_state") == "Treg")["relative_rna_content"][0]
    assert tumor_scaling > treg_scaling

    # 6. Collinearity table
    coll_df = ref.collinearity
    assert isinstance(coll_df, pl.DataFrame)
    assert "state_a" in coll_df.columns
    assert "state_b" in coll_df.columns
    assert "correlation" in coll_df.columns
    assert "is_collinear" in coll_df.columns
    assert len(coll_df) == (6 * 5) // 2  # 15 pairwise comparisons


def test_adapters_instaprism_and_bayesprism(mock_sc_anndata: ad.AnnData) -> None:
    """Verify export adapters to InstaPrism and BayesPrism."""
    res = build_deconvolution_reference(mock_sc_anndata)
    assert isinstance(res, Success)
    ref = res.unwrap()

    # InstaPrism adapter
    mat_insta, labels_insta = ref.to_instaprism()
    assert isinstance(mat_insta, np.ndarray)
    assert len(labels_insta) == 6
    assert mat_insta.shape[0] == 6

    # BayesPrism adapter
    bp_dict_res = export_to_bayesprism(ref)
    assert isinstance(bp_dict_res, Success)
    bp_dict = bp_dict_res.unwrap()
    assert "phi_cell_state" in bp_dict
    assert "phi_cell_type" in bp_dict
    assert "state_to_type_map" in bp_dict

    # BayesPrism adapter with mock bulk mixture -> Prism instance
    rng = np.random.default_rng(123)
    mock_bulk = rng.poisson(lam=10.0, size=(5, len(ref.gene_names))).astype(np.float32)
    prism_res = export_to_bayesprism(ref, mixture=mock_bulk, bulk_names=["S1", "S2", "S3", "S4", "S5"])
    assert isinstance(prism_res, Success)
    prism_obj = prism_res.unwrap()
    assert prism_obj.mixture.shape == (5, len(ref.gene_names))
    assert len(prism_obj.bulk_names) == 5


def test_export_parquet_and_anndata(mock_sc_anndata: ad.AnnData, tmp_path: Path) -> None:
    """Verify Parquet export and AnnData conversion."""
    res = build_deconvolution_reference(mock_sc_anndata)
    assert isinstance(res, Success)
    ref = res.unwrap()

    # Parquet export
    out_dir = tmp_path / "reference_parquets"
    export_res = ref.export_parquet(out_dir)
    assert isinstance(export_res, Success)
    assert (out_dir / "curated_reference_phi.parquet").exists()
    assert (out_dir / "reference_phi_cell_type.parquet").exists()
    assert (out_dir / "cell_hierarchy.parquet").exists()
    assert (out_dir / "mrna_scaling.parquet").exists()
    assert (out_dir / "reference_collinearity.parquet").exists()

    # Read back parquet with Polars
    read_phi = pl.read_parquet(out_dir / "curated_reference_phi.parquet")
    assert len(read_phi) == 6

    # AnnData export
    ref_adata = ref.to_anndata()
    assert isinstance(ref_adata, ad.AnnData)
    assert ref_adata.shape == (6, len(ref.gene_names))
    assert "mrna_scaling" in ref_adata.uns
    assert "collinearity" in ref_adata.uns


def test_auto_detect_malignant_heuristic(mock_sc_anndata: ad.AnnData) -> None:
    """Verify automatic detection of malignant cells when no explicit key is provided."""
    # Remove explicit is_malignant metadata
    adata_clean = mock_sc_anndata.copy()
    del adata_clean.obs["is_malignant"]

    # Overwrite cell_state names to not contain "Tumor" or "Malignant"
    renamed_states = [
        "Clone_X" if s == "Melanoma_Tumor" else s
        for s in adata_clean.obs["cell_state"]
    ]
    adata_clean.obs["cell_state"] = renamed_states

    config = DeconvolutionReferenceConfig(
        cell_state_key="cell_state",
        malignant_key=Nothing,
        auto_detect_malignant=True,
    )

    detected = detect_malignant_cells(adata_clean, config)
    assert isinstance(detected, np.ndarray)
    assert detected.shape == (adata_clean.n_obs,)
    # Clone_X has massive MLANA expression, should be flagged
    assert np.any(detected)


def test_missing_state_key_returns_failure(mock_sc_anndata: ad.AnnData) -> None:
    """Verify Failure result if required cell_state key is missing."""
    config = DeconvolutionReferenceConfig(cell_state_key="non_existent_column_123")
    res = build_deconvolution_reference(mock_sc_anndata, config)
    assert isinstance(res, Failure)
    assert "not found in AnnData.obs" in res.failure()
