"""Unit, property-based (Hypothesis), R XDR, and HDF5 streaming tests for T-Cell & CAF Atlases."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import pytest
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
import scipy.sparse as sp

import tme_datasets
from tme_datasets.providers.caf_atlases import (
    BUECHLER_CAF_CHECKSUM,
    BUECHLER_CAF_CLUSTER_SPECS,
    BUECHLER_CAF_DATASET_ID,
    CORDS_CAF_CLUSTER_SPECS,
    CORDS_CAF_DATASET_ID,
    CORDS_CAF_ZIP_CHECKSUM,
    CafAtlasChoice,
    GAO_CAF_CLUSTER_SPECS,
    GAO_CAF_DATASET_ID,
    GAO_CAF_EXPR_CHECKSUM,
    dispatch_caf_atlas_load,
    download_buechler_caf_atlas,
    download_gao_caf_atlas,
    harmonize_buechler_caf_obs_polars,
    harmonize_cords_caf_obs_polars,
    harmonize_gao_caf_obs_polars,
    preprocess_gao_caf_atlas,
)
from tme_datasets.providers.tcell_atlases import (
    CHU_TCELL_CLUSTER_SPECS,
    CHU_TCELL_DATASET_ID,
    CHU_TCELL_HS_MARKERS,
    SHI_PANCANCER_CHECKSUM,
    SHI_PANCANCER_DATASET_ID,
    ShiCompartmentFilter,
    TCellAtlasChoice,
    ZHENG_TCELL_CLUSTER_SPECS,
    ZHENG_TCELL_DATASET_ID,
    ZHENG_TCELL_EXPR_CHECKSUM,
    dispatch_tcell_atlas_load,
    download_chu_tcell_atlas,
    download_shi_pancancer_atlas,
    download_zheng_tcell_atlas,
    harmonize_chu_tcell_obs_polars,
    harmonize_shi_pancancer_obs_polars,
    harmonize_zheng_tcell_obs_polars,
    preprocess_shi_pancancer_atlas,
    select_shi_compartment_mask,
)
from tme_datasets.query import load_dataset
from tme_datasets.registry import (
    get_dataset_spec,
    list_registered_datasets,
    resolve_dataset_id,
)
from tme_datasets.storage import (
    read_r_dataframe,
    stream_r_object_to_staging_h5,
)
from tme_datasets.types import Modality


# =============================================================================
# 1. Registry, Alias Resolution, GSE154763 Metadata Fix & Top-Level Exports
# =============================================================================


def test_six_caf_and_tcell_atlases_registered_and_aliases() -> None:
    """Verify all 6 pan-cancer CAF and T-cell atlases and their aliases in the registry."""
    registered_ids = list_registered_datasets().ids()
    expected_specs = (
        (
            GAO_CAF_DATASET_ID,
            269899,
            "39303725",
            "10.1016/j.ccell.2024.08.020",
            Some(GAO_CAF_EXPR_CHECKSUM),
            ("CAF_Atlas", "Gao2024", "Gao_2024", "Gao_Fibroblast_Atlas", "Zenodo_13950610", "13950610"),
        ),
        (
            CORDS_CAF_DATASET_ID,
            29741,
            "38242124",
            "10.1016/j.ccell.2023.12.021",
            Some(CORDS_CAF_ZIP_CHECKSUM),
            ("Cords2024", "Cords_2024", "Cords2023", "Zenodo_7540604", "7540604"),
        ),
        (
            BUECHLER_CAF_DATASET_ID,
            10355,
            "33981032",
            "10.1038/s41586-021-03549-5",
            Some(BUECHLER_CAF_CHECKSUM),
            ("FibroXplorer", "Buechler2021", "Buechler_2021", "Human_PS_Fibro", "Human_PS_Fibro.RDS"),
        ),
        (
            ZHENG_TCELL_DATASET_ID,
            397810,
            "34914499",
            "10.1126/science.abe6474",
            Some(ZHENG_TCELL_EXPR_CHECKSUM),
            ("TCell_Atlas", "GSE156728", "Zheng2021", "Zheng_2021", "Zenodo_5461803", "5461803"),
        ),
        (
            CHU_TCELL_DATASET_ID,
            308048,
            "37248301",
            "10.1038/s41591-023-02371-y",
            Nothing,
            ("Chu2023", "Chu_2023", "TCM_Atlas", "MDAnderson_TCM", "panC_T_all"),
        ),
        (
            SHI_PANCANCER_DATASET_ID,
            1200000,
            "40437094",
            "10.1038/s41586-025-09053-4",
            Some(SHI_PANCANCER_CHECKSUM),
            ("Shi2025", "Shi_2025", "Zenodo_15554080", "15554080", "PanCancer_igt_s9_fine_counts"),
        ),
    )

    for canonical_id, n_cells, pmid, doi, expected_csum, aliases in expected_specs:
        assert canonical_id in registered_ids
        spec_maybe = get_dataset_spec(canonical_id)
        assert isinstance(spec_maybe, Some)
        spec = spec_maybe.unwrap()
        assert spec.id == canonical_id
        assert spec.modality == Modality.SINGLE_CELL
        assert spec.n_samples_or_cells == Some(n_cells)
        assert spec.publication_pmid == Some(pmid)
        assert spec.publication_doi == Some(doi)
        assert spec.checksum == expected_csum
        assert spec.qc_spec is Nothing
        assert spec.is_dysfunctional is False

        for alias in aliases:
            assert resolve_dataset_id(alias) == canonical_id
            assert get_dataset_spec(alias) == Some(spec)

    # Verify GSE154763 title correction (Cheng et al. 2021 Cell Myeloid Atlas)
    gse154763 = get_dataset_spec("GSE154763").unwrap()
    assert "Cheng et al. 2021" in gse154763.title
    assert "Myeloid" in gse154763.title
    assert gse154763.is_dysfunctional is True

    # Verify curated ontology cardinalities and top-level package exports
    assert len(GAO_CAF_CLUSTER_SPECS) == 20
    assert len(CORDS_CAF_CLUSTER_SPECS) == 11
    assert len(BUECHLER_CAF_CLUSTER_SPECS) == 6
    assert len(ZHENG_TCELL_CLUSTER_SPECS) == 41
    assert len(CHU_TCELL_CLUSTER_SPECS) == 39
    assert len(CHU_TCELL_HS_MARKERS) == 15
    for sym in (
        "load_gao_caf_atlas",
        "load_cords_caf_atlas",
        "load_buechler_caf_atlas",
        "load_zheng_tcell_atlas",
        "load_chu_tcell_atlas",
        "load_shi_pancancer_atlas",
        "read_r_dataframe",
        "stream_r_object_to_staging_h5",
    ):
        assert hasattr(tme_datasets, sym)


# =============================================================================
# 2. Hypothesis Property-Based Tests for Polars `.obs` Harmonizers
# =============================================================================


@settings(max_examples=30, deadline=8000)
@given(
    cluster_keys=st.lists(
        st.sampled_from(sorted(GAO_CAF_CLUSTER_SPECS.keys())),
        min_size=1,
        max_size=20,
    ),
    tissues=st.lists(
        st.sampled_from(["Lung", "Pancreas", "Breast", "Colon", "Liver", "Skin", "Ovary"]),
        min_size=1,
        max_size=20,
    ),
    phenos=st.lists(
        st.sampled_from(["Tumor", "Normal", "Healthy", "Inflamed", "Fetal"]),
        min_size=1,
        max_size=20,
    ),
)
def test_property_harmonize_gao_caf_obs_polars(
    cluster_keys: list[str],
    tissues: list[str],
    phenos: list[str],
) -> None:
    """Property: Gao CAF Polars harmonizer preserves rows, all raw columns, UMAP, and c01..c20."""
    n = len(cluster_keys)
    barcodes = [f"gao_cell_{i}" for i in range(n)]
    raw = pd.DataFrame(
        {
            "Cluster": cluster_keys,
            "Tissue": [tissues[i % len(tissues)] for i in range(n)],
            "Phenotype": [phenos[i % len(phenos)] for i in range(n)],
            "Dataset": [f"Study_{i % 4}" for i in range(n)],
            "Patient": [f"Pat_{i % 5}" for i in range(n)],
            "Donor_Age_Custom": [40 + (i % 35) for i in range(n)],
            "UMAP_1": [float(i) * 0.5 for i in range(n)],
            "UMAP_2": [float(-i) * 0.25 for i in range(n)],
        },
        index=barcodes,
    )
    out, umap_maybe = harmonize_gao_caf_obs_polars(raw)
    assert len(out) == n
    assert list(out.index) == barcodes
    assert isinstance(umap_maybe, Some)
    assert umap_maybe.unwrap().shape == (n, 2)
    for col in raw.columns:
        assert col in out.columns
    assert list(out["dataset_id"]) == [GAO_CAF_DATASET_ID] * n
    assert list(out["lineage_compartment"]) == ["STROMAL_CAF"] * n
    for idx, ck in enumerate(cluster_keys):
        spec = GAO_CAF_CLUSTER_SPECS[ck]
        assert out["cell_state"].iloc[idx] == spec.subtype_label
        assert out["cell_type"].iloc[idx] == spec.broad_cell_type
        assert out["lineage_branch"].iloc[idx] == spec.lineage_branch


@settings(max_examples=30, deadline=8000)
@given(
    cluster_keys=st.lists(
        st.sampled_from(sorted(CORDS_CAF_CLUSTER_SPECS.keys())),
        min_size=1,
        max_size=20,
    ),
)
def test_property_harmonize_cords_caf_obs_polars(
    cluster_keys: list[str],
) -> None:
    """Property: Cords CAF Polars harmonizer preserves informative columns and maps CAF phenotypes."""
    n = len(cluster_keys)
    raw = pd.DataFrame(
        {
            "CAFtype": cluster_keys,
            "Patient": [f"P_{i % 5}" for i in range(n)],
            "custom_author_metric": [float(i) * 0.25 for i in range(n)],
        },
        index=[f"cords_{i}" for i in range(n)],
    )
    out = harmonize_cords_caf_obs_polars(
        raw,
        cohort_id="Breast_Cords",
        cancer_code="BC",
        organ="Breast",
    )
    assert len(out) == n
    assert "CAFtype" in out.columns
    assert "custom_author_metric" in out.columns
    assert list(out["lineage_compartment"]) == ["STROMAL_CAF"] * n
    for idx, ck in enumerate(cluster_keys):
        spec = CORDS_CAF_CLUSTER_SPECS[ck]
        assert out["cell_state"].iloc[idx] == spec.subtype_label
        assert out["cell_type"].iloc[idx] == spec.broad_cell_type


@settings(max_examples=30, deadline=8000)
@given(
    cluster_keys=st.lists(
        st.sampled_from(sorted(BUECHLER_CAF_CLUSTER_SPECS.keys())),
        min_size=1,
        max_size=20,
    ),
    tissues=st.lists(
        st.sampled_from(["NSCLC", "PDAC", "IPF", "COVID", "Ulcerative Colitis"]),
        min_size=1,
        max_size=20,
    ),
)
def test_property_harmonize_buechler_caf_obs_polars(
    cluster_keys: list[str],
    tissues: list[str],
) -> None:
    """Property: Buechler FibroXplorer harmonizer preserves all columns and maps universal/perturbed states."""
    n = len(cluster_keys)
    raw = pd.DataFrame(
        {
            "ClustName": cluster_keys,
            "Tissue": [tissues[i % len(tissues)] for i in range(n)],
            "ExtraAuthorCol": [f"meta_{i}" for i in range(n)],
        },
        index=[f"S{i % 3}_buechler{i}" for i in range(n)],
    )
    out = harmonize_buechler_caf_obs_polars(raw)
    assert len(out) == n
    for col in raw.columns:
        assert col in out.columns
    assert list(out["lineage_compartment"]) == ["STROMAL_CAF"] * n
    for idx, ck in enumerate(cluster_keys):
        spec = BUECHLER_CAF_CLUSTER_SPECS[ck]
        assert out["cell_state"].iloc[idx] == spec.subtype_label
        assert out["lineage_branch"].iloc[idx] == spec.lineage_branch


@settings(max_examples=30, deadline=8000)
@given(
    cluster_keys=st.lists(
        st.sampled_from(sorted(ZHENG_TCELL_CLUSTER_SPECS.keys())),
        min_size=1,
        max_size=20,
    ),
    locs=st.lists(
        st.sampled_from(["T", "N", "P", "L"]),
        min_size=1,
        max_size=20,
    ),
)
def test_property_harmonize_zheng_tcell_obs_polars(
    cluster_keys: list[str],
    locs: list[str],
) -> None:
    """Property: Zheng 2021 T-cell harmonizer preserves TCR/pseudotime/clinical columns and maps 41 states."""
    n = len(cluster_keys)
    raw = pd.DataFrame(
        {
            "cellID.uniq": [f"zheng_cell_{i}" for i in range(n)],
            "meta.cluster": cluster_keys,
            "loc": [locs[i % len(locs)] for i in range(n)],
            "cancerType": ["MEL" if i % 2 == 0 else "NSCLC" for i in range(n)],
            "patient": [f"P316_{i % 6}" for i in range(n)],
            "clone_id": [f"clonotype_{i % 3}" for i in range(n)],
            "clone_size": [1 + (i % 12) for i in range(n)],
            "ptrt_nes": [float(i % 5) * 0.6 for i in range(n)],
            "dpt_pseudotime": [float(i) / max(1, n) for i in range(n)],
            "patient_stage": ["Stage_III" for _ in range(n)],
        },
        index=[f"zheng_cell_{i}" for i in range(n)],
    )
    out = harmonize_zheng_tcell_obs_polars(raw)
    assert len(out) == n
    for col in raw.columns:
        assert col in out.columns
    for idx, ck in enumerate(cluster_keys):
        spec = ZHENG_TCELL_CLUSTER_SPECS[ck]
        assert out["cell_state"].iloc[idx] == ck
        assert out["cell_type"].iloc[idx] == spec.broad_cell_type
        assert out["lineage_compartment"].iloc[idx] == spec.lineage_compartment
        assert out["lineage_branch"].iloc[idx] == spec.lineage_branch
        assert out["trajectory_group"].iloc[idx] == spec.trajectory_group


@settings(max_examples=30, deadline=8000)
@given(
    cluster_keys=st.lists(
        st.sampled_from(sorted(CHU_TCELL_CLUSTER_SPECS.keys())),
        min_size=1,
        max_size=20,
    ),
)
def test_property_harmonize_chu_tcell_obs_polars(
    cluster_keys: list[str],
) -> None:
    """Property: Chu 2023 T-cell harmonizer preserves columns and identifies CD8_c4_Tstr / CD4_c4_Tstr."""
    n = len(cluster_keys)
    raw = pd.DataFrame(
        {
            "cell.type": cluster_keys,
            "cancer": ["LUAD" if i % 2 == 0 else "BRCA" for i in range(n)],
            "Response": ["Responder" if i % 3 == 0 else "Non-responder" for i in range(n)],
        },
        index=[f"chu_{i}" for i in range(n)],
    )
    out = harmonize_chu_tcell_obs_polars(raw)
    assert len(out) == n
    for col in raw.columns:
        assert col in out.columns
    assert list(out["lineage_compartment"]) == ["T_CELL"] * n
    for idx, ck in enumerate(cluster_keys):
        spec = CHU_TCELL_CLUSTER_SPECS[ck]
        assert out["cell_state"].iloc[idx] == ck
        assert out["cell_type"].iloc[idx] == spec.broad_cell_type
        assert bool(out["is_stress_tstr"].iloc[idx]) == (ck in ("CD8_c4_Tstr", "CD4_c4_Tstr"))


@settings(max_examples=30, deadline=8000)
@given(
    major_labels=st.lists(
        st.sampled_from(
            [
                "CD8 T cell",
                "CD4 Treg",
                "NK cell",
                "Cancer-associated fibroblast",
                "Pericyte",
                "Macrophage",
                "B cell",
                "Epithelial malignant",
            ]
        ),
        min_size=1,
        max_size=25,
    ),
)
def test_property_shi_compartment_mask_boolean_algebra(
    major_labels: list[str],
) -> None:
    """Property: Shi compartment mask satisfies `T_AND_CAF == (T_NK_ONLY | STROMAL_CAF_ONLY)` and `ALL_CELLS == ones`."""
    n = len(major_labels)
    raw = pd.DataFrame(
        {
            "cell_type_major": major_labels,
            "cell_type_fine": [f"fine_{lbl}" for lbl in major_labels],
            "cancer_type": ["CRC" for _ in range(n)],
            "cm_module": [f"CM{(i % 12) + 1:02d}" for i in range(n)],
        },
        index=[f"shi_{i}" for i in range(n)],
    )
    harm = harmonize_shi_pancancer_obs_polars(raw)
    assert len(harm) == n
    for col in raw.columns:
        assert col in harm.columns

    mask_t = select_shi_compartment_mask(harm, ShiCompartmentFilter.T_NK_ONLY)
    mask_caf = select_shi_compartment_mask(harm, ShiCompartmentFilter.STROMAL_CAF_ONLY)
    mask_both = select_shi_compartment_mask(harm, ShiCompartmentFilter.T_AND_CAF)
    mask_all = select_shi_compartment_mask(harm, ShiCompartmentFilter.ALL_CELLS)

    np.testing.assert_array_equal(mask_both, np.logical_or(mask_t, mask_caf))
    assert bool(np.all(mask_all))
    assert not bool(np.any(np.logical_and(mask_t, mask_caf)))


# =============================================================================
# 3. Out-of-Core HDF5 CSR Slicing, R `.rda` / `SingleCellExperiment` & Guard Tests
# =============================================================================


def test_shi_pancancer_out_of_core_csr_slice_and_load_dataset(tmp_path: Path) -> None:
    """Verify out-of-core HDF5 CSR row slicing on a synthetic Shi_PanCancer_Atlas `.h5ad`."""
    raw_dir = tmp_path / "data" / "raw" / SHI_PANCANCER_DATASET_ID
    prep_dir = tmp_path / "data" / "preprocessed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    prep_dir.mkdir(parents=True, exist_ok=True)
    (tmp_path / "pyproject.toml").write_text("[project]\nname = 'mock_repo'\n", encoding="utf-8")

    counts_dense = np.array(
        [
            [10.0, 0.0, 5.0, 0.0],   # 0: CD8 Tex -> keep
            [0.0, 20.0, 0.0, 0.0],   # 1: Macrophage -> skip
            [4.0, 0.0, 12.0, 2.0],   # 2: myCAF -> keep
            [0.0, 0.0, 0.0, 15.0],   # 3: Epithelial -> skip
            [8.0, 2.0, 0.0, 0.0],    # 4: CD4 Treg -> keep
            [1.0, 0.0, 9.0, 0.0],    # 5: iCAF -> keep
        ],
        dtype=np.float32,
    )
    obs_df = pd.DataFrame(
        {
            "cell_type_major": ["CD8 T cell", "Macrophage", "Fibroblast", "Epithelial", "CD4 T cell", "Fibroblast"],
            "cell_type_fine": ["CD8_Tex", "SPP1_TAM", "LRRC15_myCAF", "Tumor_Epi", "FOXP3_Treg", "CXCL12_iCAF"],
            "cancer_type": ["NSCLC", "NSCLC", "PDAC", "PDAC", "MEL", "BRCA"],
            "donor_id": ["D1", "D1", "D2", "D2", "D3", "D4"],
            "cm_ecotype": ["CM03", "CM05", "CM08", "CM01", "CM03", "CM07"],
        },
        index=[f"shi_cell_{i}" for i in range(6)],
    )
    var_df = pd.DataFrame(
        {"gene_name": ["PDCD1", "SPP1", "LRRC15", "EPCAM"]},
        index=["PDCD1", "SPP1", "LRRC15", "EPCAM"],
    )
    umap_coords = np.arange(12, dtype=np.float32).reshape(6, 2)

    raw_h5ad = raw_dir / "PanCancer_igt_s9_fine_counts.h5ad"
    raw_adata = ad.AnnData(
        X=sp.csr_matrix(counts_dense),
        obs=obs_df,
        var=var_df,
        obsm={"X_umap": umap_coords},
    )
    raw_adata.write_h5ad(raw_h5ad)

    out_h5ad = prep_dir / f"{SHI_PANCANCER_DATASET_ID}.h5ad"
    prep_res = preprocess_shi_pancancer_atlas(
        h5ad_path_or_dir=raw_h5ad,
        output_h5ad=Some(out_h5ad),
        compartment_filter=ShiCompartmentFilter.T_AND_CAF,
        normalize_ensembl=False,
        batch_size=2,
    )
    assert isinstance(prep_res, Success)
    assert out_h5ad.is_file()

    load_res = load_dataset(
        "Shi2025",
        base_dir=tmp_path,
        auto_download=False,
        normalize_ensembl=False,
        include_all_layers=True,
    )
    assert isinstance(load_res, Success)
    sliced = load_res.unwrap()
    assert sliced.n_obs == 4
    assert list(sliced.obs_names) == ["shi_cell_0", "shi_cell_2", "shi_cell_4", "shi_cell_5"]
    assert "cm_ecotype" in sliced.obs.columns
    assert set(sliced.obs["lineage_compartment"]) == {"T_CELL", "STROMAL_CAF"}
    assert "counts" in sliced.layers
    assert "log1p_norm" in sliced.layers
    np.testing.assert_allclose(
        sliced.layers["counts"].toarray(),
        counts_dense[[0, 2, 4, 5], :],
        rtol=1e-6,
    )
    np.testing.assert_allclose(sliced.obsm["X_umap"], umap_coords[[0, 2, 4, 5], :], rtol=1e-6)


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript required for R .rda / SCE round-trip test")
def test_r_rda_and_sce_dgtmatrix_streaming_round_trip(tmp_path: Path) -> None:
    """Verify R `.rda` (`Gao_CAF_Atlas` layout) and `SingleCellExperiment` (`dgTMatrix`) streaming."""
    raw_gao = tmp_path / "Gao_CAF_Atlas"
    raw_gao.mkdir(parents=True, exist_ok=True)
    expr_rda = raw_gao / "Fib_Atlas_expression.rda"
    meta_rda = raw_gao / "Fib_Atlas_metadata.rda"
    sce_rds = tmp_path / "mock_zheng_cohort.sce.rds"

    r_code = f"""
    library(Matrix)
    cells <- c("c1", "c2", "c3", "c4")
    genes <- c("FAP", "LRRC15", "PI16")
    mat_dense <- matrix(c(5, 0, 2, 0, 10, 0, 1, 3, 0, 4, 0, 8), nrow = 3, ncol = 4,
                        dimnames = list(genes, cells))
    expression <- as(mat_dense, "dgCMatrix")
    save(expression, file = "{expr_rda.as_posix()}")

    metadata <- data.frame(
        Cluster = c("c04", "c01", "c19", "c10"),
        Tissue = c("Lung", "Pancreas", "Breast", "Colon"),
        Phenotype = c("Tumor", "Tumor", "Normal", "Normal"),
        Dataset = c("S1", "S1", "S2", "S2"),
        UMAP_1 = c("NA", "NA", "NA", "NA"),
        UMAP_2 = c("NA", "NA", "NA", "NA"),
        row.names = cells,
        stringsAsFactors = FALSE
    )
    save(metadata, file = "{meta_rda.as_posix()}")

    dgt <- as(mat_dense, "dgTMatrix")
    env_assays <- new.env(parent = emptyenv())
    assign(".->data", list(counts = dgt), envir = env_assays)
    setClass("MockSimpleListAssays", slots = c(.xData = "environment"))
    setClass("MockDFrame", slots = c(rownames = "character", nrows = "integer", listData = "list"))
    setClass("MockSCE", slots = c(
        assays = "MockSimpleListAssays",
        colData = "MockDFrame",
        int_colData = "MockDFrame"
    ))
    cd <- new("MockDFrame", rownames = cells, nrows = 4L, listData = list(
        meta.cluster = c("CD8.c12.Tex.CXCL13", "CD8.c01.Tn.MAL", "CD4.c20.Treg.TNFRSF9", "CD8.c10.Trm.ZNF683"),
        loc = c("T", "P", "T", "N")
    ))
    icd <- new("MockDFrame", rownames = cells, nrows = 4L, listData = list(reducedDims = list()))
    sce_obj <- new("MockSCE",
        assays = new("MockSimpleListAssays", .xData = env_assays),
        colData = cd,
        int_colData = icd
    )
    saveRDS(sce_obj, file = "{sce_rds.as_posix()}")
    """
    subprocess.run(["Rscript", "-e", r_code], check=True, capture_output=True, text=True)

    df_res = read_r_dataframe(meta_rda)
    assert isinstance(df_res, Success)
    assert list(df_res.unwrap().index) == ["c1", "c2", "c3", "c4"]

    out_gao = tmp_path / "Gao_CAF_Atlas.h5ad"
    gao_res = preprocess_gao_caf_atlas(
        raw_dir=raw_gao,
        output_h5ad=Some(out_gao),
        normalize_ensembl=False,
        batch_size=2,
    )
    assert isinstance(gao_res, Success)
    gao_adata = ad.read_h5ad(out_gao)
    assert gao_adata.shape == (4, 3)
    assert "Dataset" in gao_adata.obs.columns
    assert list(gao_adata.obs["lineage_compartment"]) == ["STROMAL_CAF"] * 4

    stage_h5 = tmp_path / "sce_stage.h5"
    sce_res = stream_r_object_to_staging_h5(sce_rds, stage_h5)
    assert isinstance(sce_res, Success)
    man = sce_res.unwrap()
    assert len(man.obs_df) == 4
    assert "counts" in man.assays
    assert man.assays["counts"].n_cells == 4
    assert list(man.assays["counts"].gene_names) == ["FAP", "LRRC15", "PI16"]
    assert list(man.obs_df.columns) == ["meta.cluster", "loc"]


def test_offline_and_cscs_guards_for_new_atlases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify `auto_download=False` actionable failures and CSCS Eiger download blocks."""
    (tmp_path / "pyproject.toml").write_text("[project]\nname = 'mock_repo'\n", encoding="utf-8")
    empty_raw = tmp_path / "empty_raw"
    empty_raw.mkdir(parents=True, exist_ok=True)

    for choice_caf in (CafAtlasChoice.GAO, CafAtlasChoice.CORDS, CafAtlasChoice.BUECHLER):
        res_caf = dispatch_caf_atlas_load(
            choice_caf,
            raw_dir=Some(empty_raw),
            auto_download=False,
            repo_root=Some(tmp_path),
        )
        assert isinstance(res_caf, Failure)

    for choice_t in (TCellAtlasChoice.ZHENG, TCellAtlasChoice.CHU, TCellAtlasChoice.SHI):
        res_t = dispatch_tcell_atlas_load(
            choice_t,
            raw_dir=Some(empty_raw),
            auto_download=False,
            repo_root=Some(tmp_path),
        )
        assert isinstance(res_t, Failure)

    chu_msg = dispatch_tcell_atlas_load(
        TCellAtlasChoice.CHU,
        raw_dir=Some(empty_raw),
        auto_download=False,
        repo_root=Some(tmp_path),
    ).failure()
    assert "10.1038/s41591-023-02371-y" in chu_msg
    assert "Chu_TCell_Atlas" in chu_msg

    monkeypatch.setenv("SLURM_CLUSTER_NAME", "eiger")
    assert isinstance(download_gao_caf_atlas(raw_dir=Some(empty_raw)), Failure)
    assert isinstance(download_buechler_caf_atlas(raw_dir=Some(empty_raw), repo_root=Some(tmp_path)), Failure)
    assert isinstance(download_zheng_tcell_atlas(raw_dir=Some(empty_raw), repo_root=Some(tmp_path)), Failure)
    assert isinstance(download_chu_tcell_atlas(raw_dir=Some(empty_raw), repo_root=Some(tmp_path)), Failure)
    assert isinstance(download_shi_pancancer_atlas(raw_dir=Some(empty_raw), repo_root=Some(tmp_path)), Failure)
