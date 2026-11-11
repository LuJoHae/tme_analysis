"""Unit, property-based (Hypothesis), and R XDR round-trip tests for Coulton_TAM_Atlas."""

from __future__ import annotations

import hashlib
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
import shutil
import subprocess
import threading
from typing import Any

import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import pytest
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
import scipy.sparse as sp

import tme_datasets
from tme_datasets.download.fetcher import download_single_file
from tme_datasets.models import ChecksumSpec
from tme_datasets.preprocessing.gene_normalization import _project_sparse_chunk_to_ensembl
from tme_datasets.preprocessing.normalization import standardize_dual_layers
from tme_datasets.providers.tam_atlas import (
    COULTON_TAM_CHECKSUM,
    COULTON_TAM_CLUSTER_LABELS,
    COULTON_TAM_DATASET_ID,
    _build_secondary_assay_projection,
    harmonize_tam_atlas_obs,
    load_coulton_tam_atlas,
    preprocess_coulton_tam_atlas,
)
from tme_datasets.query import load_dataset
from tme_datasets.registry import (
    get_dataset_spec,
    list_preprocessed_datasets,
    list_registered_datasets,
    query_preprocessed_datasets,
    resolve_dataset_id,
)
from tme_datasets.types import Modality


# =============================================================================
# 1. Registry & Alias Resolution Tests
# =============================================================================


def test_coulton_tam_atlas_registered_and_aliases() -> None:
    """Verify Coulton_TAM_Atlas registration metadata and all aliases."""
    spec_maybe = get_dataset_spec("Coulton_TAM_Atlas")
    assert isinstance(spec_maybe, Some)
    spec = spec_maybe.unwrap()
    assert spec.id == COULTON_TAM_DATASET_ID
    assert spec.modality == Modality.SINGLE_CELL
    assert spec.n_samples_or_cells == Some(363315)
    assert spec.publication_pmid == Some("38969631")
    assert spec.publication_doi == Some("10.1038/s41467-024-49885-8")
    assert spec.checksum == Some(COULTON_TAM_CHECKSUM)
    assert spec.qc_spec is Nothing
    assert spec.is_dysfunctional is False

    assert "Coulton_TAM_Atlas" in list_registered_datasets().ids()
    assert len(COULTON_TAM_CLUSTER_LABELS) == 24

    for alias in (
        "TAM_Atlas",
        "Coulton2024",
        "Coulton_2024",
        "Coulton",
        "Zenodo_11222158",
        "11222158",
        "mac.atlas.zenodo.200524",
        "mac.atlas.zenodo.200524.rds",
    ):
        assert resolve_dataset_id(alias) == "Coulton_TAM_Atlas"
        assert get_dataset_spec(alias) == Some(spec)

    # Top-level API exports
    assert hasattr(tme_datasets, "download_coulton_tam_atlas")
    assert hasattr(tme_datasets, "preprocess_coulton_tam_atlas")
    assert hasattr(tme_datasets, "load_coulton_tam_atlas")
    assert hasattr(tme_datasets, "COULTON_TAM_CLUSTER_LABELS")


# =============================================================================
# 2. Hypothesis Property-Based Tests
# =============================================================================


@settings(max_examples=40, deadline=None)
@given(
    n_cells=st.integers(min_value=1, max_value=12),
    n_raw_genes=st.integers(min_value=2, max_value=10),
    seed=st.integers(min_value=0, max_value=10_000),
)
def test_property_log1p_sum_ensembl_projection(
    n_cells: int,
    n_raw_genes: int,
    seed: int,
) -> None:
    """Property: `log1p_sum` on `log1p(C)` equals `log1p(C_kept @ M)` for all non-negative matrices."""
    rng = np.random.default_rng(seed)
    linear_counts = rng.uniform(0.0, 25.0, size=(n_cells, n_raw_genes)).astype(np.float64)
    zero_mask = rng.random(size=(n_cells, n_raw_genes)) < 0.5
    linear_counts[zero_mask] = 0.0
    log_input = sp.csr_matrix(np.log1p(linear_counts).astype(np.float32))

    keep_mask = rng.random(size=n_raw_genes) >= 0.2
    if not bool(keep_mask.any()):
        keep_mask[0] = True

    n_kept = int(keep_mask.sum())
    n_targets = max(1, n_kept // 2)
    target_cols = rng.integers(0, n_targets, size=n_kept, dtype=np.int32)
    M = sp.csr_matrix(
        (np.ones(n_kept, dtype=np.float32), (np.arange(n_kept, dtype=np.int32), target_cols)),
        shape=(n_kept, n_targets),
        dtype=np.float32,
    )

    projected = _project_sparse_chunk_to_ensembl(log_input, keep_mask, M, "log1p_sum")
    expected_linear = np.expm1(log_input[:, keep_mask].toarray().astype(np.float64)) @ M.toarray().astype(
        np.float64
    )
    expected_log = np.log1p(expected_linear).astype(np.float32)

    assert projected.shape == (n_cells, n_targets)
    assert projected.dtype == np.float32
    np.testing.assert_allclose(projected.toarray(), expected_log, rtol=1e-5, atol=1e-5)


@settings(max_examples=35, deadline=None)
@given(
    cluster_ids=st.lists(st.integers(min_value=0, max_value=23), min_size=1, max_size=25),
    cancers=st.lists(
        st.sampled_from(["thyroid", "melanoma", "breast", "liver", "ESO", "LUAD", "CRC", "PDAC"]),
        min_size=1,
        max_size=25,
    ),
    prim_mets=st.lists(
        st.sampled_from(["pri", "primary/laparoscopy", "met", "metastatic", "norm", "NA"]),
        min_size=1,
        max_size=25,
    ),
)
def test_property_harmonize_tam_atlas_obs_invariants(
    cluster_ids: list[int],
    cancers: list[str],
    prim_mets: list[str],
) -> None:
    """Property: Polars metadata harmonization preserves rows and maps all 24 TAM states & cancer codes."""
    n = len(cluster_ids)
    cancer_col = [cancers[i % len(cancers)] for i in range(n)]
    pmet_col = [prim_mets[i % len(prim_mets)] for i in range(n)]
    barcodes = [f"cell_{i}" for i in range(n)]

    raw_df = pd.DataFrame(
        {
            "study": [f"Study_{i % 3}" for i in range(n)],
            "patient": [f"Pat_{i}" for i in range(n)],
            "sample": ["NA" if i % 2 == 0 else f"Samp_{i}" for i in range(n)],
            "cancer": cancer_col,
            "tissue": ["NA" if i % 2 == 0 else "Lung" for i in range(n)],
            "primary_met2": pmet_col,
            "integrated_snn_res.0.525": [str(k) for k in cluster_ids],
        },
        index=barcodes,
    )

    out = harmonize_tam_atlas_obs(raw_df)
    assert len(out) == n
    assert list(out.index) == barcodes

    for idx, cid in enumerate(cluster_ids):
        expected_short = COULTON_TAM_CLUSTER_LABELS[str(cid)]
        assert out["short.label"].iloc[idx] == expected_short
        assert out["cell_state"].iloc[idx] == expected_short
        if expected_short == "19_ClassMono":
            assert out["cell_type"].iloc[idx] == "Monocyte"
        elif expected_short == "20_TDoub":
            assert out["cell_type"].iloc[idx] == "T/Myeloid Doublet"
        else:
            assert out["cell_type"].iloc[idx] == "Macrophage"

        if idx % 2 == 0:
            # sample == "NA" falls back to patient (`mac.util.R` rule)
            assert out["sample"].iloc[idx] == f"Pat_{idx}"
            assert out["tissue"].iloc[idx] == "Tumor"

    assert set(out["cancer_code"]).issubset({"THCA", "MEL", "BC", "LIHC", "ESCA", "LUAD", "CRC", "PDAC"})
    assert set(out["primary_met_harmonized"]).issubset({"primary", "metastasis", "normal", "unknown"})


@settings(max_examples=30, deadline=None)
@given(
    target_size=st.integers(min_value=3, max_value=20),
    sec_indices=st.lists(st.integers(min_value=0, max_value=24), min_size=0, max_size=15),
)
def test_property_secondary_assay_projection_alignment(
    target_size: int,
    sec_indices: list[int],
) -> None:
    """Property: Secondary assay (`SCT`) sparse projection aligns columns 1-to-1 with target `.var_names`."""
    target_ids = [f"ENSG{i:011d}" for i in range(target_size)]
    raw_to_target = {f"SYM_{i}": f"ENSG{i:011d}" for i in range(target_size)}
    sec_genes = [f"SYM_{idx}" for idx in sec_indices]

    spec = _build_secondary_assay_projection(sec_genes, target_ids, raw_to_target)
    expected_covered = {idx for idx in sec_indices if idx < target_size}
    actual_covered = set(np.flatnonzero(spec.covered_var_mask).tolist())
    assert actual_covered == expected_covered
    assert spec.projection_matrix.shape == (int(spec.keep_mask.sum()), target_size)


# =============================================================================
# 3. Pre-Normalized Dual-Layer Safeguard & Resumable HTTP 206 Download Tests
# =============================================================================


def test_standardize_dual_layers_preserves_log_normalized_without_counts() -> None:
    """Verify `standardize_dual_layers` does not fabricate `counts` or double-log1p pre-normalized data."""
    X_log = sp.csr_matrix([[0.0, 1.3862944], [2.1972246, 0.0]], dtype=np.float32)
    adata = ad.AnnData(
        X=X_log.copy(),
        uns={"is_raw_counts": False, "expression_type": "log_normalized"},
    )
    res = standardize_dual_layers(adata)
    assert isinstance(res, Success)
    out = res.unwrap()
    assert "counts" not in out.layers
    assert "log1p_norm" in out.layers
    np.testing.assert_allclose(out.X.toarray(), X_log.toarray(), rtol=1e-6)
    np.testing.assert_allclose(out.layers["log1p_norm"].toarray(), X_log.toarray(), rtol=1e-6)


def test_download_single_file_http_206_resume_and_md5(tmp_path: Path) -> None:
    """Verify HTTP `206 Partial Content` resume from `.part` and MD5 verification."""
    payload = b"Coulton_TAM_Atlas_Zenodo_Stream_" * 4096  # 131,072 bytes
    expected_md5 = hashlib.md5(payload).hexdigest()

    class _RangeHandler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802
            range_hdr = self.headers.get("Range")
            if range_hdr and range_hdr.startswith("bytes="):
                offset = int(range_hdr.split("=", 1)[1].split("-", 1)[0])
                body = payload[offset:]
                self.send_response(206)
                self.send_header("Content-Type", "application/octet-stream")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Content-Range", f"bytes {offset}-{len(payload) - 1}/{len(payload)}")
                self.end_headers()
                self.wfile.write(body)
            else:
                self.send_response(200)
                self.send_header("Content-Type", "application/octet-stream")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

    server = HTTPServer(("127.0.0.1", 0), _RangeHandler)
    port = server.server_port
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        dest = tmp_path / "sample.rds"
        part = tmp_path / "sample.rds.part"
        # Pre-populate first 50,000 bytes in `.part`
        part.write_bytes(payload[:50_000])

        good_spec = ChecksumSpec(algorithm="md5", expected_hash=expected_md5, size_bytes=Some(len(payload)))
        res = download_single_file(
            f"http://127.0.0.1:{port}/sample.rds",
            dest,
            resume=True,
            checksum=Some(good_spec),
        )
        assert isinstance(res, Success)
        assert dest.read_bytes() == payload
        assert not part.exists()

        # Verify wrong MD5 fails cleanly
        bad_dest = tmp_path / "bad.rds"
        bad_spec = ChecksumSpec(algorithm="md5", expected_hash="0" * 32, size_bytes=Some(len(payload)))
        bad_res = download_single_file(
            f"http://127.0.0.1:{port}/sample.rds",
            bad_dest,
            resume=False,
            max_retries=1,
            checksum=Some(bad_spec),
        )
        assert isinstance(bad_res, Failure)
        assert not bad_dest.exists()
    finally:
        server.shutdown()
        server.server_close()


# =============================================================================
# 4. End-to-End R 4.4 Seurat S4 `.rds` -> `.h5ad` -> `tme_datasets` API Test
# =============================================================================


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript required for S4 RDS round-trip test")
@pytest.mark.parametrize("integrated_is_sparse", [True, False])
def test_end_to_end_r_seurat_rds_to_h5ad_and_api(
    tmp_path: Path,
    integrated_is_sparse: bool,
) -> None:
    """Create a multi-assay Seurat S4 `.rds` in R, stream to `.h5ad`, and query via `tme_datasets`."""
    raw_dir = tmp_path / "data" / "raw" / "Coulton_TAM_Atlas"
    prep_dir = tmp_path / "data" / "preprocessed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    prep_dir.mkdir(parents=True, exist_ok=True)
    (tmp_path / "pyproject.toml").write_text("[project]\nname = 'mock_repo'\n", encoding="utf-8")

    # Symlink local Ensembl cache if available so pyensembl never hits network
    real_ensembl = Path("/Users/halu/Code/tme_analysis/data/ensembl")
    use_ensembl = real_ensembl.is_dir()
    if use_ensembl:
        (tmp_path / "data" / "ensembl").symlink_to(real_ensembl)

    rds_path = raw_dir / "mac.atlas.zenodo.200524.rds"
    int_expr = (
        'as(m_int_dense, "dgCMatrix")'
        if integrated_is_sparse
        else "m_int_dense"
    )

    r_script = f"""
    library(Matrix)
    library(compiler)

    setClass("MockSCTModel", slots = c(model = "ANY", cell.attributes = "data.frame", fn = "ANY"))
    setClass("MockAssay", slots = c(
        counts = "dgCMatrix",
        data = "ANY",
        scale.data = "matrix",
        var.features = "character",
        SCTModel.list = "list"
    ))
    setClass("MockDimReduc", slots = c(cell.embeddings = "matrix", key = "character"))
    setClass("MockSeurat", slots = c(
        assays = "list",
        meta.data = "data.frame",
        active.assay = "character",
        active.ident = "factor",
        graphs = "list",
        commands = "list",
        reductions = "list",
        project.name = "character"
    ))

    genes_rna <- c("CD68", "C1QA", "SPP1", "APOE", "CXCL9", "FCN1", "UNMAPPED_SYNTH_GENE_99")
    cells <- c("Cell_1", "Cell_2", "Cell_3", "Cell_4", "Cell_5")

    # 7 genes x 5 cells log-normalized RNA@data (counts empty!)
    rna_vals <- matrix(c(
        1.5, 0.0, 2.2, 0.0, 0.5, 0.0, 1.1,
        0.0, 1.8, 0.0, 2.4, 0.0, 0.0, 0.0,
        1.2, 0.9, 0.0, 0.0, 3.1, 0.0, 0.0,
        0.0, 0.0, 2.9, 1.7, 0.0, 0.0, 0.4,
        0.2, 0.0, 0.0, 0.0, 0.0, 2.8, 0.0
    ), nrow = 7, ncol = 5, dimnames = list(genes_rna, cells))
    m_rna <- as(rna_vals, "dgCMatrix")
    m_empty <- new("dgCMatrix", Dim = c(7L, 5L), p = rep(0L, 6L), Dimnames = list(genes_rna, cells))

    # SCT@data on 4-gene subset + complex SCTModel.list
    genes_sct <- c("CD68", "SPP1", "APOE", "CXCL9")
    sct_vals <- rna_vals[genes_sct, ] * 0.9
    m_sct <- as(sct_vals, "dgCMatrix")
    env1 <- new.env()
    assign("alpha", 1:50, envir = env1)
    cmp_fn <- cmpfun(function(u) u * 2 + 1)
    sct_mod <- new(
        "MockSCTModel",
        model = as.formula("y ~ log_umi", env = env1),
        cell.attributes = data.frame(umi = 1:5),
        fn = list(cmp_fn, getClass("MockAssay"))
    )

    # integrated@data on 3 anchor genes x 5 cells
    genes_int <- c("SPP1", "APOE", "CXCL9")
    m_int_dense <- matrix(
        seq(from = -0.5, by = 0.25, length.out = 15),
        nrow = 3, ncol = 5,
        dimnames = list(genes_int, cells)
    )

    a_rna <- new("MockAssay", counts = m_empty, data = m_rna, scale.data = matrix(0, 0, 0),
                 var.features = c("CD68", "SPP1", "CXCL9"), SCTModel.list = list())
    a_sct <- new("MockAssay", counts = m_sct, data = m_sct, scale.data = matrix(rnorm(100), 10, 10),
                 var.features = c("SPP1", "APOE"), SCTModel.list = list(s1 = sct_mod))
    a_int <- new("MockAssay", counts = new("dgCMatrix", Dim = c(3L, 5L), p = rep(0L, 6L)),
                 data = {int_expr}, scale.data = matrix(rnorm(100), 10, 10),
                 var.features = genes_int, SCTModel.list = list())

    emb_pca <- matrix(as.numeric(1:20), nrow = 5, ncol = 4,
                      dimnames = list(cells, c("PC_1", "PC_2", "PC_3", "PC_4")))
    emb_umap <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, -1.1, -1.2, -1.3, -1.4, -1.5),
                       nrow = 5, ncol = 2, dimnames = list(cells, c("UMAP_1", "UMAP_2")))

    meta <- data.frame(
        study = factor(c("Azizi2018", "Pelka2021", "Leader2021", "Biermann2022", "Pu2021")),
        patient = c("P1", "P2", "P3", "P4", "P5"),
        sample = c("NA", "S2", "S3", "NA", "S5"),
        cancer = c("breast", "CRC", "LUAD", "melanoma", "thyroid"),
        tissue = c("NA", "Colon", "Lung", "Brain", "Thyroid"),
        primary_met2 = c("pri", "primary", "pri", "met", "norm"),
        pre.or.post.treatment = c("Pre", "Pre", "Post", "Pre", "NA"),
        nCount_RNA = c(1200, 2300, 3100, 1800, 950),
        nFeature_RNA = c(450, 810, 920, 640, 390),
        integrated_snn_res.0.525 = factor(c("0", "3", "4", "9", "19")),
        row.names = cells,
        check.names = FALSE
    )

    seu <- new(
        "MockSeurat",
        assays = list(RNA = a_rna, SCT = a_sct, integrated = a_int),
        meta.data = meta,
        active.assay = "integrated",
        active.ident = factor(c("0", "3", "4", "9", "19")),
        graphs = list(integrated_snn = m_rna),
        commands = list(RunUMAP = list(fn = cmp_fn)),
        reductions = list(pca = new("MockDimReduc", cell.embeddings = emb_pca, key = "PC_"),
                          umap = new("MockDimReduc", cell.embeddings = emb_umap, key = "UMAP_")),
        project.name = "Coulton_TAM_Mock"
    )
    saveRDS(seu, file = "{rds_path}")
    """
    subprocess.run(["Rscript", "-e", r_script], check=True)
    assert rds_path.is_file()

    target_h5ad = prep_dir / "Coulton_TAM_Atlas.h5ad"
    prep_res = preprocess_coulton_tam_atlas(
        rds_path_or_dir=rds_path,
        output_h5ad=Some(target_h5ad),
        normalize_ensembl=use_ensembl,
        batch_size=2,  # Force 3 streaming batches across 5 cells
    )
    assert isinstance(prep_res, Success), f"Preprocessing failed: {prep_res}"
    assert target_h5ad.is_file()

    # Verify discovery through `list_preprocessed_datasets` and `query_preprocessed_datasets`
    prep_list = list_preprocessed_datasets(modality="single_cell", repo_root=tmp_path)
    assert "Coulton_TAM_Atlas" in prep_list.ids()

    queried = query_preprocessed_datasets(["TAM_Atlas"], repo_root=tmp_path)
    assert len(queried) == 1
    assert queried[0].id == "Coulton_TAM_Atlas"

    # 1. Default `load_dataset("TAM_Atlas")` (selective low-RAM read)
    res_default = load_dataset("TAM_Atlas", base_dir=tmp_path, auto_download=False)
    assert isinstance(res_default, Success)
    adata_def = res_default.unwrap()
    assert adata_def.n_obs == 5
    assert adata_def.n_vars == (6 if use_ensembl else 7)
    if use_ensembl:
        assert all(str(v).startswith("ENSG") for v in adata_def.var_names)
    assert "counts" not in adata_def.layers
    assert adata_def.uns["is_raw_counts"] is False
    assert adata_def.uns["expression_type"] == "log_normalized"
    assert "X_pca" in adata_def.obsm and adata_def.obsm["X_pca"].shape == (5, 4)
    assert "X_umap" in adata_def.obsm and adata_def.obsm["X_umap"].shape == (5, 2)
    assert list(adata_def.obs["short.label"]) == [
        "0_AlvMac",
        "3_ICIMac1",
        "4_ICIMac2",
        "9_AngioMac",
        "19_ClassMono",
    ]
    assert list(adata_def.obs["cell_type"]) == [
        "Macrophage",
        "Macrophage",
        "Macrophage",
        "Macrophage",
        "Monocyte",
    ]
    assert list(adata_def.obs["cancer_code"]) == ["BC", "CRC", "LUAD", "MEL", "THCA"]

    # 2. Full multi-layer read (`include_all_layers=True`)
    res_full = load_dataset(
        "Coulton_TAM_Atlas",
        base_dir=tmp_path,
        auto_download=False,
        include_all_layers=True,
    )
    assert isinstance(res_full, Success)
    adata_full = res_full.unwrap()
    assert "log1p_norm" in adata_full.layers
    assert "sct_data" in adata_full.layers
    assert "X_rpca_integrated" in adata_full.obsm
    assert adata_full.obsm["X_rpca_integrated"].shape == (5, 3)
    assert list(adata_full.uns["rpca_integrated_var_names"]) == ["SPP1", "APOE", "CXCL9"]
    np.testing.assert_allclose(
        adata_full.X.toarray(),
        adata_full.layers["log1p_norm"].toarray(),
        rtol=1e-6,
    )
    # Verify exact UMAP column-major -> C-order reconstruction
    np.testing.assert_allclose(
        adata_full.obsm["X_umap"][:, 0],
        np.array([0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float32),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        adata_full.obsm["X_umap"][:, 1],
        np.array([-1.1, -1.2, -1.3, -1.4, -1.5], dtype=np.float32),
        rtol=1e-6,
    )

    # 3. Direct assay selection (`assay="SCT"` and `assay="integrated"`)
    res_sct = load_coulton_tam_atlas(
        output_h5ad=Some(target_h5ad),
        assay="SCT",
        repo_root=Some(tmp_path),
    )
    assert isinstance(res_sct, Success)
    adata_sct = res_sct.unwrap()
    assert adata_sct.shape == adata_full.shape
    np.testing.assert_allclose(
        adata_sct.X.toarray(),
        adata_full.layers["sct_data"].toarray(),
        rtol=1e-6,
    )

    res_int = load_coulton_tam_atlas(
        output_h5ad=Some(target_h5ad),
        assay="integrated",
        repo_root=Some(tmp_path),
    )
    assert isinstance(res_int, Success)
    adata_int = res_int.unwrap()
    assert adata_int.shape == (5, 3)
    assert list(adata_int.var_names) == ["SPP1", "APOE", "CXCL9"]
    np.testing.assert_allclose(
        np.asarray(adata_int.X),
        adata_full.obsm["X_rpca_integrated"],
        rtol=1e-6,
    )
