"""Unit tests for the individual dataset processing pipeline (02_convert_to_sparse_h5.py, 03_preprocess_datasets.py, 04_plot_dataset_embeddings.py)."""

from pathlib import Path
import pytest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad
from returns.result import Success, Failure

import importlib
raw_to_sparse = importlib.import_module("scripts.02_convert_to_sparse_h5")
preprocess_h5 = importlib.import_module("scripts.03_preprocess_datasets")
plot_h5 = importlib.import_module("scripts.04_plot_dataset_embeddings")

convert_df_to_sparse_adata = raw_to_sparse.convert_df_to_sparse_adata
preprocess_single_dataset = preprocess_h5.preprocess_single_dataset
harmonize_dataset_metadata = preprocess_h5.harmonize_dataset_metadata
compute_dataset_embeddings = plot_h5.compute_dataset_embeddings
export_dataset_svgs = plot_h5.export_dataset_svgs


def test_convert_df_to_sparse_adata_individual() -> None:
    """Tests raw DataFrame matrix to sparse CSR AnnData conversion with per-dataset gene normalization."""
    df_raw = pd.DataFrame(
        np.random.randint(10, 500, size=(50, 20)),
        index=[f"Cell_{i}" for i in range(50)],
        columns=[f"ENSG00000{j:06d}" for j in range(20)],
    )

    res = convert_df_to_sparse_adata(df_raw, "AziziSingleCellMapDiverse2018Adata", min_cell_counts=10)
    assert isinstance(res, Success)
    adata = res.unwrap()

    assert isinstance(adata.X, csr_matrix)
    assert adata.n_obs == 50
    assert adata.n_vars == 20
    assert "original.barcode" in adata.obs.columns
    assert (adata.obs["dataset"] == "AziziSingleCellMapDiverse2018Adata").all()


def test_harmonize_dataset_metadata() -> None:
    """Tests metadata column harmonization (dataset, cancer_code, staging, sequencing_tech)."""
    X = csr_matrix(np.random.randint(10, 100, size=(10, 5)))
    obs = pd.DataFrame({"geo_id": ["GEO1"] * 10, "patient": ["P1"] * 10, "tissue": ["Breast"] * 10})
    var = pd.DataFrame(index=[f"G{i}" for i in range(5)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    adata_harm = harmonize_dataset_metadata(adata, "AziziSingleCellMapDiverse2018Adata")
    
    assert adata_harm.obs["cancer_code"].iloc[0] == "BRCA"
    assert adata_harm.obs["staging"].iloc[0] == "primary tumor"
    assert adata_harm.obs["sequencing_tech"].iloc[0] == "10x_UMI"
    assert "geo_id" in adata_harm.obs.columns  # Preserves original obs fields


def test_preprocess_single_dataset_individual(tmp_path: Path) -> None:
    """Tests individual dataset QC filtering, metadata harmonization, and preprocessed .h5ad export."""
    X = csr_matrix(np.random.randint(50, 1000, size=(60, 20), dtype=np.int64))
    obs = pd.DataFrame({"original.barcode": [f"Cell_{i}" for i in range(60)]})
    var = pd.DataFrame({"contig": ["1"] * 18 + ["MT"] * 2}, index=[f"ENSG00000{i:06d}" for i in range(20)])
    
    sparse_h5 = tmp_path / "Azizi_sparse.h5ad"
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(sparse_h5)

    out_h5 = tmp_path / "Azizi_processed.h5ad"
    qc_spec = preprocess_h5.DatasetQCSpec(min_genes=5, max_genes=100, min_counts=50, max_counts=50000, max_mt_content=50.0)

    res = preprocess_single_dataset(sparse_h5, qc_spec, out_h5)
    assert isinstance(res, Success)
    assert out_h5.exists()

    loaded = ad.read_h5ad(out_h5)
    assert loaded.n_obs > 0
    assert loaded.obs["cancer_code"].iloc[0] == "BRCA"
    assert "MT" not in loaded.var["contig"].values


def test_plot_dataset_embeddings_individual(tmp_path: Path) -> None:
    """Tests per-dataset PCA/UMAP embedding computation and SVG scatter plot export."""
    X = csr_matrix(np.random.randint(10, 500, size=(40, 15), dtype=np.int64))
    obs = pd.DataFrame({
        "cell_type": ["T_cell"] * 20 + ["Myeloid"] * 20,
        "cancer_code": ["BRCA"] * 40,
        "sequencing_tech": ["10x_UMI"] * 40,
    })
    var = pd.DataFrame(index=[f"G{i}" for i in range(15)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    res_emb = compute_dataset_embeddings(adata, n_hvg=10, n_pcs=5)
    assert isinstance(res_emb, Success)
    adata_emb = res_emb.unwrap()

    assert "X_pca" in adata_emb.obsm
    assert "X_umap" in adata_emb.obsm

    res_svg = export_dataset_svgs(adata_emb, "Azizi2018", tmp_path)
    assert isinstance(res_svg, Success)
    paths = res_svg.unwrap()

    assert len(paths) == 6
    for p in paths:
        assert p.exists()
        assert p.suffix == ".svg"
