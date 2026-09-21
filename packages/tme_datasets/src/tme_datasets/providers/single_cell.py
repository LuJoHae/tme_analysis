"""Single-cell reference and ICB response dataset providers."""

from __future__ import annotations

import gzip
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from ..models import DatasetSpec, QualityControlSpec
from ..preprocessing.gene_filtering import filter_confounding_genes
from ..preprocessing.metadata import harmonize_obs_metadata
from ..preprocessing.normalization import normalize_total_counts
from ..types import Modality


def load_sade_feldman(
    raw_or_scratch_dir: Path,
    auto_download: bool = True,
) -> Result[ad.AnnData, str]:
    """Load and preprocess Sade-Feldman et al. 2018 (GSE120575, 16,288 cells)."""
    parquet_path = raw_or_scratch_dir / "gse120575_tpm.parquet"
    meta_path = raw_or_scratch_dir / "gse120575_tpm_cell_metadata.parquet"

    # Fallback to data/raw/GSE120575
    gz_tpm = raw_or_scratch_dir / "GSE120575_tpm.txt.gz"
    gz_meta = raw_or_scratch_dir / "GSE120575_meta.txt.gz"

    # Auto-download if missing
    if auto_download and not parquet_path.exists() and not gz_tpm.exists():
        from ..download.fetcher import download_single_file

        url_tpm = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
        url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz"
        download_single_file(url_tpm, gz_tpm)
        download_single_file(url_meta, gz_meta)

    try:
        if parquet_path.exists() and meta_path.exists():
            df_tpm = pl.read_parquet(parquet_path)
            df_meta = pl.read_parquet(meta_path).to_pandas()
            genes = df_tpm["gene"].to_list()
            cell_cols = [c for c in df_tpm.columns if c != "gene"]
            X_mat = df_tpm.select(cell_cols).to_numpy().T.astype(np.float32)

            adata = ad.AnnData(
                X=sp.csr_matrix(X_mat),
                obs=df_meta.set_index(df_meta.columns[0]),
                var=pd.DataFrame(index=genes),
            )
            return Success(adata)

        elif gz_tpm.exists():
            # Robust reading of GEO TPM file with ragged line 2
            with gzip.open(gz_tpm, "rt", encoding="utf-8", errors="replace") as f:
                header = f.readline().rstrip("\r\n").split("\t")
            cell_ids = [c for c in header[1:] if c.strip()]
            usecols = [0] + list(range(1, len(cell_ids) + 1))

            df = pd.read_csv(
                gz_tpm,
                sep="\t",
                skiprows=2,
                header=None,
                usecols=usecols,
                index_col=0,
                engine="c",
            )
            df.columns = cell_ids

            adata = ad.AnnData(
                X=sp.csr_matrix(df.values.T.astype(np.float32)),
                var=pd.DataFrame(index=df.index),
            )
            adata.obs_names = list(cell_ids)

            if gz_meta.exists():
                meta = pd.read_csv(
                    gz_meta,
                    sep="\t",
                    skiprows=19,
                    encoding="latin1",
                )
                if "title" in meta.columns:
                    meta = meta.dropna(subset=["title"])
                    meta = meta.drop_duplicates(subset=["title"])
                    meta = meta.set_index("title")
                    meta_filtered = meta.reindex(adata.obs_names)
                    adata.obs = meta_filtered

            return Success(adata)
        else:
            return Failure(f"Sade-Feldman data files not found in {raw_or_scratch_dir}")
    except Exception as exc:
        return Failure(f"Failed to load Sade-Feldman: {exc}")


def load_jerby_arnon(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Load Jerby-Arnon et al. 2018 (GSE115978, 7,186 cells)."""
    tpm_path = raw_dir / "GSE115978_tpm.csv.gz"
    if not tpm_path.exists():
        tpm_path = raw_dir / "GSE115978_tpm.gz"

    if not tpm_path.exists():
        return Failure(f"GSE115978 file not found in {raw_dir}")

    try:
        df = pd.read_csv(tpm_path, index_col=0)
        adata = ad.AnnData(
            X=sp.csr_matrix(df.values.T.astype(np.float32)),
            obs=pd.DataFrame(index=df.columns),
            var=pd.DataFrame(index=df.index),
        )
        return filter_confounding_genes(adata)
    except Exception as exc:
        return Failure(f"Failed to load Jerby-Arnon: {exc}")


def load_maynard(repo_root: Path) -> Result[ad.AnnData, str]:
    """Load Maynard et al. 2020 NSCLC dataset (3,000 cells)."""
    h5ad_path = repo_root / "jupyter/data/maynard2020_3k.h5ad"
    if not h5ad_path.exists():
        return Failure(f"Maynard dataset not found at {h5ad_path}")

    try:
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load Maynard: {exc}")


def load_ma_liver(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Load Ma et al. 2019 HCC dataset (GSE125449, 5,115 cells)."""
    matrix_path = raw_dir / "GSE125449_set1_matrix.gz"
    genes_path = raw_dir / "GSE125449_set1_genes.gz"
    barcodes_path = raw_dir / "GSE125449_set1_barcodes.gz"

    if not matrix_path.exists():
        return Failure(f"GSE125449 matrix not found in {raw_dir}")

    try:
        import scipy.io as sio
        mat = sio.mmread(matrix_path).T.tocsr()
        genes = pd.read_csv(genes_path, header=None)[0].tolist()
        barcodes = pd.read_csv(barcodes_path, header=None)[0].tolist()

        adata = ad.AnnData(
            X=mat.astype(np.float32),
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=genes),
        )
        return normalize_total_counts(adata)
    except Exception as exc:
        return Failure(f"Failed to load Ma Liver: {exc}")
