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

from ..logging import get_logger
from ..models import DatasetSpec, QualityControlSpec
from ..preprocessing.gene_filtering import filter_confounding_genes
from ..preprocessing.metadata import harmonize_obs_metadata
from ..preprocessing.normalization import normalize_total_counts
from ..types import Modality

logger = get_logger("providers.single_cell")


def load_sade_feldman(
    raw_or_scratch_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load and preprocess Sade-Feldman et al. 2018 (GSE120575, 16,288 cells)."""
    parquet_path = raw_or_scratch_dir / "gse120575_tpm.parquet"
    meta_path = raw_or_scratch_dir / "gse120575_tpm_cell_metadata.parquet"

    # Fallback to data/raw/GSE120575
    gz_tpm = raw_or_scratch_dir / "GSE120575_tpm.txt.gz"
    gz_meta = raw_or_scratch_dir / "GSE120575_meta.txt.gz"

    # Auto-download if missing or force_download
    if (auto_download or force_download) and (force_download or (not parquet_path.exists() and not gz_tpm.exists())):
        from ..download.fetcher import download_single_file

        url_tpm = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
        url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz"
        if force_download:
            if gz_tpm.exists():
                gz_tpm.unlink()
            if gz_meta.exists():
                gz_meta.unlink()
        download_single_file(url_tpm, gz_tpm)
        download_single_file(url_meta, gz_meta)

    try:
        if not force_download and parquet_path.exists() and meta_path.exists():
            logger.info("Loading preprocessed Sade-Feldman parquet files from %s...", raw_or_scratch_dir)
            df_tpm = pl.read_parquet(parquet_path)
            df_meta = pl.read_parquet(meta_path).to_pandas()
            genes = df_tpm["gene"].to_list()
            cell_cols = [c for c in df_tpm.columns if c != "gene"]
            X_mat = df_tpm.select(cell_cols).to_numpy().T.astype(np.float32)

            obs_df = df_meta.set_index(df_meta.columns[0])
            obs_df.index = obs_df.index.astype(str)
            obs_df.index.name = "cell_id"

            var_df = pd.DataFrame(index=pd.Index(genes, dtype=str, name="gene_id"))

            adata = ad.AnnData(
                X=sp.csr_matrix(X_mat),
                obs=obs_df,
                var=var_df,
            )
            return Success(adata)

        elif gz_tpm.exists():
            # Robust reading of GEO TPM file with ragged line 2
            logger.info("Parsing gzipped GEO TPM expression matrix (%s)...", gz_tpm.name)
            with gzip.open(gz_tpm, "rt", encoding="utf-8", errors="replace") as f:
                header = f.readline().rstrip("\r\n").split("\t")
            cell_ids = [c for c in header[1:] if c.strip()]
            usecols = [0] + list(range(1, len(cell_ids) + 1))
            logger.debug("Identified %d cell barcode columns in header", len(cell_ids))

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
            df.index = df.index.astype(str)
            df.index.name = "gene"
            logger.info("Parsed %d genes across %d single cells from %s", len(df), len(cell_ids), gz_tpm.name)

            var_df = pd.DataFrame(index=pd.Index(df.index, dtype=str, name="gene_id"))
            adata = ad.AnnData(
                X=sp.csr_matrix(df.values.T.astype(np.float32)),
                var=var_df,
            )
            adata.obs_names = list(cell_ids)
            adata.obs_names.name = "cell_id"

            if gz_meta.exists():
                logger.info("Parsing metadata annotations from %s...", gz_meta.name)
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
                    meta.index = meta.index.astype(str)
                    meta_filtered = meta.reindex(adata.obs_names)
                    meta_filtered.index.name = "cell_id"
                    adata.obs = meta_filtered

            # Cache to parquet for near-instant future loading
            try:
                logger.info("Caching parsed GSE120575 matrix to parquet in %s...", raw_or_scratch_dir)
                pl_tpm = pl.from_pandas(df.reset_index())
                pl_tpm.write_parquet(parquet_path)
                pl_meta = pl.from_pandas(adata.obs.reset_index())
                pl_meta.write_parquet(meta_path)
                logger.info("Successfully cached GSE120575 parquet files.")
            except Exception as cache_exc:
                logger.debug("Skipped parquet caching: %s", cache_exc)

            return Success(adata)
        else:
            return Failure(f"Sade-Feldman data files not found in {raw_or_scratch_dir}")
    except Exception as exc:
        msg = f"Failed to load Sade-Feldman: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_jerby_arnon(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load Jerby-Arnon et al. 2018 (GSE115978, 7,186 cells)."""
    tpm_path = raw_dir / "GSE115978_tpm.csv.gz"
    if not tpm_path.exists() and (raw_dir / "GSE115978_tpm.gz").exists():
        tpm_path = raw_dir / "GSE115978_tpm.gz"

    meta_path = raw_dir / "GSE115978_cell.annotations.csv.gz"

    # Auto-download from NCBI GEO if missing or force_download
    if (auto_download or force_download) and (force_download or not tpm_path.exists()):
        from ..download.fetcher import download_single_file

        url_tpm = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_tpm.csv.gz"
        url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_cell.annotations.csv.gz"
        if force_download and tpm_path.exists():
            tpm_path.unlink()
        download_single_file(url_tpm, tpm_path)
        download_single_file(url_meta, meta_path)

    if not tpm_path.exists():
        return Failure(f"GSE115978 file not found in {raw_dir}")

    try:
        logger.info("Loading Jerby-Arnon TPM dataset from %s...", tpm_path.name)
        df = pd.read_csv(tpm_path, index_col=0)
        adata = ad.AnnData(
            X=sp.csr_matrix(df.values.T.astype(np.float32)),
            obs=pd.DataFrame(index=df.columns),
            var=pd.DataFrame(index=df.index),
        )
        logger.info("Filtering confounding genes for Jerby-Arnon (%d cells x %d genes)...", adata.n_obs, adata.n_vars)
        return filter_confounding_genes(adata)
    except Exception as exc:
        msg = f"Failed to load Jerby-Arnon: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_maynard(repo_root: Path) -> Result[ad.AnnData, str]:
    """Load Maynard et al. 2020 NSCLC dataset (3,000 cells)."""
    from ..paths import find_dataset_h5ad

    h5ad_maybe = find_dataset_h5ad("Maynard_NSCLC", repo_root=repo_root)
    match h5ad_maybe:
        case Some(h5ad_path):
            try:
                logger.info("Loading Maynard NSCLC H5AD from %s...", h5ad_path)
                adata = ad.read_h5ad(h5ad_path)
                return Success(adata)
            except Exception as exc:
                msg = f"Failed to load Maynard: {exc}"
                logger.error(msg)
                return Failure(msg)
        case _:
            msg = (
                "Maynard NSCLC dataset not found in candidate paths. "
                "Expected file at 'data/manual_download/Maynard_NSCLC.h5ad' (or 'jupyter/data/maynard2020_3k.h5ad'). "
                "Please run 'python scripts/setup_manual_downloads.py' (or 'make setup-manual-downloads') to copy it."
            )
            logger.error(msg)
            return Failure(msg)


def load_ma_liver(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load Ma et al. 2019 HCC dataset (GSE125449, 5,115 cells)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    matrix_candidates = [
        raw_dir / "GSE125449_Set1_matrix.mtx.gz",
        raw_dir / "GSE125449_set1_matrix.gz",
        raw_dir / "GSE125449_set1_matrix.mtx.gz",
    ]
    genes_candidates = [
        raw_dir / "GSE125449_Set1_genes.tsv.gz",
        raw_dir / "GSE125449_set1_genes.gz",
        raw_dir / "GSE125449_set1_genes.tsv.gz",
    ]
    barcodes_candidates = [
        raw_dir / "GSE125449_Set1_barcodes.tsv.gz",
        raw_dir / "GSE125449_set1_barcodes.gz",
        raw_dir / "GSE125449_set1_barcodes.tsv.gz",
    ]

    matrix_path = next((p for p in matrix_candidates if p.exists()), matrix_candidates[0])
    genes_path = next((p for p in genes_candidates if p.exists()), genes_candidates[0])
    barcodes_path = next((p for p in barcodes_candidates if p.exists()), barcodes_candidates[0])

    if (auto_download or force_download) and (force_download or not matrix_path.exists()):
        from ..download.fetcher import download_single_file

        url_mat = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_matrix.mtx.gz"
        url_genes = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_genes.tsv.gz"
        url_barcodes = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_barcodes.tsv.gz"
        if force_download and matrix_path.exists():
            matrix_path.unlink()
        download_single_file(url_mat, matrix_path)
        download_single_file(url_genes, genes_path)
        download_single_file(url_barcodes, barcodes_path)

    if not matrix_path.exists():
        return Failure(f"GSE125449 matrix file not found in {raw_dir}")

    try:
        logger.info("Loading Ma Liver HCC matrix and features from %s...", raw_dir)
        import scipy.io as sio

        mat = sio.mmread(matrix_path).T.tocsr()
        genes = pd.read_csv(genes_path, header=None, sep=r"\s+")[0].tolist()
        barcodes = pd.read_csv(barcodes_path, header=None, sep=r"\s+")[0].tolist()

        adata = ad.AnnData(
            X=mat.astype(np.float32),
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=genes),
        )
        logger.info("Normalizing total counts for Ma Liver dataset (%d cells x %d genes)...", adata.n_obs, adata.n_vars)
        return normalize_total_counts(adata)
    except Exception as exc:
        msg = f"Failed to load Ma Liver: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_yost(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
    n_cells: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load and preprocess Yost et al. 2019 BCC/SCC dataset (GSE123813)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    counts_candidates = [
        raw_dir / "GSE123813_bcc_scRNA_counts.txt.gz",
        raw_dir / "GSE123813_bcc_counts.txt.gz",
        raw_dir / "GSE123813_scc_scRNA_counts.txt.gz",
    ]
    counts_path = next((p for p in counts_candidates if p.exists()), counts_candidates[0])

    if (auto_download or force_download) and (force_download or not counts_path.exists()):
        from ..download.fetcher import download_single_file

        url_counts = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_bcc_scRNA_counts.txt.gz"
        if force_download and counts_path.exists():
            counts_path.unlink()
        download_single_file(url_counts, counts_path)

    if not counts_path.exists():
        return Failure(f"GSE123813 counts file not found in {raw_dir}")

    try:
        logger.info("Loading Yost et al. BCC counts from %s...", counts_path.name)
        genes: list[str] = []
        matrix_rows: list[np.ndarray] = []
        with gzip.open(counts_path, "rt") as f:
            header = f.readline().strip().split("\t")
            cell_names = header[1 : n_cells + 1] if n_cells else header[1:]
            num_cols = len(cell_names)
            for line in f:
                parts = line.strip().split("\t")
                genes.append(parts[0])
                arr = np.fromiter((float(x) for x in parts[1 : num_cols + 1]), dtype=np.float32, count=num_cols)
                matrix_rows.append(arr)

        raw_mat = np.vstack(matrix_rows)  # (n_genes, n_cells)
        sparse_mat = sp.csr_matrix(raw_mat.T, dtype=np.float32)

        # Parse patient IDs and response mapping
        response_map = {
            "su001": 1, "su002": 1, "su003": 1, "su004": 1, "su009": 1, "su011": 1, "su012": 1,
            "su005": 0, "su006": 0, "su007": 0, "su008": 0, "su010": 0, "su013": 0, "su014": 0,
        }
        obs_df = pd.DataFrame(index=cell_names)

        def _extract_patient(barcode: str) -> str:
            parts = barcode.replace("_", ".").split(".")
            for p in parts:
                if p.lower().startswith("su") and len(p) >= 4 and p[2:].isdigit():
                    return p.lower()
            return "unknown"

        obs_df["patient"] = [_extract_patient(b) for b in cell_names]
        obs_df["response_binary"] = obs_df["patient"].map(response_map)

        adata = ad.AnnData(
            X=sparse_mat,
            obs=obs_df,
            var=pd.DataFrame(index=genes),
        )
        logger.info("Successfully loaded Yost BCC dataset: %d cells x %d genes", adata.n_obs, adata.n_vars)
        return normalize_total_counts(adata)
    except Exception as exc:
        msg = f"Failed to load Yost GSE123813: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_gse179994(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load and preprocess Tietscher et al. Pan-Cancer T-cell Atlas (GSE179994)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    rds_gz_path = raw_dir / "GSE179994_all.Tcell.rawCounts.rds.gz"
    rds_path = raw_dir / "GSE179994_all.Tcell.rawCounts.rds"
    meta_path = raw_dir / "GSE179994_Tcell.metadata.tsv.gz"

    if (auto_download or force_download) and (force_download or (not rds_gz_path.exists() and not rds_path.exists())):
        from ..download.fetcher import download_single_file

        url_rds = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz"
        url_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_Tcell.metadata.tsv.gz"
        if force_download and rds_gz_path.exists():
            rds_gz_path.unlink()
        download_single_file(url_rds, rds_gz_path)
        download_single_file(url_meta, meta_path)

    # Decompress RDS if needed for pyreadr
    if not rds_path.exists() and rds_gz_path.exists():
        logger.info("Decompressing %s...", rds_gz_path.name)
        with gzip.open(rds_gz_path, "rb") as f_in, open(rds_path, "wb") as f_out:
            import shutil

            shutil.copyfileobj(f_in, f_out)

    if not rds_path.exists():
        return Failure(f"GSE179994 RDS count file not found in {raw_dir}")

    try:
        import pyreadr

        logger.info("Reading GSE179994 RDS raw count matrix via pyreadr from %s...", rds_path.name)
        rds_dict = pyreadr.read_r(str(rds_path))
        key = next(iter(rds_dict.keys()))
        df_counts = rds_dict[key]

        if df_counts.shape[0] > df_counts.shape[1]:  # genes x cells
            cell_names = list(df_counts.columns)
            gene_names = list(df_counts.index)
            mat_csr = sp.csr_matrix(df_counts.values.T.astype(np.float32))
        else:  # cells x genes
            cell_names = list(df_counts.index)
            gene_names = list(df_counts.columns)
            mat_csr = sp.csr_matrix(df_counts.values.astype(np.float32))

        obs_df = pd.DataFrame(index=cell_names)

        if meta_path.exists():
            df_meta = pd.read_csv(meta_path, sep="\t", index_col=0)
            obs_df = obs_df.join(df_meta, how="left")
            if "sample" in obs_df.columns:
                s = obs_df["sample"].astype(str)
                is_post = s.str.contains(r"\.post|\.tr", regex=True, case=False)
                is_pre = s.str.contains(r"\.pre|\.ut", regex=True, case=False)
                obs_df["response_treatment"] = np.where(is_post, "On-treatment", np.where(is_pre, "Pre-treatment", "Unknown"))

        adata = ad.AnnData(
            X=mat_csr,
            obs=obs_df,
            var=pd.DataFrame(index=gene_names),
        )
        logger.info("Successfully loaded GSE179994: %d cells x %d genes", adata.n_obs, adata.n_vars)
        return normalize_total_counts(adata)
    except Exception as exc:
        msg = f"Failed to load GSE179994: {exc}"
        logger.error(msg)
        return Failure(msg)

