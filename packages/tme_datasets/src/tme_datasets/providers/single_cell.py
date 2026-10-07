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
from returns.pipeline import is_successful
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
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load and preprocess Sade-Feldman et al. 2018 (GSE120575, 16,288 cells, Smart-seq2)."""
    parquet_path = raw_or_scratch_dir / "gse120575_tpm.parquet"
    meta_path = raw_or_scratch_dir / "gse120575_tpm_cell_metadata.parquet"

    tpm_candidates = [
        raw_or_scratch_dir / "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
        raw_or_scratch_dir / "GSE120575_tpm.txt.gz",
        raw_or_scratch_dir / "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",
    ]
    gz_tpm = next((p for p in tpm_candidates if p.exists()), tpm_candidates[0])

    meta_candidates = [
        raw_or_scratch_dir / "GSE120575_patient_ID_single_cells.txt.gz",
        raw_or_scratch_dir / "GSE120575_meta.txt.gz",
        raw_or_scratch_dir / "GSE120575_patient_ID_single_cells.txt",
    ]
    gz_meta = next((p for p in meta_candidates if p.exists()), meta_candidates[0])

    # Auto-download from NCBI GEO if missing or force_download
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
        adata: ad.AnnData | None = None
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
            var_df["gene_name"] = genes
            adata = ad.AnnData(X=sp.csr_matrix(X_mat), obs=obs_df, var=var_df)

        elif gz_tpm.exists():
            logger.info("Parsing gzipped GEO TPM expression matrix (%s)...", gz_tpm.name)
            open_fn = gzip.open if str(gz_tpm).endswith(".gz") else open
            with open_fn(gz_tpm, "rt", encoding="utf-8", errors="replace") as f:
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
            df.index = df.index.astype(str)
            df.index.name = "gene"

            var_df = pd.DataFrame(index=pd.Index(df.index, dtype=str, name="gene_id"))
            var_df["gene_name"] = df.index.tolist()
            adata = ad.AnnData(
                X=sp.csr_matrix(df.values.T.astype(np.float32)),
                var=var_df,
            )
            adata.obs_names = list(cell_ids)
            adata.obs_names.name = "cell_id"

            if gz_meta.exists():
                logger.info("Parsing metadata annotations from %s...", gz_meta.name)
                open_m_fn = gzip.open if str(gz_meta).endswith(".gz") else open
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

        if adata is None:
            return Failure(f"Sade-Feldman data files not found in {raw_or_scratch_dir}")

        # Extract patient ID and treatment timepoint from metadata columns
        obs = adata.obs.copy()
        patient_col_candidate = [c for c in obs.columns if "patinet" in c.lower() or "patient" in c.lower()]
        if patient_col_candidate:
            raw_pats = obs[patient_col_candidate[0]].astype(str)
            # Typically formatted as 'Pre_P1' or 'Post_P1'
            pats = []
            treats = []
            for val in raw_pats:
                parts = str(val).split("_", 1)
                if len(parts) == 2 and parts[0].lower() in {"pre", "post"}:
                    treats.append("pre-treatment" if parts[0].lower() == "pre" else "post-treatment")
                    pats.append(parts[1])
                else:
                    treats.append("pre-treatment")
                    pats.append(str(val))
            obs["patient_id"] = pats
            obs["sample_id"] = raw_pats
            obs["treatment_status"] = treats

        resp_col_candidate = [c for c in obs.columns if "response" in c.lower()]
        if resp_col_candidate:
            obs["clinical_response_raw"] = obs[resp_col_candidate[0]].astype(str)

        adata.obs = obs

        from .tier0_single_cell import _standardize_tier0_obs, _apply_subset_and_subsample

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE120575",
            indication="Melanoma",
            technology="Smart-seq2",
            cell_selection="FACS-sorted (CD45+)",
            patient_col="patient_id" if "patient_id" in adata.obs.columns else None,
            sample_col="sample_id" if "sample_id" in adata.obs.columns else None,
            response_col="clinical_response_raw" if "clinical_response_raw" in adata.obs.columns else None,
            treatment_col="treatment_status" if "treatment_status" in adata.obs.columns else None,
        )

        # Set layers: linear TPM in X and layers["tpm"], log1p TPM in layers["log1p_norm"]
        if sp.issparse(adata.X):
            tpm_mat = adata.X.tocsr().astype(np.float32)
            log1p_mat = tpm_mat.copy()
            log1p_mat.data = np.log1p(log1p_mat.data)
        else:
            tpm_mat = sp.csr_matrix(adata.X, dtype=np.float32)
            log1p_mat = sp.csr_matrix(np.log1p(np.asarray(adata.X, dtype=np.float32)))

        adata.X = tpm_mat
        adata.layers["tpm"] = tpm_mat
        adata.layers["log1p_norm"] = log1p_mat
        adata.uns["expression_type"] = "tpm"
        adata.uns["is_smartseq2"] = True

        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        msg = f"Failed to load Sade-Feldman: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_jerby_arnon(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Jerby-Arnon et al. 2018 (GSE115978, 7,186 cells, Smart-seq2)."""
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
        obs_df = pd.DataFrame(index=df.columns)
        if meta_path.exists():
            try:
                meta = pd.read_csv(meta_path, index_col=0)
                obs_df = obs_df.join(meta, how="left")
            except Exception as exc:
                logger.warning("Could not join metadata for Jerby-Arnon: %s", exc)

        var_df = pd.DataFrame(index=df.index)
        var_df["gene_name"] = df.index.tolist()

        tpm_mat = sp.csr_matrix(df.values.T.astype(np.float32))
        log1p_mat = tpm_mat.copy()
        log1p_mat.data = np.log1p(log1p_mat.data)

        adata = ad.AnnData(X=tpm_mat, obs=obs_df, var=var_df)
        adata.layers["tpm"] = tpm_mat
        adata.layers["log1p_norm"] = log1p_mat
        adata.uns["expression_type"] = "tpm"
        adata.uns["is_smartseq2"] = True

        from .tier0_single_cell import _standardize_tier0_obs, _apply_subset_and_subsample

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE115978",
            indication="Melanoma",
            technology="Smart-seq2",
            cell_selection="FACS-sorted (CD45+)",
            patient_col="patient" if "patient" in adata.obs.columns else "tumor",
            cell_type_col="cell.types" if "cell.types" in adata.obs.columns else None,
        )

        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        msg = f"Failed to load Jerby-Arnon: {exc}"
        logger.error(msg)
        return Failure(msg)


def download_and_build_maynard_full(
    raw_dir: Path,
    output_h5ad: Path | None = None,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Download full Maynard et al. 2020 NSCLC cohort from CZ Biohub and construct AnnData.

    Profiles 27,489 single cells across 45 patients with advanced NSCLC.

    Args:
        raw_dir: Directory to store downloaded raw CSV matrices.
        output_h5ad: Optional path to save constructed H5AD (e.g. data/manual_download/Maynard_NSCLC.h5ad).
        force_download: If True, re-downloads even if files already exist locally.

    Returns:
        Result[ad.AnnData, str]: Parsed AnnData with complete clinical and cell lineage metadata.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    csv_mat = raw_dir / "S01_datafinal.csv"
    csv_meta = raw_dir / "S01_metacells.csv"

    # Also check scratch/ directory for existing copies if available
    scratch_mat = raw_dir.parent.parent / "scratch" / "S01_datafinal.csv"
    scratch_meta = raw_dir.parent.parent / "scratch" / "S01_metacells.csv"
    if not csv_mat.exists() and scratch_mat.exists():
        import shutil

        shutil.copy2(scratch_mat, csv_mat)
    if not csv_meta.exists() and scratch_meta.exists():
        import shutil

        shutil.copy2(scratch_meta, csv_meta)

    from ..download.gdrive import download_gdrive_file

    # File IDs from official czbiohub-sf/scell_lung_adenocarcinoma Google Drive repository
    gdrive_mat_id = "1CFXuwvTJqtr71GAKNagb-C7shrFK-jsd"
    gdrive_meta_id = "12g84ooMTsWV58JWbp2ql8TMfgS4ct8as"

    if force_download or not csv_meta.exists():
        res_meta = download_gdrive_file(gdrive_meta_id, csv_meta)
        if not is_successful(res_meta):
            return Failure(f"Failed to download Maynard metadata: {res_meta.failure()}")

    if force_download or not csv_mat.exists():
        res_mat = download_gdrive_file(gdrive_mat_id, csv_mat)
        if not is_successful(res_mat):
            return Failure(f"Failed to download Maynard expression matrix: {res_mat.failure()}")

    try:
        logger.info("Parsing Maynard metadata from %s...", csv_meta.name)
        meta_pl = pl.read_csv(csv_meta)
        drop_cols = [c for c in meta_pl.columns if not c or c == ""]
        if drop_cols:
            meta_pl = meta_pl.drop(drop_cols)

        meta_df = meta_pl.to_pandas()
        meta_df.index = meta_df["cell_id"].astype(str)

        # Harmonize standardized clinical columns
        meta_df["patient"] = meta_df["patient_id"].astype(str)
        meta_df["sample"] = meta_df["sample_name"].astype(str)
        if "biopsy_time_status" in meta_df.columns:
            meta_df["timepoint"] = meta_df["biopsy_time_status"].astype(str)
            status_clean = meta_df["biopsy_time_status"].str.upper()
            meta_df["response_binary"] = np.where(
                status_clean.isin(["PR", "CR", "RESPONDER"]),
                1.0,
                np.where(status_clean.isin(["PD", "NON-RESPONDER"]), 0.0, np.nan),
            )

        logger.info("Streaming and parsing Maynard expression matrix (27,489 cells x ~26,577 genes)...")
        genes: list[str] = []
        row_parts: list[np.ndarray] = []
        col_parts: list[np.ndarray] = []
        data_parts: list[np.ndarray] = []

        with open(csv_mat, "r") as f:
            header_line = f.readline().strip()
            cell_names = [c.strip('"') for c in header_line.split(",")][1:]
            n_cells = len(cell_names)

            for gene_idx, line in enumerate(f):
                comma_idx = line.find(",")
                if comma_idx == -1:
                    continue
                gene_name = line[:comma_idx].strip('"')
                genes.append(gene_name)

                row_arr = np.fromstring(line[comma_idx + 1:], sep=",", dtype=np.float32)
                nz = np.flatnonzero(row_arr)
                if len(nz) > 0:
                    col_parts.append(nz)
                    row_parts.append(np.full(len(nz), gene_idx, dtype=np.int32))
                    data_parts.append(row_arr[nz])

        all_rows = np.concatenate(row_parts) if row_parts else np.empty(0, dtype=np.int32)
        all_cols = np.concatenate(col_parts) if col_parts else np.empty(0, dtype=np.int32)
        all_data = np.concatenate(data_parts) if data_parts else np.empty(0, dtype=np.float32)
        del row_parts, col_parts, data_parts

        mat_genes_cells = sp.coo_matrix((all_data, (all_rows, all_cols)), shape=(len(genes), n_cells), dtype=np.float32)
        mat_cells_genes = mat_genes_cells.T.tocsr()
        del all_rows, all_cols, all_data, mat_genes_cells

        # Align obs metadata with matrix cell order
        obs_df = meta_df.reindex(cell_names)

        # Sanitize columns for robust H5AD serialization
        for col in obs_df.columns:
            if obs_df[col].dtype == object:
                obs_df[col] = obs_df[col].fillna("").astype(str)

        adata = ad.AnnData(
            X=mat_cells_genes,
            obs=obs_df,
            var=pd.DataFrame(index=genes),
        )
        adata.obs_names.name = "cell_id"
        adata.var_names.name = "gene_name"

        # Compute TME major lineage scoring to populate cell_type_major and cell_type
        try:
            from ..genesets import get_tme_major_lineage_collection, score_geneset_zscore

            col = get_tme_major_lineage_collection()
            scores_res = score_geneset_zscore(adata, col)
            match scores_res:
                case Success(scores_df):
                    feature_cols = [c for c in scores_df.columns if c != "sample_id"]
                    if feature_cols:
                        scores_np = scores_df.select(feature_cols).to_numpy()
                        best_indices = np.argmax(scores_np, axis=1)
                        predicted_major = [feature_cols[idx] for idx in best_indices]
                        adata.obs["cell_type_major"] = predicted_major
                        adata.obs["cell_type"] = predicted_major
                case _:
                    adata.obs["cell_type_major"] = "Unknown"
                    adata.obs["cell_type"] = "Unknown"
        except Exception as scoring_exc:
            logger.debug("TME lineage scoring on Maynard cohort deferred: %s", scoring_exc)
            adata.obs["cell_type"] = "Unknown"

        # Ensure all processed single-cell datasets use canonical Ensembl IDs
        try:
            from ..preprocessing.gene_normalization import normalize_genes_to_ensembl

            logger.info("Normalizing Maynard genes to Ensembl Release 111...")
            adata = normalize_genes_to_ensembl(adata)
        except Exception as norm_exc:
            logger.warning("Could not normalize Maynard genes to Ensembl: %s", norm_exc)

        if output_h5ad:
            output_h5ad.parent.mkdir(parents=True, exist_ok=True)
            logger.info("Serializing full Maynard NSCLC cohort to %s...", output_h5ad)
            adata.write_h5ad(output_h5ad)
            size_mb = output_h5ad.stat().st_size / (1024 * 1024)
            logger.info("Successfully serialized full Maynard NSCLC cohort to %s (%.1f MB)", output_h5ad.name, size_mb)

        return Success(adata)
    except Exception as exc:
        msg = f"Failed to build full Maynard dataset: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_maynard(
    repo_root: Path,
    raw_dir: Path | None = None,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load Maynard et al. 2020 NSCLC full cohort (27,489 cells)."""
    from ..paths import find_dataset_h5ad, get_manual_download_dir, get_raw_dataset_dir

    if not force_download:
        h5ad_maybe = find_dataset_h5ad("Maynard_NSCLC", repo_root=repo_root)
        match h5ad_maybe:
            case Some(h5ad_path):
                try:
                    logger.info("Loading Maynard NSCLC H5AD from %s...", h5ad_path)
                    adata = ad.read_h5ad(h5ad_path)
                    if adata.n_obs > 10000 or not auto_download:
                        return Success(adata)
                    logger.info(
                        "Found legacy 3k Maynard dataset (%d cells). Auto-downloading full 27,489-cell cohort...",
                        adata.n_obs,
                    )
                except Exception as exc:
                    logger.warning("Failed to load existing Maynard H5AD: %s. Re-building...", exc)
            case _:
                pass

    if auto_download:
        target_raw = raw_dir or get_raw_dataset_dir("Maynard_NSCLC", repo_root=repo_root)
        target_h5ad = get_manual_download_dir(repo_root=repo_root) / "Maynard_NSCLC.h5ad"
        return download_and_build_maynard_full(
            raw_dir=target_raw,
            output_h5ad=target_h5ad,
            force_download=force_download,
        )

    expected_file = get_manual_download_dir(repo_root=repo_root) / "Maynard_NSCLC.h5ad"
    msg = (
        f"Maynard NSCLC dataset not found in candidate paths. "
        f"Expected file at '{expected_file}'. "
        f"Set auto_download=True or run download_and_build_maynard_full() to download and construct it."
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

    # 1. Prefer 10x MTX matrix if already downloaded / converted
    mtx_candidates = [raw_dir / "matrix.mtx.gz", raw_dir / "matrix.mtx"]
    mtx_file = next((p for p in mtx_candidates if p.exists()), None)
    barcodes_candidates = [raw_dir / "barcodes.tsv.gz", raw_dir / "barcodes.tsv"]
    barcodes_file = next((p for p in barcodes_candidates if p.exists()), None)
    features_candidates = [
        raw_dir / "features.tsv.gz",
        raw_dir / "features.tsv",
        raw_dir / "genes.tsv.gz",
        raw_dir / "genes.tsv",
    ]
    features_file = next((p for p in features_candidates if p.exists()), None)
    meta_candidates = [
        raw_dir / "GSE179994_Tcell.metadata.tsv.gz",
        raw_dir / "GSE179994_meta.gz",
        raw_dir / "GSE179994_meta.tsv",
    ]
    actual_meta = next((p for p in meta_candidates if p.exists()), None)

    if mtx_file and barcodes_file and features_file:
        try:
            import scipy.io as sio

            logger.info("Reading GSE179994 MTX raw count matrix from %s...", mtx_file.name)
            open_bc = gzip.open if barcodes_file.suffix == ".gz" else open
            with open_bc(barcodes_file, "rt") as f:
                cell_names = [line.strip() for line in f if line.strip()]

            open_feat = gzip.open if features_file.suffix == ".gz" else open
            with open_feat(features_file, "rt") as f:
                gene_names = [line.strip().split("\t")[0] for line in f if line.strip()]

            mat_csr = sio.mmread(str(mtx_file)).tocsr()
            if mat_csr.shape[0] == len(gene_names) and mat_csr.shape[1] == len(cell_names):
                mat_csr = mat_csr.T.tocsr()

            obs_df = pd.DataFrame(index=cell_names)
            if actual_meta:
                df_meta = pd.read_csv(actual_meta, sep="\t", index_col=0)
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
            logger.info("Successfully loaded GSE179994 via MTX: %d cells x %d genes", adata.n_obs, adata.n_vars)
            return normalize_total_counts(adata)
        except Exception as exc:
            logger.warning("Failed to load GSE179994 via MTX: %s. Falling back to RDS...", exc)

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
            # Chunk column slicing to prevent 60+ GB dense allocation
            chunk_size = 25000
            csr_parts: list[sp.csr_matrix] = []
            for start in range(0, len(cell_names), chunk_size):
                end = min(start + chunk_size, len(cell_names))
                sub_arr = df_counts.iloc[:, start:end].to_numpy(dtype=np.float32).T
                csr_parts.append(sp.csr_matrix(sub_arr))
                del sub_arr
            mat_csr = sp.vstack(csr_parts, format="csr")
            del csr_parts
        else:  # cells x genes
            cell_names = list(df_counts.index)
            gene_names = list(df_counts.columns)
            chunk_size = 25000
            csr_parts = []
            for start in range(0, len(cell_names), chunk_size):
                end = min(start + chunk_size, len(cell_names))
                sub_arr = df_counts.iloc[start:end, :].to_numpy(dtype=np.float32)
                csr_parts.append(sp.csr_matrix(sub_arr))
                del sub_arr
            mat_csr = sp.vstack(csr_parts, format="csr")
            del csr_parts

        del df_counts, rds_dict
        import gc
        gc.collect()

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

