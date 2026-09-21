"""Script 3: Individual Quality Control, Metadata Harmonization, and Gene Normalization per Dataset.

Processes each raw downloaded single-cell dataset from subdirectories, handles 10X mtx triplets, GEO series matrix comments,
TPM expression matrices, applies QC, metadata harmonization, gene symbol normalization, and saves preprocessed .h5ad files.
"""

import concurrent.futures
import gc
import gzip
import io
import multiprocessing
import os
from pathlib import Path
from pydantic import BaseModel, ConfigDict
import anndata as ad
import polars as pl
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix, issparse
from returns.result import Result, Success, Failure
from single_cell_immuno_datasets import (
    TIER_1_DATASETS,
    QualityControlSpec,
    preprocess_anndata,
)


class PreprocessConfig(BaseModel):
    """Immutable configuration for dataset preprocessing into HDF5 files."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    output_dir: Path
    default_qc: QualityControlSpec = QualityControlSpec()
    workers: int = 8
    intra_threads: int = 5


def read_tpm_matrix(file_path: Path) -> ad.AnnData:
    """Reads a gene-by-cell or cell-by-gene expression table into sparse AnnData."""
    df = pl.read_csv(
        file_path,
        separator="\t" if "\t" in file_path.name or "tpm.txt" in file_path.name else ",",
        truncate_ragged_lines=True,
        ignore_errors=True,
        encoding="utf8-lossy"
    )
    first_col = df.columns[0]
    genes = df[first_col].cast(pl.String).to_list()
    cell_cols = df.columns[1:]
    mat = df.select(cell_cols).to_numpy().astype(np.float32).T
    sparse_mat = csr_matrix(mat, dtype=np.float32)
    gc.collect()
    return ad.AnnData(X=sparse_mat, obs=pd.DataFrame(index=cell_cols), var=pd.DataFrame(index=genes))


def _read_single_10x_h5(h5_file: Path) -> Result[ad.AnnData, str]:
    """Helper to read a single 10X .h5 file."""
    try:
        adata_sub = sc.read_10x_h5(h5_file)
        if not issparse(adata_sub.X):
            adata_sub.X = csr_matrix(adata_sub.X, dtype=np.float32)
        adata_sub.var_names_make_unique()
        sample_id = h5_file.stem.split("_filtered")[0]
        adata_sub.obs["sample"] = sample_id
        return Success(adata_sub)
    except Exception as ex:
        return Failure(f"Error reading H5 file {h5_file.name}: {ex}")


def read_10x_h5_dir(parent_dir: Path, max_workers: int = 8) -> Result[ad.AnnData, str]:
    """Reads all 10X filtered_gene_bc_matrices .h5 files in a directory in parallel and concatenates them."""
    h5_files = sorted(
        f for f in parent_dir.glob("*.h5")
        if f.is_file() and f.stat().st_size > 10_000
    )
    if not h5_files:
        return Failure(f"No .h5 files found in {parent_dir}")

    workers = min(max_workers, max(1, len(h5_files)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_read_single_10x_h5, h5_files))

    adatas = [r.unwrap() for r in results if isinstance(r, Success)]
    for r in results:
        if isinstance(r, Failure):
            print(r.failure())

    if not adatas:
        return Failure(f"Failed to parse any .h5 files in {parent_dir}")

    combined = ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    return Success(combined)




def read_gse120575_matrix(file_path: Path) -> ad.AnnData:
    """Parses GSE120575 TPM text matrix with non-standard header lines using multithreaded Polars."""
    try:
        with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as f:
            line1 = f.readline().strip("\r\n").split("\t")
            _line2 = f.readline()
    except (gzip.BadGzipFile, UnicodeDecodeError):
        with open(file_path, "rt", encoding="utf-8", errors="replace") as f:
            line1 = f.readline().strip("\r\n").split("\t")
            _line2 = f.readline()

    raw_cell_ids = line1 if (line1[0].startswith("A") or line1[0].startswith("P")) else line1[1:]
    cell_ids = [c for c in raw_cell_ids if c.strip()]

    df = pl.read_csv(
        file_path,
        separator="\t",
        has_header=False,
        skip_rows=2,
        truncate_ragged_lines=True,
        ignore_errors=True,
        encoding="utf8-lossy",
    )
    gene_col = df.columns[0]
    genes = df[gene_col].cast(pl.String).to_list()
    data_cols = df.columns[1 : 1 + len(cell_ids)]
    matched_cell_ids = cell_ids[: len(data_cols)]
    mat = df.select(data_cols).to_numpy().astype(np.float32).T
    sparse_mat = csr_matrix(mat, dtype=np.float32)
    gc.collect()
    return ad.AnnData(X=sparse_mat, obs=pd.DataFrame(index=matched_cell_ids), var=pd.DataFrame(index=genes))


def attach_gse120575_response_metadata(adata: ad.AnnData, parent_dir: Path) -> ad.AnnData:
    """Parses GSE120575 GEO metadata template to extract cell title → response mapping.

    The metadata file has comment lines starting with # and a data section starting with
    'Sample name\ttitle\t...characteristics: response\t...' (one row per cell).
    """
    meta_file = parent_dir / "GSE120575_patient_ID_single_cells.txt.gz"
    if not meta_file.exists():
        return adata
    try:
        rows = []
        header: list[str] = []
        with gzip.open(meta_file, "rt", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = [p.strip().strip('"') for p in stripped.split("\t")]
                if parts[0] == "Sample name":
                    header = parts
                    continue
                if header:
                    # Pad short rows to header length so DataFrame constructor succeeds
                    padded = parts + [""] * max(0, len(header) - len(parts))
                    rows.append(padded[:len(header)])

        if not header or not rows:
            return adata

        meta_df = pd.DataFrame(rows, columns=header)
        title_col = "title"
        response_col = next(
            (c for c in meta_df.columns if "response" in c.lower() and "characteristics" in c.lower()),
            None
        )
        patient_col = next(
            (c for c in meta_df.columns if "patient" in c.lower() or "patinet" in c.lower()),
            None
        )
        if title_col in meta_df.columns and response_col:
            # Build title→response lookup (drop rows with empty title)
            lookup = meta_df[meta_df[title_col].str.len() > 0].set_index(title_col)
            # Use pd.Series.map — obs_names is a pandas Index; convert to Series first
            obs_series = pd.Series(adata.obs_names, index=adata.obs_names)
            adata.obs["response"] = obs_series.map(lookup[response_col].to_dict()).fillna("Unknown").values
            if patient_col and patient_col in lookup.columns:
                adata.obs["patient_geo"] = obs_series.map(lookup[patient_col].to_dict()).fillna("Unknown").values
    except Exception as ex:
        print(f"GSE120575 response metadata attach warning: {ex}")
    return adata





def read_geo_series_matrix(file_path: Path) -> ad.AnnData:
    """Parses a GEO series matrix file ignoring comment lines starting with ! into sparse AnnData."""
    lines: list[str] = []
    in_matrix = False
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_matrix = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_matrix or (not line.startswith("!") and not line.startswith("#") and not line.startswith("^")):
                lines.append(line)
    
    content = "".join(lines)
    df = pl.read_csv(
        io.StringIO(content),
        separator="\t",
        truncate_ragged_lines=True,
        ignore_errors=True
    )
    first_col = df.columns[0]
    genes = df[first_col].cast(pl.String).to_list()
    cell_cols = df.columns[1:]
    mat = df.select(cell_cols).to_numpy().astype(np.float32).T
    sparse_mat = csr_matrix(mat, dtype=np.float32)
    gc.collect()
    return ad.AnnData(X=sparse_mat, obs=pd.DataFrame(index=cell_cols), var=pd.DataFrame(index=genes))


def _read_single_gondal_csv(file_path: Path, csv_name: str) -> Result[ad.AnnData, str]:
    """Helper to read a single sub-cohort CSV from the Gondal Zip archive."""
    import zipfile
    meta_numeric = {'Unnamed: 0', 'nCount_RNA', 'nFeature_RNA', 'Var.41'}
    try:
        with zipfile.ZipFile(file_path, "r") as z:
            with z.open(csv_name) as f:
                df = pd.read_csv(f)
                cell_ids = df['cell_id'].astype(str).tolist() if 'cell_id' in df.columns else df.index.astype(str).tolist()
                num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                gene_columns = [c for c in num_cols if c not in meta_numeric]
                obs_columns = [c for c in df.columns if c not in gene_columns]

                mat = df[gene_columns].to_numpy(dtype=np.float32)
                sparse_mat = csr_matrix(mat, dtype=np.float32)
                obs_df = df[obs_columns].copy()
                obs_df.index = pd.Index(cell_ids)

                return Success(ad.AnnData(
                    X=sparse_mat,
                    obs=obs_df,
                    var=pd.DataFrame(index=gene_columns)
                ))
    except Exception as ex:
        return Failure(f"Error reading Gondal CSV {csv_name}: {ex}")


def read_gondal_zip(file_path: Path, max_workers: int = 8) -> ad.AnnData:
    """Parses Gondal2025 Zip archive containing cohort CSV tables into a combined AnnData object concurrently."""
    import zipfile
    with zipfile.ZipFile(file_path, "r") as z:
        csv_names = [name for name in z.namelist() if name.endswith(".csv")]

    if not csv_names:
        raise ValueError(f"No CSV files found in Zip archive {file_path}")

    workers = min(max_workers, len(csv_names))
    print(f"Reading {len(csv_names)} Gondal sub-cohort CSVs concurrently ({workers} threads)...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(lambda name: _read_single_gondal_csv(file_path, name), csv_names))

    adatas = [r.unwrap() for r in results if isinstance(r, Success)]
    if not adatas:
        raise ValueError(f"Failed to load any CSV files from Zip archive {file_path}")

    combined = ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    return combined


def _read_single_smartseq2_sample(gsm_file: Path) -> Result[ad.AnnData, str]:
    """Helper to read a single GSM Smart-seq2 expression table."""
    import re
    m = re.match(r"GSM\d+_([^.]+)", gsm_file.name)
    patient_id = m.group(1) if m else gsm_file.stem
    try:
        df = pl.read_csv(
            gsm_file, separator="\t", truncate_ragged_lines=True,
            ignore_errors=True, encoding="utf8-lossy"
        )
        first_col = df.columns[0]
        genes = df[first_col].cast(pl.String).to_list()
        cell_cols = df.columns[1:]
        mat = df.select(cell_cols).to_numpy().astype(np.float32).T
        sparse_mat = csr_matrix(mat, dtype=np.float32)
        obs_df = pd.DataFrame({"patient": patient_id}, index=list(cell_cols))
        return Success(ad.AnnData(X=sparse_mat, obs=obs_df, var=pd.DataFrame(index=genes)))
    except Exception as ex:
        return Failure(f"Error reading per-sample file {gsm_file.name}: {ex}")


def read_per_sample_smartseq2_dir(parent_dir: Path, max_workers: int = 8) -> Result[ad.AnnData, str]:
    """Reads GSE123139-style per-sample Smart-seq2 files in parallel where each GSM*.txt.gz is one patient."""
    gsm_files = sorted(
        f for f in parent_dir.glob("GSM*.txt.gz")
        if f.is_file() and f.stat().st_size > 1_000
    )
    if not gsm_files:
        return Failure(f"No per-sample GSM*.txt.gz files found in {parent_dir}")

    workers = min(max_workers, max(1, len(gsm_files)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_read_single_smartseq2_sample, gsm_files))

    adatas = [r.unwrap() for r in results if isinstance(r, Success)]
    for r in results:
        if isinstance(r, Failure):
            print(r.failure())

    if not adatas:
        return Failure(f"Failed to load any per-sample files from {parent_dir}")

    combined = ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    return Success(combined)




def _read_single_10x_triplet(
    mat_file: Path,
    parent_dir: Path,
    single_matrix: bool
) -> Result[ad.AnnData, str]:
    """Helper to parse a single 10X sparse matrix triplet (matrix, barcodes, genes)."""
    prefix = mat_file.name
    for suffix in ["_matrix.mtx.gz", "_matrix.gz", ".matrix.mtx.gz", "matrix.mtx.gz", ".mtx.gz"]:
        if prefix.endswith(suffix):
            prefix = prefix[:-len(suffix)]
            break

    # Find barcodes file
    barcode_candidates = [
        f for f in parent_dir.glob("*")
        if f.is_file() and ("barcode" in f.name.lower()) and (prefix in f.name or single_matrix)
    ]
    # Find genes/features file
    gene_candidates = [
        f for f in parent_dir.glob("*")
        if f.is_file() and ("gene" in f.name.lower() or "feature" in f.name.lower()) and (prefix in f.name or single_matrix)
    ]

    try:
        adata = sc.read_mtx(mat_file).T
        if barcode_candidates:
            bc_df = pd.read_csv(barcode_candidates[0], header=None)
            adata.obs_names = bc_df[0].astype(str).tolist()
            adata.obs_names_make_unique()
        if gene_candidates:
            gene_df = pd.read_csv(gene_candidates[0], header=None, sep="\t" if "\t" in gene_candidates[0].name or ".tsv" in gene_candidates[0].name else r"\s+")
            col_idx = 1 if gene_df.shape[1] > 1 else 0
            adata.var_names = gene_df[col_idx].astype(str).tolist()
            adata.var_names_make_unique()

        if not issparse(adata.X):
            adata.X = csr_matrix(adata.X, dtype=np.float32)
        # Tag each cell with its sample of origin when multiple triplets exist
        if prefix:
            adata.obs["sample"] = prefix
        return Success(adata)
    except Exception as ex:
        return Failure(f"Error reading 10X matrix triplet for {mat_file.name}: {ex}")


def read_10x_matrix_dir(parent_dir: Path, max_workers: int = 8) -> Result[ad.AnnData, str]:
    """Finds and parses 10X sparse matrix triplets (matrix, barcodes, genes/features) in parent_dir."""
    matrix_files = [
        f for f in parent_dir.glob("*")
        if f.is_file() and ("matrix" in f.name.lower() or "mtx" in f.name.lower())
        and not f.name.endswith(".txt") and not f.name.endswith(".txt.gz")
        and f.stat().st_size > 10_000
    ]
    if not matrix_files:
        return Failure(f"No 10X matrix files found in {parent_dir}")

    if len(matrix_files) == 1:
        return _read_single_10x_triplet(matrix_files[0], parent_dir, single_matrix=True)

    workers = min(max_workers, len(matrix_files))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(
            executor.map(
                lambda mf: _read_single_10x_triplet(mf, parent_dir, single_matrix=False),
                matrix_files,
            )
        )

    adatas = [r.unwrap() for r in results if isinstance(r, Success)]
    if not adatas:
        return Failure(f"Failed to parse 10X matrix triplets in {parent_dir}")

    combined = ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    return Success(combined)



def attach_metadata_files(adata: ad.AnnData, parent_dir: Path) -> ad.AnnData:
    """Finds and merges external cell metadata tables (.csv.gz, .tsv.gz, .txt.gz, .anno, .meta) onto adata.obs."""
    meta_files = [
        f for f in parent_dir.glob("*")
        if f.is_file() and any(k in f.name.lower() for k in ["meta", "anno", "patient", "cell.annotations"])
        and f.stat().st_size > 100
        and not f.name.endswith(".h5ad") and not f.name.endswith(".parquet")
        and not f.name.endswith(".h5")
    ]

    # Extended set of recognized cell-ID column names (case-insensitive match)
    _CELL_ID_COLS = frozenset({
        "cell_id", "cellid", "cell_name", "cell.name", "cell_ids",
        "barcode", "barcodes", "cell_id_single_cells", "cell", "cells", "samples",
        "Cell_ID",
    })

    for mf in meta_files:
        try:
            # CSV if extension is .csv.gz — never override based on filename keywords alone
            is_csv_ext = mf.name.endswith(".csv.gz") or mf.name.endswith(".csv")
            is_tsv = not is_csv_ext and (".tsv" in mf.name or "txt" in mf.name or "_meta" in mf.name or "_anno" in mf.name)
            sep = "\t" if is_tsv else ","
            df_meta = pd.read_csv(
                mf, sep=sep, compression="gzip" if mf.name.endswith(".gz") else None,
                low_memory=False, encoding_errors="replace"
            )

            cell_col = next(
                (col for col in df_meta.columns if col in _CELL_ID_COLS or col.lower() in _CELL_ID_COLS),
                None
            )

            if cell_col:
                df_meta[cell_col] = df_meta[cell_col].astype(str)
                df_meta_indexed = df_meta.drop_duplicates(subset=[cell_col]).set_index(cell_col)
                cols_to_add = [c for c in df_meta_indexed.columns if c not in adata.obs.columns]
                if cols_to_add:
                    adata.obs = adata.obs.join(df_meta_indexed[cols_to_add], how="left")
            elif len(df_meta) == adata.n_obs:
                cols_to_add = [c for c in df_meta.columns if c not in adata.obs.columns]
                for col in cols_to_add:
                    adata.obs[col] = df_meta[col].values
        except Exception as ex:
            print(f"Metadata attach warning for {mf.name}: {ex}")

    return adata


def _derive_patient_from_col_prefix(adata: ad.AnnData) -> ad.AnnData:
    """For datasets like GSE123813 where obs names are '<cancer>.<patient>.<pre/post>.<celltype>_<barcode>',
    extracts the patient ID token (second dot-delimited segment) and sets it as adata.obs['patient']."""
    import re
    if "patient" in adata.obs.columns and adata.obs["patient"].nunique() > 1:
        return adata
    sample = adata.obs_names[0] if adata.n_obs > 0 else ""
    # Pattern: word.patient_id.word... OR word-patient_id_barcode
    m = re.match(r"[a-z]+\.([a-z0-9]+)\.", sample, re.IGNORECASE)
    if m:
        patients = adata.obs_names.to_series().str.extract(r"[a-z]+\.([a-z0-9]+)\.", expand=False, flags=re.IGNORECASE)
        if patients.notna().sum() > 0.3 * adata.n_obs:
            adata.obs["patient"] = patients.fillna("Unknown").values
    return adata


def attach_gse123139_metadata(adata: ad.AnnData, parent_dir: Path) -> ad.AnnData:
    """Attaches true patient IDs and treatment/response status from GSE123139_matrix.txt.gz."""
    if "gse123139" not in str(parent_dir).lower():
        return adata
    matrix_files = list(parent_dir.glob("*matrix.txt.gz"))
    if not matrix_files:
        return adata
    try:
        mf = matrix_files[0]
        batches: list[str] = []
        patients: list[str] = []
        with gzip.open(mf, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_characteristics_ch1"):
                    parts = [p.strip().strip('"') for p in line.strip().split("\t")]
                    if len(parts) > 1:
                        if "amplification batch" in parts[1]:
                            batches = [p.split(":")[-1].strip() for p in parts[1:]]
                        elif "patient id" in parts[1]:
                            patients = [p.split(":")[-1].strip() for p in parts[1:]]

        if batches and patients and len(batches) == len(patients):
            batch_to_full_id = dict(zip(batches, patients))

            def get_condition(full_id: str) -> str:
                if "-N" in full_id:
                    return "Treatment-naive"
                elif "IT" in full_id:
                    return "Immunotherapy"
                elif "-T" in full_id:
                    return "Other-treatment"
                return "Unknown"

            if "patient" in adata.obs.columns:
                batch_series = adata.obs["patient"].astype(str)
                adata.obs["batch"] = batch_series
                full_ids = batch_series.map(batch_to_full_id).fillna(batch_series)
                clean_pats = full_ids.apply(lambda s: s.split("-")[0] if s.startswith("p") else s)
                adata.obs["patient"] = clean_pats
                adata.obs["response"] = full_ids.apply(get_condition)
                print(f"GSE123139: Mapped {len(batch_to_full_id)} batches to {adata.obs['patient'].nunique()} patients with response conditions: {adata.obs['response'].value_counts().to_dict()}")
    except Exception as ex:
        print(f"GSE123139 metadata parse warning: {ex}")

    return adata


def attach_gse179994_metadata(adata: ad.AnnData, parent_dir: Path) -> ad.AnnData:
    """For GSE179994: extracts Pre-treatment vs On-treatment condition from sample or cell metadata."""
    if "gse179994" not in str(parent_dir).lower():
        return adata
    try:
        if "sample" in adata.obs.columns:
            s = adata.obs["sample"].astype(str)
            is_post = s.str.contains(r"\.post|\.tr", regex=True, case=False)
            is_pre = s.str.contains(r"\.pre|\.ut", regex=True, case=False)
            cond = np.where(is_post, "On-treatment", np.where(is_pre, "Pre-treatment", "Unknown"))
            adata.obs["response"] = cond
            adata.obs["patient_biopsy"] = s
            # Use per-biopsy sample ID as patient/sample to avoid duplicate-sample collision in milopy GLM
            adata.obs["patient"] = s
            print(f"GSE179994: Mapped {adata.obs['patient'].nunique()} biopsy samples to response conditions: {pd.Series(cond).value_counts().to_dict()}")
        elif "cellid" in adata.obs.columns:
            c = adata.obs["cellid"].astype(str)
            is_post = c.str.contains(r"\.tr\.", regex=True, case=False)
            is_pre = c.str.contains(r"\.ut\.", regex=True, case=False)
            cond = np.where(is_post, "On-treatment", np.where(is_pre, "Pre-treatment", "Unknown"))
            adata.obs["response"] = cond
    except Exception as ex:
        print(f"GSE179994 metadata parse warning: {ex}")

    return adata


def load_raw_as_anndata(raw_file: Path, max_workers: int = 8) -> Result[ad.AnnData, str]:
    """Loads raw dataset file (.h5ad, .parquet, 10X .h5, 10X mtx, GEO series matrix, or count tables) into AnnData."""
    try:
        name = raw_file.name.lower()
        parent_dir = raw_file.parent

        if "gse120575" in name or "gse120575" in parent_dir.name.lower():
            adata = read_gse120575_matrix(raw_file)
            adata = attach_gse120575_response_metadata(adata, parent_dir)
            return Success(attach_metadata_files(adata, parent_dir))
        elif "gondal" in name or "gondal" in parent_dir.name.lower() or name.endswith(".zip"):
            adata = read_gondal_zip(raw_file, max_workers=max_workers)
            return Success(attach_metadata_files(adata, parent_dir))

        # Priority 1: Per-sample H5 files (e.g. GSE159115 with GSM*.h5 per patient)
        h5_res = read_10x_h5_dir(parent_dir, max_workers=max_workers)
        if isinstance(h5_res, Success):
            adata = attach_metadata_files(h5_res.unwrap(), parent_dir)
            # Promote sample tag from H5 reader to patient
            if "sample" in adata.obs.columns and ("patient" not in adata.obs.columns or adata.obs["patient"].nunique() <= 1):
                adata.obs["patient"] = adata.obs["sample"]
            return Success(adata)

        # Priority 2: Per-sample Smart-seq2 GSM*.txt.gz files (e.g. GSE123139)
        gsm_res = read_per_sample_smartseq2_dir(parent_dir, max_workers=max_workers)
        if isinstance(gsm_res, Success):
            adata = attach_metadata_files(gsm_res.unwrap(), parent_dir)
            adata = attach_gse123139_metadata(adata, parent_dir)
            return Success(adata)

        # Priority 3: 10X sparse matrix triplets
        mtx_res = read_10x_matrix_dir(parent_dir, max_workers=max_workers)
        if isinstance(mtx_res, Success):
            adata = mtx_res.unwrap()
            adata.obs_names_make_unique()
            adata = attach_metadata_files(adata, parent_dir)
            # Derive patient from obs index prefix pattern (e.g. GSE123813 bcc.su001.pre...)
            adata = _derive_patient_from_col_prefix(adata)
            adata = attach_gse179994_metadata(adata, parent_dir)
            return Success(adata)

        if name.endswith(".h5ad"):
            adata = ad.read_h5ad(raw_file)
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X, dtype=np.float32)
            return Success(attach_metadata_files(adata, parent_dir))
        elif name.endswith(".parquet"):
            df = pl.read_parquet(raw_file)
            first_col = df.columns[0]
            genes = df[first_col].cast(pl.String).to_list()
            cell_cols = df.columns[1:]
            mat = df.select(cell_cols).to_numpy().astype(np.float32).T
            sparse_mat = csr_matrix(mat, dtype=np.float32)
            gc.collect()
            adata = ad.AnnData(X=sparse_mat, obs=pd.DataFrame(index=cell_cols), var=pd.DataFrame(index=genes))
            return Success(attach_metadata_files(adata, parent_dir))
        elif "series_matrix" in name or "matrix.txt" in name:
            if raw_file.stat().st_size < 100_000:
                return Failure(f"Skipping GEO soft header file: {raw_file.name}")
            adata = read_geo_series_matrix(raw_file)
            return Success(attach_metadata_files(adata, parent_dir))
        elif "mtx" in name:
            try:
                adata = sc.read_10x_mtx(parent_dir, var_names="gene_symbols", make_unique=True)
                if not issparse(adata.X):
                    adata.X = csr_matrix(adata.X, dtype=np.float32)
                return Success(attach_metadata_files(adata, parent_dir))
            except Exception:
                adata = read_tpm_matrix(raw_file)
                return Success(attach_metadata_files(adata, parent_dir))
        elif any(ext in name for ext in ["tpm", "counts", "matrix", "txt", "csv"]):
            try:
                adata = read_tpm_matrix(raw_file)
                return Success(attach_metadata_files(adata, parent_dir))
            except Exception:
                adata = sc.read_text(raw_file).T
                if not issparse(adata.X):
                    adata.X = csr_matrix(adata.X, dtype=np.float32)
                return Success(attach_metadata_files(adata, parent_dir))
        else:
            return Failure(f"Unsupported format: {raw_file.name}")
    except Exception as e:
        return Failure(f"Failed to load raw dataset {raw_file.name}: {str(e)}")


def preprocess_single_raw_dataset(
    raw_dir_or_file: Path,
    out_dir: Path,
    qc_spec: QualityControlSpec,
    max_workers: int = 8,
) -> Result[Path, str]:
    """Preprocesses a raw dataset directory or file matching a registered Tier 1 specification."""
    accession = raw_dir_or_file.parent.name if raw_dir_or_file.is_file() else raw_dir_or_file.name
    spec_matches = [s for s in TIER_1_DATASETS if s.accession.lower() in accession.lower() or accession.lower() in s.accession.lower()]
    if not spec_matches:
        spec_matches = [TIER_1_DATASETS[0]]
    
    spec = spec_matches[0]
    out_h5ad = out_dir / f"{spec.accession}_processed.h5ad"

    if out_h5ad.exists() and out_h5ad.stat().st_size > 5_000_000:
        print(f"Skipping '{spec.accession}' - already preprocessed: {out_h5ad.name} ({out_h5ad.stat().st_size // 1024} KB)")
        return Success(out_h5ad)

    print(f"Preprocessing raw dataset '{spec.accession}' from {raw_dir_or_file.name}...")
    
    match load_raw_as_anndata(raw_dir_or_file, max_workers=max_workers):
        case Success(adata):
            res = preprocess_anndata(adata, spec, qc_spec, out_h5ad)
            gc.collect()
            return res
        case Failure(err):
            return Failure(err)


def _worker_preprocess_dataset(
    args: tuple[Path, Path, QualityControlSpec, int]
) -> Result[Path, str]:
    """Top-level worker function executed inside ProcessPoolExecutor processes."""
    raw_file, out_dir, qc_spec, intra_threads = args
    os.environ["POLARS_MAX_THREADS"] = str(intra_threads)
    os.environ["OMP_NUM_THREADS"] = str(intra_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(intra_threads)
    os.environ["MKL_NUM_THREADS"] = str(intra_threads)
    return preprocess_single_raw_dataset(raw_file, out_dir, qc_spec, max_workers=intra_threads)


def run_dataset_preprocessing_pipeline(config: PreprocessConfig) -> Result[list[Path], str]:
    """Scans raw input directory recursively and preprocesses all single-cell datasets into .h5ad files in parallel."""
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    raw_files = [
        f for f in config.input_dir.rglob("*")
        if f.is_file() and any(f.name.endswith(ext) for ext in [".h5ad", ".txt.gz", ".csv.gz", ".gz", ".parquet", ".zip"])
        and not f.name.startswith(".") and not f.name.endswith(".txt")
    ]

    if not raw_files:
        return Failure(f"No raw expression dataset files found in input directory: {config.input_dir}")

    # Group files by parent directory, sorting each group so primary expression files come first
    _EXPR_KEYWORDS = frozenset({"tpm", "counts", "matrix", "mtx", "zenodo", ".h5", ".zip", "expression", "raw.tar"})
    _META_KEYWORDS = frozenset({"meta", "anno", "annotation", "patient", "barcode", "features", "feature", "genes", "_gene", "sra"})

    def _file_priority(f: Path) -> int:
        """Lower number = higher priority = picked as the representative file for its directory."""
        n = f.name.lower()
        if any(k in n for k in _EXPR_KEYWORDS):
            return 0
        if any(k in n for k in _META_KEYWORDS):
            return 2
        return 1

    from itertools import groupby
    from operator import attrgetter

    grouped = {}
    for f in raw_files:
        grouped.setdefault(f.parent, []).append(f)

    # Pick the highest-priority file from each directory group as the representative
    unique_candidates: list[Path] = []
    seen_dirs: set[Path] = set()
    for parent_dir, files_in_dir in grouped.items():
        if parent_dir not in seen_dirs:
            seen_dirs.add(parent_dir)
            best_file = min(files_in_dir, key=_file_priority)
            unique_candidates.append(best_file)

    print(f"Found {len(unique_candidates)} unique dataset candidate directories in {config.input_dir}.")

    if config.workers <= 1:
        processed_paths: list[Path] = []
        for candidate in unique_candidates:
            match preprocess_single_raw_dataset(candidate, config.output_dir, config.default_qc, max_workers=config.intra_threads):
                case Success(p):
                    processed_paths.append(p)
                    print(f"Saved preprocessed dataset -> {p}")
                case Failure(err):
                    print(f"Warning preprocessing {candidate.name}: {err}")
        return Success(processed_paths)

    worker_args = [
        (cand, config.output_dir, config.default_qc, config.intra_threads)
        for cand in unique_candidates
    ]
    effective_workers = min(config.workers, max(1, len(worker_args)))
    print(f"Starting parallel preprocessing with {effective_workers} worker processes (intra_threads={config.intra_threads})...")

    processed_paths: list[Path] = []
    ctx = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(max_workers=effective_workers, mp_context=ctx) as executor:
        future_to_cand = {
            executor.submit(_worker_preprocess_dataset, arg): arg[0]
            for arg in worker_args
        }
        for future in concurrent.futures.as_completed(future_to_cand):
            cand = future_to_cand[future]
            try:
                res = future.result()
                match res:
                    case Success(p):
                        processed_paths.append(p)
                        print(f"Successfully processed {cand.name} -> {p}")
                    case Failure(err):
                        print(f"Warning preprocessing {cand.name}: {err}")
            except Exception as exc:
                print(f"Worker crashed on {cand.name}: {exc}")

    return Success(processed_paths)


def main() -> None:
    import argparse
    default_cpu = os.cpu_count() or 4
    default_workers = min(8, default_cpu)
    default_intra_threads = max(1, default_cpu // default_workers)

    parser = argparse.ArgumentParser(
        description="Preprocess single-cell datasets into individual HDF5 (.h5ad) files with multi-core parallelism."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="/storage/halu/data/raw",
        help="Input directory with raw dataset files",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="/storage/halu/data/preprocessed",
        help="Output directory for preprocessed .h5ad files",
    )
    parser.add_argument(
        "--workers",
        "-j",
        type=int,
        default=default_workers,
        help=f"Number of parallel dataset worker processes (default: {default_workers})",
    )
    parser.add_argument(
        "--intra-threads",
        type=int,
        default=default_intra_threads,
        help=f"Number of threads per worker for intra-dataset loading & Polars (default: {default_intra_threads})",
    )
    args = parser.parse_args()

    config = PreprocessConfig(
        input_dir=Path(args.input_dir).resolve(),
        output_dir=Path(args.out_dir).resolve(),
        workers=max(1, args.workers),
        intra_threads=max(1, args.intra_threads),
    )
    print(f"Starting single-cell dataset preprocessing pipeline...")
    print(f"Input dir: {config.input_dir}")
    print(f"Output dir: {config.output_dir}")
    print(f"Dataset worker processes: {config.workers}")
    print(f"Threads per worker: {config.intra_threads}")
    print(f"Total core allocation: {config.workers * config.intra_threads} cores")

    match run_dataset_preprocessing_pipeline(config):
        case Success(paths):
            print(f"\nSuccessfully preprocessed {len(paths)} datasets into .h5ad files:")
            for p in paths:
                print(f" - {p}")
        case Failure(err):
            print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()


