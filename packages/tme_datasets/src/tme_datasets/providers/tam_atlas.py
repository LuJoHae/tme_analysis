"""Coulton et al. (2024) Pan-Cancer Tumor-Associated Macrophage (TAM) Single-Cell Atlas provider.

Publication:
    Coulton et al., "Using a pan-cancer atlas to investigate tumour associated macrophages
    as regulators of immunotherapy response", Nature Communications 15, 5665 (2024).
    DOI: 10.1038/s41467-024-49885-8 | PMID: 38969631

Zenodo Deposit:
    Record 11222158 (DOI: 10.5281/zenodo.11222158)
    File: mac.atlas.zenodo.200524.rds (19,748,833,115 bytes, MD5: 49d3bea7addd7fcd07d8e3843ca901c2)

Assay Architecture Note:
    This Zenodo atlas intentionally omits raw integer UMI counts (`RNA@counts` has 0 non-zeros;
    raw counts reside in the 32 constituent source studies) and provides three processed assays:
    1. `RNA@data`        -> `.X` and `.layers["log1p_norm"]` (sparse CSR `float32`, log-normalized)
    2. `SCT@data`        -> `.layers["sct_data"]` (sparse CSR `float32`, SCTransform log-corrected)
    3. `integrated@data` -> `.obsm["X_rpca_integrated"]` (`float32` 2,000-anchor RPCA matrix,
                            with feature names in `.uns["rpca_integrated_var_names"]`)
"""

from __future__ import annotations

from contextlib import ExitStack
import gc
from pathlib import Path
from typing import Any, Final, Literal, Mapping, Sequence, assert_never

import anndata as ad  # type: ignore[import-untyped]
import anndata.io  # type: ignore[import-untyped]
import h5py  # type: ignore[import-untyped]
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import scipy.sparse as sp

from ..download.fetcher import download_single_file
from ..logging import get_logger
from ..models import ChecksumSpec
from ..paths import (
    find_dataset_h5ad,
    find_repo_root,
    get_ensembl_dir,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_scratch_dataset_dir,
)
from ..preprocessing.gene_normalization import (
    _build_ensembl_projection_operator,
    _load_or_update_mapping_cache,
    _project_sparse_chunk_to_ensembl,
    ensure_ensembl_release_installed,
)
from ..storage.incremental_writer import H5ADSparseIncrementalWriter
from ..storage.rds_stream import (
    StagedAssaySpec,
    stream_seurat_rds_to_staging_h5,
)
from .atlas_single_cell import _apply_subset_and_subsample

logger = get_logger("providers.tam_atlas")

COULTON_TAM_DATASET_ID: Final[str] = "Coulton_TAM_Atlas"
COULTON_TAM_ZENODO_RECORD: Final[str] = "11222158"
COULTON_TAM_ZENODO_DOI: Final[str] = "10.5281/zenodo.11222158"
COULTON_TAM_PUBLICATION_DOI: Final[str] = "10.1038/s41467-024-49885-8"
COULTON_TAM_PUBLICATION_PMID: Final[str] = "38969631"
COULTON_TAM_RDS_FILENAME: Final[str] = "mac.atlas.zenodo.200524.rds"
COULTON_TAM_ZENODO_URL: Final[str] = (
    "https://zenodo.org/records/11222158/files/mac.atlas.zenodo.200524.rds?download=1"
)
COULTON_TAM_EXPECTED_MD5: Final[str] = "49d3bea7addd7fcd07d8e3843ca901c2"
COULTON_TAM_EXPECTED_SIZE_BYTES: Final[int] = 19_748_833_115

COULTON_TAM_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash=COULTON_TAM_EXPECTED_MD5,
    size_bytes=Some(COULTON_TAM_EXPECTED_SIZE_BYTES),
)

# Canonical 24 cluster labels (`0`..`23`) from Coulton et al. 2024 (`figures/mac.util.R`, res=0.525)
COULTON_TAM_CLUSTER_LABELS: Final[Mapping[str, str]] = {
    "0": "0_AlvMac",
    "1": "1_MetM2Mac",
    "2": "2_C3Mac",
    "3": "3_ICIMac1",
    "4": "4_ICIMac2",
    "5": "5_StressMac",
    "6": "6_SPP1AREGMac",
    "7": "7_IFNMac",
    "8": "8_IFNGMac",
    "9": "9_AngioMac",
    "10": "10_InflamMac",
    "11": "11_MetalloMac",
    "12": "12_MBMMac",
    "13": "13_CalciumMac",
    "14": "14_ProliMac",
    "15": "15_LYZMac",
    "16": "16_ECMHomeoMac",
    "17": "17_IFNMac3",
    "18": "18_ECMMac",
    "19": "19_ClassMono",
    "20": "20_TDoub",
    "21": "21_HemeMac",
    "22": "22_IFNMac4",
    "23": "23_NA",
}

_CANCER_CODE_NORMALIZATION: Final[Mapping[str, str]] = {
    "thyroid": "THCA",
    "melanoma": "MEL",
    "Melanoma": "MEL",
    "MEL": "MEL",
    "breast": "BC",
    "Breast": "BC",
    "BRCA": "BC",
    "BC": "BC",
    "liver": "LIHC",
    "Liver": "LIHC",
    "HCC": "LIHC",
    "LIHC": "LIHC",
    "ESO": "ESCA",
    "ESCA": "ESCA",
}

_CANCER_CODE_TO_FULL_NAME: Final[Mapping[str, str]] = {
    "BC": "Breast Cancer",
    "CRC": "Colorectal Cancer",
    "MEL": "Melanoma",
    "LIHC": "Hepatocellular Carcinoma",
    "LUAD": "Lung Adenocarcinoma",
    "LUSC": "Lung Squamous Cell Carcinoma",
    "NSCLC": "Non-Small Cell Lung",
    "SCLC": "Small Cell Lung Cancer",
    "LC": "Lung Cancer",
    "PDAC": "Pancreatic Ductal Adenocarcinoma",
    "PAAD": "Pancreatic Adenocarcinoma",
    "OV": "Ovarian Cancer",
    "THCA": "Thyroid Carcinoma",
    "ESCA": "Esophageal Carcinoma",
    "GC": "Gastric Cancer",
    "STAD": "Gastric Adenocarcinoma",
    "NPC": "Nasopharyngeal Carcinoma",
    "HNSCC": "Head and Neck Squamous Cell Carcinoma",
    "HNSC": "Head and Neck Squamous Cell Carcinoma",
    "RCC": "Renal Cell Carcinoma",
    "ccRCC": "Clear Cell Renal Cell Carcinoma",
    "KIRC": "Clear Cell Renal Cell Carcinoma",
    "GBM": "Glioblastoma",
    "BLCA": "Bladder Urothelial Carcinoma",
    "PRAD": "Prostate Adenocarcinoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma",
    "CHOL": "Cholangiocarcinoma",
    "ICC": "Intrahepatic Cholangiocarcinoma",
    "UVM": "Uveal Melanoma",
}

_CANCER_CODE_TO_ORGAN: Final[Mapping[str, str]] = {
    "BC": "Breast",
    "CRC": "Colon / Rectum",
    "MEL": "Skin",
    "LIHC": "Liver",
    "LUAD": "Lung",
    "LUSC": "Lung",
    "NSCLC": "Lung",
    "SCLC": "Lung",
    "LC": "Lung",
    "PDAC": "Pancreas",
    "PAAD": "Pancreas",
    "OV": "Ovary",
    "THCA": "Thyroid",
    "ESCA": "Esophagus",
    "GC": "Stomach",
    "STAD": "Stomach",
    "NPC": "Nasopharynx",
    "HNSCC": "Head and Neck",
    "HNSC": "Head and Neck",
    "RCC": "Kidney",
    "ccRCC": "Kidney",
    "KIRC": "Kidney",
    "GBM": "Brain",
    "BLCA": "Bladder",
    "PRAD": "Prostate",
    "UCEC": "Uterus",
    "CHOL": "Liver / Bile Duct",
    "ICC": "Liver / Bile Duct",
    "UVM": "Eye",
}

TamAssayChoice = Literal["RNA", "SCT", "integrated"]


class _AlignedProjectionSpec(BaseModel):
    """Immutable column-alignment operator mapping raw assay features onto canonical `.var`."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    keep_mask: np.ndarray[Any, Any]
    projection_matrix: sp.csr_matrix
    covered_var_mask: np.ndarray[Any, Any]


# -----------------------------------------------------------------------------
# Pure Functional Core: Polars Metadata & Feature Alignment
# -----------------------------------------------------------------------------


def _normalize_cluster_token(raw_token: str) -> str:
    """Normalize a raw Seurat cluster identifier (e.g., `'9'`, `'9.0'`, `'9_AngioMac'`) to `'0'`..`'23'`."""
    cleaned = raw_token.strip()
    if not cleaned or cleaned in {"NA", "nan", "None"}:
        return "NA"
    prefix = cleaned.split("_", 1)[0]
    if prefix.endswith(".0") and prefix[:-2].isdigit():
        return prefix[:-2]
    return prefix if prefix in COULTON_TAM_CLUSTER_LABELS else cleaned


def _cluster_id_to_short_label(cluster_id: str) -> str:
    """Map normalized cluster ID to Coulton et al. 2024 `short.label`."""
    if cluster_id in COULTON_TAM_CLUSTER_LABELS:
        return COULTON_TAM_CLUSTER_LABELS[cluster_id]
    return cluster_id if "_" in cluster_id else f"{cluster_id}_Unmapped"


def _short_label_to_cell_type(short_label: str) -> str:
    """Map Coulton `short.label` to canonical broad lineage (`Macrophage`, `Monocyte`, or doublet)."""
    match short_label:
        case "19_ClassMono":
            return "Monocyte"
        case "20_TDoub":
            return "T/Myeloid Doublet"
        case _:
            return "Macrophage"


def _normalize_cancer_code(raw_cancer: str) -> str:
    """Apply Coulton `mac.util.R` cancer-code harmonization."""
    cleaned = raw_cancer.strip()
    if not cleaned or cleaned in {"NA", "nan", "None"}:
        return "Pan-Cancer"
    return _CANCER_CODE_NORMALIZATION.get(cleaned, cleaned)


def _normalize_primary_met(raw_val: str) -> str:
    """Harmonize `primary_met2` / `primary_met` annotations per `mac.util.R`."""
    low = raw_val.strip().lower()
    match low:
        case "pri" | "primary" | "primary/laparoscopy" | "primary_tumor":
            return "primary"
        case "met" | "metastasis" | "metastatic":
            return "metastasis"
        case "norm" | "normal" | "adjacent" | "adjacent_normal":
            return "normal"
        case _:
            return "unknown"


def harmonize_tam_atlas_obs(
    raw_obs: pd.DataFrame,
    active_ident: Sequence[str] = (),
    cell_names_fallback: Sequence[str] = (),
) -> pd.DataFrame:
    """Harmonize Seurat `@meta.data` into canonical `tme_datasets` `.obs` using pure Polars transforms.

    Preserves all original `@meta.data` columns and enriches canonical columns:
    `cell_id`, `dataset`, `dataset_id`, `source_study`, `patient`, `patient_id`,
    `sample`, `sample_id`, `organ`, `cancer_code`, `cancer_type`, `cell_type`,
    `cell_state`, `cell_type_author`, `tam_cluster_id`, `short.label`,
    `primary_met_harmonized`, and `treatment_status`.
    """
    n_rows = len(raw_obs) if len(raw_obs) > 0 else len(cell_names_fallback)
    if n_rows == 0:
        return pd.DataFrame(index=pd.Index([], name="cell_id"))

    barcodes = (
        [str(idx) for idx in raw_obs.index]
        if len(raw_obs) == n_rows and len(raw_obs.index) == n_rows
        else [str(c) for c in cell_names_fallback[:n_rows]]
    )

    records: dict[str, Any] = {"cell_id": barcodes}
    for col in raw_obs.columns:
        col_str = str(col).strip()
        if not col_str or col_str == "cell_id":
            continue
        series = raw_obs[col]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            records[col_str] = series.to_numpy()
        elif pd.api.types.is_bool_dtype(series):
            records[col_str] = series.to_numpy(dtype=bool)
        else:
            records[col_str] = [
                "NA" if pd.isna(v) or str(v).strip() == "" else str(v)
                for v in series.tolist()
            ]

    if len(active_ident) == n_rows:
        records["_active_ident"] = [str(v) for v in active_ident]
    else:
        records["_active_ident"] = ["NA"] * n_rows

    pl_df = pl.DataFrame(records)

    # Resolve cluster source column in priority order
    cluster_candidates = (
        "short.label",
        "integrated_snn_res.0.525",
        "seurat_clusters",
        "cluster",
        "_active_ident",
    )
    available_cluster_cols = [c for c in cluster_candidates if c in pl_df.columns]
    raw_cluster_expr = (
        pl.coalesce(
            [
                pl.when(
                    pl.col(c).cast(pl.String).is_not_null()
                    & ~pl.col(c).cast(pl.String).is_in(["", "NA", "nan", "None"])
                )
                .then(pl.col(c).cast(pl.String))
                .otherwise(None)
                for c in available_cluster_cols
            ]
        ).fill_null("NA")
        if available_cluster_cols
        else pl.lit("NA")
    )

    study_expr = (
        pl.col("study").cast(pl.String).fill_null("Unknown")
        if "study" in pl_df.columns
        else pl.lit("Coulton2024")
    )
    patient_raw_expr = (
        pl.col("patient").cast(pl.String).fill_null("NA")
        if "patient" in pl_df.columns
        else pl.lit("NA")
    )
    sample_raw_expr = (
        pl.col("sample").cast(pl.String).fill_null("NA")
        if "sample" in pl_df.columns
        else pl.lit("NA")
    )
    cancer_raw_expr = (
        pl.col("cancer").cast(pl.String).fill_null("Pan-Cancer")
        if "cancer" in pl_df.columns
        else pl.lit("Pan-Cancer")
    )
    tissue_raw_expr = (
        pl.col("tissue").cast(pl.String).fill_null("Tumor")
        if "tissue" in pl_df.columns
        else pl.lit("Tumor")
    )
    prim_met_raw_expr = (
        pl.col("primary_met2").cast(pl.String).fill_null("NA")
        if "primary_met2" in pl_df.columns
        else (
            pl.col("primary_met").cast(pl.String).fill_null("NA")
            if "primary_met" in pl_df.columns
            else pl.lit("NA")
        )
    )
    treatment_raw_expr = (
        pl.col("pre.or.post.treatment").cast(pl.String).fill_null("NA")
        if "pre.or.post.treatment" in pl_df.columns
        else pl.lit("NA")
    )
    author_ct_raw_expr = (
        pl.col("celltype").cast(pl.String).fill_null("NA")
        if "celltype" in pl_df.columns
        else pl.lit("NA")
    )

    enriched = (
        pl_df.lazy()
        .with_columns(
            raw_cluster_expr.map_elements(_normalize_cluster_token, return_dtype=pl.String).alias("tam_cluster_id"),
            cancer_raw_expr.map_elements(_normalize_cancer_code, return_dtype=pl.String).alias("cancer_code"),
            study_expr.alias("source_study"),
            pl.when(patient_raw_expr.is_in(["", "NA", "nan", "None"]))
            .then(study_expr)
            .otherwise(patient_raw_expr)
            .alias("_clean_patient"),
            pl.when(tissue_raw_expr.is_in(["", "NA", "nan", "None"]))
            .then(pl.lit("Tumor"))
            .otherwise(tissue_raw_expr)
            .alias("tissue"),
            prim_met_raw_expr.map_elements(_normalize_primary_met, return_dtype=pl.String).alias(
                "primary_met_harmonized"
            ),
            pl.when(treatment_raw_expr.is_in(["", "NA", "nan", "None"]))
            .then(pl.lit("unknown"))
            .otherwise(treatment_raw_expr)
            .alias("treatment_status"),
            author_ct_raw_expr.alias("_author_ct_raw"),
        )
        .with_columns(
            pl.col("tam_cluster_id")
            .map_elements(_cluster_id_to_short_label, return_dtype=pl.String)
            .alias("short.label"),
            pl.col("cancer_code")
            .replace_strict(_CANCER_CODE_TO_FULL_NAME, default=pl.col("cancer_code"))
            .alias("cancer_type"),
            pl.col("cancer_code")
            .replace_strict(_CANCER_CODE_TO_ORGAN, default=pl.lit("Multi-organ"))
            .alias("organ"),
            pl.col("_clean_patient").alias("patient"),
            pl.col("_clean_patient").alias("patient_id"),
            pl.when(sample_raw_expr.is_in(["", "NA", "nan", "None"]))
            .then(pl.col("_clean_patient"))
            .otherwise(sample_raw_expr)
            .alias("_clean_sample"),
        )
        .with_columns(
            pl.col("_clean_sample").alias("sample"),
            pl.col("_clean_sample").alias("sample_id"),
            pl.col("short.label")
            .map_elements(_short_label_to_cell_type, return_dtype=pl.String)
            .alias("cell_type"),
            pl.col("short.label").alias("cell_state"),
            pl.when(pl.col("_author_ct_raw").is_in(["", "NA", "nan", "None"]))
            .then(pl.col("short.label"))
            .otherwise(pl.col("_author_ct_raw"))
            .alias("cell_type_author"),
            pl.lit(COULTON_TAM_DATASET_ID).alias("dataset"),
            pl.lit(COULTON_TAM_DATASET_ID).alias("dataset_id"),
        )
        .drop(["_active_ident", "_clean_patient", "_clean_sample", "_author_ct_raw"])
        .collect()
    )

    out_pdf = enriched.to_pandas()
    out_pdf = out_pdf.set_index("cell_id", drop=False)
    out_pdf.index.name = None
    return H5ADSparseIncrementalWriter._sanitize_dataframe(out_pdf, default_index_name="cell_id")


def _build_secondary_assay_projection(
    secondary_genes: Sequence[str],
    target_var_ids: Sequence[str],
    raw_symbol_to_target_id: Mapping[str, str],
) -> _AlignedProjectionSpec:
    """Build a sparse CSR column-projection operator aligning `secondary_genes` onto `target_var_ids`."""
    n_target = len(target_var_ids)
    target_id_to_col = {gid: idx for idx, gid in enumerate(target_var_ids)}

    mapped_cols: list[int] = []
    keep_flags: list[bool] = []
    covered_mask = np.zeros(n_target, dtype=bool)

    for raw_g in secondary_genes:
        target_id = raw_symbol_to_target_id.get(raw_g, raw_g)
        if target_id in target_id_to_col:
            c_idx = target_id_to_col[target_id]
            keep_flags.append(True)
            mapped_cols.append(c_idx)
            covered_mask[c_idx] = True
        else:
            keep_flags.append(False)

    keep_mask = np.asarray(keep_flags, dtype=bool)
    n_kept = int(keep_mask.sum())
    if n_kept == 0:
        empty_m = sp.csr_matrix(
            (np.empty(0, dtype=np.float32), (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32))),
            shape=(0, n_target),
            dtype=np.float32,
        )
        return _AlignedProjectionSpec(
            keep_mask=keep_mask,
            projection_matrix=empty_m,
            covered_var_mask=covered_mask,
        )

    rows = np.arange(n_kept, dtype=np.int32)
    cols = np.asarray(mapped_cols, dtype=np.int32)
    weights = np.ones(n_kept, dtype=np.float32)
    proj_m = sp.csr_matrix((weights, (rows, cols)), shape=(n_kept, n_target), dtype=np.float32)
    return _AlignedProjectionSpec(
        keep_mask=keep_mask,
        projection_matrix=proj_m,
        covered_var_mask=covered_mask,
    )


# -----------------------------------------------------------------------------
# Imperative Shell: Streaming HDF5 Helpers, Downloader, Preprocessor & Loader
# -----------------------------------------------------------------------------


def _slice_staged_csr_batch(
    assay_grp: h5py.Group,
    start_row: int,
    end_row: int,
    n_cols: int,
) -> sp.csr_matrix:
    """Read a `[start_row:end_row]` CSR batch directly from staging HDF5 in O(batch_nnz) RAM."""
    indptr_slice = np.asarray(assay_grp["indptr"][start_row : end_row + 1], dtype=np.int64)
    p0 = int(indptr_slice[0])
    p1 = int(indptr_slice[-1])
    if p1 <= p0:
        return sp.csr_matrix(
            (np.empty(0, dtype=np.float32), (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32))),
            shape=(end_row - start_row, n_cols),
            dtype=np.float32,
        )

    data_slice = np.asarray(assay_grp["data"][p0:p1], dtype=np.float32)
    indices_slice = np.asarray(assay_grp["indices"][p0:p1], dtype=np.int32)
    local_indptr = indptr_slice - p0
    return sp.csr_matrix(
        (data_slice, indices_slice, local_indptr),
        shape=(end_row - start_row, n_cols),
        dtype=np.float32,
    )


def _init_extendable_h5_csr_group(
    parent: h5py.File | h5py.Group,
    name: str,
    n_obs: int,
    n_vars: int,
    chunk_nnz: int = 262_144,
    chunk_indptr: int = 32_768,
) -> h5py.Group:
    """Create a gzip-compressed extendable AnnData CSR matrix group (`X` or `layers/<name>`)."""
    grp = parent.create_group(name)
    grp.attrs["encoding-type"] = "csr_matrix"
    grp.attrs["encoding-version"] = "0.1.0"
    grp.attrs["shape"] = np.array([n_obs, n_vars], dtype=np.int64)
    grp.create_dataset(
        "data",
        shape=(0,),
        maxshape=(None,),
        dtype=np.float32,
        chunks=(chunk_nnz,),
        compression="gzip",
        compression_opts=4,
    )
    grp.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(None,),
        dtype=np.int32,
        chunks=(chunk_nnz,),
        compression="gzip",
        compression_opts=4,
    )
    grp.create_dataset(
        "indptr",
        shape=(1,),
        maxshape=(None,),
        dtype=np.int64,
        chunks=(min(max(1, n_obs + 1), chunk_indptr),),
        data=np.array([0], dtype=np.int64),
    )
    return grp


def _append_csr_to_h5_group(grp: h5py.Group, batch_csr: sp.csr_matrix) -> None:
    """Append a projected CSR cell batch to an extendable HDF5 CSR group."""
    ds_data = grp["data"]
    ds_indices = grp["indices"]
    ds_indptr = grp["indptr"]

    cur_nnz = int(ds_data.shape[0])
    cur_indptr_len = int(ds_indptr.shape[0])
    batch_nnz = int(batch_csr.nnz)
    n_rows = int(batch_csr.shape[0])

    if batch_nnz > 0:
        ds_data.resize((cur_nnz + batch_nnz,))
        ds_data[cur_nnz:] = batch_csr.data.astype(np.float32, copy=False)
        ds_indices.resize((cur_nnz + batch_nnz,))
        ds_indices[cur_nnz:] = batch_csr.indices.astype(np.int32, copy=False)

    ds_indptr.resize((cur_indptr_len + n_rows,))
    ds_indptr[cur_indptr_len:] = batch_csr.indptr[1:].astype(np.int64) + cur_nnz


def download_coulton_tam_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
) -> Result[Path, str]:
    """Download `mac.atlas.zenodo.200524.rds` (19.75 GB) from Zenodo record 11222158 with HTTP 206 resume.

    Args:
        raw_dir: Optional raw target directory wrapped in `Maybe` (defaults to `data/raw/Coulton_TAM_Atlas`).
        force_download: Re-download from scratch even if a verified local copy exists.
        verify_md5: Verify file byte length (`19,748,833,115`) and MD5 (`49d3bea7addd7fcd07d8e3843ca901c2`).

    Returns:
        `Success(rds_path)` or `Failure(error_message)`.
    """
    target_dir = (
        raw_dir.unwrap()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(COULTON_TAM_DATASET_ID)
    )
    target_dir.mkdir(parents=True, exist_ok=True)
    dest_rds = target_dir / COULTON_TAM_RDS_FILENAME
    marker_file = target_dir / f".{COULTON_TAM_RDS_FILENAME}.md5.ok"

    if force_download:
        dest_rds.unlink(missing_ok=True)
        marker_file.unlink(missing_ok=True)
        dest_rds.with_name(dest_rds.name + ".part").unlink(missing_ok=True)

    # Fast sidecar check to avoid re-hashing 19.75 GB if size & mtime_ns match previous verification
    if verify_md5 and dest_rds.is_file() and marker_file.is_file():
        st = dest_rds.stat()
        expected_stamp = f"{st.st_size}:{st.st_mtime_ns}:{COULTON_TAM_EXPECTED_MD5}"
        try:
            if (
                st.st_size == COULTON_TAM_EXPECTED_SIZE_BYTES
                and marker_file.read_text(encoding="utf-8").strip() == expected_stamp
            ):
                logger.info(
                    "Verified cached %s via MD5 sidecar (%d bytes)",
                    dest_rds.name,
                    st.st_size,
                )
                return Success(dest_rds)
        except OSError:
            pass

    res = download_single_file(
        url=COULTON_TAM_ZENODO_URL,
        dest_path=dest_rds,
        timeout_sec=180.0,
        max_retries=10,
        resume=True,
        checksum=Some(COULTON_TAM_CHECKSUM) if verify_md5 else Nothing,
    )
    match res:
        case Success(verified_path):
            if verify_md5 and verified_path.is_file():
                try:
                    st = verified_path.stat()
                    marker_file.write_text(
                        f"{st.st_size}:{st.st_mtime_ns}:{COULTON_TAM_EXPECTED_MD5}\n",
                        encoding="utf-8",
                    )
                except OSError:
                    pass
            return Success(verified_path)
        case Failure(err):
            return Failure(err)
        case _:
            return Failure("Unexpected download result state")


def preprocess_coulton_tam_atlas(
    rds_path_or_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 25_000,
) -> Result[Path, str]:
    """Stream `mac.atlas.zenodo.200524.rds` into a multi-assay Ensembl-normalized `.h5ad` (< 1.5 GB RAM).

    Writes:
    - `.X` and `.layers["log1p_norm"]`: `RNA@data` sparse CSR (`float32`, `log1p_sum` Ensembl aggregation)
    - `.layers["sct_data"]`: `SCT@data` sparse CSR (`float32`, aligned to `.var_names`)
    - `.obsm["X_rpca_integrated"]`: `integrated@data` (`float32` 2D gzip-compressed array)
    - `.obsm["X_pca"]`, `.obsm["X_umap"]`: Seurat dimensional reductions (`float32`)
    - `.obs`: Polars-harmonized metadata + all original Seurat `@meta.data` columns
    - `.var`: Ensembl Release 111 metadata + `in_sct_assay`, `in_rna_hvg`, `in_sct_hvg`, `is_rpca_anchor`
    - `.uns`: Provenance, `is_raw_counts=False`, `expression_type="log_normalized"`, `rpca_integrated_var_names`
    """
    rds_path = (
        rds_path_or_dir / COULTON_TAM_RDS_FILENAME
        if rds_path_or_dir.is_dir()
        else rds_path_or_dir
    )
    if not rds_path.is_file():
        return Failure(f"Coulton TAM Atlas RDS file not found at {rds_path}")

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(COULTON_TAM_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    scratch_dir = get_scratch_dataset_dir(COULTON_TAM_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    staging_h5 = scratch_dir / "coulton_tam_staging.h5"

    try:
        # 1. Stream R XDR v3 Seurat S4 `.rds` -> staging HDF5
        stage_res = stream_seurat_rds_to_staging_h5(rds_path=rds_path, staging_h5_path=staging_h5)
        match stage_res:
            case Failure(err):
                return Failure(err)
            case Success(manifest):
                pass
            case _:
                return Failure("Unexpected staging result state")

        if "RNA" not in manifest.assays:
            return Failure(
                f"Expected 'RNA' assay in {rds_path.name}, found assays: {sorted(manifest.assays.keys())}"
            )

        rna_spec: StagedAssaySpec = manifest.assays["RNA"]
        n_cells = rna_spec.n_cells
        raw_rna_genes = list(rna_spec.gene_names)
        if not raw_rna_genes:
            raw_rna_genes = [f"GENE_{i}" for i in range(rna_spec.n_genes)]

        # 2. Harmonize `.obs` with Polars
        obs_clean = harmonize_tam_atlas_obs(
            raw_obs=manifest.obs_df,
            active_ident=manifest.active_ident,
            cell_names_fallback=rna_spec.cell_names,
        )
        if len(obs_clean) != n_cells:
            return Failure(
                f"Cell count mismatch between RNA@data ({n_cells}) and @meta.data ({len(obs_clean)})"
            )

        # 3. Build `.var` and sparse column projection operators for RNA and SCT
        raw_rna_var = pd.DataFrame(index=pd.Index(raw_rna_genes, name="original_id"))
        unmapped_genes: list[str] = []
        raw_symbol_to_target_id: dict[str, str] = {}

        if normalize_ensembl:
            ensembl_dir = get_ensembl_dir()
            ensembl = ensure_ensembl_release_installed(
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
            )
            cache_path = Path(ensembl_dir) / f"gene_mapping_cache_release_{ensembl_release}.parquet"

            # Pre-warm mapping cache across RNA + SCT gene unions in one pass
            all_query_genes = list(dict.fromkeys([
                *raw_rna_genes,
                *(manifest.assays["SCT"].gene_names if "SCT" in manifest.assays else ()),
            ]))
            mapping_df = _load_or_update_mapping_cache(all_query_genes, ensembl, cache_path)
            for q_sym, gid_val, m_status in zip(
                mapping_df["query"].tolist(),
                mapping_df["gene_id"].tolist(),
                mapping_df["mapping_status"].tolist(),
                strict=False,
            ):
                if gid_val and str(m_status) != "unmapped":
                    raw_symbol_to_target_id[str(q_sym)] = str(gid_val)
                elif not drop_unmapped:
                    raw_symbol_to_target_id[str(q_sym)] = str(q_sym)

            var_final, rna_keep_mask, rna_M, unmapped_genes = _build_ensembl_projection_operator(
                raw_var_names=raw_rna_genes,
                raw_var=raw_rna_var,
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
                drop_unmapped=drop_unmapped,
                aggregation="sum",
            )
        else:
            unique_idx = ad.utils.make_index_unique(pd.Index(raw_rna_genes, name="gene_id"))
            var_final = pd.DataFrame(
                {"gene_name": raw_rna_genes, "original_id": raw_rna_genes},
                index=unique_idx,
            )
            rna_keep_mask = np.ones(len(raw_rna_genes), dtype=bool)
            rna_M = sp.eye(len(raw_rna_genes), dtype=np.float32, format="csr")
            raw_symbol_to_target_id = {
                raw_g: str(uid) for raw_g, uid in zip(raw_rna_genes, unique_idx.tolist(), strict=False)
            }

        target_var_ids = [str(v) for v in var_final.index.tolist()]
        n_final_vars = len(target_var_ids)

        sct_proj_maybe: Maybe[_AlignedProjectionSpec] = Nothing
        if "SCT" in manifest.assays:
            sct_spec = manifest.assays["SCT"]
            sct_genes = (
                list(sct_spec.gene_names)
                if sct_spec.gene_names
                else [f"SCT_GENE_{i}" for i in range(sct_spec.n_genes)]
            )
            sct_proj = _build_secondary_assay_projection(
                secondary_genes=sct_genes,
                target_var_ids=target_var_ids,
                raw_symbol_to_target_id=raw_symbol_to_target_id,
            )
            sct_proj_maybe = Some(sct_proj)
            var_final["in_sct_assay"] = sct_proj.covered_var_mask
        else:
            var_final["in_sct_assay"] = np.zeros(n_final_vars, dtype=bool)

        rna_hvg_target_ids = {
            raw_symbol_to_target_id.get(g, g) for g in rna_spec.var_features
        }
        sct_hvg_target_ids = (
            {raw_symbol_to_target_id.get(g, g) for g in manifest.assays["SCT"].var_features}
            if "SCT" in manifest.assays
            else set()
        )
        rpca_genes = (
            list(manifest.assays["integrated"].gene_names)
            if "integrated" in manifest.assays
            else []
        )
        rpca_target_ids = {raw_symbol_to_target_id.get(g, g) for g in rpca_genes}

        var_final["in_rna_hvg"] = np.asarray([vid in rna_hvg_target_ids for vid in target_var_ids], dtype=bool)
        var_final["in_sct_hvg"] = np.asarray([vid in sct_hvg_target_ids for vid in target_var_ids], dtype=bool)
        var_final["is_rpca_anchor"] = np.asarray([vid in rpca_target_ids for vid in target_var_ids], dtype=bool)
        clean_var = H5ADSparseIncrementalWriter._sanitize_dataframe(var_final, default_index_name="gene_id")
        if clean_var.index.name in clean_var.columns:
            clean_var.index.name = None

        if obs_clean.index.name in obs_clean.columns:
            obs_clean.index.name = None

        # 4. Stream batches from staging HDF5 -> final `.h5ad.tmp`
        tmp_h5ad.unlink(missing_ok=True)
        rna_agg = "log1p_sum" if normalize_ensembl else "sum"

        with ExitStack() as stack:
            src_h5 = stack.enter_context(h5py.File(staging_h5, "r"))
            dst_h5 = stack.enter_context(h5py.File(tmp_h5ad, "w"))

            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"

            anndata.io.write_elem(dst_h5, "obs", obs_clean)
            anndata.io.write_elem(dst_h5, "var", clean_var)

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            # Pass 4a: Stream `RNA@data` -> `/X` and `/layers/log1p_norm`
            logger.info(
                "Streaming RNA@data (%d cells x %d -> %d genes) to /X and /layers/log1p_norm...",
                n_cells,
                rna_spec.n_genes,
                n_final_vars,
            )
            grp_x = _init_extendable_h5_csr_group(dst_h5, "X", n_obs=n_cells, n_vars=n_final_vars)
            grp_log1p = _init_extendable_h5_csr_group(
                layers_grp,
                "log1p_norm",
                n_obs=n_cells,
                n_vars=n_final_vars,
            )
            src_rna_grp = src_h5[rna_spec.h5_group]

            for start in range(0, n_cells, batch_size):
                end = min(start + batch_size, n_cells)
                if rna_spec.is_sparse:
                    raw_chunk = _slice_staged_csr_batch(src_rna_grp, start, end, rna_spec.n_genes)
                else:
                    raw_chunk = sp.csr_matrix(src_rna_grp["dense"][start:end, :], dtype=np.float32)
                proj_chunk = _project_sparse_chunk_to_ensembl(
                    raw_chunk,
                    rna_keep_mask,
                    rna_M,
                    rna_agg,
                )
                _append_csr_to_h5_group(grp_x, proj_chunk)
                _append_csr_to_h5_group(grp_log1p, proj_chunk)
                gc.collect()

            # Pass 4b: Stream `SCT@data` -> `/layers/sct_data`
            match sct_proj_maybe:
                case Some(sct_proj):
                    sct_spec_cur = manifest.assays["SCT"]
                    logger.info(
                        "Streaming SCT@data (%d cells x %d -> %d genes) to /layers/sct_data...",
                        n_cells,
                        sct_spec_cur.n_genes,
                        n_final_vars,
                    )
                    grp_sct = _init_extendable_h5_csr_group(
                        layers_grp,
                        "sct_data",
                        n_obs=n_cells,
                        n_vars=n_final_vars,
                    )
                    src_sct_grp = src_h5[sct_spec_cur.h5_group]
                    for start in range(0, n_cells, batch_size):
                        end = min(start + batch_size, n_cells)
                        if sct_spec_cur.is_sparse:
                            raw_sct = _slice_staged_csr_batch(
                                src_sct_grp,
                                start,
                                end,
                                sct_spec_cur.n_genes,
                            )
                        else:
                            raw_sct = sp.csr_matrix(
                                src_sct_grp["dense"][start:end, :],
                                dtype=np.float32,
                            )
                        proj_sct = _project_sparse_chunk_to_ensembl(
                            raw_sct,
                            sct_proj.keep_mask,
                            sct_proj.projection_matrix,
                            rna_agg,
                        )
                        _append_csr_to_h5_group(grp_sct, proj_sct)
                        gc.collect()
                case _:
                    pass

            # Pass 4c: Stream `integrated@data` & reductions -> `/obsm`
            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"

            if "integrated" in manifest.assays:
                int_spec = manifest.assays["integrated"]
                n_int_genes = int_spec.n_genes
                logger.info(
                    "Streaming integrated@data (%d cells x %d anchor genes) to /obsm/X_rpca_integrated...",
                    n_cells,
                    n_int_genes,
                )
                int_row_chunk = min(max(1, n_cells), 2048)
                int_col_chunk = min(max(1, n_int_genes), 2000)
                ds_rpca = obsm_grp.create_dataset(
                    "X_rpca_integrated",
                    shape=(n_cells, n_int_genes),
                    dtype=np.float32,
                    chunks=(int_row_chunk, int_col_chunk),
                    compression="gzip",
                    compression_opts=4,
                )
                ds_rpca.attrs["encoding-type"] = "array"
                ds_rpca.attrs["encoding-version"] = "0.2.0"

                src_int_grp = src_h5[int_spec.h5_group]
                dense_step = 8_192
                for start in range(0, n_cells, dense_step):
                    end = min(start + dense_step, n_cells)
                    if int_spec.is_sparse:
                        csr_int = _slice_staged_csr_batch(src_int_grp, start, end, n_int_genes)
                        ds_rpca[start:end, :] = csr_int.toarray()
                    else:
                        ds_rpca[start:end, :] = np.asarray(
                            src_int_grp["dense"][start:end, :],
                            dtype=np.float32,
                        )
                    gc.collect()

            for red_key, red_spec in manifest.reductions.items():
                if red_spec.n_cells == n_cells:
                    obsm_name = f"X_{red_key.lower()}" if not red_key.startswith("X_") else red_key
                    emb_arr = np.asarray(src_h5[red_spec.h5_dataset][:], dtype=np.float32)
                    anndata.io.write_elem(obsm_grp, obsm_name, emb_arr)

            # Pass 4d: Write `/uns` provenance & expression metadata
            uns_payload: dict[str, Any] = {
                "dataset_id": COULTON_TAM_DATASET_ID,
                "title": "Coulton et al. 2024 Pan-Cancer Tumor-Associated Macrophage (TAM) Atlas",
                "publication_doi": COULTON_TAM_PUBLICATION_DOI,
                "publication_pmid": COULTON_TAM_PUBLICATION_PMID,
                "zenodo_record": COULTON_TAM_ZENODO_RECORD,
                "zenodo_doi": COULTON_TAM_ZENODO_DOI,
                "source_file": COULTON_TAM_RDS_FILENAME,
                "is_raw_counts": False,
                "has_raw_counts": False,
                "raw_counts_note": (
                    "Zenodo record 11222158 omits raw UMI counts (obtainable from original source studies) "
                    "and provides log-normalized RNA@data (.X / .layers['log1p_norm']), "
                    "SCTransformed SCT@data (.layers['sct_data']), and RPCA integrated@data "
                    "(.obsm['X_rpca_integrated'])."
                ),
                "expression_type": "log_normalized",
                "assays_available": sorted(manifest.assays.keys()),
                "seurat_active_assay": manifest.active_assay,
                "seurat_project_name": manifest.project_name,
                "tam_cluster_labels": dict(COULTON_TAM_CLUSTER_LABELS),
                "rpca_integrated_var_names": rpca_genes,
                "rna_var_features": list(rna_spec.var_features),
                "sct_var_features": (
                    list(manifest.assays["SCT"].var_features) if "SCT" in manifest.assays else []
                ),
                "n_unmapped_genes": len(unmapped_genes),
            }
            if unmapped_genes:
                uns_payload["unmapped_genes"] = unmapped_genes

            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        # Atomic commit
        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        size_mb = target_h5ad.stat().st_size / (1024.0 * 1024.0)
        logger.info(
            "Successfully wrote Coulton TAM Atlas H5AD to %s (%.1f MB, %d cells x %d genes)",
            target_h5ad,
            size_mb,
            n_cells,
            n_final_vars,
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        msg = f"Failed to preprocess Coulton TAM Atlas from {rds_path}: {exc}"
        logger.error(msg)
        return Failure(msg)
    finally:
        staging_h5.unlink(missing_ok=True)


def _read_coulton_h5ad_by_assay(
    h5ad_path: Path,
    assay: TamAssayChoice = "RNA",
    include_all_layers: bool = False,
) -> ad.AnnData:
    """Load `Coulton_TAM_Atlas.h5ad` with low-RAM selective assay reading or full multi-layer reading."""
    if assay == "RNA" and include_all_layers:
        return ad.read_h5ad(h5ad_path)

    with h5py.File(h5ad_path, "r") as f:
        obs = anndata.io.read_elem(f["obs"])
        uns = anndata.io.read_elem(f["uns"]) if "uns" in f else {}
        obsm_dict: dict[str, Any] = {}
        if "obsm" in f:
            for k in f["obsm"].keys():
                if include_all_layers or k != "X_rpca_integrated":
                    obsm_dict[k] = anndata.io.read_elem(f["obsm"][k])

        match assay:
            case "RNA":
                var = anndata.io.read_elem(f["var"])
                X_mat = anndata.io.read_elem(f["X"])
                return ad.AnnData(X=X_mat, obs=obs, var=var, uns=uns, obsm=obsm_dict)

            case "SCT":
                var = anndata.io.read_elem(f["var"])
                if "layers" not in f or "sct_data" not in f["layers"]:
                    raise KeyError(f"'layers/sct_data' not present in {h5ad_path}")
                sct_mat = anndata.io.read_elem(f["layers"]["sct_data"])
                layers_out = (
                    {lk: anndata.io.read_elem(f["layers"][lk]) for lk in f["layers"].keys()}
                    if include_all_layers
                    else {"sct_data": sct_mat}
                )
                uns_sct = dict(uns)
                uns_sct["active_loaded_assay"] = "SCT"
                return ad.AnnData(
                    X=sct_mat,
                    obs=obs,
                    var=var,
                    layers=layers_out,
                    uns=uns_sct,
                    obsm=obsm_dict,
                )

            case "integrated":
                if "obsm" not in f or "X_rpca_integrated" not in f["obsm"]:
                    raise KeyError(f"'obsm/X_rpca_integrated' not present in {h5ad_path}")
                rpca_mat = anndata.io.read_elem(f["obsm"]["X_rpca_integrated"])
                rpca_genes = [str(g) for g in uns.get("rpca_integrated_var_names", [])]
                if len(rpca_genes) != rpca_mat.shape[1]:
                    rpca_genes = [f"ANCHOR_{i}" for i in range(rpca_mat.shape[1])]
                rpca_var = pd.DataFrame(
                    {"gene_name": rpca_genes},
                    index=ad.utils.make_index_unique(pd.Index(rpca_genes, name="gene_symbol")),
                )
                uns_int = dict(uns)
                uns_int["active_loaded_assay"] = "integrated"
                uns_int["expression_type"] = "rpca_integrated"
                return ad.AnnData(
                    X=rpca_mat,
                    obs=obs,
                    var=rpca_var,
                    uns=uns_int,
                    obsm=obsm_dict,
                )

            case _ as unreachable:
                assert_never(unreachable)


def load_coulton_tam_atlas(
    raw_dir: Maybe[Path] = Nothing,
    output_h5ad: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    force_recompute: bool = False,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    include_all_layers: bool = False,
    assay: TamAssayChoice = "RNA",
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Load the Coulton et al. (2024) Pan-Cancer TAM Atlas (`Coulton_TAM_Atlas`).

    Prioritizes loading from `data/preprocessed/Coulton_TAM_Atlas.h5ad`. On cold start,
    streams `mac.atlas.zenodo.200524.rds` (auto-downloading with HTTP 206 resume if enabled)
    directly to `Coulton_TAM_Atlas.h5ad` in < 1.5 GB RAM.

    Args:
        raw_dir: Optional directory containing `mac.atlas.zenodo.200524.rds`.
        output_h5ad: Optional explicit `.h5ad` output path.
        auto_download: Automatically download `mac.atlas.zenodo.200524.rds` from Zenodo if absent.
        force_download: Re-download the `.rds` file from Zenodo.
        force_recompute: Rebuild the `.h5ad` file from the raw `.rds` file.
        normalize_ensembl: Map gene symbols to canonical Ensembl Release 111 IDs (`log1p_sum`).
        ensembl_release: Target Ensembl release (default 111).
        include_all_layers: Load `.layers['log1p_norm']`, `.layers['sct_data']`, and
            `.obsm['X_rpca_integrated']` simultaneously into memory alongside `.X`.
        assay: Which assay to expose in `.X` (`'RNA'`, `'SCT'`, or `'integrated'`).
        subset: Optional column-value mapping to subset observations.
        subsample_n: Optional cell count to subsample.
        repo_root: Optional repository root override.
    """
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(COULTON_TAM_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(COULTON_TAM_DATASET_ID, repo_root=root)
    )

    h5ad_to_read: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        if target_h5ad.is_file() and target_h5ad.stat().st_size > 0:
            h5ad_to_read = Some(target_h5ad)
        else:
            h5ad_to_read = find_dataset_h5ad(COULTON_TAM_DATASET_ID, repo_root=root)

    match h5ad_to_read:
        case Some(existing_h5ad):
            try:
                adata = _read_coulton_h5ad_by_assay(
                    existing_h5ad,
                    assay=assay,
                    include_all_layers=include_all_layers,
                )
                if subset or subsample_n is not None:
                    adata = _apply_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                return Success(adata)
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read cached {existing_h5ad}: {exc}")
        case _:
            pass

    rds_candidate = resolved_raw_dir / COULTON_TAM_RDS_FILENAME
    if force_download or not rds_candidate.is_file():
        if not auto_download:
            return Failure(
                f"Coulton TAM Atlas raw RDS not found at {rds_candidate} and auto_download=False. "
                f"Download from {COULTON_TAM_ZENODO_URL} or set auto_download=True."
            )
        dl_res = download_coulton_tam_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
        )
        match dl_res:
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(downloaded_path):
                rds_candidate = downloaded_path
            case _:
                return Failure("Unexpected download result state")

    prep_res = preprocess_coulton_tam_atlas(
        rds_path_or_dir=rds_candidate,
        output_h5ad=Some(target_h5ad),
        normalize_ensembl=normalize_ensembl,
        ensembl_release=ensembl_release,
    )
    match prep_res:
        case Failure(prep_err):
            return Failure(prep_err)
        case Success(written_h5ad):
            try:
                adata = _read_coulton_h5ad_by_assay(
                    written_h5ad,
                    assay=assay,
                    include_all_layers=include_all_layers,
                )
                if subset or subsample_n is not None:
                    adata = _apply_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                return Success(adata)
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Preprocessed {written_h5ad} succeeded, but reading failed: {exc}")
        case _:
            return Failure("Unexpected preprocessing result state")


__all__ = [
    "COULTON_TAM_CHECKSUM",
    "COULTON_TAM_CLUSTER_LABELS",
    "COULTON_TAM_DATASET_ID",
    "COULTON_TAM_EXPECTED_MD5",
    "COULTON_TAM_EXPECTED_SIZE_BYTES",
    "COULTON_TAM_PUBLICATION_DOI",
    "COULTON_TAM_PUBLICATION_PMID",
    "COULTON_TAM_RDS_FILENAME",
    "COULTON_TAM_ZENODO_DOI",
    "COULTON_TAM_ZENODO_RECORD",
    "COULTON_TAM_ZENODO_URL",
    "TamAssayChoice",
    "download_coulton_tam_atlas",
    "harmonize_tam_atlas_obs",
    "load_coulton_tam_atlas",
    "preprocess_coulton_tam_atlas",
]
