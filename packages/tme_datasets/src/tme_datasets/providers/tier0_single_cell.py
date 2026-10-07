"""Pure functional loaders for Tier 0 Premier Benchmark Core single-cell RNA-seq cohorts.

Covers all 22 high-powered clinical immunotherapy trial cohorts across 9 human solid tumor indications:
- Breast: GSE246613, GSE300475, GSE212707
- CRC: GSE236581, GSE299651, CELLxGENE_829a3cd1
- Gastric: GSE270680
- HCC: GSE313642, GSE245906
- HNSCC: GSE301741, GSE287301, GSE200996
- Melanoma: CELLxGENE_7b20c613, GSE218429, GSE344166
- NSCLC: GSE207422, GSE317309, GSE243013
- PDAC: GSE311789, GSE316195
- ccRCC: GSE210038, GSE314072

All loaders:
- Return Result[ad.AnnData, str] (pure functional Railway Oriented Programming).
- Produce cells x genes CSR sparse matrix in adata.X with raw non-negative integer counts.
- Enforce unique gene identifiers via adata.var_names_make_unique().
- Standardize clinical observation metadata in adata.obs (patient_id, sample_id, treatment_status, clinical_response, indication, technology, cell_selection_strategy).
- Tag expression metadata via tag_expression_metadata().
"""

from __future__ import annotations

import gc
import gzip
from pathlib import Path
import re
import shutil
import tarfile
from typing import Callable, Mapping
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from returns.result import Failure, Result, Success
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

from ..logging import get_logger
from ..preprocessing.matrix_inspection import tag_expression_metadata

logger = get_logger("providers.tier0_single_cell")


# =========================================================================
# Shared Helper Functions & DRY Utilities
# =========================================================================
def _to_csr(mat: object) -> sp.csr_matrix:
    """Converts a sparse or dense matrix (genes x cells) into cells x genes CSR matrix."""
    if sp.issparse(mat):
        return mat.T.tocsr()  # type: ignore[attr-defined]
    return sp.csr_matrix(np.asarray(mat).T, dtype=np.float32)


def _harmonize_single_response(val: str | None) -> tuple[str, str]:
    """Purely maps a raw response label to canonical (clinical_response, clinical_response_raw)."""
    if val is None or pd.isna(val):
        return "not-evaluable", "NA"
    raw_str = str(val).strip()
    v_clean = raw_str.lower().replace("_", " ").replace("-", " ")

    # Check non-responders first to avoid substring matching
    # (e.g. 'unfavourable' contains 'favourable', 'non response' contains 'response')
    if any(tok in v_clean for tok in ["unfavourable", "unfavorable", "progressive disease", "residual disease", "non mpr", "npcr", "nmpr", "non responder", "non response"]):
        return "non-responder", raw_str
    if v_clean in {"pd", "non responder", "non response", "nr", "rd", "low"}:
        return "non-responder", raw_str

    if any(tok in v_clean for tok in ["favourable", "favorable", "complete response", "partial response", "pcr", "mpr"]):
        return "responder", raw_str
    if v_clean in {"cr", "pr", "responder", "response", "r", "high"}:
        return "responder", raw_str

    if any(tok in v_clean for tok in ["stable disease", "stable"]):
        return "stable", raw_str
    if v_clean in {"sd", "medium"}:
        return "stable", raw_str

    return "not-evaluable", raw_str


def _standardize_tier0_obs(
    adata: ad.AnnData,
    accession: str,
    indication: str,
    technology: str,
    cell_selection: str,
    patient_col: str | None = None,
    sample_col: str | None = None,
    response_col: str | None = None,
    response_map: Mapping[str, str] | None = None,
    treatment_col: str | None = None,
    cell_type_col: str | None = None,
) -> ad.AnnData:
    """Standardize .obs columns while preserving all existing columns purely."""
    obs = adata.obs.copy()
    obs["cell_id"] = adata.obs_names.astype(str)
    obs["dataset_id"] = accession
    obs["indication"] = indication
    obs["cancer_type"] = indication
    obs["technology"] = technology
    obs["cell_selection_strategy"] = cell_selection

    # Map patient/donor
    if patient_col and patient_col in obs.columns:
        obs["patient_id"] = obs[patient_col].astype(str)
    elif "patient_id" not in obs.columns:
        for c in ("PMID_donor_id", "patient", "patientID", "donor", "donor_id", "subject", "sample"):
            if c in obs.columns:
                obs["patient_id"] = obs[c].astype(str)
                break
        else:
            obs["patient_id"] = f"{accession}_donor_unspecified"

    # Map sample
    if sample_col and sample_col in obs.columns:
        obs["sample_id"] = obs[sample_col].astype(str)
    elif "sample_id" not in obs.columns:
        obs["sample_id"] = obs["patient_id"]

    # Map response
    raw_response_series = None
    if response_col and response_col in obs.columns:
        raw_response_series = obs[response_col]
    elif "clinical_response_raw" in obs.columns:
        raw_response_series = obs["clinical_response_raw"]
    elif "clinical_response" in obs.columns:
        raw_response_series = obs["clinical_response"]
    else:
        for c in ("response", "Response", "RECIST", "benefit", "responder", "outcome", "Combined_outcome", "pCR_status", "Pathologic_Response", "pathological_response", "Path_response"):
            if c in obs.columns:
                raw_response_series = obs[c]
                break

    if raw_response_series is not None:
        obs["clinical_response_raw"] = raw_response_series.astype(str)
        if response_map:
            obs["clinical_response"] = obs["clinical_response_raw"].map(response_map).fillna("not-evaluable")
        else:
            harmonized = [_harmonize_single_response(v) for v in raw_response_series]
            obs["clinical_response"] = [h[0] for h in harmonized]
            obs["clinical_response_raw"] = [h[1] for h in harmonized]
    else:
        obs["clinical_response"] = "not-evaluable"
        obs["clinical_response_raw"] = "Unspecified"

    # Map treatment
    if treatment_col and treatment_col in obs.columns:
        obs["treatment_status"] = obs[treatment_col].astype(str)
    elif "treatment_status" not in obs.columns:
        for c in ("treatment", "treatment_status", "timepoint", "Timepoint", "cycle"):
            if c in obs.columns:
                obs["treatment_status"] = obs[c].astype(str)
                break
        else:
            obs["treatment_status"] = "Pre-treatment"

    # Map cell type
    if cell_type_col and cell_type_col in obs.columns:
        obs["cell_type"] = obs[cell_type_col].astype(str)
    elif "cell_type" not in obs.columns:
        for c in ("cell_type", "celltype", "CellType", "major_cell_type", "cluster"):
            if c in obs.columns:
                obs["cell_type"] = obs[c].astype(str)
                break
        else:
            obs["cell_type"] = "Unspecified"

    # Clean all object columns to prevent HDF5 serialization failures with mixed types
    for c in obs.columns:
        if obs[c].dtype == "object":
            obs[c] = obs[c].fillna("NA").astype(str)

    new_adata = adata.copy()
    new_adata.obs = obs
    new_adata.obs_names = obs["cell_id"]
    new_adata.obs_names.name = "cell_id"
    new_adata.var_names_make_unique()
    return new_adata


def _extract_tar_safely(tar_path: Path, extract_dir: Path) -> Result[Path, str]:
    """Safely extracts a tar archive avoiding path traversal vulnerabilities."""
    if not tar_path.exists():
        return Failure(f"Tar archive not found at {tar_path}")
    try:
        extract_dir.mkdir(parents=True, exist_ok=True)
        # Check if already extracted
        if any(extract_dir.iterdir()):
            return Success(extract_dir)

        with tarfile.open(tar_path, "r:*") as archive:
            for member in archive.getmembers():
                dest_path = (extract_dir / member.name).resolve()
                if not str(dest_path).startswith(str(extract_dir.resolve())):
                    return Failure(f"Path traversal detected in archive: {member.name}")
            archive.extractall(extract_dir)
        return Success(extract_dir)
    except Exception as exc:
        return Failure(f"Failed to extract {tar_path.name}: {exc}")


def _parse_single_geo_series_matrix_file(matrix_file: Path) -> Result[pd.DataFrame, str]:
    """Parses a single NCBI GEO Series Matrix file into a DataFrame."""
    try:
        with gzip.open(matrix_file, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        gsm_list: list[str] = []
        title_list: list[str] = []
        characteristics: dict[str, list[str]] = {}

        for line in lines:
            line_str = line.strip()
            if line_str.startswith("!Sample_geo_accession"):
                gsm_list = [t.strip().strip('"') for t in line_str.split("\t")[1:]]
            elif line_str.startswith("!Sample_title"):
                title_list = [t.strip().strip('"') for t in line_str.split("\t")[1:]]
            elif line_str.startswith("!Sample_characteristics_ch1"):
                tokens = [t.strip().strip('"') for t in line_str.split("\t")[1:]]
                for i, tok in enumerate(tokens):
                    if ":" in tok:
                        k, v = tok.split(":", 1)
                        k = k.strip()
                        v = v.strip()
                        if k not in characteristics:
                            characteristics[k] = ["" for _ in range(len(tokens))]
                        characteristics[k][i] = v

        if not gsm_list:
            return Failure(f"No samples found in series matrix {matrix_file.name}")

        meta_df = pd.DataFrame(index=gsm_list)
        meta_df["Sample_geo_accession"] = gsm_list
        if title_list and len(title_list) == len(gsm_list):
            meta_df["Sample_title"] = title_list

        for k, v in characteristics.items():
            if len(v) == len(gsm_list):
                meta_df[k] = v

        return Success(meta_df)
    except Exception as exc:
        return Failure(f"Failed to parse series matrix {matrix_file.name}: {exc}")


def _parse_geo_series_matrix(raw_dir: Path, gse_id: str) -> Result[pd.DataFrame, str]:
    """Parses NCBI GEO Series Matrix sample annotations into a DataFrame indexed by GSM and title."""
    matrix_files = sorted(raw_dir.glob(f"{gse_id}*series_matrix.txt.gz"))
    if not matrix_files:
        matrix_file = raw_dir / f"{gse_id}_series_matrix.txt.gz"
        # Fallback to web fetch if matrix file is missing
        prefix = gse_id[:-3]
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Python tme_datasets)"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                matrix_file.write_bytes(resp.read())
            matrix_files = [matrix_file]
        except Exception as exc:
            return Failure(f"Could not locate or fetch GEO series matrix for {gse_id}: {exc}")

    dfs: list[pd.DataFrame] = []
    for mf in matrix_files:
        res = _parse_single_geo_series_matrix_file(mf)
        if isinstance(res, Success):
            dfs.append(res.unwrap())
        else:
            logger.warning("Could not parse %s: %s", mf.name, res.failure())

    if not dfs:
        return Failure(f"No valid series matrix parsed for {gse_id}")

    combined_df = dfs[0] if len(dfs) == 1 else pd.concat(dfs, axis=0)
    return Success(combined_df)


def _load_10x_triplets_or_h5(
    extract_dir: Path,
    min_counts: int = 10,
    min_genes: int = 1,
) -> Result[list[ad.AnnData], str]:
    """Discovers and parses 10x .h5 files or (matrix.mtx, barcodes, features) triplets."""
    try:
        # 1. Look for .h5.gz and decompress
        for hg in extract_dir.glob("*.h5.gz"):
            decomp = hg.with_suffix("")
            if not decomp.exists():
                with gzip.open(hg, "rb") as f_in, open(decomp, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        # 2. Check for .h5 files
        h5_files = sorted(list(extract_dir.glob("*.h5"))) or sorted(list(extract_dir.glob("**/*.h5")))
        adatas: list[ad.AnnData] = []
        if h5_files:
            for hf in h5_files:
                sample_name = hf.name.replace(".h5", "").replace("filtered_feature_bc_matrix_", "").replace("raw_feature_bc_matrix_", "")
                try:
                    sub_a = sc.read_10x_h5(hf)
                    sub_a.obs_names = [f"{sample_name}_{b}" for b in sub_a.obs_names]
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = sample_name
                    patient = sample_name.split("_")[0]
                    sub_a.obs["patient_id"] = patient
                    adatas.append(sub_a)
                except Exception as exc:
                    logger.warning("Failed to parse h5 file %s: %s", hf.name, exc)
            if adatas:
                return Success(adatas)

        # 3. Check for MTX triplets
        mtx_candidates = (
            list(extract_dir.glob("*matrix.mtx.gz"))
            + list(extract_dir.glob("**/*matrix.mtx.gz"))
            + list(extract_dir.glob("*.mtx.gz"))
            + list(extract_dir.glob("**/*.mtx.gz"))
            + list(extract_dir.glob("*.mtx"))
            + list(extract_dir.glob("**/*.mtx"))
        )
        mtx_files = sorted(list(dict.fromkeys(mtx_candidates)))
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "").replace(".matrix.mtx.gz", "").replace("matrix.mtx.gz", "").replace(".mtx.gz", "").replace(".mtx", "")
            bc_p = mp.parent / f"{prefix}_barcodes.tsv.gz"
            feat_p = mp.parent / f"{prefix}_features.tsv.gz"
            if not bc_p.exists():
                bc_candidates = (
                    list(mp.parent.glob(f"*{prefix}*barcode*"))
                    + list(mp.parent.glob("*barcode*"))
                    + list(extract_dir.glob(f"*{prefix}*barcode*"))
                    + list(extract_dir.glob("*barcode*"))
                )
                if bc_candidates:
                    bc_p = bc_candidates[0]
            if not feat_p.exists():
                feat_candidates = (
                    list(mp.parent.glob(f"*{prefix}*feature*"))
                    + list(mp.parent.glob(f"*{prefix}*gene*"))
                    + list(mp.parent.glob("*feature*"))
                    + list(mp.parent.glob("*gene*"))
                    + list(extract_dir.glob(f"*{prefix}*feature*"))
                    + list(extract_dir.glob(f"*{prefix}*gene*"))
                    + list(extract_dir.glob("*feature*"))
                    + list(extract_dir.glob("*gene*"))
                )
                if feat_candidates:
                    feat_p = feat_candidates[0]

            if not (bc_p.exists() and feat_p.exists()):
                continue

            mat = _to_csr(sio.mmread(mp))
            try:
                bcs = pd.read_csv(bc_p, header=None, sep=None, engine="python")[0].astype(str).tolist()
                feats = pd.read_csv(feat_p, header=None, sep=None, engine="python")
            except Exception:
                bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
                feats = pd.read_csv(feat_p, header=None, sep="\t")

            # Check orientation (10x default is genes x cells)
            if mat.shape[0] != len(bcs) and mat.shape[1] == len(bcs):
                mat = mat.T

            # Check feature type column if 3+ cols (e.g. Gene Expression vs Antibody Capture)
            if len(feats.columns) >= 3:
                rna_mask = feats[2].astype(str) == "Gene Expression"
                if np.any(rna_mask):
                    feats = feats[rna_mask].reset_index(drop=True)
                    mat = mat[:, rna_mask.values]

            # Pre-filter empty droplets
            counts_per_cell = np.asarray(mat.sum(axis=1)).ravel()
            genes_per_cell = np.asarray((mat > 0).sum(axis=1)).ravel()
            valid_cells = (counts_per_cell >= min_counts) & (genes_per_cell >= min_genes)
            if not np.any(valid_cells):
                continue
            mat = mat[valid_cells]
            bcs = [bcs[i] for i, v in enumerate(valid_cells) if v]

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            sample_id = prefix.split("_")[0]
            sub_a = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_id
            sub_a.obs["patient_id"] = sample_id
            adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable 10x sample matrices found in {extract_dir}")
        return Success(adatas)
    except Exception as exc:
        return Failure(f"Failed to load 10x matrices from {extract_dir}: {exc}")


def _load_cellxgene_filtered_h5ad(
    h5ad_path: Path,
    cell_filter_fn: Callable[[pd.DataFrame], pd.Series | np.ndarray] | None = None,
) -> Result[ad.AnnData, str]:
    """Memory-safe loader for large CELLxGENE files using backed slicing before RAM allocation."""
    if not h5ad_path.exists():
        return Failure(f"H5AD file not found at {h5ad_path}")
    try:
        backed_adata = ad.read_h5ad(h5ad_path, backed="r")
        if cell_filter_fn is not None:
            mask = cell_filter_fn(backed_adata.obs)
            if not np.any(mask):
                return Failure("Cell filter yielded 0 matching cells")
            # Slice in memory
            sub_adata = backed_adata[mask].to_memory()
        else:
            sub_adata = backed_adata.to_memory()
        return Success(sub_adata)
    except Exception as exc:
        return Failure(f"Failed to load CELLxGENE H5AD {h5ad_path.name}: {exc}")


def _apply_subset_and_subsample(
    adata: ad.AnnData,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
    seed: int = 42,
) -> ad.AnnData:
    """Purely filters and/or downsamples an AnnData object."""
    res_adata = adata

    if subset:
        mask = np.ones(res_adata.n_obs, dtype=bool)
        for col, val in subset.items():
            if col in res_adata.obs.columns:
                mask = mask & (res_adata.obs[col].astype(str).str.lower() == str(val).lower())
            else:
                logger.warning("Subset column '%s' not present in .obs (available: %s)", col, list(res_adata.obs.columns))
        if np.any(mask):
            res_adata = res_adata[mask].copy()

    if subsample_n and subsample_n < res_adata.n_obs:
        rng = np.random.default_rng(seed)
        indices = np.sort(rng.choice(res_adata.n_obs, size=subsample_n, replace=False))
        res_adata = res_adata[indices].copy()

    return res_adata


def _restore_raw_anndata(adata: ad.AnnData) -> ad.AnnData:
    """Safely restore AnnData from adata.raw if present, preserving obs, obsm, and uns."""
    if adata.raw is not None and adata.raw.X is not None:
        try:
            obs = adata.obs.copy()
            obsm = dict(adata.obsm)
            uns = dict(adata.uns)
            var = adata.raw.var.copy() if (hasattr(adata.raw, "var") and len(adata.raw.var) > 0) else adata.var.copy()
            raw_X = adata.raw.X
            new_adata = ad.AnnData(X=raw_X, obs=obs, var=var, obsm=obsm, uns=uns)
            return new_adata
        except Exception:
            return adata
    return adata



# =========================================================================
# 1. GSE246613 (Breast — Triple Negative, Pembrolizumab + RT)
# =========================================================================
def load_gse246613_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE246613 Breast Cancer atlas (~532k cells, 266 patients)."""
    imm_gz = raw_dir / "GSE246613_PembroRT_immune_R100_final.h5ad.gz"
    non_imm_gz = raw_dir / "GSE246613_PembroRT_non_immune_cells.h5ad.gz"
    comb_gz = raw_dir / "GSE246613_combined_RTPDv4_scvi_celltypist.h5ad.gz"

    scratch_dir = raw_dir / "decompressed"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    try:
        adatas: list[ad.AnnData] = []
        if imm_gz.exists() and non_imm_gz.exists():
            for gz_path in (imm_gz, non_imm_gz):
                decomp_path = scratch_dir / gz_path.name.replace(".gz", "")
                if not decomp_path.exists():
                    logger.info("Decompressing %s...", gz_path.name)
                    with gzip.open(gz_path, "rb") as f_in, open(decomp_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                sub_a = ad.read_h5ad(decomp_path)
                sub_a.var_names_make_unique()
                adatas.append(sub_a)
            adata = ad.concat(adatas, join="outer")
        elif comb_gz.exists():
            decomp_path = scratch_dir / comb_gz.name.replace(".gz", "")
            if not decomp_path.exists():
                with gzip.open(comb_gz, "rb") as f_in, open(decomp_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            adata = ad.read_h5ad(decomp_path)
        else:
            # Fallback to any .h5ad
            h5ad_files = list(raw_dir.glob("*.h5ad"))
            if not h5ad_files:
                return Failure(f"GSE246613 H5AD files not found in {raw_dir}")
            adata = ad.read_h5ad(h5ad_files[0])

        adata = _restore_raw_anndata(adata)

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE246613",
            indication="Breast",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="patient",
            response_col="responder",
            treatment_col="cycle",
            cell_type_col="cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE246613: {exc}")


# =========================================================================
# 2. GSE300475 (Breast — HR+ Immunotherapy Longitudinal)
# =========================================================================
def load_gse300475_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE300475 Breast Cancer longitudinal atlas (~64k cells, 32 patients)."""
    tar_path = raw_dir / "GSE300475_RAW.tar"
    xlsx_path = raw_dir / "GSE300475_feature_ref.xlsx"
    extract_dir = raw_dir / "unpacked_gse300475"

    if not tar_path.exists():
        return Failure(f"GSE300475 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        # Look for 10x MTX files
        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
        if not mtx_files:
            mtx_files = sorted(list(extract_dir.glob("**/*matrix.mtx.gz")))

        meta_df: pd.DataFrame | None = None
        if xlsx_path.exists():
            try:
                meta_df = pd.read_excel(xlsx_path)
            except Exception:
                pass

        adatas: list[ad.AnnData] = []
        for mtx_p in mtx_files:
            prefix = mtx_p.name.replace("_matrix.mtx.gz", "")
            bc_p = mtx_p.parent / f"{prefix}_barcodes.tsv.gz"
            feat_p = mtx_p.parent / f"{prefix}_features.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                feat_p = mtx_p.parent / f"{prefix}_genes.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                continue

            mat = _to_csr(sio.mmread(mtx_p))
            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            gsm_id = prefix.split("_")[0]
            sub_adata = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{gsm_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_adata.var_names_make_unique()
            sub_adata.obs["sample_id"] = gsm_id
            sub_adata.obs["patient_id"] = gsm_id
            adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid sample MTX triplets parsed in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        if meta_df is not None and "sample_id" in meta_df.columns:
            meta_indexed = meta_df.set_index("sample_id")
            combined.obs = combined.obs.join(meta_indexed, on="sample_id", how="left")

        combined = _standardize_tier0_obs(
            combined,
            accession="GSE300475",
            indication="Breast",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE300475: {exc}")


# =========================================================================
# 3. GSE212707 (Breast — Multiome CAF & TME Architecture)
# =========================================================================
def load_gse212707_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE212707 Breast Multiome CAF atlas (~42k cells, 21 patients)."""
    tar_path = raw_dir / "GSE212707_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse212707"

    if not tar_path.exists():
        return Failure(f"GSE212707 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        # Look for sub-archives (*.tar.gz) e.g. GSM6543819_C1.tar.gz
        sub_tars = sorted(list(extract_dir.glob("*.tar.gz")))
        if sub_tars:
            for st in sub_tars:
                sub_dir = extract_dir / st.name.replace(".tar.gz", "")
                if not sub_dir.exists():
                    sub_dir.mkdir(parents=True, exist_ok=True)
                    try:
                        with tarfile.open(st, "r:*") as sub_arch:
                            sub_arch.extractall(sub_dir)
                    except Exception:
                        pass

        # Look for .h5 files or MTX files recursively
        h5_files = sorted(list(extract_dir.glob("**/*.h5")))
        adatas: list[ad.AnnData] = []

        if h5_files:
            for hf in h5_files:
                sample_name = hf.name.split("_")[0]
                try:
                    sub_a = sc.read_10x_h5(hf, gex_only=True)
                    sub_a.obs_names = [f"{sample_name}_{b}" for b in sub_a.obs_names]
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = sample_name
                    sub_a.obs["patient_id"] = sample_name
                    adatas.append(sub_a)
                except Exception:
                    pass
        else:
            dirs_with_mtx = [d for d in extract_dir.glob("**/") if list(d.glob("*matrix.mtx*"))]
            for d in dirs_with_mtx:
                try:
                    sub_a = sc.read_10x_mtx(d, gex_only=True)
                    sample_name = d.parent.name if d.name in ("C1", "C2", "C3", "S1", "S2", "S3", "L1", "L2", "L4") else d.name
                    sub_a.obs_names = [f"{sample_name}_{b}" for b in sub_a.obs_names]
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = sample_name
                    sub_a.obs["patient_id"] = sample_name
                    adatas.append(sub_a)
                except Exception:
                    pass

            if not adatas:
                mtx_files = sorted(list(extract_dir.glob("**/*matrix.mtx*")))
                for mp in mtx_files:
                    bc_candidates = list(mp.parent.glob("*barcodes.tsv*"))
                    feat_candidates = list(mp.parent.glob("*features.tsv*"))
                    if not (bc_candidates and feat_candidates):
                        continue
                    mat = _to_csr(sio.mmread(mp))
                    bcs = pd.read_csv(bc_candidates[0], header=None, sep="\t")[0].astype(str).tolist()
                    feats = pd.read_csv(feat_candidates[0], header=None, sep="\t")
                    var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
                    var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values
                    sample_name = mp.parent.name
                    sub_a = ad.AnnData(
                        X=mat.astype(np.float32),
                        obs=pd.DataFrame(index=[f"{sample_name}_{b}" for b in bcs]),
                        var=var_df,
                    )
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = sample_name
                    sub_a.obs["patient_id"] = sample_name
                    adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable multiome files in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE212707",
            indication="Breast",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE212707: {exc}")


# =========================================================================
# 4. GSE236581 (CRC — Longitudinal Anti-PD-1 Spatiotemporal Dynamics)
# =========================================================================
def load_gse236581_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE236581 Colorectal Cancer ICB atlas (~338k cells, 169 patients)."""
    mtx_p = raw_dir / "GSE236581_counts.mtx.gz"
    bc_p = raw_dir / "GSE236581_barcodes.tsv.gz"
    feat_p = raw_dir / "GSE236581_features.tsv.gz"
    meta_p = raw_dir / "GSE236581_CRC-ICB_metadata.txt.gz"

    if not mtx_p.exists():
        return Failure(f"GSE236581 counts matrix not found at {mtx_p}")

    try:
        mat = _to_csr(sio.mmread(mtx_p))
        if bc_p.exists():
            barcodes = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
        else:
            barcodes = [f"cell_{i}" for i in range(mat.shape[0])]

        if feat_p.exists():
            feats = pd.read_csv(feat_p, header=None, sep="\t")
            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values
        else:
            var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(mat.shape[1])])
            var_df["gene_name"] = var_df.index

        obs_df = pd.DataFrame(index=barcodes)
        if meta_p.exists():
            try:
                meta = pl.read_csv(meta_p, separator="\t").to_pandas()
                first_col = meta.columns[0]
                meta = meta.set_index(first_col)
                obs_df = obs_df.join(meta, how="left")
            except Exception:
                pass

        adata = ad.AnnData(X=mat.astype(np.float32), obs=obs_df, var=var_df)
        adata = _standardize_tier0_obs(
            adata,
            accession="GSE236581",
            indication="CRC",
            technology="High-Throughput scRNA-seq",
            cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
            patient_col="Patient" if "Patient" in adata.obs.columns else "patient",
            response_col="Response" if "Response" in adata.obs.columns else "response",
            treatment_col="Timepoint" if "Timepoint" in adata.obs.columns else "timepoint",
            cell_type_col="CellType" if "CellType" in adata.obs.columns else "cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE236581: {exc}")


# =========================================================================
# 5. GSE299651 (CRC — Pooled Screening & ICB Response)
# =========================================================================
def load_gse299651_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE299651 Colorectal Cancer pool atlas (~40k cells, 20 patients)."""
    h5_files = sorted(list(raw_dir.glob("*.h5")))
    if not h5_files:
        return Failure(f"GSE299651 .h5 files not found in {raw_dir}")

    try:
        adatas: list[ad.AnnData] = []
        for hf in h5_files:
            sub_a = sc.read_10x_h5(hf, gex_only=True)
            pool_tag = hf.name.split("_")[2] if len(hf.name.split("_")) > 2 else hf.stem
            sub_a.obs_names = [f"{pool_tag}_{b}" for b in sub_a.obs_names]
            sub_a.var_names_make_unique()
            sub_a.obs["pool_id"] = pool_tag
            sub_a.obs["patient_id"] = pool_tag
            sub_a.obs["sample_id"] = pool_tag
            adatas.append(sub_a)

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE299651",
            indication="CRC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE299651: {exc}")


# =========================================================================
# 6. CELLxGENE_829a3cd1 (CRC — Metastatic Plasticity & ICB Cohort)
# =========================================================================
def load_cellxgene_829a3cd1_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load CELLxGENE 829a3cd1 Colorectal Cancer atlas (~47k cells, 29 patients)."""
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"CELLxGENE 829a3cd1 H5AD not found in {raw_dir}")

    try:
        adata = ad.read_h5ad(h5ad_files[0])
        adata = _restore_raw_anndata(adata)

        adata = _standardize_tier0_obs(
            adata,
            accession="CELLxGENE_829a3cd1",
            indication="CRC",
            technology="10x Chromium 3' v3/v3.1",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="donor_id",
            cell_type_col="cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load CELLxGENE 829a3cd1: {exc}")


# =========================================================================
# 7. GSE270680 (Gastric — Spatially Resolved Immunotherapy Atlas)
# =========================================================================
def load_gse270680_gastric(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE270680 Gastric Cancer atlas (~154k cells, 77 patients)."""
    tar_path = raw_dir / "GSE270680_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse270680"

    if not tar_path.exists():
        return Failure(f"GSE270680 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "")
            bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
            feat_p = extract_dir / f"{prefix}_features.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                continue
            mat = _to_csr(sio.mmread(mp))
            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")
            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values
            sample_id = prefix.split("_")[0]
            sub_a = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_id
            sub_a.obs["patient_id"] = sample_id
            adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable sample matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE270680",
            indication="Gastric",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE270680: {exc}")


# =========================================================================
# 8. GSE313642 (HCC — Phase II Sorafenib + Nivolumab Trial)
# =========================================================================
def load_gse313642_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE313642 HCC CITE-seq trial atlas (~388k cells, 194 patients)."""
    tar_path = raw_dir / "GSE313642_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse313642"

    if not tar_path.exists():
        return Failure(f"GSE313642 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "")
            bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
            feat_p = extract_dir / f"{prefix}_features.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            # Check feature type column if 3+ cols (e.g. Gene Expression vs Antibody Capture)
            if len(feats.columns) >= 3:
                rna_mask = feats[2].astype(str) == "Gene Expression"
                if np.any(rna_mask):
                    feats_rna = feats[rna_mask]
                    mat = mat[:, rna_mask.values]
                else:
                    feats_rna = feats
            else:
                feats_rna = feats

            var_df = pd.DataFrame(index=feats_rna[0].astype(str).tolist())
            var_df["gene_name"] = (feats_rna[1] if len(feats_rna.columns) > 1 else feats_rna[0]).astype(str).values

            sample_id = prefix.split("_")[0]
            sub_a = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_id
            sub_a.obs["patient_id"] = sample_id
            adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable CITE-seq matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE313642",
            indication="HCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE313642: {exc}")


# =========================================================================
# 9. GSE245906 (HCC — Innate Myeloid & Monocyte Response)
# =========================================================================
def load_gse245906_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE245906 HCC innate myeloid atlas (~40k cells, 20 patients)."""
    count_p = raw_dir / "GSE245906_Giraud_Chalopin_innate_HCC_processed_count_data_mat.tsv.gz"
    meta_p = raw_dir / "GSE245906_Giraud_Chalopin_innate_HCC_metadata.tsv.gz"

    if not count_p.exists():
        return Failure(f"GSE245906 count matrix not found at {count_p}")

    try:
        with gzip.open(count_p, "rt") as f:
            header_line = f.readline().strip().split("\t")

        schema_cols = ["gene_name"] + header_line
        overrides = {c: pl.Float32 for c in header_line}
        overrides["gene_name"] = pl.String

        df_counts = pl.read_csv(
            count_p,
            separator="\t",
            has_header=False,
            skip_rows=1,
            new_columns=schema_cols,
            schema_overrides=overrides,
        )
        genes = df_counts["gene_name"].to_list()
        barcodes = header_line
        mat = df_counts.select(barcodes).to_numpy().astype(np.float32).T
        sparse_mat = sp.csr_matrix(mat, dtype=np.float32)
        del df_counts, mat

        obs_df = pd.DataFrame(index=barcodes)
        if meta_p.exists():
            try:
                meta = pl.read_csv(meta_p, separator="\t").to_pandas()
                m_first = meta.columns[0]
                meta = meta.set_index(m_first)
                obs_df = obs_df.join(meta, how="left")
            except Exception:
                pass

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes
        adata = ad.AnnData(X=sparse_mat, obs=obs_df, var=var_df)

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE245906",
            indication="HCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="patient" if "patient" in adata.obs.columns else "donor",
            response_col="response" if "response" in adata.obs.columns else None,
            cell_type_col="cell_type" if "cell_type" in adata.obs.columns else "celltype",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE245906: {exc}")


# =========================================================================
# 10. GSE301741 (HNSCC — Immunotherapy Response & Tertiary Lymphoids)
# =========================================================================
def load_gse301741_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE301741 HNSCC atlas (~116k cells, 58 patients)."""
    tar_path = raw_dir / "GSE301741_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse301741"

    if not tar_path.exists():
        return Failure(f"GSE301741 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        h5_files = sorted(list(extract_dir.glob("*.h5"))) or sorted(list(extract_dir.glob("**/*.h5")))
        adatas: list[ad.AnnData] = []

        if h5_files:
            for hf in h5_files:
                sample_name = hf.name.split("_")[0]
                try:
                    sub_a = sc.read_10x_h5(hf)
                    sub_a.obs_names = [f"{sample_name}_{b}" for b in sub_a.obs_names]
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = sample_name
                    patient = hf.name.split("_")[1] if len(hf.name.split("_")) > 1 else sample_name
                    sub_a.obs["patient_id"] = patient.split("-")[0]
                    adatas.append(sub_a)
                except Exception:
                    pass
        else:
            mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
            for mp in mtx_files:
                prefix = mp.name.replace("_matrix.mtx.gz", "")
                bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
                feat_p = extract_dir / f"{prefix}_features.tsv.gz"
                if not (bc_p.exists() and feat_p.exists()):
                    continue

                mat = _to_csr(sio.mmread(mp))
                bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
                feats = pd.read_csv(feat_p, header=None, sep="\t")

                var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
                var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

                sample_id = prefix.split("_")[0]
                sub_a = ad.AnnData(
                    X=mat.astype(np.float32),
                    obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                    var=var_df,
                )
                sub_a.var_names_make_unique()
                sub_a.obs["sample_id"] = sample_id
                sub_a.obs["patient_id"] = sample_id
                adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable sample matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE301741",
            indication="HNSCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE301741: {exc}")


# =========================================================================
# 11. GSE287301 (HNSCC — Infiltrating Architecture & CosMx Spatial)
# =========================================================================
def load_gse287301_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE287301 HNSCC atlas (~96k cells, 48 patients)."""
    tar_path = raw_dir / "GSE287301_filtered_feature_bc_matrix.tar.gz"
    extract_dir = raw_dir / "unpacked_gse287301"

    if not tar_path.exists():
        return Failure(f"GSE287301 feature matrix archive not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        # Look for 10x folders or MTX files
        dirs = [d for d in extract_dir.glob("**/") if list(d.glob("*matrix.mtx*"))]
        adatas: list[ad.AnnData] = []
        for d in dirs:
            try:
                sub_a = sc.read_10x_mtx(d, gex_only=True)
                counts_per_cell = np.asarray(sub_a.X.sum(axis=1)).ravel()
                genes_per_cell = np.asarray((sub_a.X > 0).sum(axis=1)).ravel()
                valid_cells = (counts_per_cell >= 100) & (genes_per_cell >= 50)
                sub_a = sub_a[valid_cells].copy()
                sub_a.obs_names = [f"{d.name}_{b}" for b in sub_a.obs_names]
                sub_a.var_names_make_unique()
                sub_a.obs["sample_id"] = d.name
                sub_a.obs["patient_id"] = d.name
                adatas.append(sub_a)
            except Exception:
                pass

        if not adatas:
            # Fallback to direct mtx search
            mtx_files = sorted(list(extract_dir.glob("**/*matrix.mtx*")))
            for mp in mtx_files:
                prefix = mp.name.replace("_matrix.mtx.gz", "").replace("matrix.mtx.gz", "")
                bc_p = list(mp.parent.glob("*barcodes.tsv*"))[0]
                feat_p = list(mp.parent.glob("*features.tsv*"))[0]
                mat = _to_csr(sio.mmread(mp))
                bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
                feats = pd.read_csv(feat_p, header=None, sep="\t")
                var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
                var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values
                sample_id = mp.parent.name
                sub_a = ad.AnnData(
                    X=mat.astype(np.float32),
                    obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                    var=var_df,
                )
                sub_a.var_names_make_unique()
                sub_a.obs["sample_id"] = sample_id
                sub_a.obs["patient_id"] = sample_id
                adatas.append(sub_a)

        if not adatas:
            return Failure(f"No 10x feature matrices parsed from {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE287301",
            indication="HNSCC",
            technology="Subcellular Spatial Transcriptomics (CosMx/Xenium)",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE287301: {exc}")


# =========================================================================
# 12. GSE200996 (HNSCC — Tissue-Resident Memory T Cells in Anti-PD-1)
# =========================================================================
def load_gse200996_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE200996 HNSCC TIL atlas (~408k cells, 204 patients) with verified Path_response."""
    tar_path = raw_dir / "GSE200996_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse200996"

    ext_res = _extract_tar_safely(tar_path, extract_dir)
    match ext_res:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    try:
        # Load and combine available metadata tables (tumor and PBMC)
        meta_dfs: list[pd.DataFrame] = []
        for mf in sorted(list(raw_dir.glob("*meta*.txt.gz")) + list(raw_dir.glob("*meta*.tsv.gz"))):
            try:
                df = pd.read_csv(mf, sep="\t")
                b_col = df.columns[0]
                df = df.set_index(b_col)
                meta_dfs.append(df)
            except Exception as exc:
                logger.warning("Could not read %s: %s", mf.name, exc)

        if meta_dfs:
            combined_meta = pd.concat(meta_dfs)
            combined_meta = combined_meta[~combined_meta.index.duplicated(keep="first")]
            lower_meta_map = {idx.lower(): idx for idx in combined_meta.index}
        else:
            combined_meta = pd.DataFrame()
            lower_meta_map = {}

        h5_files = sorted(list(extract_dir.glob("*.h5"))) or sorted(list(extract_dir.glob("**/*.h5")))
        mtx_files = sorted(list(extract_dir.glob("**/*matrix.mtx*"))) if not h5_files else []
        if not h5_files and not mtx_files:
            return Failure(f"No .h5 or .mtx files found in {extract_dir}")

        adatas: list[ad.AnnData] = []
        if not h5_files and mtx_files:
            for mp in mtx_files:
                try:
                    fname = mp.name.replace("_matrix.mtx.gz", "").replace("matrix.mtx.gz", "")
                    tokens = fname.split("_")
                    patient = tokens[1] if len(tokens) > 1 else (tokens[0] if tokens else "unknown")
                    stage = tokens[2] if len(tokens) > 2 else "unspecified"
                    bc_files = list(mp.parent.glob("*barcodes.tsv*"))
                    feat_files = list(mp.parent.glob("*features.tsv*"))
                    if not bc_files or not feat_files:
                        continue
                    mat = _to_csr(sio.mmread(mp))
                    bcs = pd.read_csv(bc_files[0], header=None, sep="\t")[0].astype(str).tolist()
                    feats = pd.read_csv(feat_files[0], header=None, sep="\t")
                    gene_names = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values
                    var_df = pd.DataFrame(index=gene_names)
                    sub_a = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
                    sub_a.var_names_make_unique()
                    sub_a.obs["sample_id"] = f"{patient}_{stage}"
                    sub_a.obs["patient_id"] = patient
                    sub_a.obs["stage"] = stage
                    adatas.append(sub_a)
                except Exception as exc:
                    logger.warning("Could not read mtx %s: %s", mp.name, exc)
        else:
            for hf in h5_files:
                try:
                    fname = hf.name.replace(".h5", "")
                    if "feature_bc_matrix_" in fname:
                        specimen_part = fname.split("feature_bc_matrix_")[-1]
                    else:
                        specimen_part = fname
                    tokens = specimen_part.split("_")
                    patient = tokens[0] if len(tokens) > 0 else "unknown"
                    stage = tokens[1] if len(tokens) > 1 else "unspecified"

                    sub_a = sc.read_10x_h5(hf)
                    sub_a.var_names_make_unique()

                    clean_bcs = [b.split("-")[0] for b in sub_a.obs_names]
                    candidate_keys = [f"{bc}_{patient}_{stage}".lower() for bc in clean_bcs]

                    if lower_meta_map:
                        matched_indices = [i for i, k in enumerate(candidate_keys) if k in lower_meta_map]
                        if matched_indices:
                            sub_a = sub_a[matched_indices].copy()
                            matched_keys = [candidate_keys[i] for i in matched_indices]
                            canonical_bcs = [lower_meta_map[k] for k in matched_keys]
                            sub_a.obs_names = canonical_bcs
                        else:
                            counts = np.asarray(sub_a.X.sum(axis=1)).ravel()
                            keep = counts >= 200
                            if np.any(keep):
                                sub_a = sub_a[keep].copy()
                                sub_a.obs_names = [f"{patient}_{stage}_{b}" for b in sub_a.obs_names]
                            else:
                                continue
                    else:
                        counts = np.asarray(sub_a.X.sum(axis=1)).ravel()
                        keep = counts >= 200
                        if np.any(keep):
                            sub_a = sub_a[keep].copy()
                            sub_a.obs_names = [f"{patient}_{stage}_{b}" for b in sub_a.obs_names]
                        else:
                            continue

                    sub_a.obs["sample_id"] = f"{patient}_{stage}"
                    sub_a.obs["patient_id"] = patient
                    sub_a.obs["stage"] = stage
                    adatas.append(sub_a)
                except Exception as exc:
                    logger.warning("Could not read h5 %s: %s", hf.name, exc)

        if not adatas:
            return Failure(f"Failed to parse any single-cell data from {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        if not combined_meta.empty:
            obs_df = combined.obs.join(combined_meta, how="left")
            combined.obs = obs_df

        combined = _standardize_tier0_obs(
            combined,
            accession="GSE200996",
            indication="HNSCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
            sample_col="sample_id",
            patient_col="Patient_ID" if "Patient_ID" in combined.obs.columns else "patient_id",
            response_col="Path_response" if "Path_response" in combined.obs.columns else None,
            treatment_col="Stage" if "Stage" in combined.obs.columns else "stage",
            cell_type_col="CellType_ID" if "CellType_ID" in combined.obs.columns else "cell_type",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE200996: {exc}")


# =========================================================================
# 13. CELLxGENE_7b20c613 (Melanoma — Multi-Cohort ICB Meta-Atlas)
# =========================================================================
def load_cellxgene_7b20c613_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load CELLxGENE 7b20c613 Melanoma meta-atlas (~355k cells, 167 patients)."""
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"CELLxGENE 7b20c613 H5AD not found in {raw_dir}")

    try:
        adata = ad.read_h5ad(h5ad_files[0])
        adata = _restore_raw_anndata(adata)

        # Standardize patient identifier to prevent inter-study collisions
        if "PMID_donor_id" in adata.obs.columns:
            adata.obs["patient_id"] = adata.obs["PMID_donor_id"].astype(str)
        elif "donor_id" in adata.obs.columns:
            adata.obs["patient_id"] = adata.obs["donor_id"].astype(str)

        adata = _standardize_tier0_obs(
            adata,
            accession="CELLxGENE_7b20c613",
            indication="Melanoma",
            technology="10x Chromium 3' v3/v3.1",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="patient_id",
            response_col="Combined_outcome" if "Combined_outcome" in adata.obs.columns else "outcome",
            cell_type_col="cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load CELLxGENE 7b20c613: {exc}")


# =========================================================================
# 14. GSE218429 (Melanoma — KEAP1/NRF2 Modulated Resistance)
# =========================================================================
def load_gse218429_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE218429 Melanoma KEAP1 atlas (~70k cells, 35 patients)."""
    count_p = raw_dir / "GSE218429_counts.csv.gz"
    if not count_p.exists():
        return Failure(f"GSE218429 counts CSV not found at {count_p}")

    try:
        df = pl.read_csv(count_p)
        first_col = df.columns[0]
        genes = df[first_col].cast(pl.String).to_list()
        barcodes = df.columns[1:]
        mat = df.select(barcodes).to_numpy().astype(np.float32).T
        sparse_mat = sp.csr_matrix(mat, dtype=np.float32)

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes
        adata = ad.AnnData(X=sparse_mat, obs=pd.DataFrame(index=barcodes), var=var_df)

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE218429",
            indication="Melanoma",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE218429: {exc}")


# =========================================================================
# 15. GSE344166 (Melanoma — Acral Melanoma Anti-PD-1 Cohort)
# =========================================================================
def load_gse344166_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE344166 Acral Melanoma atlas (~232k cells, 116 patients)."""
    norm_p = raw_dir / "GSE344166_Acral_GeoMX_norm.csv.gz"
    sum_p = raw_dir / "GSE344166_Experiment_Summary.csv.gz"

    if not norm_p.exists():
        return Failure(f"GSE344166 matrix not found at {norm_p}")

    try:
        df = pl.read_csv(norm_p)
        first_col = df.columns[0]
        genes = df[first_col].cast(pl.String).to_list()
        barcodes = df.columns[1:]
        mat = df.select(barcodes).to_numpy().astype(np.float32).T
        sparse_mat = sp.csr_matrix(mat, dtype=np.float32)

        obs_df = pd.DataFrame(index=barcodes)
        if sum_p.exists():
            try:
                summary = pl.read_csv(sum_p).to_pandas()
                s_first = summary.columns[0]
                summary = summary.set_index(s_first)
                obs_df = obs_df.join(summary, how="left")
            except Exception:
                pass

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes
        adata = ad.AnnData(X=sparse_mat, obs=obs_df, var=var_df)

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE344166",
            indication="Melanoma",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="patient" if "patient" in adata.obs.columns else "Sample_ID",
            response_col="response" if "response" in adata.obs.columns else "Best_Response",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE344166: {exc}")


def _stream_genes_cells_tsv_to_csr(
    tsv_gz_path: Path,
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Stream a large genes x cells gzipped TSV into a cells x genes sparse CSR matrix with minimal RAM."""
    with gzip.open(tsv_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
        barcodes = [b.strip() for b in header[1:] if b.strip()]
        n_cells = len(barcodes)

        row_indices: list[np.ndarray] = []
        col_indices: list[np.ndarray] = []
        data_values: list[np.ndarray] = []
        genes: list[str] = []

        for g_idx, line in enumerate(f):
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split("\t", 1)
            genes.append(parts[0])
            if len(parts) > 1 and parts[1]:
                vals = np.fromstring(parts[1], dtype=np.float32, sep="\t")
                nz = np.nonzero(vals)[0]
                if len(nz) > 0:
                    row_indices.append(nz.astype(np.int32))
                    col_indices.append(np.full(len(nz), g_idx, dtype=np.int32))
                    data_values.append(vals[nz])

        if not row_indices:
            empty_mat = sp.csr_matrix((n_cells, len(genes)), dtype=np.float32)
            return empty_mat, barcodes, genes

        all_rows = np.concatenate(row_indices)
        all_cols = np.concatenate(col_indices)
        all_data = np.concatenate(data_values)

        coo = sp.coo_matrix(
            (all_data, (all_rows, all_cols)),
            shape=(n_cells, len(genes)),
            dtype=np.float32,
        )
        return coo.tocsr(), barcodes, genes


# =========================================================================
# 16. GSE207422 (NSCLC — Neoadjuvant Chemo-Immunotherapy Remodeling)
# =========================================================================
def load_gse207422_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE207422 NSCLC Chemo-IO atlas (~78k cells, 39 patients) with verified Pathologic Response."""
    umi_p = raw_dir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_p = raw_dir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"

    if not umi_p.exists():
        return Failure(f"GSE207422 UMI matrix not found at {umi_p}")

    try:
        sparse_mat, barcodes, genes = _stream_genes_cells_tsv_to_csr(umi_p)

        sample_ids = [b.rsplit("_", 1)[0] for b in barcodes]
        obs_df = pd.DataFrame({"sample_id": sample_ids}, index=barcodes)
        if meta_p.exists():
            try:
                meta = pd.read_excel(meta_p)
                # Strip whitespace on column headers
                meta.columns = [str(c).strip() for c in meta.columns]
                sample_col = next((c for c in meta.columns if c.lower() == "sample"), meta.columns[0])
                meta = meta.set_index(sample_col)
                obs_df = obs_df.join(meta, on="sample_id", how="left")
            except Exception as exc:
                logger.warning("Could not read Excel metadata for GSE207422: %s", exc)

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes
        adata = ad.AnnData(X=sparse_mat, obs=obs_df, var=var_df)

        resp_col = "Pathologic Response" if "Pathologic Response" in adata.obs.columns else ("RECIST" if "RECIST" in adata.obs.columns else "response")
        adata = _standardize_tier0_obs(
            adata,
            accession="GSE207422",
            indication="NSCLC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="Patient" if "Patient" in adata.obs.columns else "sample_id",
            response_col=resp_col,
            treatment_col="Chemotherapy" if "Chemotherapy" in adata.obs.columns else "treatment",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE207422: {exc}")


# =========================================================================
# 17. GSE317309 (NSCLC — Pembrolizumab + CCL21-DC Vaccine Phase I Trial)
# =========================================================================
def load_gse317309_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE317309 NSCLC Phase I vaccine trial (~128k cells, 64 patients)."""
    tar_path = raw_dir / "GSE317309_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse317309"

    if not tar_path.exists():
        return Failure(f"GSE317309 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "")
            bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
            feat_p = extract_dir / f"{prefix}_features.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            # Check feature type column if 3+ cols (e.g. Gene Expression vs Antibody Capture)
            if len(feats.columns) >= 3:
                rna_mask = feats[2].astype(str) == "Gene Expression"
                if np.any(rna_mask):
                    feats = feats[rna_mask].reset_index(drop=True)
                    mat = mat[:, rna_mask.values]

            # Pre-filter empty droplets to prevent 41 million cell RAM explosion
            counts_per_cell = np.asarray(mat.sum(axis=1)).ravel()
            genes_per_cell = np.asarray((mat > 0).sum(axis=1)).ravel()
            valid_cells = (counts_per_cell >= 100) & (genes_per_cell >= 50)
            if not np.any(valid_cells):
                continue
            mat = mat[valid_cells]
            bcs = [bcs[i] for i, v in enumerate(valid_cells) if v]

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            sample_id = prefix.split("_")[0]
            sub_a = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_id
            sub_a.obs["patient_id"] = sample_id
            adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable sample matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE317309",
            indication="NSCLC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE317309: {exc}")


# =========================================================================
# 18. GSE243013 (NSCLC — Immune Heterogeneity & Treatment Response)
# =========================================================================
def load_gse243013_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE243013 NSCLC immune atlas (~486k cells, 243 patients) with verified pathological response."""
    mtx_p = raw_dir / "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
    bc_p = raw_dir / "GSE243013_barcodes.csv.gz"
    gene_p = raw_dir / "GSE243013_genes.csv.gz"
    meta_p = raw_dir / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"

    if not mtx_p.exists():
        return Failure(f"GSE243013 count matrix not found at {mtx_p}")

    try:
        if bc_p.exists():
            df_bc = pl.read_csv(bc_p, has_header=True)
            barcodes = df_bc[df_bc.columns[0]].cast(pl.String).to_list()
        else:
            barcodes = None

        if gene_p.exists():
            df_genes = pl.read_csv(gene_p, has_header=True)
            genes = df_genes[df_genes.columns[0]].cast(pl.String).to_list()
        else:
            genes = None

        try:
            import fast_matrix_market as fmm
            fmm_available = True
        except ImportError:
            fmm_available = False

        if fmm_available:
            with gzip.open(mtx_p, "rb") as f_gz:
                raw_mat = fmm.mmread(f_gz)
        else:
            raw_mat = sio.mmread(mtx_p)

        n_cells = len(barcodes) if barcodes is not None else (raw_mat.shape[0] if raw_mat.shape[0] > raw_mat.shape[1] else raw_mat.shape[1])
        n_feats = len(genes) if genes is not None else (raw_mat.shape[1] if raw_mat.shape[0] > raw_mat.shape[1] else raw_mat.shape[0])

        if barcodes is None:
            barcodes = [f"cell_{i}" for i in range(n_cells)]
        if genes is None:
            genes = [f"gene_{i}" for i in range(n_feats)]

        if sp.issparse(raw_mat):
            if raw_mat.shape[0] == len(barcodes):
                mat = raw_mat.tocsr().astype(np.float32)
            elif raw_mat.shape[1] == len(barcodes):
                mat = raw_mat.T.tocsr().astype(np.float32)
            elif raw_mat.shape[0] == len(genes):
                mat = raw_mat.T.tocsr().astype(np.float32)
            else:
                mat = raw_mat.tocsr().astype(np.float32)
        else:
            arr = np.asarray(raw_mat, dtype=np.float32)
            if arr.shape[0] == len(barcodes):
                mat = sp.csr_matrix(arr)
            else:
                mat = sp.csr_matrix(arr.T)
        del raw_mat

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes
        obs_df = pd.DataFrame(index=barcodes)

        if meta_p.exists():
            try:
                meta = pl.read_csv(meta_p, infer_schema_length=0).to_pandas()
                join_col = "cellID" if "cellID" in meta.columns else ("barcode" if "barcode" in meta.columns else meta.columns[0])
                meta = meta.set_index(join_col)
                obs_df = obs_df.join(meta, how="left")
            except Exception as exc:
                logger.warning("Could not join metadata for GSE243013: %s", exc)

        adata = ad.AnnData(X=mat, obs=obs_df, var=var_df)
        del mat, obs_df, var_df
        gc.collect()
        adata = _standardize_tier0_obs(
            adata,
            accession="GSE243013",
            indication="NSCLC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="sampleID" if "sampleID" in adata.obs.columns else "patient_id",
            sample_col="sampleID" if "sampleID" in adata.obs.columns else "sample_id",
            response_col="pathological_response" if "pathological_response" in adata.obs.columns else "radiological_response",
            cell_type_col="major_cell_type" if "major_cell_type" in adata.obs.columns else "cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE243013: {exc}")


# =========================================================================
# 19. GSE311789 (PDAC — Fibroblast DeCAF Immunotherapy Response)
# =========================================================================
def load_gse311789_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE311789 PDAC DeCAF atlas (~284k cells, 142 patients)."""
    tar_path = raw_dir / "GSE311789_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse311789"

    if not tar_path.exists():
        return Failure(f"GSE311789 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "")
            bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
            feat_p = extract_dir / f"{prefix}_features.tsv.gz"
            if not (bc_p.exists() and feat_p.exists()):
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            sample_id = prefix.split("_")[0]
            sub_a = ad.AnnData(
                X=mat.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                var=var_df,
            )
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_id
            sub_a.obs["patient_id"] = sample_id
            adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable sample matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE311789",
            indication="PDAC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE311789: {exc}")


# =========================================================================
# 20. GSE316195 (PDAC — Phase 1 Clinical Trial Single-Nucleus Atlas)
# =========================================================================
def load_gse316195_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE316195 PDAC single-nucleus atlas (~44k nuclei, 22 patients) with verified response."""
    tar_path = raw_dir / "GSE316195_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse316195"

    ext_res = _extract_tar_safely(tar_path, extract_dir)
    match ext_res:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    load_res = _load_10x_triplets_or_h5(extract_dir)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(adatas):
            combined = ad.concat(adatas, join="outer")

    try:
        # Join NCBI GEO Series Matrix response characteristics
        series_res = _parse_geo_series_matrix(raw_dir, "GSE316195")
        if isinstance(series_res, Success):
            meta_df = series_res.unwrap()
            resp_mapping: dict[str, str] = {}
            for gsm, row in meta_df.iterrows():
                title = str(row.get("Sample_title", "")).strip()
                resp = str(row.get("response", "")).strip()
                if resp:
                    resp_mapping[str(gsm)] = resp
                    if title:
                        resp_mapping[title] = resp

            # Match sample_id or patient_id against mapping
            for col in ["sample_id", "patient_id"]:
                if col in combined.obs.columns:
                    mapped = combined.obs[col].map(resp_mapping)
                    if mapped.notna().any():
                        combined.obs["clinical_response_raw"] = mapped.fillna("NA")
                        break

        combined.uns["is_single_nucleus"] = True
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE316195",
            indication="PDAC",
            technology="Single-Nucleus RNA-seq",
            cell_selection="Nuclei Isolation (snRNA-seq)",
            sample_col="sample_id",
            patient_col="patient_id",
            response_col="clinical_response_raw" if "clinical_response_raw" in combined.obs.columns else None,
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE316195: {exc}")


# =========================================================================
# 21. GSE210038 (ccRCC — Mesenchymal Tumor Cells & Nivolumab)
# =========================================================================
def load_gse210038_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE210038 ccRCC Nivolumab trial (~18k cells, 9 patients)."""
    tar_path = raw_dir / "GSE210038_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse210038"

    if not tar_path.exists():
        return Failure(f"GSE210038 RAW.tar not found in {raw_dir}")

    try:
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(tar_path, "r:*") as archive:
                archive.extractall(extract_dir)

        tsv_files = sorted(list(extract_dir.glob("*raw_counts.tsv.gz"))) or sorted(list(extract_dir.glob("*.tsv.gz")))
        adatas: list[ad.AnnData] = []

        if tsv_files:
            for tp in tsv_files:
                sample_name = tp.name.split(".")[0].replace("_raw_counts", "")
                pldf = pl.read_csv(tp, separator="\t")
                gene_col = pldf.columns[0]
                genes = pldf[gene_col].to_list()
                cell_ids = pldf.columns[1:]
                mat_data = pldf.select(cell_ids).to_numpy().T
                mat = sp.csr_matrix(mat_data, dtype=np.float32)
                del pldf, mat_data

                sub_a = ad.AnnData(
                    X=mat,
                    obs=pd.DataFrame(index=[f"{sample_name}_{cid}" for cid in cell_ids]),
                    var=pd.DataFrame(index=genes),
                )
                sub_a.var_names_make_unique()
                sub_a.obs["sample_id"] = sample_name
                sub_a.obs["patient_id"] = sample_name.split("_")[1] if "_" in sample_name else sample_name
                adatas.append(sub_a)
        else:
            mtx_files = sorted(list(extract_dir.glob("*matrix.mtx.gz")))
            for mp in mtx_files:
                prefix = mp.name.replace("_matrix.mtx.gz", "")
                bc_p = extract_dir / f"{prefix}_barcodes.tsv.gz"
                feat_p = extract_dir / f"{prefix}_features.tsv.gz"
                if not (bc_p.exists() and feat_p.exists()):
                    continue

                mat = _to_csr(sio.mmread(mp))
                bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
                feats = pd.read_csv(feat_p, header=None, sep="\t")

                var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
                var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

                sample_id = prefix.split("_")[0]
                sub_a = ad.AnnData(
                    X=mat.astype(np.float32),
                    obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in bcs]),
                    var=var_df,
                )
                sub_a.var_names_make_unique()
                sub_a.obs["sample_id"] = sample_id
                sub_a.obs["patient_id"] = sample_id
                adatas.append(sub_a)

        if not adatas:
            return Failure(f"No parseable sample matrices in {extract_dir}")

        combined = ad.concat(adatas, join="outer")
        combined = _standardize_tier0_obs(
            combined,
            accession="GSE210038",
            indication="ccRCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE210038: {exc}")


# =========================================================================
# 22. GSE314072 (ccRCC — Intratumoral CD4+CD8+ T-Cell States)
# =========================================================================
def load_gse314072_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE314072 ccRCC CD4/CD8 multiplexed atlas (~48k cells, 24 patients)."""
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"GSE314072 H5AD not found in {raw_dir}")

    try:
        adata = ad.read_h5ad(h5ad_files[0])
        adata = _restore_raw_anndata(adata)

        adata = _standardize_tier0_obs(
            adata,
            accession="GSE314072",
            indication="ccRCC",
            technology="High-Throughput scRNA-seq",
            cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
            patient_col="patient" if "patient" in adata.obs.columns else "donor",
            response_col="response" if "response" in adata.obs.columns else "Response",
            cell_type_col="cell_type" if "cell_type" in adata.obs.columns else "CellType",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE314072: {exc}")


# =========================================================================
# 23. CELLxGENE_05a8c945 (CRC — Global Core Atlas, Marteau et al. 2026)
# =========================================================================
def load_cellxgene_05a8c945_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load CELLxGENE 05a8c945 Colorectal Cancer Atlas (Marteau et al., Cancer Cell 2026).

    Uses memory-safe backed filtering to load only ICB-treated patients with response annotations.
    """
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"CELLxGENE 05a8c945 H5AD not found in {raw_dir}")

    try:
        backed_adata = ad.read_h5ad(h5ad_files[0], backed="r")
        obs_all = backed_adata.obs
        has_resp = (
            obs_all["treatment_response"].notna()
            & (obs_all["treatment_response"] != "")
            & (obs_all["treatment_response"] != "nan")
            if "treatment_response" in obs_all.columns
            else pd.Series(False, index=obs_all.index)
        )
        has_recist = (
            obs_all["RECIST"].notna()
            & (obs_all["RECIST"] != "")
            & (obs_all["RECIST"] != "nan")
            if "RECIST" in obs_all.columns
            else pd.Series(False, index=obs_all.index)
        )
        resp_mask = (has_resp | has_recist).values
        if not np.any(resp_mask):
            resp_mask = np.ones(len(obs_all), dtype=bool)

        # Pre-filter cells using existing QC covariates to avoid allocating dead droplets into RAM
        if "total_counts" in obs_all.columns and "n_genes_by_counts" in obs_all.columns and "pct_counts_mito" in obs_all.columns:
            tc = obs_all["total_counts"].to_numpy()
            ng = obs_all["n_genes_by_counts"].to_numpy()
            pm = obs_all["pct_counts_mito"].to_numpy()
            qc_pass = (tc >= 500) & (ng >= 200) & (ng <= 9000) & (pm <= 20.0)
            mask = resp_mask & qc_pass
            if not np.any(mask):
                mask = resp_mask
        else:
            mask = resp_mask

        obs = obs_all.iloc[mask].copy()
        if "pct_counts_mito" in obs.columns and "pct_counts_mt" not in obs.columns:
            obs["pct_counts_mt"] = obs["pct_counts_mito"]

        if backed_adata.raw is not None and backed_adata.raw.X is not None:
            var = backed_adata.raw.var.copy() if hasattr(backed_adata.raw, "var") and len(backed_adata.raw.var) > 0 else backed_adata.var.copy()
            raw_X = backed_adata.raw.X[mask]
        else:
            var = backed_adata.var.copy()
            raw_X = backed_adata.X[mask]
        backed_adata.file.close()

        adata = ad.AnnData(X=raw_X, obs=obs, var=var)
        del raw_X, obs, var, obs_all
        gc.collect()

        resp_col = "treatment_response" if "treatment_response" in adata.obs.columns else ("RECIST" if "RECIST" in adata.obs.columns else None)
        adata = _standardize_tier0_obs(
            adata,
            accession="CELLxGENE_05a8c945",
            indication="Colorectal Cancer",
            technology="10x Chromium 3' / 5'",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="donor_id" if "donor_id" in adata.obs.columns else "patient_id",
            response_col=resp_col,
            treatment_col="treatment_status_before_resection" if "treatment_status_before_resection" in adata.obs.columns else "treatment_status",
            cell_type_col="cell_type" if "cell_type" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load CELLxGENE 05a8c945: {exc}")


# =========================================================================
# 24. CELLxGENE_6f9de485 (TNBC — ARTEMIS Trial, Navin et al. Nature 2026)
# =========================================================================
def load_cellxgene_6f9de485_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load CELLxGENE 6f9de485 Triple-Negative Breast Cancer ARTEMIS trial (Navin et al., Nature 2026)."""
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"CELLxGENE 6f9de485 H5AD not found in {raw_dir}")

    try:
        backed_adata = ad.read_h5ad(h5ad_files[0], backed="r")
        obs = backed_adata.obs.copy()
        if backed_adata.raw is not None and backed_adata.raw.X is not None:
            var = backed_adata.raw.var.copy() if hasattr(backed_adata.raw, "var") and len(backed_adata.raw.var) > 0 else backed_adata.var.copy()
            raw_X = backed_adata.raw.X[:]
        else:
            var = backed_adata.var.copy()
            raw_X = backed_adata.X[:]
        backed_adata.file.close()

        adata = ad.AnnData(X=raw_X, obs=obs, var=var)
        del raw_X, obs, var
        gc.collect()

        adata = _standardize_tier0_obs(
            adata,
            accession="CELLxGENE_6f9de485",
            indication="Breast Cancer",
            technology="10x Chromium 3' v3",
            cell_selection="Unselected / Total Single-Cell Suspension",
            patient_col="donor_id" if "donor_id" in adata.obs.columns else "patient_id",
            response_col="pCR_status" if "pCR_status" in adata.obs.columns else "response",
            cell_type_col="cell_type" if "cell_type" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load CELLxGENE 6f9de485: {exc}")


# =========================================================================
# 25. GSE233203 (NSCLC — Combination ABCP Immunotherapy Trial, Hong et al. 2025)
# =========================================================================
def load_gse233203_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load GSE233203 NSCLC Combination ABCP Immunotherapy Trial (Hong et al., MedComm 2025)."""
    tar_path = raw_dir / "GSE233203_RAW.tar"
    extract_dir = raw_dir / "unpacked_gse233203"

    ext_res = _extract_tar_safely(tar_path, extract_dir)
    match ext_res:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    load_res = _load_10x_triplets_or_h5(extract_dir)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(adatas):
            combined = ad.concat(adatas, join="outer")

    try:
        # Join NCBI GEO Series Matrix response characteristics
        series_res = _parse_geo_series_matrix(raw_dir, "GSE233203")
        if isinstance(series_res, Success):
            meta_df = series_res.unwrap()
            resp_mapping: dict[str, str] = {}
            for gsm, row in meta_df.iterrows():
                title = str(row.get("Sample_title", "")).strip()
                resp = str(row.get("therapeutic response", "")).strip()
                if resp:
                    resp_mapping[str(gsm)] = resp
                    if title:
                        resp_mapping[title] = resp

            # Match sample_id or patient_id against mapping
            for col in ["sample_id", "patient_id"]:
                if col in combined.obs.columns:
                    mapped = combined.obs[col].map(resp_mapping)
                    if mapped.notna().any():
                        combined.obs["clinical_response_raw"] = mapped.fillna("NA")
                        break

        combined = _standardize_tier0_obs(
            combined,
            accession="GSE233203",
            indication="NSCLC",
            technology="10x Chromium 5'",
            cell_selection="Unselected / Total Single-Cell Suspension",
            sample_col="sample_id",
            patient_col="patient_id",
            response_col="clinical_response_raw" if "clinical_response_raw" in combined.obs.columns else None,
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to load GSE233203: {exc}")


__all__ = [
    "load_gse246613_breast",
    "load_gse300475_breast",
    "load_gse212707_breast",
    "load_gse236581_crc",
    "load_gse299651_crc",
    "load_cellxgene_829a3cd1_crc",
    "load_gse270680_gastric",
    "load_gse313642_hcc",
    "load_gse245906_hcc",
    "load_gse301741_hnscc",
    "load_gse287301_hnscc",
    "load_gse200996_hnscc",
    "load_cellxgene_7b20c613_melanoma",
    "load_cellxgene_05a8c945_crc",
    "load_cellxgene_6f9de485_breast",
    "load_gse218429_melanoma",
    "load_gse344166_melanoma",
    "load_gse207422_nsclc",
    "load_gse317309_nsclc",
    "load_gse243013_nsclc",
    "load_gse233203_nsclc",
    "load_gse311789_pdac",
    "load_gse316195_pdac",
    "load_gse210038_ccrcc",
    "load_gse314072_ccrcc",
]
