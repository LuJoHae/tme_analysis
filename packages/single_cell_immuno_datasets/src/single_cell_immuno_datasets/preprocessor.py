"""Functional preprocessing, metadata harmonization, and QC pipeline for ICB single-cell datasets."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from returns.result import Result, Success, Failure
from gene_utils import norm_genes
from .config import DatasetSpec, QualityControlSpec, DataDirectories


def check_transform_state(adata: ad.AnnData) -> tuple[ad.AnnData, dict[str, bool | float]]:
    """Detects if expression matrix is log-transformed and/or total-sum normalized.

    Returns updated AnnData and detection dictionary.
    """
    has_uns_log = "log1p" in adata.uns or "log_transformed" in adata.uns
    has_uns_norm = "normalized" in adata.uns
    
    X_mat = adata.X
    if hasattr(X_mat, "toarray"):
        sample_vals = X_mat[:min(100, adata.n_obs)].toarray()
    else:
        sample_vals = X_mat[:min(100, adata.n_obs)]

    max_val = float(sample_vals.max()) if sample_vals.size > 0 else 0.0
    is_non_integer = not np.all(np.equal(np.mod(sample_vals, 1), 0))
    is_log = has_uns_log or (max_val <= 35.0 and is_non_integer and max_val > 0.0)

    row_sums = sample_vals.sum(axis=1) if sample_vals.ndim > 1 else sample_vals
    is_norm = has_uns_norm or np.allclose(row_sums, 1e4, rtol=1e-2) or np.allclose(row_sums, 1e6, rtol=1e-2)

    status: dict[str, bool | float] = {
        "is_log1p": is_log,
        "is_normalized": is_norm,
        "max_value": max_val,
    }
    adata.uns["transform_state"] = status

    if is_log:
        if hasattr(adata.X, "data"):
            adata.X.data = np.expm1(adata.X.data)
        else:
            adata.X = np.expm1(adata.X)

    return adata, status


def harmonize_metadata(adata: ad.AnnData, spec: DatasetSpec) -> ad.AnnData:
    """Harmonizes obs metadata columns to standard schema while retaining original fields."""
    obs_df = adata.obs

    if "original.barcode" not in obs_df.columns:
        obs_df["original.barcode"] = obs_df.index.astype(str)

    obs_df["dataset"] = spec.accession
    obs_df["cancer_code"] = spec.cancer_type
    obs_df["sequencing_tech"] = spec.sequencing_tech
    obs_df["therapy"] = spec.therapy

    # Apply specific obs column renames from DatasetSpec
    for old_col, new_col in spec.obs_mapping.items():
        if old_col in obs_df.columns and new_col not in obs_df.columns:
            obs_df[new_col] = obs_df[old_col]

    if "cell_type" not in obs_df.columns:
        for col in ["cell_type_main", "celltype", "cell_types", "CellType", "major_cell_type"]:
            if col in obs_df.columns:
                obs_df["cell_type"] = obs_df[col]
                break
        else:
            obs_df["cell_type"] = "Unknown"

    if "patient" not in obs_df.columns or obs_df["patient"].nunique() <= 1:
        for col in ["patient_id", "Patient_ID", "patient_geo", "sample", "Sample", "subject", "Subject", "samples", "Cohort", "geo_accession"]:
            if col in obs_df.columns and obs_df[col].nunique() > 1:
                obs_df["patient"] = obs_df[col]
                break
        else:
            # Try underscore/dash prefix: e.g. 'GSM3496325_AB1889', 'BCC01_barcode'
            prefixes = obs_df.index.to_series().str.extract(r"^([A-Za-z0-9]+[_-][A-Za-z0-9]+)[_-]")[0]
            if prefixes.notna().sum() > 0.3 * len(obs_df) and prefixes.nunique() > 1:
                obs_df["patient"] = prefixes.fillna("Unknown").values
            else:
                # Try dot-separated prefix: e.g. 'bcc.su001.pre.tcell_barcode' -> 'su001'
                dot_prefixes = obs_df.index.to_series().str.extract(r"^[A-Za-z]+\.([A-Za-z0-9]+)\.", expand=False)
                if dot_prefixes.notna().sum() > 0.3 * len(obs_df) and dot_prefixes.nunique() > 1:
                    obs_df["patient"] = dot_prefixes.fillna("Unknown").values
                else:
                    obs_df["patient"] = "Unknown"

    if "response" not in obs_df.columns or obs_df["response"].nunique() <= 1:
        for col in ["characteristics: response", "Response", "RECIST", "response_status", "therapy_response",
                    "ICB_response", "ICB_Response", "outcome", "Outcome", "Combined_outcome",
                    "treatment.group", "treatment_group", "Treatment_group",
                    "group", "Group", "cohort", "Cohort"]:
            if col in obs_df.columns and obs_df[col].nunique() > 1:
                obs_df["response"] = obs_df[col]
                break

        else:
            # Do NOT use patient as response — it creates N unique conditions and breaks milopy
            if "response" not in obs_df.columns:
                obs_df["response"] = "Unknown"

    # Remap 'nan' string (produced by str() of actual NaN) to 'Unknown' in response
    if "response" in obs_df.columns:
        obs_df["response"] = obs_df["response"].astype(str).replace({"nan": "Unknown", "None": "Unknown", "NaN": "Unknown"})

    # Apply hardcoded patient → response map from DatasetSpec when response is still uniform
    if spec.response_map and ("response" not in obs_df.columns or obs_df["response"].nunique() <= 1):
        if "patient" in obs_df.columns:
            obs_df["response"] = obs_df["patient"].map(spec.response_map).fillna("Unknown")

    adata.obs = obs_df
    return adata



def preprocess_anndata(
    adata: ad.AnnData,
    spec: DatasetSpec,
    qc_spec: QualityControlSpec,
    out_h5ad_path: Path,
) -> Result[Path, str]:
    """Applies transformation state check, metadata harmonization, gene normalization, QC, and writes .h5ad."""
    try:
        # 1. Detect prior transformation
        adata, status = check_transform_state(adata)

        # 2. Harmonize metadata
        adata = harmonize_metadata(adata, spec)

        # 3. Gene symbol normalization
        try:
            adata = norm_genes(adata, pre_id_transform="auto")
        except Exception:
            pass

        # 4. Compute QC metrics
        if "contig" in adata.var.columns:
            adata.var["mt"] = (adata.var["contig"] == "MT").values
        else:
            adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # 5. Filter cells in a single boolean slice (handling TPM datasets)
        max_total = adata.obs["total_counts"].max() if "total_counts" in adata.obs.columns else 0.0
        is_tpm = "tpm" in spec.sequencing_tech.lower() or "smart-seq" in spec.sequencing_tech.lower() or max_total > 100_000

        valid_counts = (
            (adata.obs["total_counts"] >= qc_spec.min_counts) if is_tpm else
            ((adata.obs["total_counts"] >= qc_spec.min_counts) & (adata.obs["total_counts"] <= qc_spec.max_counts))
        )

        # Smart-seq2 / TPM: skip max_genes upper bound (full-transcript capture → more genes detected)
        valid_genes = (
            (adata.obs["n_genes_by_counts"] >= qc_spec.min_genes) if is_tpm else
            ((adata.obs["n_genes_by_counts"] >= qc_spec.min_genes) & (adata.obs["n_genes_by_counts"] <= qc_spec.max_genes))
        )

        valid_cells = valid_genes & valid_counts

        # Skip pct_counts_mt filter for TPM/Smart-seq2: high MT% is a feature of full-transcript capture, not damage
        if "pct_counts_mt" in adata.obs.columns and not is_tpm:
            valid_cells = valid_cells & (adata.obs["pct_counts_mt"] < qc_spec.max_mt_content)

        adata = adata[valid_cells, :].copy()


        # Clean obs column data types for HDF5 string compatibility
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object or str(adata.obs[col].dtype) == "category":
                adata.obs[col] = adata.obs[col].astype(str).fillna("Unknown")

        # 6. Save preprocessed HDF5 (.h5ad)
        out_h5ad_path.parent.mkdir(parents=True, exist_ok=True)
        ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(out_h5ad_path, compression="gzip")

        return Success(out_h5ad_path)
    except Exception as e:
        return Failure(f"Preprocessing failed for dataset {spec.accession}: {str(e)}")
