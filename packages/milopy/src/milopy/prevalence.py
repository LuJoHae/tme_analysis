"""Biological Replicate Prevalence and Multi-Cohort Recurrence Filtering for Milopy.

Enforces rigorous donor prevalence gates and cross-study recurrence filters
to eliminate false discoveries driven by single-patient clonal expansions
and unintegrated private tumor variants.

Adheres strictly to functional Python, immutable Pydantic models, and returns.result.
"""

from __future__ import annotations

from typing import Sequence
import anndata as ad
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field
from returns.result import Failure, Result, Success
from scipy import sparse as sp


class ReplicatePrevalenceConfig(BaseModel):
    """Immutable configuration for biological replicate prevalence and recurrence filtering."""

    model_config = ConfigDict(frozen=True)

    min_replicates: int = Field(default=2, ge=1, description="Minimum independent biological donors in enriched arm")
    min_prevalence_frac: float = Field(default=0.05, ge=0.0, le=1.0, description="Minimum fraction of donors in enriched arm")
    min_datasets: int = Field(default=2, ge=1, description="Minimum cohorts for cross-study conserved hits")
    dataset_purity_threshold: float = Field(default=0.85, ge=0.0, le=1.0, description="Purity ceiling for conserved hits")
    fdr_threshold: float = Field(default=0.10, ge=0.0, le=1.0, description="Nominal FDR threshold")
    require_permutation_significance: bool = Field(default=False, description="Whether permutation FDR must also pass")
    perm_fdr_threshold: float = Field(default=0.15, ge=0.0, le=1.0, description="Permutation FDR threshold")


class PrevalenceSummary(BaseModel):
    """Structured summary of prevalence and recurrence filtering."""

    model_config = ConfigDict(frozen=True)

    total_nhoods: int
    n_sig_up_replicated: int
    n_sig_down_replicated: int
    n_private_spikes_blocked: int
    n_conserved_hits: int
    median_patients_per_nhood: float


def evaluate_replicate_prevalence(
    adata: ad.AnnData,
    res_df: pd.DataFrame,
    patient_col: str,
    design_col: str,
    dataset_col: str | None = None,
    permutation_pvalues: Sequence[float] | None = None,
    permutation_fdr: Sequence[float] | None = None,
    config: ReplicatePrevalenceConfig = ReplicatePrevalenceConfig(),
) -> Result[pd.DataFrame, str]:
    """Annotates neighborhoods with donor prevalence and evaluates recurrence filters.

    Args:
        adata: AnnData with neighborhood incidence matrix in `adata.obsm['nhoods']`.
        res_df: Neighborhood test results DataFrame (must contain 'logFC' and 'FDR').
        patient_col: Column name for patient / donor ID in `adata.obs`.
        design_col: Column name for clinical response in `adata.obs`.
        dataset_col: Optional column name for study / cohort ID in `adata.obs`.
        permutation_pvalues: Optional sequence of permutation empirical p-values.
        permutation_fdr: Optional sequence of permutation Benjamini-Hochberg FDR values.
        config: ReplicatePrevalenceConfig.

    Returns:
        Success(pd.DataFrame) with added prevalence columns and status, or Failure(error).
    """
    if "nhoods" not in adata.obsm:
        return Failure("Neighborhood incidence matrix 'nhoods' not found in adata.obsm")

    if patient_col not in adata.obs.columns:
        return Failure(f"Patient column '{patient_col}' not found in adata.obs")

    if design_col not in adata.obs.columns:
        return Failure(f"Design column '{design_col}' not found in adata.obs")

    df = res_df.copy()
    nhoods_mat = adata.obsm["nhoods"]
    if not sp.isspmatrix_csc(nhoods_mat):
        nhoods_mat = sp.csc_matrix(nhoods_mat)

    n_cells, n_nhoods = nhoods_mat.shape
    if len(df) != n_nhoods:
        return Failure(f"Result DataFrame has {len(df)} rows, but nhoods matrix has {n_nhoods} columns")

    patients = adata.obs[patient_col].astype(str).values
    responses = adata.obs[design_col].astype(str).str.lower().values

    # Total unique donors in cohort by response
    unique_resp_pts = set(patients[responses == "responder"])
    unique_non_resp_pts = set(patients[responses == "non-responder"])
    n_total_resp_pts = max(1, len(unique_resp_pts))
    n_total_non_resp_pts = max(1, len(unique_non_resp_pts))

    has_dataset = dataset_col is not None and dataset_col in adata.obs.columns
    datasets = adata.obs[dataset_col].astype(str).values if has_dataset else None

    # Track metrics per neighborhood
    n_pts_total: list[int] = []
    n_pts_resp: list[int] = []
    n_pts_non_resp: list[int] = []
    frac_pts_resp: list[float] = []
    frac_pts_non_resp: list[float] = []
    n_ds_total: list[int] = []

    for j in range(n_nhoods):
        idx = nhoods_mat[:, j].nonzero()[0]
        if len(idx) == 0:
            n_pts_total.append(0)
            n_pts_resp.append(0)
            n_pts_non_resp.append(0)
            frac_pts_resp.append(0.0)
            frac_pts_non_resp.append(0.0)
            if has_dataset:
                n_ds_total.append(0)
            continue

        pts_in_nh = patients[idx]
        u_pts = np.unique(pts_in_nh)
        n_pts_total.append(len(u_pts))

        resp_mask = responses[idx] == "responder"
        non_resp_mask = responses[idx] == "non-responder"

        nr_resp = len(np.unique(pts_in_nh[resp_mask]))
        nr_non_resp = len(np.unique(pts_in_nh[non_resp_mask]))

        n_pts_resp.append(nr_resp)
        n_pts_non_resp.append(nr_non_resp)
        frac_pts_resp.append(nr_resp / n_total_resp_pts)
        frac_pts_non_resp.append(nr_non_resp / n_total_non_resp_pts)

        if has_dataset and datasets is not None:
            n_ds_total.append(len(np.unique(datasets[idx])))

    df["n_patients_total"] = n_pts_total
    df["n_patients_responder"] = n_pts_resp
    df["n_patients_non_responder"] = n_pts_non_resp
    df["prevalence_frac_responder"] = frac_pts_resp
    df["prevalence_frac_non_responder"] = frac_pts_non_resp

    if has_dataset:
        df["n_datasets_total"] = n_ds_total

    # Add permutation statistics if provided
    has_perm = permutation_fdr is not None and len(permutation_fdr) == n_nhoods
    if has_perm and permutation_fdr is not None and permutation_pvalues is not None:
        df["PValue_perm"] = list(permutation_pvalues)
        df["FDR_perm"] = list(permutation_fdr)
        perm_sig_mask = np.array(permutation_fdr) < config.perm_fdr_threshold
    else:
        perm_sig_mask = np.ones(n_nhoods, dtype=bool)

    # Statistical significance under parametric GLM
    df["neg_log10_fdr"] = -np.log10(df["FDR"].clip(lower=1e-300))
    is_glm_sig = df["FDR"] < config.fdr_threshold

    # Biological Replicate Prevalence Gate:
    # Responder enriched: >= min_replicates AND >= min_prevalence_frac of all responders
    # Non-responder enriched: >= min_replicates AND >= min_prevalence_frac of all non-responders
    logfc = df["logFC"].values
    resp_prev_ok = (np.array(n_pts_resp) >= config.min_replicates) & (np.array(frac_pts_resp) >= config.min_prevalence_frac)
    non_resp_prev_ok = (np.array(n_pts_non_resp) >= config.min_replicates) & (np.array(frac_pts_non_resp) >= config.min_prevalence_frac)

    has_prevalence = ((logfc > 0) & resp_prev_ok) | ((logfc < 0) & non_resp_prev_ok)

    # Combined significance filter
    if config.require_permutation_significance:
        df["is_significant"] = is_glm_sig & has_prevalence & perm_sig_mask
    else:
        df["is_significant"] = is_glm_sig & has_prevalence

    # Multi-cohort recurrence categorization
    is_multi_dataset = has_dataset and "dataset_id_purity" in df.columns
    status_list: list[str] = []
    hit_type_list: list[str] = []

    for j in range(n_nhoods):
        fdr_ok = is_glm_sig.iloc[j] if hasattr(is_glm_sig, "iloc") else is_glm_sig[j]
        prev_ok = has_prevalence[j]
        perm_ok = perm_sig_mask[j] if config.require_permutation_significance else True
        lfc = logfc[j]

        if fdr_ok and prev_ok and perm_ok:
            # Verified replicated hit
            base_status = "Enriched in Responders" if lfc > 0 else "Enriched in Non-Responders"
            status_list.append(base_status)

            if is_multi_dataset:
                n_ds = n_ds_total[j]
                purity = df.loc[df.index[j], "dataset_id_purity"] if "dataset_id_purity" in df.columns else 1.0
                if n_ds >= config.min_datasets and purity < config.dataset_purity_threshold:
                    hit_type_list.append("Conserved Recurrent Hit")
                else:
                    hit_type_list.append("Cohort-Private Replicated Hit")
            else:
                hit_type_list.append("Replicated Hit")
        elif fdr_ok and (not prev_ok or not perm_ok):
            # Candidate hit that collapsed due to single-donor spike or permutation test
            status_list.append("Private Clonal Spike")
            hit_type_list.append("Private Clonal Spike")
        else:
            status_list.append("Not Significant")
            hit_type_list.append("Not Significant")

    df["status"] = status_list
    df["hit_type"] = hit_type_list
    return Success(df)
