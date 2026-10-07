"""Pure functional projection of neighborhood differential abundance metrics to single cells.
"""

from __future__ import annotations

from typing import Any, Sequence
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse as sp


def project_nhoods_to_cells(
    nhoods_mat: sp.spmatrix,
    res_df: pd.DataFrame | pl.DataFrame,
    cell_ids: Sequence[str],
    patient_ids: Sequence[str],
    clinical_responses: Sequence[str],
    cell_types: Sequence[str] | None = None,
    treatment_statuses: Sequence[str] | None = None,
    umap_coords: np.ndarray | None = None,
    fdr_threshold: float = 0.10,
) -> pl.DataFrame:
    """Projects neighborhood log2FC and significance status back to single cells via sparse matrix multiplication.

    Args:
        nhoods_mat: Sparse matrix of shape (n_cells, n_nhoods).
        res_df: DataFrame with neighborhood test results (must contain 'logFC' and either 'is_significant' or 'FDR').
        cell_ids: Sequence of cell identifiers.
        patient_ids: Sequence of donor/patient identifiers.
        clinical_responses: Sequence of response strings ('responder', 'non-responder').
        cell_types: Optional sequence of cell type annotations.
        treatment_statuses: Optional sequence of biopsy timepoints ('pre-treatment', etc.).
        umap_coords: Optional (n_cells, 2) array of UMAP coordinates.
        fdr_threshold: Nominal FDR significance threshold.

    Returns:
        polars.DataFrame containing cell-level scores and mapped DA statuses.
    """
    if not sp.isspmatrix_csc(nhoods_mat):
        nhoods_mat = sp.csc_matrix(nhoods_mat)

    n_cells, n_nhoods = nhoods_mat.shape

    # Extract logFC vector
    if isinstance(res_df, pl.DataFrame):
        logfc_vec = np.asarray(res_df["logFC"].fill_null(0.0).to_numpy(), dtype=np.float64)
    else:
        logfc_vec = np.asarray(res_df["logFC"].fillna(0.0).values, dtype=np.float64)

    # 1. Project continuous log2FC
    cell_nhood_counts = np.array(nhoods_mat.sum(axis=1)).flatten()
    cell_logfc_sum = nhoods_mat.dot(logfc_vec)

    cell_logfc = np.zeros(n_cells, dtype=np.float64)
    valid_mask = cell_nhood_counts > 0
    cell_logfc[valid_mask] = cell_logfc_sum[valid_mask] / cell_nhood_counts[valid_mask]

    # 2. Determine discrete DA status (DA+, DA-, Not Significant)
    v_resp = np.zeros(n_nhoods, dtype=np.float32)
    v_non_resp = np.zeros(n_nhoods, dtype=np.float32)

    if isinstance(res_df, pl.DataFrame):
        has_is_sig = "is_significant" in res_df.columns
        if has_is_sig:
            sig_r_idx = res_df.filter((pl.col("is_significant")) & (pl.col("logFC") > 0)).select(pl.first()).to_series().to_list()
            sig_nr_idx = res_df.filter((pl.col("is_significant")) & (pl.col("logFC") < 0)).select(pl.first()).to_series().to_list()
        else:
            sig_r_idx = res_df.filter((pl.col("FDR") < fdr_threshold) & (pl.col("logFC") > 0)).select(pl.first()).to_series().to_list()
            sig_nr_idx = res_df.filter((pl.col("FDR") < fdr_threshold) & (pl.col("logFC") < 0)).select(pl.first()).to_series().to_list()
    else:
        has_is_sig = "is_significant" in res_df.columns
        if has_is_sig:
            sig_r_idx = list(np.where(res_df["is_significant"] & (res_df["logFC"] > 0))[0])
            sig_nr_idx = list(np.where(res_df["is_significant"] & (res_df["logFC"] < 0))[0])
        else:
            sig_r_idx = list(np.where((res_df["FDR"] < fdr_threshold) & (res_df["logFC"] > 0))[0])
            sig_nr_idx = list(np.where((res_df["FDR"] < fdr_threshold) & (res_df["logFC"] < 0))[0])

    if sig_r_idx:
        v_resp[sig_r_idx] = 1.0
    if sig_nr_idx:
        v_non_resp[sig_nr_idx] = 1.0

    c_resp = np.asarray(nhoods_mat.dot(v_resp)).flatten()
    c_non_resp = np.asarray(nhoods_mat.dot(v_non_resp)).flatten()

    in_resp = c_resp > 0
    in_non_resp = c_non_resp > 0

    da_status: list[str] = []
    for i in range(n_cells):
        ir = in_resp[i]
        inr = in_non_resp[i]
        if ir and not inr:
            da_status.append("DA+ (Responder Enriched)")
        elif inr and not ir:
            da_status.append("DA- (Non-Responder Enriched)")
        elif ir and inr:
            # Overlapping opposite neighborhoods: assign by net sign of cell_logfc
            da_status.append("DA+ (Responder Enriched)" if cell_logfc[i] > 0 else "DA- (Non-Responder Enriched)")
        else:
            da_status.append("Not Significant")

    # Build Polars DataFrame
    cell_dict: dict[str, Any] = {
        "cell_id": list(cell_ids),
        "patient_id": list(patient_ids),
        "clinical_response": list(clinical_responses),
        "treatment_status": list(treatment_statuses) if treatment_statuses is not None else ["Unknown"] * n_cells,
        "cell_type": list(cell_types) if cell_types is not None else ["Unknown"] * n_cells,
        "cell_logfc": [float(x) for x in cell_logfc],
        "nhood_count": [int(x) for x in cell_nhood_counts],
        "da_status": da_status,
    }

    if umap_coords is not None and umap_coords.shape[0] == n_cells:
        cell_dict["UMAP1"] = [float(x) for x in umap_coords[:, 0]]
        cell_dict["UMAP2"] = [float(x) for x in umap_coords[:, 1]]

    return pl.DataFrame(cell_dict)
