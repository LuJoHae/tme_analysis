"""Quality control calculation, adaptive thresholding, and filtering for single-cell RNA-seq."""

from __future__ import annotations

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import QualityControlSpec

logger = get_logger("preprocessing.qc")


def compute_qc_covariates(adata: ad.AnnData) -> ad.AnnData:
    """Compute cellular quality control metrics: count depth, gene detection, mito/ribo/hb percentage.

    Adds to adata.obs:
        - total_counts: Total UMI count depth per cell.
        - n_genes_by_counts: Number of non-zero detected genes.
        - pct_counts_mt: Percentage of counts belonging to mitochondrial genes.
        - pct_counts_ribo: Percentage of counts belonging to ribosomal genes.
        - pct_counts_hb: Percentage of counts belonging to hemoglobin genes.

    Args:
        adata: Input AnnData with raw counts in .X.

    Returns:
        AnnData with computed QC metrics in .obs.
    """
    new_adata = adata
    new_adata.var_names_make_unique()

    # 1. Identify mitochondrial genes
    mt_mask = new_adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
    if not np.any(mt_mask) and "gene_name" in new_adata.var.columns:
        mt_mask = new_adata.var["gene_name"].astype(str).str.startswith(("MT-", "mt-", "Mt-")).to_numpy()
    if not np.any(mt_mask) and "contig" in new_adata.var.columns:
        mt_mask = new_adata.var["contig"].astype(str).isin(["MT", "M", "chrM", "chrMT"]).to_numpy()

    # 2. Identify ribosomal genes
    ribo_mask = new_adata.var_names.str.startswith(("RPS", "RPL", "rps", "rpl"))
    if not np.any(ribo_mask) and "gene_name" in new_adata.var.columns:
        ribo_mask = new_adata.var["gene_name"].astype(str).str.startswith(("RPS", "RPL", "rps", "rpl")).to_numpy()

    # 3. Identify hemoglobin genes
    hb_mask = new_adata.var_names.str.startswith(("HBA", "HBB", "HBD", "HBE", "HBG", "hba", "hbb"))
    if not np.any(hb_mask) and "gene_name" in new_adata.var.columns:
        hb_mask = new_adata.var["gene_name"].astype(str).str.startswith(("HBA", "HBB", "HBD", "HBE", "HBG")).to_numpy()

    X = new_adata.X
    if sp.isspmatrix_csr(X):
        total_counts = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        genes_detected = np.diff(X.indptr).astype(np.int64)

        if np.any(mt_mask) or np.any(ribo_mask) or np.any(hb_mask):
            row_idx = np.repeat(np.arange(new_adata.n_obs), np.diff(X.indptr))

            if np.any(mt_mask):
                is_mt = np.zeros(new_adata.n_vars, dtype=bool)
                is_mt[mt_mask] = True
                m = is_mt[X.indices]
                mt_counts = np.bincount(row_idx[m], weights=X.data[m], minlength=new_adata.n_obs)
            else:
                mt_counts = np.zeros(new_adata.n_obs, dtype=np.float64)

            if np.any(ribo_mask):
                is_rb = np.zeros(new_adata.n_vars, dtype=bool)
                is_rb[ribo_mask] = True
                m = is_rb[X.indices]
                ribo_counts = np.bincount(row_idx[m], weights=X.data[m], minlength=new_adata.n_obs)
            else:
                ribo_counts = np.zeros(new_adata.n_obs, dtype=np.float64)

            if np.any(hb_mask):
                is_hb = np.zeros(new_adata.n_vars, dtype=bool)
                is_hb[hb_mask] = True
                m = is_hb[X.indices]
                hb_counts = np.bincount(row_idx[m], weights=X.data[m], minlength=new_adata.n_obs)
            else:
                hb_counts = np.zeros(new_adata.n_obs, dtype=np.float64)
        else:
            mt_counts = np.zeros(new_adata.n_obs, dtype=np.float64)
            ribo_counts = np.zeros(new_adata.n_obs, dtype=np.float64)
            hb_counts = np.zeros(new_adata.n_obs, dtype=np.float64)
    elif sp.issparse(X):
        total_counts = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        genes_detected = np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int64)
        mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).ravel() if np.any(mt_mask) else np.zeros(new_adata.n_obs)
        ribo_counts = np.asarray(X[:, ribo_mask].sum(axis=1)).ravel() if np.any(ribo_mask) else np.zeros(new_adata.n_obs)
        hb_counts = np.asarray(X[:, hb_mask].sum(axis=1)).ravel() if np.any(hb_mask) else np.zeros(new_adata.n_obs)
    else:
        arr = np.asarray(X, dtype=np.float64)
        total_counts = np.sum(arr, axis=1)
        genes_detected = np.sum(arr > 0, axis=1)

        mt_counts = np.sum(arr[:, mt_mask], axis=1) if np.any(mt_mask) else np.zeros(new_adata.n_obs)
        ribo_counts = np.sum(arr[:, ribo_mask], axis=1) if np.any(ribo_mask) else np.zeros(new_adata.n_obs)
        hb_counts = np.sum(arr[:, hb_mask], axis=1) if np.any(hb_mask) else np.zeros(new_adata.n_obs)

    safe_total = np.maximum(total_counts, 1.0)
    pct_mt = (mt_counts / safe_total) * 100.0
    pct_ribo = (ribo_counts / safe_total) * 100.0
    pct_hb = (hb_counts / safe_total) * 100.0

    new_adata.obs["total_counts"] = total_counts
    new_adata.obs["n_genes_by_counts"] = genes_detected
    new_adata.obs["pct_counts_mt"] = pct_mt
    new_adata.obs["pct_counts_ribo"] = pct_ribo
    new_adata.obs["pct_counts_hb"] = pct_hb

    # Annotate .var with feature classification
    new_adata.var["is_mitochondrial"] = mt_mask
    new_adata.var["is_ribosomal"] = ribo_mask
    new_adata.var["is_hemoglobin"] = hb_mask

    return new_adata


def calculate_adaptive_thresholds(
    adata: ad.AnnData,
    n_mads: float = 3.0,
    base_spec: QualityControlSpec | None = None,
) -> QualityControlSpec:
    """Determine data-driven outlier filtering thresholds using Median Absolute Deviation (MAD).

    Calculates:
        y = log10(metric)
        MAD = median(|y - median(y)|) * 1.4826
        lower = 10^(median(y) - n_mads * MAD)
        upper = 10^(median(y) + n_mads * MAD)

    Boundaries are safeguarded with sensible biological floor/ceiling limits.

    Args:
        adata: AnnData with total_counts, n_genes_by_counts, and pct_counts_mt in .obs.
        n_mads: Number of Median Absolute Deviations (standard is 3.0 or 5.0).
        base_spec: Optional QualityControlSpec to override specific fields.

    Returns:
        Resolved QualityControlSpec with adaptive thresholds.
    """
    total_counts = adata.obs["total_counts"].to_numpy()
    n_genes = adata.obs["n_genes_by_counts"].to_numpy()
    pct_mt = adata.obs["pct_counts_mt"].to_numpy()

    # Total counts MAD (log10 space)
    log_counts = np.log10(np.maximum(total_counts, 1.0))
    med_c = np.median(log_counts)
    mad_c = np.median(np.abs(log_counts - med_c)) * 1.4826
    adapt_min_counts = int(np.clip(10.0 ** (med_c - n_mads * mad_c), 200.0, 2000.0))
    adapt_max_counts = int(10.0 ** (med_c + n_mads * mad_c)) if (med_c + n_mads * mad_c) < 6.0 else None

    # Genes detected MAD (log10 space)
    log_genes = np.log10(np.maximum(n_genes, 1.0))
    med_g = np.median(log_genes)
    mad_g = np.median(np.abs(log_genes - med_g)) * 1.4826
    adapt_min_genes = int(np.clip(10.0 ** (med_g - n_mads * mad_g), 100.0, 1000.0))
    adapt_max_genes = int(np.clip(10.0 ** (med_g + n_mads * mad_g), 4000.0, 15000.0))

    # Mitochondrial percentage MAD (linear space)
    med_mt = np.median(pct_mt)
    mad_mt = np.median(np.abs(pct_mt - med_mt)) * 1.4826
    adapt_max_mt = float(np.clip(med_mt + n_mads * mad_mt, 5.0, 20.0))

    # Combine with base_spec overrides if provided
    spec = base_spec or QualityControlSpec()
    resolved_min_counts = max(spec.min_counts_per_cell, adapt_min_counts)
    resolved_min_genes = max(spec.min_genes_per_cell, adapt_min_genes)
    resolved_max_genes = min(spec.max_genes_per_cell, adapt_max_genes)
    resolved_max_mt = min(spec.max_pct_mitochondrial, adapt_max_mt)
    resolved_max_counts = spec.max_counts_per_cell if isinstance(spec.max_counts_per_cell, Some) else (
        Some(adapt_max_counts) if adapt_max_counts is not None and adapt_max_counts > resolved_min_counts else Nothing
    )

    return QualityControlSpec(
        min_counts_per_cell=resolved_min_counts,
        max_counts_per_cell=resolved_max_counts,
        min_genes_per_cell=resolved_min_genes,
        max_genes_per_cell=resolved_max_genes,
        max_pct_mitochondrial=resolved_max_mt,
        min_cells_per_gene=spec.min_cells_per_gene,
        filter_confounding_genes=spec.filter_confounding_genes,
    )


def extract_qc_metrics_dataframe(adata: ad.AnnData) -> pl.DataFrame:
    """Extract cell-level QC covariates from adata.obs into a typed Polars DataFrame."""
    cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    for c in ["pct_counts_ribo", "pct_counts_hb", "patient", "sample", "cell_type_author"]:
        if c in adata.obs.columns:
            cols.append(c)

    data = {c: adata.obs[c].to_numpy() for c in cols}
    data["cell_id"] = adata.obs_names.to_numpy().astype(str)
    return pl.DataFrame(data)


def apply_quality_control(
    adata: ad.AnnData,
    qc_spec: QualityControlSpec | None = None,
    use_adaptive_qc: bool = False,
    n_mads: float = 3.0,
    min_cells_per_gene: int = 3,
) -> Result[ad.AnnData, str]:
    """Filter low-quality droplets and non-informative genes per Luecken & Theis (2019).

    Args:
        adata: Input AnnData object with raw counts in X.
        qc_spec: Quality control threshold specification.
        use_adaptive_qc: If True, derives adaptive MAD thresholds from data distributions.
        n_mads: Number of MADs for adaptive thresholding.
        min_cells_per_gene: Minimum number of cells a gene must be expressed in.

    Returns:
        Success(filtered_adata) with QC metadata in .obs and .uns, or Failure(err).
    """
    try:
        if "total_counts" in adata.obs.columns and "n_genes_by_counts" in adata.obs.columns and "pct_counts_mt" in adata.obs.columns:
            new_adata = adata
        else:
            new_adata = compute_qc_covariates(adata)

        if use_adaptive_qc:
            resolved_spec = calculate_adaptive_thresholds(new_adata, n_mads=n_mads, base_spec=qc_spec)
            logger.info(
                "Derived adaptive QC thresholds: min_counts=%d, min_genes=%d, max_genes=%d, max_mt=%.1f%%",
                resolved_spec.min_counts_per_cell,
                resolved_spec.min_genes_per_cell,
                resolved_spec.max_genes_per_cell,
                resolved_spec.max_pct_mitochondrial,
            )
        else:
            resolved_spec = qc_spec or QualityControlSpec(min_cells_per_gene=min_cells_per_gene)

        total_counts = new_adata.obs["total_counts"].to_numpy()
        genes_detected = new_adata.obs["n_genes_by_counts"].to_numpy()
        pct_mt = new_adata.obs["pct_counts_mt"].to_numpy()

        # Cell-level filtering mask
        cell_mask = (
            (genes_detected >= resolved_spec.min_genes_per_cell)
            & (genes_detected <= resolved_spec.max_genes_per_cell)
            & (total_counts >= resolved_spec.min_counts_per_cell)
            & (pct_mt <= resolved_spec.max_pct_mitochondrial)
        )
        if isinstance(resolved_spec.max_counts_per_cell, Some):
            cell_mask &= (total_counts <= resolved_spec.max_counts_per_cell.unwrap())

        new_adata.obs["is_retained_qc"] = cell_mask

        # Apply cell filtering
        if np.all(cell_mask):
            filtered = new_adata
        else:
            filtered = new_adata[cell_mask].copy()

        # Gene-level filtering: keep genes expressed in at least min_cells_per_gene
        min_cells = resolved_spec.min_cells_per_gene
        if min_cells > 0:
            X_filt = filtered.X
            if sp.isspmatrix_csr(X_filt):
                cells_per_gene = np.bincount(X_filt.indices, minlength=X_filt.shape[1])
            elif sp.issparse(X_filt):
                cells_per_gene = np.asarray((X_filt > 0).sum(axis=0)).ravel()
            else:
                cells_per_gene = np.asarray((np.asarray(X_filt) > 0).sum(axis=0)).ravel()

            gene_mask = cells_per_gene >= min_cells
            filtered.var["cells_per_gene"] = cells_per_gene
            if not np.all(gene_mask):
                if sp.isspmatrix_csr(X_filt):
                    col_map = np.full(X_filt.shape[1], -1, dtype=np.int32)
                    col_map[gene_mask] = np.arange(gene_mask.sum(), dtype=np.int32)
                    kept = col_map[X_filt.indices] >= 0
                    new_data = X_filt.data[kept]
                    new_indices = col_map[X_filt.indices[kept]]
                    row_idx = np.repeat(np.arange(X_filt.shape[0]), np.diff(X_filt.indptr))
                    new_row_counts = np.bincount(row_idx[kept], minlength=X_filt.shape[0])
                    new_indptr = np.zeros(X_filt.shape[0] + 1, dtype=X_filt.indptr.dtype)
                    new_indptr[1:] = np.cumsum(new_row_counts)
                    new_X = sp.csr_matrix((new_data, new_indices, new_indptr), shape=(X_filt.shape[0], int(gene_mask.sum())))
                    filtered = ad.AnnData(X=new_X, obs=filtered.obs, var=filtered.var.iloc[gene_mask].copy())
                else:
                    filtered = filtered[:, gene_mask].copy()
            logger.info(
                "Filtered genes expressed in <%d cells: %d -> %d genes retained",
                min_cells,
                len(gene_mask),
                filtered.n_vars,
            )

        n_pre = new_adata.n_obs
        n_post = filtered.n_obs
        retention_pct = (n_post / max(n_pre, 1)) * 100.0

        filtered.uns["qc_spec"] = {
            "min_genes_per_cell": resolved_spec.min_genes_per_cell,
            "max_genes_per_cell": resolved_spec.max_genes_per_cell,
            "min_counts_per_cell": resolved_spec.min_counts_per_cell,
            "max_counts_per_cell": resolved_spec.max_counts_per_cell.value_or(None),
            "max_pct_mitochondrial": resolved_spec.max_pct_mitochondrial,
            "min_cells_per_gene": resolved_spec.min_cells_per_gene,
            "n_cells_pre_qc": n_pre,
            "n_cells_post_qc": n_post,
            "pct_cells_retained": retention_pct,
            "n_genes_pre_qc": new_adata.n_vars,
            "n_genes_post_qc": filtered.n_vars,
        }

        logger.info(
            "Quality control filtering complete: %d -> %d cells (%.1f%% retained), %d -> %d genes",
            n_pre,
            n_post,
            retention_pct,
            new_adata.n_vars,
            filtered.n_vars,
        )
        return Success(filtered)
    except Exception as exc:
        msg = f"Failed to apply quality control: {exc}"
        logger.error(msg)
        return Failure(msg)


__all__ = [
    "apply_quality_control",
    "calculate_adaptive_thresholds",
    "compute_qc_covariates",
    "extract_qc_metrics_dataframe",
]
