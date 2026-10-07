"""Malignant (cancer) cell detection and annotation heuristics for deconvolution references."""

from __future__ import annotations

import re
from typing import Sequence
import anndata as ad
import numpy as np
import scipy.sparse as sp

from .models import DeconvolutionReferenceConfig

# Known tumor and malignant marker gene signatures
TUMOR_MARKER_SIGNATURES: dict[str, tuple[str, ...]] = {
    "melanoma": ("MLANA", "PMEL", "MITF", "TYR", "S100B"),
    "carcinoma": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
    "glioma": ("SOX2", "EGFR", "NES", "OLIG2"),
    "pan_cancer_proliferation": ("MKI67", "TOP2A", "PCNA"),
}

MALIGNANT_PATTERN = re.compile(
    r"(malignant|tumor|tumour|neoplastic|cancer|carcinoma|melanoma_tumor)",
    flags=re.IGNORECASE,
)


def detect_malignant_cells(
    adata: ad.AnnData,
    config: DeconvolutionReferenceConfig,
) -> np.ndarray:
    """Detect malignant cells across an AnnData dataset using metadata or gene expression heuristics.

    Returns:
        Boolean numpy array of shape (n_obs,) where True denotes a malignant cell.
    """
    n_obs = adata.n_obs

    # 1. Check explicit user-specified malignant_key
    match config.malignant_key:
        case config.malignant_key if config.malignant_key.value_or(None) is not None:
            col_name = config.malignant_key.unwrap()
            if col_name in adata.obs.columns:
                col_vals = adata.obs[col_name]
                if col_vals.dtype == bool or col_vals.dtype == "boolean":
                    return np.asarray(col_vals, dtype=bool)
                return np.asarray([
                    bool(MALIGNANT_PATTERN.search(str(v)))
                    or str(v).lower() == config.malignant_label.lower()
                    for v in col_vals
                ], dtype=bool)
        case _:
            pass

    # 2. Check standard metadata columns for explicit malignant annotations
    candidate_cols = [
        "is_malignant",
        "malignant",
        config.cell_state_key,
        config.cell_type_key.value_or("cell_type"),
        "cell_type",
        "cell_subtype",
        "lineage",
    ]

    for col in candidate_cols:
        if col in adata.obs.columns:
            vals = adata.obs[col]
            if col in ("is_malignant", "malignant") and (vals.dtype == bool or vals.dtype == "boolean"):
                return np.asarray(vals, dtype=bool)
            flagged = np.asarray([bool(MALIGNANT_PATTERN.search(str(v))) for v in vals], dtype=bool)
            if np.any(flagged):
                return flagged

    # 3. If auto_detect_malignant is enabled, check tumor lineage marker expression
    if config.auto_detect_malignant:
        symbol_col = config.gene_symbol_key.value_or(None)
        if symbol_col is not None and symbol_col in adata.var.columns:
            var_names = [str(g).upper() for g in adata.var[symbol_col]]
        else:
            var_names = [str(g).upper() for g in adata.var_names]
        gene_to_idx = {g: i for i, g in enumerate(var_names)}

        # Search for matched tumor markers
        detected_tumor_genes = []
        for sig_genes in TUMOR_MARKER_SIGNATURES.values():
            for g in sig_genes:
                if g in gene_to_idx:
                    detected_tumor_genes.append(gene_to_idx[g])

        if len(detected_tumor_genes) >= 2:
            X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)
            sub_expr = X[:, detected_tumor_genes]
            if sp.issparse(sub_expr):
                mean_marker = np.asarray(sub_expr.mean(axis=1)).flatten()
            else:
                mean_marker = np.mean(sub_expr, axis=1)

            # Cells with substantial expression in the top 90th percentile above background
            threshold = float(np.percentile(mean_marker, 90))
            if threshold > 0.5:
                return mean_marker >= threshold

    return np.zeros(n_obs, dtype=bool)


def calculate_cnv_proxy_scores(
    adata: ad.AnnData,
    contig_key: str = "contig",
) -> np.ndarray | None:
    """Calculate chromosomal expression variance across genomic contigs as a CNV proxy score.

    Returns:
        np.ndarray of shape (n_obs,) containing variance scores, or None if contig metadata is unavailable.
    """
    if contig_key not in adata.var.columns:
        return None

    chroms = [str(i) for i in range(1, 23)] + ["X", "chr1", "chr2"]
    chrom_means = []

    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)

    for chrom in chroms:
        mask = np.asarray(adata.var[contig_key].astype(str) == chrom)
        if np.sum(mask) >= 10:
            sub = X[:, mask]
            if sp.issparse(sub):
                m = np.asarray(sub.mean(axis=1)).flatten()
            else:
                m = np.mean(sub, axis=1)
            chrom_means.append(m)

    if len(chrom_means) < 5:
        return None

    # Variance across chromosomes per cell: high variance indicates aneuploid CNV events
    chrom_mat = np.vstack(chrom_means)
    return np.var(chrom_mat, axis=0)
