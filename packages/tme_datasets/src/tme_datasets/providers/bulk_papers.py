"""Bulk ICI response dataset providers directly from primary publication deposits."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..preprocessing.metadata import binarize_response, standardize_recist

logger = get_logger("providers.bulk_papers")

PAPER_DATASETS = (
    "Auslander",
    "Chen-CTLA4",
    "Chen-PD1",
    "Freeman",
    "Gide",
    "Hugo",
    "Lauss",
    "Liu",
    "Prat",
    "Ravi",
    "Riaz",
    "Rose",
    "Snyder",
    "VanAllen",
)


def load_paper_h5ad(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Load an authentic paper H5AD dataset from publication archives or dataset_papers/."""
    if not h5ad_path.exists():
        msg = f"H5AD file not found: {h5ad_path}"
        logger.error(msg)
        return Failure(msg)

    try:
        logger.info("Reading paper H5AD file from %s...", h5ad_path.name)
        adata = ad.read_h5ad(h5ad_path)
        # Harmonize response column if present
        resp_cols = [c for c in adata.obs.columns if any(k in c.lower() for k in ("response", "recist", "dcb"))]
        if resp_cols and "response_binary" not in adata.obs.columns:
            resp_col = resp_cols[0]
            logger.debug("Standardizing response column '%s' for %s", resp_col, h5ad_path.name)
            adata.obs["response_binary"] = adata.obs[resp_col].apply(binarize_response)
            adata.obs["response_recist"] = adata.obs[resp_col].apply(standardize_recist)

        logger.info("Successfully parsed %s: %d samples x %d genes", h5ad_path.name, adata.n_obs, adata.n_vars)
        return Success(adata)
    except Exception as exc:
        msg = f"Failed to load paper H5AD from {h5ad_path}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_genentech_egad(align_dir: Path) -> Result[ad.AnnData, str]:
    """Load Genentech IMvigor210 raw RNA-seq alignment read counts (EGAD00001006631)."""
    if not align_dir.exists():
        msg = f"EGAD alignment directory not found: {align_dir}"
        logger.error(msg)
        return Failure(msg)

    try:
        logger.info("Scanning patient subdirectories in %s...", align_dir)
        all_counts = []
        for dirpath in align_dir.iterdir():
            if not dirpath.is_dir():
                continue
            tab_file = dirpath / "ReadsPerGene.out.tab"
            if tab_file.exists():
                patient_id = dirpath.name
                counts = pd.read_csv(tab_file, sep="\t", header=None, skiprows=4, index_col=0)
                # Select stranded_reverse (col 3)
                series = counts[3]
                series.name = patient_id
                all_counts.append(series)

        if not all_counts:
            msg = f"No ReadsPerGene.out.tab files found in {align_dir}"
            logger.error(msg)
            return Failure(msg)

        logger.info("Concatenating %d patient alignment count profiles...", len(all_counts))
        df_counts = pd.concat(all_counts, axis=1).T
        adata = ad.AnnData(
            X=df_counts.values.astype(np.float32),
            obs=pd.DataFrame(index=df_counts.index),
            var=pd.DataFrame(index=df_counts.columns),
        )
        return Success(adata)
    except Exception as exc:
        msg = f"Failed to load Genentech EGAD dataset: {exc}"
        logger.error(msg)
        return Failure(msg)
