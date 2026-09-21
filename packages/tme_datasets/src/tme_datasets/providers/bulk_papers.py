"""Bulk ICI response dataset providers directly from primary publication deposits."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Failure, Result, Success

from ..preprocessing.metadata import binarize_response, standardize_recist

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
        return Failure(f"H5AD file not found: {h5ad_path}")

    try:
        adata = ad.read_h5ad(h5ad_path)
        # Harmonize response column if present
        resp_cols = [c for c in adata.obs.columns if any(k in c.lower() for k in ("response", "recist", "dcb"))]
        if resp_cols and "response_binary" not in adata.obs.columns:
            resp_col = resp_cols[0]
            adata.obs["response_binary"] = adata.obs[resp_col].apply(binarize_response)
            adata.obs["response_recist"] = adata.obs[resp_col].apply(standardize_recist)

        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load paper H5AD from {h5ad_path}: {exc}")


def load_genentech_egad(align_dir: Path) -> Result[ad.AnnData, str]:
    """Load Genentech IMvigor210 raw RNA-seq alignment read counts (EGAD00001006631)."""
    if not align_dir.exists():
        return Failure(f"EGAD alignment directory not found: {align_dir}")

    try:
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
            return Failure(f"No ReadsPerGene.out.tab files found in {align_dir}")

        df_counts = pd.concat(all_counts, axis=1).T
        adata = ad.AnnData(
            X=df_counts.values.astype(np.float32),
            obs=pd.DataFrame(index=df_counts.index),
            var=pd.DataFrame(index=df_counts.columns),
        )
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load Genentech EGAD dataset: {exc}")
