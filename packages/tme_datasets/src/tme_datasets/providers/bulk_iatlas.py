"""Bulk ICI validation cohort providers via cBioPortal / iAtlas."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Failure, Result, Success

from ..preprocessing.metadata import binarize_response, standardize_recist

IATLAS_COHORTS = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
)


def load_iatlas_cohort(cohort_dir: Path) -> Result[ad.AnnData, str]:
    """Load a single cBioPortal / iAtlas bulk cohort directory into an AnnData object."""
    if not cohort_dir.is_dir():
        return Failure(f"Directory does not exist: {cohort_dir}")

    try:
        # 1. Identify mRNA expression file
        expr_files = [
            "data_mrna_seq_expression.txt",
            "data_mrna_seq_tpm.txt",
            "data_mrna_seq_rpkm.txt",
            "data_mrna_seq_fpkm.txt",
        ]
        found_expr = next((cohort_dir / f for f in expr_files if (cohort_dir / f).exists()), None)
        if not found_expr:
            return Failure(f"No mRNA expression file found in {cohort_dir}")

        df_expr = pd.read_csv(found_expr, sep="\t", index_col=0)
        # Drop redundant Hugo symbol or Entrez column if present
        if "Entrez_Gene_Id" in df_expr.columns:
            df_expr = df_expr.drop(columns=["Entrez_Gene_Id"])

        # Samples are columns, genes are index
        sample_ids = list(df_expr.columns)
        genes = list(df_expr.index)

        # 2. Load clinical annotations
        clinical_sample_file = cohort_dir / "data_clinical_sample.txt"
        df_clinical = pd.DataFrame(index=sample_ids)

        if clinical_sample_file.exists():
            # Skip comment lines beginning with '#'
            df_raw_clin = pd.read_csv(clinical_sample_file, sep="\t", comment="#")
            id_col = next((c for c in df_raw_clin.columns if any(k in c.upper() for k in ("SAMPLE_ID", "SAMPLEID"))), df_raw_clin.columns[0])
            df_raw_clin = df_raw_clin.set_index(id_col)
            # Reindex to match expression samples
            common_samples = df_clinical.index.intersection(df_raw_clin.index)
            df_clinical.loc[common_samples, df_raw_clin.columns] = df_raw_clin.loc[common_samples]

        # 3. Standardize response
        resp_col = next((c for c in df_clinical.columns if any(k in c.upper() for k in ("RESPONSE", "RECIST", "CLINICAL_BENEFIT"))), None)
        if resp_col:
            df_clinical["response_binary"] = df_clinical[resp_col].apply(binarize_response)
            df_clinical["response_recist"] = df_clinical[resp_col].apply(standardize_recist)

        X_mat = df_expr.values.T.astype(np.float32)
        adata = ad.AnnData(
            X=X_mat,
            obs=df_clinical,
            var=pd.DataFrame(index=genes),
        )
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load iAtlas cohort from {cohort_dir}: {exc}")
