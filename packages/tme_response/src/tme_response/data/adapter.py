"""Data adapter bridging tme_datasets AnnData objects into MultiOmicCohort models."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from tme_datasets.paths import find_repo_root, get_raw_dataset_dir
from tme_datasets.preprocessing.metadata import harmonize_obs_metadata
from tme_datasets.query import load_dataset
from tme_datasets.registry import get_dataset_spec, resolve_dataset_id

from ..schemas import MultiOmicCohort
from .filters import filter_baseline_samples
from .genomics import parse_cna_file, parse_maf_mutations


def load_multi_omic_cohort(
    cohort_id: str,
    repo_root: Path | None = None,
    filter_baseline: bool = True,
) -> Result[MultiOmicCohort, str]:
    """Load a synchronized bulk multi-omic cohort from tme_datasets.

    Ingests:
      1. Bulk RNA expression matrix (TPM)
      2. Harmonized clinical annotations (RECIST, binary response, timepoint)
      3. Genomic somatic mutations / TMB (if available in raw cache)
      4. Copy number alterations / CNA (if available in raw cache)
    """
    root = repo_root or find_repo_root()
    canonical_id = resolve_dataset_id(cohort_id)

    # 1. Load AnnData from tme_datasets (preserve HUGO gene symbols for signatures)
    match load_dataset(cohort_id, base_dir=root, normalize_ensembl=False):
        case Failure(err):
            return Failure(f"Failed to load dataset '{cohort_id}': {err}")
        case Success(raw_adata):
            pass

    # 2. Harmonize clinical observations
    match harmonize_obs_metadata(raw_adata.obs):
        case Failure(err):
            return Failure(f"Failed to harmonize metadata for '{cohort_id}': {err}")
        case Success(harmonized_obs):
            pass

    harmonized_obs = harmonized_obs.with_columns(
        pl.Series("sample_id", [str(x) for x in raw_adata.obs_names])
    )

    # 3. Filter baseline if requested
    if filter_baseline:
        match filter_baseline_samples(harmonized_obs):
            case Failure(err):
                return Failure(f"Baseline filtering failed for '{cohort_id}': {err}")
            case Success(filtered_clinical):
                pass
    else:
        filtered_clinical = harmonized_obs

    valid_samples = set(filtered_clinical["sample_id"].to_list())
    sample_indices = [
        idx for idx, name in enumerate(raw_adata.obs_names)
        if str(name) in valid_samples
    ]

    if not sample_indices:
        return Failure(f"No matching samples after filtering for cohort {cohort_id}")

    # 4. Extract expression matrix (samples x genes)
    sample_names = tuple(str(raw_adata.obs_names[i]) for i in sample_indices)

    # Prefer gene_name or hugo_symbol or original_id over Ensembl IDs so signatures (CYT, GEP, etc.) work seamlessly
    symbol_col = next(
        (c for c in ("gene_name", "hugo_symbol", "original_id") if c in raw_adata.var.columns),
        None,
    )
    if symbol_col:
        gene_names = tuple(str(g) for g in raw_adata.var[symbol_col])
    else:
        gene_names = tuple(str(g) for g in raw_adata.var_names)

    X_sub = raw_adata.X[sample_indices, :]
    X_dense = X_sub.toarray() if hasattr(X_sub, "toarray") else np.asarray(X_sub)

    # Construct expression Polars DataFrame
    expr_dict: dict[str, list[object]] = {"sample_id": list(sample_names)}
    for g_idx, g_name in enumerate(gene_names):
        expr_dict[g_name] = [float(v) for v in X_dense[:, g_idx]]

    expr_df = pl.DataFrame(expr_dict)

    # Cast any categorical columns in filtered_clinical to String to prevent concatenation errors
    filtered_clinical = filtered_clinical.with_columns([
        pl.col(c).cast(pl.String)
        for c, dtype in zip(filtered_clinical.columns, filtered_clinical.dtypes)
        if dtype == pl.Categorical
    ])


    # 5. Extract genomic data from raw download cache if present
    raw_dir = get_raw_dataset_dir(canonical_id, repo_root=root)
    tmb_scores: Maybe[pl.DataFrame] = Nothing
    driver_mutations: Maybe[pl.DataFrame] = Nothing
    cna_scores: Maybe[pl.DataFrame] = Nothing

    if raw_dir.is_dir():
        # Look for mutation files (data_mutations.txt, *.maf)
        mut_files = list(raw_dir.rglob("data_mutations.txt")) + list(raw_dir.rglob("*.maf"))
        if mut_files:
            match parse_maf_mutations(mut_files[0]):
                case Success((tmb_df, drv_df)):
                    tmb_scores = Some(tmb_df.filter(pl.col("sample_id").is_in(list(sample_names))))
                    driver_mutations = Some(drv_df.filter(pl.col("sample_id").is_in(list(sample_names))))
                case Failure(_):
                    pass

        # Look for CNA files
        cna_files = list(raw_dir.rglob("data_cna.txt")) + list(raw_dir.rglob("data_linear_cna.txt"))
        if cna_files:
            match parse_cna_file(cna_files[0]):
                case Success(cna_df):
                    cna_scores = Some(cna_df.filter(pl.col("sample_id").is_in(list(sample_names))))
                case Failure(_):
                    pass

    # Resolve cancer type from spec
    spec_res = get_dataset_spec(cohort_id)
    cancer_type = spec_res.map(lambda s: s.cancer_type).value_or("Pan-Cancer")

    return Success(
        MultiOmicCohort(
            cohort_id=cohort_id,
            cancer_type=cancer_type,
            sample_ids=sample_names,
            expression_tpm=expr_df,
            clinical_annotations=filtered_clinical,
            tmb_scores=tmb_scores,
            driver_mutations=driver_mutations,
            cna_scores=cna_scores,
        )
    )
