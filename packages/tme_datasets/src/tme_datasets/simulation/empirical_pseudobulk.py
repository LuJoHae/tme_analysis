"""Pure functional empirical sample/patient-grouped pseudobulk aggregation and Parquet caching."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Final, assert_never, cast

import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
import scipy.sparse as sp
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import EmpiricalPseudobulkConfig, EmpiricalPseudobulkResult
from ..paths import (
    find_cached_pseudobulk,
    get_pseudobulk_cache_path,
    validate_no_store_io,
)
from ..preprocessing.metadata import binarize_response, standardize_timepoint
from ..registry import resolve_dataset_id

logger = get_logger("simulation.empirical_pseudobulk")

SAMPLE_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "sample_id",
    "patient_id",
    "patient",
    "donor_id",
    "donor",
    "sample",
    "characteristics: patinet ID (Pre=baseline; Post= on treatment)",
    "characteristics: patient id",
    "characteristics: patient ID",
    "Sample name",
)

RESPONSE_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "response",
    "clinical_response",
    "treatment_response",
    "RECIST",
    "characteristics: response",
    "characteristics: therapeutic response",
)

TREATMENT_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "treatment_status",
    "biopsy_timepoint",
    "timepoint",
    "time_point",
    "biopsy",
)


# ------------------------------------------------------------------------------
# Pure Functional Core
# ------------------------------------------------------------------------------


def resolve_adata_gene_symbols(adata: ad.AnnData) -> tuple[str, ...]:
    """Extract upper-case gene symbols or identifiers from AnnData .var metadata purely."""
    var_cols = tuple(str(c) for c in adata.var.columns)
    for col in ("gene_name", "gene_symbol", "gene_symbols", "GeneSymbol", "feature_name", "symbol"):
        if col in var_cols:
            return tuple(str(g).strip().upper() for g in adata.var[col])

    var_names_list = [str(g).strip().upper() for g in adata.var_names]
    if (
        "original_id" in var_cols
        and len(var_names_list) > 0
        and not var_names_list[0].startswith("ENSG")
    ):
        return tuple(var_names_list)
    if "original_id" in var_cols:
        return tuple(str(g).strip().upper() for g in adata.var["original_id"])
    return tuple(var_names_list)


def compute_empirical_cell_fractions(df_cells: pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Compute exact ground-truth cell-state proportions per sample from single-cell annotations."""
    if df_cells.is_empty():
        return Failure("Cannot compute empirical cell fractions on an empty DataFrame.")

    col_mapping: dict[str, str] = {}
    if "patient_id" in df_cells.columns and "sample_id" not in df_cells.columns:
        col_mapping["patient_id"] = "sample_id"
    if "clinical_response" in df_cells.columns and "response" not in df_cells.columns:
        col_mapping["clinical_response"] = "response"
    if "cell_type" in df_cells.columns and "cell_state" not in df_cells.columns:
        col_mapping["cell_type"] = "cell_state"

    std_df = df_cells.rename(col_mapping) if col_mapping else df_cells
    if "sample_id" not in std_df.columns or "cell_state" not in std_df.columns:
        return Failure(
            f"compute_empirical_cell_fractions requires 'sample_id' and 'cell_state' columns; "
            f"got {list(std_df.columns)}"
        )

    with_resp = (
        std_df
        if "response" in std_df.columns
        else std_df.with_columns(pl.lit("Unknown").alias("response"))
    )
    with_meta = (
        with_resp
        if "treatment_status" in with_resp.columns
        else with_resp.with_columns(pl.lit("Combined").alias("treatment_status"))
    ).with_columns(
        [
            pl.col("sample_id").cast(pl.String),
            pl.col("cell_state").cast(pl.String),
            pl.col("response").cast(pl.String),
            pl.col("treatment_status").cast(pl.String),
        ]
    )

    sample_meta = with_meta.select(["sample_id", "response", "treatment_status"]).unique(
        subset=["sample_id"], keep="first", maintain_order=True
    )
    sample_totals = with_meta.group_by("sample_id").agg(pl.len().alias("n_total"))
    state_counts = with_meta.group_by(["sample_id", "cell_state"]).agg(pl.len().alias("n_cells"))

    unique_samples = sample_meta["sample_id"].to_list()
    unique_states = with_meta["cell_state"].unique().sort().to_list()

    grid = pl.DataFrame(
        {
            "sample_id": [s for s in unique_samples for _ in unique_states],
            "cell_state": [st for _ in unique_samples for st in unique_states],
        }
    )

    joined = (
        grid.join(state_counts, on=["sample_id", "cell_state"], how="left")
        .with_columns(pl.col("n_cells").fill_null(0))
        .join(sample_totals, on="sample_id", how="left")
        .join(sample_meta, on="sample_id", how="left")
        .with_columns((pl.col("n_cells") / pl.col("n_total")).alias("true_fraction"))
    )
    return Success(joined)


def _find_obs_column(
    obs_cols: tuple[str, ...],
    explicit_col: Maybe[str],
    candidates: tuple[str, ...],
    keyword_substrings: tuple[str, ...],
) -> Maybe[str]:
    """Purely resolve an AnnData .obs column name from explicit override or fallback candidates."""
    match explicit_col:
        case Some(col_name):
            return Some(col_name) if col_name in obs_cols else Nothing
        case _:
            found_exact = next((c for c in candidates if c in obs_cols), "")
            if len(found_exact) > 0:
                return Some(found_exact)
            found_kw = next(
                (c for c in obs_cols if any(k in c.lower() for k in keyword_substrings)),
                "",
            )
            return Some(found_kw) if len(found_kw) > 0 else Nothing


def _standardize_response_label(raw_val: object) -> tuple[str, int]:
    """Map heterogeneous clinical response string/number to ('responder'|'non-responder'|'unknown', 1|0|-1)."""
    bin_val = binarize_response(raw_val)
    if not np.isnan(bin_val):
        return ("responder", 1) if bin_val >= 0.5 else ("non-responder", 0)

    text = str(raw_val).strip().lower().replace("_", " ").replace("-", " ")
    if any(
        tok in text
        for tok in ("unfavourable", "progressive", "residual", "non mpr", "npcr", "nmpr", "non responder")
    ):
        return ("non-responder", 0)
    if any(
        tok in text
        for tok in ("favourable", "complete response", "partial response", "pcr", "mpr", "responder")
    ):
        return ("responder", 1)
    return ("unknown", -1)


def resolve_sample_column_name(
    obs_df: pd.DataFrame,
    sample_col: str = "sample_id",
) -> Result[str, str]:
    """Resolve the biological sample/patient grouping column in AnnData .obs."""
    obs_cols = tuple(str(c) for c in obs_df.columns)
    if sample_col in obs_cols:
        return Success(sample_col)
    if sample_col == "sample_id":
        match _find_obs_column(
            obs_cols=obs_cols,
            explicit_col=Nothing,
            candidates=SAMPLE_COL_CANDIDATES,
            keyword_substrings=("sample", "patient", "donor", "patinet"),
        ):
            case Some(resolved):
                return Success(resolved)
            case _:
                pass
    return Failure(
        f"Sample grouping column '{sample_col}' not found in AnnData.obs. "
        f"Available columns: {list(obs_cols)}"
    )


def _build_sample_metadata_df(
    obs_df: pd.DataFrame,
    sample_labels: npt.NDArray[np.str_],
    unique_samples: tuple[str, ...],
    sample_counts: npt.NDArray[np.int64],
    config: EmpiricalPseudobulkConfig,
    sample_metadata_override: Maybe[pl.DataFrame],
) -> pl.DataFrame:
    """Construct harmonized per-sample clinical metadata Polars DataFrame purely."""
    obs_cols = tuple(str(c) for c in obs_df.columns)
    resp_col_maybe = _find_obs_column(
        obs_cols=obs_cols,
        explicit_col=config.response_col,
        candidates=RESPONSE_COL_CANDIDATES,
        keyword_substrings=("response", "recist", "clinical_benefit"),
    )
    tp_col_maybe = _find_obs_column(
        obs_cols=obs_cols,
        explicit_col=config.treatment_col,
        candidates=TREATMENT_COL_CANDIDATES,
        keyword_substrings=("timepoint", "time_point", "biopsy", "treatment_status"),
    )

    first_idx_by_sample: dict[str, int] = {}
    for row_idx, s_val in enumerate(sample_labels.tolist()):
        if s_val not in first_idx_by_sample:
            first_idx_by_sample[s_val] = row_idx

    resp_labels: list[str] = []
    resp_binaries: list[int] = []
    treatment_labels: list[str] = []

    for s_id in unique_samples:
        r_idx = first_idx_by_sample[s_id]
        match resp_col_maybe:
            case Some(r_col):
                r_lbl, r_bin = _standardize_response_label(obs_df[r_col].iloc[r_idx])
            case _:
                r_lbl, r_bin = ("unknown", -1)
        resp_labels.append(r_lbl)
        resp_binaries.append(r_bin)

        match tp_col_maybe:
            case Some(t_col):
                raw_tp = str(obs_df[t_col].iloc[r_idx]).strip()
                std_tp = standardize_timepoint(raw_tp)
                treatment_labels.append(raw_tp if std_tp == "Unknown" and len(raw_tp) > 0 else std_tp)
            case _:
                # Infer Pre/Post from sample_id prefix when formatted like Pre_P1 / Post_P1
                s_lower = s_id.lower()
                inferred_tp = (
                    "Pre"
                    if s_lower.startswith("pre_")
                    else ("Post" if s_lower.startswith("post_") else "Combined")
                )
                treatment_labels.append(inferred_tp)

    base_meta = pl.DataFrame(
        {
            "sample_id": list(unique_samples),
            "n_cells": [int(c) for c in sample_counts.tolist()],
            "response": resp_labels,
            "response_binary": resp_binaries,
            "treatment_status": treatment_labels,
        }
    )

    match sample_metadata_override:
        case Some(override_raw) if not override_raw.is_empty() and "sample_id" in override_raw.columns:
            ov = override_raw.with_columns(pl.col("sample_id").cast(pl.String)).unique(
                subset=["sample_id"], keep="first", maintain_order=True
            )
            ov_map_resp: dict[str, tuple[str, int]] = {}
            if "response" in ov.columns:
                for r in ov.select(["sample_id", "response"]).iter_rows(named=True):
                    raw_r = str(r["response"]).strip()
                    lbl, bval = _standardize_response_label(raw_r)
                    ov_map_resp[str(r["sample_id"])] = (
                        raw_r if len(raw_r) > 0 else lbl,
                        bval,
                    )
            ov_map_tp: dict[str, str] = (
                {
                    str(r["sample_id"]): str(r["treatment_status"])
                    for r in ov.select(["sample_id", "treatment_status"]).iter_rows(named=True)
                }
                if "treatment_status" in ov.columns
                else {}
            )
            merged_rows: list[dict[str, object]] = []
            for row in base_meta.iter_rows(named=True):
                sid = str(row["sample_id"])
                r_str, r_bin = ov_map_resp.get(
                    sid, (str(row["response"]), int(row["response_binary"]))
                )
                t_str = ov_map_tp.get(sid, str(row["treatment_status"]))
                merged_rows.append(
                    {
                        "sample_id": sid,
                        "n_cells": int(row["n_cells"]),
                        "response": r_str,
                        "response_binary": r_bin,
                        "treatment_status": t_str,
                    }
                )
            return pl.DataFrame(merged_rows)
        case _:
            return base_meta


def _select_valid_genes(
    all_genes: tuple[str, ...],
    target_genes: Maybe[frozenset[str]],
) -> tuple[tuple[int, ...], tuple[str, ...]]:
    """Deduplicate valid gene identifiers in .var order with optional target_genes filtering."""
    valid_indices: list[int] = []
    genes: list[str] = []
    seen: set[str] = set()

    has_filter = isinstance(target_genes, Some)
    allowed: frozenset[str] = target_genes.unwrap() if isinstance(target_genes, Some) else frozenset()

    for idx, g in enumerate(all_genes):
        if has_filter and g not in allowed:
            continue
        if len(g) > 0 and g not in {"NONE", "NAN"} and g not in seen:
            valid_indices.append(idx)
            genes.append(g)
            seen.add(g)

    return tuple(valid_indices), tuple(genes)


def _extract_linear_expression(
    adata: ad.AnnData,
    cell_indices: npt.NDArray[np.int32],
    gene_indices: tuple[int, ...],
    normalize_cells_cp10k: bool,
) -> sp.csr_matrix | npt.NDArray[np.float64]:
    """Extract linear-scale single-cell expression slice without mutating input AnnData."""
    base_layer = (
        adata.layers["counts"]
        if normalize_cells_cp10k and "counts" in adata.layers
        else adata.X
    )
    sub_x = base_layer[cell_indices, :][:, list(gene_indices)]

    if sp.issparse(sub_x):
        x_csr: sp.csr_matrix = sp.csr_matrix(cast(Any, sub_x), dtype=np.float64)
        n_rows = int(x_csr.shape[0])
        row_sums_check = (
            np.asarray(x_csr[: min(200, n_rows)].sum(axis=1)).flatten()
            if n_rows > 0
            else np.array([0.0], dtype=np.float64)
        )
        med_row_sum = float(np.median(row_sums_check)) if len(row_sums_check) > 0 else 0.0
        max_val = float(x_csr.data.max()) if n_rows > 0 and len(x_csr.data) > 0 else 0.0

        if max_val < 30.0 and med_row_sum < 500.0:
            x_lin_sp: sp.csr_matrix = x_csr.copy()
            x_lin_sp.data = np.expm1(np.clip(x_lin_sp.data, 0.0, 14.0))
        else:
            x_lin_sp = x_csr

        if not normalize_cells_cp10k:
            return x_lin_sp

        cell_sums = np.asarray(x_lin_sp.sum(axis=1)).flatten()
        safe_sums = np.where(cell_sums > 0.0, cell_sums, 1.0)
        scaled_sp: sp.csr_matrix = sp.csr_matrix(sp.diags(1e4 / safe_sums) @ x_lin_sp, dtype=np.float64)
        return scaled_sp

    x_dense = np.asarray(sub_x, dtype=np.float64)
    n_rows_d = int(x_dense.shape[0])
    row_sums_d = (
        np.asarray(x_dense[: min(200, n_rows_d)].sum(axis=1)).flatten()
        if n_rows_d > 0
        else np.array([0.0], dtype=np.float64)
    )
    med_row_d = float(np.median(row_sums_d)) if len(row_sums_d) > 0 else 0.0
    max_val_d = float(x_dense.max()) if x_dense.size > 0 else 0.0

    x_lin_d = (
        np.asarray(np.expm1(np.clip(x_dense, 0.0, 14.0)), dtype=np.float64)
        if max_val_d < 30.0 and med_row_d < 500.0
        else x_dense
    )
    if not normalize_cells_cp10k:
        return x_lin_d

    dense_sums = x_lin_d.sum(axis=1, keepdims=True)
    safe_dense_sums = np.where(dense_sums > 0.0, dense_sums, 1.0)
    return np.asarray((x_lin_d / safe_dense_sums) * 1e4, dtype=np.float64)


def aggregate_empirical_pseudobulk(
    adata: ad.AnnData,
    config: EmpiricalPseudobulkConfig = EmpiricalPseudobulkConfig(),
    dataset_id: str = "in_memory",
    sample_metadata_override: Maybe[pl.DataFrame] = Nothing,
) -> Result[EmpiricalPseudobulkResult, str]:
    """Aggregate single-cell expression into empirical sample/patient-grouped pseudobulks purely."""
    if adata.n_obs == 0 or adata.n_vars == 0:
        return Failure(
            f"Cannot aggregate empirical pseudobulk on empty AnnData (shape={adata.shape})."
        )

    match resolve_sample_column_name(adata.obs, config.sample_col):
        case Failure(err):
            return Failure(err)
        case Success(resolved_sample_col):
            sample_col_used = resolved_sample_col

    raw_samples = np.asarray(
        [str(s).strip() for s in adata.obs[sample_col_used].to_numpy()],
        dtype=np.str_,
    )
    valid_cell_mask = np.array(
        [len(s) > 0 and s.lower() not in {"nan", "none", "unknown"} for s in raw_samples.tolist()],
        dtype=bool,
    )
    if not np.any(valid_cell_mask):
        return Failure(
            f"No valid non-empty sample labels found in column '{sample_col_used}' for '{dataset_id}'."
        )

    min_cells = max(1, config.min_cells_per_sample)
    uniq_vals, uniq_counts = np.unique(raw_samples[valid_cell_mask], return_counts=True)
    kept_pairs = sorted(
        (str(s), int(c))
        for s, c in zip(uniq_vals.tolist(), uniq_counts.tolist(), strict=True)
        if int(c) >= min_cells
    )
    if not kept_pairs:
        return Failure(
            f"Zero samples met min_cells_per_sample={min_cells} in '{dataset_id}'."
        )

    unique_samples = tuple(p[0] for p in kept_pairs)
    sample_counts = np.array([p[1] for p in kept_pairs], dtype=np.int64)
    sample_to_idx = {s: i for i, s in enumerate(unique_samples)}

    kept_cell_indices = np.array(
        [i for i, s in enumerate(raw_samples.tolist()) if s in sample_to_idx],
        dtype=np.int32,
    )
    kept_cell_samples = raw_samples[kept_cell_indices]

    all_genes = resolve_adata_gene_symbols(adata)
    valid_gene_indices, genes = _select_valid_genes(all_genes, config.target_genes)
    if not genes:
        return Failure(f"Zero valid or matching target genes found for '{dataset_id}'.")

    x_lin = _extract_linear_expression(
        adata=adata,
        cell_indices=kept_cell_indices,
        gene_indices=valid_gene_indices,
        normalize_cells_cp10k=config.normalize_cells_cp10k,
    )

    n_kept = len(kept_cell_indices)
    row_ind = np.array([sample_to_idx[str(s)] for s in kept_cell_samples.tolist()], dtype=np.int32)
    col_ind = np.arange(n_kept, dtype=np.int32)
    indicator_mat = sp.csr_matrix(
        (np.ones(n_kept, dtype=np.float64), (row_ind, col_ind)),
        shape=(len(unique_samples), n_kept),
    )

    raw_bulk = indicator_mat @ x_lin
    bulk_sums: npt.NDArray[np.float64] = (
        np.asarray(raw_bulk.toarray(), dtype=np.float64)
        if sp.issparse(raw_bulk)
        else np.asarray(raw_bulk, dtype=np.float64)
    )

    match config.aggregation:
        case "sum":
            bulk_matrix = bulk_sums
        case "mean":
            bulk_matrix = bulk_sums / np.maximum(sample_counts[:, np.newaxis].astype(np.float64), 1.0)
        case _ as unreachable:
            assert_never(unreachable)

    col_dict: dict[str, object] = {"sample_id": list(unique_samples)}
    for g_idx, g_name in enumerate(genes):
        col_dict[g_name] = bulk_matrix[:, g_idx]
    pseudobulk_df = pl.DataFrame(col_dict)

    kept_obs_df = adata.obs.iloc[kept_cell_indices]
    sample_metadata_df = _build_sample_metadata_df(
        obs_df=kept_obs_df,
        sample_labels=kept_cell_samples,
        unique_samples=unique_samples,
        sample_counts=sample_counts,
        config=config,
        sample_metadata_override=sample_metadata_override,
    )

    true_fractions_maybe: Maybe[pl.DataFrame] = Nothing
    match config.cell_state_col:
        case Some(state_col) if state_col in kept_obs_df.columns:
            meta_lookup = {
                str(r["sample_id"]): (str(r["response"]), str(r["treatment_status"]))
                for r in sample_metadata_df.iter_rows(named=True)
            }
            df_cells = pl.DataFrame(
                {
                    "sample_id": kept_cell_samples.tolist(),
                    "cell_state": [str(v) for v in kept_obs_df[state_col].to_numpy().tolist()],
                    "response": [meta_lookup[str(s)][0] for s in kept_cell_samples.tolist()],
                    "treatment_status": [meta_lookup[str(s)][1] for s in kept_cell_samples.tolist()],
                }
            )
            match compute_empirical_cell_fractions(df_cells):
                case Success(frac_df):
                    true_fractions_maybe = Some(frac_df)
                case Failure(frac_err):
                    return Failure(frac_err)
        case Some(missing_col):
            return Failure(
                f"Configured cell_state_col '{missing_col}' not found in AnnData.obs. "
                f"Available columns: {list(kept_obs_df.columns)}"
            )
        case _:
            pass

    return Success(
        EmpiricalPseudobulkResult(
            dataset_id=dataset_id,
            sample_ids=unique_samples,
            genes=genes,
            bulk_matrix=bulk_matrix,
            pseudobulk_df=pseudobulk_df,
            sample_metadata_df=sample_metadata_df,
            true_fractions_df=true_fractions_maybe,
            cache_path=Nothing,
        )
    )


def empirical_pseudobulk_to_anndata(result: EmpiricalPseudobulkResult) -> ad.AnnData:
    """Convert an EmpiricalPseudobulkResult into a sample-by-gene AnnData object purely."""
    obs_pd = result.sample_metadata_df.to_pandas().set_index("sample_id")
    obs_aligned = obs_pd.reindex(list(result.sample_ids))
    var_pd = pd.DataFrame({"gene_symbol": list(result.genes)}, index=list(result.genes))
    return ad.AnnData(
        X=result.bulk_matrix.copy(),
        obs=obs_aligned,
        var=var_pd,
    )


# ------------------------------------------------------------------------------
# Imperative Shell (Atomic Parquet Caching & Dataset API Loading)
# ------------------------------------------------------------------------------


def _companion_meta_path(parquet_path: Path) -> Path:
    """Return companion sample-metadata Parquet path adjacent to pseudobulk Parquet."""
    return parquet_path.with_name(f"{parquet_path.stem}.sample_meta.parquet")


def _companion_fractions_path(parquet_path: Path) -> Path:
    """Return companion ground-truth cell-fractions Parquet path adjacent to pseudobulk Parquet."""
    return parquet_path.with_name(f"{parquet_path.stem}.true_fractions.parquet")


def _atomic_write_parquet(df: pl.DataFrame, target_path: Path) -> Result[Path, str]:
    """Atomically write a Polars DataFrame to Parquet via temporary sibling + os.replace."""
    try:
        validated = validate_no_store_io(target_path)
        validated.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = validated.with_name(f".{validated.stem}.{os.getpid()}.tmp.parquet")
        df.write_parquet(tmp_path)
        os.replace(tmp_path, validated)
        return Success(validated)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed atomic Parquet write to '{target_path}': {exc}")


def cache_empirical_pseudobulk(
    result: EmpiricalPseudobulkResult,
    cache_dir: Maybe[Path] = Nothing,
    cache_tag: Maybe[str] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[EmpiricalPseudobulkResult, str]:
    """Persist an EmpiricalPseudobulkResult atomically to the configured pseudobulk cache."""
    canonical_id = resolve_dataset_id(result.dataset_id)
    try:
        target_parquet = get_pseudobulk_cache_path(
            dataset_id=canonical_id,
            cache_tag=cache_tag,
            repo_root=repo_root,
            cache_dir=cache_dir,
        )
    except ValueError as store_err:
        return Failure(str(store_err))

    match _atomic_write_parquet(result.pseudobulk_df, target_parquet):
        case Failure(err):
            return Failure(err)
        case Success(written_pb_path):
            pass

    meta_target = _companion_meta_path(written_pb_path)
    match _atomic_write_parquet(result.sample_metadata_df, meta_target):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    match result.true_fractions_df:
        case Some(frac_df):
            frac_target = _companion_fractions_path(written_pb_path)
            match _atomic_write_parquet(frac_df, frac_target):
                case Failure(err):
                    return Failure(err)
                case Success(_):
                    pass
        case _:
            pass

    logger.info(
        "Cached empirical pseudobulk for '%s' (%d samples x %d genes) -> %s",
        canonical_id,
        len(result.sample_ids),
        len(result.genes),
        written_pb_path,
    )
    return Success(
        result.model_copy(
            update={
                "dataset_id": canonical_id,
                "cache_path": Some(written_pb_path),
            }
        )
    )


def _read_cached_pseudobulk_result(
    dataset_id: str,
    parquet_path: Path,
    config: EmpiricalPseudobulkConfig,
    adata: Maybe[ad.AnnData],
    sample_metadata_override: Maybe[pl.DataFrame],
) -> Maybe[Result[EmpiricalPseudobulkResult, str]]:
    """Attempt to satisfy load_empirical_pseudobulk from an existing cached Parquet.

    Returns:
        Some(Success(result)) when cache satisfies all requested samples and genes,
        Some(Failure(err)) if standalone read fails gene intersection,
        Nothing when an in-memory AnnData requires genes/samples absent from cache (triggering recompute).
    """
    try:
        validated = validate_no_store_io(parquet_path)
        pb_df = pl.read_parquet(validated)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Corrupt or unreadable pseudobulk cache '%s' (%s); recomputing.", parquet_path, exc)
        return Nothing

    if "sample_id" not in pb_df.columns:
        return Nothing

    cached_samples = tuple(pb_df["sample_id"].cast(pl.String).to_list())
    cached_genes = tuple(c for c in pb_df.columns if c != "sample_id")
    cached_gene_set = set(cached_genes)

    # When an in-memory AnnData is supplied, verify exact sample set & gene coverage before reusing cache
    match adata:
        case Some(in_mem_adata):
            match resolve_sample_column_name(in_mem_adata.obs, config.sample_col):
                case Failure(_):
                    return Nothing
                case Success(s_col):
                    raw_s = [str(v).strip() for v in in_mem_adata.obs[s_col].to_numpy().tolist()]
                    valid_s = [v for v in raw_s if len(v) > 0 and v.lower() not in {"nan", "none", "unknown"}]
                    u_vals, u_cnts = np.unique(np.asarray(valid_s, dtype=np.str_), return_counts=True)
                    expected_samples = tuple(
                        sorted(
                            str(s)
                            for s, c in zip(u_vals.tolist(), u_cnts.tolist(), strict=True)
                            if int(c) >= max(1, config.min_cells_per_sample)
                        )
                    )
            if expected_samples != cached_samples:
                return Nothing

            _, expected_genes = _select_valid_genes(
                resolve_adata_gene_symbols(in_mem_adata),
                config.target_genes,
            )
            if not set(expected_genes).issubset(cached_gene_set):
                return Nothing
            selected_genes = expected_genes
        case _:
            match config.target_genes:
                case Some(req_genes):
                    selected_genes = tuple(g for g in cached_genes if g in req_genes)
                    if not selected_genes:
                        return Some(
                            Failure(
                                f"Zero requested target_genes found in cached pseudobulk '{validated}'."
                            )
                        )
                case _:
                    selected_genes = cached_genes

    sliced_pb_df = pb_df.select(["sample_id", *selected_genes])
    bulk_mat = sliced_pb_df.select(list(selected_genes)).to_numpy().astype(np.float64)

    meta_file = _companion_meta_path(validated)
    if meta_file.exists() and meta_file.stat().st_size > 0:
        sample_meta_df = pl.read_parquet(meta_file)
    else:
        sample_meta_df = pl.DataFrame(
            {
                "sample_id": list(cached_samples),
                "n_cells": [0 for _ in cached_samples],
                "response": ["unknown" for _ in cached_samples],
                "response_binary": [-1 for _ in cached_samples],
                "treatment_status": ["Combined" for _ in cached_samples],
            }
        )

    match sample_metadata_override:
        case Some(ov_df) if not ov_df.is_empty() and "sample_id" in ov_df.columns:
            # Refresh metadata if caller provided explicit sample metadata override
            ov_clean = ov_df.with_columns(pl.col("sample_id").cast(pl.String)).unique(
                subset=["sample_id"], keep="first", maintain_order=True
            )
            if "response" in ov_clean.columns:
                resp_lookup = {
                    str(r["sample_id"]): _standardize_response_label(r["response"])
                    for r in ov_clean.iter_rows(named=True)
                }
                raw_resp_lookup = {
                    str(r["sample_id"]): str(r["response"])
                    for r in ov_clean.iter_rows(named=True)
                }
                tp_lookup = (
                    {
                        str(r["sample_id"]): str(r["treatment_status"])
                        for r in ov_clean.iter_rows(named=True)
                    }
                    if "treatment_status" in ov_clean.columns
                    else {}
                )
                updated_rows: list[dict[str, object]] = []
                for row in sample_meta_df.iter_rows(named=True):
                    sid = str(row["sample_id"])
                    r_bin = resp_lookup[sid][1] if sid in resp_lookup else int(row.get("response_binary", -1))
                    r_str = raw_resp_lookup.get(sid, str(row.get("response", "unknown")))
                    t_str = tp_lookup.get(sid, str(row.get("treatment_status", "Combined")))
                    updated_rows.append(
                        {
                            "sample_id": sid,
                            "n_cells": int(row.get("n_cells", 0)),
                            "response": r_str,
                            "response_binary": r_bin,
                            "treatment_status": t_str,
                        }
                    )
                sample_meta_df = pl.DataFrame(updated_rows)
        case _:
            pass

    frac_file = _companion_fractions_path(validated)
    frac_maybe: Maybe[pl.DataFrame] = (
        Some(pl.read_parquet(frac_file))
        if frac_file.exists() and frac_file.stat().st_size > 0
        else Nothing
    )

    return Some(
        Success(
            EmpiricalPseudobulkResult(
                dataset_id=dataset_id,
                sample_ids=cached_samples,
                genes=selected_genes,
                bulk_matrix=bulk_mat,
                pseudobulk_df=sliced_pb_df,
                sample_metadata_df=sample_meta_df,
                true_fractions_df=frac_maybe,
                cache_path=Some(validated),
            )
        )
    )


def load_empirical_pseudobulk(
    dataset_id: str,
    config: EmpiricalPseudobulkConfig = EmpiricalPseudobulkConfig(),
    *,
    cache: bool = True,
    force_recompute: bool = False,
    cache_dir: Maybe[Path] = Nothing,
    cache_tag: Maybe[str] = Nothing,
    adata: Maybe[ad.AnnData] = Nothing,
    sample_metadata_override: Maybe[pl.DataFrame] = Nothing,
    base_dir: Maybe[Path] = Nothing,
    auto_download: bool = False,
) -> Result[EmpiricalPseudobulkResult, str]:
    """Load cached empirical sample/patient-grouped pseudobulk or generate and cache from AnnData."""
    canonical_id = resolve_dataset_id(dataset_id)

    try:
        if not force_recompute:
            found_cache = find_cached_pseudobulk(
                dataset_id=canonical_id,
                cache_tag=cache_tag,
                repo_root=base_dir,
                cache_dir=cache_dir,
            )
            match found_cache:
                case Some(cached_parquet):
                    cached_attempt = _read_cached_pseudobulk_result(
                        dataset_id=canonical_id,
                        parquet_path=cached_parquet,
                        config=config,
                        adata=adata,
                        sample_metadata_override=sample_metadata_override,
                    )
                    match cached_attempt:
                        case Some(hit_result):
                            return cast(Result[EmpiricalPseudobulkResult, str], hit_result)
                        case _:
                            pass
                case _:
                    pass
    except ValueError as store_err:
        return Failure(str(store_err))

    # Cache miss or force_recompute: obtain single-cell AnnData
    match adata:
        case Some(provided_adata):
            working_adata = provided_adata
        case _:
            from ..query import load_dataset

            base_dir_arg = base_dir.unwrap() if isinstance(base_dir, Some) else None
            match load_dataset(
                canonical_id,
                base_dir=base_dir_arg,
                auto_download=auto_download,
                normalize_ensembl=False,
            ):
                case Failure(load_err):
                    return Failure(
                        f"Failed to load single-cell dataset '{canonical_id}' for empirical pseudobulk: {load_err}"
                    )
                case Success(loaded_adata):
                    working_adata = loaded_adata

    match aggregate_empirical_pseudobulk(
        adata=working_adata,
        config=config,
        dataset_id=canonical_id,
        sample_metadata_override=sample_metadata_override,
    ):
        case Failure(agg_err):
            return Failure(agg_err)
        case Success(computed):
            if not cache:
                return Success(computed)
            return cache_empirical_pseudobulk(
                result=computed,
                cache_dir=cache_dir,
                cache_tag=cache_tag,
                repo_root=base_dir,
            )
    return Failure("Unreachable")
