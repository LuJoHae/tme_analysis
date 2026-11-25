"""Pure functional gene signature scoring algorithms for transcriptomic data."""

from __future__ import annotations

from typing import Any, Literal, Sequence, cast
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import stats
import scipy.sparse as sp

from ..simulation.empirical_pseudobulk import resolve_adata_gene_symbols
from .models import GeneSetCollection

RankMethod = Literal["average", "min", "max", "dense", "ordinal"]


def score_geneset_zscore(
    adata: ad.AnnData,
    collection: GeneSetCollection,
) -> Result[pl.DataFrame, str]:
    """Calculate mean standardized Z-scores for each gene set across cells/samples."""
    if adata.n_obs == 0 or adata.n_vars == 0:
        return Failure("Cannot compute Z-score signatures on an empty AnnData container.")
    var_dict = {str(name): idx for idx, name in enumerate(adata.var_names)}
    X = np.asarray(cast(Any, adata.X).toarray() if sp.issparse(adata.X) else adata.X, dtype=np.float64)

    gene_means = np.mean(X, axis=0)
    gene_stds = np.std(X, axis=0)
    safe_stds = np.where(gene_stds > 0.0, gene_stds, 1.0)
    X_z = (X - gene_means) / safe_stds

    scores_dict: dict[str, list[str] | list[float]] = {
        "sample_id": [str(x) for x in adata.obs_names],
        **{
            set_id: (
                [float(v) for v in np.mean(X_z[:, valid_indices], axis=1)]
                if (valid_indices := [var_dict[g] for g in gs.genes if g in var_dict])
                else [0.0] * adata.n_obs
            )
            for set_id, gs in collection.gene_sets.items()
        },
    }
    return Success(pl.DataFrame(scores_dict))


def score_geneset_auc(
    adata: ad.AnnData,
    collection: GeneSetCollection,
    top_fraction: float = 0.20,
) -> Result[pl.DataFrame, str]:
    """Calculate rank-based AUCell enrichment scores for gene sets in each cell/sample."""
    try:
        var_names = list(adata.var_names)
        var_dict = {name: idx for idx, name in enumerate(var_names)}
        X = cast(Any, adata.X).toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

        n_cells, n_genes = X.shape
        k_threshold = max(1, int(n_genes * top_fraction))

        # Compute ranks per cell (descending: highest expression gets rank 1)
        ranks = np.argsort(np.argsort(-X, axis=1), axis=1) + 1

        scores_dict: dict[str, list[str] | list[float]] = {
            "sample_id": [str(x) for x in adata.obs_names]
        }

        for set_id, gs in collection.gene_sets.items():
            valid_indices = [var_dict[g] for g in gs.genes if g in var_dict]
            if not valid_indices:
                scores_dict[set_id] = [0.0] * n_cells
                continue

            set_ranks = ranks[:, valid_indices]
            # AUC over ranks within top threshold
            recovered_in_top = np.sum(set_ranks <= k_threshold, axis=1) / len(valid_indices)
            scores_dict[set_id] = [float(v) for v in recovered_in_top]

        return Success(pl.DataFrame(scores_dict))
    except Exception as exc:
        return Failure(f"Failed to compute rank AUC signatures: {exc}")


def compute_ssgsea_matrix(
    expr_matrix: npt.NDArray[np.float64],
    sample_ids: Sequence[str],
    gene_symbols: Sequence[str],
    collection: GeneSetCollection,
    alpha: float = 0.25,
    rank_method: RankMethod = "max",
) -> Result[pl.DataFrame, str]:
    """Compute exact Bagaev et al. (2021) / Barbie et al. (2009) ssGSEA enrichment scores purely.

    Formula per sample for gene set S within expressed gene universe G:
        E(S) = sum_{g in S}(r_g^(1 + alpha)) / sum_{g in S}(r_g^alpha) - (|G| - |S| + 1) / 2
    where r_g in {1, ..., |G|} are ascending within-sample expression ranks (`method=rank_method`).
    Duplicate Hugo symbols (e.g. from Ensembl-to-symbol projection) are collapsed by taking the
    per-sample maximum expression prior to ranking, matching BostonGene MFP preprocessing.
    """
    if expr_matrix.ndim != 2:
        return Failure(f"Expected 2D expression matrix, got shape {expr_matrix.shape}")
    n_samples, n_cols = expr_matrix.shape
    if n_samples == 0 or n_cols == 0:
        return Failure(f"Cannot compute ssGSEA on empty matrix of shape {expr_matrix.shape}")
    if len(sample_ids) != n_samples or len(gene_symbols) != n_cols:
        return Failure(
            f"Dimension mismatch: matrix={expr_matrix.shape}, samples={len(sample_ids)}, genes={len(gene_symbols)}"
        )

    # Deduplicate Hugo symbols deterministically (sort unique symbols; take column-wise max on ties)
    clean_symbols = [str(g).strip().upper() for g in gene_symbols]
    sym_to_cols: dict[str, list[int]] = {}
    for col_idx, sym in enumerate(clean_symbols):
        if sym and sym != "NAN":
            sym_to_cols.setdefault(sym, []).append(col_idx)

    if not sym_to_cols:
        return Failure("No valid non-empty gene symbols provided for ssGSEA.")

    unique_syms = sorted(sym_to_cols.keys())
    has_dups = any(len(idxs) > 1 for idxs in sym_to_cols.values())
    if has_dups:
        dedup_cols = [
            expr_matrix[:, idxs[0]]
            if len(idxs) == 1
            else np.nanmax(expr_matrix[:, idxs], axis=1)
            for idxs in (sym_to_cols[s] for s in unique_syms)
        ]
        work_mat = np.column_stack(dedup_cols).astype(np.float64, copy=False)
    elif len(unique_syms) == n_cols and unique_syms == clean_symbols:
        work_mat = np.asarray(expr_matrix, dtype=np.float64)
    else:
        col_order = [sym_to_cols[s][0] for s in unique_syms]
        work_mat = np.asarray(expr_matrix[:, col_order], dtype=np.float64)

    # Replace NaNs with -inf so missing values rank at bottom (matching pandas na_option='bottom')
    safe_mat = np.where(np.isnan(work_mat), -np.inf, work_mat)
    ranks = stats.rankdata(safe_mat, method=rank_method, axis=1).astype(np.float64)
    ranks_num = np.power(ranks, 1.0 + alpha)
    ranks_den = np.power(ranks, alpha)

    sym_lookup = {s: idx for idx, s in enumerate(unique_syms)}
    n_total_genes = len(unique_syms)
    out_cols: dict[str, list[str] | list[float]] = {
        "sample_id": [str(sid) for sid in sample_ids],
    }

    for set_id, gs in collection.gene_sets.items():
        set_upper = sorted({str(g).strip().upper() for g in gs.genes})
        valid_idx = [sym_lookup[g] for g in set_upper if g in sym_lookup]
        if not valid_idx:
            out_cols[set_id] = [0.0 for _ in range(n_samples)]
        else:
            num = np.sum(ranks_num[:, valid_idx], axis=1)
            den = np.sum(ranks_den[:, valid_idx], axis=1)
            baseline = 0.5 * float(n_total_genes - len(valid_idx) + 1)
            scores = np.where(den > 0.0, (num / den) - baseline, 0.0)
            out_cols[set_id] = [float(v) for v in scores]

    return Success(pl.DataFrame(out_cols))


def score_geneset_ssgsea(
    adata: ad.AnnData,
    collection: GeneSetCollection,
    alpha: float = 0.25,
    rank_method: RankMethod = "max",
) -> Result[pl.DataFrame, str]:
    """Compute exact Bagaev et al. (2021) ssGSEA scores on an AnnData container."""
    try:
        symbols = resolve_adata_gene_symbols(adata)
        raw_x = cast(Any, adata.X).toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        mat = np.asarray(raw_x, dtype=np.float64)
        sample_ids = tuple(str(x) for x in adata.obs_names)
        return compute_ssgsea_matrix(
            expr_matrix=mat,
            sample_ids=sample_ids,
            gene_symbols=symbols,
            collection=collection,
            alpha=alpha,
            rank_method=rank_method,
        )
    except Exception as exc:
        return Failure(f"Failed to compute ssGSEA signatures: {exc}")


def _median_mad_expr(col: str, group_col: Maybe[str], clip: Maybe[float]) -> pl.Expr:
    """Build pure Polars expression for BostonGene MFP median/MAD scaling."""
    val = pl.col(col).cast(pl.Float64)
    match group_col:
        case Some(g_col):
            grp_expr = pl.col(g_col).cast(pl.String).fill_null("None")
            med = val.median().over(grp_expr)
            mu = val.mean().over(grp_expr)
            mad = (val - mu).abs().mean().over(grp_expr)
        case _:
            med = val.median()
            mu = val.mean()
            mad = (val - mu).abs().mean()
    scaled = pl.when(mad.is_null() | (mad <= 1e-12)).then(0.0).otherwise((val - med) / mad)
    match clip:
        case Some(bound):
            b = float(bound)
            return scaled.clip(-b, b).alias(col)
        case _:
            return scaled.alias(col)


def median_mad_scale_signatures(
    scores_df: pl.DataFrame,
    signature_cols: Sequence[str],
    group_col: Maybe[str] = Nothing,
    clip: Maybe[float] = Some(4.0),
) -> Result[pl.DataFrame, str]:
    """Apply BostonGene/Bagaev median-centering and MAD-scaling globally or within cohort groups."""
    if scores_df.is_empty():
        return Failure("Cannot scale empty signature DataFrame.")
    missing = [c for c in signature_cols if c not in scores_df.columns]
    if missing:
        return Failure(f"Missing signature columns in scores_df: {missing}")

    match group_col:
        case Some(g_col) if g_col not in scores_df.columns:
            return Failure(f"Grouping column '{g_col}' not found in scores_df.")
        case _:
            pass

    exprs = [_median_mad_expr(col, group_col, clip) for col in signature_cols]
    return Success(scores_df.with_columns(exprs))
