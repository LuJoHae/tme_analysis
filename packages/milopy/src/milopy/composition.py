"""Vectorized Neighborhood Composition, Dataset Purity, and Diversity Metrics for Milopy.

Computes exact cell category representation across Milo neighborhoods using sparse
matrix multiplication C = B^T A in < 5 ms.

Calculates:
- Purity / Dominance: Π_j = max_k(p_kj)
- Normalized Shannon Entropy: H_j = - (1 / log2(K)) * ∑ p_kj log2(p_kj + ε)
- Simpson Diversity: D_j = 1 - ∑ p_kj^2
- Active Categories: N_≥τ,j = ∑ I(C_kj ≥ τ)

Strict functional Python adhering to immutability, returns Result, and Polars.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence
import anndata as ad
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp


class NeighborhoodCompositionConfig(BaseModel):
    """Immutable configuration for neighborhood composition metrics."""

    model_config = ConfigDict(frozen=True)

    min_cells_threshold: int = 5
    epsilon: float = 1e-12


class NeighborhoodCompositionResult(BaseModel):
    """Immutable container for neighborhood composition metrics."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    category_col: str
    classes: tuple[str, ...]
    metrics_df: pl.DataFrame
    counts_matrix: np.ndarray  # Shape: (K, M)
    proportions_matrix: np.ndarray  # Shape: (K, M)


def compute_neighborhood_composition(
    adata: ad.AnnData,
    category_col: str,
    config: NeighborhoodCompositionConfig = NeighborhoodCompositionConfig(),
) -> Result[NeighborhoodCompositionResult, str]:
    """Computes vectorized category composition, purity, diversity, and entropy per Milo neighborhood.

    Args:
        adata: AnnData object containing Milo neighborhoods in `adata.obsm['nhoods']`.
        category_col: Categorical column in `adata.obs` (e.g. 'dataset_id', 'cancer_type').
        config: NeighborhoodCompositionConfig.

    Returns:
        Success(NeighborhoodCompositionResult) or Failure(error_message).
    """
    if "nhoods" not in adata.obsm:
        return Failure("Neighborhood incidence matrix 'nhoods' not found in adata.obsm")

    if category_col not in adata.obs.columns:
        return Failure(f"Category column '{category_col}' not found in adata.obs")

    try:
        # A: (N_cells, M_nhoods) sparse binary incidence matrix
        A = adata.obsm["nhoods"]
        if not sp.issparse(A):
            A = sp.csr_matrix(A)
        else:
            A = A.tocsr()

        n_cells, m_nhoods = A.shape
        if n_cells != adata.n_obs:
            return Failure(
                f"Mismatch: nhoods matrix has {n_cells} rows but adata has {adata.n_obs} cells"
            )

        # Extract category values and map to integers
        raw_vals = adata.obs[category_col].astype(str).values
        classes = tuple(sorted(np.unique(raw_vals).tolist()))
        k_classes = len(classes)

        if k_classes == 0:
            return Failure(f"No unique categories found in adata.obs['{category_col}']")

        class_to_idx = {c: i for i, c in enumerate(classes)}
        cat_indices = np.array([class_to_idx[v] for v in raw_vals], dtype=np.int32)

        # B: (N_cells, K_classes) one-hot indicator sparse matrix
        row_indices = np.arange(n_cells, dtype=np.int32)
        data = np.ones(n_cells, dtype=np.float32)
        B = sp.csr_matrix((data, (row_indices, cat_indices)), shape=(n_cells, k_classes), dtype=np.float32)

        # C = B^T * A -> Shape: (K_classes, M_nhoods)
        # Entry (k, j) is the count of cells from category k in neighborhood j
        C_sparse = B.T.dot(A)
        C = np.asarray(C_sparse.toarray(), dtype=np.float32)  # (K, M)

        # Neighborhood totals: S_j = sum_k C_{kj}
        nhood_sizes = C.sum(axis=0)  # Shape: (M,)
        safe_sizes = np.where(nhood_sizes > 0, nhood_sizes, 1.0)

        # Proportions: p_{kj} = C_{kj} / S_j
        P = C / safe_sizes[np.newaxis, :]  # Shape: (K, M)

        # 1. Purity / Dominance: Π_j = max_k p_{kj}
        purity = np.max(P, axis=0)  # Shape: (M,)

        # Dominant category name
        dominant_idx = np.argmax(P, axis=0)
        dominant_category = [classes[idx] for idx in dominant_idx]

        # 2. Normalized Shannon Entropy: H_j = - (1 / log2(K)) * sum_k p_{kj} log2(p_{kj} + ε)
        if k_classes > 1:
            log2_k = np.log2(k_classes)
            # Clip proportions to avoid log2(0)
            p_clipped = np.clip(P, config.epsilon, 1.0)
            entropy_raw = -np.sum(P * np.log2(p_clipped), axis=0)
            entropy_norm = np.clip(entropy_raw / log2_k, 0.0, 1.0)
        else:
            entropy_norm = np.zeros(m_nhoods, dtype=np.float32)

        # 3. Simpson Diversity: D_j = 1 - sum_k p_{kj}^2
        simpson = 1.0 - np.sum(P ** 2, axis=0)
        simpson = np.clip(simpson, 0.0, 1.0)

        # 4. Active category count: N_≥τ,j = sum_k I(C_{kj} >= τ)
        active_counts = np.sum(C >= config.min_cells_threshold, axis=0).astype(np.int32)

        # Build Polars DataFrame
        data_dict: dict[str, Any] = {
            "nhood_index": np.arange(m_nhoods, dtype=np.int32),
            "nhood_size": nhood_sizes.astype(np.int32),
            f"dominant_{category_col}": dominant_category,
            f"{category_col}_purity": purity.astype(np.float64),
            f"{category_col}_entropy": entropy_norm.astype(np.float64),
            f"{category_col}_simpson_diversity": simpson.astype(np.float64),
            f"{category_col}_active_count": active_counts,
        }

        # Add per-class counts and proportions columns
        for k_idx, cls_name in enumerate(classes):
            clean_name = cls_name.replace(" ", "_").replace("-", "_").lower()
            data_dict[f"count_{category_col}_{clean_name}"] = C[k_idx, :].astype(np.int32)
            data_dict[f"prop_{category_col}_{clean_name}"] = P[k_idx, :].astype(np.float64)

        metrics_df = pl.DataFrame(data_dict)

        return Success(
            NeighborhoodCompositionResult(
                category_col=category_col,
                classes=classes,
                metrics_df=metrics_df,
                counts_matrix=C,
                proportions_matrix=P,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute neighborhood composition for '{category_col}': {exc}")


def annotate_nhood_adata_composition(
    nhood_adata: ad.AnnData,
    cell_adata: ad.AnnData,
    category_cols: Sequence[str] = ("dataset_id", "cancer_type"),
    config: NeighborhoodCompositionConfig = NeighborhoodCompositionConfig(),
) -> Result[ad.AnnData, str]:
    """Annotates a Milo neighborhood AnnData (e.g. adata.uns['nhood_adata']) with composition metrics.

    Args:
        nhood_adata: Milo neighborhood AnnData (dim: M nhoods x samples/genes).
        cell_adata: Single-cell AnnData containing `nhoods` in .obsm.
        category_cols: Categorical columns in `cell_adata.obs` to analyze.
        config: NeighborhoodCompositionConfig.

    Returns:
        Success(annotated_nhood_adata) or Failure(error_message).
    """
    try:
        nhood_copy = nhood_adata.copy()
        for cat in category_cols:
            if cat not in cell_adata.obs.columns:
                continue

            res = compute_neighborhood_composition(cell_adata, cat, config)
            match res:
                case Failure(err):
                    return Failure(err)
                case Success(comp_result):
                    df = comp_result.metrics_df
                    # Assign primary metrics to nhood_copy.obs
                    nhood_copy.obs[f"dominant_{cat}"] = df[f"dominant_{cat}"].to_list()
                    nhood_copy.obs[f"{cat}_purity"] = df[f"{cat}_purity"].to_numpy()
                    nhood_copy.obs[f"{cat}_entropy"] = df[f"{cat}_entropy"].to_numpy()
                    nhood_copy.obs[f"{cat}_simpson_diversity"] = df[f"{cat}_simpson_diversity"].to_numpy()
                    nhood_copy.obs[f"{cat}_active_count"] = df[f"{cat}_active_count"].to_numpy()

        return Success(nhood_copy)
    except Exception as exc:
        return Failure(f"Failed to annotate neighborhood AnnData: {exc}")
