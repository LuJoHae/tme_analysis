from typing import Union
import numpy as np
import scipy.sparse as sp
import torch
import polars as pl
from scipy import stats
from returns.result import Result, Success, Failure
from bayesprism.models import MatrixType


def norm_to_one(ref: torch.Tensor, pseudo_min: float = 1e-8) -> torch.Tensor:
    """
    Normalize reference expression matrix such that each row sums to 1.
    If entries are zero, apply pseudo_min adjustment:
    phi = (ref / row_sums) * (1 - pseudo_min * G) + pseudo_min
    For rows where all genes have non-zero expression, simple normalization is used.
    """
    if not isinstance(ref, torch.Tensor):
        ref_t = torch.as_tensor(ref, dtype=torch.float32)
    else:
        ref_t = ref.to(dtype=torch.float32)

    G = ref_t.shape[1]
    if G == 0:
        return ref_t

    row_sums = ref_t.sum(dim=1, keepdim=True)
    row_sums = torch.where(row_sums == 0, torch.ones_like(row_sums), row_sums)

    norm_matrix = (ref_t / row_sums) * (1.0 - pseudo_min * G) + pseudo_min

    # Check rows with min > 0
    min_vals, _ = ref_t.min(dim=1, keepdim=True)
    nonzero_min_mask = min_vals > 0

    if torch.any(nonzero_min_mask):
        direct_norm = ref_t / row_sums
        norm_matrix = torch.where(nonzero_min_mask, direct_norm, norm_matrix)

    return norm_matrix


def collapse(
    ref: MatrixType,
    labels: tuple[str, ...],
) -> tuple[torch.Tensor, tuple[str, ...]]:
    """
    Sum read counts for each unique label level.
    Returns collapsed count matrix (K x G) and unique label tuple.
    """
    if isinstance(ref, sp.spmatrix):
        ref_dense = torch.from_numpy(ref.toarray()).to(dtype=torch.float32)
    elif isinstance(ref, np.ndarray):
        ref_dense = torch.from_numpy(ref).to(dtype=torch.float32)
    else:
        ref_dense = ref.to(dtype=torch.float32)

    unique_labels = tuple(dict.fromkeys(labels))
    K = len(unique_labels)
    G = ref_dense.shape[1]

    collapsed_mat = torch.zeros((K, G), dtype=torch.float32, device=ref_dense.device)

    for i, label in enumerate(unique_labels):
        mask = torch.tensor([l == label for l in labels], device=ref_dense.device)
        collapsed_mat[i, :] = ref_dense[mask, :].sum(dim=0)

    return collapsed_mat, unique_labels


def filter_bulk_outlier(
    mixture: torch.Tensor,
    outlier_cut: float = 0.01,
    outlier_fraction: float = 0.1,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Filter outlier genes in bulk mixture data whose expression fraction
    exceeds outlier_cut in more than outlier_fraction of bulk samples.
    Returns (filtered_mixture, keep_gene_mask).
    """
    row_sums = mixture.sum(dim=1, keepdim=True)
    row_sums = torch.where(row_sums == 0, torch.ones_like(row_sums), row_sums)
    mixture_norm = mixture / row_sums

    outlier_mask = (mixture_norm > outlier_cut).float().sum(dim=0) / mixture.shape[0] > outlier_fraction
    keep_mask = ~outlier_mask

    if not torch.any(keep_mask):
        # If all genes were flagged as outliers (e.g. small mock dataset), retain all genes
        keep_mask = torch.ones(mixture.shape[1], dtype=torch.bool, device=mixture.device)

    filtered_mixture = mixture[:, keep_mask]
    return filtered_mixture, keep_mask


def compute_specificity(
    input_matrix: torch.Tensor,
    pseudo_min: float = 1e-8,
) -> torch.Tensor:
    """Compute maximum specificity score for each gene across cell states/types."""
    ref_ct = norm_to_one(input_matrix, pseudo_min=pseudo_min)
    exp_spec = ref_ct.t() / ref_ct.sum(dim=0, keepdim=True)
    max_spec, _ = exp_spec.max(dim=1)
    return max_spec
