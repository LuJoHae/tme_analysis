"""Pure functional subsampling and feature randomization routines.

Implements:
- Complementary Pairs subsampling (Shah & Samworth 2013)
- Standard subsampling without replacement (Meinshausen & Bühlmann 2010)
- Randomized Lasso feature weighting (W_j ~ Uniform(alpha, 1))
"""

from typing import Tuple
import numpy as np


def generate_complementary_pairs(
    n_samples: int, B: int, seed: int = 42
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    """Generate B complementary pairs of disjoint subsamples of size floor(n / 2).

    Parameters
    ----------
    n_samples : int
        Total number of observations n.
    B : int
        Number of complementary pairs to generate (total subsamples = 2 * B).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    tuple of (subsample_A, subsample_B) index arrays.
    """
    rng = np.random.default_rng(seed)
    m = n_samples // 2
    pairs = []

    for _ in range(B):
        perm = rng.permutation(n_samples)
        pair_a = np.sort(perm[:m])
        pair_b = np.sort(perm[m : 2 * m])
        pairs.append((pair_a, pair_b))

    return tuple(pairs)


def generate_subsamples(
    n_samples: int, B: int, seed: int = 42
) -> tuple[np.ndarray, ...]:
    """Generate B standard subsamples of size floor(n / 2) without replacement.

    Parameters
    ----------
    n_samples : int
        Total number of observations n.
    B : int
        Number of subsamples to generate.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    tuple of subsample index arrays.
    """
    rng = np.random.default_rng(seed)
    m = n_samples // 2
    subsamples = []

    for _ in range(B):
        sample = np.sort(rng.choice(n_samples, size=m, replace=False))
        subsamples.append(sample)

    return tuple(subsamples)


def apply_randomized_weights(
    X: np.ndarray, weakness: float, rng: np.random.Generator
) -> np.ndarray:
    """Randomly reweight feature columns for Randomized Lasso.

    Scales column j by W_j ~ Uniform(weakness, 1.0).
    Under Lasso, minimizing ||y - X beta||^2 + lambda * sum(|beta_j| / W_j)
    is equivalent to scaling X_j by W_j.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix of shape (n_samples, n_features).
    weakness : float
        Weakness parameter alpha in (0, 1]. When weakness=1.0, no randomization occurs.
    rng : np.random.Generator
        NumPy random number generator.

    Returns
    -------
    np.ndarray
        New scaled feature matrix of shape (n_samples, n_features).
    """
    if weakness >= 1.0:
        return X.copy()

    p = X.shape[1]
    weights = rng.uniform(low=weakness, high=1.0, size=p)
    return X * weights


def generate_stratified_complementary_pairs(
    strata: np.ndarray, B: int, seed: int = 42
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    """Generate B complementary pairs stratified by cohort/group labels.

    Within each unique stratum, observations are permuted and split into
    two disjoint halves of size floor(n_stratum / 2). These are aggregated
    across all strata to form complementary pairs (pair_A, pair_B).

    Guarantees:
    - pair_A and pair_B are strictly disjoint: pair_A intersect pair_B = empty.
    - Each stratum is represented proportionally in every subsample.

    Parameters
    ----------
    strata : np.ndarray of shape (n_samples,)
        1D array of categorical labels (e.g. cohort IDs or classes).
    B : int
        Number of complementary pairs to generate.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    tuple of (subsample_A, subsample_B) index arrays.
    """
    rng = np.random.default_rng(seed)
    unique_strata = np.unique(strata)
    strata_indices = {
        label: np.where(strata == label)[0] for label in unique_strata
    }

    pairs = []
    for _ in range(B):
        pair_a_parts = []
        pair_b_parts = []
        for label, idxs in strata_indices.items():
            n_c = len(idxs)
            m_c = n_c // 2
            if m_c > 0:
                perm = rng.permutation(idxs)
                pair_a_parts.append(perm[:m_c])
                pair_b_parts.append(perm[m_c : 2 * m_c])

        pair_a = np.sort(np.concatenate(pair_a_parts)) if pair_a_parts else np.array([], dtype=np.int64)
        pair_b = np.sort(np.concatenate(pair_b_parts)) if pair_b_parts else np.array([], dtype=np.int64)
        pairs.append((pair_a, pair_b))

    return tuple(pairs)


def generate_stratified_subsamples(
    strata: np.ndarray, B: int, seed: int = 42
) -> tuple[np.ndarray, ...]:
    """Generate B stratified subsamples of size sum_c floor(n_c / 2)."""
    rng = np.random.default_rng(seed)
    unique_strata = np.unique(strata)
    strata_indices = {
        label: np.where(strata == label)[0] for label in unique_strata
    }

    subsamples = []
    for _ in range(B):
        sample_parts = []
        for label, idxs in strata_indices.items():
            n_c = len(idxs)
            m_c = n_c // 2
            if m_c > 0:
                sample_parts.append(rng.choice(idxs, size=m_c, replace=False))
        sample = np.sort(np.concatenate(sample_parts)) if sample_parts else np.array([], dtype=np.int64)
        subsamples.append(sample)

    return tuple(subsamples)

