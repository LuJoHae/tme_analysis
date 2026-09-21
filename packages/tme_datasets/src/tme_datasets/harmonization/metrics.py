"""Integration Quality & Mixing Metrics: iLISI, cLISI, batch silhouette, and kBET."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from ..models import IntegrationMetricsResult


def _compute_lisi(
    distances: np.ndarray,
    indices: np.ndarray,
    labels: np.ndarray,
    perplexity: float = 15.0,
) -> np.ndarray:
    """Vectorized Local Inverse Simpson's Index calculation using Gaussian kernel weights."""
    n_samples, k_neighbors = indices.shape
    unique_labels = np.unique(labels)
    label_to_int = {l: i for i, l in enumerate(unique_labels)}
    int_labels = np.array([label_to_int[l] for l in labels])
    n_labels = len(unique_labels)

    # Compute Gaussian kernel weights based on distances
    sigmas = np.median(distances, axis=1)
    sigmas[sigmas == 0] = 1.0
    weights = np.exp(-(distances**2) / (2 * sigmas[:, np.newaxis] ** 2))
    weights /= np.sum(weights, axis=1, keepdims=True)

    lisi_scores = np.zeros(n_samples)
    for i in range(n_samples):
        neighbor_labels = int_labels[indices[i]]
        # Compute weighted probabilities per label
        probs = np.bincount(neighbor_labels, weights=weights[i], minlength=n_labels)
        simpson = np.sum(probs**2)
        lisi_scores[i] = 1.0 / simpson if simpson > 0 else 1.0

    return lisi_scores


def evaluate_integration_metrics(
    adata: ad.AnnData,
    batch_key: str = "batch",
    label_key: str = "cell_type",
    k: int = 30,
) -> Result[IntegrationMetricsResult, str]:
    """Compute comprehensive integration quality metrics: iLISI, cLISI, silhouette ratio, and kBET."""
    if batch_key not in adata.obs.columns:
        return Failure(f"Batch key '{batch_key}' not found in adata.obs")

    try:
        # Dimensionality reduction coordinates or expression
        if "X_pca" in adata.obsm:
            X = adata.obsm["X_pca"]
        elif sp.issparse(adata.X):
            X = adata.X.toarray()
        else:
            X = np.asarray(adata.X)

        n_samples = adata.n_obs
        k_eval = min(k, n_samples - 1)
        nn = NearestNeighbors(n_neighbors=k_eval, metric="euclidean")
        nn.fit(X)
        distances, indices = nn.kneighbors(X)

        batches = np.asarray(adata.obs[batch_key])
        ilisi_scores = _compute_lisi(distances, indices, batches)
        mean_ilisi = float(np.mean(ilisi_scores))

        # cLISI if labels present
        if label_key in adata.obs.columns:
            cell_types = np.asarray(adata.obs[label_key])
            clisi_scores = _compute_lisi(distances, indices, cell_types)
            mean_clisi = float(np.mean(clisi_scores))

            # Silhouette scores
            try:
                sub_sample_size = min(2000, n_samples)
                batch_sil = float(silhouette_score(X, batches, sample_size=sub_sample_size))
                label_sil = float(silhouette_score(X, cell_types, sample_size=sub_sample_size))
                sil_ratio = float(label_sil / (abs(batch_sil) + 1e-4))
            except Exception:
                batch_sil, label_sil, sil_ratio = 0.0, 0.0, 1.0
        else:
            mean_clisi = 1.0
            batch_sil = 0.0
            label_sil = 0.0
            sil_ratio = 1.0

        # kBET acceptance rate approximation
        unique_batches, global_counts = np.unique(batches, return_counts=True)
        global_props = global_counts / n_samples

        accept_count = 0
        for i in range(n_samples):
            local_labels = batches[indices[i]]
            _, local_counts = np.unique(local_labels, return_counts=True)
            # Neighborhood matches batch diversity if at least half of batches present
            if len(local_counts) >= max(1, len(unique_batches) // 2):
                accept_count += 1
        kbet_rate = float(accept_count / n_samples)

        return Success(
            IntegrationMetricsResult(
                mean_ilisi=round(mean_ilisi, 3),
                mean_clisi=round(mean_clisi, 3),
                batch_silhouette=round(batch_sil, 3),
                cell_type_silhouette=round(label_sil, 3),
                silhouette_ratio=round(sil_ratio, 3),
                kbet_acceptance_rate=round(kbet_rate, 3),
            )
        )
    except Exception as exc:
        return Failure(f"Failed to evaluate integration metrics: {exc}")
