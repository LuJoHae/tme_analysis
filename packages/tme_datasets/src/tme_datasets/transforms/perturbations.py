"""Stochastic perturbation models: Negative Binomial counts, dropout, and Gaussian jitter."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..models import NegativeBinomialConfig
from .nb_inference import compute_size_factors, infer_dataset_nb_parameters


def randomize_negative_binomial(
    adata: ad.AnnData,
    config: NegativeBinomialConfig,
) -> Result[ad.AnnData, str]:
    """Sample randomized single-cell count data according to a Negative Binomial distribution.

    Two modes supported:
    1. Simple entry-wise mode (default, when config.estimation_method is Nothing):
       Uses each entry X_{c,g} as expected rate mu_{c,g} with fixed dispersion alpha.
       E[Y_{c,g}] = mu_{c,g}, Var(Y_{c,g}) = mu_{c,g} + alpha * mu_{c,g}^2.
    2. Parameter inference mode (when config.estimation_method is Some):
       Infers gene-level (mu_g, alpha_g) parameters across the dataset or per cluster
       using Moments, MLE, or Empirical Bayes, then resamples counts:
       lambda_{i,g} ~ Gamma(shape=1/alpha_g, scale=alpha_g * (s_i * mu_g + baseline))
       counts_{i,g} ~ Poisson(lambda_{i,g}).
    """
    try:
        rng_seed = config.seed.value_or(None) if isinstance(config.seed, Some) else None
        rng = np.random.default_rng(rng_seed)

        new_adata = adata.copy()
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=np.float32)

        # Store unperturbed counts in layers
        new_adata.layers["raw_counts"] = adata.X.copy()

        # Check if inference mode is requested
        if isinstance(config.estimation_method, Some):
            method = config.estimation_method.unwrap()
            params_res = infer_dataset_nb_parameters(
                adata,
                method=method,
                cluster_key=config.cluster_key,
                min_dispersion=config.min_dispersion,
                max_dispersion=config.max_dispersion,
            )
            if not isinstance(params_res, Success):
                return Failure(params_res.failure())

            param_dict = params_res.unwrap()
            N, G = X.shape
            sampled_counts = np.zeros((N, G), dtype=np.float32)

            cluster_col = (
                config.cluster_key.value_or(None)
                if isinstance(config.cluster_key, Some)
                else None
            )

            if cluster_col is not None and cluster_col in adata.obs:
                clusters = adata.obs[cluster_col].astype(str).values
                for c, (mu, alpha) in param_dict.items():
                    idx = np.where(clusters == c)[0]
                    if len(idx) == 0:
                        continue
                    sub_X = X[idx]
                    sf = (
                        compute_size_factors(sub_X).reshape(-1, 1)
                        if config.library_size_scaling
                        else np.ones((len(idx), 1), dtype=np.float32)
                    )
                    expected_rates = sf * mu.reshape(1, -1) + config.baseline_rate
                    # Shape = 1 / alpha, scale = alpha * expected_rate
                    alpha_mat = np.maximum(alpha.reshape(1, -1), config.min_dispersion)
                    shape = 1.0 / alpha_mat
                    scale = alpha_mat * expected_rates
                    lam = rng.gamma(shape=shape, scale=np.maximum(scale, 1e-8))
                    lam = np.where(expected_rates > 0, lam, 0.0)
                    sampled_counts[idx] = rng.poisson(lam).astype(np.float32)
            else:
                mu, alpha = param_dict["global"]
                sf = (
                    compute_size_factors(X).reshape(-1, 1)
                    if config.library_size_scaling
                    else np.ones((N, 1), dtype=np.float32)
                )
                expected_rates = sf * mu.reshape(1, -1) + config.baseline_rate
                alpha_mat = np.maximum(alpha.reshape(1, -1), config.min_dispersion)
                shape = 1.0 / alpha_mat
                scale = alpha_mat * expected_rates
                lam = rng.gamma(shape=shape, scale=np.maximum(scale, 1e-8))
                lam = np.where(expected_rates > 0, lam, 0.0)
                sampled_counts = rng.poisson(lam).astype(np.float32)

            # Store inferred parameters in varm
            if "global" in param_dict:
                new_adata.varm["nb_means"] = param_dict["global"][0].reshape(-1, 1)
                new_adata.varm["nb_dispersions"] = param_dict["global"][1].reshape(-1, 1)
        else:
            # Simple Mode (Entry-wise Gamma-Poisson mixture)
            mu = np.clip(X, 0.0, None) + config.baseline_rate
            alpha = max(config.min_dispersion, config.dispersion)
            shape = 1.0 / alpha

            scale = alpha * mu
            lam = rng.gamma(shape=shape, scale=np.maximum(scale, 1e-8))
            lam = np.where(mu > 0, lam, 0.0)
            sampled_counts = rng.poisson(lam).astype(np.float32)

        if sp.issparse(adata.X):
            new_adata.X = sp.csr_matrix(sampled_counts)
        else:
            new_adata.X = sampled_counts

        new_adata.layers["randomized_nb"] = new_adata.X.copy()
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to randomize Negative Binomial counts: {exc}")


def simulate_dropout(
    adata: ad.AnnData,
    rate: float = 0.05,
    seed: int | None = None,
) -> Result[ad.AnnData, str]:
    """Simulate random technical capture dropouts by setting non-zero values to 0 with probability `rate`."""
    if not 0.0 <= rate <= 1.0:
        return Failure(f"Dropout rate must be in [0, 1], got {rate}")

    try:
        rng = np.random.default_rng(seed)
        new_adata = adata.copy()
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X).copy()

        mask = rng.uniform(0.0, 1.0, size=X.shape) < rate
        X[mask] = 0.0

        new_adata.X = sp.csr_matrix(X) if sp.issparse(adata.X) else X
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to simulate dropout: {exc}")


def add_expression_jitter(
    adata: ad.AnnData,
    sigma: float = 0.1,
    seed: int | None = None,
) -> Result[ad.AnnData, str]:
    """Add multiplicative log-normal expression noise to AnnData: X' = X * exp(N(0, sigma^2))."""
    if sigma < 0:
        return Failure(f"Jitter sigma must be non-negative, got {sigma}")

    try:
        rng = np.random.default_rng(seed)
        new_adata = adata.copy()
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X).copy()

        noise = np.exp(rng.normal(0.0, sigma, size=X.shape))
        X_perturbed = X * noise

        new_adata.X = sp.csr_matrix(X_perturbed) if sp.issparse(adata.X) else X_perturbed
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to add expression jitter: {exc}")
