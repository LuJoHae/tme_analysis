"""Variance-stabilizing transformation via regularized Negative Binomial regression and Analytic Pearson Residuals.

References:
    - Hafemeister, C., & Satija, R. (2019).
      Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.
      Genome Biology, 20(1), 296.
    - Lause, J., Berens, P., & Kobak, D. (2021).
      Analytic Pearson residuals for normalization and variable gene selection in single-cell RNA-seq.
      Genome Biology, 22(1), 258.
"""

from __future__ import annotations

import time
import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import SCTransformConfig
from ..types import SCTransformFlavor

logger = get_logger("preprocessing.sctransform")


def _kernel_smooth(x: np.ndarray, y: np.ndarray, bandwidth: float = 0.5) -> np.ndarray:
    """Gaussian kernel regression to smooth parameter y over independent variable x."""
    # Distance matrix
    diff = x[:, None] - x[None, :]
    weights = np.exp(-0.5 * (diff / max(bandwidth, 1e-4)) ** 2)
    # Normalize rows
    row_sums = np.sum(weights, axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-12)
    weights /= row_sums
    return np.dot(weights, y)


def _normalize_analytic_pearson(
    adata: ad.AnnData,
    config: SCTransformConfig,
) -> Result[ad.AnnData, str]:
    """Compute Analytic Pearson Residuals using closed-form offset Negative Binomial (Lause et al., 2021)."""
    try:
        start_time = time.time()
        logger.info(
            "Computing Analytic Pearson Residuals for %d cells x %d genes (theta=%.2f)...",
            adata.n_obs,
            adata.n_vars,
            config.theta,
        )
        new_adata = adata.copy()

        # Preserve unnormalized counts in layers
        new_adata.layers["raw_counts"] = adata.X.copy()

        # Optional HVG selection via Pearson residuals flavor
        if isinstance(config.n_top_genes, Some):
            n_top = min(config.n_top_genes.unwrap(), new_adata.n_vars)
            logger.debug("Calculating highly variable genes (top %d)...", n_top)
            sc.experimental.pp.highly_variable_genes(
                new_adata,
                flavor="pearson_residuals",
                n_top_genes=n_top,
                subset=False,
            )

        # Compute Pearson residuals
        # scanpy.experimental.pp.normalize_pearson_residuals modifies in-place or into layer
        theta = config.theta
        clip = (
            float(np.sqrt(new_adata.n_obs))
            if config.clip_residuals and not isinstance(config.max_residual, Some)
            else (config.max_residual.value_or(None) if isinstance(config.max_residual, Some) else None)
        )

        sc.experimental.pp.normalize_pearson_residuals(
            new_adata,
            theta=theta,
            clip=clip,
            check_values=True,
        )

        # Store residuals in layers
        residuals = new_adata.X.copy()
        new_adata.layers["pearson_residuals"] = residuals

        if not config.use_layer_as_x:
            new_adata.X = new_adata.layers["raw_counts"].copy()

        elapsed = max(0.01, time.time() - start_time)
        logger.info(
            "Successfully completed Analytic Pearson Residuals normalization in %.2fs",
            elapsed,
        )
        return Success(new_adata)
    except Exception as exc:
        msg = f"Analytic Pearson Residuals normalization failed: {exc}"
        logger.error(msg)
        return Failure(msg)


def _normalize_regularized_glm(
    adata: ad.AnnData,
    config: SCTransformConfig,
) -> Result[ad.AnnData, str]:
    """Compute regularized Negative Binomial regression residuals (Hafemeister & Satija, 2019).

    1. Fits log-linear regression of gene counts on cell depth log10(N_c).
    2. Estimates dispersion parameter theta = 1 / alpha.
    3. Smooths parameters across genes via kernel regression over log10(mean expression).
    4. Computes clipped Pearson residuals z_{c,g} = (n_{c,g} - mu_{c,g}) / sigma_{c,g}.
    """
    try:
        start_time = time.time()
        logger.info(
            "Fitting regularized Negative Binomial GLM for %d cells x %d genes (theta=%.2f)...",
            adata.n_obs,
            adata.n_vars,
            config.theta,
        )
        new_adata = adata.copy()
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=np.float64)

        N_cells, G_genes = X.shape
        new_adata.layers["raw_counts"] = adata.X.copy()

        # Cell sequencing depth (covariate)
        depths = np.sum(X, axis=1)
        depths = np.maximum(depths, 1.0)
        log10_depth = np.log10(depths)
        # Center predictor
        x_pred = log10_depth - np.mean(log10_depth)
        X_design = np.column_stack([np.ones(N_cells), x_pred])

        # Gene mean expression
        gene_means = np.mean(X, axis=0)
        log10_means = np.log10(np.maximum(gene_means, 1e-4))

        # 1. Fit per-gene OLS / log-linear regression
        # log(y + 1) ~ beta_0 + beta_1 * x_pred
        logger.debug("Fitting per-gene linear coefficients on cell depths...")
        log_y = np.log(X + 1.0)
        # Solve beta for all genes: (X^T X)^-1 X^T log_y
        beta, _, _, _ = np.linalg.lstsq(X_design, log_y, rcond=None)
        beta_0_raw = beta[0]
        beta_1_raw = beta[1]

        # 2. Kernel smoothing across genes over log10(gene_means)
        logger.debug("Applying Gaussian kernel regression smoothing over log10 gene means...")
        sort_idx = np.argsort(log10_means)
        sorted_log_means = log10_means[sort_idx]

        beta_0_smooth = np.empty_like(beta_0_raw)
        beta_1_smooth = np.empty_like(beta_1_raw)

        beta_0_smooth[sort_idx] = _kernel_smooth(sorted_log_means, beta_0_raw[sort_idx], bandwidth=0.3)
        beta_1_smooth[sort_idx] = _kernel_smooth(sorted_log_means, beta_1_raw[sort_idx], bandwidth=0.3)

        # 3. Predict expected counts: mu_{c,g} = exp(X_design @ beta_smooth)
        beta_smooth_mat = np.vstack([beta_0_smooth, beta_1_smooth])
        log_mu = np.dot(X_design, beta_smooth_mat)
        mu = np.exp(np.clip(log_mu, -20.0, 20.0))

        # 4. Variance with overdispersion theta: Var = mu + mu^2 / theta
        theta = max(1.0, config.theta)
        var_mat = mu + (mu**2) / theta
        sigma = np.sqrt(np.maximum(var_mat, 1e-12))

        # Pearson residuals
        residuals = (X - mu) / sigma

        # Clip residuals
        clip_val = (
            float(np.sqrt(N_cells))
            if config.clip_residuals and not isinstance(config.max_residual, Some)
            else (config.max_residual.value_or(None) if isinstance(config.max_residual, Some) else None)
        )
        if clip_val is not None:
            residuals = np.clip(residuals, -clip_val, clip_val)

        residuals = residuals.astype(np.float32)
        new_adata.layers["pearson_residuals"] = sp.csr_matrix(residuals) if sp.issparse(adata.X) else residuals

        if config.use_layer_as_x:
            new_adata.X = residuals

        # Record gene statistics
        new_adata.var["sct_mean"] = gene_means
        new_adata.var["sct_beta0"] = beta_0_smooth
        new_adata.var["sct_beta1"] = beta_1_smooth

        elapsed = max(0.01, time.time() - start_time)
        logger.info(
            "Successfully completed regularized GLM SCTransform in %.2fs",
            elapsed,
        )
        return Success(new_adata)
    except Exception as exc:
        msg = f"Regularized GLM SCTransform normalization failed: {exc}"
        logger.error(msg)
        return Failure(msg)


def normalize_sctransform(
    adata: ad.AnnData,
    config: SCTransformConfig = SCTransformConfig(),
) -> Result[ad.AnnData, str]:
    """Normalize single-cell count data via variance-stabilizing transformation.

    Supports two flavors:
    1. `SCTransformFlavor.ANALYTIC` (Default):
       Analytic Pearson Residuals (Lause et al., 2021). Fast, closed-form, native Scanpy implementation.
    2. `SCTransformFlavor.REGULARIZED_GLM`:
       Regularized Negative Binomial regression with parameter kernel smoothing (Hafemeister & Satija, 2019).
    """
    logger.info(
        "Starting SCTransform normalization (flavor=%s, %d cells x %d genes)...",
        config.flavor.value,
        adata.n_obs,
        adata.n_vars,
    )
    match config.flavor:
        case SCTransformFlavor.ANALYTIC:
            return _normalize_analytic_pearson(adata, config)
        case SCTransformFlavor.REGULARIZED_GLM:
            return _normalize_regularized_glm(adata, config)
