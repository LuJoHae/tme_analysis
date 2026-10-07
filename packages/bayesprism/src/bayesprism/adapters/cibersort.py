"""
Adapter for CIBERSORT (linear nu-SVR) Deconvolution.

Provides native functional integration of CIBERSORT (Newman et al., Nature Methods 2015):
- Linear nu-Support Vector Regression (nu-SVR) across candidate margin parameters nu in {0.25, 0.50, 0.75}
- Non-negative projection and proportion simplex normalization
- Mixture reconstruction RMSE and Pearson correlation metrics
- Optional Monte Carlo gene permutation significance testing (empirical p-values)
- Multi-threaded batch execution across bulk samples

Adheres to strict functional programming standards, Pydantic (frozen=True), and Returns (Result/Maybe).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Some, Nothing
from returns.result import Result, Success, Failure
from sklearn.svm import NuSVR  # type: ignore


class CibersortConfig(BaseModel):
    """Configuration options for CIBERSORT deconvolution."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    nu_values: tuple[float, ...] = (0.25, 0.50, 0.75)
    c_param: float = 1.0
    tol: float = 1e-4
    max_iter: int = 2000
    n_perm: int = 0
    n_cpus: Maybe[int] = Field(default=Nothing)


class CibersortDeconvResult(BaseModel):
    """Container holding inferred cell proportions and diagnostics from CIBERSORT."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    proportions: pl.DataFrame
    cell_types: tuple[str, ...]
    sample_names: tuple[str, ...]
    rmse: pl.Series
    correlation: pl.Series
    best_nu: pl.Series
    p_values: Maybe[pl.Series] = Field(default=Nothing)


def _deconvolve_single_sample(
    y_vec: np.ndarray,
    s_std: np.ndarray,
    nu_values: tuple[float, ...],
    c_param: float,
    tol: float,
    max_iter: int,
    n_perm: int,
    seed: int = 42,
) -> tuple[np.ndarray, float, float, float, float | None]:
    """Pure helper function solving nu-SVR for a single mixture sample."""
    K = s_std.shape[1]
    y_std_val = float(np.std(y_vec))
    if y_std_val < 1e-12:
        return np.ones(K) / K, 0.0, 0.0, nu_values[0], None

    y_std = (y_vec - np.mean(y_vec)) / y_std_val

    best_rmse = float("inf")
    best_w = np.ones(K) / K
    best_nu = nu_values[0]
    best_corr = 0.0

    for nu in nu_values:
        svr = NuSVR(nu=nu, C=c_param, kernel="linear", tol=tol, max_iter=max_iter)
        svr.fit(s_std, y_std)
        raw_w = svr.coef_[0]
        pos_w = np.maximum(0.0, raw_w)
        s_sum = float(np.sum(pos_w))
        w = pos_w / s_sum if s_sum > 1e-12 else np.ones(K) / K

        recon = s_std @ w
        rmse = float(np.sqrt(np.mean((recon - y_std) ** 2)))
        recon_std = float(np.std(recon))
        corr = float(np.corrcoef(recon, y_std)[0, 1]) if recon_std > 1e-12 else 0.0

        if rmse < best_rmse:
            best_rmse = rmse
            best_w = w
            best_nu = nu
            best_corr = corr

    # Optional Monte Carlo gene label permutation test
    p_val: float | None = None
    if n_perm > 0:
        rng = np.random.default_rng(seed)
        perm_corrs: list[float] = []
        for _ in range(n_perm):
            y_perm = rng.permutation(y_std)
            svr_p = NuSVR(nu=best_nu, C=c_param, kernel="linear", tol=tol, max_iter=max_iter)
            svr_p.fit(s_std, y_perm)
            w_p = np.maximum(0.0, svr_p.coef_[0])
            s_p = float(np.sum(w_p))
            w_p_norm = w_p / s_p if s_p > 1e-12 else np.ones(K) / K
            r_p = s_std @ w_p_norm
            c_p = float(np.corrcoef(r_p, y_perm)[0, 1]) if np.std(r_p) > 1e-12 else 0.0
            perm_corrs.append(c_p)

        n_better = sum(1 for cp in perm_corrs if cp >= best_corr)
        p_val = float(n_better / n_perm)

    return best_w, best_rmse, best_corr, best_nu, p_val


def deconvolve_cibersort(
    mixture: np.ndarray | pl.DataFrame,
    reference: np.ndarray | pl.DataFrame,
    state_names: Sequence[str],
    sample_names: Sequence[str],
    gene_names: Sequence[str],
    config: Maybe[CibersortConfig] = Nothing,
) -> Result[CibersortDeconvResult, str]:
    """
    Execute CIBERSORT (linear nu-SVR) deconvolution on bulk mixtures.

    Parameters:
    -----------
    mixture: Array or DataFrame of shape (N samples, G genes)
    reference: Array or DataFrame of shape (K states, G genes) or (G genes, K states)
    state_names: Sequence of K cell state/type names
    sample_names: Sequence of N sample names
    gene_names: Sequence of G gene names
    config: Optional CibersortConfig

    Returns:
    --------
    Success(CibersortDeconvResult) or Failure(error_message)
    """
    try:
        # Convert mixture to numpy
        match mixture:
            case pl.DataFrame() as pldf:
                num_cols = [c for c in pldf.columns if c not in ("sample_id", "sample")]
                mix_mat = pldf.select(num_cols).to_numpy().astype(np.float64)
            case np.ndarray() as arr:
                mix_mat = arr.astype(np.float64)
            case _:
                return Failure(f"Unsupported mixture type: {type(mixture)}")

        # Convert reference to numpy
        match reference:
            case pl.DataFrame() as pldf:
                num_cols = [c for c in pldf.columns if c not in ("cluster", "cell_type", "state", "cell_state")]
                ref_mat = pldf.select(num_cols).to_numpy().astype(np.float64)
            case np.ndarray() as arr:
                ref_mat = arr.astype(np.float64)
            case _:
                return Failure(f"Unsupported reference type: {type(reference)}")

        N, G_mix = mix_mat.shape
        K = len(state_names)
        G = len(gene_names)

        if len(sample_names) != N:
            return Failure(f"Sample names count ({len(sample_names)}) != mixture rows ({N}).")
        if G_mix != G:
            return Failure(f"Mixture gene count ({G_mix}) != gene names ({G}).")

        # Determine reference orientation: (K, G) vs (G, K)
        if ref_mat.shape == (K, G):
            S = ref_mat.T  # (G, K)
        elif ref_mat.shape == (G, K):
            S = ref_mat
        else:
            return Failure(
                f"Reference shape {ref_mat.shape} does not match (K={K}, G={G}) or (G={G}, K={K})."
            )

        # Standardize signature matrix S over all elements (as in CIBERSORT.R)
        s_std_val = float(np.std(S))
        if s_std_val < 1e-12:
            return Failure("Reference signature matrix has zero variance.")
        S_std = (S - np.mean(S)) / s_std_val

        cfg = config.value_or(CibersortConfig())
        n_workers = cfg.n_cpus.value_or(1)

        def worker_fn(idx: int) -> tuple[np.ndarray, float, float, float, float | None]:
            return _deconvolve_single_sample(
                y_vec=mix_mat[idx],
                s_std=S_std,
                nu_values=cfg.nu_values,
                c_param=cfg.c_param,
                tol=cfg.tol,
                max_iter=cfg.max_iter,
                n_perm=cfg.n_perm,
                seed=42 + idx,
            )

        if n_workers > 1:
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                results = list(executor.map(worker_fn, range(N)))
        else:
            results = [worker_fn(i) for i in range(N)]

        weights_arr = np.vstack([r[0] for r in results])
        rmse_arr = np.array([r[1] for r in results], dtype=np.float64)
        corr_arr = np.array([r[2] for r in results], dtype=np.float64)
        nu_arr = np.array([r[3] for r in results], dtype=np.float64)

        # Build Polars DataFrame
        data_dict: dict[str, object] = {"sample_id": list(sample_names)}
        for idx, s_name in enumerate(state_names):
            data_dict[s_name] = weights_arr[:, idx]

        proportions_df = pl.DataFrame(data_dict)

        p_vals_series: Maybe[pl.Series] = Nothing
        if cfg.n_perm > 0:
            p_arr = np.array([r[4] if r[4] is not None else 1.0 for r in results], dtype=np.float64)
            p_vals_series = Some(pl.Series("p_value", p_arr))

        return Success(
            CibersortDeconvResult(
                proportions=proportions_df,
                cell_types=tuple(state_names),
                sample_names=tuple(sample_names),
                rmse=pl.Series("rmse", rmse_arr),
                correlation=pl.Series("correlation", corr_arr),
                best_nu=pl.Series("best_nu", nu_arr),
                p_values=p_vals_series,
            )
        )

    except Exception as exc:  # pylint: disable=broad-except
        return Failure(f"CIBERSORT deconvolution encountered an error: {exc}")
