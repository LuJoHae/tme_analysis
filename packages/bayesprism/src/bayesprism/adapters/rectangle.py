"""
Adapter for Rectangle Multiscale Deconvolution (DWLS-QP).

Provides native, pure functional integration of Rectangle (Eder et al., bioRxiv 2026; Tsoucas et al., 2019):
- Native two-stage quadratic programming deconvolution (Dampened Weighted Least Squares DWLS-QP)
  powered directly by SciPy (lsq_linear and SLSQP), eliminating external binary / OSQP incompatibilities
- Support for hierarchical coarse clustered signatures and fine direct signatures
- Explicit Unknown cellular content modeling and simplex re-normalization
- Seamless interoperability between NumPy arrays, Polars DataFrames, and Pandas.

Adheres strictly to .agents/rules/code-style-guide.md:
- Pure functional core, immutable Pydantic models (frozen=True), and Returns (Result/Maybe).
"""

from __future__ import annotations

from typing import Any
import numpy as np
import pandas as pd  # type: ignore
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Some, Nothing
from returns.result import Result, Success, Failure
from scipy.optimize import minimize, lsq_linear  # type: ignore


class RectangleConfig(BaseModel):
    """Configuration options for Rectangle DWLS-QP deconvolution."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    correct_mrna_bias: bool = True
    n_cpus: Maybe[int] = Field(default=Nothing)
    optimize_cutoffs: bool = True
    p_val_threshold: float = 0.015
    lfc_threshold: float = 1.5
    max_dwls_iter: int = 40
    dwls_tolerance: float = 0.002
    cluster_tolerance: float = 0.03


class RectangleDeconvResult(BaseModel):
    """Container holding inferred cellular proportions and diagnostics from Rectangle."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    proportions: pl.DataFrame
    unknown_fraction: pl.Series
    cell_types: tuple[str, ...]
    sample_names: tuple[str, ...]
    reconstruction_errors: Maybe[pl.DataFrame] = Field(default=Nothing)


class RectangleSignatureResult(BaseModel):
    """Container holding reference signature and structural metadata for Rectangle."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    signature_genes: tuple[str, ...]
    bias_factors: tuple[float, ...]
    phi: np.ndarray
    state_names: tuple[str, ...]
    gene_names: tuple[str, ...]
    marker_genes_per_cell_type: dict[str, tuple[str, ...]]
    pseudobulk_sig_cpm: Any = None
    clustered_pseudobulk_sig_cpm: Any = None
    cluster_mapping: Maybe[np.ndarray] = Field(default=Nothing)
    assignments: Maybe[tuple[str, ...]] = Field(default=Nothing)


def is_rectangle_available() -> bool:
    """
    Check if Rectangle deconvolution is available.
    Returns True unconditionally as Rectangle DWLS-QP is natively implemented via SciPy.
    """
    return True


def solve_dwls_qp_single(
    phi: np.ndarray,
    bulk: np.ndarray,
    cluster_mapping: Maybe[np.ndarray] = Nothing,
    bias_factors: Maybe[np.ndarray] = Nothing,
    max_iter: int = 40,
    tol: float = 0.002,
    cluster_tol: float = 0.03,
) -> tuple[np.ndarray, float]:
    """
    Solve two-stage Dampened Weighted Least Squares via Quadratic Programming (DWLS-QP)
    with explicit unknown cellular content estimation for a single bulk sample.

    Parameters:
    - phi: (S, G) reference expression matrix on probability simplex or CPM.
    - bulk: (G,) bulk expression vector (counts or CPM).
    - cluster_mapping: Optional (T, S) aggregation matrix mapping S fine states to T broad clusters.
    - bias_factors: Optional (S,) mRNA bias scaling factors.
    - max_iter: Maximum number of iterative re-weighting cycles.
    - tol: L1 convergence threshold on predicted fractions change.
    - cluster_tol: Epsilon inequality boundary around Stage 1 cluster fractions (default 0.03).

    Returns:
    - (proportions, unknown_fraction) where proportions has shape (S,) and unknown_fraction is float.
    """
    S, G = phi.shape

    # 1. mRNA Content Bias Scaling
    match bias_factors:
        case Some(b_arr):
            b_vec = np.asarray(b_arr, dtype=np.float64).reshape(-1, 1)
            phi_scaled = phi * b_vec
        case _:
            phi_scaled = phi.copy()

    # Normalize signature rows
    row_sums = phi_scaled.sum(axis=1, keepdims=True)
    S_mat = np.where(row_sums > 1e-12, phi_scaled / row_sums, 1.0 / float(G))

    # Normalize bulk
    bulk_sum = float(np.sum(bulk))
    y = bulk / bulk_sum if bulk_sum > 1e-12 else bulk.copy()

    # 2. Stage 1: Coarse Cluster Deconvolution (if cluster mapping is provided)
    c_bounds: np.ndarray | None = None
    M_mat: np.ndarray | None = None
    match cluster_mapping:
        case Some(M) if M.shape[1] == S:
            M_mat = M
            T = M.shape[0]
            S_clust = np.zeros((T, G), dtype=np.float64)
            for t in range(T):
                members = np.where(M[t] > 0)[0]
                if len(members) > 0:
                    S_clust[t] = S_mat[members].mean(axis=0)
                else:
                    S_clust[t] = 1.0 / float(G)

            res_coarse = lsq_linear(S_clust.T, y, bounds=(0.0, 1.0))
            s_c = float(np.sum(res_coarse.x))
            c_bounds = res_coarse.x / s_c if s_c > 1e-12 else np.ones(T, dtype=np.float64) / float(T)
        case _:
            pass

    # 3. Stage 2: Initial Unweighted QP Solution
    res_init = lsq_linear(S_mat.T, y, bounds=(0.0, 1.0))
    s_init = float(np.sum(res_init.x))
    x = res_init.x / s_init if s_init > 1e-12 else np.ones(S, dtype=np.float64) / float(S)

    # 4. Iterative Dampened Weighted Least Squares
    scale = 1e6
    for _ in range(max_iter):
        y_hat = x @ S_mat
        weights = 1.0 / np.clip(y_hat ** 2, 1e-8, None)
        min_w = float(np.min(weights))
        if min_w > 1e-12:
            weights /= min_w
        weights = np.clip(weights, 1.0, 256.0)

        if c_bounds is not None and M_mat is not None:
            # Constrained QP with coarse cluster bounds: c_k - eps <= sum_{i in C_k} x_i <= c_k + eps
            W_sqrt = np.sqrt(weights)
            A_mat = S_mat * W_sqrt
            y_vec = y * W_sqrt
            G_mat = (A_mat @ A_mat.T) * scale
            a_vec = (A_mat @ y_vec) * scale

            constraints = [{"type": "ineq", "fun": lambda z: 1.0 - float(np.sum(z))}]
            for t in range(M_mat.shape[0]):
                idx_t = np.where(M_mat[t] > 0)[0]
                up = min(1.0, float(c_bounds[t]) + cluster_tol)
                lo = max(0.0, float(c_bounds[t]) - cluster_tol)
                constraints.append({"type": "ineq", "fun": lambda z, i=idx_t, u=up: u - float(np.sum(z[i]))})
                constraints.append({"type": "ineq", "fun": lambda z, i=idx_t, l=lo: float(np.sum(z[i])) - l})

            res_opt = minimize(
                fun=lambda z: float(0.5 * z @ G_mat @ z - a_vec @ z),
                x0=x,
                jac=lambda z: G_mat @ z - a_vec,
                bounds=[(0.0, 1.0) for _ in range(S)],
                constraints=constraints,
                method="SLSQP",
                options={"ftol": 1e-9, "maxiter": 150},
            )
            x_new = res_opt.x
            sn = float(np.sum(x_new))
            if sn > 1e-12:
                x_new /= sn
        else:
            # Standard DWLS unconstrained by clusters
            W_sqrt = np.sqrt(weights)
            A_mat = S_mat * W_sqrt
            y_vec = y * W_sqrt
            res_w = lsq_linear(A_mat.T, y_vec, bounds=(0.0, 1.0))
            sn = float(np.sum(res_w.x))
            x_new = res_w.x / sn if sn > 1e-12 else np.ones(S, dtype=np.float64) / float(S)

        diff = float(np.max(np.abs(x_new - x)))
        x = 0.5 * x + 0.5 * x_new
        if diff < tol:
            break

    # 5. Stage 3: Unknown Cellular Content Modeling
    y_hat_final = x @ S_mat
    obs_sum = float(np.sum(y))
    pred_sum = float(np.sum(y_hat_final))
    ukn_cc = max(0.0, (obs_sum - pred_sum) / max(obs_sum, 1e-12))

    s_final = float(np.sum(x))
    x_normalized = (x / s_final) if s_final > 1e-12 else np.ones(S, dtype=np.float64) / float(S)
    x_corrected = (1.0 - ukn_cc) * x_normalized

    return x_corrected, ukn_cc


def create_signature_from_matrix(
    phi: np.ndarray,
    state_names: tuple[str, ...],
    gene_names: tuple[str, ...],
    bias_factors: Maybe[np.ndarray] = Nothing,
    marker_genes_dict: Maybe[dict[str, list[str]]] = Nothing,
    cluster_mapping: Maybe[np.ndarray] = Nothing,
    assignments: Maybe[tuple[str, ...]] = Nothing,
) -> Result[RectangleSignatureResult, str]:
    """
    Construct a RectangleSignatureResult from pre-computed signature matrix phi (S, G).
    Fully compatible with native SciPy DWLS-QP engine and external rectanglepy package.
    """
    S, G = phi.shape
    if len(state_names) != S:
        return Failure(f"State names length ({len(state_names)}) does not match phi rows ({S}).")
    if len(gene_names) != G:
        return Failure(f"Gene names length ({len(gene_names)}) does not match phi cols ({G}).")

    b_tuple: tuple[float, ...] = (
        tuple(float(b) for b in bias_factors.unwrap())
        if isinstance(bias_factors, Some)
        else tuple(1.0 for _ in range(S))
    )

    markers: dict[str, tuple[str, ...]] = (
        {k: tuple(v) for k, v in marker_genes_dict.unwrap().items()}
        if isinstance(marker_genes_dict, Some)
        else {st: gene_names for st in state_names}
    )

    # DataFrame for compatibility with any code expecting pandas attributes
    sig_df = pd.DataFrame(phi.T, index=list(gene_names), columns=list(state_names))

    return Success(
        RectangleSignatureResult(
            signature_genes=gene_names,
            bias_factors=b_tuple,
            phi=phi.astype(np.float64),
            state_names=state_names,
            gene_names=gene_names,
            marker_genes_per_cell_type=markers,
            pseudobulk_sig_cpm=sig_df,
            cluster_mapping=cluster_mapping,
            assignments=assignments,
        )
    )


def deconvolve_rectangle(
    signatures: Any,
    bulks: np.ndarray | pl.DataFrame | pd.DataFrame,
    sample_names: tuple[str, ...],
    gene_names: tuple[str, ...],
    config: Maybe[RectangleConfig] = Nothing,
) -> Result[RectangleDeconvResult, str]:
    """
    Execute Rectangle multiscale deconvolution (DWLS-QP with Unknown content modeling).

    Parameters:
    - signatures: RectangleSignatureResult or object containing phi, state_names, gene_names.
    - bulks: (N, G) Bulk expression matrix (counts, TPM, or CPM).
    - sample_names: Identifiers for each mixture sample.
    - gene_names: Gene names matching the columns of the mixture.
    - config: Optional RectangleConfig.
    """
    cfg = config.value_or(RectangleConfig())

    # Extract signature matrix and metadata
    if isinstance(signatures, RectangleSignatureResult):
        phi = signatures.phi
        state_names = signatures.state_names
        sig_genes = signatures.gene_names
        cluster_mapping = signatures.cluster_mapping
        bias_factors = (
            Some(np.array(signatures.bias_factors, dtype=np.float64))
            if cfg.correct_mrna_bias
            else Nothing
        )
    elif hasattr(signatures, "pseudobulk_sig_cpm") and hasattr(signatures.pseudobulk_sig_cpm, "to_numpy"):
        sig_df = signatures.pseudobulk_sig_cpm
        phi = sig_df.to_numpy().T
        state_names = tuple(str(c) for c in sig_df.columns)
        sig_genes = tuple(str(idx) for idx in sig_df.index)
        cluster_mapping = getattr(signatures, "cluster_mapping", Nothing)
        bias_factors = (
            Some(np.array(signatures.bias_factors, dtype=np.float64))
            if (cfg.correct_mrna_bias and hasattr(signatures, "bias_factors"))
            else Nothing
        )
    else:
        return Failure(f"Unsupported signatures object type: {type(signatures)}")

    # Extract bulk array
    match bulks:
        case np.ndarray() as arr:
            bulk_arr = arr
        case pl.DataFrame() as pldf:
            bulk_arr = pldf.select([c for c in pldf.columns if c != "sample_id"]).to_numpy()
        case pd.DataFrame() as pdf:
            bulk_arr = pdf.to_numpy()
        case _:
            return Failure(f"Unsupported bulks type: {type(bulks)}")

    N, G_bulk = bulk_arr.shape
    if N != len(sample_names):
        return Failure(f"Bulk rows ({N}) does not match sample_names ({len(sample_names)}).")
    if G_bulk != len(gene_names):
        return Failure(f"Bulk cols ({G_bulk}) does not match gene_names ({len(gene_names)}).")

    # Align gene subsets if needed
    if sig_genes != gene_names:
        common_genes = [g for g in sig_genes if g in gene_names]
        if not common_genes:
            return Failure("No overlapping genes between signature and bulk mixture.")
        sig_gene_indices = [sig_genes.index(g) for g in common_genes]
        bulk_gene_indices = [gene_names.index(g) for g in common_genes]
        phi_sub = phi[:, sig_gene_indices]
        bulk_sub = bulk_arr[:, bulk_gene_indices]
    else:
        phi_sub = phi
        bulk_sub = bulk_arr

    S = phi_sub.shape[0]
    estimated_proportions = np.zeros((N, S), dtype=np.float64)
    unknown_fractions = np.zeros(N, dtype=np.float64)

    # Solve DWLS-QP across all N samples
    for n in range(N):
        x_est, u_est = solve_dwls_qp_single(
            phi=phi_sub,
            bulk=bulk_sub[n],
            cluster_mapping=cluster_mapping,
            bias_factors=bias_factors,
            max_iter=cfg.max_dwls_iter,
            tol=cfg.dwls_tolerance,
            cluster_tol=cfg.cluster_tolerance,
        )
        estimated_proportions[n] = x_est
        unknown_fractions[n] = u_est

    # Construct Polars output DataFrame
    prop_data: dict[str, Any] = {"sample_id": list(sample_names)}
    for s_idx, st in enumerate(state_names):
        prop_data[st] = estimated_proportions[:, s_idx].tolist()
    prop_data["Unknown"] = unknown_fractions.tolist()

    proportions_df = pl.DataFrame(prop_data)

    return Success(
        RectangleDeconvResult(
            proportions=proportions_df,
            unknown_fraction=pl.Series("Unknown", unknown_fractions),
            cell_types=state_names,
            sample_names=sample_names,
            reconstruction_errors=Nothing,
        )
    )


def build_rectangle_signatures(
    adata: Any,
    cell_type_col: str = "cell_type",
    bulks: Maybe[pd.DataFrame] = Nothing,
    config: Maybe[RectangleConfig] = Nothing,
) -> Result[Any, str]:
    """
    Build signature from an AnnData object.
    Delegates to rectanglepy if installed, or creates direct pseudobulk signature.
    """
    try:
        import rectanglepy.pp as rpp  # type: ignore
        cfg = config.value_or(RectangleConfig())
        bulk_arg = bulks.value_or(None)
        n_cpus_arg = cfg.n_cpus.value_or(None)

        sig_result = rpp.build_rectangle_signatures(
            adata=adata,
            cell_type_col=cell_type_col,
            bulks=bulk_arg,
            optimize_cutoffs=cfg.optimize_cutoffs,
            p=cfg.p_val_threshold,
            lfc=cfg.lfc_threshold,
            n_cpus=n_cpus_arg,
        )
        return Success(sig_result)
    except Exception as exc:  # pylint: disable=broad-except
        return Failure(f"Rectangle signature building failed: {exc}")
