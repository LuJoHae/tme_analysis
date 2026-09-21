#!/usr/bin/env python3
"""
Step 3: In Silico Pseudo-bulk LOPO Validation.
Generates synthetic pseudo-bulk mixtures from scRNA-seq patient profiles,
deconvolves them against the reference, and benchmarks recovery accuracy of true fractions
and response-associated effect sizes under zero-noise ground-truth conditions.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad  # type: ignore
import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import polars as pl
from scipy import optimize, sparse, stats  # type: ignore


class ValidationConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_path: Path
    patient_col: str
    cluster_col: str
    response_col: str
    reference_path: Path
    hierarchy_path: Path
    sc_effects_path: Path
    min_cells_per_patient: int
    out_dir: Path


def parse_args() -> ValidationConfig:
    parser = argparse.ArgumentParser(
        description="Run in silico pseudo-bulk validation of deconvolution identifiability."
    )
    parser.add_argument("--adata", type=Path, required=True, help="Path to single-cell AnnData")
    parser.add_argument("--patient-col", type=str, default="patient", help="Patient identifier column")
    parser.add_argument("--cluster-col", type=str, default="celltypist_leiden_0.5", help="Cluster column")
    parser.add_argument("--response-col", type=str, default="response", help="Response column")
    parser.add_argument("--reference", type=Path, required=True, help="Path to curated_reference_phi.parquet")
    parser.add_argument("--hierarchy", type=Path, required=True, help="Path to cell_hierarchy.parquet")
    parser.add_argument("--sc-effects", type=Path, required=True, help="Path to sc_patient_response_effects.parquet")
    parser.add_argument("--min-cells", type=int, default=50, help="Minimum cells per patient for pseudo-bulk")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    args = parser.parse_args()
    return ValidationConfig(
        adata_path=args.adata,
        patient_col=args.patient_col,
        cluster_col=args.cluster_col,
        response_col=args.response_col,
        reference_path=args.reference,
        hierarchy_path=args.hierarchy,
        sc_effects_path=args.sc_effects,
        min_cells_per_patient=args.min_cells,
        out_dir=args.out_dir,
    )


def generate_pseudobulk(
    adata: ad.AnnData,
    patient_col: str,
    cluster_col: str,
    response_col: str,
    reference_genes: list[str],
    min_cells: int,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str], list[str]]:
    """Synthesize pseudo-bulk counts and ground-truth proportions per patient."""
    obs_df = pl.from_pandas(adata.obs.reset_index())
    counts = adata.raw.X if adata.raw is not None else adata.X
    all_genes = list(adata.raw.var_names if adata.raw is not None else adata.var_names)

    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    shared_genes = [g for g in reference_genes if g in gene_to_idx]
    sub_indices = [gene_to_idx[g] for g in shared_genes]

    sub_counts = counts[:, sub_indices]
    if sparse.issparse(sub_counts):
        sub_counts = sub_counts.tocsr()

    # Filter patients with >= min_cells
    patient_counts = (
        obs_df.group_by(patient_col)
        .agg([pl.len().alias("n_cells")])
        .filter(pl.col("n_cells") >= min_cells)
    )
    valid_patients = sorted(patient_counts[patient_col].to_list())

    unique_states = sorted(obs_df[cluster_col].unique().to_list())
    n_patients = len(valid_patients)
    n_states = len(unique_states)
    n_genes = len(shared_genes)

    pseudobulk_mat = np.zeros((n_patients, n_genes), dtype=np.float64)
    true_props_mat = np.zeros((n_patients, n_states), dtype=np.float64)
    patient_responses: list[str] = []

    state_to_idx = {s: i for i, s in enumerate(unique_states)}
    obs_patients = obs_df[patient_col].to_numpy()
    obs_clusters = obs_df[cluster_col].to_numpy()

    for p_idx, patient_id in enumerate(valid_patients):
        p_mask = obs_patients == patient_id
        # Sum expression across cells
        pseudobulk_mat[p_idx, :] = np.array(sub_counts[p_mask].sum(axis=0)).flatten()

        # Compute ground truth proportions
        p_clusters = obs_clusters[p_mask]
        total_p_cells = len(p_clusters)
        for s, count in zip(*np.unique(p_clusters, return_counts=True), strict=False):
            if s in state_to_idx:
                true_props_mat[p_idx, state_to_idx[s]] = count / total_p_cells

        # Extract patient response
        resp_vals = obs_df.filter(pl.col(patient_col) == patient_id)[response_col].drop_nulls()
        resp_str = str(resp_vals[0]) if len(resp_vals) > 0 else "Unknown"
        patient_responses.append(resp_str)

    return (
        pseudobulk_mat,
        true_props_mat,
        valid_patients,
        unique_states,
        shared_genes,
        patient_responses,
    )


def fast_nnls_deconvolve(
    pseudobulk_mat: np.ndarray,
    reference_phi: np.ndarray,
) -> np.ndarray:
    """Non-negative least squares deconvolution per pseudo-bulk sample."""
    n_samples = pseudobulk_mat.shape[0]
    n_states = reference_phi.shape[0]
    inferred_theta = np.zeros((n_samples, n_states), dtype=np.float64)

    # reference_phi: states x genes (transpose to genes x states for Ax = b)
    A = reference_phi.T

    for i in range(n_samples):
        b = pseudobulk_mat[i, :]
        sol, _ = optimize.nnls(A, b)
        total = sol.sum()
        inferred_theta[i, :] = sol / total if total > 0 else 1.0 / n_states

    return inferred_theta


def evaluate_insilico_performance(
    true_props: np.ndarray,
    inferred_theta: np.ndarray,
    states: list[str],
    sc_effects_df: pl.DataFrame,
    patient_responses: list[str],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Compute per-state recovery metrics and compare in silico response effects against scRNA-seq."""
    # Binarize response
    y = np.array([
        1 if str(r).lower() in ["responder", "r", "cr", "pr", "yes", "1"] else 0
        for r in patient_responses
    ])

    metrics: list[dict[str, object]] = []

    for idx, state in enumerate(states):
        true_col = true_props[:, idx]
        inf_col = inferred_theta[:, idx]

        # Pearson r and RMSE
        if np.std(true_col) > 1e-6 and np.std(inf_col) > 1e-6:
            r_val, _ = stats.pearsonr(true_col, inf_col)
        else:
            r_val = 0.0

        rmse = float(np.sqrt(np.mean((true_col - inf_col) ** 2)))

        # In silico response association
        inf_std = (inf_col - np.mean(inf_col)) / (np.std(inf_col) + 1e-9)
        pb_r, pb_pval = stats.pointbiserialr(y, inf_std)
        pb_r_clip = np.clip(pb_r, -0.999, 0.999)
        beta_pseudobulk = float(2.0 * pb_r_clip / np.sqrt(1.0 - pb_r_clip**2 + 1e-12))

        metrics.append({
            "cell_state": state,
            "pearson_r_recovery": float(r_val),
            "rmse_recovery": rmse,
            "beta_pseudobulk": beta_pseudobulk,
            "p_val_pseudobulk": float(pb_pval),
            "is_identifiable": bool(r_val >= 0.5),
        })

    metrics_df = pl.DataFrame(metrics)

    # Join with single-cell effects
    comparison_df = metrics_df.join(sc_effects_df, on="cell_state", how="inner")
    return metrics_df, comparison_df


def run_pipeline(config: ValidationConfig) -> Result[None, str]:
    """Execute Step 3 pure in silico validation."""
    if not config.adata_path.exists():
        return Failure(f"AnnData file not found: {config.adata_path}")
    if not config.reference_path.exists():
        return Failure(f"Reference file not found: {config.reference_path}")
    if not config.sc_effects_path.exists():
        return Failure(f"Single-cell effects file not found: {config.sc_effects_path}")

    try:
        adata = ad.read_h5ad(config.adata_path, backed=False)
        phi_df = pl.read_parquet(config.reference_path)
        sc_effects_df = pl.read_parquet(config.sc_effects_path)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to read input files: {exc}")

    ref_genes = [c for c in phi_df.columns if c != "cell_state"]
    ref_states = phi_df["cell_state"].to_list()

    # Generate pseudo-bulk
    (
        pb_mat,
        true_props,
        valid_patients,
        states,
        shared_genes,
        patient_responses,
    ) = generate_pseudobulk(
        adata,
        config.patient_col,
        config.cluster_col,
        config.response_col,
        ref_genes,
        config.min_cells_per_patient,
    )

    if len(valid_patients) < 3:
        return Failure(f"Too few patients with >= {config.min_cells_per_patient} cells: {len(valid_patients)}")

    # Align reference states and genes
    phi_sub = phi_df.filter(pl.col("cell_state").is_in(states)).select(["cell_state"] + shared_genes)
    ref_phi_mat = phi_sub.select(shared_genes).to_numpy().astype(np.float64)

    # Deconvolve
    inferred_theta = fast_nnls_deconvolve(pb_mat, ref_phi_mat)

    # Evaluate
    metrics_df, comparison_df = evaluate_insilico_performance(
        true_props,
        inferred_theta,
        states,
        sc_effects_df,
        patient_responses,
    )

    config.out_dir.mkdir(parents=True, exist_ok=True)
    metrics_df.write_parquet(config.out_dir / "insilico_validation_metrics.parquet")
    comparison_df.write_parquet(config.out_dir / "insilico_recovery_comparison.parquet")

    # Output fractions
    frac_dict: dict[str, list[object]] = {
        "patient": valid_patients,
        "response": patient_responses,
    }
    for s_idx, state_name in enumerate(states):
        frac_dict[f"true_{state_name}"] = true_props[:, s_idx].tolist()
        frac_dict[f"inferred_{state_name}"] = inferred_theta[:, s_idx].tolist()
    pl.DataFrame(frac_dict).write_parquet(config.out_dir / "insilico_pseudobulk_fractions.parquet")

    # Compute overall effect correlation
    if comparison_df.height >= 3:
        sc_b = comparison_df["beta_sc"].to_numpy()
        pb_b = comparison_df["beta_pseudobulk"].to_numpy()
        overall_r, p_val = stats.spearmanr(sc_b, pb_b)
        print(
            f"[INFO] In Silico Pseudo-bulk Validation Complete ({len(valid_patients)} patients). "
            f"Effect size rank correlation: Spearman rho = {overall_r:.3f} (p = {p_val:.3e})."
        )

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] In silico validation results written to {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
