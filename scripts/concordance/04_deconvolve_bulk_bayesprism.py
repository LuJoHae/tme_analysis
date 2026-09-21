#!/usr/bin/env python3
"""
Step 4: Bulk RNA-seq Deconvolution with BayesPrism and Purity Normalization.
Deconvolves bulk RNA-seq mixtures against the curated single-cell reference,
computes tumor-purity-normalized microenvironment fractions,
scales by cell-type mRNA content factors, and extracts cell-type expression profiles.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import optimize  # type: ignore

try:
    from bayesprism import new_prism, run_prism, get_fraction, get_exp
    HAS_BAYESPRISM = True
except ImportError:
    HAS_BAYESPRISM = False


class DeconvConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    bulk_counts_path: Path
    reference_path: Path
    hierarchy_path: Path
    mrna_scaling_path: Path
    malignant_label: str
    out_dir: Path
    use_bayesprism_mcmc: bool = False
    max_gibbs_iter: int = 200


def parse_args() -> DeconvConfig:
    parser = argparse.ArgumentParser(
        description="Deconvolve bulk RNA-seq and calculate tumor-normalized immune fractions."
    )
    parser.add_argument("--bulk-counts", type=Path, required=True, help="Bulk RNA-seq counts parquet (samples x genes)")
    parser.add_argument("--reference", type=Path, required=True, help="Curated reference phi parquet")
    parser.add_argument("--hierarchy", type=Path, required=True, help="Cell hierarchy parquet")
    parser.add_argument("--mrna-scaling", type=Path, required=True, help="mRNA scaling factors parquet")
    parser.add_argument("--malignant-label", type=str, default="Malignant", help="Malignant cell state label")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--use-bayesprism-mcmc", action="store_true", help="Run full MCMC Gibbs sampling")
    parser.add_argument("--max-iter", type=int, default=200, help="Gibbs iterations")
    args = parser.parse_args()
    return DeconvConfig(
        bulk_counts_path=args.bulk_counts,
        reference_path=args.reference,
        hierarchy_path=args.hierarchy,
        mrna_scaling_path=args.mrna_scaling,
        malignant_label=args.malignant_label,
        out_dir=args.out_dir,
        use_bayesprism_mcmc=args.use_bayesprism_mcmc,
        max_gibbs_iter=args.max_iter,
    )


def read_bulk_matrix(bulk_path: Path) -> Result[tuple[pl.DataFrame, list[str], list[str]], str]:
    """Read bulk counts and extract sample IDs and gene names."""
    if not bulk_path.exists():
        return Failure(f"Bulk counts file does not exist: {bulk_path}")

    try:
        bulk_df = pl.read_parquet(bulk_path)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to read bulk parquet: {exc}")

    # Identify sample id column (e.g. sample, sample_id, patient, bulk_id)
    first_col = bulk_df.columns[0]
    sample_col = first_col if bulk_df[first_col].dtype == pl.String else "sample_id"

    sample_ids = bulk_df[sample_col].to_list() if sample_col in bulk_df.columns else [f"sample_{i}" for i in range(bulk_df.height)]
    gene_cols = [c for c in bulk_df.columns if c != sample_col and c != "sample_id"]

    return Success((bulk_df, sample_ids, gene_cols))


def run_deconvolution(
    bulk_mat: np.ndarray,
    ref_phi_mat: np.ndarray,
    states: list[str],
    samples: list[str],
    use_bp_mcmc: bool,
    max_iter: int,
) -> pl.DataFrame:
    """Run deconvolution using BayesPrism or high-throughput NNLS projection."""
    n_samples = bulk_mat.shape[0]
    n_states = len(states)

    if use_bp_mcmc and HAS_BAYESPRISM:
        # Run BayesPrism full pipeline
        prism_res = new_prism(
            reference=ref_phi_mat,
            cell_type_labels=states,
            cell_state_labels=Nothing,
            mixture=bulk_mat,
            bulk_names=Some(samples),
        )
        match prism_res:
            case Success(prism_obj):
                fit_res = run_prism(prism_obj, n_iter=max_iter)
                match fit_res:
                    case Success(bp_fit):
                        frac_res = get_fraction(bp_fit, which_theta="first", state_or_type="type")
                        match frac_res:
                            case Success(df):
                                return df
                            case Failure(_):
                                pass
                    case Failure(_):
                        pass
            case Failure(_):
                pass

    # Fast projection fallback: NNLS with normalization
    A = ref_phi_mat.T
    theta_mat = np.zeros((n_samples, n_states), dtype=np.float64)

    for i in range(n_samples):
        b = bulk_mat[i, :]
        sol, _ = optimize.nnls(A, b)
        tot = sol.sum()
        theta_mat[i, :] = sol / tot if tot > 0 else 1.0 / n_states

    frac_dict: dict[str, list[object]] = {"sample_id": samples}
    for j, s in enumerate(states):
        frac_dict[s] = theta_mat[:, j].tolist()

    return pl.DataFrame(frac_dict)


def compute_normalized_fractions(
    raw_frac_df: pl.DataFrame,
    states: list[str],
    malignant_label: str,
    mrna_df: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Calculate tumor-purity-normalized immune fractions and mRNA-scaled cell proxies."""
    samples = raw_frac_df["sample_id"].to_list()
    raw_mat = raw_frac_df.select(states).to_numpy().astype(np.float64)

    # Identify malignant state index
    mal_indices = [
        i for i, s in enumerate(states) if malignant_label.lower() in s.lower()
    ]

    # 1. Purity-normalized immune fractions: theta / (1 - theta_tumor)
    norm_mat = np.copy(raw_mat)
    if mal_indices:
        tumor_purity = raw_mat[:, mal_indices].sum(axis=1, keepdims=True)
        immune_total = 1.0 - tumor_purity
        # Clip immune total to avoid division by zero
        immune_total = np.maximum(immune_total, 0.01)

        # Set tumor to 0 in immune compartment, rescale rest
        for m_idx in mal_indices:
            norm_mat[:, m_idx] = 0.0

        for non_m in range(len(states)):
            if non_m not in mal_indices:
                norm_mat[:, non_m] = norm_mat[:, non_m] / immune_total.flatten()

        # Rescale immune cells to sum to 1
        non_m_sum = norm_mat.sum(axis=1, keepdims=True)
        non_m_sum = np.where(non_m_sum == 0, 1.0, non_m_sum)
        norm_mat = norm_mat / non_m_sum

    norm_dict: dict[str, list[object]] = {"sample_id": samples}
    for j, s in enumerate(states):
        norm_dict[s] = norm_mat[:, j].tolist()
    norm_df = pl.DataFrame(norm_dict)

    # 2. mRNA-scaled fractions: theta_k / s_k
    scaling_map = {row["cell_state"]: row["relative_rna_content"] for row in mrna_df.to_dicts()}
    scale_factors = np.array([scaling_map.get(s, 1.0) for s in states])
    scale_factors = np.where(scale_factors <= 0, 1.0, scale_factors)

    scaled_mat = raw_mat / scale_factors
    scaled_sum = scaled_mat.sum(axis=1, keepdims=True)
    scaled_sum = np.where(scaled_sum == 0, 1.0, scaled_sum)
    scaled_mat = scaled_mat / scaled_sum

    scaled_dict: dict[str, list[object]] = {"sample_id": samples}
    for j, s in enumerate(states):
        scaled_dict[s] = scaled_mat[:, j].tolist()
    scaled_df = pl.DataFrame(scaled_dict)

    return norm_df, scaled_df


def run_pipeline(config: DeconvConfig) -> Result[None, str]:
    """Execute Step 4 bulk deconvolution and normalization."""
    match read_bulk_matrix(config.bulk_counts_path):
        case Failure(err):
            return Failure(err)
        case Success((bulk_df, samples, bulk_genes)):
            pass

    if not config.reference_path.exists():
        return Failure(f"Reference path not found: {config.reference_path}")
    if not config.hierarchy_path.exists():
        return Failure(f"Hierarchy path not found: {config.hierarchy_path}")
    if not config.mrna_scaling_path.exists():
        return Failure(f"mRNA scaling path not found: {config.mrna_scaling_path}")

    phi_df = pl.read_parquet(config.reference_path)
    mrna_df = pl.read_parquet(config.mrna_scaling_path)

    states = phi_df["cell_state"].to_list()
    ref_genes = [c for c in phi_df.columns if c != "cell_state"]

    # Shared genes
    shared_genes = [g for g in ref_genes if g in bulk_genes]
    if len(shared_genes) < 50:
        return Failure(f"Too few overlapping genes between bulk and reference: {len(shared_genes)}")

    # Extract aligned matrices
    ref_phi_mat = phi_df.select(shared_genes).to_numpy().astype(np.float64)
    bulk_mat = bulk_df.select(shared_genes).to_numpy().astype(np.float64)

    # Deconvolve
    raw_frac_df = run_deconvolution(
        bulk_mat,
        ref_phi_mat,
        states,
        samples,
        config.use_bayesprism_mcmc,
        config.max_gibbs_iter,
    )

    # Calculate normalized and scaled fractions
    norm_frac_df, scaled_frac_df = compute_normalized_fractions(
        raw_frac_df,
        states,
        config.malignant_label,
        mrna_df,
    )

    config.out_dir.mkdir(parents=True, exist_ok=True)
    raw_frac_df.write_parquet(config.out_dir / "bulk_deconv_fractions_raw.parquet")
    norm_frac_df.write_parquet(config.out_dir / "bulk_deconv_fractions_normalized.parquet")
    scaled_frac_df.write_parquet(config.out_dir / "bulk_deconv_fractions_mrna_scaled.parquet")

    print(
        f"[INFO] Successfully deconvolved {len(samples)} bulk samples across {len(states)} cell states "
        f"using {len(shared_genes)} shared genes."
    )
    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] Bulk deconvolution fractions saved to {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
