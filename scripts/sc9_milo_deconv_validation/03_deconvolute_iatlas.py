#!/usr/bin/env python3
"""Step 2: Deconvolute iAtlas Bulk RNA-seq Datasets using Multi-Resolution References.
Uses instaprism deconvolution to compute cell state proportions across immunotherapy cohorts.
Auto-discovers reference matrices and outputs deconv_fractions_{tag}.parquet.

Follows functional programming standards: pure functions, immutable data structures,
and explicit error handling using Returns Result monads.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import polars as pl
from joblib import Parallel, delayed  # type: ignore[import-untyped]
from pydantic import BaseModel, ConfigDict, Field
from returns.result import Failure, Result, Success

# Ensure repo packages are discoverable
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "tme_datasets" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "tme_datasets" / "src"))
if str(REPO_ROOT / "packages" / "instaprism" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "instaprism" / "src"))

import instaprism  # type: ignore[import-untyped]
from tme_datasets import load_dataset

DEFAULT_COHORTS: Final[tuple[str, ...]] = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
)

COHORT_CANCER_MAP: Final[dict[str, str]] = {
    "Hugo-iAtlas": "Melanoma",
    "Riaz-iAtlas": "Melanoma",
    "Liu-iAtlas": "Melanoma",
    "Gide-iAtlas": "Melanoma",
    "Rosenberg-iAtlas": "Bladder",
    "Padron-iAtlas": "Pancreatic",
    "Anders-iAtlas": "Breast",
    "McDermott-iAtlas": "Renal Cell",
    "Choueiri-iAtlas": "Renal Cell",
}


class DeconvConfig(BaseModel):
    """Immutable configuration for multi-cohort bulk deconvolution."""

    model_config = ConfigDict(frozen=True)

    ref_dir: Path
    out_dir: Path
    cohorts: tuple[str, ...] = DEFAULT_COHORTS
    n_iter: int = Field(default=50, ge=10)
    n_jobs: int = Field(default=-1)


def load_reference(ref_path: Path) -> Result[tuple[tuple[str, ...], tuple[str, ...], np.ndarray], str]:
    """Pure boundary to load deconvolution reference matrix from parquet."""
    try:
        df_ref = pl.read_parquet(ref_path)
        clusters = tuple(df_ref["cluster"].to_list())
        gene_cols = tuple(c for c in df_ref.columns if c != "cluster")
        ref_matrix = df_ref.select(gene_cols).to_numpy().astype(np.float64)
        return Success((clusters, gene_cols, ref_matrix))
    except Exception as exc:
        return Failure(f"Failed to load reference from {ref_path}: {exc}")


def load_cohort_bulk_expression(cohort_name: str) -> Result[pd.DataFrame, str]:
    """Load bulk RNA-seq expression (genes x samples, TPM) via tme_datasets."""
    match load_dataset(cohort_name, auto_download=False, normalize_ensembl=False):
        case Success(adata):
            try:
                x_raw = (
                    adata.X.toarray()
                    if hasattr(adata.X, "toarray")
                    else np.asarray(adata.X, dtype=np.float64)
                )
                x_mat = np.nan_to_num(
                    np.asarray(x_raw, dtype=np.float64),
                    nan=0.0,
                    posinf=0.0,
                    neginf=0.0,
                )
                raw_genes = (
                    adata.var["gene_name"].astype(str).str.strip().str.upper().tolist()
                    if "gene_name" in adata.var.columns
                    else [str(g).strip().upper() for g in adata.var_names]
                )
                sample_ids = [str(s) for s in adata.obs_names]
                df = pd.DataFrame(x_mat.T, index=raw_genes, columns=sample_ids)
                df = df[~df.index.isin(["", "NAN", "NONE", "<NA>", "UNKNOWN"])]
                df = df.groupby(level=0).mean()

                if float(df.to_numpy().max()) < 50.0:
                    df = np.maximum(np.power(2.0, np.clip(df, 0.0, 32.0)) - 1.0, 0.0)
                else:
                    df = np.maximum(df, 0.0)

                col_sums = df.sum(axis=0).replace(0.0, 1.0)
                df = df.div(col_sums, axis=1) * 1e6
                return Success(df)
            except Exception as exc:
                return Failure(f"Failed to convert bulk AnnData for {cohort_name}: {exc}")
        case Failure(err):
            return Failure(f"Could not load bulk cohort {cohort_name}: {err}")


def deconvolute_single_sample(
    bulk_vector: np.ndarray,
    norm_reference: np.ndarray,
    n_iter: int,
) -> np.ndarray:
    """Run instaprism on a single sample."""
    res = instaprism.insta_prism(
        bulk=bulk_vector,
        reference=norm_reference,
        n_iter=n_iter,
    )
    # res is (cell_state_gene_expression, probability_matrix, cell_fractions)
    cell_fracs = res[2] if len(res) >= 3 else res[0]
    return np.asarray(cell_fracs, dtype=np.float64)


def process_cohort_deconvolution(
    cohort: str,
    bulk_df: pd.DataFrame,
    clusters: tuple[str, ...],
    ref_genes: tuple[str, ...],
    ref_matrix: np.ndarray,
    n_iter: int,
    n_jobs: int,
) -> Result[pl.DataFrame, str]:
    """Pure wrapper executing InstaPrism on aligned bulk expression matrix."""
    try:
        ref_gene_set = set(ref_genes)
        common_genes = sorted(list(set(bulk_df.index).intersection(ref_gene_set)))
        if len(common_genes) < 50:
            return Failure(f"Too few overlapping genes ({len(common_genes)}) for {cohort}")

        bulk_sub = bulk_df.loc[common_genes].values.astype(np.float64)
        gene_indices = [ref_genes.index(g) for g in common_genes]
        phi_sub = ref_matrix[:, gene_indices]

        # Simplex re-normalization
        row_sums = phi_sub.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        norm_ref = phi_sub / row_sums

        sample_ids = list(bulk_df.columns)
        results: list[np.ndarray] = Parallel(n_jobs=n_jobs)(
            delayed(deconvolute_single_sample)(bulk_sub[:, s_idx], norm_ref, n_iter)
            for s_idx in range(len(sample_ids))
        )
        fraction_matrix = np.vstack(results)  # (n_samples, n_clusters)

        records: list[dict[str, object]] = []
        for s_idx, sample_id in enumerate(sample_ids):
            row_dict: dict[str, object] = {
                "sample_id": sample_id,
                "cohort": cohort,
                "cancer_type": COHORT_CANCER_MAP.get(cohort, "Other"),
            }
            for c_idx, cluster_name in enumerate(clusters):
                row_dict[cluster_name] = float(fraction_matrix[s_idx, c_idx])
            records.append(row_dict)

        return Success(pl.DataFrame(records))
    except Exception as exc:
        return Failure(f"Deconvolution failed for {cohort}: {exc}")


def run_deconvolution_pipeline(config: DeconvConfig) -> Result[Path, str]:
    """Discovers all reference parquets and deconvolutes all bulk cohorts."""
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Cache bulk expression
    print(f"Loading bulk RNA-seq cohorts ({len(config.cohorts)} cohorts)...")
    bulk_cache: dict[str, pd.DataFrame] = {}
    for c in config.cohorts:
        match load_cohort_bulk_expression(c):
            case Success(b_df):
                bulk_cache[c] = b_df
                print(f"  Loaded {c} ({b_df.shape[1]} samples, {b_df.shape[0]} genes)")
            case Failure(err):
                print(f"  Warning: {err}")

    if not bulk_cache:
        return Failure("No bulk cohorts could be loaded.")

    # 2. Discover reference parquet files
    all_ref_files = sorted(config.ref_dir.glob("reference_phi_*_res*.parquet"))
    ref_targets: list[tuple[str, Path, str]] = []

    for rf in all_ref_files:
        if "tidy" in rf.name:
            continue
        m = re.match(r"^reference_phi_(.+)_res([0-9\.]+)\.parquet$", rf.name)
        if m:
            d_id, res_str = m.group(1), m.group(2)
            tag = f"{d_id}_res{res_str}"
            out_file = f"deconv_fractions_{tag}.parquet"
            ref_targets.append((tag, rf, out_file))

    if not ref_targets:
        return Failure(f"No reference_phi_*_res*.parquet files found in {config.ref_dir}")

    print(f"\nDiscovered {len(ref_targets)} reference matrices for benchmarking.")
    last_out: Path | None = None

    for tag, ref_p, out_name in ref_targets:
        out_p = config.out_dir / out_name
        if out_p.exists() and out_p.stat().st_size > 1000:
            print(f"[{tag}] Cached deconvolution exists at {out_name}. Skipping.")
            last_out = out_p
            continue

        print(f"\n--- Deconvoluting against reference: {tag} ({ref_p.name}) ---")
        ref_res = load_reference(ref_p)
        match ref_res:
            case Failure(err):
                print(f"Error reading {ref_p}: {err}", file=sys.stderr)
                continue
            case Success((clusters, ref_genes, ref_matrix)):
                pass

        cohort_dfs: list[pl.DataFrame] = []
        for cohort, bulk_df in bulk_cache.items():
            d_res = process_cohort_deconvolution(
                cohort=cohort,
                bulk_df=bulk_df,
                clusters=clusters,
                ref_genes=ref_genes,
                ref_matrix=ref_matrix,
                n_iter=config.n_iter,
                n_jobs=config.n_jobs,
            )
            match d_res:
                case Failure(err):
                    print(f"  Warning: {cohort} failed for {tag}: {err}")
                case Success(df_cohort):
                    cohort_dfs.append(df_cohort)

        if cohort_dfs:
            master_df = pl.concat(cohort_dfs, how="vertical")
            master_df.write_parquet(out_p)
            print(f"[{tag}] Saved deconvolution fractions ({master_df.height} samples) -> {out_p}")
            last_out = out_p

    if last_out is not None:
        return Success(last_out)
    return Failure("Deconvolution failed across all reference targets.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 2: Deconvolute iAtlas bulk RNA-seq cohorts across multi-resolution references."
    )
    parser.add_argument(
        "--ref-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/01_references",
        help="Directory containing generated reference parquets",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/02_deconvolution",
        help="Directory to save deconvolution fractions",
    )
    parser.add_argument(
        "--n-iter",
        type=int,
        default=50,
        help="Number of InstaPrism iterations",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Parallel worker count (-1 = all available cores)",
    )
    args = parser.parse_args()

    config = DeconvConfig(
        ref_dir=Path(args.ref_dir),
        out_dir=Path(args.out_dir),
        n_iter=args.n_iter,
        n_jobs=args.n_jobs,
    )

    match run_deconvolution_pipeline(config):
        case Success(p):
            print(f"\nDeconvolution pipeline completed successfully: {p}")
            sys.exit(0)
        case Failure(err):
            print(f"\nDeconvolution pipeline failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
