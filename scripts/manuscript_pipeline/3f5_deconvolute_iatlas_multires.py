#!/usr/bin/env python3
"""
Step 2: Deconvolute iAtlas Bulk RNA-seq Datasets using Sade-Feldman Reference.
Uses instaprism deconvolution to compute cell state proportions across immunotherapy cohorts.
Outputs deconv_fractions.parquet.
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
from joblib import Parallel, delayed  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

# Local repo packages
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "packages"))
from scripts.manuscript_pipeline._helpers.paths import (  # noqa: E402
    MULTIRES_REF_STAGES,
    discover_upstream_glob,
    resolve_upstream_artifact,
)
import instaprism  # type: ignore
from tme_datasets import load_dataset  # noqa: E402


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
    model_config = ConfigDict(frozen=True)
    reference_path: Path
    lair_dir: Path | None = None
    cohorts: tuple[str, ...]
    out_dir: Path
    output_name: str = "deconv_fractions.parquet"
    benchmark_mode: bool = False
    n_iter: int = 50
    n_jobs: int = -1


def load_reference(ref_path: Path) -> Result[tuple[tuple[str, ...], tuple[str, ...], np.ndarray], str]:
    """Pure boundary to load deconvolution reference matrix from parquet."""
    try:
        df_ref = pl.read_parquet(ref_path)
        clusters = tuple(df_ref["cluster"].to_list())
        gene_cols = tuple(c for c in df_ref.columns if c != "cluster")
        # Extract matrix: (n_clusters, n_genes)
        ref_matrix = df_ref.select(gene_cols).to_numpy().astype(np.float64)
        return Success((clusters, gene_cols, ref_matrix))
    except Exception as exc:
        return Failure(f"Failed to load reference from {ref_path}: {exc}")


def load_cohort_bulk_expression(
    cohort_name: str,
    lair_dir: Path | None = None,
) -> Result[pd.DataFrame, str]:
    """Load bulk RNA-seq expression (genes x samples, 1e6 TPM) via tme_datasets H5AD."""
    match load_dataset(cohort_name, auto_download=False):
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
                return Failure(f"Error formatting tme_datasets H5AD for {cohort_name}: {exc}")
        case Failure(tme_err):
            _ = lair_dir
            return Failure(f"Failed to load {cohort_name} via tme_datasets.load_dataset(): {tme_err}")


def deconvolute_single_sample(
    bulk_vector: np.ndarray,
    norm_reference: np.ndarray,
    n_iter: int,
) -> np.ndarray:
    """Run instaprism on a single sample."""
    # instaprism returns (probability_matrix, cell_state_gene_expression, cell_fractions, ...)
    _, _, cell_fracs, _ = instaprism.insta_prism(
        bulk=bulk_vector,
        reference=norm_reference,
        n_iter=n_iter,
    )
    return cell_fracs


def process_cohort_deconvolution(
    cohort: str,
    bulk_df: pd.DataFrame,
    clusters: tuple[str, ...],
    ref_genes: tuple[str, ...],
    ref_matrix: np.ndarray,
    n_iter: int,
    n_jobs: int,
) -> Result[pl.DataFrame, str]:
    """Deconvolute a single cohort and return a Polars DataFrame of fractions."""
    common_genes = sorted(list(set(bulk_df.index).intersection(set(ref_genes))))
    if len(common_genes) < 100:
        return Failure(
            f"Insufficient common genes between {cohort} and reference ({len(common_genes)} < 100)"
        )

    print(f"[{cohort}] Deconvoluting {bulk_df.shape[1]} samples with {len(common_genes)} common genes...")

    # Gene index map for reference
    gene_to_idx = {g: i for i, g in enumerate(ref_genes)}
    ref_gene_indices = [gene_to_idx[g] for g in common_genes]

    # Subset reference and normalize rows to 1
    sub_ref = ref_matrix[:, ref_gene_indices].copy()
    row_sums = sub_ref.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    norm_ref = sub_ref / row_sums

    # Subset bulk
    sub_bulk = bulk_df.loc[common_genes].values.astype(np.float64)
    sample_ids = list(bulk_df.columns)

    # Parallel deconvolution
    results: list[np.ndarray] = Parallel(n_jobs=n_jobs)(
        delayed(deconvolute_single_sample)(sub_bulk[:, s_idx], norm_ref, n_iter)
        for s_idx in range(len(sample_ids))
    )

    fraction_matrix = np.vstack(results)  # (n_samples, n_clusters)
    cancer_type = COHORT_CANCER_MAP.get(cohort, "Other")

    records: list[dict[str, object]] = []
    for s_idx, sid in enumerate(sample_ids):
        rec: dict[str, object] = {
            "sample_id": str(sid),
            "cohort": cohort,
            "cancer_type": cancer_type,
        }
        for c_idx, cl in enumerate(clusters):
            rec[cl] = float(fraction_matrix[s_idx, c_idx])
        records.append(rec)

    return Success(pl.DataFrame(records))


def run_deconv_pipeline(config: DeconvConfig) -> Result[Path, str]:
    """Orchestrates deconvolution across cohorts, supporting single reference or multi-resolution benchmark mode."""
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # Preload bulk expression for all valid cohorts once to maximize efficiency
    print("Loading bulk RNA-seq cohorts via tme_datasets...")
    bulk_cache: dict[str, pd.DataFrame] = {}
    for cohort in config.cohorts:
        bulk_res = load_cohort_bulk_expression(cohort, config.lair_dir)
        match bulk_res:
            case Failure(err):
                print(f"Warning: Skipping {cohort}: {err}")
            case Success(bulk_df):
                bulk_cache[cohort] = bulk_df
                print(f"Cached {cohort} bulk expression ({bulk_df.shape[1]} samples, {bulk_df.shape[0]} genes).")

    if not bulk_cache:
        return Failure("No bulk cohorts could be loaded.")

    if config.benchmark_mode:
        print("\n=== RUNNING MULTI-RESOLUTION BENCHMARK DECONVOLUTION ===")
        # Dynamically identify all generated multi-resolution reference parquet files across 3f1..3f4
        ref_targets: list[tuple[str, Path, str]] = []
        explicit_ref_dirs = (config.out_dir, config.reference_path.parent)
        fallback_ref_dirs = (
            Path("output/output/sade_feldman_deconv_validation"),
            Path("output/sade_feldman_deconv_validation"),
        )

        # 1. Standard Sade-Feldman & Integrated Multi-Atlas
        for res in (0.5, 1.0, 1.5, 2.0):
            match resolve_upstream_artifact(
                f"reference_phi_res{res}.parquet",
                MULTIRES_REF_STAGES,
                explicit_dirs=explicit_ref_dirs,
                fallback_dirs=fallback_ref_dirs,
            ):
                case Some(sf_file):
                    ref_targets.append((f"sf_res{res}", sf_file, f"deconv_fractions_sf_res{res}.parquet"))
                case Nothing:
                    pass
            match resolve_upstream_artifact(
                f"integrated_reference_phi_res{res}.parquet",
                MULTIRES_REF_STAGES,
                explicit_dirs=explicit_ref_dirs,
                fallback_dirs=fallback_ref_dirs,
            ):
                case Some(comb_file):
                    ref_targets.append((f"comb_res{res}", comb_file, f"deconv_fractions_comb_res{res}.parquet"))
                case Nothing:
                    pass

        # 2. Any additional dataset/combined references: reference_phi_{dataset_id}_res{res}.parquet
        for phi_file in discover_upstream_glob(
            "reference_phi_*_res*.parquet",
            MULTIRES_REF_STAGES,
            explicit_dirs=explicit_ref_dirs,
            fallback_dirs=fallback_ref_dirs,
        ):
            fname = phi_file.name
            # parse {dataset_id} and {res}, explicitly ignoring tidy files
            m = re.match(r"^reference_phi_(?!tidy_)([a-zA-Z0-9_\.-]+)_res([0-9\.]+)\.parquet$", fname)
            if m:
                d_id, r_val = m.group(1), m.group(2)
                tag = f"{d_id}_res{r_val}"
                out_name = f"deconv_fractions_{d_id}_res{r_val}.parquet"
                if not any(t[0] == tag for t in ref_targets):
                    ref_targets.append((tag, phi_file, out_name))

        if not ref_targets:
            # Check fallback to primary reference_phi.parquet
            if config.reference_path.exists():
                ref_targets.append(("primary", config.reference_path, "deconv_fractions.parquet"))
            else:
                return Failure("No reference parquet files found in out-dir for benchmarking.")

        print(f"Found {len(ref_targets)} reference configurations for deconvolution benchmarking.")
        last_out: Path = config.out_dir / "deconv_fractions.parquet"
        for tag, ref_path, out_name in ref_targets:
            out_p = config.out_dir / out_name
            if out_p.exists() and out_p.stat().st_size > 1000:
                print(f"[{tag}] Fractions already exist at {out_p.name} (cached). Skipping deconvolution.")
                last_out = out_p
                continue

            print(f"\n--- Deconvoluting against reference: {tag} ({ref_path.name}) ---")
            ref_res = load_reference(ref_path)
            match ref_res:
                case Failure(err):
                    print(f"Error loading {ref_path}: {err}")
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
                        print(f"Warning: {cohort} failed for {tag}: {err}")
                    case Success(df_cohort):
                        cohort_dfs.append(df_cohort)

            if cohort_dfs:
                master_df = pl.concat(cohort_dfs, how="vertical")
                out_p = config.out_dir / out_name
                master_df.write_parquet(out_p)
                print(f"[{tag}] Saved deconvolution fractions ({master_df.height} samples) to: {out_p}")
                last_out = out_p

                # Backward compatibility
                if tag == "sf_res0.5":
                    master_df.write_parquet(config.out_dir / "deconv_fractions.parquet")
                    print("Saved default deconv_fractions.parquet (from sf_res0.5).")

        return Success(last_out)

    else:
        # Single reference execution
        ref_res = load_reference(config.reference_path)
        match ref_res:
            case Failure(err):
                return Failure(err)
            case Success((clusters, ref_genes, ref_matrix)):
                pass

        print(f"Reference loaded: {len(clusters)} clusters x {len(ref_genes)} genes.")
        cohort_dfs = []
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
                    print(f"Warning: {cohort} failed deconvolution: {err}")
                case Success(df_cohort):
                    cohort_dfs.append(df_cohort)
                    print(f"[{cohort}] Successfully deconvoluted {df_cohort.height} samples.")

        if not cohort_dfs:
            return Failure("No cohorts were successfully deconvoluted.")

        master_df = pl.concat(cohort_dfs, how="vertical")
        out_path = config.out_dir / config.output_name
        master_df.write_parquet(out_path)
        print(f"\nCompleted deconvolution across {len(cohort_dfs)} cohorts.")
        print(f"Total samples: {master_df.height}. Output written to: {out_path}")
        return Success(out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 2: Deconvolute iAtlas cohorts across single or multi-resolution references."
    )
    parser.add_argument(
        "--reference",
        type=str,
        default="output/manuscript_pipeline/3f1_build_sf_reference/reference_phi.parquet",
        help="Path to reference_phi.parquet",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="",
        help="Unused deprecated CLI flag retained for backward compatibility",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="Hugo-iAtlas,Riaz-iAtlas,Liu-iAtlas,Gide-iAtlas,Rosenberg-iAtlas,Padron-iAtlas,Anders-iAtlas,McDermott-iAtlas,Choueiri-iAtlas",
        help="Comma-separated cohort names",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/manuscript_pipeline/3f5_deconvolute_iatlas_multires",
        help="Directory to save deconv fractions",
    )
    parser.add_argument(
        "--output-name",
        type=str,
        default="deconv_fractions.parquet",
        help="Filename for single output deconvolution parquet",
    )
    parser.add_argument(
        "--benchmark-mode",
        action="store_true",
        help="Execute multi-resolution deconvolution across all 8 reference variants",
    )
    parser.add_argument(
        "--n-iter",
        type=int,
        default=50,
        help="Number of instaprism iterations",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Number of parallel worker processes",
    )
    args = parser.parse_args()

    cohort_list = tuple(c.strip() for c in args.cohorts.split(",") if c.strip())
    config = DeconvConfig(
        reference_path=Path(args.reference),
        lair_dir=Path(args.lair_dir) if args.lair_dir else None,
        cohorts=cohort_list,
        out_dir=Path(args.out_dir),
        output_name=args.output_name,
        benchmark_mode=args.benchmark_mode,
        n_iter=args.n_iter,
        n_jobs=args.n_jobs,
    )

    match run_deconv_pipeline(config):
        case Success(out_file):
            print(f"Step 2 finished successfully: {out_file}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 2 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
