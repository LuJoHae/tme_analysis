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
from returns.result import Failure, Result, Success

# Local repo packages
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "packages"))
import datalair  # type: ignore
import ici_datasets  # type: ignore
import instaprism  # type: ignore


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
    lair_dir: Path
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


def load_cohort_bulk_expression(cohort_name: str, lair: datalair.Lair) -> Result[pd.DataFrame, str]:
    """Load bulk RNA-seq expression for a cohort from cBioPortal via datalair."""
    try:
        ds_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
        ds = ds_class(name=cohort_name)
        lair.safe_derive(ds)
        filepaths = lair.get_dataset_filepaths(ds)

        unpacked_key = next(f for f in filepaths.keys() if not f.endswith(".tar.gz"))
        p_dir = filepaths[unpacked_key] / filepaths[unpacked_key].name

        dispatch_map = {
            "data_mrna_seq_expression.txt": True,
            "data_mrna_seq_tpm.txt": True,
            "data_mrna_seq_rpkm.txt": True,
        }
        existing = [f for f in dispatch_map if (p_dir / f).is_file()]
        if not existing:
            return Failure(f"No bulk mRNA expression file found in {p_dir} for {cohort_name}")

        target = existing[0]
        df = pd.read_csv(p_dir / target, sep="\t", index_col=0)

        # Normalize index to uppercase gene symbols
        df.index = df.index.astype(str).str.upper()
        # Deduplicate multiple rows with same gene by averaging
        df = df.groupby(level=0).mean()

        # If log-transformed (max value < 50), linearize
        if df.max().max() < 50:
            df = np.power(2.0, df) - 1.0

        # Size-factor normalize to sum to 1e6 (TPM-like)
        col_sums = df.sum(axis=0)
        col_sums[col_sums == 0] = 1.0
        df = df.div(col_sums, axis=1) * 1e6

        return Success(df)
    except Exception as exc:
        return Failure(f"Error loading bulk data for {cohort_name}: {exc}")


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
    lair = datalair.Lair(str(config.lair_dir))

    # Preload bulk expression for all valid cohorts once to maximize efficiency
    print("Loading bulk RNA-seq cohorts from cBioPortal / datalair...")
    bulk_cache: dict[str, pd.DataFrame] = {}
    for cohort in config.cohorts:
        bulk_res = load_cohort_bulk_expression(cohort, lair)
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
        # Dynamically identify all generated multi-resolution reference parquet files
        ref_targets: list[tuple[str, Path, str]] = []
        
        # 1. Standard Sade-Feldman & Integrated Multi-Atlas
        for res in (0.5, 1.0, 1.5, 2.0):
            sf_file = config.out_dir / f"reference_phi_res{res}.parquet"
            if sf_file.exists():
                ref_targets.append((f"sf_res{res}", sf_file, f"deconv_fractions_sf_res{res}.parquet"))
            comb_file = config.out_dir / f"integrated_reference_phi_res{res}.parquet"
            if comb_file.exists():
                ref_targets.append((f"comb_res{res}", comb_file, f"deconv_fractions_comb_res{res}.parquet"))

        # 2. Any additional dataset references: reference_phi_{dataset_id}_res{res}.parquet
        for phi_file in sorted(config.out_dir.glob("reference_phi_*_res*.parquet")):
            fname = phi_file.name
            # parse {dataset_id} and {res}, explicitly ignoring tidy files
            m = re.match(r"^reference_phi_(?!tidy_)([a-zA-Z0-9_\.-]+)_res([0-9\.]+)\.parquet$", fname)
            if m:
                d_id, r_val = m.group(1), m.group(2)
                tag = f"{d_id}_res{r_val}"
                out_name = f"deconv_fractions_{d_id}_res{r_val}.parquet"
                if (tag, phi_file, out_name) not in ref_targets:
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
        default="output/sade_feldman_deconv_validation/reference_phi.parquet",
        help="Path to reference_phi.parquet",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="scratch/lair" if Path("scratch/lair").exists() else "/storage/halu/lair",
        help="Path to datalair directory",
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
        default="output/sade_feldman_deconv_validation",
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
        lair_dir=Path(args.lair_dir),
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
