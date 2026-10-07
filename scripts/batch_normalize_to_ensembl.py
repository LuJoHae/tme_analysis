#!/usr/bin/env python3
"""Batch Ensembl Normalization Pipeline for Preprocessed Single-Cell Datasets.

Scans cached .h5ad files, detects cohorts whose var_names are not yet canonical
Ensembl gene IDs (ENSG...), converts them to Ensembl Release 111 with persistent
Parquet caching and atomic disk writes.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gc
import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

from tme_datasets import Modality, list_datasets
from tme_datasets.logging import get_logger
from tme_datasets.paths import find_dataset_h5ad, get_data_paths
from tme_datasets.preprocessing.gene_normalization import (
    ensure_ensembl_release_installed,
    normalize_dataset_to_ensembl,
)

logger = get_logger("batch_normalize_ensembl")


class NormalizationResult(BaseModel):
    """Immutable record of cohort normalization outcome."""

    model_config = ConfigDict(frozen=True)

    dataset_id: str
    status: str
    n_obs: int
    n_vars_before: int
    n_vars_after: int
    elapsed_sec: float
    error: str = ""


def check_is_ensembl_indexed(h5ad_path: Path) -> bool:
    """Fast inspection via h5py to determine if var_names are already Ensembl IDs."""
    if not h5ad_path.is_file() or h5ad_path.stat().st_size == 0:
        return False
    try:
        with h5py.File(h5ad_path, "r") as f:
            if "var" not in f:
                return False
            var_grp = f["var"]
            idx_name = var_grp.attrs.get("_index", "_index")
            if idx_name in var_grp:
                idx_ds = var_grp[idx_name]
            elif "_index" in var_grp:
                idx_ds = var_grp["_index"]
            elif isinstance(var_grp, h5py.Dataset):
                idx_ds = var_grp
            else:
                return False
            sample_slice = idx_ds[:min(20, len(idx_ds))]
            names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in sample_slice]
            return len(names) > 0 and any(n.startswith("ENSG") for n in names)
    except Exception:
        return False


def normalize_single_h5ad(
    h5ad_path: Path,
    release: int = 111,
    species: str = "human",
    ensembl_dir: Path | None = None,
) -> NormalizationResult:
    """Normalize a single .h5ad file in-place using atomic staging."""
    dataset_id = h5ad_path.stem
    t0 = time.perf_counter()

    try:
        # 1. Read existing AnnData
        adata = ad.read_h5ad(h5ad_path)
        n_obs = adata.n_obs
        n_vars_before = adata.n_vars

        if n_vars_before == 0:
            return NormalizationResult(
                dataset_id=dataset_id,
                status="SKIPPED_EMPTY",
                n_obs=n_obs,
                n_vars_before=0,
                n_vars_after=0,
                elapsed_sec=round(time.perf_counter() - t0, 2),
                error="H5AD has 0 genes before normalization; skipping.",
            )

        # 2. Normalize to canonical Ensembl
        norm_res = normalize_dataset_to_ensembl(
            adata=adata,
            release=release,
            species=species,
            ensembl_dir=ensembl_dir,
            drop_unmapped=True,
        )

        match norm_res:
            case Failure(err):
                return NormalizationResult(
                    dataset_id=dataset_id,
                    status="FAILED",
                    n_obs=n_obs,
                    n_vars_before=n_vars_before,
                    n_vars_after=0,
                    elapsed_sec=round(time.perf_counter() - t0, 2),
                    error=str(err),
                )
            case Success(norm_adata):
                # Safety check: never write an empty feature matrix over an existing dataset!
                if norm_adata.n_vars == 0 or (n_vars_before >= 100 and norm_adata.n_vars < 100):
                    return NormalizationResult(
                        dataset_id=dataset_id,
                        status="FAILED_ZERO_GENES",
                        n_obs=n_obs,
                        n_vars_before=n_vars_before,
                        n_vars_after=norm_adata.n_vars,
                        elapsed_sec=round(time.perf_counter() - t0, 2),
                        error=f"Mapped to {norm_adata.n_vars} genes from {n_vars_before}. Aborting to protect data.",
                    )

                # 3. Sanitize index names for robust h5py serialization
                if norm_adata.obs_names.name is not None and norm_adata.obs_names.name in norm_adata.obs.columns:
                    norm_adata.obs_names.name = None
                if norm_adata.var_names.name is not None and norm_adata.var_names.name in norm_adata.var.columns:
                    norm_adata.var_names.name = None

                # 4. Atomic staging write
                staging_path = h5ad_path.with_suffix(".h5ad.tmp")
                if staging_path.exists():
                    staging_path.unlink()

                norm_adata.write_h5ad(staging_path)
                staging_path.replace(h5ad_path)

                n_vars_after = norm_adata.n_vars
                del norm_adata, adata
                gc.collect()

                elapsed = round(time.perf_counter() - t0, 2)
                logger.info(
                    "[%s] Successfully normalized: %d cells, %d -> %d genes (took %.2fs)",
                    dataset_id,
                    n_obs,
                    n_vars_before,
                    n_vars_after,
                    elapsed,
                )
                return NormalizationResult(
                    dataset_id=dataset_id,
                    status="SUCCESS",
                    n_obs=n_obs,
                    n_vars_before=n_vars_before,
                    n_vars_after=n_vars_after,
                    elapsed_sec=elapsed,
                )
    except Exception as exc:
        gc.collect()
        return NormalizationResult(
            dataset_id=dataset_id,
            status="FAILED",
            n_obs=0,
            n_vars_before=0,
            n_vars_after=0,
            elapsed_sec=round(time.perf_counter() - t0, 2),
            error=str(exc),
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch Ensembl normalization for preprocessed H5AD caches.")
    parser.add_argument("--workers", "-w", type=int, default=4, help="Concurrent workers (default: 4).")
    parser.add_argument("--release", "-r", type=int, default=111, help="Ensembl release (default: 111).")
    parser.add_argument("--force", "-f", action="store_true", help="Force re-normalization of all files.")
    parser.add_argument("--limit", "-l", type=int, default=None, help="Optional limit on number of cohorts to process.")
    parser.add_argument("--single-cell-only", action="store_true", default=True, help="Only process single-cell datasets.")
    parser.add_argument("--all-cohorts", action="store_true", help="Process all modalities including bulk.")
    args = parser.parse_args()

    paths_cfg = get_data_paths()
    preprocessed_dir = paths_cfg.preprocessed_dir
    ensembl_dir = paths_cfg.ensembl_dir

    print("=" * 70)
    print("Batch Ensembl Normalization Pipeline")
    print(f"Preprocessed Directory: {preprocessed_dir}")
    print(f"Ensembl Directory:      {ensembl_dir} (Release {args.release})")
    print(f"Workers:                {args.workers}")
    print("=" * 70)

    # 1. Warm-up Ensembl release and mapping cache
    print("--> Verifying Ensembl Release database...")
    ensure_ensembl_release_installed(release=args.release, ensembl_dir=ensembl_dir)
    print("--> Ensembl database verified.\n")

    # 2. Collect candidate H5AD files
    if args.all_cohorts:
        all_h5ad = sorted(preprocessed_dir.glob("*.h5ad"))
    else:
        # Default: single-cell datasets in registry
        sc_cohorts = list_datasets(modality=Modality.SINGLE_CELL)
        all_h5ad = []
        for ds in sc_cohorts:
            found = find_dataset_h5ad(ds.id)
            if found:
                all_h5ad.append(found.unwrap())

    candidates: list[Path] = []
    already_done = 0

    for p in all_h5ad:
        if p.name.endswith(".tmp"):
            continue
        if check_is_ensembl_indexed(p):
            already_done += 1
            continue
        try:
            h_backed = ad.read_h5ad(p, backed="r")
            if h_backed.n_vars > 0:
                candidates.append(p)
            else:
                print(f"Skipping {p.name}: 0 genes in cached file.")
        except Exception:
            pass

    if args.limit is not None:
        candidates = candidates[: args.limit]

    print(f"Total Target Cohorts:          {len(all_h5ad)}")
    print(f"Already Ensembl Normalized:    {already_done}")
    print(f"Pending Normalization:         {len(candidates)}\n")

    if not candidates:
        print("[DONE] All target preprocessed datasets already use canonical Ensembl gene IDs!")
        return 0

    # 3. Execute concurrently with ProcessPoolExecutor
    t_start = time.perf_counter()
    completed = 0
    failed = 0
    total = len(candidates)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_path = {
            executor.submit(
                normalize_single_h5ad,
                h5ad_path=p,
                release=args.release,
                ensembl_dir=ensembl_dir,
            ): p
            for p in candidates
        }

        for future in concurrent.futures.as_completed(future_to_path):
            p = future_to_path[future]
            try:
                res = future.result()
                if res.status == "SUCCESS":
                    completed += 1
                    print(
                        f"[{completed + failed}/{total}] SUCCESS: {res.dataset_id} "
                        f"({res.n_obs:,} cells, {res.n_vars_before} -> {res.n_vars_after} genes, {res.elapsed_sec:.1f}s)"
                    )
                else:
                    failed += 1
                    print(f"[{completed + failed}/{total}] {res.status}: {res.dataset_id} - Error: {res.error}")
            except Exception as exc:
                failed += 1
                print(f"[{completed + failed}/{total}] FAILED:  {p.stem} - Exception: {exc}")

    total_time = time.perf_counter() - t_start
    print("\n" + "=" * 70)
    print(f"Normalization Complete in {total_time / 60:.2f} minutes")
    print(f"Successfully Normalized: {completed}/{total}")
    print(f"Failed/Skipped:          {failed}/{total}")
    print("=" * 70)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
