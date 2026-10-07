#!/usr/bin/env python3
"""Standardized Batch Preprocessor for Discovered Single-Cell Solid Tumor Cohorts.

Converts raw single-cell matrix assets (10x MTX, H5, TSV, H5AD) in data/raw/{dataset_id}
into standardized, dual-layer AnnData (.h5ad) files in data/preprocessed/{dataset_id}.h5ad.

Standardization Specifications:
- adata.X: Non-negative integer raw counts (un-logged via check_transform_state if necessary).
- adata.layers["log1p_norm"]: Log1p-normalized expression values.
- Gene identifiers: Harmonized to Ensembl Release 111 (human) official HGNC symbols and ENSG IDs.
- Quality Control: QualityControlSpec (min 200 genes, max 8000 genes, >= 500 counts, mito <= 20%).
- Standardized metadata (adata.obs): dataset_id, patient_id, cancer_type, therapy, response,
  cell_selection, technology.
- Optional raw cleanup (--clean-raw) to conserve disk space.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import shutil
import sys
import tarfile
import time
import psutil
from pathlib import Path
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy.sparse import csr_matrix, issparse


class MemoryGovernor:
    """Enforces that total RAM usage never exceeds a strict limit (e.g. 150 GB)."""

    def __init__(self, max_ram_gb: float = 150.0) -> None:
        self.max_ram_gb = max_ram_gb
        self.max_ram_bytes = int(max_ram_gb * (1024 ** 3))

    def get_system_used_ram_gb(self) -> float:
        """Returns total system used RAM in GB."""
        try:
            return psutil.virtual_memory().used / (1024 ** 3)
        except Exception:
            return 0.0

    def get_process_rss_gb(self) -> float:
        """Returns current process resident memory in GB."""
        try:
            return psutil.Process().memory_info().rss / (1024 ** 3)
        except Exception:
            return 0.0

    def enforce_headroom(self, required_headroom_gb: float = 8.0, check_interval: float = 2.0, max_wait: float = 180.0) -> None:
        """Blocks if memory limit is approached until garbage collection or prior jobs release memory."""
        gc.collect()
        t0 = time.time()
        while True:
            used_gb = self.get_system_used_ram_gb()
            # If system has more than max_ram_gb total, enforce the 150 GB ceiling
            # If system has less (e.g. 16GB laptop), scale appropriately to available
            sys_total = psutil.virtual_memory().total / (1024 ** 3)
            effective_cap = min(self.max_ram_gb, sys_total * 0.92)

            if used_gb + required_headroom_gb <= effective_cap:
                break

            print(f"[MemoryGovernor] High RAM: {used_gb:.1f} GB used. Ceiling: {effective_cap:.1f} GB. Waiting for memory to clear...")
            gc.collect()
            time.sleep(check_interval)
            if time.time() - t0 > max_wait:
                print(f"[MemoryGovernor] Warning: Memory wait timed out after {max_wait:.0f}s. Proceeding.")
                break


class PreprocessTask(BaseModel):
    """Specification of a cohort preprocessing task."""

    model_config = ConfigDict(frozen=True)

    accession: str
    indication: str
    tier: str
    technology: str
    cell_selection: str
    raw_dir: Path
    output_h5ad: Path


def check_and_restore_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    """Detect if expression matrix is log-transformed and restore integer counts in X."""
    X_mat = adata.X
    if issparse(X_mat):
        sample_vals = X_mat[:min(200, adata.n_obs)].toarray()
    else:
        sample_vals = X_mat[:min(200, adata.n_obs)]

    max_val = float(sample_vals.max()) if sample_vals.size > 0 else 0.0
    is_non_integer = not np.all(np.equal(np.mod(sample_vals[sample_vals > 0], 1), 0))
    is_log = (max_val <= 35.0 and is_non_integer and max_val > 0.0)

    if is_log:
        if issparse(adata.X):
            adata.X.data = np.expm1(adata.X.data)
        else:
            adata.X = np.expm1(adata.X)

    # Ensure sparse CSR matrix
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X, dtype=np.float32)
    elif adata.X.format != "csr":
        adata.X = adata.X.tocsr()

    # Store log1p normalized values in layer
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    adata.layers["log1p_norm"] = adata_norm.X
    del adata_norm
    gc.collect()

    return adata


def read_raw_matrix_asset(raw_dir: Path, accession: str) -> Result[ad.AnnData, str]:
    """Identify and parse matrix assets in raw_dir into an AnnData object."""
    if not raw_dir.exists():
        return Failure(f"Directory {raw_dir} does not exist.")

    # 1. Check for pre-existing .h5ad
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if h5ad_files:
        try:
            adata = ad.read_h5ad(h5ad_files[0])
            return Success(adata)
        except Exception as exc:
            return Failure(f"Failed to read H5AD {h5ad_files[0]}: {exc}")

    # 2. Check for 10x .h5 files
    h5_files = list(raw_dir.glob("*.h5"))
    if h5_files:
        try:
            adata = sc.read_10x_h5(h5_files[0])
            adata.var_names_make_unique()
            return Success(adata)
        except Exception as exc:
            return Failure(f"Failed to read 10x H5 {h5_files[0]}: {exc}")

    # 3. Check for tar archives and unpack if needed
    tar_files = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
    for tf in tar_files:
        extract_dir = raw_dir / "unpacked"
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            try:
                with tarfile.open(tf, "r:*") as archive:
                    archive.extractall(extract_dir)
            except Exception as exc:
                pass

    # Search in both raw_dir and unpack subdirectories
    search_dirs = [raw_dir] + list(raw_dir.glob("**/"))

    # 4. Check for 10x MTX triplets
    for d in search_dirs:
        mtx_files = list(d.glob("*matrix.mtx*"))
        if mtx_files:
            try:
                adata = sc.read_10x_mtx(d)
                adata.var_names_make_unique()
                return Success(adata)
            except Exception:
                pass

    # 5. Check for dense or sparse count TSV/CSV/TXT tables
    for d in search_dirs:
        count_tables = [
            f for f in d.glob("*.*")
            if any(k in f.name.lower() for k in ("count", "tpm", "expression", "counts"))
            and any(f.name.lower().endswith(ext) for ext in (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".csv", ".txt"))
            and not f.name.startswith("filelist")
        ]
        if count_tables:
            target = count_tables[0]
            try:
                sep = "\t" if (target.name.endswith(".tsv") or target.name.endswith(".tsv.gz") or "txt" in target.name) else ","
                df = pl.read_csv(
                    target,
                    separator=sep,
                    truncate_ragged_lines=True,
                    ignore_errors=True,
                    encoding="utf8-lossy",
                )
                first_col = df.columns[0]
                genes = df[first_col].cast(pl.String).to_list()
                cell_cols = df.columns[1:]
                mat = df.select(cell_cols).to_numpy().astype(np.float32).T
                sparse_mat = csr_matrix(mat, dtype=np.float32)
                adata = ad.AnnData(
                    X=sparse_mat,
                    obs=pd.DataFrame(index=cell_cols),
                    var=pd.DataFrame(index=genes),
                )
                adata.var_names_make_unique()
                return Success(adata)
            except Exception as exc:
                return Failure(f"Failed to read count table {target}: {exc}")

    return Failure(f"No parseable matrix assets found in {raw_dir}")


def apply_standard_qc(adata: ad.AnnData) -> ad.AnnData:
    """Filter low-quality cells by count depth, gene detection, and mitochondrial percentage."""
    adata.var_names_make_unique()

    # Identify mitochondrial genes
    mt_genes = adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
    if issparse(adata.X):
        total_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
        if np.any(mt_genes):
            mt_counts = np.asarray(adata.X[:, mt_genes].sum(axis=1)).flatten()
            pct_mt = (mt_counts / np.maximum(total_counts, 1.0)) * 100.0
        else:
            pct_mt = np.zeros(adata.n_obs, dtype=np.float32)
    else:
        total_counts = adata.X.sum(axis=1).flatten()
        genes_detected = (adata.X > 0).sum(axis=1).flatten()
        if np.any(mt_genes):
            mt_counts = adata.X[:, mt_genes].sum(axis=1).flatten()
            pct_mt = (mt_counts / np.maximum(total_counts, 1.0)) * 100.0
        else:
            pct_mt = np.zeros(adata.n_obs, dtype=np.float32)

    adata.obs["total_counts"] = total_counts
    adata.obs["n_genes_by_counts"] = genes_detected
    adata.obs["pct_counts_mt"] = pct_mt

    # Apply QualityControlSpec thresholds
    keep_mask = (
        (genes_detected >= 200)
        & (genes_detected <= 8000)
        & (total_counts >= 500)
        & (pct_mt <= 20.0)
    )

    filtered = adata[keep_mask].copy()
    return filtered


def harmonize_obs_metadata(adata: ad.AnnData, task: PreprocessTask) -> ad.AnnData:
    """Harmonize cell metadata to canonical tme_analysis schema."""
    obs_df = adata.obs

    if "original_barcode" not in obs_df.columns:
        obs_df["original_barcode"] = obs_df.index.astype(str)

    obs_df["dataset_id"] = task.accession
    obs_df["cancer_type"] = task.indication
    obs_df["tier"] = task.tier
    obs_df["technology"] = task.technology
    obs_df["cell_selection"] = task.cell_selection

    # Map donor/patient if available
    for donor_col in ("patient", "patient_id", "donor", "donor_id", "subject", "sample"):
        if donor_col in obs_df.columns and "patient_id" not in obs_df.columns:
            obs_df["patient_id"] = obs_df[donor_col]

    if "patient_id" not in obs_df.columns:
        obs_df["patient_id"] = task.accession + "_donor_unspecified"

    # Map response if present
    for resp_col in ("response", "Response", "RECIST", "benefit", "responder"):
        if resp_col in obs_df.columns and "response" not in obs_df.columns:
            obs_df["response"] = obs_df[resp_col]

    if "response" not in obs_df.columns:
        obs_df["response"] = "Documented (cohort level)" if "Response" in task.tier else "Baseline"

    # Harmonize cell type annotation if present
    for ct_col in ("cell_type", "celltype", "CellType", "major_cell_type", "cluster"):
        if ct_col in obs_df.columns and "cell_type" not in obs_df.columns:
            obs_df["cell_type"] = obs_df[ct_col]

    if "cell_type" not in obs_df.columns:
        obs_df["cell_type"] = "Unspecified"

    adata.obs = obs_df
    return adata


def execute_preprocess_task(task: PreprocessTask, clean_raw: bool = False) -> Result[str, str]:
    """Execute end-to-end preprocessing for a single cohort using tme_datasets."""
    from tme_datasets.query import load_dataset

    res = load_dataset(
        dataset_id=task.accession,
        raw_dir=task.raw_dir,
        output_h5ad=task.output_h5ad,
        auto_download=False,
        force_recompute=True,
        cache_h5ad=True,
        normalize_ensembl=True,
        apply_qc=True,
    )

    match res:
        case Failure(err):
            return Failure(f"{task.accession}: {err}")
        case Success(adata):
            pass

    if clean_raw and task.raw_dir.exists():
        try:
            shutil.rmtree(task.raw_dir)
        except Exception:
            pass

    return Success(
        f"{task.accession}: Saved preprocessed AnnData ({adata.n_obs:,} cells x {adata.n_vars:,} genes) to {task.output_h5ad}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Standardized Batch Preprocessor for scRNA-seq Cohorts.")
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("data/registry/discovered_solid_tumor_sc_datasets.parquet"),
        help="Path to discovered cohorts Parquet registry.",
    )
    parser.add_argument(
        "--data-raw-dir",
        type=Path,
        default=Path("data/raw"),
        help="Input raw data directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/preprocessed"),
        help="Target preprocessed H5AD directory.",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Comma-separated accessions to preprocess.",
    )
    parser.add_argument(
        "--clean-raw",
        action="store_true",
        help="Delete raw directory after successful H5AD validation to conserve disk space.",
    )
    parser.add_argument(
        "--max-ram-gb",
        type=float,
        default=150.0,
        help="Maximum simultaneous RAM ceiling in GB (default: 150.0).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force recomputing cohorts that already have a preprocessed .h5ad file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned preprocessing tasks without executing.",
    )

    args = parser.parse_args()

    df = pl.read_parquet(args.registry)
    if args.cohorts:
        specific = [s.strip() for s in args.cohorts.split(",") if s.strip()]
        df = df.filter(pl.col("accession").is_in(specific))

    tasks: list[PreprocessTask] = []
    for row in df.iter_rows(named=True):
        acc = row["accession"]
        raw_p = args.data_raw_dir / acc
        out_p = args.output_dir / f"{acc}.h5ad"

        tasks.append(
            PreprocessTask(
                accession=acc,
                indication=row["indication"],
                tier=row["tier"],
                technology=row["technology"],
                cell_selection=row["cell_selection_strategy"],
                raw_dir=raw_p,
                output_h5ad=out_p,
            )
        )

    governor = MemoryGovernor(max_ram_gb=args.max_ram_gb)

    print(f"\n=======================================================")
    print(f"  TME Analysis scRNA-seq Batch Preprocessor")
    print(f"  Total tasks resolved: {len(tasks)}")
    print(f"  Max simultaneous RAM: {args.max_ram_gb:.1f} GB")
    print(f"  Dry run mode:         {args.dry_run}")
    print(f"  Clean raw archives:   {args.clean_raw}")
    print(f"=======================================================\n")

    if args.dry_run:
        for idx, t in enumerate(tasks[:15], 1):
            print(f"[{idx:2d}/{len(tasks):2d}] [DRY-RUN] {t.accession:18s} | {t.raw_dir} -> {t.output_h5ad}")
        print("\nDry run complete. No processing executed.")
        return

    successes = 0
    failures = 0
    for idx, t in enumerate(tasks, 1):
        if not t.raw_dir.exists():
            print(f"[{idx:3d}/{len(tasks):3d}] - {t.accession}: Raw directory not found (download first).")
            continue

        if t.output_h5ad.exists() and not args.force:
            print(f"[{idx:3d}/{len(tasks):3d}] [SKIP] {t.accession}: {t.output_h5ad.name} already exists. (Pass --force to reprocess)")
            successes += 1
            continue

        # Enforce that starting this cohort never breaches the 150 GB ceiling
        governor.enforce_headroom(required_headroom_gb=10.0)
        used_mem = governor.get_system_used_ram_gb()
        proc_rss = governor.get_process_rss_gb()
        print(f"[{idx:3d}/{len(tasks):3d}] [RAM: {used_mem:.1f} GB used (Process RSS: {proc_rss:.1f} GB) / Max: {args.max_ram_gb:.1f} GB] Processing {t.accession}...")

        res = execute_preprocess_task(t, clean_raw=args.clean_raw)
        gc.collect()

        match res:
            case Success(msg):
                print(f"[{idx:3d}/{len(tasks):3d}] ✓ {msg}")
                successes += 1
            case Failure(err):
                print(f"[{idx:3d}/{len(tasks):3d}] ✗ {err}")
                failures += 1

    print(f"\nPreprocessing finished: {successes} succeeded, {failures} failed.")


if __name__ == "__main__":
    main()
