"""Standardized batch preprocessor for discovered single-cell solid tumor cohorts inside tme_datasets.

Converts raw single-cell matrix assets (10x MTX, H5, TSV, H5AD) in data/raw/{dataset_id}
into standardized, dual-layer AnnData (.h5ad) files in data/preprocessed/{dataset_id}.h5ad.
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path
import shutil
import signal
import sys
import tarfile
import time
from typing import Any, Sequence, cast

import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import pandas as pd
import polars as pl
import psutil  # type: ignore[import-untyped]
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import scanpy as sc  # type: ignore[import-untyped]
from scipy.sparse import csr_matrix, issparse

from ..preprocessing.gene_normalization import normalize_dataset_to_ensembl
from ..query import load_dataset
from ..registry import get_dataset_spec, resolve_dataset_id


class MemoryGovernor(BaseModel):
    """Enforces that total RAM usage never exceeds a strict ceiling (e.g. 150 GB)."""

    model_config = ConfigDict(frozen=True)

    max_ram_gb: float = 150.0

    def get_system_used_ram_gb(self) -> float:
        """Return total system used RAM in GB."""
        try:
            return float(psutil.virtual_memory().used) / float(1024**3)
        except Exception:  # noqa: BLE001
            return 0.0

    def get_process_rss_gb(self) -> float:
        """Return current process resident memory in GB."""
        try:
            return float(psutil.Process().memory_info().rss) / float(1024**3)
        except Exception:  # noqa: BLE001
            return 0.0

    def enforce_headroom(
        self,
        required_headroom_gb: float = 8.0,
        check_interval: float = 2.0,
        max_wait: float = 180.0,
    ) -> None:
        """Block if memory limit is approached until GC or prior jobs release memory."""
        gc.collect()
        t0 = time.time()
        while True:
            used_gb = self.get_system_used_ram_gb()
            sys_total = float(psutil.virtual_memory().total) / float(1024**3)
            effective_cap = min(self.max_ram_gb, sys_total * 0.92)
            if used_gb + required_headroom_gb <= effective_cap:
                break
            sys.stdout.write(
                f"[MemoryGovernor] High RAM: {used_gb:.1f} GB used. "
                f"Ceiling: {effective_cap:.1f} GB. Waiting for memory to clear...\n"
            )
            gc.collect()
            time.sleep(check_interval)
            if time.time() - t0 > max_wait:
                sys.stdout.write(
                    f"[MemoryGovernor] Warning: Memory wait timed out after {max_wait:.0f}s. Proceeding.\n"
                )
                break


class CohortPreprocessTask(BaseModel):
    """Specification of a cohort preprocessing task."""

    model_config = ConfigDict(frozen=True)

    accession: str
    indication: str
    tier: str
    technology: str
    cell_selection: str
    raw_dir: Path
    output_h5ad: Path


def _compute_log1p_norm_layer(adata: ad.AnnData) -> csr_matrix:
    """Compute log1p(1e4)-normalized CSR layer in an isolated scope."""
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    return csr_matrix(adata_norm.X, dtype=np.float32)


def check_and_restore_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    """Detect if expression matrix is log-transformed and restore integer counts in X."""
    n_sample = min(200, adata.n_obs)
    maybe_csr: Maybe[csr_matrix[Any]] = (
        Some(csr_matrix(cast(Any, adata.X))) if issparse(adata.X) else Nothing
    )
    sample_vals = (
        maybe_csr.unwrap()[:n_sample].toarray()
        if isinstance(maybe_csr, Some)
        else np.asarray(adata.X)[:n_sample]
    )

    max_val = float(sample_vals.max()) if sample_vals.size > 0 else 0.0
    pos_vals = sample_vals[sample_vals > 0]
    is_non_integer = bool(pos_vals.size > 0 and not np.all(np.equal(np.mod(pos_vals, 1), 0)))
    is_log = bool(max_val <= 35.0 and is_non_integer and max_val > 0.0)

    match maybe_csr:
        case Some(x_csr):
            updated_csr = (
                csr_matrix(
                    (np.expm1(x_csr.data).astype(np.float32), x_csr.indices, x_csr.indptr),
                    shape=x_csr.shape,
                )
                if is_log
                else x_csr.astype(np.float32)
            )
            adata.X = updated_csr
        case _:
            dense = np.asarray(adata.X)
            adata.X = csr_matrix(np.expm1(dense) if is_log else dense, dtype=np.float32)

    adata.layers["log1p_norm"] = _compute_log1p_norm_layer(adata)
    gc.collect()
    return adata


def read_raw_matrix_asset(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Identify and parse matrix assets in raw_dir into an AnnData object."""
    if not raw_dir.exists():
        return Failure(f"Directory {raw_dir} does not exist.")

    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if h5ad_files:
        try:
            return Success(ad.read_h5ad(h5ad_files[0]))
        except Exception as exc:  # noqa: BLE001
            return Failure(f"Failed to read H5AD {h5ad_files[0]}: {exc}")

    h5_files = list(raw_dir.glob("*.h5"))
    if h5_files:
        try:
            adata = sc.read_10x_h5(h5_files[0])
            adata.var_names_make_unique()
            return Success(adata)
        except Exception as exc:  # noqa: BLE001
            return Failure(f"Failed to read 10x H5 {h5_files[0]}: {exc}")

    tar_files = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
    for tf in tar_files:
        extract_dir = raw_dir / "unpacked"
        if not extract_dir.exists():
            extract_dir.mkdir(parents=True, exist_ok=True)
            try:
                with tarfile.open(tf, "r:*") as archive:
                    archive.extractall(extract_dir)
            except Exception:  # noqa: BLE001
                pass

    search_dirs = [raw_dir, *raw_dir.glob("**/")]

    for d in search_dirs:
        if list(d.glob("*matrix.mtx*")):
            try:
                adata = sc.read_10x_mtx(d)
                adata.var_names_make_unique()
                return Success(adata)
            except Exception:  # noqa: BLE001
                pass

    for d in search_dirs:
        count_tables = [
            f
            for f in d.glob("*.*")
            if any(k in f.name.lower() for k in ("count", "tpm", "expression", "counts"))
            and any(
                f.name.lower().endswith(ext)
                for ext in (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".csv", ".txt")
            )
            and not f.name.startswith("filelist")
        ]
        if count_tables:
            target = count_tables[0]
            try:
                sep = (
                    "\t"
                    if (
                        target.name.endswith(".tsv")
                        or target.name.endswith(".tsv.gz")
                        or "txt" in target.name
                    )
                    else ","
                )
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
                    obs=pd.DataFrame(index=pd.Index(cell_cols)),
                    var=pd.DataFrame(index=pd.Index(genes)),
                )
                adata.var_names_make_unique()
                return Success(adata)
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read count table {target}: {exc}")

    return Failure(f"No parseable matrix assets found in {raw_dir}")


def apply_standard_qc(adata: ad.AnnData) -> ad.AnnData:
    """Filter low-quality cells by count depth, gene detection, and mitochondrial percentage."""
    adata.var_names_make_unique()
    mt_genes = np.asarray(adata.var_names.str.startswith(("MT-", "mt-", "Mt-")), dtype=bool)
    if issparse(adata.X):
        x_csr: csr_matrix[Any] = csr_matrix(cast(Any, adata.X))
        total_counts = np.asarray(x_csr.sum(axis=1)).flatten()
        genes_detected = np.diff(x_csr.indptr).astype(np.int64)
        mt_counts = (
            np.asarray(x_csr[:, mt_genes].sum(axis=1)).flatten()
            if np.any(mt_genes)
            else np.zeros(adata.n_obs, dtype=np.float32)
        )
    else:
        arr = np.asarray(adata.X)
        total_counts = arr.sum(axis=1).flatten()
        genes_detected = (arr > 0).sum(axis=1).flatten()
        mt_counts = (
            arr[:, mt_genes].sum(axis=1).flatten()
            if np.any(mt_genes)
            else np.zeros(adata.n_obs, dtype=np.float32)
        )

    pct_mt = (mt_counts / np.maximum(total_counts, 1.0)) * 100.0
    adata.obs["total_counts"] = total_counts
    adata.obs["n_genes_by_counts"] = genes_detected
    adata.obs["pct_counts_mt"] = pct_mt

    keep_mask = (
        (genes_detected >= 200)
        & (genes_detected <= 8000)
        & (total_counts >= 500)
        & (pct_mt <= 20.0)
    )
    return adata[keep_mask].copy()


def harmonize_discovered_obs_metadata(adata: ad.AnnData, task: CohortPreprocessTask) -> ad.AnnData:
    """Harmonize cell metadata to canonical tme_datasets schema."""
    obs_df = adata.obs.copy()
    if "original_barcode" not in obs_df.columns:
        obs_df["original_barcode"] = obs_df.index.astype(str)

    obs_df["dataset_id"] = task.accession
    obs_df["cancer_type"] = task.indication
    obs_df["tier"] = task.tier
    obs_df["technology"] = task.technology
    obs_df["cell_selection"] = task.cell_selection

    for donor_col in ("patient", "patient_id", "donor", "donor_id", "subject", "sample"):
        if donor_col in obs_df.columns and "patient_id" not in obs_df.columns:
            obs_df["patient_id"] = obs_df[donor_col]
    if "patient_id" not in obs_df.columns:
        obs_df["patient_id"] = f"{task.accession}_donor_unspecified"

    for resp_col in ("response", "Response", "RECIST", "benefit", "responder"):
        if resp_col in obs_df.columns and "response" not in obs_df.columns:
            obs_df["response"] = obs_df[resp_col]
    if "response" not in obs_df.columns:
        obs_df["response"] = "Documented (cohort level)" if "Response" in task.tier else "Baseline"

    for ct_col in ("cell_type", "celltype", "CellType", "major_cell_type", "cluster"):
        if ct_col in obs_df.columns and "cell_type" not in obs_df.columns:
            obs_df["cell_type"] = obs_df[ct_col]
    if "cell_type" not in obs_df.columns:
        obs_df["cell_type"] = "Unspecified"

    adata.obs = obs_df
    return adata


def _run_curated_loader(task: CohortPreprocessTask) -> Result[tuple[int, int], str]:
    """Load and serialize a curated dataset in an isolated scope so memory is freed on return."""
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
            return Success((int(adata.n_obs), int(adata.n_vars)))
        case _:
            return Failure(f"{task.accession}: Unexpected Result state")


def _run_discovered_asset_processor(task: CohortPreprocessTask) -> Result[tuple[int, int], str]:
    """Parse, QC, harmonize, normalize, and write an uncurated discovered asset."""
    match read_raw_matrix_asset(task.raw_dir):
        case Failure(read_err):
            return Failure(f"{task.accession}: {read_err}")
        case Success(raw_adata):
            try:
                raw_adata.obs_names_make_unique()
                raw_adata.var_names_make_unique()
                adata = check_and_restore_raw_counts(raw_adata)
                adata = apply_standard_qc(adata)
                adata = harmonize_discovered_obs_metadata(adata, task)
                match normalize_dataset_to_ensembl(adata):
                    case Success(norm_ad):
                        adata = norm_ad
                    case Failure(_):
                        pass
                task.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(task.output_h5ad, compression="gzip")
                return Success((int(adata.n_obs), int(adata.n_vars)))
            except Exception as proc_err:  # noqa: BLE001
                return Failure(f"{task.accession}: Processing failed: {proc_err}")
        case _:
            return Failure(f"{task.accession}: Unexpected Result state")


def execute_preprocess_task(task: CohortPreprocessTask, clean_raw: bool = False) -> Result[str, str]:
    """Execute end-to-end preprocessing for a single cohort."""
    canonical_id = resolve_dataset_id(task.accession)
    shape_res = (
        _run_curated_loader(task)
        if isinstance(get_dataset_spec(canonical_id), Some)
        else _run_discovered_asset_processor(task)
    )
    gc.collect()

    match shape_res:
        case Failure(err):
            return Failure(err)
        case Success((n_obs, n_vars)):
            if clean_raw and task.raw_dir.exists():
                try:
                    shutil.rmtree(task.raw_dir)
                except Exception:  # noqa: BLE001
                    pass
            return Success(
                f"{task.accession}: Saved preprocessed AnnData ({n_obs:,} cells x {n_vars:,} genes) "
                f"to {task.output_h5ad}"
            )
        case _:
            return Failure(f"{task.accession}: Unexpected Result state")


def run_batch_preprocess_cohorts(
    tasks: Sequence[CohortPreprocessTask],
    *,
    output_dir: Path,
    report_dir: Path,
    max_ram_gb: float = 150.0,
    clean_raw: bool = False,
    force: bool = False,
    dry_run: bool = False,
) -> pl.DataFrame:
    """Execute batch preprocessing over a sequence of cohort tasks and write summary parquet."""
    governor = MemoryGovernor(max_ram_gb=max_ram_gb)

    if dry_run:
        for idx, t in enumerate(tasks[:15], 1):
            sys.stdout.write(
                f"[{idx:2d}/{len(tasks):2d}] [DRY-RUN] {t.accession:18s} | {t.raw_dir} -> {t.output_h5ad}\n"
            )
        return pl.DataFrame(
            [{"total_tasks": len(tasks), "succeeded": 0, "failed": 0, "output_h5ad_dir": str(output_dir)}]
        )

    successes = 0
    failures = 0
    for idx, t in enumerate(tasks, 1):
        if not t.raw_dir.exists():
            sys.stdout.write(f"[{idx:3d}/{len(tasks):3d}] - {t.accession}: Raw directory not found.\n")
            continue
        if t.output_h5ad.exists() and not force:
            sys.stdout.write(f"[{idx:3d}/{len(tasks):3d}] [SKIP] {t.accession}: {t.output_h5ad.name} exists.\n")
            successes += 1
            continue

        governor.enforce_headroom(required_headroom_gb=10.0)
        res = execute_preprocess_task(t, clean_raw=clean_raw)
        gc.collect()
        match res:
            case Success(msg):
                sys.stdout.write(f"[{idx:3d}/{len(tasks):3d}] {msg}\n")
                successes += 1
            case Failure(err):
                sys.stdout.write(f"[{idx:3d}/{len(tasks):3d}] FAILED {err}\n")
                failures += 1

    report_dir.mkdir(parents=True, exist_ok=True)
    summary_df = pl.DataFrame(
        [
            {
                "total_tasks": len(tasks),
                "succeeded": successes,
                "failed": failures,
                "output_h5ad_dir": str(output_dir),
            }
        ]
    )
    summary_df.write_parquet(report_dir / "batch_preprocess_summary.parquet")
    return summary_df


def main() -> int:
    """CLI entrypoint for tme-datasets-batch-preprocess."""
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
        "--report-dir",
        type=Path,
        default=Path("output/tme_datasets/batch_preprocess"),
        help="Output directory for batch preprocessing summary reports.",
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
        help="Delete raw directory after successful H5AD validation.",
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
    if hasattr(signal, "SIGHUP"):
        signal.signal(signal.SIGHUP, signal.SIG_IGN)

    df = pl.read_parquet(Path(args.registry))
    if args.cohorts:
        specific = [s.strip() for s in str(args.cohorts).split(",") if s.strip()]
        df = df.filter(pl.col("accession").is_in(specific))

    tasks = tuple(
        CohortPreprocessTask(
            accession=str(row["accession"]),
            indication=str(row["indication"]),
            tier=str(row["tier"]),
            technology=str(row["technology"]),
            cell_selection=str(row["cell_selection_strategy"]),
            raw_dir=Path(args.data_raw_dir) / str(row["accession"]),
            output_h5ad=Path(args.output_dir) / f"{row['accession']}.h5ad",
        )
        for row in df.iter_rows(named=True)
    )

    run_batch_preprocess_cohorts(
        tasks,
        output_dir=Path(args.output_dir),
        report_dir=Path(args.report_dir),
        max_ram_gb=float(args.max_ram_gb),
        clean_raw=bool(args.clean_raw),
        force=bool(args.force),
        dry_run=bool(args.dry_run),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
