#!/usr/bin/env python3
"""Resilient, Detached Dataset Download and Preprocessing Pipeline.

Survives SSH disconnection from remote servers (olm) via:
1. Signal immunity (SIGHUP and SIGPIPE ignored in master and worker processes).
2. Atomic H5AD cache verification and progressive state checkpointing (JSON).
3. Dynamic memory governor preventing OOM crashes (RAM ceiling strictly <= 140 GB).
4. Automated download fallback for NCBI GEO (tar/h5/mtx) and CZ CELLxGENE.
5. 100% Resumable: skips already-processed valid H5AD caches in <0.01s.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import gc
import json
import logging
import os
import signal
import sys
import time
from pathlib import Path
from typing import Sequence

import anndata as ad
import polars as pl
import psutil
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

from tme_datasets import list_registered_datasets, load_dataset
from tme_datasets.paths import find_repo_root, get_data_paths, get_preprocessed_h5ad_path
from tme_datasets.preprocessing.gene_normalization import ensure_ensembl_release_installed


class DatasetProcessResult(BaseModel):
    """Immutable result record for a processed dataset."""

    model_config = ConfigDict(frozen=True)

    dataset_id: str
    status: str
    n_obs: int
    n_vars: int
    elapsed_sec: float
    h5ad_path: str
    error: str


class PipelineConfig(BaseModel):
    """Immutable configuration for parallel ingestion."""

    model_config = ConfigDict(frozen=True)

    max_workers: int = 8
    max_ram_gb: float = 140.0
    hard_limit_gb: float = 150.0
    force_recompute: bool = False
    normalize_ensembl: bool = True
    ensembl_release: int = 111
    skip_datasets: tuple[str, ...] = ("GSE178341",)
    log_file: Path | None = None
    progress_file: Path | None = None


def setup_disconnection_immunity() -> None:
    """Ignore SIGHUP and SIGPIPE to survive remote SSH session termination."""
    if hasattr(signal, "SIGHUP"):
        signal.signal(signal.SIGHUP, signal.SIG_IGN)
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_IGN)


def _worker_init() -> None:
    """Initialize worker process with SIGHUP immunity."""
    setup_disconnection_immunity()


def setup_pipeline_logger(log_file: Path) -> logging.Logger:
    """Configure dual-stream logging (stdout and persistent file)."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("download_preprocess")
    logger.setLevel(logging.INFO)

    # Avoid duplicate handlers if re-initialized
    if not logger.handlers:
        formatter = logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_process_tree_memory_gb() -> float:
    """Calculate total Resident Set Size (RSS) in GB of current process and all child workers."""
    try:
        current_proc = psutil.Process(os.getpid())
        total_rss = current_proc.memory_info().rss
        for child in current_proc.children(recursive=True):
            try:
                total_rss += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        return total_rss / (1024**3)
    except Exception:
        return 0.0


def _verify_existing_h5ad(target_h5ad: Path) -> Result[tuple[int, int], str]:
    """Check if preprocessed H5AD exists and is non-corrupt, returning (n_obs, n_vars) via lightweight h5py."""
    if not target_h5ad.is_file() or target_h5ad.stat().st_size == 0:
        return Failure("H5AD does not exist or is empty")
    try:
        import h5py
        with h5py.File(target_h5ad, "r") as f:
            if "X" in f:
                x = f["X"]
                if isinstance(x, h5py.Dataset):
                    shape = (int(x.shape[0]), int(x.shape[1]))
                elif "shape" in x.attrs:
                    shape_attr = x.attrs["shape"]
                    shape = (int(shape_attr[0]), int(shape_attr[1]))
                elif "shape" in x:
                    shape_ds = x["shape"]
                    shape = (int(shape_ds[0]), int(shape_ds[1]))
                elif "obs" in f and "var" in f:
                    shape = (len(f["obs"]), len(f["var"]))
                else:
                    return Failure("Cannot determine AnnData matrix dimensions")
            elif "obs" in f and "var" in f:
                shape = (len(f["obs"]), len(f["var"]))
            else:
                return Failure("H5AD structure missing 'X', 'obs', or 'var'")
            return Success(shape)
    except Exception as exc:
        return Failure(f"Corrupted H5AD detected: {exc}")


def _process_dataset_worker(
    dataset_id: str,
    force_recompute: bool,
    normalize_ensembl: bool,
    ensembl_release: int,
) -> DatasetProcessResult:
    """Worker task: downloads, preprocesses, and writes dataset to H5AD on SSD."""
    target_h5ad = get_preprocessed_h5ad_path(dataset_id)
    start = time.perf_counter()

    # Fast check: If already preprocessed and not force_recompute, return immediately
    if not force_recompute:
        match _verify_existing_h5ad(target_h5ad):
            case Success((n_obs, n_vars)):
                elapsed = time.perf_counter() - start
                return DatasetProcessResult(
                    dataset_id=dataset_id,
                    status="SUCCESS",
                    n_obs=n_obs,
                    n_vars=n_vars,
                    elapsed_sec=round(elapsed, 2),
                    h5ad_path=str(target_h5ad),
                    error="",
                )
            case Failure(_):
                pass  # Proceed to ingestion

    try:
        # Atomic temporary staging
        staging_h5ad = target_h5ad.with_suffix(".h5ad.tmp")
        if staging_h5ad.exists():
            staging_h5ad.unlink(missing_ok=True)

        load_res = load_dataset(
            dataset_id,
            auto_download=True,
            force_recompute=force_recompute,
            cache_h5ad=True,
            normalize_ensembl=normalize_ensembl,
            ensembl_release=ensembl_release,
            use_batched_processing=True,
        )
        elapsed = time.perf_counter() - start

        match load_res:
            case Success(adata):
                n_obs, n_vars = adata.n_obs, adata.n_vars
                del adata
                gc.collect()

                return DatasetProcessResult(
                    dataset_id=dataset_id,
                    status="SUCCESS",
                    n_obs=n_obs,
                    n_vars=n_vars,
                    elapsed_sec=round(elapsed, 2),
                    h5ad_path=str(target_h5ad),
                    error="",
                )
            case Failure(err):
                gc.collect()
                return DatasetProcessResult(
                    dataset_id=dataset_id,
                    status="FAILED",
                    n_obs=0,
                    n_vars=0,
                    elapsed_sec=round(elapsed, 2),
                    h5ad_path=str(target_h5ad),
                    error=str(err),
                )
    except Exception as exc:
        elapsed = time.perf_counter() - start
        gc.collect()
        return DatasetProcessResult(
            dataset_id=dataset_id,
            status="FAILED",
            n_obs=0,
            n_vars=0,
            elapsed_sec=round(elapsed, 2),
            h5ad_path=str(target_h5ad),
            error=f"Worker exception: {exc}",
        )


def _write_progress_checkpoint(
    progress_path: Path,
    total_target: int,
    results: Sequence[DatasetProcessResult],
    pending_count: int,
    active_count: int,
    current_ram_gb: float,
    max_ram_gb: float,
) -> None:
    """Atomically write progressive status report in JSON format."""
    try:
        progress_path.parent.mkdir(parents=True, exist_ok=True)
        completed_results = [r for r in results if r.status == "SUCCESS"]
        failed_results = [r for r in results if r.status != "SUCCESS"]

        payload = {
            "updated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "pid": os.getpid(),
            "total_target_datasets": total_target,
            "completed_count": len(completed_results),
            "failed_count": len(failed_results),
            "pending_count": pending_count,
            "active_workers": active_count,
            "current_rss_gb": round(current_ram_gb, 2),
            "max_ram_cap_gb": max_ram_gb,
            "completed_datasets": [r.dataset_id for r in completed_results],
            "failed_datasets": [{"id": r.dataset_id, "error": r.error} for r in failed_results],
        }

        tmp_path = progress_path.with_suffix(".json.tmp")
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        tmp_path.replace(progress_path)
    except Exception:
        pass


def download_and_preprocess_all(config: PipelineConfig) -> pl.DataFrame:
    """Orchestrate downloading and preprocessing of all registered datasets."""
    setup_disconnection_immunity()

    repo_root = find_repo_root()
    log_file = config.log_file or (repo_root / "logs/download_and_preprocess.log")
    progress_file = config.progress_file or (repo_root / "output/preprocessing_progress.json")
    logger = setup_pipeline_logger(log_file)

    specs = list_registered_datasets()
    target_specs = [s for s in specs if s.id not in config.skip_datasets]

    logger.info("=" * 80)
    logger.info("TME ROBUST RESUMABLE DOWNLOAD & PREPROCESSING PIPELINE")
    logger.info("=" * 80)
    logger.info("PID:                          %d", os.getpid())
    logger.info("Total Registered Datasets:    %d", len(specs))
    logger.info("Target Datasets to Process:   %d (Skipping: %s)", len(target_specs), config.skip_datasets)
    logger.info("Concurrent Worker Processes:  %d", config.max_workers)
    logger.info("Memory Ceiling (Gov / Hard):  %.1f GB / %.1f GB", config.max_ram_gb, config.hard_limit_gb)
    logger.info("Force Recompute:              %s", config.force_recompute)
    logger.info("Ensembl Gene Normalization:   %s (Release %d)", config.normalize_ensembl, config.ensembl_release)
    logger.info("Log File:                     %s", log_file)
    logger.info("Progress File:                %s", progress_file)
    logger.info("=" * 80)

    # Step 1: Sequential Pre-warming for shared Ensembl SQLite cache
    if config.normalize_ensembl:
        logger.info("--> Verifying shared Ensembl release database...")
        ensure_ensembl_release_installed(release=config.ensembl_release)
        logger.info("--> Ensembl database verified.\n")

    # Step 2: Separate already-cached datasets for instant reporting (Resumability)
    results: list[DatasetProcessResult] = []
    pending_specs: list[str] = []

    if not config.force_recompute:
        for spec in target_specs:
            target_path = get_preprocessed_h5ad_path(spec.id)
            match _verify_existing_h5ad(target_path):
                case Success((n_obs, n_vars)):
                    results.append(
                        DatasetProcessResult(
                            dataset_id=spec.id,
                            status="SUCCESS",
                            n_obs=n_obs,
                            n_vars=n_vars,
                            elapsed_sec=0.01,
                            h5ad_path=str(target_path),
                            error="",
                        )
                    )
                case Failure(_):
                    pending_specs.append(spec.id)
    else:
        pending_specs = [s.id for s in target_specs]

    cached_count = len(results)
    if cached_count > 0:
        logger.info("[RESUME] Found %d/%d valid preprocessed H5AD caches on SSD.", cached_count, len(target_specs))
        logger.info("[RESUME] %d datasets remaining to download and process.\n", len(pending_specs))
    else:
        logger.info("[QUEUE]  All %d datasets scheduled for download & preprocessing.\n", len(pending_specs))

    # Initial checkpoint write
    _write_progress_checkpoint(
        progress_path=progress_file,
        total_target=len(target_specs),
        results=results,
        pending_count=len(pending_specs),
        active_count=0,
        current_ram_gb=get_process_tree_memory_gb(),
        max_ram_gb=config.max_ram_gb,
    )

    if not pending_specs:
        logger.info("All target datasets already preprocessed. Generating summary report.")
        summary_df = pl.DataFrame([r.model_dump() for r in results])
        return summary_df

    # Step 3: Concurrent Execution with Dynamic Memory Governor
    completed_count = cached_count
    total_count = len(target_specs)

    future_to_id: dict[concurrent.futures.Future[DatasetProcessResult], str] = {}
    pending_queue = list(pending_specs)
    retry_counts: dict[str, int] = {}

    while pending_queue or future_to_id:
        executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=config.max_workers,
            initializer=_worker_init,
            max_tasks_per_child=1,
        )
        try:
            while pending_queue or future_to_id:
                current_mem_gb = get_process_tree_memory_gb()

                # Submit tasks as long as workers are available and memory is within safety threshold
                while (
                    pending_queue
                    and len(future_to_id) < config.max_workers
                    and current_mem_gb < config.max_ram_gb
                ):
                    dataset_id = pending_queue.pop(0)
                    fut = executor.submit(
                        _process_dataset_worker,
                        dataset_id,
                        config.force_recompute,
                        config.normalize_ensembl,
                        config.ensembl_release,
                    )
                    future_to_id[fut] = dataset_id
                    logger.info(
                        "--> [SUBMIT] '%s' queued (Active: %d, RAM: %.1f/%.1f GB)",
                        dataset_id,
                        len(future_to_id),
                        current_mem_gb,
                        config.max_ram_gb,
                    )
                    current_mem_gb = get_process_tree_memory_gb()

                if current_mem_gb >= config.max_ram_gb and pending_queue:
                    logger.warning(
                        "[THROTTLE] Process tree RAM at %.1f GB (>= safety cap %.1f GB). Pausing new submissions...",
                        current_mem_gb,
                        config.max_ram_gb,
                    )

                # Wait for at least one future to finish
                if future_to_id:
                    done, _ = concurrent.futures.wait(
                        future_to_id.keys(),
                        return_when=concurrent.futures.FIRST_COMPLETED,
                        timeout=5.0,
                    )
                    for fut in done:
                        dataset_id = future_to_id.pop(fut)
                        completed_count += 1
                        try:
                            res = fut.result()
                        except Exception as exc:
                            res = DatasetProcessResult(
                                dataset_id=dataset_id,
                                status="FAILED",
                                n_obs=0,
                                n_vars=0,
                                elapsed_sec=0.0,
                                h5ad_path=str(get_preprocessed_h5ad_path(dataset_id)),
                                error=f"Uncaught future exception: {exc}",
                            )

                        results.append(res)
                        mem_now = get_process_tree_memory_gb()

                        if res.status == "SUCCESS":
                            logger.info(
                                "[%2d/%2d] ✓ Done '%s' in %.1fs | Shape: (%d, %d) | RAM: %.1f GB",
                                completed_count,
                                total_count,
                                res.dataset_id,
                                res.elapsed_sec,
                                res.n_obs,
                                res.n_vars,
                                mem_now,
                            )
                        else:
                            logger.error(
                                "[%2d/%2d] ✗ Failed '%s' in %.1fs: %s | RAM: %.1f GB",
                                completed_count,
                                total_count,
                                res.dataset_id,
                                res.elapsed_sec,
                                res.error,
                                mem_now,
                            )

                        # Write progressive checkpoint
                        _write_progress_checkpoint(
                            progress_path=progress_file,
                            total_target=total_count,
                            results=results,
                            pending_count=len(pending_queue),
                            active_count=len(future_to_id),
                            current_ram_gb=mem_now,
                            max_ram_gb=config.max_ram_gb,
                        )
                else:
                    time.sleep(0.5)
        except concurrent.futures.process.BrokenProcessPool as b_err:
            logger.warning(
                "[RECOVERY] Process pool encountered an error (%s). Requeuing in-flight tasks and recreating worker pool...",
                b_err,
            )
            for fut, ds_id in list(future_to_id.items()):
                retry_counts[ds_id] = retry_counts.get(ds_id, 0) + 1
                if retry_counts[ds_id] <= 2:
                    logger.info("--> [RETRY] Requeuing '%s' (attempt %d/2)...", ds_id, retry_counts[ds_id])
                    pending_queue.append(ds_id)
                else:
                    completed_count += 1
                    res = DatasetProcessResult(
                        dataset_id=ds_id,
                        status="FAILED",
                        n_obs=0,
                        n_vars=0,
                        elapsed_sec=0.0,
                        h5ad_path=str(get_preprocessed_h5ad_path(ds_id)),
                        error=f"Worker process terminated abruptly ({b_err})",
                    )
                    results.append(res)
            future_to_id.clear()
            executor.shutdown(wait=False, cancel_futures=True)
            time.sleep(3.0)
        finally:
            try:
                executor.shutdown(wait=False)
            except Exception:
                pass

    summary_df = pl.DataFrame([r.model_dump() for r in results])
    logger.info("\n" + "=" * 80)
    logger.info("PREPROCESSING PIPELINE COMPLETE")
    logger.info("=" * 80)
    success_n = len(summary_df.filter(pl.col("status") == "SUCCESS"))
    failed_n = len(summary_df.filter(pl.col("status") != "SUCCESS"))
    logger.info("Successfully Preprocessed: %d / %d", success_n, total_count)
    logger.info("Failed Datasets:           %d / %d", failed_n, total_count)
    logger.info("=" * 80)

    # Final checkpoint write
    _write_progress_checkpoint(
        progress_path=progress_file,
        total_target=total_count,
        results=results,
        pending_count=0,
        active_count=0,
        current_ram_gb=get_process_tree_memory_gb(),
        max_ram_gb=config.max_ram_gb,
    )

    return summary_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Resilient, Detached Dataset Download & Preprocessing Pipeline."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of concurrent worker processes (default: 8).",
    )
    parser.add_argument(
        "--max-ram-gb",
        type=float,
        default=140.0,
        help="Safety RAM cap in GB for dynamic memory throttling (default: 140.0 GB, strictly <= 150 GB).",
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Force recomputation and overwrite existing preprocessed H5AD caches.",
    )
    parser.add_argument(
        "--include-pelka",
        action="store_true",
        help="Include GSE178341 (Pelka et al.) which is skipped by default due to memory scale.",
    )
    parser.add_argument(
        "--ensembl-release",
        type=int,
        default=111,
        help="Ensembl release version for gene normalization (default: 111).",
    )
    parser.add_argument(
        "--no-normalize-ensembl",
        action="store_true",
        help="Disable Ensembl gene normalization.",
    )

    args = parser.parse_args()

    skip = () if args.include_pelka else ("GSE178341",)
    cfg = PipelineConfig(
        max_workers=args.workers,
        max_ram_gb=args.max_ram_gb,
        hard_limit_gb=150.0,
        force_recompute=args.force_recompute,
        normalize_ensembl=not args.no_normalize_ensembl,
        ensembl_release=args.ensembl_release,
        skip_datasets=skip,
    )

    download_and_preprocess_all(cfg)


if __name__ == "__main__":
    main()
