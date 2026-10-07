#!/usr/bin/env python3
"""Automated, Phased Downloader for Discovered Solid Tumor scRNA-seq Cohorts.

Downloads raw supplementary matrices for NCBI GEO accessions and streams curated
.h5ad AnnData files from CZ CELLxGENE into the repository data/raw directory structure.

Supports:
- Phased download queues: --tier tier1_response (Phase 1), tier1_all, tier2, all.
- Filtering by cancer indication (--indications Melanoma,NSCLC).
- Target accession lists (--cohorts GSE120575,GSE123139).
- Multi-threaded concurrent transfers.
- Dry-run planning mode (--dry-run).
- Integration with tme_datasets download engine.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import shutil
import sys
import time
import urllib.request
from pathlib import Path
from typing import Sequence

import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from tme_datasets.download.geo import download_geo_supplementary


class DownloadTask(BaseModel):
    """Specification of a download task for a single cohort."""

    model_config = ConfigDict(frozen=True)

    accession: str
    indication: str
    tier: str
    repository: str
    download_urls: tuple[str, ...]
    matrix_files: tuple[str, ...]
    destination_dir: Path


def resolve_tasks_from_registry(
    registry_path: Path,
    tier_filter: str,
    indications: Sequence[str],
    specific_cohorts: Sequence[str],
    data_raw_dir: Path,
    limit: int | None = None,
) -> tuple[DownloadTask, ...]:
    """Parse Parquet registry and resolve list of target DownloadTasks."""
    df = pl.read_parquet(registry_path)

    # Filter by specific cohorts if provided
    if specific_cohorts:
        df = df.filter(pl.col("accession").is_in(list(specific_cohorts)))
    else:
        # Filter by tier
        match tier_filter:
            case "tier0_benchmark":
                df = df.filter(pl.col("tier") == "Tier 0 (Benchmark Core)")
            case "tier1_response":
                df = df.filter(pl.col("tier").is_in(["Tier 0 (Benchmark Core)", "Tier 1 (ICB Response)"]))
            case "tier1_all":
                df = df.filter(pl.col("tier").str.starts_with("Tier 0") | pl.col("tier").str.starts_with("Tier 1"))
            case "tier1_only":
                df = df.filter(pl.col("tier").str.starts_with("Tier 1"))
            case "tier2":
                df = df.filter(pl.col("tier") == "Tier 2 (Baseline Atlas)")
            case "all":
                pass
            case _:
                df = df.filter(pl.col("tier") == tier_filter)

        # Filter by indication
        if indications:
            df = df.filter(pl.col("indication").is_in(list(indications)))

    if limit is not None and limit > 0:
        df = df.head(limit)

    tasks: list[DownloadTask] = []
    for row in df.iter_rows(named=True):
        acc = row["accession"]
        repo = row["repository"]
        ind = row["indication"]
        tier = row["tier"]

        urls_raw = row.get("download_urls", "")
        urls = tuple(u.strip() for u in urls_raw.split(";") if u.strip())

        mat_raw = row.get("matrix_files", "")
        mat_files = tuple(m.strip() for m in mat_raw.split(";") if m.strip())

        dest_dir = data_raw_dir / acc

        tasks.append(
            DownloadTask(
                accession=acc,
                indication=ind,
                tier=tier,
                repository=repo,
                download_urls=urls,
                matrix_files=mat_files,
                destination_dir=dest_dir,
            )
        )

    return tuple(tasks)


def stream_download_file(url: str, dest_file: Path, timeout: int = 120) -> Result[Path, str]:
    """Stream download a file from an HTTPS URL to a local destination path."""
    try:
        dest_file.parent.mkdir(parents=True, exist_ok=True)
        temp_dest = dest_file.with_suffix(dest_file.suffix + ".part")
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "TME-Analysis-Downloader/1.0 (academic research)"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(temp_dest, "wb") as f_out:
            shutil.copyfileobj(resp, f_out)
        temp_dest.replace(dest_file)
        return Success(dest_file)
    except Exception as exc:
        if temp_dest.exists():
            temp_dest.unlink(missing_ok=True)
        return Failure(f"Failed to download {url}: {exc}")


def execute_download_task(task: DownloadTask) -> Result[str, str]:
    """Download data assets for a single cohort with timing and size metrics."""
    t0 = time.time()
    task.destination_dir.mkdir(parents=True, exist_ok=True)

    if task.repository == "NCBI GEO":
        # Download GEO supplementary files with intra-cohort file parallelism
        expected = [f for f in task.matrix_files if not f.startswith("filelist") and not f.endswith(".txt")]
        match download_geo_supplementary(
            gse_id=task.accession,
            dest_dir=task.destination_dir,
            expected_files=expected if expected else None,
            max_workers=4,
        ):
            case Success(files):
                elapsed = max(time.time() - t0, 0.01)
                total_bytes = sum(f.stat().st_size for f in files if f.exists())
                mb = total_bytes / (1024 * 1024)
                speed = mb / elapsed
                return Success(f"{task.accession}: Downloaded {len(files)} files ({mb:.1f} MB in {elapsed:.1f}s, {speed:.2f} MB/s)")
            case Failure(err):
                # Fallback to direct download URLs if available
                if task.download_urls:
                    downloaded: list[Path] = []
                    for url in task.download_urls:
                        fname = url.split("/")[-1]
                        target_p = task.destination_dir / fname
                        if target_p.exists() and target_p.stat().st_size > 0:
                            downloaded.append(target_p)
                            continue
                        match stream_download_file(url, target_p):
                            case Success(p):
                                downloaded.append(p)
                            case Failure(_):
                                pass
                    if downloaded:
                        elapsed = max(time.time() - t0, 0.01)
                        total_bytes = sum(f.stat().st_size for f in downloaded if f.exists())
                        mb = total_bytes / (1024 * 1024)
                        speed = mb / elapsed
                        return Success(f"{task.accession}: Fallback downloaded {len(downloaded)} files ({mb:.1f} MB in {elapsed:.1f}s, {speed:.2f} MB/s)")
                return Failure(f"{task.accession}: GEO download failed: {err}")

    elif task.repository == "CZ CELLxGENE":
        # Download CELLxGENE .h5ad
        if not task.download_urls:
            return Failure(f"{task.accession}: No CELLxGENE download URLs available")

        h5ad_url = task.download_urls[0]
        fname = task.matrix_files[0] if task.matrix_files else f"{task.accession}.h5ad"
        target_file = task.destination_dir / fname

        if target_file.exists() and target_file.stat().st_size > 1024:
            mb = target_file.stat().st_size / (1024 * 1024)
            return Success(f"{task.accession}: Cached on disk ({mb:.1f} MB)")

        match stream_download_file(h5ad_url, target_file):
            case Success(p):
                elapsed = max(time.time() - t0, 0.01)
                mb = p.stat().st_size / (1024 * 1024)
                speed = mb / elapsed
                return Success(f"{task.accession}: Streamed CELLxGENE .h5ad ({mb:.1f} MB in {elapsed:.1f}s, {speed:.2f} MB/s)")
            case Failure(err):
                return Failure(f"{task.accession}: CELLxGENE stream failed: {err}")

    return Failure(f"{task.accession}: Unsupported repository {task.repository}")


def run_downloader(
    tasks: Sequence[DownloadTask],
    max_workers: int = 8,
    dry_run: bool = False,
) -> None:
    """Orchestrate downloading across multiple cohorts concurrently using a thread pool."""
    print(f"\n=======================================================")
    print(f"  TME Analysis Parallel scRNA-seq Downloader")
    print(f"  Total cohorts to download: {len(tasks)}")
    print(f"  Concurrent cohort workers: {max_workers}")
    print(f"  Intra-cohort file workers: 4")
    print(f"  Dry run mode:              {dry_run}")
    print(f"=======================================================\n")

    if dry_run:
        for idx, t in enumerate(tasks, 1):
            print(f"[{idx:3d}/{len(tasks):3d}] [DRY-RUN] {t.accession:18s} | {t.indication:10s} | {t.tier:22s} | {t.repository:12s} -> {t.destination_dir}")
        print("\nDry run complete. No network connections established.")
        return

    success_count = 0
    failure_count = 0
    total = len(tasks)
    start_all = time.time()

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(execute_download_task, t): t for t in tasks
        }

        print(f"Dispatched {len(tasks)} cohort download tasks across {max_workers} worker threads.\n")

        for idx, future in enumerate(concurrent.futures.as_completed(future_to_task), 1):
            task = future_to_task[future]
            now_str = time.strftime("%H:%M:%S")
            try:
                result = future.result()
                match result:
                    case Success(msg):
                        print(f"[{now_str}] [{idx:3d}/{total:3d}] ✓ {msg}")
                        success_count += 1
                    case Failure(err):
                        print(f"[{now_str}] [{idx:3d}/{total:3d}] ✗ {err}")
                        failure_count += 1
            except Exception as exc:
                print(f"[{now_str}] [{idx:3d}/{total:3d}] ✗ {task.accession}: Unhandled Exception: {exc}")
                failure_count += 1

    total_time = max(time.time() - start_all, 0.01)
    print(f"\n=======================================================")
    print(f"  Parallel Download Summary")
    print(f"  Completed: {success_count} succeeded, {failure_count} failed in {total_time:.1f}s")
    print(f"=======================================================\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel Downloader for scRNA-seq Cohorts.")
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
        help="Target raw data directory.",
    )
    parser.add_argument(
        "--tier",
        type=str,
        default="tier0_benchmark",
        choices=["tier0_benchmark", "tier1_response", "tier1_all", "tier1_only", "tier2", "all"],
        help="Cohort tier filter (default: tier0_benchmark).",
    )
    parser.add_argument(
        "--indications",
        type=str,
        default="",
        help="Comma-separated indications to download (default: all).",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Comma-separated specific accessions to download.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum cohorts to download.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Concurrent cohort download workers (default: 8).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned downloads without executing network transfers.",
    )

    args = parser.parse_args()

    indications = tuple(s.strip() for s in args.indications.split(",") if s.strip())
    specific_cohorts = tuple(s.strip() for s in args.cohorts.split(",") if s.strip())

    tasks = resolve_tasks_from_registry(
        registry_path=args.registry,
        tier_filter=args.tier,
        indications=indications,
        specific_cohorts=specific_cohorts,
        data_raw_dir=args.data_raw_dir,
        limit=args.limit,
    )

    if not tasks:
        print("No cohorts matched the specified filters.", file=sys.stderr)
        sys.exit(1)

    run_downloader(tasks=tasks, max_workers=args.max_workers, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
