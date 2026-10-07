"""Automated URL downloading with retry and archive extraction."""

from __future__ import annotations

from pathlib import Path
import tarfile
import time
import urllib.request
from returns.result import Failure, Result, Success

from ..logging import get_logger

logger = get_logger("download")


def download_single_file(url: str, dest_path: Path) -> Result[Path, str]:
    """Download a file from a URL to dest_path if not already present or empty.

    Streams chunks to a temporary .part file and atomically renames upon completion.
    """
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        part_path = dest_path.with_name(dest_path.name + ".part")

        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0 (Python tme_datasets)"},
        )

        with urllib.request.urlopen(req) as response:
            content_len_header = response.info().get("Content-Length")
            total_bytes = int(content_len_header) if content_len_header else None

            # If dest_path exists, verify size matches Content-Length if available
            if dest_path.exists() and dest_path.stat().st_size > 0:
                if total_bytes is None or dest_path.stat().st_size == total_bytes:
                    logger.debug(
                        "File already present and verified: %s (%.1f MB)",
                        dest_path.name,
                        dest_path.stat().st_size / (1024 * 1024),
                    )
                    return Success(dest_path)
                else:
                    logger.warning(
                        "Existing file %s is incomplete (%d / %d bytes). Re-downloading...",
                        dest_path.name,
                        dest_path.stat().st_size,
                        total_bytes,
                    )
                    dest_path.unlink(missing_ok=True)

            logger.info("Downloading %s -> %s", url, dest_path)
            start_time = time.time()
            chunk_size = 1024 * 1024  # 1 MB
            downloaded = 0
            last_reported_mb = 0.0

            with open(part_path, "wb") as out_file:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    downloaded += len(chunk)

                    downloaded_mb = downloaded / (1024 * 1024)
                    # Report every 20 MB or at milestones
                    if downloaded_mb - last_reported_mb >= 20:
                        if total_bytes:
                            pct = (downloaded / total_bytes) * 100
                            total_mb = total_bytes / (1024 * 1024)
                            logger.info(
                                "Downloading %s: %.1f / %.1f MB (%.1f%%)",
                                dest_path.name,
                                downloaded_mb,
                                total_mb,
                                pct,
                            )
                        else:
                            logger.info(
                                "Downloading %s: %.1f MB downloaded",
                                dest_path.name,
                                downloaded_mb,
                            )
                        last_reported_mb = downloaded_mb

            # Verify and atomically move into place
            if total_bytes is not None and downloaded != total_bytes:
                part_path.unlink(missing_ok=True)
                return Failure(
                    f"Download incomplete for {dest_path.name}: {downloaded}/{total_bytes} bytes received"
                )

            part_path.replace(dest_path)

        elapsed = max(0.01, time.time() - start_time)
        final_mb = downloaded / (1024 * 1024)
        rate = final_mb / elapsed
        logger.info(
            "Successfully downloaded %s (%.1f MB in %.1fs, %.1f MB/s)",
            dest_path.name,
            final_mb,
            elapsed,
            rate,
        )
        return Success(dest_path)
    except Exception as exc:
        part_path = dest_path.with_name(dest_path.name + ".part")
        if part_path.exists():
            part_path.unlink(missing_ok=True)
        msg = f"Failed to download {url} -> {dest_path}: {exc}"
        logger.error(msg)
        return Failure(msg)



def unpack_tar(tar_path: Path, extract_dir: Path) -> Result[Path, str]:
    """Extract a tar archive to the specified directory."""
    if not tar_path.is_file():
        msg = f"Tar file not found: {tar_path}"
        logger.error(msg)
        return Failure(msg)

    try:
        extract_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Extracting archive %s -> %s...", tar_path.name, extract_dir)
        start_time = time.time()

        with tarfile.open(tar_path, "r:*") as tar:
            if hasattr(tarfile, "data_filter"):
                tar.extractall(path=extract_dir, filter="data")
            else:
                tar.extractall(path=extract_dir)

        elapsed = time.time() - start_time
        logger.info(
            "Successfully extracted %s in %.1fs",
            tar_path.name,
            elapsed,
        )
        return Success(extract_dir)
    except Exception as exc:
        msg = f"Failed to unpack tar {tar_path}: {exc}"
        logger.error(msg)
        return Failure(msg)
