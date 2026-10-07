"""Google Drive download utilities using gdown."""

from __future__ import annotations

import logging
from pathlib import Path

import gdown
from returns.result import Failure, Result, Success

logger = logging.getLogger("tme_datasets")


def download_gdrive_file(
    file_id_or_url: str,
    dest_path: Path,
    quiet: bool = False,
) -> Result[Path, str]:
    """Downloads a file from Google Drive given a file ID or sharing URL.

    Args:
        file_id_or_url: Google Drive file ID or full sharing URL.
        dest_path: Target destination path on local disk.
        quiet: If True, suppresses download progress bar.

    Returns:
        Result[Path, str]: Success(dest_path) if successful, Failure(err) otherwise.
    """
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        url = (
            file_id_or_url
            if file_id_or_url.startswith("http")
            else f"https://drive.google.com/uc?id={file_id_or_url}"
        )
        logger.info("Downloading file from Google Drive (%s) -> %s", file_id_or_url, dest_path)
        output = gdown.download(url=url, output=str(dest_path), quiet=quiet)
        if output is None or not dest_path.exists() or dest_path.stat().st_size == 0:
            return Failure(f"Failed to download Google Drive file {file_id_or_url} to {dest_path}")
        logger.info(
            "Successfully downloaded %s (%.1f MB)",
            dest_path.name,
            dest_path.stat().st_size / (1024 * 1024),
        )
        return Success(dest_path)
    except Exception as exc:
        return Failure(f"Exception during Google Drive download of {file_id_or_url}: {exc}")
