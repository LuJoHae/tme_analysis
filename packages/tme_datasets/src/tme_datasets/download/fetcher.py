"""Automated URL downloading with retry and archive extraction."""

from __future__ import annotations

import tarfile
import urllib.request
from pathlib import Path
from returns.result import Failure, Result, Success


def download_single_file(url: str, dest_path: Path) -> Result[Path, str]:
    """Download a file from a URL to dest_path if not already present or empty."""
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        if dest_path.exists() and dest_path.stat().st_size > 0:
            return Success(dest_path)

        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0 (Python tme_datasets)"},
        )
        with urllib.request.urlopen(req) as response, open(dest_path, "wb") as out_file:
            out_file.write(response.read())

        return Success(dest_path)
    except Exception as exc:
        return Failure(f"Failed to download {url} -> {dest_path}: {exc}")


def unpack_tar(tar_path: Path, extract_dir: Path) -> Result[Path, str]:
    """Extract a tar archive to the specified directory."""
    if not tar_path.is_file():
        return Failure(f"Tar file not found: {tar_path}")

    try:
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path, "r:*") as tar:
            tar.extractall(path=extract_dir)
        return Success(extract_dir)
    except Exception as exc:
        return Failure(f"Failed to unpack tar {tar_path}: {exc}")
