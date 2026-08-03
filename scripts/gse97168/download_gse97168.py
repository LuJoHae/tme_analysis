"""Download script for GSE97168 dataset from GEO."""

import argparse
import sys
import urllib.request
from pathlib import Path
from typing import Final, assert_never

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

SUPPL_URLS: Final[dict[str, str]] = {
    "GSE97168_metadata.txt.gz": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97168&format=file&file=GSE97168%5Fmetadata%2Etxt%2Egz",
    "GSE97168_umitab.txt.gz": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97168&format=file&file=GSE97168%5Fumitab%2Etxt%2Egz",
}


class DownloadConfig(BaseModel):
    """Immutable configuration for GSE97168 dataset download."""

    model_config = ConfigDict(frozen=True)
    out_dir: Path


import subprocess


def download_single_file(url: str, dest_path: Path) -> Result[Path, str]:
    """Retrieve a single remote URL resource using system curl."""
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        if not dest_path.exists():
            print(f"Downloading {url} -> {dest_path}")
            res = subprocess.run(["curl", "-sSL", "-o", str(dest_path), url], capture_output=True, text=True)
            if res.returncode != 0:
                return Failure(f"curl failed for {url}: {res.stderr}")
        else:
            print(f"File already exists: {dest_path}")
        return Success(dest_path)
    except Exception as err:
        return Failure(f"Failed downloading {url}: {err}")


def fetch_gse97168(config: DownloadConfig) -> Result[tuple[Path, ...], str]:
    """Monadic workflow for downloading GSE97168 dataset components."""
    downloaded_paths: list[Path] = []

    # 1. Download Supplementary Files (umitab and metadata)
    for filename, url in SUPPL_URLS.items():
        dest = config.out_dir / filename
        match download_single_file(url, dest):
            case Success(path):
                downloaded_paths.append(path)
            case Failure(err):
                return Failure(err)
            case _ as unreachable:
                assert_never(unreachable)

    return Success(tuple(downloaded_paths))


def main() -> None:
    """CLI entry point for downloading GSE97168 data."""
    parser = argparse.ArgumentParser(description="Download GSE97168 dataset.")
    _ = parser.add_argument("--out-dir", required=True, help="Target directory for downloaded files.")

    args = parser.parse_args()
    config = DownloadConfig(out_dir=Path(args.out_dir))

    match fetch_gse97168(config):
        case Success(paths):
            print("Successfully finished downloading GSE97168 dataset:")
            for p in paths:
                print(f"  - {p}")
            sys.exit(0)
        case Failure(err):
            print(f"Download failed: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
