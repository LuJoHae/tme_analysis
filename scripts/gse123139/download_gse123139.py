"""Download script for GSE123139 single-cell melanoma dataset from GEO."""

import argparse
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Final, assert_never

from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

GEO_BASE_URL: Final[str] = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139"
SUPPL_FILES: Final[tuple[str, ...]] = (
    "GSE123139_RAW.tar",
    "GSE123139_T_cells_tcrb_v2.txt.gz",
)
MATRIX_FILE: Final[str] = "GSE123139_series_matrix.txt.gz"


class DownloadConfig(BaseModel):
    """Immutable configuration for GSE123139 downloads."""

    model_config = ConfigDict(frozen=True)
    out_dir: Path
    extract_raw: bool = True


def download_single_file(url: str, dest_path: Path) -> Result[Path, str]:
    """Pure IO action to retrieve a URL to a target filesystem location."""
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        if not dest_path.exists():
            print(f"Downloading: {url} -> {dest_path}")
            _ = urllib.request.urlretrieve(url, dest_path)
        else:
            print(f"File already exists: {dest_path}")
        return Success(dest_path)
    except Exception as err:
        return Failure(f"Failed downloading {url}: {err}")


def extract_tar_archive(tar_path: Path, extract_dir: Path) -> Result[Path, str]:
    """Extract a tar archive to the specified directory safely."""
    try:
        extract_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting {tar_path} into {extract_dir}...")
        with tarfile.open(tar_path, "r:*") as archive:
            archive.extractall(path=extract_dir)
        return Success(extract_dir)
    except Exception as err:
        return Failure(f"Failed to extract tar archive {tar_path}: {err}")


def fetch_gse123139(config: DownloadConfig) -> Result[tuple[Path, ...], str]:
    """Monadic workflow for downloading GSE123139 dataset components."""
    downloaded_paths: list[Path] = []

    # 1. Download Supplementary Files
    for filename in SUPPL_FILES:
        url = f"{GEO_BASE_URL}/suppl/{filename}"
        dest = config.out_dir / filename
        match download_single_file(url, dest):
            case Success(path):
                downloaded_paths.append(path)
            case Failure(err):
                return Failure(err)
            case _ as unreachable:
                assert_never(unreachable)

    # 2. Download Series Matrix Metadata
    matrix_url = f"{GEO_BASE_URL}/matrix/{MATRIX_FILE}"
    matrix_dest = config.out_dir / MATRIX_FILE
    match download_single_file(matrix_url, matrix_dest):
        case Success(path):
            downloaded_paths.append(path)
        case Failure(err):
            return Failure(err)
        case _ as unreachable:
            assert_never(unreachable)

    # 3. Optionally Extract Raw Files
    raw_tar = config.out_dir / "GSE123139_RAW.tar"
    if config.extract_raw and raw_tar.exists():
        raw_extract_dir = config.out_dir / "raw_counts"
        match extract_tar_archive(raw_tar, raw_extract_dir):
            case Success(path):
                downloaded_paths.append(path)
            case Failure(err):
                return Failure(err)
            case _ as unreachable:
                assert_never(unreachable)

    return Success(tuple(downloaded_paths))


def main() -> None:
    """CLI entry point for downloading GSE123139 data."""
    parser = argparse.ArgumentParser(description="Download GSE123139 scRNA-seq melanoma dataset.")
    _ = parser.add_argument("--out-dir", required=True, help="Directory to save downloaded dataset.")
    _ = parser.add_argument(
        "--no-extract",
        action="store_false",
        dest="extract_raw",
        help="Skip extraction of raw count tar file.",
    )

    args = parser.parse_args()
    config = DownloadConfig(out_dir=Path(args.out_dir), extract_raw=args.extract_raw)

    match fetch_gse123139(config):
        case Success(paths):
            print("Successfully finished downloading GSE123139 dataset:")
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
