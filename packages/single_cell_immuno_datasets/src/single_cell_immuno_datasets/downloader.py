"""Functional data downloading module for ICB single-cell datasets."""

import tarfile
import urllib.request
from pathlib import Path
from returns.result import Result, Success, Failure
from .config import DatasetSpec, DataDirectories, TIER_1_DATASETS


from tme_datasets.download import download_single_file


def download_dataset(spec: DatasetSpec, raw_dir: Path) -> Result[tuple[Path, ...], str]:
    """Downloads all raw files associated with a DatasetSpec into raw_dir / spec.accession and unpacks tar archives."""
    dataset_dir = raw_dir / spec.accession
    downloaded: list[Path] = []

    for name, url in spec.download_urls.items():
        if url.endswith(".tar"):
            ext = ".tar"
        elif url.endswith(".csv.gz"):
            ext = ".csv.gz"
        elif url.endswith(".txt.gz"):
            ext = ".txt.gz"
        elif "h5ad" in url:
            ext = ".h5ad"
        else:
            ext = ".gz"

        dest_file = dataset_dir / f"{spec.accession}_{name}{ext}"
        
        match download_single_file(url, dest_file):
            case Success(p):
                downloaded.append(p)
                if p.name.endswith(".tar"):
                    print(f"Unpacking TAR archive for {spec.accession}: {p.name}...")
                    try:
                        with tarfile.open(p, "r:*") as tar:
                            tar.extractall(path=dataset_dir)
                    except Exception as ex:
                        print(f"Warning unpacking {p.name}: {ex}")
            case Failure(err):
                return Failure(f"Dataset {spec.accession} failed on file '{name}': {err}")

    return Success(tuple(downloaded))


def download_all_tier1(dirs: DataDirectories) -> Result[dict[str, tuple[Path, ...]], str]:
    """Downloads all Tier 1 ICB single-cell datasets into remote/local raw storage directory."""
    dirs.raw_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, tuple[Path, ...]] = {}

    for spec in TIER_1_DATASETS:
        match download_dataset(spec, dirs.raw_dir):
            case Success(paths):
                results[spec.accession] = paths
                print(f"Successfully downloaded {spec.accession} ({len(paths)} files)")
            case Failure(err):
                print(f"Warning: {err}")

    return Success(results)
