"""Script 1: Download GEO Single-Cell Datasets.

Inspired by the single_cell_datasets package, this script downloads raw single-cell
RNA-seq datasets from GEO (e.g. GSE120575, GSE123139, GSE97168) into specified output directories.
"""

from pathlib import Path
import urllib.request
import gzip
from typing import NamedTuple, Protocol
from pydantic import BaseModel, ConfigDict
from returns.result import Result, Success, Failure


class GeoDatasetSpec(BaseModel):
    """Metadata specification for a GEO dataset to download."""
    model_config = ConfigDict(frozen=True)
    accession: str
    urls: dict[str, str]


class DownloadConfig(BaseModel):
    """Immutable configuration for downloading GEO datasets."""
    model_config = ConfigDict(frozen=True)
    out_dir: Path
    datasets: tuple[GeoDatasetSpec, ...] = (
        GeoDatasetSpec(
            accession="GSE120575",
            urls={
                "tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
                "meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/matrix/GSE120575_series_matrix.txt.gz",
            },
        ),
        GeoDatasetSpec(
            accession="GSE123139",
            urls={
                "tcr": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/suppl/GSE123139_T_cells_tcrb_v2.txt.gz",
                "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/matrix/GSE123139_series_matrix.txt.gz",
            },
        ),
        GeoDatasetSpec(
            accession="GSE97168",
            urls={
                "umitab": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97168/suppl/GSE97168_umitab.txt.gz",
            },
        ),
    )


def download_single_file(url: str, dest_path: Path) -> Result[Path, str]:
    """Downloads a single URL to a destination file if not already present."""
    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        if dest_path.exists() and dest_path.stat().st_size > 0:
            return Success(dest_path)

        print(f"Downloading {url} -> {dest_path}...")
        urllib.request.urlretrieve(url, dest_path)
        return Success(dest_path)
    except Exception as e:
        return Failure(f"Failed to download {url}: {str(e)}")


def download_geo_dataset(spec: GeoDatasetSpec, base_dir: Path) -> Result[tuple[Path, ...], str]:
    """Downloads all supplementary files associated with a GEO dataset specification."""
    dataset_dir = base_dir / spec.accession
    downloaded_paths = []
    
    for name, url in spec.urls.items():
        ext = ".txt.gz" if url.endswith(".txt.gz") else ".gz"
        dest_file = dataset_dir / f"{spec.accession}_{name}{ext}"
        match download_single_file(url, dest_file):
            case Success(p):
                downloaded_paths.append(p)
            case Failure(err):
                return Failure(f"Dataset {spec.accession} failed on file {name}: {err}")

    return Success(tuple(downloaded_paths))


def run_geo_download(config: DownloadConfig) -> Result[dict[str, tuple[Path, ...]], str]:
    """Downloads all configured GEO datasets."""
    results = {}
    for spec in config.datasets:
        match download_geo_dataset(spec, config.out_dir):
            case Success(paths):
                results[spec.accession] = paths
            case Failure(err):
                return Failure(err)
    return Success(results)


def main() -> None:
    """CLI entry point for downloading GEO datasets."""
    import argparse
    parser = argparse.ArgumentParser(description="Download GEO Single-Cell Datasets.")
    parser.add_argument("--out-dir", type=str, default="data/raw_geo", help="Output directory for raw downloaded files")
    args = parser.parse_args()

    config = DownloadConfig(out_dir=Path(args.out_dir).resolve())
    print(f"Starting GEO dataset downloads into: {config.out_dir}")

    match run_geo_download(config):
        case Success(res):
            print("\nSuccessfully downloaded GEO datasets:")
            for acc, paths in res.items():
                print(f"[{acc}] -> {len(paths)} files")
                for p in paths:
                    print(f"  - {p}")
        case Failure(err):
            print(f"Download failed: {err}")


if __name__ == "__main__":
    main()
