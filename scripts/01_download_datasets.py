"""Script 1: Unified Downloader for All 18 Single-Cell Datasets.

Downloads/fetches raw single-cell dataset files for:
1. GEO HTTP supplementary downloads: SadeFeldmanDefiningTCell2018Adata, YostClonalReplacementTumor2019Adata, ZhengLandscapeInfiltratingTCells2017Adata.
2. singlecellrnasignature datasets: Azizi2018, Becker2022, Biermann2022, Borcherding2021, Cheng2021, Durante2020, JerbyArnon2018, Khaliq2022, Kim2020, Leader2021, Lu2022, Pelka2021, Pu2021, Qian2020, Sharma2020.
"""

from pathlib import Path
import urllib.request
from pydantic import BaseModel, ConfigDict
from returns.result import Result, Success, Failure


SINGLECELL_SIGNATURE_DATASETS = [
    "AziziSingleCellMapDiverse2018Adata",
    "BeckerSinglecellAnalysesDefine2022Adata",
    "BiermannDissectingTreatmentnaiveEcosystem2022Adata",
    "BorcherdingMappingImmuneEnvironment2021Adata",
    "ChengPancancerSinglecellTranscriptional2021Adata",
    "DuranteSinglecellAnalysisReveals2020Adata",
    "JerbyArnonCancerCellProgram2018Adata",
    "KhaliqRefiningColorectalCancer2022Adata",
    "KimSinglecellRNASequencing2020Adata",
    "LeaderSinglecellAnalysisHuman2021Adata",
    "LuSinglecellAtlasMulticellular2022Adata",
    "PelkaSpatiallyOrganizedMulticellular2021Adata",
    "PuSinglecellTranscriptomicAnalysis2021Adata",
    "QianPancancerBlueprintHeterogeneous2020aAdata",
    "SharmaOncofetalReprogrammingEndothelial2020Adata",
]


class GeoDatasetSpec(BaseModel):
    """Metadata specification for a single-cell dataset to download."""
    model_config = ConfigDict(frozen=True)
    accession: str
    urls: dict[str, str]


class DownloadConfig(BaseModel):
    """Immutable configuration for downloading single-cell datasets."""
    model_config = ConfigDict(frozen=True)
    out_dir: Path
    datasets: tuple[GeoDatasetSpec, ...] = (
        GeoDatasetSpec(
            accession="SadeFeldmanDefiningTCell2018Adata",
            urls={
                "tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
            },
        ),
        GeoDatasetSpec(
            accession="YostClonalReplacementTumor2019Adata",
            urls={
                "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/matrix/GSE123139_series_matrix.txt.gz",
            },
        ),
        GeoDatasetSpec(
            accession="ZhengLandscapeInfiltratingTCells2017Adata",
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


def download_dataset(spec: GeoDatasetSpec, base_dir: Path) -> Result[tuple[Path, ...], str]:
    """Downloads files associated with a dataset specification."""
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


def fetch_signature_dataset(name: str, base_dir: Path) -> Result[Path, str]:
    """Fetches and downloads a singlecellrnasignature dataset into raw_dir with automatic corruption cleanup."""
    try:
        import datalair
        import singlecellrnasignature.adata

        dataset_dir = base_dir / name
        dest_h5 = dataset_dir / f"{name}.h5ad"
        if dest_h5.exists() and dest_h5.stat().st_size > 0:
            return Success(dest_h5)

        print(f"Downloading singlecellrnasignature dataset '{name}'...")
        if hasattr(singlecellrnasignature.adata, name):
            cls = getattr(singlecellrnasignature.adata, name)
            ds_instance = cls()
            lair = datalair.Lair()
            try:
                lair.safe_derive(ds_instance)
            except Exception as e:
                # Remove potentially corrupted download cache in ~/.cache/datalair and /tmp
                print(f"Retrying {name} after download error: {e}")
                cache_dir = Path.home() / ".cache" / "datalair"
                if cache_dir.exists():
                    import shutil
                    for item in cache_dir.glob(f"*{name}*"):
                        if item.is_dir():
                            shutil.rmtree(item, ignore_errors=True)
                        else:
                            item.unlink(missing_ok=True)
                
                # Clear /tmp temporary files created during download
                tmp_dir = Path("/tmp")
                if tmp_dir.exists():
                    import shutil
                    for item in tmp_dir.glob("tmp*"):
                        try:
                            if item.is_dir() and any(name in f.name or "GSE184362" in f.name for f in item.glob("*")):
                                shutil.rmtree(item, ignore_errors=True)
                        except Exception:
                            pass
                lair.safe_derive(ds_instance)

            fps = lair.get_dataset_filepaths(ds_instance)
            if fps:
                src_path = list(fps.values())[0]
                dataset_dir.mkdir(parents=True, exist_ok=True)
                import shutil
                shutil.copy2(src_path, dest_h5)
                return Success(dest_h5)
            else:
                return Failure(f"No filepaths found for {name}")
        else:
            return Failure(f"Unknown class {name} in singlecellrnasignature.adata")
    except Exception as e:
        return Failure(f"Failed to download {name}: {str(e)}")


def run_download_pipeline(config: DownloadConfig) -> Result[dict[str, tuple[Path, ...]], str]:
    """Downloads all 18 dataset files uniformly."""
    results = {}
    
    # 1. Download GEO HTTP datasets
    for spec in config.datasets:
        match download_dataset(spec, config.out_dir):
            case Success(paths):
                results[spec.accession] = paths
            case Failure(err):
                print(f"Warning downloading {spec.accession}: {err}")

    # 2. Download singlecellrnasignature datasets
    for sname in SINGLECELL_SIGNATURE_DATASETS:
        match fetch_signature_dataset(sname, config.out_dir):
            case Success(p):
                results[sname] = (p,)
            case Failure(err):
                print(f"Warning downloading {sname}: {err}")

    return Success(results)


def main() -> None:
    """CLI entry point for downloading single-cell datasets."""
    import argparse
    parser = argparse.ArgumentParser(description="Download Single-Cell Datasets.")
    parser.add_argument("--out-dir", type=str, default="data/raw_geo", help="Output directory for raw downloaded files")
    parser.add_argument("--dataset", type=str, default=None, help="Optional specific dataset accession/name to download")
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()

    if args.dataset:
        dname = args.dataset
        print(f"Downloading single dataset '{dname}' into: {out_dir}")
        
        # Check if GEO HTTP dataset
        config = DownloadConfig(out_dir=out_dir)
        geo_match = [spec for spec in config.datasets if spec.accession.lower() == dname.lower()]
        if geo_match:
            match download_dataset(geo_match[0], out_dir):
                case Success(paths):
                    print(f"Successfully downloaded GEO dataset '{dname}' -> {paths}")
                case Failure(err):
                    print(f"Error downloading GEO dataset '{dname}': {err}")
        else:
            match fetch_signature_dataset(dname, out_dir):
                case Success(p):
                    print(f"Successfully fetched singlecellrnasignature dataset '{dname}' -> {p}")
                case Failure(err):
                    print(f"Error fetching dataset '{dname}': {err}")
    else:
        config = DownloadConfig(out_dir=out_dir)
        print(f"Starting downloads for all datasets into: {config.out_dir}")

        match run_download_pipeline(config):
            case Success(res):
                print(f"\nSuccessfully processed {len(res)} single-cell datasets:")
                for acc, paths in res.items():
                    print(f"[{acc}] -> {len(paths)} files")
            case Failure(err):
                print(f"Download pipeline failed: {err}")


if __name__ == "__main__":
    main()
