"""CLI Script to download Tier 1 ICB single-cell datasets into /storage/halu/data/raw."""

import argparse
from pathlib import Path
from returns.result import Success, Failure
from single_cell_immuno_datasets import DataDirectories, download_all_tier1, download_dataset, TIER_1_DATASETS


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Tier 1 ICB Single-Cell Datasets.")
    parser.add_argument("--base-dir", type=str, default="/storage/halu/data", help="Base directory (default: /storage/halu/data)")
    parser.add_argument("--accession", type=str, default=None, help="Specific dataset accession to download")
    args = parser.parse_args()

    dirs = DataDirectories.with_base(Path(args.base_dir))

    if args.accession:
        spec = [s for s in TIER_1_DATASETS if s.accession.lower() == args.accession.lower()]
        if not spec:
            print(f"Error: Unknown accession '{args.accession}'. Available: {[s.accession for s in TIER_1_DATASETS]}")
            return
        print(f"Downloading single dataset '{spec[0].accession}' into {dirs.raw_dir}...")
        match download_dataset(spec[0], dirs.raw_dir):
            case Success(paths):
                print(f"Successfully downloaded {spec[0].accession} -> {paths}")
            case Failure(err):
                print(f"Download error: {err}")
    else:
        print(f"Starting download pipeline for all Tier 1 ICB datasets into {dirs.raw_dir}...")
        match download_all_tier1(dirs):
            case Success(res):
                print(f"Downloaded {len(res)} dataset specs into {dirs.raw_dir}")
            case Failure(err):
                print(f"Download pipeline failed: {err}")


if __name__ == "__main__":
    main()
