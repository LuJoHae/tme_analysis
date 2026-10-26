"""CLI Script to download Tier 1 ICB single-cell datasets (disabled on CSCS Eiger)."""

import argparse
import os
from pathlib import Path
import sys
from returns.result import Success, Failure
from single_cell_immuno_datasets import DataDirectories, download_all_tier1, download_dataset, TIER_1_DATASETS
from tme_datasets.paths import is_cscs_environment


def main() -> None:
    if is_cscs_environment():
        print(
            "ERROR: Internet dataset downloads are strictly disabled on CSCS Alps Eiger.\n"
            "Download raw datasets on 'olm' first and transfer via Globus to /capstor/store/cscs/empa/em25/halu/data-test.",
            file=sys.stderr,
        )
        sys.exit(1)
    parser = argparse.ArgumentParser(description="Download Tier 1 ICB Single-Cell Datasets.")
    parser.add_argument(
        "--base-dir",
        type=str,
        default=os.environ.get("TME_DATA_DIR", "/ritom/scratch/cscs/lhaeuser/data-test"),
        help="Base directory",
    )
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
