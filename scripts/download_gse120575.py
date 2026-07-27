import sys
from pathlib import Path
from single_cell_immuno_datasets.gse120575 import fetch_and_format_gse120575
from returns.result import Success, Failure
import argparse

def main():
    parser = argparse.ArgumentParser(description="Download and format GSE120575 data")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args()

    print(f"Fetching GSE120575 to {args.out_dir}...")
    result = fetch_and_format_gse120575(args.out_dir)

    match result:
        case Success(paths):
            print(f"Successfully downloaded and formatted to {paths}")
            sys.exit(0)
        case Failure(err):
            print(f"Error downloading data: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
