"""CLI script to generate comparative analysis report and SVG charts for ICB single-cell datasets."""

import argparse
from pathlib import Path
from returns.result import Success, Failure
from single_cell_immuno_datasets import DataDirectories, generate_summary_report


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate ICB single-cell dataset summary report and SVG plots.")
    parser.add_argument("--base-dir", type=str, default="/storage/halu/data", help="Base directory (default: /storage/halu/data)")
    args = parser.parse_args()

    dirs = DataDirectories.with_base(Path(args.base_dir))
    print(f"Generating dataset summary report in {dirs.reports_dir}...")

    match generate_summary_report(dirs):
        case Success((df, report_md, (svg_patients, svg_tiers))):
            print(f"Successfully generated summary report:")
            print(f" - Summary Dataframe: {df.height} datasets evaluated")
            print(f" - Markdown Report: {report_md}")
            print(f" - Patient SVG Chart: {svg_patients}")
            print(f" - Tier Summary SVG Chart: {svg_tiers}")
        case Failure(err):
            print(f"Failed to generate summary report: {err}")


if __name__ == "__main__":
    main()
