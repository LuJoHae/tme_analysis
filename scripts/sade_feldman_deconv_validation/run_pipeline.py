#!/usr/bin/env python3
"""
Master Orchestration Script for Sade-Feldman Deconvolution & Validation Pipeline.
Runs steps 1 to 6 sequentially, ensuring end-to-end reproducibility.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_command(cmd: list[str]) -> None:
    """Run a subprocess command and handle errors."""
    print(f"\n>>> Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"ERROR: Command failed with exit code {result.returncode}: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run end-to-end Sade-Feldman Deconvolution & Milopy Validation Pipeline."
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="/storage/halu/data/GSE120575/gse120575_processed.h5ad",
        help="Path to gse120575_processed.h5ad",
    )
    parser.add_argument(
        "--cluster-col",
        type=str,
        default="celltypist_leiden_0.5",
        help="Cluster column name",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="/storage/halu/lair",
        help="Path to datalair directory",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Output directory for data parquets",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="results/sade_feldman_deconv_validation",
        help="Output directory for SVG figures",
    )
    parser.add_argument(
        "--skip-deconv",
        action="store_true",
        help="Skip deconvolution step if deconv_fractions.parquet already exists",
    )
    parser.add_argument(
        "--skip-milo",
        action="store_true",
        help="Skip Milo step if milo parquet already exists",
    )
    args = parser.parse_args()

    python_bin = sys.executable
    script_dir = Path(__file__).resolve().parent

    print("======================================================================")
    print("STARTING SADE-FELDMAN DECONVOLUTION & MILOPY VALIDATION PIPELINE")
    print("======================================================================")

    # Step 1: Reference preparation
    print("\n--- STEP 1: Building Reference from Sade-Feldman Dataset ---")
    run_command([
        python_bin,
        str(script_dir / "01_build_reference.py"),
        "--adata", args.adata,
        "--cluster-col", args.cluster_col,
        "--out-dir", args.out_dir,
    ])

    # Step 2: Deconvolution of iAtlas datasets
    if not (args.skip_deconv and (Path(args.out_dir) / "deconv_fractions.parquet").exists()):
        print("\n--- STEP 2: Deconvoluting iAtlas Cohorts with instaprism ---")
        run_command([
            python_bin,
            str(script_dir / "02_deconvolute_iatlas.py"),
            "--reference", str(Path(args.out_dir) / "reference_phi.parquet"),
            "--lair-dir", args.lair_dir,
            "--out-dir", args.out_dir,
        ])
    else:
        print("\n--- STEP 2: Skipping deconvolution (cached) ---")

    # Step 3: Logistic regression of fractions vs response
    print("\n--- STEP 3: Logistic Regression on Deconvoluted Fractions vs. Response ---")
    run_command([
        python_bin,
        str(script_dir / "03_logistic_regression.py"),
        "--fractions", str(Path(args.out_dir) / "deconv_fractions.parquet"),
        "--lair-dir", args.lair_dir,
        "--out-dir", args.out_dir,
    ])

    # Step 4: Milopy differential abundance analysis
    if not (args.skip_milo and (Path(args.out_dir) / "milopy_cell_state_da.parquet").exists()):
        print("\n--- STEP 4: Milopy Differential Abundance Analysis on GSE120575 ---")
        run_command([
            python_bin,
            str(script_dir / "04_analyze_milopy.py"),
            "--adata", args.adata,
            "--cluster-col", args.cluster_col,
            "--out-dir", args.out_dir,
        ])
    else:
        print("\n--- STEP 4: Skipping Milopy analysis (cached) ---")

    # Step 5: Concordance comparison
    print("\n--- STEP 5: Concordance Analysis (Bulk Deconv vs. Single-Cell Milo) ---")
    run_command([
        python_bin,
        str(script_dir / "05_compare_concordance.py"),
        "--logistic-results", str(Path(args.out_dir) / "logistic_regression_results.parquet"),
        "--milo-results", str(Path(args.out_dir) / "milopy_cell_state_da.parquet"),
        "--stratum", "Melanoma",
        "--out-dir", args.out_dir,
    ])

    # Step 6: Visualization
    print("\n--- STEP 6: Generating Publication Figures (Altair Vector SVGs) ---")
    run_command([
        python_bin,
        str(script_dir / "06_plot_figures.py"),
        "--data-dir", args.out_dir,
        "--results-dir", args.results_dir,
        "--stratum", "Melanoma",
    ])

    print("\n======================================================================")
    print("PIPELINE EXECUTION COMPLETE!")
    print(f"Data Parquets: {args.out_dir}")
    print(f"Visualizations: {args.results_dir}")
    print("======================================================================")


if __name__ == "__main__":
    main()
