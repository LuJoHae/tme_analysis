#!/usr/bin/env python3
"""Master Orchestration CLI for the 9-Cohort Multi-Resolution Deconvolution & Milo Validation Pipeline.
Orchestrates steps 01 through 07 with support for selective step execution and caching.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Final

DEFAULT_SC9_COHORTS: Final[tuple[str, ...]] = (
    "GSE120575",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE207422",
    "GSE243013",
    "GSE233203",
    "GSE200996",
    "GSE316195",
)


def run_command(cmd: list[str]) -> None:
    """Run a subprocess command and handle errors cleanly."""
    print(f"\n>>> Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"ERROR: Command failed with exit code {result.returncode}: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run end-to-end 9-Cohort Multi-Resolution Single-Cell Reference Deconvolution & Milo Validation Pipeline."
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default=",".join(DEFAULT_SC9_COHORTS),
        help="Comma-separated cohort accessions to include in reference construction",
    )
    parser.add_argument(
        "--pipeline-out-root",
        type=str,
        default="output/sc9_milo_deconv_validation",
        help="Root directory for intermediate pipeline outputs",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="results/sc9_milo_deconv_validation",
        help="Directory to save final publication figures",
    )
    parser.add_argument(
        "--skip-single-ref",
        action="store_true",
        help="Skip single-dataset reference building (Step 01)",
    )
    parser.add_argument(
        "--skip-comb-ref",
        action="store_true",
        help="Skip multi-dataset combined reference building (Step 02)",
    )
    parser.add_argument(
        "--skip-deconv",
        action="store_true",
        help="Skip bulk deconvolution step (Step 03)",
    )
    parser.add_argument(
        "--skip-lr",
        action="store_true",
        help="Skip logistic regression step (Step 04)",
    )
    parser.add_argument(
        "--skip-milo",
        action="store_true",
        help="Skip Milo DA step (Step 05)",
    )
    parser.add_argument(
        "--skip-conc",
        action="store_true",
        help="Skip concordance and benchmark table compilation (Step 06)",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip publication figure generation (Step 07)",
    )
    args = parser.parse_args()

    python_bin = sys.executable
    script_dir = Path(__file__).resolve().parent
    out_root = Path(args.pipeline_out_root)
    results_dir = Path(args.results_dir)

    dir_01_ref = out_root / "01_references"
    dir_02_deconv = out_root / "02_deconvolution"
    dir_03_lr = out_root / "03_logistic_regression"
    dir_04_milo = out_root / "04_milopy"
    dir_05_conc = out_root / "05_concordance"

    cohort_list = [c.strip() for c in args.cohorts.split(",") if c.strip()]

    # Step 01: Single-dataset reference construction
    if not args.skip_single_ref:
        print("\n========================================================")
        print("STEP 01: Building Single-Dataset References")
        print("========================================================")
        for cohort in cohort_list:
            run_command([
                python_bin,
                str(script_dir / "01_build_single_reference.py"),
                "--cohort", cohort,
                "--out-dir", str(dir_01_ref),
            ])

    # Step 02: Multi-dataset combined references
    if not args.skip_comb_ref:
        print("\n========================================================")
        print("STEP 02: Building Combined References via Harmony")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "02_build_combined_references.py"),
            "--combo", "all",
            "--out-dir", str(dir_01_ref),
        ])

    # Step 03: Bulk deconvolution with InstaPrism
    if not args.skip_deconv:
        print("\n========================================================")
        print("STEP 03: Deconvoluting iAtlas Bulk Cohorts")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "03_deconvolute_iatlas.py"),
            "--ref-dir", str(dir_01_ref),
            "--out-dir", str(dir_02_deconv),
        ])

    # Step 04: Logistic regression and predictive performance evaluation
    if not args.skip_lr:
        print("\n========================================================")
        print("STEP 04: Fitting Logistic Regression Models")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "04_logistic_regression.py"),
            "--deconv-dir", str(dir_02_deconv),
            "--out-dir", str(dir_03_lr),
        ])

    # Step 05: Single-cell Milo DA aggregation
    if not args.skip_milo:
        print("\n========================================================")
        print("STEP 05: Aggregating Milo DA across Resolutions")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "05_analyze_milopy.py"),
            "--out-dir", str(dir_04_milo),
        ])

    # Step 06: Concordance and unified benchmark summary table
    if not args.skip_conc:
        print("\n========================================================")
        print("STEP 06: Compiling Master Benchmark Summary Table")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "06_compare_concordance.py"),
            "--ref-dir", str(dir_01_ref),
            "--logistic-dir", str(dir_03_lr),
            "--milo-dir", str(dir_04_milo),
            "--out-dir", str(dir_05_conc),
        ])

    # Step 07: Publication figure plotting
    if not args.skip_plots:
        print("\n========================================================")
        print("STEP 07: Generating Publication Figures")
        print("========================================================")
        run_command([
            python_bin,
            str(script_dir / "07_plot_figures.py"),
            "--data-dir", str(dir_05_conc),
            "--results-dir", str(results_dir),
        ])

    print("\n========================================================")
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print(f"Intermediate outputs: {out_root}")
    print(f"Publication figures:  {results_dir}")
    print("========================================================")


if __name__ == "__main__":
    main()
