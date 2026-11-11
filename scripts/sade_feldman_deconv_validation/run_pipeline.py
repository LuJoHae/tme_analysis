#!/usr/bin/env python3
"""
Master Orchestration Wrapper for Multi-Resolution Sade-Feldman Deconvolution & Validation Pipeline.
Delegates steps 3f1..3f8 and 4d1..4d2 to scripts/manuscript_pipeline/ with isolated per-script output directories.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def run_command(cmd: list[str]) -> None:
    """Run a subprocess command and handle errors."""
    print(f"\n>>> Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"ERROR: Command failed with exit code {result.returncode}: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run end-to-end Multi-Resolution Sade-Feldman Deconvolution & Milopy Validation Pipeline."
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="",
        help="Optional explicit path to Sade-Feldman .h5ad (defaults to tme_datasets GSE120575)",
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
        default="",
        help="Deprecated optional path to legacy datalair directory",
    )
    parser.add_argument(
        "--pipeline-out-root",
        type=str,
        default="output/manuscript_pipeline",
        help="Common manuscript_pipeline output root directory",
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
    repo_root = Path(__file__).resolve().parents[2]
    mp_dir = repo_root / "scripts" / "manuscript_pipeline"
    out_root = Path(args.pipeline_out_root)

    out_3f1 = out_root / "3f1_build_sf_reference"
    out_3f2 = out_root / "3f2_build_integrated_reference"
    out_3f5 = out_root / "3f5_deconvolute_iatlas_multires"
    out_3f6 = out_root / "3f6_logistic_regression_multires"
    out_3f7 = out_root / "3f7_analyze_milopy_multires"
    out_3f8 = out_root / "3f8_compare_concordance_multires"
    out_4d1 = out_root / "4d1_plot_sf_multires_validation_figures"
    out_4d2 = out_root / "4d2_plot_cohort_auc_reference_selection"

    adata_flags = ["--adata", args.adata] if args.adata else []
    adata_sf_flags = ["--adata-sf", args.adata] if args.adata else []
    lair_flags = ["--lair-dir", args.lair_dir] if args.lair_dir else []

    run_command([
        python_bin,
        str(mp_dir / "3f1_build_sf_reference.py"),
        *adata_flags,
        "--out-dir", str(out_3f1),
    ])

    run_command([
        python_bin,
        str(mp_dir / "3f2_build_integrated_reference.py"),
        *adata_sf_flags,
        *lair_flags,
        "--out-dir", str(out_3f2),
    ])

    if not (args.skip_deconv and (out_3f5 / "deconv_fractions.parquet").exists()):
        run_command([
            python_bin,
            str(mp_dir / "3f5_deconvolute_iatlas_multires.py"),
            "--reference", str(out_3f1 / "reference_phi.parquet"),
            *lair_flags,
            "--out-dir", str(out_3f5),
            "--benchmark-mode",
        ])

    run_command([
        python_bin,
        str(mp_dir / "3f6_logistic_regression_multires.py"),
        "--fractions", str(out_3f5 / "deconv_fractions.parquet"),
        *lair_flags,
        "--out-dir", str(out_3f6),
        "--benchmark-mode",
    ])

    if not (args.skip_milo and (out_3f7 / "milopy_cell_state_da.parquet").exists()):
        run_command([
            python_bin,
            str(mp_dir / "3f7_analyze_milopy_multires.py"),
            *adata_flags,
            "--cluster-col", args.cluster_col,
            "--out-dir", str(out_3f7),
        ])

    run_command([
        python_bin,
        str(mp_dir / "3f8_compare_concordance_multires.py"),
        "--logistic-results", str(out_3f6 / "logistic_regression_results.parquet"),
        "--milo-dir", str(out_3f7),
        "--out-dir", str(out_3f8),
    ])

    run_command([
        python_bin,
        str(mp_dir / "4d1_plot_sf_multires_validation_figures.py"),
        "--data-dir", str(out_3f8),
        "--results-dir", str(out_4d1),
        "--stratum", "Melanoma",
    ])

    run_command([
        python_bin,
        str(mp_dir / "4d2_plot_cohort_auc_reference_selection.py"),
        "--data-dir", str(out_3f8),
        "--results-dir", str(out_4d2),
    ])


if __name__ == "__main__":
    main()
