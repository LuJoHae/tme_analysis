#!/usr/bin/env python3
"""Unified Publication Figures Generator for Stability Selection.

Orchestrates all 4 publication figure pipelines:
1. Individual cohort stability paths (MB vs SS-CPSS) with regularization budget cutoffs
2. 12-cohort composite figure (2x6 Nature Methods wireframe grid)
3. Cross-cohort gene selection overlap (UpSet + Dot-matrix + 4-set Euler partition)
4. Multi-fitter gene selection matrix (13 fitters, 2-tier capsule rack, cardinality)

Usage:
    python scripts/stability_selection_iatlas/plotting/run_all_plots.py --all
    python scripts/stability_selection_iatlas/plotting/run_all_plots.py --composite --fitter-overlap
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


def run_cmd(cmd: list[str], description: str) -> bool:
    """Execute a plotting command and report outcome."""
    print(f"\n{'='*70}")
    print(f" [Pipeline] {description}")
    print(f" Command:   {' '.join(cmd)}")
    print(f"{'='*70}")
    res = subprocess.run(cmd)
    if res.returncode == 0:
        print(f"[✓] {description} succeeded.")
        return True
    else:
        print(f"[✗] {description} failed with return code {res.returncode}.", file=sys.stderr)
        return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Unified Publication Figures Generator for Stability Selection"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate all 4 publication figure suites",
    )
    parser.add_argument(
        "--paths",
        action="store_true",
        help="Generate individual cohort stability paths (MB vs SS-CPSS)",
    )
    parser.add_argument(
        "--composite",
        action="store_true",
        help="Generate 12-cohort composite figure (2x6 Nature Methods wireframe)",
    )
    parser.add_argument(
        "--cohort-overlap",
        action="store_true",
        help="Generate cross-cohort gene selection overlap figure (UpSet + Euler)",
    )
    parser.add_argument(
        "--fitter-overlap",
        action="store_true",
        help="Generate 13-fitter 2-tier cross-cohort selection matrix figure",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Custom base output directory for figures",
    )

    args = parser.parse_args()

    # Default to --all if no specific flags were specified
    if not (args.all or args.paths or args.composite or args.cohort_overlap or args.fitter_overlap):
        args.all = True

    script_dir = Path(__file__).resolve().parent
    python_bin = sys.executable

    successes: list[str] = []
    failures: list[str] = []

    # 1. Individual Cohort Stability Paths
    if args.all or args.paths:
        cmd = [python_bin, str(script_dir / "plot_stability_selection_iatlas.py")]
        if args.output_dir is not None:
            cmd.extend(["--output-dir", str(args.output_dir)])
        ok = run_cmd(cmd, "Individual Cohort Stability Paths")
        (successes if ok else failures).append("Stability Paths")

    # 2. Composite 12-Cohort Figure
    if args.all or args.composite:
        cmd = [python_bin, str(script_dir / "combine_all_stability_paths.py")]
        if args.output_dir is not None:
            cmd.extend([
                "--output-svg", str(args.output_dir / "all_cohorts_stability_paths_composite.svg"),
                "--output-png", str(args.output_dir / "all_cohorts_stability_paths_composite.png"),
            ])
        ok = run_cmd(cmd, "12-Cohort Composite Figure")
        (successes if ok else failures).append("12-Cohort Composite")

    # 3. Cross-Cohort Overlap (UpSet + Euler)
    if args.all or args.cohort_overlap:
        cmd = [python_bin, str(script_dir / "plot_cohort_gene_overlap.py")]
        if args.output_dir is not None:
            cmd.extend([
                "--output-svg", str(args.output_dir / "cohort_gene_overlap_nature.svg"),
                "--output-png", str(args.output_dir / "cohort_gene_overlap_nature.png"),
            ])
        ok = run_cmd(cmd, "Cross-Cohort Gene Selection Overlap")
        (successes if ok else failures).append("Cohort Overlap")

    # 4. Multi-Fitter Gene Selection Matrix
    if args.all or args.fitter_overlap:
        cmd = [python_bin, str(script_dir / "plot_fitter_gene_overlap.py")]
        if args.output_dir is not None:
            cmd.extend([
                "--output-svg", str(args.output_dir / "fitter_gene_overlap_nature.svg"),
                "--output-png", str(args.output_dir / "fitter_gene_overlap_nature.png"),
            ])
        ok = run_cmd(cmd, "13-Fitter Gene Selection Overlap")
        (successes if ok else failures).append("Fitter Overlap")

    print("\n" + "=" * 70)
    print(" Figure Pipeline Summary")
    print(f" Succeeded: {len(successes)} ({', '.join(successes)})")
    if failures:
        print(f" Failed:    {len(failures)} ({', '.join(failures)})")
    print("=" * 70)

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
