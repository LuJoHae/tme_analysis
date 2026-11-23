#!/usr/bin/env python3
"""3e1b_benchmark_iatlas_perturbations.py: Perturbation & robustness benchmarks on iAtlas cohorts.

Runs systematic data perturbation sweeps (jitter, dropout, dilution, label noise)
and computes Predictor Perturbation Resilience Index (PRI) metrics across
cBioPortal iAtlas cohorts, exporting `.parquet` and `.csv` tables.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import polars as pl
from returns.result import Failure, Success

from tme_datasets.paths import find_repo_root
from tme_response import (
    compute_perturbation_resilience,
    load_multi_omic_cohort,
    run_cohort_perturbation_sweep,
)

DEFAULT_IATLAS_COHORTS = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "McDermott-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "Choueiri-iAtlas",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run perturbation & resilience benchmarks on iAtlas ICB cohorts."
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default=",".join(DEFAULT_IATLAS_COHORTS),
        help="Comma-separated cohorts to stress-test under perturbations.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3e1_benchmark_iatlas_icb_predictors"),
        help="Output directory for perturbation .parquet and .csv tables.",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=200,
        help="Number of bootstrap iterations per perturbation point (default: 200).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force recomputation even if cached perturbation evaluations exist.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = find_repo_root()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    pert_cohort_ids = [c.strip() for c in args.cohorts.split(",") if c.strip()]
    print("=== Running Systematic Data Perturbation & Robustness Benchmarks ===")
    print(f"Target Cohorts ({len(pert_cohort_ids)}): {', '.join(pert_cohort_ids)}")
    print(f"Output Directory: {output_dir}")

    pert_pq = output_dir / "perturbation_benchmark_metrics.parquet"
    pert_csv = output_dir / "perturbation_benchmark_metrics.csv"

    cached_records: list[dict[str, object]] = []
    already_computed: set[str] = set()
    if pert_pq.exists() and not args.force:
        cached_df = pl.read_parquet(pert_pq)
        cached_records = cached_df.to_dicts()
        already_computed = set(cached_df["cohort_id"].unique().to_list())
        print(
            f"Loaded {len(cached_df)} existing perturbation evaluations across "
            f"{len(already_computed)} cohorts: {', '.join(sorted(already_computed))}"
        )

    new_pert_records: list[dict[str, object]] = []
    for p_cid in pert_cohort_ids:
        if p_cid in already_computed:
            print(f"\n--- Cohort {p_cid}: Already evaluated in cache; skipping recomputation ---")
            continue

        print(f"\n--- Perturbation Sweep on: {p_cid} ---")
        match load_multi_omic_cohort(p_cid, repo_root=repo_root, filter_baseline=True):
            case Failure(err):
                print(f"  [SKIPPED] Failed loading {p_cid}: {err}")
                continue
            case Success(p_cohort):
                match run_cohort_perturbation_sweep(
                    p_cohort, n_bootstrap=args.n_bootstrap
                ):
                    case Success(pert_df):
                        print(f"  Generated {len(pert_df)} perturbation evaluations on {p_cid}.")
                        new_pert_records.extend(pert_df.to_dicts())
                    case Failure(err):
                        print(f"  Perturbation sweep failed on {p_cid}: {err}")

    all_pert_records = cached_records + new_pert_records
    if not all_pert_records:
        print("\n[WARNING] No perturbation evaluations generated.")
        return

    pert_benchmark_df = pl.DataFrame(all_pert_records)
    pert_benchmark_df.write_parquet(pert_pq)
    pert_benchmark_df.write_csv(pert_csv)
    print(f"Saved: {pert_pq}")
    print(f"Saved: {pert_csv}")

    match compute_perturbation_resilience(pert_benchmark_df):
        case Success(resilience_df):
            resil_pq = output_dir / "perturbation_resilience_metrics.parquet"
            resil_csv = output_dir / "perturbation_resilience_metrics.csv"
            resilience_df.write_parquet(resil_pq)
            resilience_df.write_csv(resil_csv)
            print(f"Saved: {resil_pq}")
            print(f"Saved: {resil_csv}")
        case Failure(err):
            print(f"Failed to compute resilience index: {err}")

    print("\n=== Perturbation Benchmark Finished Successfully ===")


if __name__ == "__main__":
    main()
