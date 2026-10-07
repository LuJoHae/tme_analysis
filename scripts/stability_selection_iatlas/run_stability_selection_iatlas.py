#!/usr/bin/env python3
"""Unified Stability Selection Runner (Single-Fitter & Multi-Fitter Matrix Benchmark).

Executes finite-sample stability selection comparing:
- Meinshausen & Bühlmann (2010) [MB]
- Shah & Samworth Complementary Pairs Stability Selection (2013) [SS-CPSS]

Key Features:
- Comprehensive 13-fitter benchmark matrix across cohorts (by default or via --fitters all).
- Flexible single-fitter or custom-fitter subsets (e.g. --fitters lasso, or --fitters oscar slope).
- Curated immunotherapy panel (101 genes) or transcriptome-wide HVG selection (--feature-mode hvg).
- Granular per-(cohort, fitter) caching in `.cache/` to eliminate redundant recomputations.
- In-memory AnnData caching to prevent repeated disk I/O.
- Enforces regularization budget q_budget to prevent PFER explosion and false discovery flooding.
- Columnar Parquet serialization with Polars.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
from typing import Final, Sequence
import anndata as ad  # type: ignore[import-untyped]
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

# Add parent directory to sys.path to enable clean imports when run as script
sys.path.append(str(Path(__file__).resolve().parent))

from shared import (
    COMBINED_COHORTS,
    DEFAULT_COHORTS,
    FITTERS_FOR_COMBINED,
    FITTERS_FOR_SINGLE,
    CohortSelectionOutput,
    assemble_parquet_tables,
    load_cohort_adata,
    run_cohort_analysis,
)
from tme_datasets import (  # type: ignore[import-untyped]
    IMMUNE_CHECKPOINT_GENES,
    IMMUNOTHERAPY_GENE_PANEL,
)

ALL_AVAILABLE_FITTERS: Final[tuple[str, ...]] = (
    "lasso",
    "elastic_net",
    "logistic",
    "rf",
    "cohort_adjusted",
    "group_lasso",
    "merf",
    "multitask_logistic",
    "meta_analysis",
    "multistudy_invariant",
    "glmm_lasso",
    "oscar",
    "slope",
)

STRATIFIED_FITTERS: Final[set[str]] = {
    "cohort_adjusted",
    "group_lasso",
    "merf",
    "multitask_logistic",
    "meta_analysis",
    "multistudy_invariant",
    "glmm_lasso",
    "oscar",
    "slope",
}


def resolve_fitters_for_cohort(
    cohort_id: str,
    requested_fitters: Sequence[str],
) -> tuple[tuple[str, bool], ...]:
    """Resolve (fitter_name, is_stratified) pairs for a given cohort."""
    is_combined = cohort_id in COMBINED_COHORTS
    if "all" in requested_fitters:
        return FITTERS_FOR_COMBINED if is_combined else FITTERS_FOR_SINGLE

    resolved: list[tuple[str, bool]] = []
    for f in requested_fitters:
        f_clean = f.lower()
        if f_clean not in ALL_AVAILABLE_FITTERS:
            print(f"[-] Warning: Unknown fitter '{f_clean}'. Skipping.")
            continue
        is_strat = (f_clean in STRATIFIED_FITTERS) if is_combined else False
        resolved.append((f_clean, is_strat))

    return tuple(resolved)


def run_pipeline(
    cohorts: Sequence[str],
    output_dir: Path,
    fitters: Sequence[str] = ("all",),
    feature_mode: str = "immunotherapy",
    n_top_genes: int = 500,
    pfer: float = 1.0,
    cutoff: float = 0.75,
    q_budget: float = 20.0,
    B: int = 50,
    seed: int = 42,
    l1_ratio: float = 0.7,
    kappa: float = 0.5,
    q_fdr: float = 0.1,
    force_stratified: bool | None = None,
    overwrite: bool = False,
) -> Result[tuple[Path, Path, Path], str]:
    """Execute stability selection runs across cohorts and fitters with caching."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    scores_dfs: list[pl.DataFrame] = []
    paths_dfs: list[pl.DataFrame] = []
    summary_dfs: list[pl.DataFrame] = []

    budget_q: float | None = q_budget if (q_budget is not None and q_budget > 0) else None

    print("==================================================================")
    print(" Unified Stability Selection Pipeline")
    print(f" Feature Mode:     {feature_mode.upper()} ({'101 IO genes' if feature_mode == 'immunotherapy' else f'Top {n_top_genes} HVGs'})")
    print(f" Cohorts ({len(cohorts)}):     {list(cohorts)}")
    print(f" Fitters:          {list(fitters)}")
    print(f" Parameters:       pi_thr={cutoff}, PFER<={pfer}, q_budget<={budget_q}, B={B}, seed={seed}")
    print(f" Fitter Hyperparams: l1_ratio={l1_ratio}, kappa={kappa}, q_fdr={q_fdr}")
    print(f" Target Directory: {output_dir}")
    print(f" Overwrite Cache:  {'YES' if overwrite else 'NO'}")
    print("==================================================================")

    start_time = time.time()
    total_runs = 0
    adata_cache: dict[str, ad.AnnData] = {}

    for cohort_id in cohorts:
        is_combined = cohort_id in COMBINED_COHORTS
        fitters_to_run = resolve_fitters_for_cohort(cohort_id, fitters)
        if not fitters_to_run:
            print(f"\n[-] No valid fitters configured for cohort '{cohort_id}'. Skipping.")
            continue

        print(f"\n[+] Processing cohort: '{cohort_id}' (is_combined={is_combined}, {len(fitters_to_run)} fitters)...")

        for fitter_name, default_stratified in fitters_to_run:
            total_runs += 1
            is_stratified = default_stratified if force_stratified is None else force_stratified

            prefix = f"{cohort_id}__{fitter_name}__{feature_mode}"
            c_scores_path = cache_dir / f"{prefix}__scores.parquet"
            c_paths_path = cache_dir / f"{prefix}__paths.parquet"
            c_summary_path = cache_dir / f"{prefix}__summary.parquet"

            legacy_prefix = f"{cohort_id}__{fitter_name}"
            l_scores_path = cache_dir / f"{legacy_prefix}__scores.parquet"
            l_paths_path = cache_dir / f"{legacy_prefix}__paths.parquet"
            l_summary_path = cache_dir / f"{legacy_prefix}__summary.parquet"

            target_scores = c_scores_path if c_scores_path.exists() else l_scores_path
            target_paths = c_paths_path if c_paths_path.exists() else l_paths_path
            target_summary = c_summary_path if c_summary_path.exists() else l_summary_path

            if (
                not overwrite
                and target_scores.exists()
                and target_paths.exists()
                and target_summary.exists()
            ):
                print(f"    --> Fitter '{fitter_name.upper()}': [LOADED FROM CACHE]")
                scores_dfs.append(pl.read_parquet(target_scores))
                paths_dfs.append(pl.read_parquet(target_paths))
                summary_dfs.append(pl.read_parquet(target_summary))
                continue

            # Load AnnData if not yet cached in memory for this cohort
            if cohort_id not in adata_cache:
                print(f"    [+] Loading AnnData for '{cohort_id}' via tme_datasets...")
                load_res = load_cohort_adata(cohort_id)
                match load_res:
                    case Failure(err):
                        print(f"    [-] Skipping '{cohort_id}': Failed to load via tme_datasets: {err}")
                        break
                    case Success(loaded_ad):
                        adata_cache[cohort_id] = loaded_ad
                        print(f"    Loaded {loaded_ad.n_obs} samples x {loaded_ad.n_vars} genes.")

            adata_obj = adata_cache[cohort_id]
            fitter_start = time.time()
            print(f"    --> Running Fitter '{fitter_name.upper()}' (stratified={is_stratified})...", end="", flush=True)

            analysis_res = run_cohort_analysis(
                cohort_name=cohort_id,
                adata=adata_obj,
                feature_mode=feature_mode,
                n_top_genes=n_top_genes,
                pfer=pfer,
                cutoff=cutoff,
                B=B,
                seed=seed,
                max_expected_q=budget_q,
                fitter_name=fitter_name,
                l1_ratio=l1_ratio,
                stratified=is_stratified,
                kappa=kappa,
                q_fdr=q_fdr,
            )

            elapsed = time.time() - fitter_start
            match analysis_res:
                case Failure(err):
                    print(f" [FAILED in {elapsed:.1f}s: {err}]")
                case Success(out):
                    mb_sel = len(out.mb_result.selected_features)
                    ss_sel = len(out.ss_result.selected_features)
                    print(f" [OK in {elapsed:.1f}s | MB sel: {mb_sel}, SS sel: {ss_sel}]")

                    cur_paths, cur_scores, cur_summary = assemble_parquet_tables([out])
                    cur_scores.write_parquet(c_scores_path)
                    cur_paths.write_parquet(c_paths_path)
                    cur_summary.write_parquet(c_summary_path)

                    scores_dfs.append(cur_scores)
                    paths_dfs.append(cur_paths)
                    summary_dfs.append(cur_summary)

    if not scores_dfs:
        return Failure("No stability selection pipeline runs succeeded.")

    print(f"\n[+] Assembled {len(scores_dfs)} of {total_runs} pipeline runs in {time.time() - start_time:.1f}s.")
    print("[+] Concatenating and writing consolidated Parquet datasets...")

    final_scores_df = pl.concat(scores_dfs)
    final_paths_df = pl.concat(paths_dfs)
    final_summary_df = pl.concat(summary_dfs)

    paths_file = output_dir / "stability_paths.parquet"
    scores_file = output_dir / "stability_scores.parquet"
    summary_file = output_dir / "cohort_fitter_summary.parquet"
    summary_legacy_file = output_dir / "cohort_summary.parquet"

    final_paths_df.write_parquet(paths_file)
    final_scores_df.write_parquet(scores_file)
    final_summary_df.write_parquet(summary_file)
    final_summary_df.write_parquet(summary_legacy_file)

    print(f"    [x] Saved scores:   {scores_file} ({final_scores_df.shape[0]:,} rows)")
    print(f"    [x] Saved paths:    {paths_file} ({final_paths_df.shape[0]:,} rows)")
    print(f"    [x] Saved summary:  {summary_file} ({final_summary_df.shape[0]:,} rows)")
    print(f"    [x] Saved summary*: {summary_legacy_file}")

    return Success((scores_file, paths_file, summary_file))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified Stability Selection Runner (Single-Fitter & Multi-Fitter Matrix Benchmark)"
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=list(DEFAULT_COHORTS),
        help=f"List of cohort IDs or group names (default: {list(DEFAULT_COHORTS)})",
    )
    parser.add_argument(
        "--fitters",
        "--fitter",
        nargs="+",
        default=["all"],
        help="One or more fitters to run, or 'all' to evaluate full benchmark matrix (default: ['all'])",
    )
    parser.add_argument(
        "--feature-mode",
        choices=["immunotherapy", "hvg"],
        default="immunotherapy",
        help="Feature selection mode: 'immunotherapy' for curated IO panel, 'hvg' for top variable genes (default: immunotherapy)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=500,
        help="Number of highly variable genes to select if feature-mode is 'hvg' (default: 500)",
    )
    parser.add_argument(
        "--pfer",
        type=float,
        default=1.0,
        help="Target Per-Family Error Rate upper bound (default: 1.0)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.75,
        help="Stability selection threshold pi_thr in [0.5, 1.0] (default: 0.75)",
    )
    parser.add_argument(
        "--q-budget",
        type=float,
        default=20.0,
        help="Maximum expected model size E[|S(lambda)|] along the path (default: 20.0). Set <= 0 to use theoretical analytical budget.",
    )
    parser.add_argument(
        "--B",
        type=int,
        default=50,
        help="Number of subsample complementary pairs (default: 50)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    parser.add_argument(
        "--l1-ratio",
        type=float,
        default=0.7,
        help="L1 ratio mixing parameter for Elastic Net fitter (default: 0.7)",
    )
    parser.add_argument(
        "--kappa",
        type=float,
        default=0.5,
        help="Clustering parameter for OSCAR fitter in [0, 1] (default: 0.5)",
    )
    parser.add_argument(
        "--q-fdr",
        type=float,
        default=0.1,
        help="Target FDR parameter for SLOPE fitter in (0, 1) (default: 0.1)",
    )
    parser.add_argument(
        "--stratified",
        action="store_true",
        default=None,
        help="Explicitly enable stratified subsampling (default: automatic per fitter)",
    )
    parser.add_argument(
        "--force-overwrite",
        action="store_true",
        default=False,
        help="Force recomputation of all selected cohorts and fitters ignoring cache",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/stability_selection_fitters_benchmark"),
        help="Target output directory (default: output/stability_selection_fitters_benchmark)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    res = run_pipeline(
        cohorts=args.cohorts,
        output_dir=args.output_dir,
        fitters=args.fitters,
        feature_mode=args.feature_mode,
        n_top_genes=args.n_top_genes,
        pfer=args.pfer,
        cutoff=args.cutoff,
        q_budget=args.q_budget,
        B=args.B,
        seed=args.seed,
        l1_ratio=args.l1_ratio,
        kappa=args.kappa,
        q_fdr=args.q_fdr,
        force_stratified=args.stratified,
        overwrite=args.force_overwrite,
    )

    match res:
        case Failure(err):
            print(f"\n[-] Error: {err}", file=sys.stderr)
            return 1
        case Success((scores_path, paths_path, summary_path)):
            print("\n==================================================================")
            print(f"[✓] Stability Selection pipeline completed successfully.")
            print(f"    Scores:  {scores_path}")
            print(f"    Paths:   {paths_path}")
            print(f"    Summary: {summary_path}")
            print("==================================================================")
            return 0


if __name__ == "__main__":
    sys.exit(main())
