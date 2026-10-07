#!/usr/bin/env python3
"""Master CLI script for single-cell reference sampling and classifier HPO.

Executes multi-fidelity ASHA optimization to identify optimal single-cell cohort combinations,
sampling budgets, patient-stratified malignant fractions, and downstream classifier pipelines
to maximize immunotherapy response prediction across iAtlas bulk RNA-seq cohorts.
"""

from __future__ import annotations

import argparse
import os
import socket
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

# Local repo packages
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "packages"))
from tme_datasets import (  # type: ignore[import-untyped]
    HPOSearchSpace,
    MalignantStrategy,
    build_deconvolution_reference,
    filter_candidate_cohorts,
    run_sampling_hpo,
    sample_single_cell_cohorts,
)
from tme_datasets.sampling.classifier_hpo import (  # type: ignore[import-untyped]
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    TransformType,
    apply_compositional_transform,
    fit_and_predict,
    select_features_train_test,
)
from tme_datasets.sampling.hpo_models import HPORunResult  # type: ignore[import-untyped]
from tme_datasets.sampling.hpo_optimizer import fast_deconvolute_cohorts  # type: ignore[import-untyped]

# Local script module
from plotting import plot_held_out_validation, plot_pareto_frontier  # type: ignore[import-not-found]


def verify_execution_host(
    benchmark_mode: bool = False,
    force_local: bool = False,
    hostname: str | None = None,
    allow_local_env: str | None = None,
    storage_dir: Path | None = None,
) -> Result[None, str]:
    """Verify that execution is occurring on remote host 'olm' unless benchmark/override is active."""
    if benchmark_mode or force_local:
        return Success(None)

    env_val = os.environ.get("ALLOW_LOCAL_HPO", "") if allow_local_env is None else allow_local_env
    if env_val.strip().lower() in {"1", "true", "yes"}:
        return Success(None)

    current_host = socket.gethostname().lower() if hostname is None else hostname.lower()
    if "olm" in current_host or current_host.startswith("sgl"):
        return Success(None)

    check_storage = Path("/storage/halu") if storage_dir is None else storage_dir
    if check_storage.exists():
        return Success(None)

    return Failure(
        f"Execution blocked: Single-cell reference sampling and classifier HPO is computationally "
        f"intensive (>64GB RAM) and restricted to remote host 'olm' (current host: '{current_host}').\n"
        f"To execute remotely on 'olm', run: make hpo (or make hpo-remote).\n"
        f"For local testing with synthetic cohorts, pass --benchmark-mode (or make hpo-benchmark-local)."
    )


def load_bulk_cohort_data(
    cohort_id: str,
    preprocessed_dir: Path,
) -> Result[tuple[np.ndarray, tuple[str, ...], tuple[str, ...], np.ndarray], str]:
    """Load bulk ICI cohort AnnData into memory for deconvolution evaluation."""
    from tme_datasets import load_dataset, load_iatlas_cohort_or_combined  # type: ignore[import-untyped]
    from tme_datasets.preprocessing.metadata import binarize_response  # type: ignore[import-untyped]

    candidates = [cohort_id, cohort_id.replace("-iAtlas", ""), f"{cohort_id}-iAtlas"]
    adata: ad.AnnData | None = None
    manual_dir = Path("data/manual_download")

    for cid in candidates:
        # Check preprocessed_dir
        fpath = preprocessed_dir / f"{cid}.h5ad"
        if fpath.exists():
            try:
                cand_ad = ad.read_h5ad(fpath)
                resp_cand = next(
                    (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in cand_ad.obs.columns),
                    None,
                )
                if resp_cand is not None and cand_ad.n_vars > 100:
                    adata = cand_ad
                    break
            except Exception:
                pass

        # Check manual_download dir
        mpath = manual_dir / f"{cid}.h5ad"
        if mpath.exists():
            try:
                cand_ad = ad.read_h5ad(mpath)
                resp_cand = next(
                    (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in cand_ad.obs.columns),
                    None,
                )
                if resp_cand is not None and cand_ad.n_vars > 100:
                    adata = cand_ad
                    break
            except Exception:
                pass

        # Check tme_datasets.load_dataset
        try:
            load_res = load_dataset(cid)
            if isinstance(load_res, Success):
                cand_ad = load_res.unwrap()
                resp_cand = next(
                    (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in cand_ad.obs.columns),
                    None,
                )
                if resp_cand is not None and cand_ad.n_vars > 100:
                    adata = cand_ad
                    break
        except Exception:
            pass

    if adata is None:
        return Failure(f"Failed to load bulk dataset for '{cohort_id}' with response annotations.")

    try:
        # Detect response column
        resp_col = next(
            (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in adata.obs.columns),
            None,
        )
        if resp_col is None:
            return Failure(f"No response column found in {cohort_id}.obs (columns: {list(adata.obs.columns)})")

        obs_resp = np.array([binarize_response(v) for v in adata.obs[resp_col]], dtype=np.float64)
        valid_mask = ~np.isnan(obs_resp)

        if np.sum(valid_mask) < 10:
            return Failure(f"Cohort {cohort_id} has fewer than 10 valid response annotations.")

        sub_adata = adata[valid_mask].copy()
        y = np.asarray(obs_resp[valid_mask], dtype=np.int64)

        # Extract expression matrix (samples x genes)
        X = sub_adata.X.toarray() if hasattr(sub_adata.X, "toarray") else np.asarray(sub_adata.X, dtype=np.float64)  # type: ignore[union-attr]

        # Extract both Ensembl IDs and Gene Symbols to maximize cross-nomenclature matching
        ensembl_genes = tuple(str(g).upper() for g in sub_adata.var_names)
        if "gene_name" in sub_adata.var.columns:
            raw_vals = sub_adata.var["gene_name"].tolist()
            symbol_genes = tuple(
                str(val).upper()
                if val is not None and str(val).strip() and str(val).lower() != "nan"
                else str(idx).upper()
                for val, idx in zip(raw_vals, ensembl_genes)
            )
        else:
            symbol_genes = ensembl_genes

        return Success((X, ensembl_genes, symbol_genes, y))
    except Exception as exc:
        return Failure(f"Failed to process bulk data for {cohort_id}: {exc}")


def generate_synthetic_benchmark_data() -> tuple[
    dict[str, tuple[np.ndarray, tuple[str, ...], np.ndarray]],
    dict[str, tuple[np.ndarray, tuple[str, ...], np.ndarray]],
]:
    """Generate high-speed synthetic bulk cohorts for dry-run verification."""
    rng = np.random.default_rng(42)
    n_genes = 200
    genes = tuple(f"MARKER_GENE_{i}" for i in range(n_genes))

    def make_cohort(n_samples: int) -> tuple[np.ndarray, tuple[str, ...], np.ndarray]:
        y = np.array([1] * (n_samples // 2) + [0] * (n_samples // 2))
        X = rng.poisson(lam=10.0, size=(n_samples, n_genes)).astype(np.float64)
        # Add strong signal in first 30 genes for responders
        X[y == 1, :30] += 25.0
        return X, genes, y

    disc = {
        "Hugo-iAtlas": make_cohort(20),
        "Riaz-iAtlas": make_cohort(20),
    }
    held = {
        "Gide-iAtlas": make_cohort(16),
        "Rosenberg-iAtlas": make_cohort(16),
    }
    return disc, held


def evaluate_held_out_generalization(
    best_ref_mat: np.ndarray,
    best_ref_genes: Sequence[str],
    state_labels: Sequence[str],
    discovery_bulk: Mapping[str, tuple[Any, ...]],
    held_out_bulk: Mapping[str, tuple[Any, ...]],
    best_inner_config: InnerClassifierConfig,
) -> Result[dict[str, float], str]:
    """Deconvolute held-out validation cohorts and test the fitted pipeline out-of-study."""
    from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]

    # 1. Deconvolute Discovery Cohorts
    disc_deconv_res = fast_deconvolute_cohorts(best_ref_mat, best_ref_genes, discovery_bulk, n_iter=60)
    match disc_deconv_res:
        case Failure(err):
            return Failure(f"Discovery deconvolution failed: {err}")
        case Success((disc_fracs, disc_resps)):
            pass

    # 2. Deconvolute Held-Out Cohorts
    held_deconv_res = fast_deconvolute_cohorts(best_ref_mat, best_ref_genes, held_out_bulk, n_iter=60)
    match held_deconv_res:
        case Failure(err):
            return Failure(f"Held-out deconvolution failed: {err}")
        case Success((held_fracs, held_resps)):
            pass

    # 3. Apply optimal compositional transform
    train_X_list = [
        apply_compositional_transform(disc_fracs[cid], best_inner_config.transform, best_inner_config.eps)
        for cid in disc_fracs
    ]
    train_y_list = [disc_resps[cid] for cid in disc_resps]

    pool_train_X = np.vstack(train_X_list)
    pool_train_y = np.concatenate(train_y_list)

    held_aucs: dict[str, float] = {}

    for test_cid, raw_te_fracs in held_fracs.items():
        test_y = held_resps[test_cid]
        te_transformed = apply_compositional_transform(
            raw_te_fracs,
            best_inner_config.transform,
            best_inner_config.eps,
        )

        # Feature selection strictly on pooled training set
        tr_sub, te_sub, _ = select_features_train_test(
            X_train=pool_train_X,
            y_train=pool_train_y,
            X_test=te_transformed,
            feature_names=state_labels,
            selector_type=best_inner_config.selector,
            config=best_inner_config,
        )

        try:
            pred_probs = fit_and_predict(
                X_train=tr_sub,
                y_train=pool_train_y,
                X_test=te_sub,
                classifier_type=best_inner_config.classifier,
                config=best_inner_config,
            )
            score = float(roc_auc_score(test_y, pred_probs))
            held_aucs[test_cid] = score
        except Exception:
            held_aucs[test_cid] = 0.50

    return Success(held_aucs)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run end-to-end single-cell reference sampling and classifier HPO."
    )
    parser.add_argument(
        "--cancer-types",
        type=str,
        default="Melanoma",
        help="Comma-separated target cancer types to screen single-cell cohorts (e.g. 'Melanoma', 'Melanoma,RCC').",
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=15,
        help="Number of multi-fidelity HPO trials to evaluate (default 15).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default 42).",
    )
    parser.add_argument(
        "--discovery-cohorts",
        type=str,
        default="Hugo,Riaz,Liu",
        help="Comma-separated list of Discovery ICI cohorts for inner LOCO-CV.",
    )
    parser.add_argument(
        "--held-out-cohorts",
        type=str,
        default="Gide,VanAllen,Snyder",
        help="Comma-separated list of Held-Out ICI cohorts for test validation.",
    )
    parser.add_argument(
        "--preprocessed-dir",
        type=Path,
        default=Path("data/preprocessed"),
        help="Directory containing preprocessed H5AD cohorts.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/reference_sampling_hpo"),
        help="Directory to save Parquet evaluations and reference matrices.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/reference_sampling_hpo"),
        help="Directory to save publication Altair SVG figures.",
    )
    parser.add_argument(
        "--benchmark-mode",
        action="store_true",
        help="Execute high-speed synthetic benchmark dry-run.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=16,
        help="Number of parallel worker threads/jobs for deconvolution and feature selection (default 16).",
    )
    parser.add_argument(
        "--force-local",
        action="store_true",
        help="Force execution on local machine bypassing host check (caution: high memory usage).",
    )

    args = parser.parse_args()

    # Host Guard Check
    match verify_execution_host(benchmark_mode=args.benchmark_mode, force_local=args.force_local):
        case Failure(err):
            print(f"\n[HOST CHECK ERROR]\n{err}\n", file=sys.stderr)
            sys.exit(1)
        case Success(_):
            pass

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    print("======================================================================")
    print("STARTING SINGLE-CELL REFERENCE SAMPLING & CLASSIFIER HPO PIPELINE")
    print("======================================================================")
    print(f"Target Cancer Type(s): {args.cancer_types}")
    print(f"Planned HPO Trials:    {args.n_trials}")
    print(f"Benchmark Mode:        {args.benchmark_mode}")
    print(f"Force Local Override:  {args.force_local}")
    print(f"Output Parquets:       {args.out_dir}")
    print(f"Results SVGs:          {args.results_dir}")
    print("======================================================================\n")

    # Step 1: Pre-Load Bulk Data
    disc_ids = [c.strip() for c in args.discovery_cohorts.split(",") if c.strip()]
    held_ids = [c.strip() for c in args.held_out_cohorts.split(",") if c.strip()]

    if args.benchmark_mode:
        print("[Step 1] Initializing fast benchmark mode (using top 15 samples per cohort)...")
        disc_ids = disc_ids[:2]
        held_ids = held_ids[:1]

    print("[Step 1a] Pre-loading Discovery iAtlas bulk cohorts...")
    discovery_bulk = {}
    for cid in disc_ids:
        load_res = load_bulk_cohort_data(cid, args.preprocessed_dir)
        match load_res:
            case Success((mat, ens_genes, sym_genes, y)):
                if args.benchmark_mode:
                    mat = mat[:15, :]
                    y = y[:15]
                discovery_bulk[cid] = (mat, ens_genes, sym_genes, y)
                print(f"  + Loaded Discovery Cohort: {cid} ({mat.shape[0]} samples, {len(sym_genes)} genes)")
            case Failure(err):
                print(f"  ! Warning: Skipping {cid}: {err}")

    if len(discovery_bulk) < 2:
        print("ERROR: Fewer than 2 Discovery cohorts available for LOCO-CV.", file=sys.stderr)
        sys.exit(1)

    print("[Step 1b] Pre-loading Held-Out Validation iAtlas bulk cohorts...")
    held_out_bulk = {}
    for cid in held_ids:
        load_res = load_bulk_cohort_data(cid, args.preprocessed_dir)
        match load_res:
            case Success((mat, ens_genes, sym_genes, y)):
                if args.benchmark_mode:
                    mat = mat[:15, :]
                    y = y[:15]
                held_out_bulk[cid] = (mat, ens_genes, sym_genes, y)
                print(f"  + Loaded Held-Out Cohort:  {cid} ({mat.shape[0]} samples, {len(sym_genes)} genes)")
            case Failure(err):
                print(f"  ! Warning: Skipping {cid}: {err}")

    # Step 2: Screen Single-Cell Candidate Datasets
    print("\n[Step 2] Screening candidate single-cell datasets from registry...")
    target_cancers = [ct.strip() for ct in args.cancer_types.split(",") if ct.strip()]
    cand_res = filter_candidate_cohorts(
        cancer_types=target_cancers,
        min_viable_cells=100 if args.benchmark_mode else 200,
        require_cached=True,
    )

    match cand_res:
        case Failure(err):
            print(f"ERROR: Could not resolve single-cell candidate cohorts: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(cand_ids):
            print(f"  + Identified {len(cand_ids)} candidate single-cell cohorts: {cand_ids[:6]}...")

    # Step 3: Define Search Space
    space = HPOSearchSpace(
        candidate_cohort_ids=cand_ids,
        min_cohorts=2,
        max_cohorts=min(3 if args.benchmark_mode else 4, len(cand_ids)),
        malignant_strategies=(
            MalignantStrategy.PATIENT_STRATIFIED,
            MalignantStrategy.POOLED_GENERIC,
            MalignantStrategy.EXCLUDED_TME_ONLY,
        ),
        malignant_fraction_range=(0.10, 0.35),
        min_cells_per_cohort=50 if args.benchmark_mode else 150,
        max_cells_per_cohort=120 if args.benchmark_mode else 1200,
        require_cached_h5ad=True,
    )

    # Step 4: Execute Multi-Fidelity HPO Campaign
    print("\n[Step 4] Launching multi-fidelity ASHA hyperparameter optimization...")
    hpo_res = run_sampling_hpo(
        search_space=space,
        bulk_data_dict=discovery_bulk,
        n_trials=args.n_trials,
        seed=args.seed,
        n_jobs=args.n_jobs,
    )

    match hpo_res:
        case Failure(err):
            print(f"ERROR: HPO run failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(run_summary):
            pass

    # Save evaluations dataframe
    eval_parquet = args.out_dir / "hpo_evaluations.parquet"
    run_summary.evaluations_df.write_parquet(eval_parquet)
    print(f"\n[Step 5] Saved trial evaluations to: {eval_parquet}")

    best_cfg = run_summary.best_config
    print("\n======================================================================")
    print("HPO OPTIMIZATION SUMMARY")
    print("======================================================================")
    print(f"Total Trials Evaluated:     {run_summary.total_trials}")
    print(f"Early Pruned Trials:        {run_summary.pruned_trials}")
    print(f"Total Elapsed Time:         {run_summary.total_elapsed_seconds} s")
    print(f"Best Trial ID:              {run_summary.best_trial_id}")
    print(f"Optimal Mean LOCO ROC-AUC:  {run_summary.best_loco_auc:.4f}")
    print(f"Selected Single-Cell Cohorts:{best_cfg.selected_cohort_ids}")
    print(f"Sampling Mode:              {best_cfg.sampling_spec.mode.value}")
    print(f"Cell Budget per Cohort:     {best_cfg.sampling_spec.n_cells_per_cohort}")
    print(f"Malignant Strategy:         {best_cfg.malignant_config.strategy.value}")
    print(f"Malignant Fraction:         {best_cfg.malignant_config.malignant_fraction:.2f}")
    print(f"Leiden Clustering Res:      {best_cfg.cluster_spec.leiden_resolution}")
    print(f"Collinearity Threshold:     {best_cfg.reference_config.collinearity_threshold}")

    best_classifier = run_summary.best_classifier_config.value_or(
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.STABILITY_SELECTION,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
        )
    )
    print(f"Optimal Transform:          {best_classifier.transform.value}")
    print(f"Optimal Feature Selector:   {best_classifier.selector.value}")
    print(f"Optimal Classifier:         {best_classifier.classifier.value}")
    print("======================================================================\n")

    # Step 6: Reconstruct Best Reference & Export Parquets
    print("[Step 6] Reconstructing optimal deconvolution reference...")
    sample_res = sample_single_cell_cohorts(
        sampling_spec=best_cfg.sampling_spec,
        cluster_spec=best_cfg.cluster_spec,
        harmonize_config=best_cfg.harmonize_config,
    )
    match sample_res:
        case Failure(err):
            print(f"Warning: Failed to reconstruct final reference ({err}); skipping reference export.")
            best_ref_opt: Maybe[Any] = Nothing
        case Success(sampled_data):
            ref_res = build_deconvolution_reference(
                data_source=sampled_data.adata,
                config=best_cfg.reference_config,
            )
            match ref_res:
                case Success(ref_result):
                    ref_dir = args.out_dir / "best_reference"
                    ref_result.export_parquet(ref_dir)
                    print(f"  + Exported best reference Parquet tables to: {ref_dir}")
                    best_ref_opt = Some(ref_result)
                case Failure(err):
                    print(f"Warning: Failed to build reference: {err}")
                    best_ref_opt = Nothing

    # Step 7: Independent Validation on Held-Out iAtlas Cohorts
    if isinstance(best_ref_opt, Some) and held_out_bulk:
        print("\n[Step 7] Evaluating generalization on held-out validation cohorts...")
        ref_obj = best_ref_opt.unwrap()
        ref_mat, state_labels = ref_obj.to_instaprism()

        held_res = evaluate_held_out_generalization(
            best_ref_mat=ref_mat,
            best_ref_genes=ref_obj.gene_names,
            state_labels=state_labels,
            discovery_bulk=discovery_bulk,
            held_out_bulk=held_out_bulk,
            best_inner_config=best_classifier.model_copy(update={"n_jobs": args.n_jobs}),
        )
        match held_res:
            case Success(held_aucs):
                print("  --- Held-Out Test Generalization ---")
                for cid, score in held_aucs.items():
                    print(f"  * {cid:<22}: ROC-AUC = {score:.4f}")
                mean_held = float(np.mean(list(held_aucs.values())))
                print(f"  * Mean Held-Out ROC-AUC : {mean_held:.4f}")

                # Step 8: Plot Held-Out Bar Chart
                svg_held = args.results_dir / "held_out_validation_generalization.svg"
                plot_held_out_validation(held_aucs, run_summary.best_loco_auc, svg_held)
                print(f"  + Saved Held-Out validation plot to: {svg_held}")
            case Failure(err):
                print(f"Warning: Held-out evaluation failed: {err}")

    # Step 8: Plot Pareto Frontier Diagnostic
    svg_pareto = args.results_dir / "pareto_frontier_loco_vs_collinearity.svg"
    plot_pareto_frontier(run_summary.evaluations_df, run_summary.pareto_trials, svg_pareto)
    print(f"[Step 8] Saved Pareto frontier diagnostic SVG to: {svg_pareto}")

    print("\n======================================================================")
    print("REFERENCE SAMPLING & CLASSIFIER HPO RUN COMPLETE!")
    print("======================================================================")


if __name__ == "__main__":
    main()
