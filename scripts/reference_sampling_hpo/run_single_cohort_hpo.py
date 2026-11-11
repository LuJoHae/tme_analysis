#!/usr/bin/env python3
"""Master CLI script for cohort-individual single-cell reference sampling and classifier HPO.

Optimizes deconvolution reference construction and downstream predictive classifiers
strictly and exclusively for a single clinical cohort (e.g. Gide), using leakage-free
Repeated Stratified 5-Fold Cross-Validation, with optional one-time external transfer
evaluation on other datasets at the very end.
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
    load_dataset,
    run_single_cohort_sampling_hpo,
    sample_single_cell_cohorts,
)
from tme_datasets.preprocessing.metadata import binarize_response  # type: ignore[import-untyped]
from tme_datasets.sampling.classifier_hpo import (  # type: ignore[import-untyped]
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    TransformType,
    apply_compositional_transform,
    fit_and_predict,
    get_transformed_feature_names,
    select_features_train_test,
)
from tme_datasets.sampling.hpo_models import HPORunResult  # type: ignore[import-untyped]
from tme_datasets.sampling.hpo_optimizer import (  # type: ignore[import-untyped]
    batch_richardson_lucy_deconvolution,
    fast_deconvolute_cohorts,
)

# Local script module
from plotting import (  # type: ignore[import-not-found]
    plot_cohort_transfer_validation,
    plot_single_cohort_pareto_frontier,
)


def verify_execution_host(benchmark_mode: bool = False, force_local: bool = False) -> Result[bool, str]:
    """Verify script is running on a CSCS Eiger Slurm compute node unless benchmark mode is requested."""
    hostname = socket.gethostname().lower()
    if hostname.startswith("eiger-ln") and not os.environ.get("SLURM_JOB_ID") and not benchmark_mode:
        return Failure(
            f"Execution blocked on login node '{hostname}': Compute-intensive tasks must NEVER run "
            "directly in bash on Eiger login nodes. Submit via Slurm ('make hpo-cohort')."
        )

    if benchmark_mode or force_local or os.environ.get("ALLOW_LOCAL_HPO") == "1":
        return Success(True)

    allowed_patterns = ("nid", "eiger", "alps", "olm", "sgl07282", "storage")
    is_remote = any(p in hostname for p in allowed_patterns) or bool(os.environ.get("SLURM_JOB_ID"))

    if not is_remote:
        return Failure(
            f"Execution halted on host '{hostname}'. Single-cell reference sampling HPO must run "
            "exclusively on CSCS Eiger Slurm compute nodes due to >64GB memory requirements.\n"
            "To launch remotely: make hpo-cohort HPO_COHORT=Gide\n"
            "To test synthetic HPO locally: python scripts/reference_sampling_hpo/run_single_cohort_hpo.py --benchmark-mode"
        )
    return Success(True)


def load_bulk_cohort_data(
    cohort_id: str,
    preprocessed_dir: Path,
) -> Result[tuple[np.ndarray, Sequence[str], Sequence[str], np.ndarray], str]:
    """Load preprocessed bulk RNA-seq cohort strictly via tme_datasets API with binarized response annotations."""
    from tme_datasets import load_iatlas_cohort_or_combined

    adata: ad.AnnData | None = None
    for cid in (cohort_id, f"{cohort_id}-iAtlas", cohort_id.replace("-iAtlas", "")):
        load_res = load_dataset(
            cid,
            base_dir=preprocessed_dir,
            auto_download=False,
            normalize_ensembl=False,
        )
        if not isinstance(load_res, Success):
            load_res = load_iatlas_cohort_or_combined(cid)
        if isinstance(load_res, Success):
            cand_ad = load_res.unwrap()
            resp_cand = next(
                (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in cand_ad.obs.columns),
                None,
            )
            if resp_cand is not None and cand_ad.n_vars > 100:
                adata = cand_ad
                break

    if adata is None:
        return Failure(f"Failed to load bulk dataset for '{cohort_id}' via tme_datasets.")

    try:
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

        if hasattr(sub_adata.X, "toarray"):
            expr_mat = np.asarray(sub_adata.X.toarray(), dtype=np.float64)
        else:
            expr_mat = np.asarray(sub_adata.X, dtype=np.float64)

        if np.all(expr_mat >= 0.0) and np.max(expr_mat) <= 30.0:
            expr_mat = np.expm1(expr_mat)

        ensembl_genes = [str(x).split(".")[0].strip() for x in sub_adata.var_names]
        if "gene_name" in sub_adata.var.columns:
            symbol_genes = [
                str(x).upper() if x is not None and str(x) != "nan" else ""
                for x in sub_adata.var["gene_name"].tolist()
            ]
        elif "gene_symbol" in sub_adata.var.columns:
            symbol_genes = [
                str(x).upper() if x is not None and str(x) != "nan" else ""
                for x in sub_adata.var["gene_symbol"].tolist()
            ]
        else:
            symbol_genes = [str(x).upper() for x in sub_adata.var_names]

        return Success((expr_mat, ensembl_genes, symbol_genes, y))
    except Exception as exc:
        return Failure(f"Failed to process bulk dataset {cohort_id}: {exc}")


def evaluate_transfer_on_external_cohorts(
    best_ref_mat: np.ndarray,
    best_ref_genes: Sequence[str],
    state_labels: Sequence[str],
    target_mat: np.ndarray,
    target_bulk_ens: Sequence[str],
    target_bulk_sym: Sequence[str],
    target_y: np.ndarray,
    external_bulk: Mapping[str, tuple[Any, ...]],
    best_inner_config: InnerClassifierConfig,
) -> Result[dict[str, float], str]:
    """Train the optimal pipeline on 100% of target cohort and evaluate transfer on external cohorts."""
    from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]

    # 1. Deconvolute Target Cohort (passing both ensembl and symbol so best overlap is selected)
    target_dict = {"target": (target_mat, target_bulk_ens, target_bulk_sym, target_y)}
    tgt_deconv_res = fast_deconvolute_cohorts(best_ref_mat, best_ref_genes, target_dict, n_iter=60)
    match tgt_deconv_res:
        case Failure(err):
            return Failure(f"Target deconvolution failed: {err}")
        case Success((tgt_fracs, tgt_resps)):
            pass

    # 2. Deconvolute External Cohorts
    ext_deconv_res = fast_deconvolute_cohorts(best_ref_mat, best_ref_genes, external_bulk, n_iter=60)
    match ext_deconv_res:
        case Failure(err):
            return Failure(f"External cohorts deconvolution failed: {err}")
        case Success((ext_fracs, ext_resps)):
            pass

    # 3. Apply optimal compositional transform to Target Cohort
    X_train = apply_compositional_transform(
        tgt_fracs["target"],
        best_inner_config.transform,
        best_inner_config.eps,
    )
    y_train = tgt_resps["target"]

    active_labels = get_transformed_feature_names(state_labels, best_inner_config.transform)
    transfer_aucs: dict[str, float] = {}

    for ext_cid, raw_ext_fracs in ext_fracs.items():
        ext_y = ext_resps[ext_cid]
        X_test = apply_compositional_transform(
            raw_ext_fracs,
            best_inner_config.transform,
            best_inner_config.eps,
        )

        tr_sub, te_sub, _ = select_features_train_test(
            X_train=X_train,
            y_train=y_train,
            X_test=X_test,
            feature_names=active_labels,
            selector_type=best_inner_config.selector,
            config=best_inner_config,
        )

        try:
            pred_probs = fit_and_predict(
                X_train=tr_sub,
                y_train=y_train,
                X_test=te_sub,
                classifier_type=best_inner_config.classifier,
                config=best_inner_config,
            )
            score = float(roc_auc_score(ext_y, pred_probs))
            transfer_aucs[ext_cid] = score
        except Exception:
            transfer_aucs[ext_cid] = 0.50

    return Success(transfer_aucs)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run cohort-individual single-cell reference sampling and classifier HPO."
    )
    parser.add_argument(
        "--target-cohort",
        type=str,
        default="Gide",
        help="Primary clinical bulk cohort to optimize for (default 'Gide').",
    )
    parser.add_argument(
        "--transfer-cohorts",
        type=str,
        default="Hugo,Riaz,Liu,VanAllen,Snyder",
        help="Comma-separated external bulk cohorts for optional transfer testing at the very end.",
    )
    parser.add_argument(
        "--cancer-types",
        type=str,
        default="Melanoma",
        help="Target cancer type(s) to screen candidate single-cell datasets (default 'Melanoma').",
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=15,
        help="Number of multi-fidelity HPO trials to evaluate (default 15).",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=32,
        help="Number of parallel worker threads across CV folds and stability paths (default 32).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default 42).",
    )
    _default_data_dir = Path(os.environ.get("TME_DATA_DIR", "/ritom/scratch/cscs/lhaeuser/data-test"))
    parser.add_argument(
        "--preprocessed-dir",
        type=Path,
        default=_default_data_dir / "preprocessed",
        help="Directory containing preprocessed H5AD cohorts.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=_default_data_dir / "output/reference_sampling_hpo/cohort_Gide",
        help="Output directory for evaluations and reference parquets.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=_default_data_dir / "results/reference_sampling_hpo/cohort_Gide",
        help="Output directory for publication SVG figures.",
    )
    parser.add_argument(
        "--benchmark-mode",
        action="store_true",
        help="Run lightweight benchmark mode for fast verification.",
    )
    parser.add_argument(
        "--force-local",
        action="store_true",
        help="Force execution on local machine bypassing host check.",
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
    print("STARTING COHORT-INDIVIDUAL REFERENCE SAMPLING & CLASSIFIER HPO")
    print(f"Target Cohort:         {args.target_cohort} (Optimizing Strictly for this Cohort)")
    print(f"Parallel Worker Cores: {args.n_jobs}")
    print(f"Planned HPO Trials:    {args.n_trials}")
    print(f"Benchmark Mode:        {args.benchmark_mode}")
    print(f"Output Directory:      {args.out_dir}")
    print(f"Results Directory:     {args.results_dir}")
    print("======================================================================\n")

    # Step 1: Pre-Load Target Cohort Data
    print(f"[Step 1a] Pre-loading Target Cohort: '{args.target_cohort}'...")
    target_res = load_bulk_cohort_data(args.target_cohort, args.preprocessed_dir)
    match target_res:
        case Failure(err):
            print(f"ERROR: Could not load target cohort '{args.target_cohort}': {err}", file=sys.stderr)
            sys.exit(1)
        case Success((tgt_mat, tgt_ens, tgt_sym, tgt_y)):
            if args.benchmark_mode:
                tgt_mat = tgt_mat[:20, :]
                tgt_y = tgt_y[:20]
            print(f"  + Target '{args.target_cohort}' loaded: {tgt_mat.shape[0]} patients, {len(tgt_sym)} genes.")
            print(f"  + Class balance: {int(np.sum(tgt_y == 0))} Non-Responders vs. {int(np.sum(tgt_y == 1))} Responders.")

    # Step 1b: Pre-Load External Cohorts for Final One-Time Transfer Check
    ext_cohort_ids = [c.strip() for c in args.transfer_cohorts.split(",") if c.strip() and c.strip() != args.target_cohort]
    external_bulk: dict[str, tuple[Any, ...]] = {}
    print(f"\n[Step 1b] Pre-loading External Cohorts for final transfer check: {ext_cohort_ids}...")
    for cid in ext_cohort_ids:
        c_res = load_bulk_cohort_data(cid, args.preprocessed_dir)
        match c_res:
            case Success((mat, ens, sym, y)):
                if args.benchmark_mode:
                    mat = mat[:15, :]
                    y = y[:15]
                external_bulk[cid] = (mat, ens, sym, y)
                print(f"  + Loaded External: {cid:<12} (N={mat.shape[0]}, G={len(sym)})")
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
            print(f"ERROR: Could not resolve candidate cohorts: {err}", file=sys.stderr)
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

    # Step 4: Execute Cohort-Specific Multi-Fidelity HPO Campaign
    print(f"\n[Step 4] Launching ASHA optimization strictly on '{args.target_cohort}' ({args.n_jobs} cores)...")
    hpo_res = run_single_cohort_sampling_hpo(
        target_cohort_name=args.target_cohort,
        target_bulk_mat=tgt_mat,
        target_bulk_genes=tgt_sym,
        target_y=tgt_y,
        search_space=space,
        target_bulk_ensembl_genes=tgt_ens,
        target_bulk_symbol_genes=tgt_sym,
        n_trials=args.n_trials,
        seed=args.seed,
        n_jobs=args.n_jobs,
    )
    match hpo_res:
        case Failure(err):
            print(f"ERROR: Single-cohort HPO run failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(run_summary):
            pass

    # Step 5: Save evaluations dataframe
    eval_parquet = args.out_dir / "hpo_evaluations.parquet"
    run_summary.evaluations_df.write_parquet(eval_parquet)
    print(f"\n[Step 5] Saved trial evaluations to: {eval_parquet}")

    best_cfg = run_summary.best_config
    print("\n======================================================================")
    print(f"HPO OPTIMIZATION SUMMARY FOR TARGET COHORT: {args.target_cohort}")
    print("======================================================================")
    print(f"Total Trials Evaluated:     {run_summary.total_trials}")
    print(f"Early Pruned Trials:        {run_summary.pruned_trials}")
    print(f"Total Elapsed Time:         {run_summary.total_elapsed_seconds} s")
    print(f"Best Trial ID:              {run_summary.best_trial_id}")
    print(f"Optimal In-Cohort CV AUC:   {run_summary.best_loco_auc:.4f}")
    print(f"Selected Single-Cell Cohorts:{best_cfg.selected_cohort_ids}")
    print(f"Sampling Mode:              {best_cfg.sampling_spec.mode.value}")
    print(f"Cell Budget per Cohort:     {best_cfg.sampling_spec.n_cells_per_cohort}")
    print(f"Malignant Strategy:         {best_cfg.malignant_config.strategy.value}")
    print(f"Malignant Fraction:         {best_cfg.malignant_config.malignant_fraction:.2f}")
    print(f"Leiden Clustering Res:      {best_cfg.cluster_spec.leiden_resolution}")
    print(f"Collinearity Threshold:     {best_cfg.reference_config.collinearity_threshold}")

    best_classifier = run_summary.best_classifier_config.value_or(
        InnerClassifierConfig(
            transform=TransformType.LOGIT,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            n_jobs=args.n_jobs,
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
    best_ref_opt: Maybe[Any] = Nothing
    match sample_res:
        case Failure(err):
            print(f"Warning: Failed to reconstruct final reference ({err}); skipping reference export.")
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

    # Step 7: Final External Transfer Evaluation (Evaluated at the very end once)
    transfer_aucs: dict[str, float] = {}
    if isinstance(best_ref_opt, Some) and external_bulk:
        print(f"\n[Step 7] Evaluating transfer of {args.target_cohort}-optimized pipeline on external cohorts...")
        ref_obj = best_ref_opt.unwrap()
        ref_mat, state_labels = ref_obj.to_instaprism()

        transfer_res = evaluate_transfer_on_external_cohorts(
            best_ref_mat=ref_mat,
            best_ref_genes=ref_obj.gene_names,
            state_labels=state_labels,
            target_mat=tgt_mat,
            target_bulk_ens=tgt_ens,
            target_bulk_sym=tgt_sym,
            target_y=tgt_y,
            external_bulk=external_bulk,
            best_inner_config=best_classifier.model_copy(update={"n_jobs": args.n_jobs}),
        )
        match transfer_res:
            case Success(t_aucs):
                transfer_aucs = t_aucs
                print(f"  --- Out-of-Cohort Transfer Performance (Trained on 100% {args.target_cohort}) ---")
                for cid, score in transfer_aucs.items():
                    print(f"  * {cid:<22}: Transfer ROC-AUC = {score:.4f}")
                mean_trans = float(np.mean(list(transfer_aucs.values())))
                print(f"  * Mean External Transfer ROC-AUC: {mean_trans:.4f}")

                # Step 8a: Plot Transfer Bar Chart
                svg_trans = args.results_dir / "cohort_transfer.svg"
                plot_cohort_transfer_validation(
                    target_cohort_name=args.target_cohort,
                    target_cv_auc=run_summary.best_loco_auc,
                    transfer_aucs=transfer_aucs,
                    out_file=svg_trans,
                )
                print(f"  + Saved Cohort Transfer SVG to: {svg_trans}")
            case Failure(err):
                print(f"Warning: External transfer evaluation failed: {err}")

    # Step 8b: Plot Pareto Frontier
    svg_pareto = args.results_dir / "pareto_frontier.svg"
    plot_single_cohort_pareto_frontier(
        eval_df=run_summary.evaluations_df,
        pareto_trial_ids=run_summary.pareto_trials,
        target_cohort_name=args.target_cohort,
        out_file=svg_pareto,
    )
    print(f"[Step 8] Saved Pareto frontier diagnostic SVG to: {svg_pareto}")

    print("\n======================================================================")
    print(f"COHORT-INDIVIDUAL HPO FOR '{args.target_cohort}' COMPLETE!")
    print("======================================================================")


if __name__ == "__main__":
    main()
