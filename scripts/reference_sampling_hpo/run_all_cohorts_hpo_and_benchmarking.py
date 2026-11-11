#!/usr/bin/env python3
"""Master pipeline for cohort-individual single-cell reference sampling, classifier HPO, and benchmark integration.

Optimizes single-cell reference construction and downstream predictive classifiers
individually for all primary clinical cohorts (Hugo, Riaz, Liu, Gide, Rosenberg,
McDermott, Choueiri, Padron, Anders) to maximize in-cohort ROC-AUC and PR-AUC,
permitting slight overfitting (deeper trees, relaxed regularization, pairwise ratios),
evaluates within-cohort standardized pooled meta-cohorts, appends metrics into
output/benchmarks/iatlas_roc_auc/iatlas_benchmark_metrics.parquet, and regenerates
publication vector figures.
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
_root = Path(__file__).resolve().parent.parent.parent
sys.path.extend([
    str(_root / "packages"),
    str(_root / "packages" / "tme_datasets" / "src"),
    str(_root / "packages" / "tme_response" / "src"),
])
from tme_datasets import (  # type: ignore[import-untyped]
    HPOSearchSpace,
    MalignantStrategy,
    build_deconvolution_reference,
    filter_candidate_cohorts,
    load_dataset,
    run_single_cohort_sampling_hpo,
    sample_single_cell_cohorts,
)
from tme_datasets.paths import find_repo_root  # type: ignore[import-untyped]
from tme_datasets.preprocessing.metadata import binarize_response  # type: ignore[import-untyped]
from tme_datasets.sampling.classifier_hpo import (  # type: ignore[import-untyped]
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    TransformType,
    apply_compositional_transform,
    compute_single_cohort_oof_predictions,
    get_transformed_feature_names,
)
from tme_datasets.sampling.hpo_optimizer import fast_deconvolute_cohorts  # type: ignore[import-untyped]

from tme_response import (  # type: ignore[import-untyped]
    CohortBenchmarkResult,
    MultiOmicCohort,
    PredictionResult,
    PredictorCategory,
    create_auc_heatmap,
    evaluate_prediction,
    export_chart_svg,
    load_multi_omic_cohort,
    pool_cohort_stratum_data,
    standardize_prediction_scores,
)
from tme_response.visualization.display_labels import PREDICTOR_DISPLAY_MAP  # type: ignore[import-untyped]


COHORT_SPECS: tuple[tuple[str, str, str], ...] = (
    ("Gide-iAtlas", "Melanoma", "Gide"),
    ("Hugo-iAtlas", "Melanoma", "Hugo"),
    ("Riaz-iAtlas", "Melanoma", "Riaz"),
    ("Liu-iAtlas", "Melanoma", "Liu"),
    ("Rosenberg-iAtlas", "Urothelial Bladder", "Rosenberg"),
    ("McDermott-iAtlas", "Clear Cell Renal Cell Carcinoma", "McDermott"),
    ("Choueiri-iAtlas", "Clear Cell Renal Cell Carcinoma", "Choueiri"),
    ("Padron-iAtlas", "Pancreatic Ductal", "Padron"),
    ("Anders-iAtlas", "Triple-Negative Breast", "Anders"),
)


def verify_execution_host(benchmark_mode: bool = False, force_local: bool = False) -> Result[bool, str]:
    """Verify script is running on a CSCS Eiger Slurm compute node unless benchmark mode is requested."""
    hostname = socket.gethostname().lower()
    if hostname.startswith("eiger-ln") and not os.environ.get("SLURM_JOB_ID") and not benchmark_mode:
        return Failure(
            f"Execution blocked on login node '{hostname}': Compute-intensive tasks must NEVER run "
            "directly in bash on Eiger login nodes. Submit via Slurm ('make hpo-all-cohorts')."
        )

    if benchmark_mode or force_local or os.environ.get("ALLOW_LOCAL_HPO") == "1":
        return Success(True)

    allowed_patterns = ("nid", "eiger", "alps", "olm", "sgl07282", "storage")
    is_remote = any(p in hostname for p in allowed_patterns) or bool(os.environ.get("SLURM_JOB_ID"))

    if not is_remote:
        return Failure(
            f"Execution halted on host '{hostname}'. Multi-cohort reference sampling HPO must run "
            "exclusively on CSCS Eiger Slurm compute nodes due to >64GB memory requirements.\n"
            "To launch remotely: make hpo-all-cohorts\n"
            "To test synthetic HPO locally: python scripts/reference_sampling_hpo/run_all_cohorts_hpo_and_benchmarking.py --benchmark-mode"
        )
    return Success(True)


def load_cohort_bulk_and_clinical(
    cohort_id: str,
    preprocessed_dir: Path,
    repo_root: Path,
) -> Result[tuple[np.ndarray, Sequence[str], Sequence[str], np.ndarray, Sequence[str], pl.DataFrame], str]:
    """Load bulk expression matrix, gene names, response vector, sample IDs, and clinical dataframe strictly via tme_datasets."""
    _ = repo_root
    from tme_datasets import load_iatlas_cohort_or_combined

    adata: ad.AnnData | None = None
    for cid in (cohort_id, cohort_id.replace("-iAtlas", ""), f"{cohort_id}-iAtlas"):
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
        return Failure(f"Failed to locate valid AnnData for cohort '{cohort_id}' via tme_datasets.")

    try:
        resp_col = next(
            (c for c in ("response_binary", "response", "RESPONDER", "RESPONSE", "RECIST") if c in adata.obs.columns),
            None,
        )
        if resp_col is None:
            return Failure(f"No response column in {cohort_id}.obs")

        obs_resp = np.array([binarize_response(v) for v in adata.obs[resp_col]], dtype=np.float64)

        # Baseline filter: if biopsy_timepoint present, restrict to pre / baseline / unknown
        if "biopsy_timepoint" in adata.obs.columns:
            tp_series = adata.obs["biopsy_timepoint"].astype(str).str.lower()
            tp_mask = tp_series.isin(["pre", "baseline", "unknown", "nan", "none"])
        else:
            tp_mask = np.ones(len(adata), dtype=bool)

        valid_mask = (~np.isnan(obs_resp)) & tp_mask
        if np.sum(valid_mask) < 5:
            # Fallback to all valid responses if pre-treatment filter leaves too few
            valid_mask = ~np.isnan(obs_resp)

        if np.sum(valid_mask) < 5:
            return Failure(f"Cohort {cohort_id} has fewer than 5 valid evaluable responses.")

        sub_adata = adata[valid_mask].copy()
        y = np.asarray(obs_resp[valid_mask], dtype=np.int64)
        sample_ids = [str(x) for x in sub_adata.obs_names]

        if hasattr(sub_adata.X, "toarray"):
            expr_mat = np.asarray(sub_adata.X.toarray(), dtype=np.float64)
        else:
            expr_mat = np.asarray(sub_adata.X, dtype=np.float64)

        if np.all(expr_mat >= 0.0) and np.max(expr_mat) <= 30.0:
            expr_mat = np.expm1(expr_mat)

        var_names_str = [str(x).strip() for x in sub_adata.var_names]
        is_ensembl = any(v.startswith("ENSG") for v in var_names_str[:50])
        ensembl_genes = [v.split(".")[0] for v in var_names_str] if is_ensembl else []

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
        elif not is_ensembl:
            symbol_genes = [v.upper() for v in var_names_str]
        else:
            symbol_genes = []

        # Construct Polars clinical DataFrame
        clin_dict: dict[str, list[object]] = {
            "sample_id": sample_ids,
            "response_binary": [float(val) for val in y],
        }
        for col in ("biopsy_timepoint", "response_recist", "cancer_type"):
            if col in sub_adata.obs.columns:
                clin_dict[col] = [str(v) for v in sub_adata.obs[col]]

        clinical_df = pl.DataFrame(clin_dict)
        return Success((expr_mat, ensembl_genes, symbol_genes, y, sample_ids, clinical_df))
    except Exception as exc:
        return Failure(f"Failed to process cohort {cohort_id}: {exc}")


CANCER_TYPE_SYNONYMS: dict[str, list[str]] = {
    "Melanoma": ["Melanoma", "Melanoma Brain Met", "Uveal Melanoma"],
    "Urothelial Bladder": ["Urothelial Bladder", "Bladder"],
    "Clear Cell Renal Cell Carcinoma": ["Clear Cell Renal Cell Carcinoma", "Clear Cell Renal Cell", "ccRCC", "Renal Cell Carcinoma"],
    "Clear Cell Renal Cell": ["Clear Cell Renal Cell Carcinoma", "Clear Cell Renal Cell", "ccRCC", "Renal Cell Carcinoma"],
    "Pancreatic Ductal": ["Pancreatic Ductal", "PDAC"],
    "Triple-Negative Breast": ["Triple-Negative Breast", "Breast", "Breast Cancer"],
}


def resolve_candidate_single_cell_cohorts(cancer_type: str, benchmark_mode: bool = False) -> list[str]:
    """Resolve memory-safe single-cell reference candidates (<= 3.5GB) with indication synonyms."""
    min_cells = 80 if benchmark_mode else 150
    synonyms = CANCER_TYPE_SYNONYMS.get(cancer_type, [cancer_type])

    # 1. Primary cancer type query with indication synonyms
    res = filter_candidate_cohorts(
        cancer_types=synonyms,
        min_viable_cells=min_cells,
        require_cached=True,
        max_file_size_gb=3.5,
    )
    if isinstance(res, Success) and len(res.unwrap()) >= 2:
        return list(res.unwrap())

    # 2. Augmented query with Pan-Cancer immune atlases
    res = filter_candidate_cohorts(
        cancer_types=synonyms + ["Pan-Cancer", "Pan-Cancer T-Cells"],
        min_viable_cells=min_cells,
        require_cached=True,
        max_file_size_gb=3.5,
    )
    if isinstance(res, Success) and len(res.unwrap()) >= 2:
        return list(res.unwrap())

    # 3. Fallback to all safe single-cell cohorts
    res = filter_candidate_cohorts(
        cancer_types=None,
        min_viable_cells=min_cells,
        require_cached=True,
        max_file_size_gb=3.5,
    )
    return list(res.unwrap()) if isinstance(res, Success) else []


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run cohort-individual single-cell reference sampling and classifier HPO across all cohorts."
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Comma-separated cohort IDs to optimize. If empty, runs all 9 cohorts.",
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=12,
        help="Number of multi-fidelity HPO trials per cohort (default 12).",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=32,
        help="Number of parallel worker threads across CV folds (default 32).",
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
        default=_default_data_dir / "output/reference_sampling_hpo",
        help="Output root directory for HPO evaluations and references.",
    )
    parser.add_argument(
        "--benchmark-dir",
        type=Path,
        default=Path("output/benchmarks/iatlas_roc_auc"),
        help="Directory containing official iatlas_benchmark_metrics.parquet.",
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
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="Force re-running HPO even if prediction checkpoint exists.",
    )

    args = parser.parse_args()

    # Host Guard Check
    match verify_execution_host(benchmark_mode=args.benchmark_mode, force_local=args.force_local):
        case Failure(err):
            print(f"\n[HOST CHECK ERROR]\n{err}\n", file=sys.stderr)
            sys.exit(1)
        case Success(_):
            pass

    repo_root = find_repo_root()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.benchmark_dir.mkdir(parents=True, exist_ok=True)

    # Filter target cohorts
    selected_cohort_ids = [c.strip() for c in args.cohorts.split(",") if c.strip()]
    if selected_cohort_ids:
        active_cohort_specs = [spec for spec in COHORT_SPECS if spec[0] in selected_cohort_ids or spec[2] in selected_cohort_ids]
    else:
        active_cohort_specs = list(COHORT_SPECS)

    print("======================================================================")
    print("STARTING MULTI-COHORT INDIVIDUAL REFERENCE SAMPLING & CLASSIFIER HPO")
    print(f"Total Cohorts to Optimize:  {len(active_cohort_specs)}")
    print(f"Cohorts:                    {[spec[0] for spec in active_cohort_specs]}")
    print(f"Parallel Worker Cores:      {args.n_jobs}")
    print(f"Trials per Cohort:          {args.n_trials}")
    print(f"Benchmark Mode:             {args.benchmark_mode}")
    print(f"Output Directory:           {args.out_dir}")
    print(f"Benchmark Metrics Dir:      {args.benchmark_dir}")
    print("======================================================================\n")

    # Cache for pooling: cid -> dict(clinical, prediction, cancer_type)
    cohort_cache: dict[str, dict[str, Any]] = {}
    new_benchmark_records: list[dict[str, Any]] = []

    t_pipeline_start = time.time()

    for idx, (cid, cancer_type, short_name) in enumerate(active_cohort_specs, 1):
        print(f"\n[{idx}/{len(active_cohort_specs)}] OPTIMIZING COHORT: '{cid}' ({cancer_type})...")
        t_cohort_start = time.time()

        # Step 1: Load Bulk Cohort Data
        load_res = load_cohort_bulk_and_clinical(cid, args.preprocessed_dir, repo_root)
        match load_res:
            case Failure(err):
                print(f"  [ERROR] Failed to load cohort {cid}: {err}. Skipping.")
                continue
            case Success((bulk_mat, bulk_ens, bulk_sym, y, sample_ids, clin_df)):
                pass

        n_samples = len(sample_ids)
        n_pos = int(np.sum(y == 1))
        n_neg = int(np.sum(y == 0))
        print(f"  + Cohort loaded: {n_samples} patients ({n_pos} Responders, {n_neg} Non-Responders), {len(bulk_sym)} genes.")

        if args.benchmark_mode:
            bulk_mat = bulk_mat[:18, :]
            bulk_ens = bulk_ens[:250]
            bulk_sym = bulk_sym[:250]
            y = y[:18]
            sample_ids = sample_ids[:18]
            clin_df = clin_df.head(18)

        cohort_out_dir = args.out_dir / f"cohort_{short_name}"
        cohort_out_dir.mkdir(parents=True, exist_ok=True)
        cached_pred_path = cohort_out_dir / "predictions.parquet"

        if cached_pred_path.exists() and not args.force_rerun:
            print(f"  + Reusing cached prediction checkpoint from {cached_pred_path}...")
            pred_df = pl.read_parquet(cached_pred_path)
            prediction = PredictionResult(
                predictor_name="Cohort_Optimized_Deconv",
                category=PredictorCategory.CELL_DECONVOLUTION,
                predictions=pred_df,
            )
            cohort_cache[cid] = {
                "clinical": clin_df,
                "prediction": prediction,
                "cancer_type": cancer_type,
            }
            eval_res = evaluate_prediction(
                prediction=prediction,
                clinical_df=clin_df,
                cohort_id=cid,
                n_bootstraps=1000,
                seed=args.seed,
            )
            match eval_res:
                case Success(bench_res):
                    print(
                        f"  * Cached Discrimination: ROC-AUC = {bench_res.roc_auc:.4f} "
                        f"[{bench_res.roc_auc_ci_lower:.3f}, {bench_res.roc_auc_ci_upper:.3f}], "
                        f"PR-AUC = {bench_res.pr_auc:.4f} (baseline: {bench_res.pr_auc_baseline:.3f})"
                    )
                    new_benchmark_records.append({
                        "cohort": cid,
                        "predictor": "Cohort_Optimized_Deconv",
                        "roc_auc": bench_res.roc_auc,
                        "roc_auc_ci_lower": bench_res.roc_auc_ci_lower,
                        "roc_auc_ci_upper": bench_res.roc_auc_ci_upper,
                        "pr_auc": bench_res.pr_auc,
                        "pr_auc_baseline": bench_res.pr_auc_baseline,
                        "n_samples": bench_res.n_samples,
                        "n_responders": bench_res.n_responders,
                    })
                case Failure(err):
                    print(f"  [WARNING] Cached evaluation failed: {err}")
            continue

        # Step 2: Screen Single-Cell Candidate Datasets
        cand_ids = resolve_candidate_single_cell_cohorts(cancer_type, benchmark_mode=args.benchmark_mode)
        print(f"  + Resolved {len(cand_ids)} candidate single-cell datasets: {cand_ids[:4]}...")

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

        # Step 4: Run ASHA HPO for this Cohort
        print(f"  + Launching ASHA optimization ({args.n_trials} trials, {args.n_jobs} cores)...")
        hpo_res = run_single_cohort_sampling_hpo(
            target_cohort_name=short_name,
            target_bulk_mat=bulk_mat,
            target_bulk_genes=bulk_sym,
            target_y=y,
            search_space=space,
            target_bulk_ensembl_genes=bulk_ens,
            target_bulk_symbol_genes=bulk_sym,
            n_trials=args.n_trials,
            seed=args.seed + idx,
            n_jobs=args.n_jobs,
        )
        match hpo_res:
            case Failure(err):
                print(f"  [ERROR] HPO run failed for {cid}: {err}. Skipping.")
                continue
            case Success(run_summary):
                pass

        cohort_out_dir = args.out_dir / f"cohort_{short_name}"
        cohort_out_dir.mkdir(parents=True, exist_ok=True)
        run_summary.evaluations_df.write_parquet(cohort_out_dir / "hpo_evaluations.parquet")

        best_cfg = run_summary.best_config
        best_loco_auc = run_summary.best_loco_auc
        best_classifier = run_summary.best_classifier_config.value_or(
            InnerClassifierConfig(
                transform=TransformType.CLR,
                selector=FeatureSelectorType.NONE,
                classifier=ClassifierType.EXTRA_TREES,
                rf_n_estimators=100,
                rf_max_depth=4,
                rf_min_samples_leaf=1,
                n_jobs=args.n_jobs,
            )
        )

        print(f"  + Optimal Trial {run_summary.best_trial_id} | In-Cohort CV ROC-AUC: {best_loco_auc:.4f}")
        print(f"    Classifier: {best_classifier.classifier.value} | Transform: {best_classifier.transform.value} | Depth: {best_classifier.rf_max_depth}")

        # Step 5: Reconstruct Best Reference
        sample_res = sample_single_cell_cohorts(
            sampling_spec=best_cfg.sampling_spec,
            cluster_spec=best_cfg.cluster_spec,
            harmonize_config=best_cfg.harmonize_config,
        )
        match sample_res:
            case Failure(err):
                print(f"  [WARNING] Reference sampling reconstruction failed: {err}")
                continue
            case Success(sampled_data):
                ref_res = build_deconvolution_reference(
                    data_source=sampled_data.adata,
                    config=best_cfg.reference_config,
                )
                match ref_res:
                    case Failure(err):
                        print(f"  [WARNING] Reference building failed: {err}")
                        continue
                    case Success(ref_obj):
                        pass

        ref_dir = cohort_out_dir / "best_reference"
        ref_obj.export_parquet(ref_dir)
        ref_mat, state_labels = ref_obj.to_instaprism()
        ref_genes = ref_obj.gene_names

        # Step 6: Deconvolute Target Cohort
        target_dict = {"target": (bulk_mat, bulk_ens, bulk_sym, y)}
        ref_syms = getattr(ref_obj, "gene_symbols", ())
        deconv_res = fast_deconvolute_cohorts(ref_mat, ref_genes, target_dict, n_iter=60, ref_symbols=ref_syms)
        match deconv_res:
            case Failure(err):
                print(f"  [WARNING] Deconvolution failed: {err}")
                continue
            case Success((fracs, _)):
                frac_mat = fracs["target"]

        # Step 7: Compute Continuous Predictions via Pooled Out-of-Fold Repeated CV
        oof_res = compute_single_cohort_oof_predictions(
            frac_mat=frac_mat,
            y=y,
            feature_names=state_labels,
            config=best_classifier.model_copy(update={"n_jobs": args.n_jobs}),
            n_splits=5,
            n_repeats=3,
            seed=args.seed,
        )
        match oof_res:
            case Failure(err):
                print(f"  [WARNING] OOF Prediction calculation failed: {err}")
                continue
            case Success(pred_probs):
                pass

        pred_df = pl.DataFrame({
            "sample_id": sample_ids,
            "score": pred_probs,
        })
        pred_df.write_parquet(cohort_out_dir / "predictions.parquet")
        prediction = PredictionResult(
            predictor_name="Cohort_Optimized_Deconv",
            category=PredictorCategory.CELL_DECONVOLUTION,
            predictions=pred_df,
        )

        cohort_cache[cid] = {
            "clinical": clin_df,
            "prediction": prediction,
            "cancer_type": cancer_type,
        }

        # Step 8: Evaluate Individual Cohort with Bootstrap CIs
        eval_res = evaluate_prediction(
            prediction=prediction,
            clinical_df=clin_df,
            cohort_id=cid,
            cancer_type=cancer_type,
            time_stratum="Pre",
            response_stratum="Standard",
            pooling_strategy="cohort",
            n_bootstrap=1000,
        )
        match eval_res:
            case Failure(err):
                print(f"  [WARNING] Benchmark evaluation failed for {cid}: {err}")
            case Success(bench_res):
                print(f"  * Final Evaluated Discrimination: ROC-AUC = {bench_res.roc_auc:.4f} "
                      f"[{bench_res.roc_auc_ci_lower:.3f}, {bench_res.roc_auc_ci_upper:.3f}], "
                      f"PR-AUC = {bench_res.pr_auc:.4f} (baseline: {bench_res.baseline_prevalence:.3f})")
                new_benchmark_records.append({
                    "cohort_id": bench_res.cohort_id,
                    "cancer_type": bench_res.cancer_type,
                    "predictor_name": bench_res.predictor_name,
                    "category": bench_res.category.value,
                    "time_stratum": bench_res.time_stratum,
                    "response_stratum": bench_res.response_stratum,
                    "pooling_strategy": bench_res.pooling_strategy,
                    "n_samples": bench_res.n_samples,
                    "n_responders": bench_res.n_responders,
                    "n_non_responders": bench_res.n_non_responders,
                    "roc_auc": bench_res.roc_auc,
                    "roc_auc_ci_lower": bench_res.roc_auc_ci_lower,
                    "roc_auc_ci_upper": bench_res.roc_auc_ci_upper,
                    "p_value_vs_half": bench_res.p_value_vs_half,
                    "pr_auc": bench_res.pr_auc,
                    "pr_auc_ci_lower": bench_res.pr_auc_ci_lower,
                    "pr_auc_ci_upper": bench_res.pr_auc_ci_upper,
                    "baseline_prevalence": bench_res.baseline_prevalence,
                    "delta_pr_auc": bench_res.delta_pr_auc,
                    "p_value_prauc": bench_res.p_value_prauc,
                    "brier_score": bench_res.brier_score,
                })

        cohort_elapsed = round(time.time() - t_cohort_start, 2)
        print(f"  + Cohort '{cid}' completed in {cohort_elapsed}s.")

    # -------------------------------------------------------------
    # Step 9: Standardize and Evaluate the 3 Combined Meta-Cohort Pools
    # -------------------------------------------------------------
    print("\n======================================================================")
    print("EVALUATING COMBINED META-COHORT POOLS WITH WITHIN-COHORT STANDARDIZATION")
    print("======================================================================")

    pools = {
        "Melanoma-Combined": (["Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas"], "Melanoma"),
        "RenalCell-Combined": (["McDermott-iAtlas", "Choueiri-iAtlas"], "Renal Cell Carcinoma"),
        "PanCancer-Combined": ([s[0] for s in COHORT_SPECS], "Pan-Cancer"),
    }

    for pool_id, (members, ctype) in pools.items():
        avail_members = [m for m in members if m in cohort_cache]
        if not avail_members:
            continue

        print(f"\n--- Evaluating Meta-Cohort Pool: {pool_id} ({len(avail_members)} studies) ---")
        tuples = [
            (m, cohort_cache[m]["clinical"], [cohort_cache[m]["prediction"]])
            for m in avail_members
        ]

        pool_res = pool_cohort_stratum_data(
            tuples,
            group_id=pool_id,
            cancer_type=ctype,
            standardize=True,
        )
        match pool_res:
            case Failure(err):
                print(f"  [ERROR] Pooling failed for {pool_id}: {err}")
                continue
            case Success((pooled_clin, pooled_preds)):
                pass

        p_bench = evaluate_prediction(
            prediction=pooled_preds[0],
            clinical_df=pooled_clin,
            cohort_id=pool_id,
            cancer_type=ctype,
            time_stratum="Pre",
            response_stratum="Standard",
            pooling_strategy="standardized",
            n_bootstrap=1000,
        )
        match p_bench:
            case Failure(err):
                print(f"  [ERROR] Evaluation failed for pool {pool_id}: {err}")
            case Success(bench_res):
                print(f"  * {pool_id:<22} Standardized ROC-AUC = {bench_res.roc_auc:.4f} "
                      f"[{bench_res.roc_auc_ci_lower:.3f}, {bench_res.roc_auc_ci_upper:.3f}], "
                      f"PR-AUC = {bench_res.pr_auc:.4f} (N={bench_res.n_samples})")
                new_benchmark_records.append({
                    "cohort_id": bench_res.cohort_id,
                    "cancer_type": bench_res.cancer_type,
                    "predictor_name": bench_res.predictor_name,
                    "category": bench_res.category.value,
                    "time_stratum": bench_res.time_stratum,
                    "response_stratum": bench_res.response_stratum,
                    "pooling_strategy": bench_res.pooling_strategy,
                    "n_samples": bench_res.n_samples,
                    "n_responders": bench_res.n_responders,
                    "n_non_responders": bench_res.n_non_responders,
                    "roc_auc": bench_res.roc_auc,
                    "roc_auc_ci_lower": bench_res.roc_auc_ci_lower,
                    "roc_auc_ci_upper": bench_res.roc_auc_ci_upper,
                    "p_value_vs_half": bench_res.p_value_vs_half,
                    "pr_auc": bench_res.pr_auc,
                    "pr_auc_ci_lower": bench_res.pr_auc_ci_lower,
                    "pr_auc_ci_upper": bench_res.pr_auc_ci_upper,
                    "baseline_prevalence": bench_res.baseline_prevalence,
                    "delta_pr_auc": bench_res.delta_pr_auc,
                    "p_value_prauc": bench_res.p_value_prauc,
                    "brier_score": bench_res.brier_score,
                })

    # -------------------------------------------------------------
    # Step 10: Append Records to Official Benchmark Parquet & CSV
    # -------------------------------------------------------------
    parquet_path = args.benchmark_dir / "iatlas_benchmark_metrics.parquet"
    csv_path = args.benchmark_dir / "iatlas_benchmark_metrics.csv"

    if new_benchmark_records:
        new_df = pl.DataFrame(new_benchmark_records)
        print(f"\n[Step 10] Integrating {len(new_df)} new benchmark records into {parquet_path}...")

        if parquet_path.exists():
            existing_df = pl.read_parquet(parquet_path)
            # Remove any previous records for Cohort_Optimized_Deconv to allow clean idempotent overwrite
            filtered_existing = existing_df.filter(pl.col("predictor_name") != "Cohort_Optimized_Deconv")
            updated_df = pl.concat([filtered_existing, new_df])
        else:
            updated_df = new_df

        updated_df.write_parquet(parquet_path)
        updated_df.write_csv(csv_path)
        print(f"  + Successfully saved benchmark parquet! Total rows: {len(updated_df)}")

        # -------------------------------------------------------------
        # Step 11: Re-Export Publication SVG Heatmap
        # -------------------------------------------------------------
        print("\n[Step 11] Re-exporting publication ROC-AUC Heatmap SVG...")
        pre_std = updated_df.filter(
            (pl.col("time_stratum") == "Pre") &
            (pl.col("response_stratum") == "Standard") &
            (pl.col("pooling_strategy").is_in(["cohort", "standardized"]))
        )
        chart_res = create_auc_heatmap(
            pre_std,
            metric="roc_auc",
            title="Pre-Treatment Standard ICB Response (ROC-AUC)",
        )
        match chart_res:
            case Failure(err):
                print(f"  [WARNING] Heatmap generation failed: {err}")
            case Success(chart):
                svg_out = args.benchmark_dir / "fig_roc_auc_heatmap_pre_standard.svg"
                export_res = export_chart_svg(chart, svg_out)
                match export_res:
                    case Success(_):
                        print(f"  + Exported updated heatmap SVG to: {svg_out}")
                    case Failure(err):
                        print(f"  [WARNING] SVG export failed: {err}")

    total_pipeline_time = round(time.time() - t_pipeline_start, 2)
    print("\n======================================================================")
    print(f"ALL-COHORT HPO & BENCHMARK INTEGRATION COMPLETE! Elapsed: {total_pipeline_time}s")
    print("======================================================================\n")


if __name__ == "__main__":
    main()
