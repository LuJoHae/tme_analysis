#!/usr/bin/env python3
"""Multi-Cohort Sampling, Perturbation Modeling, Harmony Integration, and Milopy Analysis.

Deterministically samples from the 9 premier single-cell immunotherapy response cohorts:
- All 9 cohorts (Pan-cancer pooled)
- By Cancer Type (Melanoma, NSCLC, Colorectal, Breast, HNSCC, PDAC)
- By Treatment / Mechanism (Anti-PD-1 monotherapy, Anti-PD-1 + CTLA-4, Chemo-immunotherapy, Combination)
- By Platform (10x Chromium 3', 10x Chromium 5', Smart-seq2, BD Rhapsody, Single Nucleus)
- Custom selection of cohort accessions

Enforces that SAMPLING, PERTURBATIONS, and HARMONIZATION are strictly and exclusively
executed via the tme_datasets API:
- tme_datasets.resolve_milopy_cohorts
- tme_datasets.load_dataset
- tme_datasets.transforms (dropout, jitter, knockout, overexpression, negative_binomial, subsample)
- tme_datasets.harmonization (align_and_concatenate, run_harmony_batch_integration)

Calculates and stores exhaustive neighborhood compositions, differential abundance metrics,
replicate prevalence gates, and triggers publication-grade diagnostic UMAP dashboards.

Strict functional Python adhering to immutability, returns Result, Polars, and Altair vector SVG.
"""

from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
import scanpy as sc

# Strict API imports
import milopy
from milopy import (
    NeighborhoodCompositionConfig,
    PermutationConfig,
    ReplicatePrevalenceConfig,
    compute_neighborhood_composition,
    compute_score_permutation_null,
    evaluate_replicate_prevalence,
    release_system_memory,
)
from tme_datasets import (
    HarmonizeConfig,
    HarmonizeMode,
    NegativeBinomialConfig,
    SubsampleSpec,
    add_expression_jitter,
    align_and_concatenate,
    in_silico_knockout,
    in_silico_overexpression,
    load_dataset,
    randomize_negative_binomial,
    resolve_milopy_cohorts,
    run_harmony_batch_integration,
    simulate_dropout,
    subsample_cells,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.milopy_analysis.plot_neighborhood_diagnostics import render_diagnostic_dashboard

# Canonical Cancer Type Lookup for the 9 Cohorts
COHORT_CANCER_MAP: dict[str, str] = {
    "GSE120575": "Melanoma",
    "CELLxGENE_7b20c613": "Melanoma",
    "CELLxGENE_05a8c945": "Colorectal",
    "CELLxGENE_6f9de485": "Breast",
    "GSE207422": "NSCLC",
    "GSE243013": "NSCLC",
    "GSE233203": "NSCLC",
    "GSE200996": "HNSCC",
    "GSE316195": "PDAC",
}


class SampledRunConfig(BaseModel):
    """Immutable configuration for multi-cohort sampling, perturbation, and Milopy analysis."""

    model_config = ConfigDict(frozen=True)

    sample_by: str = Field(..., description="Sampling axis: all, cancer_type, treatment, platform, custom")
    value: str | None = Field(default=None, description="Filter value matching the sampling axis")
    cohort_ids: tuple[str, ...] = Field(..., description="Resolved cohort accession IDs")
    cells_per_dataset: int = Field(default=20000, description="Max cells per dataset prior to concatenation (0 = unconstrained)")
    base_dir: Path = Field(default=Path("/ritom/scratch/cscs/lhaeuser/data-test"), description="Base directory containing preprocessed H5AD files")
    out_dir: Path = Field(default=Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts"), description="Results output directory")
    reports_dir: Path = Field(default=Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts"), description="Reports output directory")

    # In-silico perturbation
    perturbation: str = Field(default="none", description="Perturbation type: none, dropout, jitter, knockout, overexpression, negative_binomial")
    dropout_rate: float = Field(default=0.10, ge=0.0, le=1.0)
    jitter_sigma: float = Field(default=0.15, ge=0.0)
    knockout_genes: tuple[str, ...] = Field(default=("CD8A", "PDCD1"))
    knockout_efficiency: float = Field(default=1.0, ge=0.0, le=1.0)
    overexpression_genes: tuple[str, ...] = Field(default=("CXCL13", "IFNG"))
    overexpression_fold_change: float = Field(default=2.0, ge=0.0)
    nb_dispersion: float = Field(default=0.5, ge=0.0)

    # Harmonization
    harmonize: bool = Field(default=True, description="Enable latent space batch integration via harmonypy")
    harmony_key: str = Field(default="patient_id", description="Batch key for harmony: patient_id, dataset_id, or both")
    harmony_max_iter: int = Field(default=15, ge=1)

    # Milopy differential abundance parameters
    prop: float = Field(default=0.10, ge=0.01, le=0.50, description="Neighborhood sampling proportion")
    k: int = Field(default=30, ge=5, description="k-nearest neighbors")
    d: int = Field(default=30, ge=2, description="Number of principal components")
    fdr_threshold: float = Field(default=0.10, ge=0.001, le=0.50)
    min_replicates: int = Field(default=2, ge=1, description="Minimum independent biological donors in enriched arm")
    min_prevalence_frac: float = Field(default=0.05, ge=0.0, le=1.0, description="Minimum donor prevalence fraction in enriched arm")

    # Permutation testing
    run_permutation: bool = Field(default=True, description="Vectorized Rao Score permutation null calibration")
    n_permutations: int = Field(default=100, ge=10)
    require_permutation_significance: bool = Field(
        default=False,
        description="Require empirical permutation FDR < threshold for final significance call",
    )

    # Visualization
    plot: bool = Field(default=True, description="Generate diagnostic UMAP figure dashboard")
    title: str | None = Field(default=None, description="Custom title for diagnostic figures")
    seed: int = Field(default=42)


def parse_args() -> SampledRunConfig:
    """Parses and validates CLI arguments into an immutable SampledRunConfig."""
    parser = argparse.ArgumentParser(
        description="Multi-cohort single-cell sampling, perturbation modeling, harmony integration, and Milopy analysis."
    )
    parser.add_argument(
        "--sample-by",
        type=str,
        default="all",
        choices=["all", "cancer_type", "treatment", "platform", "custom"],
        help="Sampling axis across the 9 premier ICB cohorts",
    )
    parser.add_argument(
        "--value",
        type=str,
        default=None,
        help="Filter value for sample-by (e.g. 'melanoma', 'anti-pd1_monotherapy', '10x_chromium', or comma-separated IDs)",
    )
    parser.add_argument(
        "--cells-per-dataset",
        type=int,
        default=20000,
        help="Max cells to include per dataset before concatenation (0 = unconstrained)",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/ritom/scratch/cscs/lhaeuser/data-test"),
        help="Base directory containing preprocessed H5AD files",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts"),
        help="Output directory for publication figures",
    )

    # Perturbations
    parser.add_argument(
        "--perturbation",
        type=str,
        default="none",
        choices=["none", "dropout", "jitter", "knockout", "overexpression", "negative_binomial"],
        help="In-silico perturbation model applied strictly via tme_datasets.transforms",
    )
    parser.add_argument("--dropout-rate", type=float, default=0.10, help="Technical capture dropout rate (default: 0.10)")
    parser.add_argument("--jitter-sigma", type=float, default=0.15, help="Multiplicative log-normal jitter sigma (default: 0.15)")
    parser.add_argument("--knockout-genes", type=str, default="CD8A,PDCD1", help="Comma-separated genes for in-silico knockout")
    parser.add_argument("--knockout-efficiency", type=float, default=1.0, help="In-silico knockout efficiency (0.0 to 1.0)")
    parser.add_argument("--overexpression-genes", type=str, default="CXCL13,IFNG", help="Comma-separated genes for overexpression")
    parser.add_argument("--overexpression-fold-change", type=float, default=2.0, help="Overexpression fold change multiplier")
    parser.add_argument("--nb-dispersion", type=float, default=0.5, help="Dispersion parameter for Negative Binomial resampling")

    # Harmonization
    parser.add_argument(
        "--harmonize",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable Harmony latent batch integration (default: True)",
    )
    parser.add_argument(
        "--harmony-key",
        type=str,
        default="patient_id",
        choices=["patient_id", "dataset_id", "both"],
        help="Batch key for Harmony integration (default: patient_id)",
    )
    parser.add_argument("--harmony-max-iter", type=int, default=15, help="Max iterations for Harmony (default: 15)")

    # Milopy
    parser.add_argument("--prop", type=float, default=0.10, help="Proportion of cells sampled as neighborhood centroids")
    parser.add_argument("--k", type=int, default=30, help="k-nearest neighbors")
    parser.add_argument("--d", type=int, default=30, help="Number of principal components")
    parser.add_argument("--fdr-threshold", type=float, default=0.10, help="Nominal FDR threshold")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum independent biological donors in enriched arm")
    parser.add_argument("--min-prevalence-frac", type=float, default=0.05, help="Minimum donor prevalence fraction")

    # Permutation
    parser.add_argument(
        "--run-permutation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run Rao Score permutation testing (default: True)",
    )
    parser.add_argument("--n-permutations", type=int, default=100, help="Number of permutations for Rao score null")
    parser.add_argument(
        "--require-permutation-significance",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Require empirical permutation FDR < threshold for final significance call (default: False)",
    )

    # Plotting
    parser.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render diagnostic UMAP dashboard (default: True)",
    )
    parser.add_argument("--title", type=str, default=None, help="Custom title for diagnostic figures")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Resolve cohort IDs strictly via tme_datasets API
    cohort_res = resolve_milopy_cohorts(args.sample_by, args.value)
    match cohort_res:
        case Failure(err):
            print(f"Error resolving cohorts: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(cohort_ids):
            resolved = cohort_ids

    ko_genes = tuple(g.strip() for g in args.knockout_genes.split(",") if g.strip())
    oe_genes = tuple(g.strip() for g in args.overexpression_genes.split(",") if g.strip())

    return SampledRunConfig(
        sample_by=args.sample_by,
        value=args.value,
        cohort_ids=resolved,
        cells_per_dataset=args.cells_per_dataset,
        base_dir=args.base_dir,
        out_dir=args.out_dir,
        reports_dir=args.reports_dir,
        perturbation=args.perturbation,
        dropout_rate=args.dropout_rate,
        jitter_sigma=args.jitter_sigma,
        knockout_genes=ko_genes,
        knockout_efficiency=args.knockout_efficiency,
        overexpression_genes=oe_genes,
        overexpression_fold_change=args.overexpression_fold_change,
        nb_dispersion=args.nb_dispersion,
        harmonize=args.harmonize,
        harmony_key=args.harmony_key,
        harmony_max_iter=args.harmony_max_iter,
        prop=args.prop,
        k=args.k,
        d=args.d,
        fdr_threshold=args.fdr_threshold,
        min_replicates=args.min_replicates,
        min_prevalence_frac=args.min_prevalence_frac,
        run_permutation=args.run_permutation,
        n_permutations=args.n_permutations,
        require_permutation_significance=args.require_permutation_significance,
        plot=args.plot,
        title=args.title,
        seed=args.seed,
    )


def harmonize_cohort_metadata(obs: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Standardizes obs metadata columns across cohorts: dataset_id, clinical_response, patient_id, cancer_type."""
    obs = obs.copy()
    obs["dataset_id"] = accession
    obs["cancer_type"] = COHORT_CANCER_MAP.get(accession, "Other")

    # Response harmonization
    candidates = [
        "clinical_response",
        "response",
        "treatment_response",
        "RECIST",
        "characteristics: response",
        "characteristics: therapeutic response",
    ]
    resp_col = next((c for c in candidates if c in obs.columns), None)
    if resp_col is None:
        for c in obs.columns:
            if "response" in c.lower() or "recist" in c.lower():
                resp_col = c
                break

    def map_resp(val: Any) -> str:
        if val is None or pd.isna(val):
            return "not-evaluable"
        v = str(val).strip().lower().replace("_", " ").replace("-", " ")
        if v in {"responder", "non responder", "non-responder"}:
            return v.replace(" ", "-")
        if any(w in v for w in ["unfavourable", "progressive", "residual", "non mpr", "npcr", "nmpr", "pd", "nr"]):
            return "non-responder"
        if any(w in v for w in ["favourable", "complete response", "partial response", "pcr", "mpr", "cr", "pr", "r"]):
            return "responder"
        return "not-evaluable"

    if resp_col is not None:
        obs["clinical_response"] = obs[resp_col].map(map_resp)
    else:
        obs["clinical_response"] = "not-evaluable"

    # Patient ID harmonization
    patient_cand = [
        "patient_id",
        "patient",
        "donor",
        "donor_id",
        "sample_id",
        "sample",
        "characteristics: patient id",
        "characteristics: patient ID",
        "characteristics: patinet ID (Pre=baseline; Post= on treatment)",
    ]
    pt_col = next((c for c in patient_cand if c in obs.columns), None)
    if pt_col is None:
        for c in obs.columns:
            if any(k in c.lower() for k in ("patient", "donor", "patinet")):
                pt_col = c
                break

    if pt_col is not None:
        obs["patient_id"] = f"{accession}_" + obs[pt_col].astype(str).str.replace(r"^(?:Pre|Post)_", "", regex=True)
    else:
        obs["patient_id"] = f"{accession}_pt_unknown"

    return obs


def apply_cohort_perturbation(
    adata: ad.AnnData,
    config: SampledRunConfig,
) -> Result[ad.AnnData, str]:
    """Applies in-silico perturbation strictly via tme_datasets.transforms API."""
    match config.perturbation:
        case "none":
            return Success(adata)
        case "dropout":
            return simulate_dropout(adata, rate=config.dropout_rate, seed=config.seed)
        case "jitter":
            return add_expression_jitter(adata, sigma=config.jitter_sigma, seed=config.seed)
        case "knockout":
            return in_silico_knockout(
                adata,
                genes=config.knockout_genes,
                efficiency=config.knockout_efficiency,
            )
        case "overexpression":
            return in_silico_overexpression(
                adata,
                genes=config.overexpression_genes,
                fold_change=config.overexpression_fold_change,
            )
        case "negative_binomial":
            nb_cfg = NegativeBinomialConfig(
                dispersion=config.nb_dispersion,
                seed=Some(config.seed),
            )
            return randomize_negative_binomial(adata, nb_cfg)
        case other:
            return Failure(f"Unknown perturbation type: {other}")


def load_and_prepare_cohort(
    accession: str,
    config: SampledRunConfig,
) -> Result[ad.AnnData, str]:
    """Loads, standardizes, perturbs, and subsamples a cohort strictly via tme_datasets API."""
    print(f"[{accession}] Loading via tme_datasets.load_dataset...")
    load_res = load_dataset(accession, base_dir=config.base_dir, auto_download=False)
    match load_res:
        case Failure(err):
            return Failure(f"Could not load {accession}: {err}")
        case Success(loaded):
            raw_adata = loaded

    # Standardize metadata
    raw_adata.obs = harmonize_cohort_metadata(raw_adata.obs, accession)

    # Filter binary response
    valid_mask = raw_adata.obs["clinical_response"].isin(["responder", "non-responder"])
    n_valid = int(valid_mask.sum())
    if n_valid < 30:
        return Failure(f"{accession} has insufficient binary response cells ({n_valid} found)")

    # Subsample indices directly BEFORE slicing raw_adata to prevent memory exhaustion
    if 0 < config.cells_per_dataset < n_valid:
        print(f"[{accession}] Subsampling cells to budget {config.cells_per_dataset:,} (total valid: {n_valid:,})...")
        valid_indices = np.where(valid_mask.values)[0]
        rng = np.random.default_rng(config.seed)
        sub_resp = raw_adata.obs["clinical_response"].iloc[valid_indices].values

        chosen_indices: list[int] = []
        for resp_label in ["responder", "non-responder"]:
            resp_idx = valid_indices[sub_resp == resp_label]
            if len(resp_idx) == 0:
                continue
            k_sel = max(1, int(round(len(resp_idx) / len(valid_indices) * config.cells_per_dataset)))
            k_sel = min(len(resp_idx), k_sel)
            chosen_indices.extend(rng.choice(resp_idx, size=k_sel, replace=False).tolist())

        chosen_arr = np.sort(np.array(chosen_indices))
        if len(chosen_arr) > config.cells_per_dataset:
            chosen_arr = np.sort(rng.choice(chosen_arr, size=config.cells_per_dataset, replace=False))

        sub_adata = raw_adata[chosen_arr].copy()
    else:
        sub_adata = raw_adata[valid_mask].copy()

    del raw_adata
    release_system_memory()

    # Apply in-silico perturbation if requested
    if config.perturbation != "none":
        print(f"[{accession}] Applying perturbation '{config.perturbation}' via tme_datasets.transforms...")
        perturb_res = apply_cohort_perturbation(sub_adata, config)
        match perturb_res:
            case Failure(err):
                return Failure(f"Perturbation failed on {accession}: {err}")
            case Success(perturbed):
                sub_adata = perturbed

    # Prefix cell barcodes to guarantee global uniqueness
    sub_adata.obs_names = [f"{accession}_{name}" for name in sub_adata.obs_names]
    sub_adata.obs_names_make_unique()
    print(f"[{accession}] Retained {sub_adata.n_obs:,} cells across {sub_adata.n_vars:,} genes.")
    return Success(sub_adata)


def compute_nhood_centroids(adata: ad.AnnData) -> tuple[np.ndarray, np.ndarray]:
    """Computes exact centroid UMAP coordinates for all neighborhoods via sparse projection."""
    nhoods_mat = adata.obsm["nhoods"]
    if not sp.isspmatrix_csc(nhoods_mat):
        nhoods_mat = sp.csc_matrix(nhoods_mat)

    cell_umap = adata.obsm["X_umap"]
    nhood_sizes = np.asarray(nhoods_mat.sum(axis=0)).ravel()
    safe_sizes = np.where(nhood_sizes > 0, nhood_sizes, 1.0)

    # nhoods_mat: (N_cells, M_nhoods), cell_umap: (N_cells, 2)
    # (nhoods_mat.T @ cell_umap) -> (M_nhoods, 2)
    nhood_coords = np.asarray(nhoods_mat.T.dot(cell_umap)) / safe_sizes[:, np.newaxis]
    return nhood_coords[:, 0].astype(np.float32), nhood_coords[:, 1].astype(np.float32)


def project_single_cell_metrics(
    adata: ad.AnnData,
    res_df: pd.DataFrame,
    fdr_thresh: float = 0.10,
) -> pl.DataFrame:
    """Projects neighborhood logFC and FDR onto single cells via sparse weighted projection."""
    nhoods_mat = adata.obsm["nhoods"]
    if not sp.isspmatrix_csr(nhoods_mat):
        nhoods_csr = sp.csr_matrix(nhoods_mat)
    else:
        nhoods_csr = nhoods_mat

    n_cells = adata.n_obs
    logfc_vec = res_df["logFC"].fillna(0.0).values
    fdr_vec = res_df["FDR"].fillna(1.0).values

    cell_nhood_counts = np.asarray(nhoods_csr.sum(axis=1)).ravel()
    safe_counts = np.where(cell_nhood_counts > 0, cell_nhood_counts, 1.0)

    # Mean logFC across covering neighborhoods
    cell_logfc = np.asarray(nhoods_csr.dot(logfc_vec)).ravel() / safe_counts
    cell_logfc = np.where(cell_nhood_counts > 0, cell_logfc, 0.0)

    # Minimum FDR across covering neighborhoods
    cell_min_fdr = np.ones(n_cells, dtype=np.float32)
    rows, cols = nhoods_csr.nonzero()
    for r, c in zip(rows, cols):
        if fdr_vec[c] < cell_min_fdr[r]:
            cell_min_fdr[r] = fdr_vec[c]

    is_sig = (cell_min_fdr < fdr_thresh) & (cell_nhood_counts > 0)
    umap_coords = adata.obsm["X_umap"]

    return pl.DataFrame({
        "cell_id": adata.obs_names.tolist(),
        "dataset_id": adata.obs["dataset_id"].astype(str).tolist(),
        "patient_id": adata.obs["patient_id"].astype(str).tolist(),
        "clinical_response": adata.obs["clinical_response"].astype(str).tolist(),
        "cancer_type": adata.obs.get("cancer_type", pd.Series(["Unknown"] * n_cells)).astype(str).tolist(),
        "UMAP1": umap_coords[:, 0].astype(np.float32),
        "UMAP2": umap_coords[:, 1].astype(np.float32),
        "cell_logfc": cell_logfc.astype(np.float32),
        "min_fdr": cell_min_fdr.astype(np.float32),
        "is_significant": is_sig,
        "nhood_count": cell_nhood_counts.astype(np.int32),
    })


def run_sampling_pipeline(config: SampledRunConfig) -> Result[Path, str]:
    """Executes deterministic multi-cohort sampling, perturbation, harmonization, and Milopy analysis."""
    start_time = time.time()
    run_name = (
        f"{config.sample_by}_{config.value}_{config.perturbation}"
        if config.value
        else f"{config.sample_by}_{config.perturbation}"
    ).replace("/", "_").replace(" ", "_")

    run_out_dir = config.out_dir / run_name
    run_reports_dir = config.reports_dir / run_name
    run_out_dir.mkdir(parents=True, exist_ok=True)
    run_reports_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"MULTI-COHORT SAMPLING AND MILOPY ANALYSIS: {run_name.upper()}")
    print("=" * 80)
    print(f"Sampling Axis     : {config.sample_by} (filter: {config.value})")
    print(f"Resolved Cohorts  : {config.cohort_ids}")
    print(f"Max Cells/Cohort  : {config.cells_per_dataset}")
    print(f"Perturbation      : {config.perturbation}")
    print(f"Harmony Enabled   : {config.harmonize} (key: {config.harmony_key})")
    print(f"Output Directory  : {run_out_dir}")

    # Step 1: Load, perturb, and subsample candidate cohorts
    loaded_adatas: list[ad.AnnData] = []
    loaded_ids: list[str] = []

    for acc in config.cohort_ids:
        load_res = load_and_prepare_cohort(acc, config)
        match load_res:
            case Failure(err):
                print(f"[{acc}] Warning: Skipped cohort: {err}", file=sys.stderr)
            case Success(c_adata):
                loaded_adatas.append(c_adata)
                loaded_ids.append(acc)

    if not loaded_adatas:
        return Failure(f"Failed to load any cohorts for run {run_name}")

    # Step 2: Feature alignment & concatenation via tme_datasets API
    if len(loaded_adatas) > 1:
        print("\n" + "-" * 70)
        print(f"HARMONIZING {len(loaded_adatas)} COHORTS VIA tme_datasets.align_and_concatenate...")
        print("-" * 70)
        min_genes_present = min(a.n_vars for a in loaded_adatas)
        harmonize_config = HarmonizeConfig(
            mode=HarmonizeMode.INTERSECTION,
            batch_key="dataset_id",
            reconcile_genes=False,
            min_shared_genes=min(50, min_genes_present),
        )
        concat_res = align_and_concatenate(loaded_adatas, loaded_ids, harmonize_config)
        match concat_res:
            case Failure(err):
                return Failure(f"align_and_concatenate failed: {err}")
            case Success(adata_combined):
                combined_adata = adata_combined
    else:
        print(f"\nSingle cohort loaded ({loaded_ids[0]}); proceeding without multi-dataset concatenation.")
        combined_adata = loaded_adatas[0]

    del loaded_adatas
    release_system_memory()

    n_cells = combined_adata.n_obs
    n_genes = combined_adata.n_vars
    print(f"[PREPARED] Combined single-cell matrix: {n_cells:,} cells x {n_genes:,} genes.")

    # Step 3: Expression normalization, HVG selection, and joint PCA
    print("\nPreparing feature matrix and computing PCA...")
    if not sp.issparse(combined_adata.X):
        combined_adata.X = sp.csr_matrix(combined_adata.X)

    max_expr = combined_adata.X.data.max() if len(combined_adata.X.data) > 0 else 0
    if max_expr > 50.0:
        print("Expression values appear unlogged; applying log1p...")
        sc.pp.log1p(combined_adata)

    n_top = min(2000, n_genes)
    if len(loaded_ids) > 1 and "dataset_id" in combined_adata.obs.columns:
        try:
            sc.pp.highly_variable_genes(
                combined_adata,
                n_top_genes=n_top,
                flavor="seurat",
                batch_key="dataset_id",
                subset=False,
            )
        except Exception:
            sc.pp.highly_variable_genes(combined_adata, n_top_genes=n_top, flavor="seurat", subset=False)
    else:
        sc.pp.highly_variable_genes(combined_adata, n_top_genes=n_top, flavor="seurat", subset=False)

    pca_d = min(config.d, n_genes - 1, n_cells - 1)
    sc.pp.pca(combined_adata, n_comps=pca_d, mask_var="highly_variable", random_state=config.seed)

    # Step 4: Harmony latent batch integration strictly via tme_datasets API
    use_rep = "X_pca"
    if config.harmonize:
        print("\n" + "-" * 70)
        print("RUNNING HARMONY LATENT BATCH INTEGRATION (tme_datasets.run_harmony_batch_integration)...")
        print("-" * 70)

        # Determine batch keys
        if config.harmony_key == "both" and len(loaded_ids) > 1:
            batch_keys: str | Sequence[str] = ["dataset_id", "patient_id"]
        elif config.harmony_key == "dataset_id" and len(loaded_ids) > 1:
            batch_keys = "dataset_id"
        else:
            batch_keys = "patient_id"

        harmony_res = run_harmony_batch_integration(
            combined_adata,
            batch_key=batch_keys,
            n_comps=pca_d,
            max_iter=config.harmony_max_iter,
            seed=config.seed,
        )
        match harmony_res:
            case Failure(err):
                print(f"Harmony integration warning: {err}; continuing with raw PCA.", file=sys.stderr)
            case Success(adata_harm):
                combined_adata = adata_harm
                if "X_pca_harmony" in combined_adata.obsm:
                    use_rep = "X_pca_harmony"

    # Step 5: kNN Graph and UMAP
    print(f"\nConstructing kNN graph (k={config.k}, rep={use_rep}) and UMAP embedding...")
    sc.pp.neighbors(combined_adata, n_neighbors=config.k, n_pcs=pca_d, use_rep=use_rep, random_state=config.seed)
    sc.tl.umap(combined_adata, min_dist=0.3, random_state=config.seed)

    # Step 6: Milopy Neighborhood Definition and Cell Counting
    print(f"\nDefining Milopy neighborhoods (prop={config.prop})...")
    milopy.core.make_nhoods(
        combined_adata,
        prop=config.prop,
        k=config.k,
        d=pca_d,
        refined=True,
        random_state=config.seed,
    )
    n_nhoods = int(combined_adata.obsm["nhoods"].shape[1])
    print(f"Defined {n_nhoods:,} neighborhoods across {n_cells:,} single cells.")

    combined_adata.obs["patient_id"] = pd.Categorical(combined_adata.obs["patient_id"].astype(str))
    milopy.core.count_cells(combined_adata, sample_col="patient_id")

    # Step 7: Differential Abundance GLM
    print("\nFitting Quasi-Likelihood Negative Binomial GLM...")
    nhood_counts_cols = list(combined_adata.uns["nhood_counts"].columns)
    design_df = (
        combined_adata.obs[["patient_id", "dataset_id", "clinical_response"]]
        .drop_duplicates(subset=["patient_id"])
        .set_index("patient_id")
        .loc[nhood_counts_cols]
    )

    n_datasets_in_glm = len(design_df["dataset_id"].unique())
    design_formula = "~ dataset_id + clinical_response" if n_datasets_in_glm > 1 else "~ clinical_response"
    print(f"Design Formula : {design_formula}")
    print(f"Patients in GLM: {len(design_df):,} ({sum(design_df['clinical_response'] == 'responder')} R / "
          f"{sum(design_df['clinical_response'] == 'non-responder')} NR)")

    milopy.core.test_nhoods(combined_adata, design=design_formula, design_df=design_df)
    res_df = combined_adata.uns["nhood_test_results"].copy()
    res_df.index.name = "Nhood"

    # Step 8: Vectorized Neighborhood Composition (Datasets & Patients)
    print("\nComputing exact neighborhood compositions via sparse projection C = B^T A...")
    ds_comp_res = compute_neighborhood_composition(combined_adata, category_col="dataset_id")
    match ds_comp_res:
        case Failure(err):
            print(f"Dataset composition warning: {err}", file=sys.stderr)
            ds_comp_df = pl.DataFrame({"nhood_index": np.arange(n_nhoods)})
        case Success(ds_comp):
            ds_comp_df = ds_comp.metrics_df
            ds_pdf = ds_comp.metrics_df.to_pandas().set_index("nhood_index")
            # Simpson purity S = sum p^2 = 1 - simpson_diversity
            simpson_purity = 1.0 - ds_pdf["dataset_id_simpson_diversity"].values
            res_df["dataset_id_purity"] = simpson_purity
            res_df["dataset_purity"] = simpson_purity
            res_df["dominant_dataset"] = ds_pdf["dominant_dataset_id"].values
            res_df["dataset_id_entropy"] = ds_pdf["dataset_id_entropy"].values

    pt_comp_res = compute_neighborhood_composition(combined_adata, category_col="patient_id")
    match pt_comp_res:
        case Failure(err):
            print(f"Patient composition warning: {err}", file=sys.stderr)
            pt_comp_df = pl.DataFrame({"nhood_index": np.arange(n_nhoods)})
        case Success(pt_comp):
            pt_comp_df = pt_comp.metrics_df
            pt_pdf = pt_comp.metrics_df.to_pandas().set_index("nhood_index")
            # Max single-donor fraction
            max_donor_frac = pt_pdf["patient_id_purity"].values
            res_df["patient_id_max_fraction"] = max_donor_frac
            res_df["patient_max_frac"] = max_donor_frac
            res_df["dominant_patient"] = pt_pdf["dominant_patient_id"].values
            res_df["patient_id_entropy"] = pt_pdf["patient_id_entropy"].values

    # Centroid UMAP Coordinates
    umap1_cent, umap2_cent = compute_nhood_centroids(combined_adata)
    res_df["UMAP1"] = umap1_cent
    res_df["UMAP2"] = umap2_cent
    res_df["Nhood_UMAP1"] = umap1_cent
    res_df["Nhood_UMAP2"] = umap2_cent

    # Step 9: Rao Score Permutation Testing
    perm_pvals: Sequence[float] | None = None
    perm_fdrs: Sequence[float] | None = None
    if config.run_permutation:
        print(f"\nRunning vectorized Rao Score permutation calibration (B={config.n_permutations})...")
        nhood_counts_df = combined_adata.uns["nhood_counts"]
        count_mat = nhood_counts_df.values.astype(np.float64)
        pts_in_order = list(nhood_counts_df.columns)
        resp_in_order = [str(design_df.loc[p, "clinical_response"]).lower() for p in pts_in_order]

        perm_res = compute_score_permutation_null(
            count_matrix=count_mat,
            response_labels=resp_in_order,
            library_sizes=count_mat.sum(axis=0),
            config=PermutationConfig(
                n_permutations=config.n_permutations,
                seed=config.seed,
                fdr_threshold=config.fdr_threshold,
            ),
        )
        match perm_res:
            case Success(p_out):
                perm_pvals = p_out.permutation_pvalues
                perm_fdrs = p_out.permutation_fdr
                n_perm_sig = sum(p_out.is_perm_significant)
                print(f"Permutation test complete: {n_perm_sig} neighborhoods pass empirical FDR < {config.fdr_threshold}.")
            case Failure(err):
                print(f"Warning: permutation test skipped ({err})", file=sys.stderr)

    # Step 10: Replicate Prevalence and Multi-Cohort Recurrence Filtering
    print("\nApplying biological replicate prevalence and clonal spike filtering...")
    prev_res = evaluate_replicate_prevalence(
        adata=combined_adata,
        res_df=res_df,
        patient_col="patient_id",
        design_col="clinical_response",
        dataset_col="dataset_id",
        permutation_pvalues=perm_pvals,
        permutation_fdr=perm_fdrs,
        config=ReplicatePrevalenceConfig(
            min_replicates=config.min_replicates,
            min_prevalence_frac=config.min_prevalence_frac,
            fdr_threshold=config.fdr_threshold,
            require_permutation_significance=config.require_permutation_significance,
        ),
    )
    match prev_res:
        case Success(annotated_df):
            res_df = annotated_df
        case Failure(err):
            print(f"Prevalence filtering warning: {err}", file=sys.stderr)

    # Add numeric replicate prevalence alias if needed
    if "replicate_prevalence" not in res_df.columns:
        res_df["replicate_prevalence"] = res_df["n_patients_total"].values

    n_sig_up = int(((res_df["is_significant"]) & (res_df["logFC"] > 0)).sum())
    n_sig_down = int(((res_df["is_significant"]) & (res_df["logFC"] < 0)).sum())
    n_conserved = int((res_df["hit_type"] == "Conserved Recurrent Hit").sum())
    n_private_spikes = int((res_df["hit_type"] == "Private Clonal Spike").sum())

    print(
        f"\n[DA SUMMARY] Total Nhoods: {len(res_df):,} | Responders (Up): {n_sig_up:,} | "
        f"Non-Responders (Down): {n_sig_down:,} | Conserved Recurrent: {n_conserved:,} | Clonal Spikes Blocked: {n_private_spikes:,}"
    )

    # Step 11: Export Structured Parquet Tables & Sparse Neighborhood Matrix
    print("\nExporting structured Parquet results...")
    if "nhoods" in combined_adata.obsm:
        nhood_mat_path = run_out_dir / "nhood_matrix.npz"
        sp.save_npz(nhood_mat_path, sp.csc_matrix(combined_adata.obsm["nhoods"]))
        print(f"[EXPORTED] {nhood_mat_path}")

    da_parquet_path = run_out_dir / "da_results.parquet"
    res_pl = pl.from_pandas(res_df.reset_index())
    res_pl.write_parquet(da_parquet_path)
    print(f"[EXPORTED] {da_parquet_path}")

    ds_comp_path = run_out_dir / "nhood_composition_datasets.parquet"
    ds_comp_df.write_parquet(ds_comp_path)
    print(f"[EXPORTED] {ds_comp_path}")

    pt_comp_path = run_out_dir / "nhood_composition_patients.parquet"
    pt_comp_df.write_parquet(pt_comp_path)
    print(f"[EXPORTED] {pt_comp_path}")

    # Project cell-level scores
    cell_df = project_single_cell_metrics(combined_adata, res_df, fdr_thresh=config.fdr_threshold)
    cell_parquet_path = run_out_dir / "cell_level_scores.parquet"
    cell_df.write_parquet(cell_parquet_path)
    print(f"[EXPORTED] {cell_parquet_path}")

    # Summary JSON
    duration = time.time() - start_time
    summary_data = {
        "run_name": run_name,
        "timestamp": datetime.datetime.now().isoformat(),
        "duration_seconds": round(duration, 2),
        "sampling": {
            "sample_by": config.sample_by,
            "filter_value": config.value,
            "cohort_ids": list(config.cohort_ids),
            "loaded_cohorts": loaded_ids,
            "max_cells_per_cohort": config.cells_per_dataset,
        },
        "perturbation": {
            "type": config.perturbation,
            "dropout_rate": config.dropout_rate,
            "jitter_sigma": config.jitter_sigma,
            "knockout_genes": list(config.knockout_genes),
            "overexpression_genes": list(config.overexpression_genes),
            "nb_dispersion": config.nb_dispersion,
        },
        "harmonization": {
            "harmonize": config.harmonize,
            "harmony_key": config.harmony_key,
        },
        "metrics": {
            "total_cells": n_cells,
            "total_genes": n_genes,
            "total_nhoods": len(res_df),
            "sig_responder_nhoods": n_sig_up,
            "sig_non_responder_nhoods": n_sig_down,
            "conserved_recurrent_hits": n_conserved,
            "private_clonal_spikes_gated": n_private_spikes,
        },
    }
    summary_path = run_out_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary_data, indent=2), encoding="utf-8")
    print(f"[EXPORTED] {summary_path}")

    # Step 12: Generate Diagnostic Visualizations
    if config.plot:
        print("\n" + "-" * 70)
        print("RENDERING PUBLICATION-GRADE DIAGNOSTIC DASHBOARD (plot_neighborhood_diagnostics)...")
        print("-" * 70)
        fig_title = config.title or f"Multi-Cohort Analysis: {run_name.replace('_', ' ').title()}"
        report_svg = run_reports_dir / f"{run_name}_diagnostics.svg"

        plot_res = render_diagnostic_dashboard(
            da_df=res_pl,
            cell_df=cell_df,
            output_svg=report_svg,
            title=fig_title,
            fdr_thresh=config.fdr_threshold,
            max_cells=15000,
        )
        match plot_res:
            case Success(path):
                print(f"[DASHBOARD] Diagnostic dashboard saved successfully to {path}")
            case Failure(err):
                print(f"[DASHBOARD WARNING] Failed to render dashboard: {err}", file=sys.stderr)

    print("\n" + "=" * 80)
    print(f"PIPELINE COMPLETED SUCCESSFULLY IN {duration:.1f} SECONDS")
    print("=" * 80)
    return Success(da_parquet_path)


def main() -> None:
    config = parse_args()
    res = run_sampling_pipeline(config)
    match res:
        case Success(path):
            print(f"\nResult files written to: {path.parent}")
            sys.exit(0)
        case Failure(err):
            print(f"\nPipeline failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
