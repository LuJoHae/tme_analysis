#!/usr/bin/env python3
"""Cross-Cohort Combination, Harmonization, and Milopy Differential Abundance Analysis.

Combines single-cell immunotherapy cohorts grouped by cancer type (Melanoma, NSCLC)
or pooled pan-cancer STRICTLY AND EXCLUSIVELY using the API provided by tme_datasets:
- align_and_concatenate
- sample_and_harmonize_cohorts
- evaluate_integration_metrics
- HarmonizeConfig, HarmonizeMode, GeneIDType

Computes vectorized neighborhood dataset purity, diversity, and entropy metrics via
compute_nhood_composition, and models technical batch variation in Milopy GLM:
~ dataset_id + clinical_response

Adheres strictly to functional Python, returns Result, Polars, and Nature Methods vector SVG.
"""

from __future__ import annotations

import argparse
import ctypes
import gc
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence
import xml.etree.ElementTree as ET

import altair as alt  # type: ignore
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp
import scanpy as sc
import vl_convert as vlc  # type: ignore

import milopy
from milopy import (
    NeighborhoodCompositionConfig,
    annotate_nhood_adata_composition,
    compute_neighborhood_composition,
    compute_score_permutation_null,
    evaluate_replicate_prevalence,
    PermutationConfig,
    ReplicatePrevalenceConfig,
    release_system_memory,
)
from tme_datasets import (
    HarmonizeConfig,
    HarmonizeMode,
    GeneIDType,
    align_and_concatenate,
    evaluate_integration_metrics,
    load_dataset,
)
from tme_datasets.preprocessing import standardize_timepoint
from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
)

COHORT_GROUPS: dict[str, tuple[str, ...]] = {
    "melanoma": ("GSE120575", "CELLxGENE_7b20c613"),
    "nsclc": ("GSE207422", "GSE243013", "GSE233203"),
    "solid_tumors": (
        "GSE316195",
        "CELLxGENE_05a8c945",
        "CELLxGENE_6f9de485",
        "GSE200996",
    ),
    "pancancer": (
        "GSE120575",
        "CELLxGENE_7b20c613",
        "CELLxGENE_05a8c945",
        "CELLxGENE_6f9de485",
        "GSE207422",
        "GSE243013",
        "GSE233203",
        "GSE200996",
        "GSE316195",
    ),
}

COHORT_CANCER_TYPES: dict[str, str] = {
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


class CombinedRunConfig(BaseModel):
    """Configuration for cross-cohort combination and Milopy analysis."""

    model_config = ConfigDict(frozen=True)

    group: str
    cohort_ids: tuple[str, ...]
    base_dir: Path
    out_dir: Path
    reports_dir: Path
    max_cells_per_cohort: int = 25000
    prop: float = 0.10
    k: int = 30
    d: int = 30
    fdr_threshold: float = 0.10
    min_replicates: int = 2
    min_prevalence_frac: float = 0.05
    run_permutation: bool = True
    n_permutations: int = 500
    seed: int = 42


def parse_args() -> CombinedRunConfig:
    parser = argparse.ArgumentParser(
        description="Cross-cohort harmonization and Milopy differential abundance analysis using tme_datasets."
    )
    parser.add_argument(
        "--group",
        type=str,
        required=True,
        choices=["melanoma", "nsclc", "solid_tumors", "pancancer"],
        help="Cohort combination group",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/storage/halu/data-test"),
        help="Base directory containing preprocessed H5AD files",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/results/milopy_combined_analysis"),
        help="Output directory for combined analysis results",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("output/reports"),
        help="Directory to save publication figures",
    )
    parser.add_argument(
        "--max-cells-per-cohort",
        type=int,
        default=25000,
        help="Max cells to include per cohort before concatenation (0 for unlimited)",
    )
    parser.add_argument(
        "--prop",
        type=float,
        default=0.10,
        help="Proportion of cells to sample as neighborhood representatives",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=30,
        help="k-nearest neighbors",
    )
    parser.add_argument(
        "--d",
        type=int,
        default=30,
        help="Number of principal components",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Significance threshold for FDR",
    )
    parser.add_argument(
        "--min-replicates",
        type=int,
        default=2,
        help="Minimum independent biological donors in enriched arm (default: 2)",
    )
    parser.add_argument(
        "--min-prevalence-frac",
        type=float,
        default=0.05,
        help="Minimum donor prevalence fraction in enriched arm (default: 0.05)",
    )
    parser.add_argument(
        "--run-permutation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable vectorized Rao Score permutation testing across donor labels (default: True)",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=500,
        help="Number of donor label permutations for null distribution calibration (default: 500)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )
    args = parser.parse_args()

    cohort_ids = COHORT_GROUPS[args.group]
    return CombinedRunConfig(
        group=args.group,
        cohort_ids=cohort_ids,
        base_dir=args.base_dir,
        out_dir=args.out_dir,
        reports_dir=args.reports_dir,
        max_cells_per_cohort=args.max_cells_per_cohort,
        prop=args.prop,
        k=args.k,
        d=args.d,
        fdr_threshold=args.fdr_threshold,
        min_replicates=args.min_replicates,
        min_prevalence_frac=args.min_prevalence_frac,
        run_permutation=args.run_permutation,
        n_permutations=args.n_permutations,
        seed=args.seed,
    )


def harmonize_cohort_obs(obs: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Harmonizes clinical_response, patient_id, dataset_id, and cancer_type."""
    obs = obs.copy()
    obs["dataset_id"] = accession
    obs["cancer_type"] = COHORT_CANCER_TYPES.get(accession, "Other")

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
    ]
    pt_col = next((c for c in patient_cand if c in obs.columns), None)
    if pt_col is None:
        for c in obs.columns:
            if "patient" in c.lower() or "donor" in c.lower():
                pt_col = c
                break

    if pt_col is not None:
        obs["patient_id"] = f"{accession}_" + obs[pt_col].astype(str).str.replace(r"^(?:Pre|Post)_", "", regex=True)
    else:
        obs["patient_id"] = f"{accession}_pt_unknown"

    return obs


def load_and_subsample_cohort(
    accession: str,
    config: CombinedRunConfig,
) -> Result[ad.AnnData, str]:
    """Loads a cohort strictly via tme_datasets, standardizes obs, filters binary response, and subsamples cells."""
    print(f"[{accession}] Loading via tme_datasets API...")
    load_res = load_dataset(accession, base_dir=config.base_dir, auto_download=False)
    match load_res:
        case Failure(err):
            return Failure(f"Could not load {accession}: {err}")
        case Success(loaded):
            raw_adata = loaded

    obs_h = harmonize_cohort_obs(raw_adata.obs, accession)
    raw_adata.obs = obs_h

    # Filter binary response
    valid_mask = raw_adata.obs["clinical_response"].isin(["responder", "non-responder"])
    if valid_mask.sum() < 50:
        return Failure(f"{accession} has insufficient binary response cells ({valid_mask.sum()} found)")

    sub_adata = raw_adata[valid_mask].copy()
    del raw_adata
    release_system_memory()

    # Stratified subsampling if exceeding budget
    if 0 < config.max_cells_per_cohort < sub_adata.n_obs:
        rng = np.random.default_rng(config.seed)
        n_target = config.max_cells_per_cohort
        frac = n_target / sub_adata.n_obs

        selected_keys: list[str] = []
        for _, grp in sub_adata.obs.groupby(["patient_id", "clinical_response"], observed=True):
            n_sel = max(1, int(round(len(grp) * frac)))
            if n_sel >= len(grp):
                selected_keys.extend(grp.index.tolist())
            else:
                chosen = rng.choice(grp.index.values, size=n_sel, replace=False)
                selected_keys.extend(chosen.tolist())

        if len(selected_keys) > n_target:
            selected_keys = rng.choice(selected_keys, size=n_target, replace=False).tolist()

        sub_adata = sub_adata[selected_keys].copy()
        release_system_memory()

    # Make cell names globally unique with cohort prefix
    sub_adata.obs_names = [f"{accession}_{name}" for name in sub_adata.obs_names]
    sub_adata.obs_names_make_unique()
    print(f"[{accession}] Retained {sub_adata.n_obs:,} cells across {sub_adata.n_vars:,} genes.")
    return Success(sub_adata)


def run_combined_analysis(config: CombinedRunConfig) -> Result[Path, str]:
    """Orchestrates multi-cohort harmonization via tme_datasets, Milopy DA testing, and reporting."""
    print("=" * 70)
    print(f"STARTING COMBINED COHORT ANALYSIS: GROUP = {config.group.upper()}")
    print("=" * 70)
    print(f"Target Cohorts : {config.cohort_ids}")
    print(f"Max Cells/Cohort : {config.max_cells_per_cohort}")

    loaded_adatas: list[ad.AnnData] = []
    loaded_ids: list[str] = []

    for acc in config.cohort_ids:
        res = load_and_subsample_cohort(acc, config)
        match res:
            case Failure(err):
                print(f"[{acc}] Warning: Skipped cohort: {err}", file=sys.stderr)
            case Success(adata):
                loaded_adatas.append(adata)
                loaded_ids.append(acc)

    if len(loaded_adatas) < 2 and config.group != "melanoma":
        return Failure(f"Fewer than 2 cohorts loaded for group {config.group} ({loaded_ids})")

    # Step 1: Harmonize and Concatenate strictly via tme_datasets API
    print("\n" + "-" * 70)
    print(f"HARMONIZING {len(loaded_adatas)} COHORTS VIA tme_datasets.harmonization.align_and_concatenate...")
    print("-" * 70)
    harmonize_config = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        batch_key="dataset_id",
        reconcile_genes=False,  # Cohorts are already processed with standardized features
        min_shared_genes=500,
    )

    concat_res = align_and_concatenate(loaded_adatas, loaded_ids, harmonize_config)
    match concat_res:
        case Failure(err):
            return Failure(f"align_and_concatenate failed: {err}")
        case Success(combined_adata):
            pass

    del loaded_adatas
    release_system_memory()

    n_total_cells = combined_adata.n_obs
    n_shared_genes = combined_adata.n_vars
    print(f"[HARMONIZED] Combined dataset: {n_total_cells:,} cells x {n_shared_genes:,} shared genes.")

    # Evaluate integration quality metrics via tme_datasets API
    print("\n" + "-" * 70)
    print("EVALUATING INTEGRATION QUALITY METRICS (tme_datasets.harmonization.evaluate_integration_metrics)...")
    print("-" * 70)

    # Compute joint PCA for kNN, UMAP, and metrics
    if not sp.issparse(combined_adata.X):
        combined_adata.X = sp.csr_matrix(combined_adata.X)
    max_val = combined_adata.X.data.max() if len(combined_adata.X.data) > 0 else 0
    if max_val > 50.0:
        print("Expression values appear unlogged; applying log1p...")
        sc.pp.log1p(combined_adata)

    print("Selecting 2,000 Highly Variable Genes (batch-aware)...")
    try:
        sc.pp.highly_variable_genes(
            combined_adata,
            n_top_genes=2000,
            flavor="seurat",
            batch_key="dataset_id",
            subset=False,
        )
    except Exception:
        sc.pp.highly_variable_genes(combined_adata, n_top_genes=2000, flavor="seurat", subset=False)

    print("Computing joint PCA (30 components)...")
    sc.pp.pca(combined_adata, n_comps=config.d, mask_var="highly_variable", random_state=config.seed)

    label_key = "cell_type" if "cell_type" in combined_adata.obs.columns else "clinical_response"
    metrics_res = evaluate_integration_metrics(
        combined_adata,
        batch_key="dataset_id",
        label_key=label_key,
        k=30,
    )
    match metrics_res:
        case Success(m):
            print(f"Integration Metrics: iLISI={m.mean_ilisi:.2f}, cLISI={m.mean_clisi:.2f}, "
                  f"Batch Silhouette={m.batch_silhouette:.3f}, kBET={m.kbet_acceptance_rate:.3f}")
        case Failure(err):
            print(f"Integration evaluation warning: {err}", file=sys.stderr)

    # Step 2: kNN Graph and UMAP
    print("Constructing kNN graph and UMAP embedding...")
    sc.pp.neighbors(combined_adata, n_neighbors=config.k, n_pcs=config.d, random_state=config.seed)
    sc.tl.umap(combined_adata, min_dist=0.3, random_state=config.seed)

    # Step 3: Milopy Neighborhood Definition
    print(f"\nDefining Milopy neighborhoods (prop={config.prop})...")
    milopy.core.make_nhoods(
        combined_adata,
        prop=config.prop,
        k=config.k,
        d=config.d,
        refined=True,
        random_state=config.seed,
    )
    n_nhoods = int(combined_adata.obsm["nhoods"].shape[1])
    print(f"Defined {n_nhoods:,} neighborhoods across {n_total_cells:,} cells.")

    # Step 4: Count Cells per Patient
    print("Counting cells per patient in neighborhoods...")
    combined_adata.obs["patient_id"] = pd.Categorical(combined_adata.obs["patient_id"].astype(str))
    milopy.core.count_cells(combined_adata, sample_col="patient_id")

    # Step 5: Test Differential Abundance with Dataset Batch Covariate
    print("\nFitting QL GLM with technical batch adjustment: ~ dataset_id + clinical_response...")
    nhood_counts_cols = list(combined_adata.uns["nhood_counts"].columns)
    design_df = (
        combined_adata.obs[["patient_id", "dataset_id", "clinical_response"]]
        .drop_duplicates(subset=["patient_id"])
        .set_index("patient_id")
        .loc[nhood_counts_cols]
    )

    n_datasets = len(design_df["dataset_id"].unique())
    design_formula = "~ dataset_id + clinical_response" if n_datasets > 1 else "~ clinical_response"
    print(f"Design Formula : {design_formula}")
    print(f"Patients in GLM : {len(design_df):,} ({sum(design_df['clinical_response'] == 'responder')} R / "
          f"{sum(design_df['clinical_response'] == 'non-responder')} NR)")

    milopy.core.test_nhoods(combined_adata, design=design_formula, design_df=design_df)
    res_df = combined_adata.uns["nhood_test_results"].copy()
    res_df.index.name = "Nhood"

    # Step 6: Vectorized Neighborhood Composition, Dataset Purity & Entropy
    print("\nComputing vectorized neighborhood dataset purity and entropy metrics...")
    ds_comp_res = compute_neighborhood_composition(combined_adata, category_col="dataset_id")
    match ds_comp_res:
        case Failure(err):
            print(f"Dataset composition calculation warning: {err}", file=sys.stderr)
        case Success(ds_comp):
            metrics_pdf = ds_comp.metrics_df.to_pandas().set_index("nhood_index")
            for col in ["dataset_id_purity", "dataset_id_entropy", "dataset_id_simpson_diversity", "dataset_id_active_count"]:
                if col in metrics_pdf.columns:
                    res_df[col] = metrics_pdf[col].values

    if config.group == "pancancer":
        cancer_comp_res = compute_neighborhood_composition(combined_adata, category_col="cancer_type")
        match cancer_comp_res:
            case Success(c_comp):
                c_pdf = c_comp.metrics_df.to_pandas().set_index("nhood_index")
                for col in ["cancer_type_purity", "cancer_type_entropy", "cancer_type_simpson_diversity", "cancer_type_active_count"]:
                    if col in c_pdf.columns:
                        res_df[col] = c_pdf[col].values
            case Failure(err):
                pass

    # Step 6b: Vectorized Rao Score Permutation Testing across Combined Cohorts
    perm_pvals = None
    perm_fdrs = None
    if config.run_permutation:
        print(f"\nRunning vectorized Rao Score permutation test across combined cohorts (B={config.n_permutations})...")
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
                print(f"Permutation calibration complete: {n_perm_sig} nhoods pass permutation FDR < {config.fdr_threshold}.")
            case Failure(err):
                print(f"Warning: permutation test skipped ({err})")

    # Step 6c: Replicate Prevalence and Multi-Cohort Recurrence Filtering
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
            require_permutation_significance=(perm_fdrs is not None),
        ),
    )
    match prev_res:
        case Success(annotated_df):
            res_df = annotated_df
        case Failure(err):
            print(f"Prevalence evaluation warning: {err}")

    n_up = int(((res_df["is_significant"]) & (res_df["logFC"] > 0)).sum())
    n_down = int(((res_df["is_significant"]) & (res_df["logFC"] < 0)).sum())
    n_private_spikes = int((res_df["status"] == "Private Clonal Spike").sum())
    n_conserved = int((res_df["hit_type"] == "Conserved Recurrent Hit").sum())
    print(
        f"\n[RESULTS] Total Nhoods: {len(res_df):,} | Sig Up (Resp): {n_up:,} | "
        f"Sig Down (NR): {n_down:,} | Conserved Hits: {n_conserved:,} | Private Spikes Blocked: {n_private_spikes:,}"
    )

    # Step 7: Export Parquet & Figures
    group_out_dir = config.out_dir / config.group
    group_out_dir.mkdir(parents=True, exist_ok=True)
    config.reports_dir.mkdir(parents=True, exist_ok=True)

    da_parquet_path = group_out_dir / "da_results.parquet"
    res_pl = pl.from_pandas(res_df.reset_index())
    res_pl.write_parquet(da_parquet_path)
    print(f"[EXPORTED] {da_parquet_path}")

    # Project cell scores
    nhoods_sparse = combined_adata.obsm["nhoods"].tocsr()
    logfc_vec = res_df["logFC"].values
    fdr_vec = res_df["FDR"].values
    cell_nhood_counts = np.asarray(nhoods_sparse.sum(axis=1)).ravel()
    safe_counts = np.where(cell_nhood_counts > 0, cell_nhood_counts, 1.0)

    cell_logfc = np.asarray(nhoods_sparse.dot(logfc_vec)).ravel() / safe_counts
    cell_logfc = np.where(cell_nhood_counts > 0, cell_logfc, 0.0)

    # Minimum FDR across cell's neighborhoods
    # For sparse matrix, calculate min FDR efficiently
    cell_min_fdr = np.ones(n_total_cells, dtype=np.float32)
    rows, cols = nhoods_sparse.nonzero()
    for r, c in zip(rows, cols):
        if fdr_vec[c] < cell_min_fdr[r]:
            cell_min_fdr[r] = fdr_vec[c]

    umap_coords = combined_adata.obsm["X_umap"]
    cell_df = pl.DataFrame({
        "cell_id": combined_adata.obs_names.tolist(),
        "dataset_id": combined_adata.obs["dataset_id"].astype(str).tolist(),
        "patient_id": combined_adata.obs["patient_id"].astype(str).tolist(),
        "clinical_response": combined_adata.obs["clinical_response"].astype(str).tolist(),
        "cancer_type": combined_adata.obs["cancer_type"].astype(str).tolist(),
        "UMAP1": umap_coords[:, 0].astype(np.float32),
        "UMAP2": umap_coords[:, 1].astype(np.float32),
        "cell_logfc": cell_logfc.astype(np.float32),
        "min_fdr": cell_min_fdr.astype(np.float32),
    })

    cell_scores_path = group_out_dir / "cell_level_scores.parquet"
    cell_df.write_parquet(cell_scores_path)
    print(f"[EXPORTED] {cell_scores_path}")

    # Generate Publication Figures
    render_combined_volcano_and_umap(
        res_pl=res_pl,
        cell_pl=cell_df,
        group=config.group,
        reports_dir=config.reports_dir,
        fdr_thresh=config.fdr_threshold,
    )

    return Success(da_parquet_path)


def render_combined_volcano_and_umap(
    res_pl: pl.DataFrame,
    cell_pl: pl.DataFrame,
    group: str,
    reports_dir: Path,
    fdr_thresh: float = 0.10,
) -> None:
    """Renders standalone publication-grade volcano and UMAP figures for the combined cohort."""
    pdf = res_pl.to_pandas()
    pdf["neg_log10_fdr"] = -np.log10(pdf["FDR"].clip(lower=1e-300))

    # 1. Volcano Plot
    title_group = group.upper()
    volcano = (
        alt.Chart(pdf)
        .mark_circle(size=18, opacity=0.75)
        .encode(
            x=alt.X("logFC:Q", title="Milo log2 Fold Change (Responder vs Non-Responder)"),
            y=alt.Y("neg_log10_fdr:Q", title="-log10 FDR"),
            color=alt.Color(
                "status:N",
                scale=alt.Scale(
                    domain=["Enriched in Responders", "Enriched in Non-Responders", "Not Significant"],
                    range=[OKABE_VERMILION, OKABE_BLUE, COLOR_NEUTRAL_GREY],
                ),
            ),
            tooltip=["Nhood:N", "logFC:Q", "FDR:Q", "status:N"],
        )
        .properties(
            title=alt.TitleParams(
                text=f"Differential Abundance Volcano Plot: {title_group} Harmonized Cohort",
                subtitle=f"Quasi-Likelihood Negative Binomial GLM (~ dataset_id + clinical_response; FDR < {fdr_thresh})",
                fontSize=14,
                fontWeight="bold",
                subtitleFontSize=11,
                anchor="start",
            ),
            width=500,
            height=380,
        )
    )

    rule = (
        alt.Chart(pd.DataFrame([{"y": -np.log10(fdr_thresh)}]))
        .mark_rule(strokeDash=[4, 4], color=COLOR_HAIRLINE, strokeWidth=1.0)
        .encode(y="y:Q")
    )

    volcano_chart = (volcano + rule).configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    volcano_svg_path = reports_dir / f"milopy_combined_{group}_volcano.svg"
    volcano_svg_str = vlc.vegalite_to_svg(volcano_chart.to_json())
    volcano_svg_path.write_text(volcano_svg_str, encoding="utf-8")
    vlc.svg_to_png(volcano_svg_str, scale=2.0)
    volcano_png_path = volcano_svg_path.with_suffix(".png")
    volcano_png_path.write_bytes(vlc.svg_to_png(volcano_svg_str, scale=2.0))
    print(f"[EXPORTED] {volcano_svg_path} & {volcano_png_path}")

    # 2. UMAP logFC Plot (Clipped at +/- 5)
    cell_pdf = cell_pl.to_pandas()
    if len(cell_pdf) > 20000:
        cell_pdf = cell_pdf.sample(n=20000, random_state=42)

    umap_chart = (
        alt.Chart(cell_pdf)
        .mark_circle(size=10, opacity=0.80)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color(
                "cell_logfc:Q",
                title="Milo log2FC (Responder vs Non-Responder; Clamped at ±5)",
                scale=alt.Scale(
                    domain=[-5.0, 0.0, 5.0],
                    range=[OKABE_BLUE, COLOR_LIGHT_GREY, OKABE_VERMILION],
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    titleLimit=0,
                    gradientLength=300,
                    values=[-5.0, -2.5, 0.0, 2.5, 5.0],
                ),
            ),
        )
        .properties(
            title=alt.TitleParams(
                text=f"Single-Cell Response Landscape: {title_group} Harmonized Cohort",
                subtitle="2D UMAP Embedding Colored by Continuous Milo log2 Fold Change",
                fontSize=14,
                fontWeight="bold",
                subtitleFontSize=11,
                anchor="start",
            ),
            width=500,
            height=420,
        )
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )

    umap_svg_path = reports_dir / f"milopy_combined_{group}_umap_logfc.svg"
    umap_svg_str = vlc.vegalite_to_svg(umap_chart.to_json())
    umap_svg_path.write_text(umap_svg_str, encoding="utf-8")
    umap_png_path = umap_svg_path.with_suffix(".png")
    umap_png_path.write_bytes(vlc.svg_to_png(umap_svg_str, scale=2.0))
    print(f"[EXPORTED] {umap_svg_path} & {umap_png_path}")


def main() -> None:
    config = parse_args()
    res = run_combined_analysis(config)
    match res:
        case Failure(err):
            print(f"[FATAL ERROR] {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"\n[COMPLETED] Combined analysis finished successfully. Results at: {path}")


if __name__ == "__main__":
    main()
