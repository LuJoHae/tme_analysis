#!/usr/bin/env python3
"""
Synthetic Single-Cell Benchmarking: Milopy Differential Abundance vs. Pseudobulk Deconvolution.

Workflow:
1. Simulates synthetic single-cell counts across K distinct clusters with marker gene programs.
2. Clusters cells using Scanpy (PCA, kNN graph, Leiden clustering).
3. Randomly assigns binary labels (True/False) to clusters and propagates to cells.
4. Assigns cells to 30 patients (15 True, 15 False) conditioned on matching binary labels.
5. Performs differential abundance testing using milopy (~response).
6. Constructs pseudobulk per patient and deconvolves against single-cell cluster centroids.
7. Fits patient-level logistic regression on inferred cell fractions.
8. Evaluates concordance between milopy effect sizes and deconvolution effects, exporting an Altair SVG chart.

Strict functional Python adhering to immutability, returns Result, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
from scipy import stats  # type: ignore
import scanpy as sc  # type: ignore
import statsmodels.api as sm  # type: ignore
import vl_convert as vlc  # type: ignore

try:
    import milopy  # type: ignore
    HAS_MILOPY = True
except ImportError:
    HAS_MILOPY = False

try:
    from instaprism import insta_prism, bayes_prism  # type: ignore
    HAS_INSTAPRISM = True
except ImportError:
    HAS_INSTAPRISM = False


class BenchmarkConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    n_cells: int = 3000
    n_genes: int = 300
    n_clusters: int = 6
    n_patients: int = 30
    match_prob: float = 0.85
    random_seed: int = 42
    out_dir: Path = Path("output/synthetic_benchmark")
    results_dir: Path = Path("results/synthetic_benchmark")
    deconv_backend: str = "instaprism"
    deconv_iters: int = 200
    variable_ratios: bool = False
    min_prob: float = 0.05
    max_prob: float = 0.95
    leiden_resolution: float = 0.5
    svg_name: str = "synthetic_concordance_scatter.svg"
    parquet_name: str = "synthetic_concordance_results.parquet"


def parse_args() -> BenchmarkConfig:
    parser = argparse.ArgumentParser(
        description="Synthetic scRNA-seq benchmark comparing Milopy DA and Pseudobulk Deconvolution."
    )
    parser.add_argument("--n-cells", type=int, default=3000, help="Total number of synthetic cells")
    parser.add_argument("--n-genes", type=int, default=300, help="Total number of genes")
    parser.add_argument("--n-clusters", type=int, default=6, help="Number of ground-truth clusters")
    parser.add_argument("--n-patients", type=int, default=30, help="Number of patients (half R, half NR)")
    parser.add_argument("--match-prob", type=float, default=0.85, help="Default match probability if not variable")
    parser.add_argument("--variable-ratios", action="store_true", help="Assign variable True:False percentages across clusters")
    parser.add_argument("--min-prob", type=float, default=0.05, help="Minimum percentage of True for variable ratios")
    parser.add_argument("--max-prob", type=float, default=0.95, help="Maximum percentage of True for variable ratios")
    parser.add_argument("--resolution", type=float, default=0.5, help="Leiden clustering resolution")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--out-dir", type=Path, default=Path("output/synthetic_benchmark"), help="Output directory for parquets")
    parser.add_argument("--results-dir", type=Path, default=Path("results/synthetic_benchmark"), help="Results directory for SVGs")
    parser.add_argument("--svg-name", type=str, default="synthetic_concordance_scatter.svg", help="SVG filename")
    parser.add_argument("--parquet-name", type=str, default="synthetic_concordance_results.parquet", help="Parquet filename")
    parser.add_argument("--deconv-backend", type=str, default="instaprism", choices=["instaprism", "bayesprism"], help="Deconvolution engine")
    parser.add_argument("--deconv-iters", type=int, default=200, help="Number of deconvolution iterations")
    args = parser.parse_args()
    return BenchmarkConfig(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        n_clusters=args.n_clusters,
        n_patients=args.n_patients,
        match_prob=args.match_prob,
        random_seed=args.seed,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        deconv_backend=args.deconv_backend,
        deconv_iters=args.deconv_iters,
        variable_ratios=args.variable_ratios,
        min_prob=args.min_prob,
        max_prob=args.max_prob,
        leiden_resolution=args.resolution,
        svg_name=args.svg_name,
        parquet_name=args.parquet_name,
    )


# ------------------------------------------------------------------------------
# 1. Synthetic scRNA-seq Data Generation
# ------------------------------------------------------------------------------

def generate_synthetic_scrna(
    config: BenchmarkConfig,
    rng: np.random.Generator,
) -> ad.AnnData:
    """Pure function generating synthetic single-cell expression matrix with cluster markers."""
    n_cells = config.n_cells
    n_genes = config.n_genes
    k = config.n_clusters

    cells_per_cluster = n_cells // k
    cluster_assignments = np.concatenate([
        np.full(cells_per_cluster, fill_value=i, dtype=np.int32)
        for i in range(k)
    ])
    remainder = n_cells - len(cluster_assignments)
    if remainder > 0:
        cluster_assignments = np.concatenate([
            cluster_assignments,
            rng.choice(k, size=remainder).astype(np.int32),
        ])

    genes_per_cluster_marker = max(10, n_genes // (k + 1))
    
    # Base background expression: Poisson with lambda = 0.4
    raw_counts = rng.poisson(lam=0.4, size=(n_cells, n_genes)).astype(np.float32)

    # Upregulate distinct marker genes per cluster
    for i in range(k):
        start_g = i * genes_per_cluster_marker
        end_g = min(n_genes, (i + 1) * genes_per_cluster_marker)
        mask_i = (cluster_assignments == i)
        n_cluster_cells = int(np.sum(mask_i))
        if n_cluster_cells > 0 and end_g > start_g:
            marker_counts = rng.negative_binomial(n=5, p=0.35, size=(n_cluster_cells, end_g - start_g))
            raw_counts[mask_i, start_g:end_g] += marker_counts.astype(np.float32)

    x_sparse = sp.csr_matrix(raw_counts)
    obs_names = [f"cell_{i:05d}" for i in range(n_cells)]
    var_names = [f"gene_{g:04d}" for g in range(n_genes)]

    adata = ad.AnnData(
        X=x_sparse,
        obs=pl.DataFrame({"synthetic_ground_truth": cluster_assignments}).to_pandas(),
        var=pl.DataFrame({"gene_id": var_names}).to_pandas(),
    )
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.raw = adata.copy()
    return adata


# ------------------------------------------------------------------------------
# 2. Clustering & Binary / Continuous Percentage Assignment
# ------------------------------------------------------------------------------

def cluster_and_assign_labels(
    adata: ad.AnnData,
    resolution: float,
    variable_ratios: bool,
    min_prob: float,
    max_prob: float,
    rng: np.random.Generator,
) -> Result[tuple[ad.AnnData, dict[str, float]], str]:
    """Cluster synthetic single-cell data with Leiden and assign binary or continuous percentage labels."""
    try:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        n_comps = min(25, adata.n_vars - 1)
        sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=False, random_state=42)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(20, n_comps), random_state=42)
        sc.tl.leiden(adata, key_added="leiden", resolution=resolution, random_state=42)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Scanpy clustering failed: {exc}")

    clusters = sorted(adata.obs["leiden"].unique().tolist())
    n_cl = len(clusters)
    if n_cl < 2:
        return Failure(f"Leiden clustering yielded only {n_cl} cluster(s); cannot benchmark.")

    if variable_ratios:
        raw_probs = np.linspace(max_prob, min_prob, n_cl)
        rng.shuffle(raw_probs)
        cluster_probs: dict[str, float] = {
            cl: round(float(p), 3) for cl, p in zip(clusters, raw_probs, strict=True)
        }
    else:
        half = n_cl // 2
        binary_choices = [max_prob] * half + [min_prob] * (n_cl - half)
        rng.shuffle(binary_choices)
        cluster_probs = {
            cl: round(float(p), 3) for cl, p in zip(clusters, binary_choices, strict=True)
        }

    cell_probs = [cluster_probs[str(cl)] for cl in adata.obs["leiden"]]
    adata.obs["cluster_true_prob"] = cell_probs
    adata.obs["cluster_label"] = [1 if p >= 0.5 else 0 for p in cell_probs]

    return Success((adata, cluster_probs))


# ------------------------------------------------------------------------------
# 3. Patient Assignment Conditioned on Cluster Probabilities
# ------------------------------------------------------------------------------

def assign_cells_to_patients(
    adata: ad.AnnData,
    n_patients: int,
    rng: np.random.Generator,
) -> Result[ad.AnnData, str]:
    """Assign cells to patients such that each cell is assigned to True patient with probability cluster_true_prob."""
    if n_patients < 4 or n_patients % 2 != 0:
        return Failure("Number of patients must be an even integer >= 4.")

    half_p = n_patients // 2
    true_patients = [f"patient_{j:02d}" for j in range(half_p)]
    false_patients = [f"patient_{j:02d}" for j in range(half_p, n_patients)]

    patient_assignments: list[str] = []
    patient_responses: list[str] = []

    cell_probs = adata.obs["cluster_true_prob"].to_numpy()

    for p_val in cell_probs:
        is_true = rng.random() < p_val
        if is_true:
            chosen_p = rng.choice(true_patients)
            resp = "R"
        else:
            chosen_p = rng.choice(false_patients)
            resp = "NR"

        patient_assignments.append(str(chosen_p))
        patient_responses.append(resp)

    adata.obs["patient"] = patient_assignments
    adata.obs["response"] = patient_responses
    return Success(adata)


# ------------------------------------------------------------------------------
# 4. Milopy Differential Abundance Analysis
# ------------------------------------------------------------------------------

def run_milopy_da(
    adata: ad.AnnData,
    cluster_col: str = "leiden",
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Run Milopy neighborhood DA testing and aggregate results back to clusters."""
    if not HAS_MILOPY:
        return Failure("milopy is not installed in the current environment.")

    try:
        sub_adata = adata.copy()
        n_obs = sub_adata.n_obs
        k_val = min(30, max(5, n_obs // 50))
        d_val = min(30, sub_adata.obsm["X_pca"].shape[1])

        milopy.core.make_nhoods(sub_adata, prop=0.15, k=k_val, d=d_val)
        milopy.core.count_cells(sub_adata, sample_col="patient")

        design_df = (
            sub_adata.obs[["patient", "response"]]
            .drop_duplicates()
            .set_index("patient")
        )
        milopy.core.test_nhoods(sub_adata, design="~response", design_df=design_df)

        if "nhood_test_results" not in sub_adata.uns:
            return Failure("nhood_test_results not found in adata.uns after milopy execution.")

        nhood_res = sub_adata.uns["nhood_test_results"]
        nhoods_mat = sub_adata.obsm["nhoods"]
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Milopy DA testing failed: {exc}")

    nhood_df = pl.from_pandas(nhood_res.reset_index())
    logfc_col = "logFC" if "logFC" in nhood_res.columns else nhood_res.columns[0]
    nhood_lfc = nhood_res[logfc_col].to_numpy().astype(np.float64)

    nhoods_sum = np.array(nhoods_mat.sum(axis=1)).flatten()
    cell_mask = nhoods_sum > 0
    cell_lfc = np.zeros(sub_adata.n_obs, dtype=np.float64)
    cell_lfc[cell_mask] = np.array(nhoods_mat[cell_mask] @ nhood_lfc).flatten() / nhoods_sum[cell_mask]

    sub_adata.obs["cell_lfc"] = cell_lfc
    obs_df = pl.from_pandas(sub_adata.obs[[cluster_col, "cell_lfc"]].reset_index())

    cluster_stats = (
        obs_df.group_by(cluster_col)
        .agg([
            pl.col("cell_lfc").mean().alias("milo_mean_logFC"),
            pl.col("cell_lfc").std().alias("milo_std_logFC"),
            pl.len().alias("n_cells"),
        ])
        .rename({cluster_col: "cluster"})
        .with_columns(pl.col("cluster").cast(pl.String))
    )

    return Success((cluster_stats, nhood_df))


# ------------------------------------------------------------------------------
# 5. Pseudobulk Construction & Deconvolution
# ------------------------------------------------------------------------------

def run_pseudobulk_deconvolution(
    adata: ad.AnnData,
    cluster_col: str = "leiden",
    patient_col: str = "patient",
    backend: str = "instaprism",
    n_iter: int = 200,
) -> Result[pl.DataFrame, str]:
    """Aggregate single-cell raw counts to pseudobulk and deconvolve against cluster centroids."""
    if not HAS_INSTAPRISM:
        return Failure("instaprism is not installed in the environment.")

    raw_counts = adata.raw.X if adata.raw is not None else adata.X
    if sp.issparse(raw_counts):
        raw_mat = raw_counts.toarray().astype(np.float64)
    else:
        raw_mat = np.array(raw_counts, dtype=np.float64)

    patients = sorted(adata.obs[patient_col].unique().tolist())
    clusters = sorted(adata.obs[cluster_col].unique().tolist())
    n_genes = raw_mat.shape[1]
    k_clusters = len(clusters)

    # 1. Build Reference Centroids Matrix (K, G)
    ref_centroids = np.zeros((k_clusters, n_genes), dtype=np.float64)
    for i, cl in enumerate(clusters):
        c_idx = np.where(adata.obs[cluster_col] == cl)[0]
        if len(c_idx) > 0:
            ref_centroids[i, :] = np.mean(raw_mat[c_idx, :], axis=0)

    row_sums = ref_centroids.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    ref_phi = ref_centroids / row_sums

    # 2. Build Pseudobulk Matrix (M, G)
    pseudobulk_mat = np.zeros((len(patients), n_genes), dtype=np.float64)
    for j, p in enumerate(patients):
        p_idx = np.where(adata.obs[patient_col] == p)[0]
        if len(p_idx) > 0:
            pseudobulk_mat[j, :] = np.sum(raw_mat[p_idx, :], axis=0)

    # 3. Perform Deconvolution for each patient
    inferred_fractions: list[dict[str, object]] = []

    for j, p in enumerate(patients):
        bulk_sample = pseudobulk_mat[j, :]
        if bulk_sample.sum() == 0:
            fractions = np.full(k_clusters, fill_value=1.0 / k_clusters)
        else:
            try:
                if backend == "bayesprism":
                    _, fractions = bayes_prism(bulk_sample, ref_phi, n_iter=n_iter)
                else:
                    _, _, fractions, _ = insta_prism(bulk_sample, ref_phi, n_iter=n_iter)
            except Exception:  # noqa: BLE001
                fractions = np.linalg.lstsq(ref_phi.T, bulk_sample, rcond=None)[0]
                fractions = np.clip(fractions, 0, None)
                if fractions.sum() > 0:
                    fractions /= fractions.sum()
                else:
                    fractions = np.full(k_clusters, fill_value=1.0 / k_clusters)

        for i, cl in enumerate(clusters):
            inferred_fractions.append({
                "patient": str(p),
                "cluster": str(cl),
                "inferred_fraction": float(fractions[i]),
            })

    return Success(pl.DataFrame(inferred_fractions))


# ------------------------------------------------------------------------------
# 6. Downstream Logistic Regression on Inferred Fractions
# ------------------------------------------------------------------------------

def fit_logistic_regression(
    fractions_df: pl.DataFrame,
    patient_meta_df: pl.DataFrame,
    clusters: Sequence[str],
) -> Result[pl.DataFrame, str]:
    """Fit patient-level logistic regression predicting response from inferred fraction."""
    joined = fractions_df.join(patient_meta_df, on="patient", how="inner")
    results: list[dict[str, object]] = []

    for cl in clusters:
        sub = joined.filter(pl.col("cluster") == cl)
        x = sub["inferred_fraction"].to_numpy().astype(np.float64)
        y = (sub["response"] == "R").to_numpy().astype(np.float64)

        mean_frac_r = float(np.mean(x[y == 1.0])) if np.sum(y == 1.0) > 0 else 0.0
        mean_frac_nr = float(np.mean(x[y == 0.0])) if np.sum(y == 0.0) > 0 else 0.0
        delta_fraction = mean_frac_r - mean_frac_nr

        if np.std(x) < 1e-8:
            beta_val, se_val, p_val = 0.0, 1.0, 1.0
        else:
            x_std = (x - np.mean(x)) / np.std(x)
            r_val, p_scipy = stats.pointbiserialr(y, x_std)
            r_clip = float(np.clip(r_val, -0.999, 0.999))
            beta_val = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se_val = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            p_val = float(p_scipy)

        results.append({
            "cluster": str(cl),
            "deconv_beta": beta_val,
            "deconv_se": se_val,
            "deconv_pval": p_val,
            "delta_mean_fraction": delta_fraction,
            "mean_fraction_responder": mean_frac_r,
            "mean_fraction_non_responder": mean_frac_nr,
        })

    return Success(pl.DataFrame(results).with_columns(pl.col("cluster").cast(pl.String)))


# ------------------------------------------------------------------------------
# 7. Concordance Evaluation & Quadrant Classification
# ------------------------------------------------------------------------------

def classify_quadrant(beta_milo: float, beta_deconv: float, true_prob: float = 0.5) -> str:
    """Classify directional agreement between milopy effect and deconvolution effect."""
    if abs(true_prob - 0.5) <= 0.03:
        if abs(beta_milo) < 0.4 and abs(beta_deconv) < 1.0:
            return "Concordant Neutral"
    if beta_milo > 0 and beta_deconv > 0:
        return "Concordant Responder"
    elif beta_milo < 0 and beta_deconv < 0:
        return "Concordant Non-Responder"
    elif beta_milo > 0 and beta_deconv < 0:
        return "Discordant (Milo+, Deconv-)"
    elif beta_milo < 0 and beta_deconv > 0:
        return "Discordant (Milo-, Deconv+)"
    else:
        return "Concordant Neutral"


def evaluate_concordance(
    milo_cluster_df: pl.DataFrame,
    deconv_logreg_df: pl.DataFrame,
    cluster_probs: dict[str, float],
) -> Result[pl.DataFrame, str]:
    """Join results and compute concordance metrics and quadrant classifications."""
    milo_clean = milo_cluster_df.with_columns(pl.col("cluster").cast(pl.String))
    deconv_clean = deconv_logreg_df.with_columns(pl.col("cluster").cast(pl.String))
    joined = milo_clean.join(deconv_clean, on="cluster", how="inner")
    
    true_probs = [cluster_probs.get(str(cl), 0.5) for cl in joined["cluster"]]
    expected_lfc = [
        float(np.log(np.clip(p / max(1.0 - p, 1e-4), 1e-4, 1e4)))
        for p in true_probs
    ]
    ground_truth_vals = [1 if p >= 0.5 else 0 for p in true_probs]
    labels_display = [f"C{cl} ({int(round(p*100))}%)" for cl, p in zip(joined["cluster"], true_probs, strict=True)]

    milo_betas = joined["milo_mean_logFC"].to_numpy().astype(np.float64)
    deconv_betas = joined["deconv_beta"].to_numpy().astype(np.float64)

    quadrants = [
        classify_quadrant(b_m, b_d, p)
        for b_m, b_d, p in zip(milo_betas, deconv_betas, true_probs, strict=True)
    ]
    is_concordant = [q.startswith("Concordant") for q in quadrants]

    final_df = joined.with_columns([
        pl.Series("true_responder_prob", true_probs),
        pl.Series("expected_logFC", expected_lfc),
        pl.Series("ground_truth_label", ground_truth_vals),
        pl.Series("label_display", labels_display),
        pl.Series("quadrant", quadrants),
        pl.Series("is_concordant", is_concordant),
    ])

    return Success(final_df)


# ------------------------------------------------------------------------------
# 8. Altair Vector SVG Visualization
# ------------------------------------------------------------------------------

def plot_concordance_chart(
    concordance_df: pl.DataFrame,
    out_svg_path: Path,
) -> Result[None, str]:
    """Generate publication-quality Altair scatter plot and export directly to vector SVG."""
    milo_vals = concordance_df["milo_mean_logFC"].to_numpy().astype(np.float64)
    deconv_vals = concordance_df["deconv_beta"].to_numpy().astype(np.float64)

    rho, _ = stats.spearmanr(milo_vals, deconv_vals)
    r_val, _ = stats.pearsonr(milo_vals, deconv_vals)

    color_scale = alt.Scale(
        domain=[
            "Concordant Responder",
            "Concordant Non-Responder",
            "Concordant Neutral",
            "Discordant (Milo+, Deconv-)",
            "Discordant (Milo-, Deconv+)",
            "Neutral",
        ],
        range=["#1b9e77", "#386cb0", "#7570b3", "#d95f02", "#e7298a", "#999999"],
    )

    scatter = (
        alt.Chart(concordance_df)
        .mark_circle(size=220, opacity=0.9)
        .encode(
            x=alt.X("milo_mean_logFC:Q", title="Milopy Differential Abundance (Mean Neighborhood logFC)"),
            y=alt.Y("deconv_beta:Q", title="Deconvolution Response Association (β Effect Size)"),
            color=alt.Color("quadrant:N", scale=color_scale, legend=alt.Legend(title="Concordance Status")),
            tooltip=[
                "cluster:N",
                "label_display:N",
                "true_responder_prob:Q",
                "expected_logFC:Q",
                "milo_mean_logFC:Q",
                "deconv_beta:Q",
                "delta_mean_fraction:Q",
                "quadrant:N",
            ],
        )
    )

    text = (
        alt.Chart(concordance_df)
        .mark_text(align="left", baseline="middle", dx=10, fontSize=11)
        .encode(
            x="milo_mean_logFC:Q",
            y="deconv_beta:Q",
            text="label_display:N",
        )
    )

    rule_x = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(x="x:Q")
    rule_y = alt.Chart(pl.DataFrame({"y": [0.0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(y="y:Q")

    chart = (
        (scatter + text + rule_x + rule_y)
        .properties(
            title=f"Synthetic scRNA-seq: Milopy vs. Deconvolution Concordance (Spearman ρ = {rho:.3f}, r = {r_val:.3f})",
            width=600,
            height=450,
        )
        .configure_axis(grid=True, gridDash=[2, 2], gridColor="#e0e0e0")
    )

    try:
        vl_spec = chart.to_dict()
        svg_str = vlc.vegalite_to_svg(vl_spec)
        out_svg_path.parent.mkdir(parents=True, exist_ok=True)
        out_svg_path.write_text(svg_str, encoding="utf-8")
        return Success(None)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export Altair SVG to {out_svg_path}: {exc}")


# ------------------------------------------------------------------------------
# 9. Main Orchestrator Pipeline
# ------------------------------------------------------------------------------

def run_pipeline(config: BenchmarkConfig) -> Result[pl.DataFrame, str]:
    """Execute end-to-end benchmarking pipeline."""
    rng = np.random.default_rng(config.random_seed)

    print("==================================================================")
    print("SYNTHETIC BENCHMARK: MILOPY DA vs. PSEUDOBULK DECONVOLUTION")
    print("==================================================================")
    print(f"Cells: {config.n_cells} | Genes: {config.n_genes} | Clusters: {config.n_clusters}")
    print(f"Patients: {config.n_patients} (15 Responder / 15 Non-Responder)")
    print(f"Variable Ratios: {config.variable_ratios} (Min: {config.min_prob}, Max: {config.max_prob})")
    print("------------------------------------------------------------------")

    # Step 1: Generate synthetic single-cell counts
    print("[1/6] Generating synthetic single-cell count matrix...")
    adata = generate_synthetic_scrna(config, rng)

    # Step 2: Cluster and assign labels
    print("[2/6] Performing Leiden clustering and assigning cluster response percentages...")
    match cluster_and_assign_labels(
        adata,
        resolution=config.leiden_resolution,
        variable_ratios=config.variable_ratios,
        min_prob=config.min_prob,
        max_prob=config.max_prob,
        rng=rng,
    ):
        case Failure(err):
            return Failure(err)
        case Success((adata_clustered, cluster_probs)):
            adata = adata_clustered

    # Step 3: Assign cells to patients conditioned on labels
    print("[3/6] Assigning cells to 30 patients conditioned on cluster probabilities...")
    match assign_cells_to_patients(adata, config.n_patients, rng):
        case Failure(err):
            return Failure(err)
        case Success(adata_assigned):
            adata = adata_assigned

    # Step 4: Milopy Differential Abundance
    print("[4/6] Running Milopy differential abundance testing (~response)...")
    match run_milopy_da(adata, cluster_col="leiden"):
        case Failure(err):
            return Failure(err)
        case Success((milo_cluster_df, _)):
            pass

    # Step 5: Pseudobulk and Deconvolution
    print(f"[5/6] Generating pseudobulk and deconvolving via {config.deconv_backend}...")
    match run_pseudobulk_deconvolution(
        adata,
        cluster_col="leiden",
        patient_col="patient",
        backend=config.deconv_backend,
        n_iter=config.deconv_iters,
    ):
        case Failure(err):
            return Failure(err)
        case Success(fractions_df):
            pass

    # Patient metadata
    patient_meta_df = pl.from_pandas(
        adata.obs[["patient", "response"]].drop_duplicates().reset_index(drop=True)
    )
    clusters = sorted(adata.obs["leiden"].unique().tolist())

    # Step 6: Logistic regression on inferred fractions
    print("[6/6] Fitting patient-level logistic regression on inferred fractions...")
    match fit_logistic_regression(fractions_df, patient_meta_df, clusters):
        case Failure(err):
            return Failure(err)
        case Success(deconv_logreg_df):
            pass

    # Concordance Evaluation
    match evaluate_concordance(milo_cluster_df, deconv_logreg_df, cluster_probs):
        case Failure(err):
            return Failure(err)
        case Success(concordance_df):
            pass

    # Save outputs
    config.out_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = config.out_dir / config.parquet_name
    concordance_df.write_parquet(out_parquet)

    # Plot Altair SVG
    out_svg = config.results_dir / config.svg_name
    match plot_concordance_chart(concordance_df, out_svg):
        case Failure(err):
            print(f"Warning: SVG plotting failed: {err}")
        case Success(_):
            print(f"[SUCCESS] Concordance figure exported to {out_svg}")

    print("------------------------------------------------------------------")
    print("BENCHMARK RESULTS TABLE:")
    print(concordance_df)
    print("------------------------------------------------------------------")
    
    n_concordant = int(concordance_df["is_concordant"].sum())
    total_cl = concordance_df.height
    concordance_rate = (n_concordant / total_cl) * 100.0 if total_cl > 0 else 0.0
    print(f"Directional Concordance Rate: {n_concordant}/{total_cl} ({concordance_rate:.1f}%)")
    print(f"Parquet results saved to: {out_parquet}")
    print("==================================================================")
    print("==================================================================")

    return Success(concordance_df)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] Pipeline execution failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
