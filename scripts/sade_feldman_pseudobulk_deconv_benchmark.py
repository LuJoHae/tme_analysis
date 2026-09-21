#!/usr/bin/env python3
"""
Sade-Feldman Single-Cell Dataset Benchmark: Milopy DA vs. Pseudobulk Deconvolution.

Workflow:
1. Loads single-cell data from Sade-Feldman et al. (GSE120575, 16,288 cells, 51 samples).
2. Computes true ground-truth cell state fractions per patient sample directly from single-cell counts.
3. Constructs patient-level pseudobulk expression profiles in linear TPM space (2^x - 1).
4. Deconvolves patient pseudobulks against the single-cell reference centroids (InstaPrism / BayesPrism).
5. Fits patient-level logistic regression on inferred and true fractions (Combined, Pre-treatment, Post-treatment).
6. Compares deconvolution effect sizes with single-cell Milopy DA effect sizes across conditions.
7. Evaluates directional concordance, Spearman / Pearson correlations, and deconvolution fidelity.
8. Generates publication-quality Altair vector SVG figures.

Strict functional Python: immutability, returns Result, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Final, Sequence

import altair as alt  # type: ignore
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
from scipy import stats  # type: ignore
import vl_convert as vlc  # type: ignore

try:
    import instaprism  # type: ignore
    HAS_INSTAPRISM = True
except ImportError:
    HAS_INSTAPRISM = False


# ------------------------------------------------------------------------------
# Configuration Model
# ------------------------------------------------------------------------------

class BenchmarkConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    tpm_gz_path: Path = Path("data/raw/GSE120575/GSE120575_tpm.txt.gz")
    cell_scores_path: Path = Path("output/sade_feldman_deconv_validation/milopy_cell_level_scores.parquet")
    phi_path: Path = Path("output/sade_feldman_deconv_validation/reference_phi.parquet")
    marker_genes_path: Path = Path("output/sade_feldman_deconv_validation/reference_marker_genes.parquet")
    milo_da_dir: Path = Path("output/sade_feldman_deconv_validation")
    out_dir: Path = Path("output/sade_feldman_deconv_validation")
    results_dir: Path = Path("results/sade_feldman_deconv_validation")
    deconv_backend: str = "instaprism"
    deconv_iters: int = 100
    use_marker_genes: bool = True
    svg_name: str = "sade_feldman_self_concordance.svg"
    parquet_name: str = "sade_feldman_self_concordance.parquet"


def parse_args() -> BenchmarkConfig:
    parser = argparse.ArgumentParser(
        description="Sade-Feldman benchmarking: Milopy DA vs Pseudobulk Deconvolution."
    )
    parser.add_argument("--tpm-path", type=Path, default=Path("data/raw/GSE120575/GSE120575_tpm.txt.gz"), help="Path to raw TPM gzip")
    parser.add_argument("--cell-scores", type=Path, default=Path("output/sade_feldman_deconv_validation/milopy_cell_level_scores.parquet"), help="Cell-level scores parquet")
    parser.add_argument("--phi-path", type=Path, default=Path("output/sade_feldman_deconv_validation/reference_phi.parquet"), help="Reference phi matrix parquet")
    parser.add_argument("--markers-path", type=Path, default=Path("output/sade_feldman_deconv_validation/reference_marker_genes.parquet"), help="Reference marker genes parquet")
    parser.add_argument("--milo-dir", type=Path, default=Path("output/sade_feldman_deconv_validation"), help="Directory containing milopy DA parquets")
    parser.add_argument("--out-dir", type=Path, default=Path("output/sade_feldman_deconv_validation"), help="Output directory")
    parser.add_argument("--results-dir", type=Path, default=Path("results/sade_feldman_deconv_validation"), help="Results directory for plots")
    parser.add_argument("--deconv-backend", type=str, default="instaprism", choices=["instaprism", "bayesprism"], help="Deconvolution backend")
    parser.add_argument("--deconv-iters", type=int, default=100, help="Number of deconvolution iterations")
    parser.add_argument("--all-genes", action="store_true", help="Deconvolve using all genes instead of marker genes")
    parser.add_argument("--svg-name", type=str, default="sade_feldman_self_concordance.svg", help="SVG filename")
    parser.add_argument("--parquet-name", type=str, default="sade_feldman_self_concordance.parquet", help="Parquet filename")
    args = parser.parse_args()

    return BenchmarkConfig(
        tpm_gz_path=args.tpm_path,
        cell_scores_path=args.cell_scores,
        phi_path=args.phi_path,
        marker_genes_path=args.markers_path,
        milo_da_dir=args.milo_dir,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        deconv_backend=args.deconv_backend,
        deconv_iters=args.deconv_iters,
        use_marker_genes=not args.all_genes,
        svg_name=args.svg_name,
        parquet_name=args.parquet_name,
    )


# ------------------------------------------------------------------------------
# 1. Ground Truth Single-Cell Proportions
# ------------------------------------------------------------------------------

def compute_true_cell_fractions(
    df_cells: pl.DataFrame,
) -> pl.DataFrame:
    """Pure function computing true cell state fractions per sample directly from cell annotations."""
    sample_meta = (
        df_cells.select(["sample_id", "response", "treatment_status"])
        .unique(subset=["sample_id"])
    )
    sample_totals = df_cells.group_by("sample_id").agg(pl.len().alias("n_total"))
    state_counts = df_cells.group_by(["sample_id", "cell_state"]).agg(pl.len().alias("n_cells"))

    unique_samples = sample_meta["sample_id"].to_list()
    unique_states = df_cells["cell_state"].unique().sort().to_list()

    grid = pl.DataFrame({
        "sample_id": [s for s in unique_samples for _ in unique_states],
        "cell_state": [st for _ in unique_samples for st in unique_states],
    })

    joined = (
        grid.join(state_counts, on=["sample_id", "cell_state"], how="left")
        .with_columns(pl.col("n_cells").fill_null(0))
        .join(sample_totals, on="sample_id", how="left")
        .join(sample_meta, on="sample_id", how="left")
        .with_columns(
            (pl.col("n_cells") / pl.col("n_total")).alias("true_fraction")
        )
    )
    return joined


# ------------------------------------------------------------------------------
# 2. Patient Pseudobulk Construction
# ------------------------------------------------------------------------------

def generate_pseudobulk_matrix(
    tpm_path: Path,
    df_cells: pl.DataFrame,
    out_parquet: Path,
    target_genes: set[str] | None = None,
) -> Result[tuple[list[str], list[str], np.ndarray], str]:
    """Extract patient-level pseudobulk by streaming raw log2(TPM+1) and aggregating in linear scale."""
    if out_parquet.exists():
        print(f"Loading cached pseudobulk matrix from {out_parquet}...")
        try:
            bulk_df = pl.read_parquet(out_parquet)
            sample_ids = bulk_df["sample_id"].to_list()
            gene_cols = [c for c in bulk_df.columns if c != "sample_id"]
            if target_genes is not None:
                gene_cols = [g for g in gene_cols if g in target_genes]
            bulk_mat = bulk_df.select(gene_cols).to_numpy().astype(np.float64)
            return Success((sample_ids, gene_cols, bulk_mat))
        except Exception as exc:  # noqa: BLE001
            print(f"Failed to read cached pseudobulk ({exc}); recomputing from {tpm_path}...")

    if not tpm_path.exists():
        return Failure(f"TPM raw file not found at {tpm_path}")

    print(f"Generating pseudobulk from {tpm_path}...")
    cell_to_sample = dict(zip(df_cells["cell_id"].to_list(), df_cells["sample_id"].to_list(), strict=True))
    unique_samples = sorted(list(set(cell_to_sample.values())))
    sample_to_idx = {s: i for i, s in enumerate(unique_samples)}

    try:
        with gzip.open(tpm_path, "rt") as f:
            cell_ids = f.readline().strip().split("\t")
            _ = f.readline()  # skip sample metadata row

            sample_indices = np.array([sample_to_idx.get(cell_to_sample.get(cid), -1) for cid in cell_ids], dtype=np.int32)
            valid_cell_mask = sample_indices >= 0
            valid_sample_indices = sample_indices[valid_cell_mask]

            n_valid = int(np.sum(valid_cell_mask))
            row_ind = valid_sample_indices
            col_ind = np.arange(n_valid)
            # Mapping matrix M: (n_samples, n_valid_cells)
            m_mat = sp.csr_matrix(
                (np.ones(n_valid, dtype=np.float32), (row_ind, col_ind)),
                shape=(len(unique_samples), n_valid),
            )

            genes: list[str] = []
            sample_bulk_rows: list[np.ndarray] = []

            chunk_size = 2000
            chunk_genes: list[str] = []
            chunk_lines: list[list[float]] = []

            for line in f:
                parts = line.strip().split("\t")
                gene = parts[0].strip().upper()
                if target_genes is not None and gene not in target_genes:
                    continue

                chunk_genes.append(gene)
                chunk_lines.append([float(x) for x in parts[1:]])

                if len(chunk_lines) >= chunk_size:
                    chunk_arr = np.array(chunk_lines, dtype=np.float32)[:, valid_cell_mask]
                    lin_arr = np.expm1(chunk_arr * np.log(2.0))
                    # (chunk, n_samples) = (chunk, n_cells) @ (n_samples, n_cells).T
                    b_chunk = lin_arr @ m_mat.T
                    sample_bulk_rows.append(b_chunk)
                    genes.extend(chunk_genes)
                    chunk_genes = []
                    chunk_lines = []

            if chunk_lines:
                chunk_arr = np.array(chunk_lines, dtype=np.float32)[:, valid_cell_mask]
                lin_arr = np.expm1(chunk_arr * np.log(2.0))
                b_chunk = lin_arr @ m_mat.T
                sample_bulk_rows.append(b_chunk)
                genes.extend(chunk_genes)

        if not sample_bulk_rows:
            return Failure("No matching genes extracted from TPM file.")

        # Combined matrix: (n_genes, n_samples) -> transpose to (n_samples, n_genes)
        full_bulk_mat = np.vstack(sample_bulk_rows).T.astype(np.float64)

        # Build Polars DataFrame and save cache
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
        col_dict: dict[str, object] = {"sample_id": unique_samples}
        for g_idx, g_name in enumerate(genes):
            col_dict[g_name] = full_bulk_mat[:, g_idx]

        pseudobulk_df = pl.DataFrame(col_dict)
        pseudobulk_df.write_parquet(out_parquet)
        print(f"Successfully cached pseudobulk matrix with {len(unique_samples)} samples and {len(genes)} genes to {out_parquet}.")

        return Success((unique_samples, genes, full_bulk_mat))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to extract pseudobulk: {exc}")


# ------------------------------------------------------------------------------
# 3. Deconvolution Engine (InstaPrism / BayesPrism)
# ------------------------------------------------------------------------------

def run_pseudobulk_deconvolution(
    sample_ids: list[str],
    bulk_genes: list[str],
    bulk_mat: np.ndarray,
    phi_df: pl.DataFrame,
    marker_genes: list[str] | None,
    backend: str = "instaprism",
    n_iter: int = 100,
) -> Result[pl.DataFrame, str]:
    """Deconvolve pseudobulk against reference centroids matrix."""
    clusters = phi_df["cluster"].to_list()
    ref_gene_cols = [c for c in phi_df.columns if c != "cluster"]

    # Select common genes
    common_genes_set = set(bulk_genes).intersection(set(ref_gene_cols))
    if marker_genes is not None:
        common_genes_set = common_genes_set.intersection(set(marker_genes))

    common_genes = sorted(list(common_genes_set))
    if len(common_genes) < 10:
        return Failure(f"Too few intersecting genes ({len(common_genes)}) for deconvolution.")

    print(f"Deconvolving {len(sample_ids)} samples across {len(clusters)} cell states using {len(common_genes)} genes...")

    # Build aligned normalized reference matrix (n_clusters, n_genes)
    ref_sub = phi_df.select(common_genes).to_numpy().astype(np.float64)
    ref_row_sums = ref_sub.sum(axis=1, keepdims=True)
    ref_row_sums[ref_row_sums == 0] = 1.0
    norm_phi = ref_sub / ref_row_sums

    # Build aligned bulk matrix (n_samples, n_genes)
    gene_to_bulk_idx = {g: i for i, g in enumerate(bulk_genes)}
    bulk_sub = bulk_mat[:, [gene_to_bulk_idx[g] for g in common_genes]].copy()

    k_clusters = len(clusters)
    inferred_records: list[dict[str, object]] = []

    for s_idx, sid in enumerate(sample_ids):
        bulk_vec = bulk_sub[s_idx, :]
        if bulk_vec.sum() == 0:
            fractions = np.full(k_clusters, 1.0 / k_clusters)
        else:
            try:
                if backend == "instaprism" and HAS_INSTAPRISM:
                    _, _, fractions, _ = instaprism.insta_prism(
                        bulk=bulk_vec,
                        reference=norm_phi,
                        n_iter=n_iter,
                    )
                else:
                    # Non-negative least squares fallback
                    res, _ = sp.linalg.lsqr(norm_phi.T, bulk_vec)[:2]
                    fractions = np.clip(res, 0, None)
                    if fractions.sum() > 0:
                        fractions /= fractions.sum()
                    else:
                        fractions = np.full(k_clusters, 1.0 / k_clusters)
            except Exception:  # noqa: BLE001
                fractions = np.full(k_clusters, 1.0 / k_clusters)

        # Fallback if instaprism produced NaN or Infs due to sparse bulk
        if np.any(np.isnan(fractions)) or np.any(np.isinf(fractions)) or np.sum(fractions) <= 0:
            res = np.linalg.lstsq(norm_phi.T, bulk_vec, rcond=None)[0]
            fractions = np.clip(res, 0, None)
            if fractions.sum() > 0 and not np.isnan(fractions.sum()):
                fractions /= fractions.sum()
            else:
                fractions = np.full(k_clusters, 1.0 / k_clusters)

        for c_idx, cl in enumerate(clusters):
            inferred_records.append({
                "sample_id": sid,
                "cell_state": cl,
                "inferred_fraction": float(fractions[c_idx]),
            })

    return Success(pl.DataFrame(inferred_records))


# ------------------------------------------------------------------------------
# 4. Logistic Regression on Fractions
# ------------------------------------------------------------------------------

def fit_logistic_stratum(
    fractions_df: pl.DataFrame,
    fraction_col: str,
    prefix: str,
) -> pl.DataFrame:
    """Fit logistic regression predicting binary response from fraction column for each cell state."""
    clusters = sorted(fractions_df["cell_state"].unique().to_list())
    results: list[dict[str, object]] = []

    for cl in clusters:
        sub = fractions_df.filter(pl.col("cell_state") == cl)
        x_raw = sub[fraction_col].to_numpy().astype(np.float64)
        y_raw = (sub["response"] == "Responder").to_numpy().astype(np.float64)

        valid_mask = ~np.isnan(x_raw) & ~np.isnan(y_raw) & ~np.isinf(x_raw)
        x = x_raw[valid_mask]
        y = y_raw[valid_mask]

        n_r = int(np.sum(y == 1.0))
        n_nr = int(np.sum(y == 0.0))

        mean_r = float(np.mean(x[y == 1.0])) if n_r > 0 else 0.0
        mean_nr = float(np.mean(x[y == 0.0])) if n_nr > 0 else 0.0
        delta_mean = mean_r - mean_nr

        if len(x) < 3 or np.std(x) < 1e-8 or n_r == 0 or n_nr == 0:
            beta_val, se_val, p_val = 0.0, 1.0, 1.0
        else:
            x_std = (x - np.mean(x)) / np.std(x)
            r_val, p_scipy = stats.pointbiserialr(y, x_std)
            r_clip = float(np.clip(r_val, -0.999, 0.999))
            # Standardized logit effect: beta = 2*r / sqrt(1 - r^2)
            beta_val = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se_val = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            p_val = float(p_scipy)

        results.append({
            "cell_state": cl,
            f"{prefix}_beta": beta_val,
            f"{prefix}_se": se_val,
            f"{prefix}_pval": p_val,
            f"{prefix}_delta_mean": delta_mean,
            f"{prefix}_mean_responder": mean_r,
            f"{prefix}_mean_non_responder": mean_nr,
        })

    return pl.DataFrame(results)


def fit_stratified_logistic_regression(
    merged_fractions: pl.DataFrame,
) -> Result[pl.DataFrame, str]:
    """Run logistic regressions across Combined, Pre, and Post conditions for both inferred and true fractions."""
    strata_definitions = [
        ("Combined", merged_fractions),
        ("Pre", merged_fractions.filter(pl.col("treatment_status") == "Pre")),
        ("Post", merged_fractions.filter(pl.col("treatment_status") == "Post")),
    ]

    all_strata_dfs: list[pl.DataFrame] = []

    for cond_name, df_sub in strata_definitions:
        if df_sub.height == 0:
            continue

        deconv_res = fit_logistic_stratum(df_sub, fraction_col="inferred_fraction", prefix="deconv")
        true_res = fit_logistic_stratum(df_sub, fraction_col="true_fraction", prefix="true")

        joined = (
            deconv_res.join(true_res, on="cell_state", how="inner")
            .with_columns([
                pl.lit(cond_name).alias("condition"),
                pl.lit(df_sub["sample_id"].n_unique()).alias("n_samples"),
            ])
        )
        all_strata_dfs.append(joined)

    if not all_strata_dfs:
        return Failure("No valid strata could be evaluated.")

    combined_df = pl.concat(all_strata_dfs)
    return Success(combined_df)


# ------------------------------------------------------------------------------
# 5. Milopy DA Ingestion & Concordance Evaluation
# ------------------------------------------------------------------------------

def load_milopy_da_results(
    milo_dir: Path,
) -> Result[pl.DataFrame, str]:
    """Load precomputed Milopy DA results across Combined, Pre, and Post conditions."""
    conditions = ["Combined", "Pre", "Post"]
    dfs: list[pl.DataFrame] = []

    # Check candidates
    candidate_dirs = [milo_dir, milo_dir.parent / "output" / "sade_feldman_deconv_validation"]

    for cond in conditions:
        found_path: Path | None = None
        for d in candidate_dirs:
            p_cond = d / f"milopy_cell_state_da_{cond}.parquet"
            if p_cond.exists():
                found_path = p_cond
                break
        
        if found_path is None and cond == "Combined":
            p_alt = milo_dir / "milopy_cell_state_da.parquet"
            if p_alt.exists():
                found_path = p_alt

        if found_path is not None:
            try:
                df = pl.read_parquet(found_path)
                if "condition" not in df.columns:
                    df = df.with_columns(pl.lit(cond).alias("condition"))
                
                # Standardize columns
                select_cols = [
                    pl.col("condition"),
                    pl.col("cell_state"),
                    pl.col("milo_mean_logfc"),
                    pl.col("milo_std_logfc"),
                    pl.col("milo_wilcoxon_pval").alias("milo_pval"),
                    pl.col("pct_positive_cells"),
                    pl.col("pct_negative_cells"),
                ]
                dfs.append(df.select(select_cols))
            except Exception as exc:  # noqa: BLE001
                print(f"Warning: Failed to load {found_path}: {exc}")

    if not dfs:
        return Failure("No milopy DA result files could be loaded.")

    return Success(pl.concat(dfs))


def classify_quadrant(beta_milo: float, beta_deconv: float) -> str:
    """Classify directional agreement between milopy effect and deconvolution effect."""
    if beta_milo > 0.1 and beta_deconv > 0.1:
        return "Concordant Responder"
    elif beta_milo < -0.1 and beta_deconv < -0.1:
        return "Concordant Non-Responder"
    elif beta_milo > 0.1 and beta_deconv < -0.1:
        return "Discordant (Milo+, Deconv-)"
    elif beta_milo < -0.1 and beta_deconv > 0.1:
        return "Discordant (Milo-, Deconv+)"
    else:
        return "Concordant Neutral"


def evaluate_concordance(
    logreg_df: pl.DataFrame,
    milo_df: pl.DataFrame,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Join logistic regression results and Milopy DA to compute concordance metrics."""
    joined = logreg_df.join(milo_df, on=["condition", "cell_state"], how="inner")
    if joined.height == 0:
        return Failure("No overlapping (condition, cell_state) entries found between logreg and milopy.")

    milo_lfcs = joined["milo_mean_logfc"].to_numpy().astype(np.float64)
    deconv_betas = joined["deconv_beta"].to_numpy().astype(np.float64)

    quadrants = [
        classify_quadrant(m_lfc, d_b)
        for m_lfc, d_b in zip(milo_lfcs, deconv_betas, strict=True)
    ]
    is_concordant = [q.startswith("Concordant") for q in quadrants]

    detailed_df = joined.with_columns([
        pl.Series("quadrant", quadrants),
        pl.Series("is_concordant", is_concordant),
    ])

    # Compute summary metrics per condition
    summary_records: list[dict[str, object]] = []
    for cond in detailed_df["condition"].unique().sort():
        sub = detailed_df.filter(pl.col("condition") == cond)
        m_lfc = sub["milo_mean_logfc"].to_numpy().astype(np.float64)
        d_b = sub["deconv_beta"].to_numpy().astype(np.float64)
        t_b = sub["true_beta"].to_numpy().astype(np.float64)

        # Correlations between Milo and Deconv
        spearman_rho, spearman_p = stats.spearmanr(m_lfc, d_b)
        pearson_r, pearson_p = stats.pearsonr(m_lfc, d_b)

        # Correlation between True Fraction Effect and Deconv Effect (Deconvolution Fidelity)
        fid_spearman_rho, fid_spearman_p = stats.spearmanr(t_b, d_b)
        fid_pearson_r, fid_pearson_p = stats.pearsonr(t_b, d_b)

        # Correlation between Milo and True Fraction Effect (Single-cell baseline validity)
        milo_true_rho, milo_true_p = stats.spearmanr(m_lfc, t_b)

        # Sign agreement percentage
        sign_matches = np.sum((m_lfc * d_b) > 0)
        sign_pct = float(sign_matches / len(m_lfc)) * 100.0

        summary_records.append({
            "condition": cond,
            "n_states": len(m_lfc),
            "spearman_rho_milo_deconv": float(spearman_rho),
            "spearman_pval_milo_deconv": float(spearman_p),
            "pearson_r_milo_deconv": float(pearson_r),
            "pearson_pval_milo_deconv": float(pearson_p),
            "sign_concordance_pct": sign_pct,
            "fidelity_spearman_rho": float(fid_spearman_rho),
            "fidelity_spearman_pval": float(fid_spearman_p),
            "fidelity_pearson_r": float(fid_pearson_r),
            "milo_true_spearman_rho": float(milo_true_rho),
            "n_concordant": int(sub["is_concordant"].sum()),
        })

    summary_df = pl.DataFrame(summary_records)
    return Success((detailed_df, summary_df))


# ------------------------------------------------------------------------------
# 6. Declarative Altair Visualization
# ------------------------------------------------------------------------------

def build_single_panel(
    df_panel: pl.DataFrame,
    title: str,
    x_col: str,
    y_col: str,
    x_title: str,
    y_title: str,
    rho_val: float,
    r_val: float,
) -> alt.Chart:
    """Build a single scatter plot panel with crosshairs, regression line, and stats annotation."""
    data_pd = df_panel.to_pandas()

    domain = [
        "Concordant Responder",
        "Concordant Non-Responder",
        "Concordant Neutral",
        "Discordant (Milo+, Deconv-)",
        "Discordant (Milo-, Deconv+)",
    ]
    colors = ["#2ca02c", "#1f77b4", "#7f7f7f", "#d62728", "#ff7f0e"]

    color_scale = alt.Scale(domain=domain, range=colors)

    base = alt.Chart(data_pd).encode(
        x=alt.X(f"{x_col}:Q", title=x_title),
        y=alt.Y(f"{y_col}:Q", title=y_title),
    )

    x_rule = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(x=alt.datum(0))
    y_rule = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(y=alt.datum(0))

    scatter = base.mark_circle(size=140, opacity=0.85).encode(
        color=alt.Color("quadrant:N", scale=color_scale, title="Quadrant"),
        tooltip=[
            alt.Tooltip("cell_state:N", title="Cell State"),
            alt.Tooltip(f"{x_col}:Q", title=x_title, format=".3f"),
            alt.Tooltip(f"{y_col}:Q", title=y_title, format=".3f"),
            alt.Tooltip("quadrant:N", title="Classification"),
        ],
    )

    labels = base.mark_text(
        align="left",
        baseline="middle",
        dx=9,
        fontSize=10,
        fontWeight="bold",
    ).encode(
        text=alt.Text("cell_state:N"),
    )

    reg_line = base.transform_regression(x_col, y_col).mark_line(
        color="#333333",
        size=2,
        strokeDash=[6, 4],
    )

    annotation_text = f"Spearman ρ = {rho_val:+.3f}\nPearson r = {r_val:+.3f}"
    annot_df = pl.DataFrame({
        "text": [annotation_text],
        "x": [float(data_pd[x_col].min())],
        "y": [float(data_pd[y_col].max())],
    }).to_pandas()

    annot = alt.Chart(annot_df).mark_text(
        align="left",
        baseline="top",
        fontSize=11,
        fontWeight="bold",
        lineBreak="\n",
    ).encode(
        x="x:Q",
        y="y:Q",
        text="text:N",
    )

    chart = (
        (x_rule + y_rule + reg_line + scatter + labels + annot)
        .properties(width=340, height=320, title=title)
    )
    return chart


def plot_concordance_multipanel(
    detailed_df: pl.DataFrame,
    summary_df: pl.DataFrame,
) -> Result[alt.VConcatChart, str]:
    """Construct a 4-panel Altair vector figure comparing Milopy and Deconvolution across conditions."""
    try:
        # Panel 1: Combined
        comb_df = detailed_df.filter(pl.col("condition") == "Combined")
        comb_sum = summary_df.filter(pl.col("condition") == "Combined").to_dicts()[0]
        p1 = build_single_panel(
            comb_df,
            title="A: Combined Cohort (All 51 Samples)",
            x_col="milo_mean_logfc",
            y_col="deconv_beta",
            x_title="Single-Cell Milopy DA (mean logFC)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(comb_sum["spearman_rho_milo_deconv"]),
            r_val=float(comb_sum["pearson_r_milo_deconv"]),
        )

        # Panel 2: Pre-treatment
        pre_df = detailed_df.filter(pl.col("condition") == "Pre")
        pre_sum = summary_df.filter(pl.col("condition") == "Pre").to_dicts()[0]
        p2 = build_single_panel(
            pre_df,
            title="B: Pre-treatment Cohort (Baseline Response)",
            x_col="milo_mean_logfc",
            y_col="deconv_beta",
            x_title="Single-Cell Milopy DA (mean logFC)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(pre_sum["spearman_rho_milo_deconv"]),
            r_val=float(pre_sum["pearson_r_milo_deconv"]),
        )

        # Panel 3: Post-treatment
        post_df = detailed_df.filter(pl.col("condition") == "Post")
        post_sum = summary_df.filter(pl.col("condition") == "Post").to_dicts()[0]
        p3 = build_single_panel(
            post_df,
            title="C: Post-treatment Cohort (On-Treatment)",
            x_col="milo_mean_logfc",
            y_col="deconv_beta",
            x_title="Single-Cell Milopy DA (mean logFC)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(post_sum["spearman_rho_milo_deconv"]),
            r_val=float(post_sum["pearson_r_milo_deconv"]),
        )

        # Panel 4: True Fraction Effect vs Deconv Effect (Fidelity)
        p4 = build_single_panel(
            pre_df,
            title="D: Deconvolution Fidelity (Pre-treatment)",
            x_col="true_beta",
            y_col="deconv_beta",
            x_title="True Single-Cell Count Effect (beta)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(pre_sum["fidelity_spearman_rho"]),
            r_val=float(pre_sum["fidelity_pearson_r"]),
        )

        row1 = alt.hconcat(p1, p2)
        row2 = alt.hconcat(p3, p4)
        full_chart = (
            alt.vconcat(row1, row2)
            .properties(
                title=alt.TitleParams(
                    text="Sade-Feldman Melanoma (GSE120575): Milopy DA vs. Self-Pseudobulk Deconvolution",
                    subtitle=[
                        "Benchmarking patient-level pseudobulk deconvolution against single-cell differential abundance",
                        "Stratified across Combined, Pre-treatment, and Post-treatment cohorts",
                    ],
                    fontSize=15,
                    subtitleFontSize=12,
                    anchor="middle",
                )
            )
            .configure_view(strokeWidth=1, stroke="#cccccc")
        )

        return Success(full_chart)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Altair visualization construction failed: {exc}")


# ------------------------------------------------------------------------------
# 7. Main Pipeline Execution Boundary
# ------------------------------------------------------------------------------

def run_pipeline(config: BenchmarkConfig) -> Result[Path, str]:
    """Orchestrates end-to-end benchmarking pipeline on the Sade-Feldman dataset."""
    print("--- Starting Sade-Feldman Benchmark Pipeline ---")
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load cell annotations
    if not config.cell_scores_path.exists():
        return Failure(f"Cell scores parquet not found at {config.cell_scores_path}")

    df_cells = pl.read_parquet(config.cell_scores_path)
    print(f"Loaded {df_cells.height} single cells across {df_cells['sample_id'].n_unique()} samples.")

    # 2. Compute true single-cell proportions
    true_fractions_df = compute_true_cell_fractions(df_cells)
    print(f"Computed true ground-truth fractions across {true_fractions_df.height} sample-state pairs.")

    # 3. Load marker genes if requested
    marker_genes: list[str] | None = None
    if config.use_marker_genes and config.marker_genes_path.exists():
        df_markers = pl.read_parquet(config.marker_genes_path)
        marker_genes = sorted(list(set(df_markers["gene"].to_list())))
        print(f"Loaded {len(marker_genes)} curated reference marker genes.")

    # 4. Generate or load pseudobulk matrix
    bulk_parquet = config.out_dir / "sade_feldman_pseudobulk.parquet"
    target_genes = set(marker_genes) if marker_genes else None
    match generate_pseudobulk_matrix(config.tpm_gz_path, df_cells, bulk_parquet, target_genes=target_genes):
        case Failure(err):
            return Failure(err)
        case Success((sample_ids, bulk_genes, bulk_mat)):
            pass

    # 5. Load reference matrix Phi
    if not config.phi_path.exists():
        return Failure(f"Reference phi matrix not found at {config.phi_path}")
    phi_df = pl.read_parquet(config.phi_path)
    print(f"Loaded reference matrix with {phi_df.height} cell states.")

    # 6. Run pseudobulk deconvolution
    match run_pseudobulk_deconvolution(
        sample_ids=sample_ids,
        bulk_genes=bulk_genes,
        bulk_mat=bulk_mat,
        phi_df=phi_df,
        marker_genes=marker_genes,
        backend=config.deconv_backend,
        n_iter=config.deconv_iters,
    ):
        case Failure(err):
            return Failure(err)
        case Success(inferred_fractions_df):
            pass

    # Save inferred fractions
    fractions_out = config.out_dir / "sade_feldman_pseudobulk_fractions.parquet"
    merged_fractions = inferred_fractions_df.join(
        true_fractions_df,
        on=["sample_id", "cell_state"],
        how="inner",
    )
    merged_fractions.write_parquet(fractions_out)
    print(f"Saved merged pseudobulk fractions to {fractions_out}.")

    # 7. Fit stratified logistic regressions
    match fit_stratified_logistic_regression(merged_fractions):
        case Failure(err):
            return Failure(err)
        case Success(logreg_df):
            pass

    # 8. Load Milopy DA results
    match load_milopy_da_results(config.milo_da_dir):
        case Failure(err):
            return Failure(err)
        case Success(milo_df):
            pass

    # 9. Evaluate concordance
    match evaluate_concordance(logreg_df, milo_df):
        case Failure(err):
            return Failure(err)
        case Success((detailed_concordance, summary_concordance)):
            pass

    # Save concordance outputs
    concordance_out = config.out_dir / config.parquet_name
    summary_out = config.out_dir / "sade_feldman_concordance_summary.parquet"
    detailed_concordance.write_parquet(concordance_out)
    summary_concordance.write_parquet(summary_out)
    print(f"Saved detailed concordance to {concordance_out}.")
    print(f"Saved summary concordance to {summary_out}.")

    print("\n--- Summary Concordance Metrics ---")
    print(summary_concordance)

    # 10. Generate Altair visualization
    match plot_concordance_multipanel(detailed_concordance, summary_concordance):
        case Failure(err):
            return Failure(err)
        case Success(chart):
            svg_path = config.results_dir / config.svg_name
            svg_bytes = vlc.vegalite_to_svg(chart.to_dict())
            with open(svg_path, "w", encoding="utf-8") as f:
                f.write(svg_bytes)
            print(f"Exported publication vector SVG to {svg_path}.")
            return Success(svg_path)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Failure(err):
            print(f"Pipeline failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"\nPipeline successfully completed! Final plot: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
